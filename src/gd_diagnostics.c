/*
 * gd_diagnostics.c — MCMC convergence diagnostics for Clepsydra (M8).
 *
 * Rank-normalized split-R-hat and bulk-ESS (Vehtari et al. 2021). See the
 * header for the layout contract. Only bounded scratch is allocated, outside
 * any likelihood hot path; degenerate inputs return NAN.
 */
#include "gd_diagnostics.h"

/* The public header stays dependency-free on purpose — every reported quantity
 * is a double with a documented sentinel, so a caller needs nothing but
 * <stddef.h>. Only the internal scratch helper below has a real failure to
 * report, and it reports it the way the rest of the library does.
 *
 * Every allocation here passes a NULL gd_error for the same reason: these
 * functions answer with a number, not a status, so the NAN they return on
 * failure is the whole answer. The checked helpers are still used — the point
 * is the overflow check, not the message. */
#include "gd_error.h"
#include "gd_memory.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------- inverse normal CDF */

double gd_inv_normal_cdf(double p)
{
    if (isnan(p) || p < 0.0 || p > 1.0) return NAN;
    if (p == 0.0) return -HUGE_VAL;
    if (p == 1.0) return HUGE_VAL;

    /* Acklam's rational approximation. */
    static const double a[6] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00 };
    static const double b[5] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01 };
    static const double c[6] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00 };
    static const double d[4] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00 };
    const double plow = 0.02425, phigh = 1.0 - 0.02425;
    double x, q, r;
    if (p < plow) {
        q = sqrt(-2.0 * log(p));
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (p <= phigh) {
        q = p - 0.5;
        r = q * q;
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
    } else {
        q = sqrt(-2.0 * log(1.0 - p));
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    /* One Halley step for full double precision. */
    double e = 0.5 * erfc(-x / sqrt(2.0)) - p;
    double u = e * sqrt(2.0 * M_PI) * exp(0.5 * x * x);
    x = x - u / (1.0 + 0.5 * x * u);
    return x;
}

/* ------------------------------------------------------------- layout */

/*
 * Everything below follows ArviZ 0.23 (arviz/stats/diagnostics.py) line for
 * line, so that the numbers can be golden-tested against it
 * (tests/unit/test_diagnostics.c). The estimators are those of Vehtari,
 * Gelman, Simpson, Carpenter & Bürkner (2021).
 *
 * Splitting: each chain of n draws gives two chains of h = n / 2 draws, the
 * first h and the LAST h (the middle draw of an odd chain is dropped). The
 * split chains are stored chain-major, m = 2 * n_chains of them.
 */
static double *split_chains(const double *samples, size_t n_chains,
                            size_t n_per_chain, size_t *m_out, size_t *h_out)
{
    size_t h = n_per_chain / 2, m = 2 * n_chains, tot;
    if (!gd_size_mul(m, h, &tot)) return NULL;
    double *out = (double *)gd_malloc_array_checked(tot ? tot : 1, sizeof(double), NULL);
    if (out == NULL) return NULL;
    for (size_t j = 0; j < n_chains; j++) {
        const double *c = samples + j * n_per_chain;
        for (size_t i = 0; i < h; i++) {
            out[j * h + i] = c[i];
            out[(n_chains + j) * h + i] = c[n_per_chain - h + i];
        }
    }
    *m_out = m; *h_out = h;
    return out;
}

static bool any_nan(const double *v, size_t n)
{
    for (size_t i = 0; i < n; i++) if (isnan(v[i])) return true;
    return false;
}

/* (value, original index) pair for the rank sort. */
typedef struct { double v; size_t i; } gd_diag_pair;

static int gd_diag_pair_cmp(const void *pa, const void *pb)
{
    const gd_diag_pair *a = (const gd_diag_pair *)pa;
    const gd_diag_pair *b = (const gd_diag_pair *)pb;
    if (a->v < b->v) return -1;
    if (a->v > b->v) return 1;
    return 0;
}

/* Rank-normalise in place: average ranks for ties, then the Blom normal
 * scores Phi^{-1}((r - 3/8) / (N + 1/4)). */
static int z_scale(double *v, size_t total)
{
    gd_diag_pair *pr = (gd_diag_pair *)gd_malloc_array_checked(total ? total : 1,
                                                              sizeof(gd_diag_pair), NULL);
    if (pr == NULL) return GD_E_NOMEM;
    for (size_t i = 0; i < total; i++) { pr[i].v = v[i]; pr[i].i = i; }
    qsort(pr, total, sizeof(gd_diag_pair), gd_diag_pair_cmp);
    for (size_t s = 0; s < total; ) {
        size_t e = s + 1;
        while (e < total && pr[e].v == pr[s].v) e++;
        double avg_rank = ((double)(s + 1) + (double)e) / 2.0;
        double zz = gd_inv_normal_cdf((avg_rank - 0.375) / ((double)total + 0.25));
        for (size_t k = s; k < e; k++) v[pr[k].i] = zz;
        s = e;
    }
    free(pr);
    return 0;
}

static int dbl_cmp(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

/* Type-7 (R default, numpy linear) quantile of n values; sorts a copy. */
static double quantile7(const double *v, size_t n, double p)
{
    double *c = (double *)gd_malloc_array_checked(n ? n : 1, sizeof(double), NULL);
    if (c == NULL) return NAN;
    for (size_t i = 0; i < n; i++) c[i] = v[i];
    qsort(c, n, sizeof(double), dbl_cmp);
    double hpos = (double)(n - 1) * p;
    size_t lo = (size_t)floor(hpos);
    size_t hi = (lo + 1 < n) ? lo + 1 : lo;
    double q = c[lo] + (hpos - (double)lo) * (c[hi] - c[lo]);
    free(c);
    return q;
}

/* The median as numpy computes it: the middle value, or (a + b) / 2 of the
 * two middle values. Interpolating as a + (b - a) / 2 rounds differently in
 * the last bit, which is enough to create or break a tie in the folded draws
 * and move the tail R-hat. */
static double median_of(const double *v, size_t n)
{
    double *c = (double *)gd_malloc_array_checked(n ? n : 1, sizeof(double), NULL);
    if (c == NULL || n == 0) { free(c); return NAN; }
    for (size_t i = 0; i < n; i++) c[i] = v[i];
    qsort(c, n, sizeof(double), dbl_cmp);
    double med = (n % 2) ? c[n / 2] : (c[n / 2 - 1] + c[n / 2]) / 2.0;
    free(c);
    return med;
}

/* ------------------------------------------------------------ R-hat core */

/* Classic R-hat on m chains of n draws (no splitting here). */
static double rhat_core(const double *v, size_t m, size_t n)
{
    if (m < 2 || n < 2 || any_nan(v, m * n)) return NAN;
    double grand = 0.0, W = 0.0;
    double *means = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    if (means == NULL) return NAN;
    for (size_t j = 0; j < m; j++) {
        const double *c = v + j * n;
        double mean = 0.0, var = 0.0;
        for (size_t i = 0; i < n; i++) mean += c[i];
        mean /= (double)n;
        for (size_t i = 0; i < n; i++) { double d = c[i] - mean; var += d * d; }
        W += var / (double)(n - 1);
        means[j] = mean;
        grand += mean;
    }
    W /= (double)m;
    grand /= (double)m;
    double vm = 0.0;
    for (size_t j = 0; j < m; j++) { double d = means[j] - grand; vm += d * d; }
    vm /= (double)(m - 1);
    free(means);
    if (!(W > 0.0)) return NAN;
    double B = (double)n * vm;
    return sqrt((B / W + (double)n - 1.0) / (double)n);
}

/* ------------------------------------------------------------- ESS core */

/* Biased (divide-by-n) autocovariance of every chain at lag t, averaged. */
static double mean_acov(const double *v, const double *means, size_t m,
                        size_t n, size_t t)
{
    double acc = 0.0;
    for (size_t j = 0; j < m; j++) {
        const double *c = v + j * n;
        double mu = means[j], a = 0.0;
        for (size_t i = 0; i + t < n; i++) a += (c[i] - mu) * (c[i + t] - mu);
        acc += a / (double)n;
    }
    return acc / (double)m;
}

/* Multi-chain ESS of m chains of n draws (ArviZ `_ess`). Autocovariances are
 * computed lag by lag, only as far as Geyer's sequence needs them. */
static double ess_core(const double *v, size_t m, size_t n)
{
    size_t tot = m * n;
    if (m < 1 || n < 4 || any_nan(v, tot)) return NAN;
    double vmin = v[0], vmax = v[0];
    for (size_t i = 1; i < tot; i++) { if (v[i] < vmin) vmin = v[i]; if (v[i] > vmax) vmax = v[i]; }
    if (vmax - vmin < 1e-15) return (double)tot;

    double *means = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    double *rho = (double *)gd_calloc_checked(n, sizeof(double), NULL);
    if (means == NULL || rho == NULL) { free(means); free(rho); return NAN; }
    double grand = 0.0;
    for (size_t j = 0; j < m; j++) {
        double mu = 0.0;
        for (size_t i = 0; i < n; i++) mu += v[j * n + i];
        means[j] = mu / (double)n;
        grand += means[j];
    }
    grand /= (double)m;

    double mean_var = mean_acov(v, means, m, n, 0) * (double)n / ((double)n - 1.0);
    double var_plus = mean_var * ((double)n - 1.0) / (double)n;
    if (m > 1) {
        double vm = 0.0;
        for (size_t j = 0; j < m; j++) { double d = means[j] - grand; vm += d * d; }
        var_plus += vm / (double)(m - 1);
    }

    double rho_even = 1.0;
    rho[0] = rho_even;
    double rho_odd = 1.0 - (mean_var - mean_acov(v, means, m, n, 1)) / var_plus;
    rho[1] = rho_odd;

    /* Geyer's initial positive sequence. */
    size_t t = 1;
    while (t + 3 < n && (rho_even + rho_odd) > 0.0) {
        rho_even = 1.0 - (mean_var - mean_acov(v, means, m, n, t + 1)) / var_plus;
        rho_odd  = 1.0 - (mean_var - mean_acov(v, means, m, n, t + 2)) / var_plus;
        if ((rho_even + rho_odd) >= 0.0) {
            rho[t + 1] = rho_even;
            rho[t + 2] = rho_odd;
        }
        t += 2;
    }
    /* ArviZ sets max_t = t - 2, which is -1 when the loop never ran (t == 1).
     * Unsigned arithmetic cannot hold that, so the case is carried as a flag. */
    bool neg_max = (t < 3);
    size_t max_t = neg_max ? 0 : t - 2;
    if (!neg_max && rho_even > 0.0) rho[max_t + 1] = rho_even;
    if (neg_max && rho_even > 0.0) rho[0] = rho_even;   /* index max_t + 1 = 0 */

    /* Geyer's initial monotone sequence. */
    if (!neg_max) {
        for (size_t u = 1; u + 2 <= max_t; u += 2) {
            if ((rho[u + 1] + rho[u + 2]) > (rho[u - 1] + rho[u])) {
                rho[u + 1] = (rho[u - 1] + rho[u]) / 2.0;
                rho[u + 2] = rho[u + 1];
            }
        }
    }

    double sum = 0.0;
    double extra;
    if (neg_max) {
        /* rho[:0] is empty and rho[0:1] is rho[0]. */
        extra = rho[0];
    } else {
        for (size_t u = 0; u <= max_t; u++) sum += rho[u];
        extra = (max_t + 1 < n) ? rho[max_t + 1] : 0.0;
    }
    double ess = (double)tot;
    double tau = -1.0 + 2.0 * sum + extra;
    double floor_tau = 1.0 / log10(ess);
    if (tau < floor_tau) tau = floor_tau;
    bool bad = false;
    for (size_t u = 0; u < n; u++) if (isnan(rho[u])) bad = true;
    free(means); free(rho);
    return bad ? NAN : ess / tau;
}

/* ------------------------------------------------------------ public API */

double gd_split_rhat(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 2 || n_per_chain < 4) return NAN;
    size_t m, h;
    double *sp = split_chains(samples, n_chains, n_per_chain, &m, &h);
    if (sp == NULL) return NAN;
    double r = rhat_core(sp, m, h);
    free(sp);
    return r;
}

double gd_rhat_bulk(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 2 || n_per_chain < 4) return NAN;
    if (any_nan(samples, n_chains * n_per_chain)) return NAN;
    size_t m, h;
    double *sp = split_chains(samples, n_chains, n_per_chain, &m, &h);
    if (sp == NULL) return NAN;
    double r = (z_scale(sp, m * h) == 0) ? rhat_core(sp, m, h) : NAN;
    free(sp);
    return r;
}

double gd_rhat_tail(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 2 || n_per_chain < 4) return NAN;
    if (any_nan(samples, n_chains * n_per_chain)) return NAN;
    size_t m, h;
    double *sp = split_chains(samples, n_chains, n_per_chain, &m, &h);
    if (sp == NULL) return NAN;
    double med = median_of(sp, m * h);
    for (size_t i = 0; i < m * h; i++) sp[i] = fabs(sp[i] - med);
    double r = (z_scale(sp, m * h) == 0) ? rhat_core(sp, m, h) : NAN;
    free(sp);
    return r;
}

double gd_rhat_rank(const double *samples, size_t n_chains, size_t n_per_chain)
{
    double b = gd_rhat_bulk(samples, n_chains, n_per_chain);
    double t = gd_rhat_tail(samples, n_chains, n_per_chain);
    if (isnan(b) || isnan(t)) return NAN;
    return (b > t) ? b : t;
}

double gd_ess_bulk(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 1 || n_per_chain < 4) return NAN;
    if (any_nan(samples, n_chains * n_per_chain)) return NAN;
    size_t m, h;
    double *sp = split_chains(samples, n_chains, n_per_chain, &m, &h);
    if (sp == NULL) return NAN;
    double e = (z_scale(sp, m * h) == 0) ? ess_core(sp, m, h) : NAN;
    free(sp);
    return e;
}

double gd_ess_mean(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 1 || n_per_chain < 4) return NAN;
    size_t m, h;
    double *sp = split_chains(samples, n_chains, n_per_chain, &m, &h);
    if (sp == NULL) return NAN;
    double e = ess_core(sp, m, h);
    free(sp);
    return e;
}

/* ESS of the indicator I(x <= q_p), q_p the type-7 quantile of all draws. */
static double ess_quantile(const double *samples, size_t n_chains,
                           size_t n_per_chain, double p)
{
    size_t tot = n_chains * n_per_chain;
    double q = quantile7(samples, tot, p);
    if (isnan(q)) return NAN;
    double *ind = (double *)gd_malloc_array_checked(tot, sizeof(double), NULL);
    if (ind == NULL) return NAN;
    for (size_t i = 0; i < tot; i++) ind[i] = (samples[i] <= q) ? 1.0 : 0.0;
    double e = gd_ess_mean(ind, n_chains, n_per_chain);
    free(ind);
    return e;
}

double gd_ess_tail(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 1 || n_per_chain < 4) return NAN;
    if (any_nan(samples, n_chains * n_per_chain)) return NAN;
    double lo = ess_quantile(samples, n_chains, n_per_chain, 0.05);
    double hi = ess_quantile(samples, n_chains, n_per_chain, 0.95);
    if (isnan(lo) || isnan(hi)) return NAN;
    return (lo < hi) ? lo : hi;
}

double gd_mcse_mean(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 1 || n_per_chain < 4) return NAN;
    size_t tot = n_chains * n_per_chain;
    if (any_nan(samples, tot)) return NAN;
    double ess = gd_ess_mean(samples, n_chains, n_per_chain);
    double mu = 0.0, var = 0.0;
    for (size_t i = 0; i < tot; i++) mu += samples[i];
    mu /= (double)tot;
    for (size_t i = 0; i < tot; i++) { double d = samples[i] - mu; var += d * d; }
    var /= (double)(tot - 1);
    return sqrt(var) / sqrt(ess);
}

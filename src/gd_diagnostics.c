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

/* ------------------------------------------------------- split-R-hat core */

/* Classic split-R-hat over the chain-major array v (n_chains x n_per_chain).
 * Returns NAN on degenerate input. */
static double split_rhat_array(const double *v, size_t n_chains,
                               size_t n_per_chain)
{
    if (n_chains < 2 || n_per_chain < 4) return NAN;
    size_t n = n_per_chain / 2;          /* half-chain length */
    size_t m = 2 * n_chains;             /* number of half-chains */

    double grand = 0.0;
    double *means = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    double *vars  = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    if (means == NULL || vars == NULL) { free(means); free(vars); return NAN; }

    for (size_t j = 0; j < n_chains; j++) {
        for (int half = 0; half < 2; half++) {
            const double *seg = v + j * n_per_chain + (size_t)half * n;
            double mean = 0.0;
            for (size_t i = 0; i < n; i++) mean += seg[i];
            mean /= (double)n;
            double var = 0.0;
            for (size_t i = 0; i < n; i++) {
                double dd = seg[i] - mean;
                var += dd * dd;
            }
            var /= (double)(n - 1);
            size_t hc = j * 2 + (size_t)half;
            means[hc] = mean;
            vars[hc]  = var;
            grand += mean;
        }
    }
    grand /= (double)m;

    double W = 0.0, Bvar = 0.0;
    for (size_t hc = 0; hc < m; hc++) {
        W += vars[hc];
        double dd = means[hc] - grand;
        Bvar += dd * dd;
    }
    W /= (double)m;
    Bvar /= (double)(m - 1);             /* between-chain variance of means */
    free(means);
    free(vars);
    if (!(W > 0.0)) return NAN;
    double var_plus = (double)(n - 1) / (double)n * W + Bvar;
    return sqrt(var_plus / W);
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

/* Rank-normalize the pooled draws into z (chain-major, same layout). Ties get
 * the average rank; ranks map to normal scores via z = Phi^{-1}((r-3/8)/(N-1/4)). */
static int rank_normalize(const double *samples, size_t total, double *z)
{
    gd_diag_pair *pr = (gd_diag_pair *)gd_malloc_array_checked(total, sizeof(gd_diag_pair), NULL);
    if (pr == NULL) return GD_E_NOMEM;
    for (size_t i = 0; i < total; i++) { pr[i].v = samples[i]; pr[i].i = i; }
    qsort(pr, total, sizeof(gd_diag_pair), gd_diag_pair_cmp);

    for (size_t s = 0; s < total; ) {
        size_t e = s + 1;
        while (e < total && pr[e].v == pr[s].v) e++;
        /* positions s..e-1 (0-based) => 1-based ranks s+1..e; average them. */
        double avg_rank = ((double)(s + 1) + (double)e) / 2.0;
        double zz = gd_inv_normal_cdf((avg_rank - 0.375) /
                                      ((double)total - 0.25));
        for (size_t k = s; k < e; k++) z[pr[k].i] = zz;
        s = e;
    }
    free(pr);
    return 0;
}

double gd_split_rhat(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL) return NAN;
    return split_rhat_array(samples, n_chains, n_per_chain);
}

double gd_rhat_rank(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 2 || n_per_chain < 4) return NAN;
    size_t total = n_chains * n_per_chain;
    double *z = (double *)gd_malloc_array_checked(total, sizeof(double), NULL);
    if (z == NULL) return NAN;
    if (rank_normalize(samples, total, z) != 0) { free(z); return NAN; }
    double r = split_rhat_array(z, n_chains, n_per_chain);
    free(z);
    return r;
}

/* --------------------------------------------------------------- bulk ESS */

/* Multi-chain ESS on the chain-major array v of split half-chains. */
static double ess_array(const double *v, size_t n_chains, size_t n_per_chain)
{
    if (n_chains < 1 || n_per_chain < 4) return NAN;
    size_t n = n_per_chain / 2;
    size_t m = 2 * n_chains;
    if (n < 2) return NAN;

    double *means = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    double *vars  = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    if (means == NULL || vars == NULL) { free(means); free(vars); return NAN; }

    /* Contiguous copy of the half-chains for cache-friendly autocov. The
     * element count is itself a product, so it is multiplied out under a check
     * before the allocation sees it. */
    size_t nseg;
    if (!gd_size_mul(m, n, &nseg)) { free(means); free(vars); return NAN; }
    double *seg = (double *)gd_malloc_array_checked(nseg, sizeof(double), NULL);
    if (seg == NULL) { free(means); free(vars); return NAN; }

    double grand = 0.0;
    for (size_t j = 0; j < n_chains; j++) {
        for (int half = 0; half < 2; half++) {
            size_t hc = j * 2 + (size_t)half;
            const double *src = v + j * n_per_chain + (size_t)half * n;
            double mean = 0.0;
            for (size_t i = 0; i < n; i++) { seg[hc*n+i] = src[i]; mean += src[i]; }
            mean /= (double)n;
            double var = 0.0;
            for (size_t i = 0; i < n; i++) {
                double dd = src[i] - mean;
                var += dd * dd;
            }
            means[hc] = mean;
            vars[hc]  = var / (double)(n - 1);
            grand += mean;
        }
    }
    grand /= (double)m;

    double W = 0.0, Bvar = 0.0;
    for (size_t hc = 0; hc < m; hc++) {
        W += vars[hc];
        double dd = means[hc] - grand;
        Bvar += dd * dd;
    }
    W /= (double)m;
    Bvar /= (double)(m - 1);
    if (!(W > 0.0)) { free(means); free(vars); free(seg); return NAN; }
    double var_plus = (double)(n - 1) / (double)n * W + Bvar;

    /* rho[t] = 1 - (W - mean_chain_autocov(t)) / var_plus, biased autocov. */
    double *rho = (double *)gd_malloc_array_checked(n, sizeof(double), NULL);
    if (rho == NULL) { free(means); free(vars); free(seg); return NAN; }
    rho[0] = 1.0;
    for (size_t t = 1; t < n; t++) {
        double acov = 0.0;
        for (size_t hc = 0; hc < m; hc++) {
            const double *s = seg + hc * n;
            double mean = means[hc];
            double a = 0.0;
            for (size_t i = 0; i + t < n; i++) {
                a += (s[i] - mean) * (s[i + t] - mean);
            }
            acov += a / (double)n;         /* biased (÷ n) */
        }
        acov /= (double)m;
        rho[t] = 1.0 - (W - acov) / var_plus;
    }

    /* Geyer initial monotone positive sequence on paired sums. */
    double tau = 1.0;
    double prev_pair = HUGE_VAL;
    for (size_t k = 1; 2 * k <= n - 1; k++) {
        double pair = rho[2 * k - 1] + rho[2 * k];
        if (pair <= 0.0) break;
        if (pair > prev_pair) pair = prev_pair;   /* enforce monotonicity */
        tau += 2.0 * pair;
        prev_pair = pair;
    }
    if (tau < 1.0) tau = 1.0;

    free(means); free(vars); free(seg); free(rho);
    return (double)(m * n) / tau;
}

double gd_ess_bulk(const double *samples, size_t n_chains, size_t n_per_chain)
{
    if (samples == NULL || n_chains < 1 || n_per_chain < 4) return NAN;
    size_t total = n_chains * n_per_chain;
    double *z = (double *)gd_malloc_array_checked(total, sizeof(double), NULL);
    if (z == NULL) return NAN;
    if (rank_normalize(samples, total, z) != 0) { free(z); return NAN; }
    double e = ess_array(z, n_chains, n_per_chain);
    free(z);
    return e;
}

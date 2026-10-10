/*
 * gd_loo.c — PSIS-LOO and stacking. See gd_loo.h.
 *
 * The smoothing follows the reference algorithm as the loo package and ArviZ
 * implement it: tail size M = ceil(min(S/5, 3 sqrt(S))), threshold at the
 * (M+1)-th largest log ratio, the M exceedances replaced by the fitted GPD's
 * quantiles at (z - 1/2)/M, and truncation at the largest raw ratio.
 */
#include "gd_loo.h"
#include "gd_math.h"
#include "gd_memory.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

double gd_gpd_fit(const double *x, size_t n, double *sigma)
{
    const double prior_bs = 3.0, prior_k = 10.0;
    size_t m = 30 + (size_t)sqrt((double)n);
    double *b = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    double *l = (double *)gd_malloc_array_checked(m, sizeof(double), NULL);
    if (b == NULL || l == NULL || n == 0) {
        free(b); free(l);
        *sigma = NAN;
        return NAN;
    }
    size_t q = (size_t)((double)n / 4.0 + 0.5);
    double xq = x[q > 0 ? q - 1 : 0], xmax = x[n - 1];
    for (size_t j = 0; j < m; j++) {
        b[j] = 1.0 - sqrt((double)m / ((double)j + 0.5));
        b[j] = b[j] / (prior_bs * xq) + 1.0 / xmax;
        double kj = 0.0;
        for (size_t i = 0; i < n; i++) kj += log1p(-b[j] * x[i]);
        kj /= (double)n;
        l[j] = (double)n * (log(-(b[j] / kj)) - kj - 1.0);
    }
    /* weights_j = 1 / sum_t exp(l_t - l_j); negligible ones dropped. */
    double wsum = 0.0, bpost = 0.0;
    for (size_t j = 0; j < m; j++) {
        double s = 0.0;
        for (size_t t = 0; t < m; t++) s += exp(l[t] - l[j]);
        double w = 1.0 / s;
        if (!(w >= 10.0 * DBL_EPSILON)) continue;
        wsum += w;
        bpost += w * b[j];
    }
    bpost /= wsum;
    double k = 0.0;
    for (size_t i = 0; i < n; i++) k += log1p(-bpost * x[i]);
    k /= (double)n;
    *sigma = -k / bpost;
    k = ((double)n * k + prior_k * 0.5) / ((double)n + prior_k);
    free(b); free(l);
    return k;
}

/* Inverse GPD CDF at p for shape k, scale sigma. */
static double gpd_inv(double p, double k, double sigma)
{
    if (fabs(k) < 1e-12) return -sigma * log1p(-p);
    return sigma * expm1(-k * log1p(-p)) / k;
}

typedef struct { double v; size_t i; } vi_pair;

static int vi_cmp(const void *a, const void *b)
{
    double x = ((const vi_pair *)a)->v, y = ((const vi_pair *)b)->v;
    return (x > y) - (x < y);
}

double gd_psis_smooth(double *lw, size_t S)
{
    if (S == 0) return NAN;
    double mx = -HUGE_VAL, mn = HUGE_VAL;
    for (size_t s = 0; s < S; s++) { if (lw[s] > mx) mx = lw[s]; if (lw[s] < mn) mn = lw[s]; }
    for (size_t s = 0; s < S; s++) lw[s] -= mx;
    if (mx - mn < 1e-12) {
        /* Equal ratios (a concept every draw scores the same, such as one
         * the likelihood excludes): uniform weights, and the estimate is
         * exact, so there is no tail to worry about. */
        for (size_t s = 0; s < S; s++) lw[s] = -log((double)S);
        return 0.0;
    }
    double mtail = ceil(fmin(0.2 * (double)S, 3.0 * sqrt((double)S)));
    size_t M = (size_t)mtail;
    double k = HUGE_VAL;
    if (M >= 5 && M < S) {
        vi_pair *v = (vi_pair *)gd_malloc_array_checked(S, sizeof(vi_pair), NULL);
        double *ex = (double *)gd_malloc_array_checked(M, sizeof(double), NULL);
        if (v == NULL || ex == NULL) { free(v); free(ex); return NAN; }
        for (size_t s = 0; s < S; s++) { v[s].v = lw[s]; v[s].i = s; }
        qsort(v, S, sizeof(vi_pair), vi_cmp);
        double cut = v[S - M - 1].v;
        if (cut < log(DBL_MIN)) cut = log(DBL_MIN);
        double ecut = exp(cut);
        /* Exceedances above the cut, ascending: v[S-M .. S-1]. */
        size_t t0 = S - M;
        size_t nt = 0;
        for (size_t j = t0; j < S; j++) if (v[j].v > cut) ex[nt++] = exp(v[j].v) - ecut;
        if (nt > 4) {
            double sigma;
            k = gd_gpd_fit(ex, nt, &sigma);
            if (isfinite(k)) {
                size_t first = S - nt;
                for (size_t z = 0; z < nt; z++) {
                    double p = ((double)z + 0.5) / (double)nt;
                    double sm = log(gpd_inv(p, k, sigma) + ecut);
                    lw[v[first + z].i] = sm > 0.0 ? 0.0 : sm;   /* truncate */
                }
            }
        }
        free(v); free(ex);
    }
    double lse = gd_logsumexp(lw, S);
    for (size_t s = 0; s < S; s++) lw[s] -= lse;
    return k;
}

int gd_psis_loo(const double *ll, size_t S, size_t N, double *elpd_i, double *k_i)
{
    double *lw = (double *)gd_malloc_array_checked(S ? S : 1, sizeof(double), NULL);
    double *t = (double *)gd_malloc_array_checked(S ? S : 1, sizeof(double), NULL);
    if (lw == NULL || t == NULL) { free(lw); free(t); return 1; }
    for (size_t i = 0; i < N; i++) {
        for (size_t s = 0; s < S; s++) lw[s] = -ll[s * N + i];
        double k = gd_psis_smooth(lw, S);
        if (isnan(k)) { free(lw); free(t); return 1; }
        for (size_t s = 0; s < S; s++) t[s] = lw[s] + ll[s * N + i];
        elpd_i[i] = gd_logsumexp(t, S);
        k_i[i] = k;
    }
    free(lw); free(t);
    return 0;
}

size_t gd_stacking_weights(const double *elpd, size_t M, size_t N, double *w)
{
    if (M == 0) return 0;
    /* p[m][i] = exp(elpd - max over m), so the scale cancels per concept. */
    double *p = (double *)gd_malloc_array_checked(M * (N ? N : 1), sizeof(double), NULL);
    double *g = (double *)gd_malloc_array_checked(M, sizeof(double), NULL);
    if (p == NULL || g == NULL) { free(p); free(g); return 0; }
    for (size_t i = 0; i < N; i++) {
        double mx = -HUGE_VAL;
        for (size_t m = 0; m < M; m++) if (elpd[m * N + i] > mx) mx = elpd[m * N + i];
        for (size_t m = 0; m < M; m++) p[m * N + i] = exp(elpd[m * N + i] - mx);
    }
    for (size_t m = 0; m < M; m++) w[m] = 1.0 / (double)M;
    /* EM fixed point for mixture weights with fixed components: every step
     * increases the stacking objective, and the fixed point is its maximum. */
    size_t it;
    for (it = 1; it <= 100000; it++) {
        for (size_t m = 0; m < M; m++) g[m] = 0.0;
        for (size_t i = 0; i < N; i++) {
            double den = 0.0;
            for (size_t m = 0; m < M; m++) den += w[m] * p[m * N + i];
            for (size_t m = 0; m < M; m++) g[m] += w[m] * p[m * N + i] / den;
        }
        double change = 0.0;
        for (size_t m = 0; m < M; m++) {
            double nw = g[m] / (double)N;
            change = fmax(change, fabs(nw - w[m]));
            w[m] = nw;
        }
        if (change < 1e-10) break;
    }
    free(p); free(g);
    return it;
}

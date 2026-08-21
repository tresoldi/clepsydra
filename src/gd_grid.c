/*
 * gd_grid.c — reduce a log-density on a grid to moments, quantiles and
 * intervals.
 *
 * The arithmetic here is P0's and P1's, preserved operation for operation:
 * moments weight each interval's trapezoidal mass by its midpoint time, and
 * quantiles come from a trapezoidal CDF read back by binary search. That
 * choice is not arbitrary — it is what the shipped pairwise model already
 * did, so adopting it moved no number.
 */
#include "gd_grid.h"
#include "gd_memory.h"

#include <math.h>
#include <string.h>

/* Trapezoidal CDF over the grid, normalised so cdf[n-1] == 1. */
static void build_cdf(const double *t, const double *logdens, size_t n,
                      double *cdf)
{
    double total = 0.0;
    cdf[0] = 0.0;
    for (size_t i = 1; i < n; i++) {
        double dt = t[i] - t[i - 1];
        double f0 = exp(logdens[i - 1]);
        double f1 = exp(logdens[i]);
        total += 0.5 * (f0 + f1) * dt;
        cdf[i] = total;
    }
    if (total > 0.0) {
        for (size_t i = 0; i < n; i++) {
            cdf[i] /= total;
        }
    }
}

/* Inverse CDF by binary search, linearly interpolated within the bracket. */
static double inv_cdf(const double *t, const double *cdf, size_t n,
                      double target)
{
    if (target <= 0.0) return t[0];
    if (target >= 1.0) return t[n - 1];
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (cdf[mid] < target) lo = mid;
        else hi = mid;
    }
    double c0 = cdf[lo], c1 = cdf[hi];
    if (c1 == c0) return t[lo];
    double frac = (target - c0) / (c1 - c0);
    return t[lo] + frac * (t[hi] - t[lo]);
}

/* The mass a grid point owns: half the interval on each side of it. */
static double cell_width(const double *t, size_t n, size_t k)
{
    double w = 0.0;
    if (k > 0)     w += 0.5 * (t[k] - t[k - 1]);
    if (k < n - 1) w += 0.5 * (t[k + 1] - t[k]);
    return w;
}

/*
 * Narrowest interval holding `level` of the mass: take grid points in
 * descending density until the mass is reached, and report the span of what
 * was taken. `idx` is scratch of n entries, already sorted by descending
 * density. On a grid this is an approximation, and it is the approximation
 * P0 has always used.
 */
static void hpd_interval(const double *t, const double *logdens, size_t n,
                         const size_t *idx, double level,
                         double *lo, double *hi)
{
    double total_mass = 0.0;
    for (size_t i = 0; i < n; i++) {
        size_t k = idx[i];
        total_mass += exp(logdens[k]) * cell_width(t, n, k);
    }
    if (total_mass <= 0.0) { *lo = t[0]; *hi = t[n - 1]; return; }

    double target = level * total_mass;
    double acc = 0.0;
    double min_t = t[n - 1], max_t = t[0];
    for (size_t i = 0; i < n; i++) {
        size_t k = idx[i];
        acc += exp(logdens[k]) * cell_width(t, n, k);
        if (t[k] < min_t) min_t = t[k];
        if (t[k] > max_t) max_t = t[k];
        if (acc >= target) break;
    }
    *lo = min_t;
    *hi = max_t;
}

int gd_grid_summarize(const double *t, const double *logdens, size_t n,
                      gd_grid_summary *out, gd_error *err)
{
    int rc = GD_E_INTERNAL;
    double *cdf = NULL;
    size_t *idx = NULL;
    gd_grid_summary s;

    if (t == NULL || logdens == NULL || out == NULL || n < 2) {
        return gd_error_set(err, GD_E_INVALID_ARG, "grid: bad argument");
    }

    cdf = (double *)gd_calloc_checked(n, sizeof(double), err);
    if (cdf == NULL) { rc = GD_E_NOMEM; goto done; }
    idx = (size_t *)gd_malloc_array_checked(n, sizeof(size_t), err);
    if (idx == NULL) { rc = GD_E_NOMEM; goto done; }

    memset(&s, 0, sizeof s);

    /* Mean, and the mass it is normalised by. */
    double mean = 0.0, total = 0.0;
    for (size_t i = 1; i < n; i++) {
        double dt = t[i] - t[i - 1];
        double f0 = exp(logdens[i - 1]);
        double f1 = exp(logdens[i]);
        double t_mid = 0.5 * (t[i] + t[i - 1]);
        double mass = 0.5 * (f0 + f1) * dt;
        mean += t_mid * mass;
        total += mass;
    }
    if (total > 0) mean /= total;
    s.mean = mean;

    /* Second moment, in its own pass, so the summation order matches what
     * P0 and P1 did before this module existed. */
    double m2 = 0.0;
    for (size_t i = 1; i < n; i++) {
        double dt = t[i] - t[i - 1];
        double f0 = exp(logdens[i - 1]);
        double f1 = exp(logdens[i]);
        double t_mid = 0.5 * (t[i] + t[i - 1]);
        double mass = 0.5 * (f0 + f1) * dt;
        m2 += t_mid * t_mid * mass;
    }
    if (total > 0) m2 /= total;
    s.sd = sqrt(fabs(m2 - mean * mean));

    /* Mode. */
    size_t mode_idx = 0;
    for (size_t i = 1; i < n; i++) {
        if (logdens[i] > logdens[mode_idx]) mode_idx = i;
    }
    s.mode = t[mode_idx];

    /* Equal-tailed intervals. */
    build_cdf(t, logdens, n, cdf);
    s.median   = inv_cdf(t, cdf, n, 0.5);
    s.et_lo_50 = inv_cdf(t, cdf, n, 0.25);
    s.et_hi_50 = inv_cdf(t, cdf, n, 0.75);
    s.et_lo_80 = inv_cdf(t, cdf, n, 0.10);
    s.et_hi_80 = inv_cdf(t, cdf, n, 0.90);
    s.et_lo_90 = inv_cdf(t, cdf, n, 0.05);
    s.et_hi_90 = inv_cdf(t, cdf, n, 0.95);
    s.et_lo_95 = inv_cdf(t, cdf, n, 0.025);
    s.et_hi_95 = inv_cdf(t, cdf, n, 0.975);

    /* HPD intervals share one ordering of the grid by descending density.
     * Selection sort: the grid is small (typically 256-1024), and its tie
     * order is what P0's intervals were computed against. */
    for (size_t i = 0; i < n; i++) idx[i] = i;
    for (size_t i = 0; i < n; i++) {
        size_t best = i;
        for (size_t j = i + 1; j < n; j++) {
            if (logdens[idx[j]] > logdens[idx[best]]) best = j;
        }
        size_t tmp = idx[i]; idx[i] = idx[best]; idx[best] = tmp;
    }
    hpd_interval(t, logdens, n, idx, 0.80, &s.hpd_lo_80, &s.hpd_hi_80);
    hpd_interval(t, logdens, n, idx, 0.90, &s.hpd_lo_90, &s.hpd_hi_90);
    hpd_interval(t, logdens, n, idx, 0.95, &s.hpd_lo_95, &s.hpd_hi_95);

    *out = s;
    rc = GD_E_OK;

done:
    if (rc != GD_E_OK && err != NULL && !gd_error_failed(err)) {
        rc = gd_error_set(err, rc, "grid: summarize failed");
    }
    gd_free(idx);
    gd_free(cdf);
    return rc;
}

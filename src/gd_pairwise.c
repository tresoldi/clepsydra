/*
 * gd_pairwise.c — pairwise dating engine for Clepsydra (M2).
 *
 * P0 classical baseline + deterministic posterior integration. The
 * hierarchical P1 model extends this with the log-rate hierarchy; the grid
 * integration machinery is shared.
 */
#include "gd.h"
#include "gd_memory.h"
#include "gd_pairwise.h"
#include "gd_grid.h"
#include "gd_math.h"
#include "gd_rng.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

double gd_p0_loglik_concept(double r, double d, int match, double y_c,
                            double eps, double q_c, gd_error *err)
{
    (void)err;
    if (r <= 0.0 || d < 0.0) {
        return -HUGE_VAL;
    }

    /* Path-length decay. */
    double s = exp(-r * d);       /* P(same inherited cognate) */

    /* Contamination mixture. */
    double p_match;
    if (eps > 0.0) {
        p_match = (1.0 - eps) * s + eps * q_c;
    } else {
        p_match = s;
    }

    /* Soft observation: L_c = y * p + (1-y) * (1-p). */
    double p_obs;
    if (match) {
        p_obs = p_match;
    } else {
        p_obs = 1.0 - p_match;
    }

    /* Fold in uncertain cognacy. */
    double L;
    if (y_c >= 1.0) {
        L = p_obs;
    } else if (y_c <= 0.0) {
        L = 1.0 - p_obs;
    } else {
        L = y_c * p_obs + (1.0 - y_c) * (1.0 - p_obs);
    }

    /* Guard against L <= 0 (log of non-positive). */
    if (L <= 0.0) {
        return -HUGE_VAL;
    }
    return log(L);
}

/* ----------------------------------------------------------------- posterior */

/*
 * We compute the posterior over T on a grid. For P0 with a fixed rate r,
 * the likelihood is a product over concepts. The marginal posterior (with a
 * uniform prior on T over [T_lo, T_hi]) is proportional to the likelihood.
 *
 * For the integration we use the trapezoidal rule on a grid in T, with
 * log-space summation via logsumexp. The posterior density is normalized so
 * that the integral over [T_lo, T_hi] is 1.
 *
 * Summaries (mean, median, quantiles, HPD) are computed from the grid CDF.
 */

/* Build the posterior grid: log-posterior at each grid point. */
static int build_grid(const double *rates, const int *matches,
                      const double *y_c, size_t n, double a1, double a2,
                      double eps, const double *q_c,
                      double T_lo, double T_hi, size_t grid_n,
                      double *grid_T, double *log_post)
{
    /* Evenly spaced grid in T. */
    double dT = (T_hi - T_lo) / (double)(grid_n - 1);
    for (size_t i = 0; i < grid_n; i++) {
        grid_T[i] = T_lo + (double)i * dT;
    }

    for (size_t i = 0; i < grid_n; i++) {
        double T = grid_T[i];
        double d = 2.0 * T - a1 - a2;
        if (d < 0.0) {
            log_post[i] = -HUGE_VAL;
            continue;
        }
        double ll = 0.0;
        for (size_t c = 0; c < n; c++) {
            double r = rates[c];
            double y = (y_c != NULL) ? y_c[c] : 1.0;
            double q = (q_c != NULL) ? q_c[c] : 0.5;
            double lc = gd_p0_loglik_concept(r, d, matches[c], y, eps, q, NULL);
            if (lc == -HUGE_VAL) { ll = -HUGE_VAL; break; }
            ll += lc;
        }
        log_post[i] = ll;
    }

    /* Normalize so the density integrates to 1 (trapezoidal), not just PMF.
     * First subtract max for stability, then compute the trapezoidal integral
     * and divide by it. */
    double m = -HUGE_VAL;
    for (size_t i = 0; i < grid_n; i++) {
        if (log_post[i] > m) m = log_post[i];
    }
    if (m == -HUGE_VAL) {
        return GD_E_RANGE;  /* every grid point has zero likelihood */
    }
    for (size_t i = 0; i < grid_n; i++) {
        log_post[i] -= m;  /* subtract max for stability */
    }
    double dTgrid = (T_hi - T_lo) / (double)(grid_n - 1);
    double integral = 0.0;
    for (size_t i = 1; i < grid_n; i++) {
        double f0 = exp(log_post[i - 1]);
        double f1 = exp(log_post[i]);
        integral += 0.5 * (f0 + f1) * dTgrid;
    }
    if (integral <= 0.0) {
        return GD_E_RANGE;
    }
    double log_integral = log(integral);
    for (size_t i = 0; i < grid_n; i++) {
        log_post[i] -= log_integral;  /* now log_post is normalized log-density */
    }
    return 0;
}

int gd_p0_posterior(const double *rates, const int *matches,
                    const double *y_c, size_t n_concepts,
                    double a1, double a2,
                    double T_lo, double T_hi,
                    double eps, const double *q_c,
                    size_t grid_n,
                    gd_pair_result *result, gd_error *err)
{
    if (n_concepts == 0 || rates == NULL || matches == NULL || result == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "p0_posterior: bad args");
        return GD_E_INVALID_ARG;
    }
    if (T_lo < a1 || T_lo < a2) {
        T_lo = (a1 > a2) ? a1 : a2;
    }
    if (T_hi <= T_lo) {
        gd_error_set(err, GD_E_RANGE, "p0_posterior: T_hi <= T_lo");
        return GD_E_RANGE;
    }
    if (grid_n < 16) grid_n = 256;

    double *grid_T = (double *)gd_malloc_array_checked(grid_n, sizeof(double), err);
    double *log_post = (double *)gd_malloc_array_checked(grid_n, sizeof(double), err);
    if (grid_T == NULL || log_post == NULL) {
        free(grid_T); free(log_post);
        gd_error_set(err, GD_E_NOMEM, "p0_posterior: alloc");
        return GD_E_NOMEM;
    }

    int rc = build_grid(rates, matches, y_c, n_concepts, a1, a2, eps, q_c,
                        T_lo, T_hi, grid_n, grid_T, log_post);
    if (rc != 0) {
        free(grid_T); free(log_post);
        gd_error_set(err, GD_E_RANGE, "p0_posterior: all-zero likelihood");
        return GD_E_RANGE;
    }

    gd_grid_summary gs;
    rc = gd_grid_summarize(grid_T, log_post, grid_n, &gs, err);
    if (rc != GD_E_OK) {
        free(grid_T); free(log_post);
        return rc;
    }

    memset(result, 0, sizeof(*result));
    result->mean   = gs.mean;
    result->sd     = gs.sd;
    result->median = gs.median;
    result->mode   = gs.mode;
    result->et_lo_50 = gs.et_lo_50; result->et_hi_50 = gs.et_hi_50;
    result->et_lo_80 = gs.et_lo_80; result->et_hi_80 = gs.et_hi_80;
    result->et_lo_90 = gs.et_lo_90; result->et_hi_90 = gs.et_hi_90;
    result->et_lo_95 = gs.et_lo_95; result->et_hi_95 = gs.et_hi_95;
    result->hpd_lo_80 = gs.hpd_lo_80; result->hpd_hi_80 = gs.hpd_hi_80;
    result->hpd_lo_90 = gs.hpd_lo_90; result->hpd_hi_90 = gs.hpd_hi_90;
    result->hpd_lo_95 = gs.hpd_lo_95; result->hpd_hi_95 = gs.hpd_hi_95;

    /* Store grid (caller may inspect). */
    result->grid_T = grid_T;
    {
        double *density = (double *)gd_malloc_array_checked(grid_n, sizeof(double), err);
        if (density != NULL) {
            for (size_t i = 0; i < grid_n; i++) {
                density[i] = exp(log_post[i]);
            }
        }
        result->grid_density = density;
    }
    result->grid_n = grid_n;

    /* Absolute anchoring: P0 never has an absolute anchor unless tip ages
     * are nonzero (ancient tips). For contemporaneous tips (a1=a2=0), report
     * lexical-time units only. */
    result->absolute_anchored = (a1 > 0.0 || a2 > 0.0);

    free(log_post);
    /* grid_T and grid_density are owned by the result. */
    return 0;
}

int gd_p0_simulate(gd_rng *rng, double r, double d, size_t n_concepts,
                   int *matches, gd_error *err)
{
    if (rng == NULL || r <= 0.0 || d < 0.0 || matches == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "p0_simulate: bad args");
        return GD_E_INVALID_ARG;
    }
    double s = exp(-r * d);  /* retention probability */
    for (size_t c = 0; c < n_concepts; c++) {
        double u = gd_rng_uniform_closed_open(rng);
        matches[c] = (u < s) ? 1 : 0;
    }
    return 0;
}

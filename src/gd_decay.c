/*
 * gd_decay.c — alternative pairwise lexical-decay laws (M12 extension).
 *
 * The three published laws ported from the distglot toolkit (Swadesh /
 * Starostin / Vasilyev), plus a configurable two-component form, each wrapped
 * in a proper Bayesian grid posterior over the divergence time so they can be
 * compared on the same data. See gd_decay.h.
 */
#include "gd_decay.h"
#include "gd_grid.h"
#include "gd_memory.h"
#include "gd.h"
#include "gd_math.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void gd_decay_params_default(gd_decay_params *p)
{
    if (p == NULL) return;
    p->swadesh_r        = 0.86;
    p->starostin_lambda = 0.05;
    p->vasilyev_floor   = 0.08;
    p->vasilyev_rate    = 0.45;
    p->tc_stable_base   = 0.20;
    p->tc_stable_decay  = 0.80;
    p->tc_stable_rate   = 0.14;
    p->tc_regular_r     = 0.86;
    p->tc_stable_weight = 0.50;
}

double gd_decay_similarity(gd_decay_model m, double T_years, const gd_decay_params *p)
{
    gd_decay_params def;
    if (p == NULL) { gd_decay_params_default(&def); p = &def; }
    if (!(T_years >= 0.0)) return NAN;
    double t = T_years / 1000.0;   /* kyr */
    switch (m) {
        case GD_DECAY_SWADESH:
            if (!(p->swadesh_r > 0.0 && p->swadesh_r < 1.0)) return NAN;
            return pow(p->swadesh_r, 2.0 * t);
        case GD_DECAY_STAROSTIN:
            if (!(p->starostin_lambda > 0.0)) return NAN;
            return exp(-2.0 * p->starostin_lambda * t * t);
        case GD_DECAY_VASILYEV: {
            /* Published single-curve fit: a decaying core above a fixed floor.
             * floor and rate are fixed constants of the law (not estimated
             * here); absolute dates therefore rest on those constants, not only
             * on the data — see docs/decay-laws.md. */
            double fl = p->vasilyev_floor, ra = p->vasilyev_rate;
            if (!(fl >= 0.0 && fl < 1.0) || !(ra > 0.0)) return NAN;
            return fl + (1.0 - fl) * exp(-ra * t) * (1.0 + ra * t);
        }
        case GD_DECAY_TWO_COMPONENT: {
            double ws = p->tc_stable_weight, wr = 1.0 - ws;
            if (!(ws >= 0.0 && ws <= 1.0) ||
                !(p->tc_regular_r > 0.0 && p->tc_regular_r < 1.0)) return NAN;
            double stable = p->tc_stable_base + p->tc_stable_decay * exp(-p->tc_stable_rate * t);
            double regular = pow(p->tc_regular_r, t);
            return ws * stable + wr * regular;
        }
    }
    return NAN;
}

double gd_decay_invert(gd_decay_model m, double similarity,
                       const gd_decay_params *p, double max_years)
{
    if (!(similarity > 0.0) || similarity > 1.0 || !(max_years > 0.0)) return NAN;
    double s0 = gd_decay_similarity(m, 0.0, p);
    if (isnan(s0)) return NAN;
    if (similarity >= s0) return 0.0;
    double smax = gd_decay_similarity(m, max_years, p);
    if (similarity < smax) return NAN;   /* below what the law reaches in range */

    double lo = 0.0, hi = max_years;     /* s(lo) >= sim >= s(hi); s decreasing */
    for (int it = 0; it < 100; it++) {
        double mid = 0.5 * (lo + hi);
        double sm = gd_decay_similarity(m, mid, p);
        if (sm > similarity) lo = mid; else hi = mid;
        if (hi - lo < 1e-6) break;
    }
    return 0.5 * (lo + hi);
}

/* Concept log-likelihood at similarity s (prob of a shared cognate). */
static double concept_loglik(double s, int match, double y_c)
{
    double p_obs = match ? s : (1.0 - s);
    double L = (y_c >= 1.0) ? p_obs
             : (y_c <= 0.0) ? (1.0 - p_obs)
             : y_c * p_obs + (1.0 - y_c) * (1.0 - p_obs);
    return (L > 0.0) ? log(L) : -HUGE_VAL;
}

int gd_decay_posterior(gd_decay_model m, const gd_decay_params *p,
                       const int *matches, const double *y_c, size_t n_concepts,
                       double T_lo, double T_hi, size_t grid_n,
                       gd_pair_result *result, gd_error *err)
{
    gd_decay_params def;
    if (p == NULL) { gd_decay_params_default(&def); p = &def; }
    if (result == NULL || matches == NULL || grid_n < 2 || !(T_hi > T_lo) || T_lo < 0.0) {
        gd_error_set(err, GD_E_INVALID_ARG, "decay_posterior: bad argument");
        return GD_E_INVALID_ARG;
    }
    double *T = (double *)gd_malloc_array_checked(grid_n, sizeof(double), err);
    double *logw = (double *)gd_malloc_array_checked(grid_n, sizeof(double), err);
    double *dens = (double *)gd_malloc_array_checked(grid_n, sizeof(double), err);
    if (!T || !logw || !dens) {
        free(T); free(logw); free(dens);
        gd_error_set(err, GD_E_NOMEM, "decay_posterior: alloc");
        return GD_E_NOMEM;
    }

    double step = (T_hi - T_lo) / (double)(grid_n - 1);
    double maxlw = -HUGE_VAL;
    for (size_t i = 0; i < grid_n; i++) {
        T[i] = T_lo + (double)i * step;
        double s = gd_decay_similarity(m, T[i], p);
        double ll = 0.0;
        if (isnan(s) || s <= 0.0 || s >= 1.0) {
            /* clamp the boundary similarity so log-likelihoods stay finite */
            s = (isnan(s)) ? 0.5 : (s <= 0.0 ? 1e-12 : 1.0 - 1e-12);
        }
        for (size_t c = 0; c < n_concepts; c++)
            ll += concept_loglik(s, matches[c], y_c ? y_c[c] : 1.0);
        logw[i] = ll;
        if (ll > maxlw) maxlw = ll;
    }

    /* Normalize (uniform prior on T) via the trapezoidal rule. */
    double Z = 0.0;
    for (size_t i = 0; i < grid_n; i++) dens[i] = exp(logw[i] - maxlw);
    for (size_t i = 1; i < grid_n; i++) Z += 0.5 * (dens[i - 1] + dens[i]) * step;
    if (!(Z > 0.0)) {
        free(T); free(logw); free(dens);
        gd_error_set(err, GD_E_RANGE, "decay_posterior: degenerate likelihood");
        return GD_E_RANGE;
    }
    for (size_t i = 0; i < grid_n; i++) dens[i] /= Z;

    gd_grid_summary gs;
    {
        /* dens is already the normalised density; gd_grid takes a log-density
         * and is invariant to its scale, so log(dens) is the right input. */
        for (size_t i = 0; i < grid_n; i++) logw[i] = log(dens[i]);
        int grc = gd_grid_summarize(T, logw, grid_n, &gs, err);
        if (grc != GD_E_OK) {
            free(T); free(logw); free(dens);
            return grc;
        }
    }

    memset(result, 0, sizeof *result);
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
    result->prior_mean = 0.5 * (T_lo + T_hi);
    result->absolute_anchored = false;   /* pairwise, no tip-age anchor */
    /* Prior dominance: posterior sd relative to the uniform prior sd. P1
     * reports this unclamped; before D2 the decay laws clamped it to [0,1],
     * so the same quantity meant two things depending on the model. */
    double prior_sd = (T_hi - T_lo) / sqrt(12.0);
    result->prior_dominance = (prior_sd > 0.0) ? (result->sd / prior_sd) : 1.0;
    result->grid_T = T;
    result->grid_density = dens;
    result->grid_n = grid_n;

    free(logw);
    (void)err;
    return 0;
}

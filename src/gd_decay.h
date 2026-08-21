/*
 * gd_decay.h — alternative pairwise lexical-decay laws (M12 extension).
 *
 * Clepsydra's pairwise likelihood assumes a constant replacement rate, so the
 * expected cognate similarity of two languages that diverged `T` years ago is
 * `exp(-r·2T)` (the Swadesh law). That is one hypothesis about how basic
 * vocabulary decays; it is not the only one in the literature. This module adds
 * three published alternatives — ported as *models to compare*, not to prefer,
 * in the spirit of docs/model-cards.md (flexibility is not correctness):
 *
 *   - SWADESH   s(T) = r^(2·T_kyr)                     constant rate (baseline)
 *   - STAROSTIN s(T) = exp(-2λ·T_kyr²)                 rate rises with time
 *               (Starostin 2000: retained words are increasingly stable, so
 *                replacement appears to accelerate)
 *   - VASILYEV  s(T) = 0.08 + 0.92·e^{-0.45 T_kyr}(1+0.45 T_kyr)
 *               (Vasilyev's coherent+independent two-component form, with an
 *                ultra-stable ~8% floor)
 *   - TWO_COMPONENT  w_s(a+b·e^{-k T_kyr}) + w_r·r^{T_kyr}
 *               (a configurable stable-core + regular-vocabulary mixture)
 *
 * These come from the `distglot` toolkit (Swadesh / Starostin / Vasilyev
 * models). Here each law defines a proper Bayesian posterior over the
 * divergence time (grid integration), so alternatives can be compared on the
 * same data by held-out prediction rather than by a point estimate.
 *
 * `s(T)` is the expected fraction of shared cognates between the two languages
 * at divergence time T (both lineages have evolved for T). Strict C99, libm.
 */
#ifndef GD_DECAY_H
#define GD_DECAY_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_pairwise.h"

#include <stddef.h>

typedef enum gd_decay_model {
    GD_DECAY_SWADESH = 0,
    GD_DECAY_STAROSTIN,
    GD_DECAY_VASILYEV,
    GD_DECAY_TWO_COMPONENT
} gd_decay_model;

typedef struct gd_decay_params {
    double swadesh_r;        /* per-kyr retention (default 0.86) */
    double starostin_lambda; /* acceleration (default 0.05) */
    /* Vasilyev single-curve law: floor + (1-floor)·exp(-rate·t)·(1+rate·t).
     * Constants are taken as published (not fit here); defaults 0.08 / 0.45. */
    double vasilyev_floor;   /* retention floor (default 0.08) */
    double vasilyev_rate;    /* decay rate per kyr (default 0.45) */
    /* Two-component (Vasilyev-style stable core + regular vocabulary). */
    double tc_stable_base;   /* a: stable floor (default 0.20) */
    double tc_stable_decay;  /* b: decaying part of the stable core (default 0.80) */
    double tc_stable_rate;   /* k: stable decay per kyr (default 0.14) */
    double tc_regular_r;     /* regular-vocabulary retention per kyr (default 0.86) */
    double tc_stable_weight; /* w_s in [0,1]; w_r = 1 - w_s (default 0.5) */
} gd_decay_params;

/* Fill sensible published defaults. */
void gd_decay_params_default(gd_decay_params *p);

/* Model names are parsed by gd_pair_model_from_string (gd_analysis.h), which
 * is the single place a user-facing model name is spelled. There is
 * deliberately no second parser here: two of them would drift into accepting
 * different spellings. */

/*
 * Expected shared-cognate similarity at divergence time `T_years` (>= 0) under
 * the given law. Monotone non-increasing in T; in [0,1]. Returns NAN on a bad
 * parameter.
 */
double gd_decay_similarity(gd_decay_model m, double T_years,
                           const gd_decay_params *p);

/*
 * Point divergence-time estimate: invert s(T) = observed similarity by
 * bisection over [0, max_years]. Returns the time in years, or NAN if the
 * similarity is below the law's floor (unreachable) or out of (0,1].
 */
double gd_decay_invert(gd_decay_model m, double similarity,
                       const gd_decay_params *p, double max_years);

/*
 * Bayesian posterior over the divergence time T for a pair, under a chosen
 * decay law, given per-concept match indicators `matches` (1 = shared cognate,
 * 0 = not) and optional per-concept cognacy confidences `y_c` (NULL = certain).
 * Uniform prior on T over [T_lo, T_hi]; grid integration with `grid_n` points.
 * Fills mean/median/sd, equal-tail intervals, mode and (optional) grid in
 * `result`. `result->grid_T/grid_density` are heap-allocated (free with
 * gd_free) when grid_n > 0.
 */
int gd_decay_posterior(gd_decay_model m, const gd_decay_params *p,
                       const int *matches, const double *y_c, size_t n_concepts,
                       double T_lo, double T_hi, size_t grid_n,
                       gd_pair_result *result, gd_error *err);

#endif /* GD_DECAY_H */

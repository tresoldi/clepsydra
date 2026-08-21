/*
 * gd_pairwise.h — internal pairwise dating engine for Clepsydra (M2/M4).
 */
#ifndef GD_PAIRWISE_H
#define GD_PAIRWISE_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"

#include <stdbool.h>
#include <stddef.h>

/*
 * The posterior over one pair's divergence time, as the pairwise engines
 * compute it. Internal: the public pairwise result is
 * gd_analysis_pair_result, which reports the same quantities in years BP
 * without the diagnostic detail the engines carry between themselves.
 *
 * `grid_T` / `grid_density` are borrowed from whoever ran the integration and
 * are NULL when no grid was requested.
 */
typedef struct gd_pair_result {
    double mean;
    double median;
    double sd;
    double mode;            /* NAN if unstable/unavailable */
    double et_lo_50,  et_hi_50;
    double et_lo_80,  et_hi_80;
    double et_lo_90,  et_hi_90;
    double et_lo_95,  et_hi_95;
    double hpd_lo_80, hpd_hi_80;
    double hpd_lo_90, hpd_hi_90;
    double hpd_lo_95, hpd_hi_95;
    double prior_mean;      /* for prior-to-posterior comparison */
    double prior_dominance; /* 0..1 diagnostic */
    bool   absolute_anchored; /* 1 if an absolute scale anchor was supplied */
    /* evaluated posterior grid (optional; NULL if not requested) */
    const double *grid_T;       /* length grid_n */
    const double *grid_density; /* posterior density at grid_T, length grid_n */
    size_t grid_n;
} gd_pair_result;

/* P0: classical baseline. */
double gd_p0_loglik_concept(double r, double d, int match, double y_c,
                            double eps, double q_c, gd_error *err);
int gd_p0_posterior(const double *rates, const int *matches,
                    const double *y_c, size_t n_concepts,
                    double a1, double a2,
                    double T_lo, double T_hi,
                    double eps, const double *q_c,
                    size_t grid_n,
                    gd_pair_result *result, gd_error *err);
int gd_p0_simulate(gd_rng *rng, double r, double d, size_t n_concepts,
                   int *matches, gd_error *err);

/* P1: hierarchical pairwise with log-rate marginalization. */

/* Per-concept reference prior parameters (from the loaded prior package). */
typedef struct {
    double log_rate_mean;   /* m_c^ref */
    double log_rate_sd;     /* s_c^ref */
} gd_concept_ref;

/*
 * P1 concept log-likelihood, marginalized over the lognormal concept rate.
 *   log r_c ~ N(log_mean, log_sd^2)
 *   p_c = P(match | T, r_c, eps) integrated over r_c.
 *
 * Uses Gauss-Legendre quadrature nodes (caller provides nodes/weights) or
 * a fixed grid in log-rate space. n_nodes typically 20-40.
 *
 * Returns the log marginal likelihood for one concept.
 */
double gd_p1_loglik_concept_marginal(double d, int match, double y_c,
                                     double log_rate_mean, double log_rate_sd,
                                     double eps, double q_c,
                                     const double *quad_nodes,
                                     const double *quad_weights,
                                     size_t n_nodes);

/*
 * P1 posterior over divergence time T, with per-concept rates marginalized
 * out via quadrature. Family clock β_f and contamination ε are integrated
 * on grids (n_beta, n_eps points). For a single-family pair, β_f is a
 * scalar offset applied to all concepts.
 *
 * sensitivity preset: 0=tight, 1=standard, 2=wide (controls hyperprior SDs).
 */
typedef enum {
    GD_SENSITIVITY_TIGHT = 0,
    GD_SENSITIVITY_STANDARD = 1,
    GD_SENSITIVITY_WIDE = 2,
} gd_sensitivity;

int gd_p1_posterior(const gd_concept_ref *refs, size_t n_concepts,
                    const int *matches, const double *y_c,
                    double a1, double a2,
                    double T_lo, double T_hi,
                    double beta_mean, double beta_sd,
                    double eps, double eps_a, double eps_b,
                    gd_sensitivity sens,
                    size_t grid_T, size_t grid_rate,
                    gd_pair_result *result, gd_error *err);

/*
 * Simulation-based calibration for pairwise P1.
 * Runs n_sbc replicates: draw T from prior, draw rates from hierarchy,
 * simulate, fit P1 posterior, record rank of true T.
 * Fills ranks[0..n_sbc-1] with the rank of the true T in each posterior.
 * Returns 0 on success.
 */
int gd_p1_sbc(gd_rng *rng, size_t n_sbc, size_t n_concepts,
              double T_lo, double T_hi,
              double log_rate_mean, double log_rate_sd,
              double beta_sd, double eps,
              gd_sensitivity sens,
              uint32_t *ranks, gd_error *err);

#endif /* GD_PAIRWISE_H */

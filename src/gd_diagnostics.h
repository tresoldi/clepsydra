/*
 * gd_diagnostics.h — MCMC convergence diagnostics for Clepsydra (M8).
 *
 * Rank-normalized split-R-hat and bulk effective sample size (ESS), following
 * Vehtari, Gelman, Simpson, Carpenter & Bürkner 2021 ("Rank-normalization,
 * folding, and localization: an improved R-hat"), the same estimators Stan
 * reports. See docs/inference.md §7.2.
 *
 * Samples are laid out chain-major: value of draw i in chain j is at
 * samples[j * n_per_chain + i]. All functions are pure and allocate only
 * bounded scratch. They never call exit()/abort(); a degenerate input (too few
 * draws, zero variance) returns a documented sentinel rather than crashing.
 */
#ifndef GD_DIAGNOSTICS_H
#define GD_DIAGNOSTICS_H

#include <stddef.h>

/*
 * Rank-normalized split-R-hat over n_chains chains of n_per_chain draws each.
 * Each chain is split in half (so 2*n_chains half-chains), the pooled draws
 * are rank-normalized to normal scores, and the classic R-hat is computed on
 * those. Converged chains give a value near 1.0; a stuck or offset chain
 * pushes it well above 1 (the > 1.01 rule of thumb).
 *
 * Returns NAN if n_chains < 2 or n_per_chain < 4, or if every draw is
 * identical (no variance to assess).
 */
double gd_rhat_rank(const double *samples, size_t n_chains, size_t n_per_chain);

/* The classic split-R-hat on raw values (no rank-normalization), exposed for
 * testing and for parameters already on a well-behaved scale. Same degenerate
 * conventions as gd_rhat_rank. */
double gd_split_rhat(const double *samples, size_t n_chains, size_t n_per_chain);

/*
 * Bulk effective sample size: the rank-normalized multi-chain ESS. Combines
 * within- and between-chain autocorrelation with Geyer's initial monotone
 * positive-sequence truncation. For well-mixed independent draws it approaches
 * the total draw count n_chains * n_per_chain; for autocorrelated draws it is
 * smaller.
 *
 * Returns NAN if n_chains < 1 or n_per_chain < 4 or the draws have no
 * variance.
 */
double gd_ess_bulk(const double *samples, size_t n_chains, size_t n_per_chain);

/* Inverse standard-normal CDF (quantile function), Acklam's rational
 * approximation refined by one Halley step; |error| < 1e-12 on (0,1).
 * Exposed because the rank normalization needs it and it is independently
 * testable. Returns -/+HUGE_VAL at 0/1 and NAN outside [0,1]. */
double gd_inv_normal_cdf(double p);

#endif /* GD_DIAGNOSTICS_H */

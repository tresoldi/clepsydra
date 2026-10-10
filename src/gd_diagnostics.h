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
 * Every estimator below follows ArviZ 0.23 exactly and is golden-tested
 * against it. Chains are split in half (first and last n/2 draws), so every
 * function needs n_per_chain >= 4; R-hat also needs n_chains >= 2. A
 * degenerate input returns NAN.
 *
 * gd_rhat_rank is max(bulk, tail): the bulk R-hat is the classic split R-hat
 * on rank-normalised draws, the tail R-hat the same on draws folded about
 * their median, |x - median|. The folded version catches chains that agree
 * in location but not in scale. Converged chains give < 1.01.
 */
double gd_rhat_rank(const double *samples, size_t n_chains, size_t n_per_chain);
double gd_rhat_bulk(const double *samples, size_t n_chains, size_t n_per_chain);
double gd_rhat_tail(const double *samples, size_t n_chains, size_t n_per_chain);

/* The classic split-R-hat on raw values (no rank-normalization). */
double gd_split_rhat(const double *samples, size_t n_chains, size_t n_per_chain);

/*
 * Effective sample sizes, with Geyer's initial positive and monotone
 * sequences (Stan/ArviZ). They may exceed the draw count for antithetic
 * chains.
 *   gd_ess_bulk  — on rank-normalised draws: reliability of central summaries.
 *   gd_ess_tail  — min of the ESS of I(x <= q05) and I(x <= q95): reliability
 *                  of the interval ends, which is what a reported HPD needs.
 *   gd_ess_mean  — on raw draws: the ESS for the posterior mean.
 */
double gd_ess_bulk(const double *samples, size_t n_chains, size_t n_per_chain);
double gd_ess_tail(const double *samples, size_t n_chains, size_t n_per_chain);
double gd_ess_mean(const double *samples, size_t n_chains, size_t n_per_chain);

/* Monte Carlo standard error of the posterior mean: sd / sqrt(ess_mean). */
double gd_mcse_mean(const double *samples, size_t n_chains, size_t n_per_chain);

/* Inverse standard-normal CDF (quantile function), Acklam's rational
 * approximation refined by one Halley step; |error| < 1e-12 on (0,1).
 * Exposed because the rank normalization needs it and it is independently
 * testable. Returns -/+HUGE_VAL at 0/1 and NAN outside [0,1]. */
double gd_inv_normal_cdf(double p);

#endif /* GD_DIAGNOSTICS_H */

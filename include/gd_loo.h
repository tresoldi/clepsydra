/*
 * gd_loo.h — leave-one-concept-out predictive accuracy and stacking (E2.8).
 *
 * PSIS-LOO (Vehtari, Gelman & Gabry 2017; Vehtari et al. 2024) estimates each
 * concept's leave-one-out predictive density from posterior draws of its
 * pointwise log-likelihood, smoothing the importance-ratio tail with a fitted
 * generalised Pareto distribution. The fitted shape k is the reliability
 * diagnostic: above 0.7 the estimate for that concept is not trustworthy.
 *
 * Stacking (Yao, Vehtari, Simpson & Gelman 2018) then chooses model weights
 * on the simplex maximising sum_i log sum_m w_m exp(elpd_loo[m][i]): the
 * weights of the mixture of posterior predictives that would have predicted
 * the held-out concepts best. Unlike Bayes factors it neither needs marginal
 * likelihoods nor goes to 0/1 when every model is wrong in a different way.
 *
 * Used by the `average` command; public so a host can compare its own runs.
 */
#ifndef GD_LOO_H
#define GD_LOO_H

#include <stddef.h>

/* Generalised Pareto fit to positive exceedances x[0..n-1], sorted ascending
 * (Zhang & Stephens 2009, with the weakly informative prior on k used by the
 * loo package). Returns k and writes the scale to *sigma. */
double gd_gpd_fit(const double *x, size_t n, double *sigma);

/* Pareto-smoothed log importance weights for one concept, in place.
 * lw[0..S-1] holds log ratios on entry (for LOO: -log p(y_i | theta_s)) and
 * normalised smoothed log weights (log-sum-exp 0) on exit. Returns k; the
 * value is +inf when the tail is too short to fit (the weights are then the
 * raw ones, normalised), 0 when all ratios are equal (the estimate is then
 * exact), and NaN on allocation failure. */
double gd_psis_smooth(double *lw, size_t S);

/* PSIS-LOO for an S x N matrix of pointwise log-likelihoods ll[s*N + i].
 * Writes elpd_i[N] and k_i[N]; returns 0, or nonzero on allocation failure. */
int gd_psis_loo(const double *ll, size_t S, size_t N, double *elpd_i, double *k_i);

/* Stacking weights for M models from their pointwise elpd_loo rows
 * elpd[m*N + i]. Writes w[M] (non-negative, summing to 1). Returns the
 * number of iterations used, or 0 on allocation failure. */
size_t gd_stacking_weights(const double *elpd, size_t M, size_t N, double *w);

#endif /* GD_LOO_H */

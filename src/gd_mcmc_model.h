/*
 * gd_mcmc_model.h — the seam between the sampler and a likelihood.
 *
 * Its own header because both sides need it and neither should include the
 * other: the sampler must not know which likelihoods exist, and a likelihood
 * must not depend on the sampler. Until E3 these types lived in gd_mcmc.h,
 * which then had to include gd_model_novel_state.h to declare the one adapter
 * — the compile-time dependency on T1 that docs/c99-alignment-plan.md
 * §Deferred names, and the reason gd_mcmc_model had exactly one
 * implementation for so long.
 */
#ifndef GD_MCMC_MODEL_H
#define GD_MCMC_MODEL_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"

#include <stddef.h>

/*
 * The sampler's current parameter state, handed to the likelihood callback.
 * `branch_logmult` (length tree->nnodes, indexed by child node; root entry 0)
 * and `concept_logdev` (length n_pointwise) are NULL when their block is off.
 */
typedef struct gd_mcmc_state {
    /*
     * The sampled rates, in log space. A likelihood that has one clock reads
     * rates[0]; one with more -- a binary CTMC has a gain and a loss rate --
     * reads as many as it declares.
     *
     * This was a single `double logr` until E3, which is why gd_mcmc_model
     * had exactly one implementation: T2 has been fully validated since M7 and
     * could not reach the sampler, because a scalar cannot carry two rates.
     * A vector rather than a second named field so that a third likelihood
     * needs no change here, and a plain array rather than an accessor so the
     * inner loop still just dereferences.
     */
    const double *rates;
    size_t        n_rates;
    const double *concept_logdev;   /* mu_c, or NULL (all zero) */
    const double *branch_logmult;   /* log m_b, or NULL (strict clock) */
    double        epsilon;          /* contamination weight in [0,1) */
    /* Spread of the integrated-out among-concept rate distribution (log sd);
     * 0 = homogeneous. A likelihood without such a mixture ignores it. */
    double        rate_mix_sd;
    int           rate_mix_family;  /* 0 log-normal, 1 gamma (CV = rate_mix_sd) */
    /* Per-node extra hazard on the edge above a tip with observation error,
     * -log(1 - e_v); NULL when no tip has one. */
    const double *tip_hazard;
} gd_mcmc_state;

/*
 * Likelihood callback. The sampler writes the current node ages into `tree`,
 * the other parameters into `st`, and expects the total data log-likelihood in
 * *out_total; if `pointwise != NULL`, fill n_pointwise per-concept log-liks.
 *
 * Optional parallel support: `clone` returns an independent context bound to
 * `tree` (a copy of the run's tree), sharing only read-only data with `ctx`;
 * `clone_free` releases one. With both set, gd_mcmc_run may run each chain in
 * its own thread on its own tree copy. NULL: the chains run one after another.
 */
typedef struct gd_mcmc_model {
    void  *ctx;
    size_t n_pointwise;
    int  (*loglik)(void *ctx, gd_tree *tree, const gd_mcmc_state *st,
                   double *pointwise, double *out_total, gd_error *err);
    void *(*clone)(void *ctx, gd_tree *tree, gd_error *err);
    void  (*clone_free)(void *clone);
    /* Optional: a hash of the data and settings the likelihood was built
     * from, so a checkpoint cannot be restored against different data.
     * NULL: the data are not part of the binding. */
    uint64_t (*fingerprint)(void *ctx);
} gd_mcmc_model;

#endif /* GD_MCMC_MODEL_H */

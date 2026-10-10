/*
 * gd_tree_prior.h — a proper prior on node ages for a fixed topology.
 *
 * Before this module the sampler had no density on uncalibrated ages: they were
 * flat on the ordered region inside a hidden [0, 1e15] box. That measure tilts
 * the root towards whatever bound exists by root^(m-1), and with the
 * contamination floor of T1 it makes the posterior improper
 * (docs/research/reviews/2026-10-10-codex.md, findings 1-2).
 *
 * The prior here is generative and proper:
 *
 *   rho                ~ the clock prior (log-normal, in the sampler)
 *   L = rho * a_root   ~ LogNormal(lex_logmean, lex_logsd)   "lexical root"
 *   free ages | a_root ~ uniform on the ordered region below the root
 *
 * L is the expected number of replacements per concept along a root-to-tip
 * path, which the cognates measure directly. The calendar root is L / rho, so
 * with no calendar calibration its posterior is the lexical posterior divided
 * by the clock prior, exactly (docs/roadmap.md E1.1). Calendar calibrations
 * multiply in as factors (E1.2), and the report shows their effective prior.
 *
 * "Uniform on the ordered region" is the MrBayes uniform-clock construction,
 * equivalently the lambda -> 0 limit of a birth-death prior conditioned on the
 * root. Its density is 1 / V(a_root), where V is the volume of the region. With
 * every tip at age 0, V = c * a_root^m. With dated tips, or internal nodes fixed
 * by calibrations, V is a polynomial in a_root above the oldest fixed age below
 * the root; it is built once, by a post-order recursion, and evaluated per move.
 */
#ifndef GD_TREE_PRIOR_H
#define GD_TREE_PRIOR_H

#include "gd_error.h"
#include "gd_tree.h"

#include <stdbool.h>
#include <stddef.h>

typedef struct gd_tree_prior {
    size_t root;
    bool   root_free;
    /* yule_kappa > 0 selects the Yule alternative (below); 0 the uniform. */
    double yule_kappa;
    const gd_tree *tree;      /* borrowed: read for the Yule ages */
    bool          *is_free;   /* owned copy, nnodes */
    double lex_logmean, lex_logsd;
    /* V(x) = sum_k coef[k] * (x / scale)^k, valid for x >= floor. */
    double *coef;
    size_t  ncoef;
    double  scale;
    double  floor;
} gd_tree_prior;

/*
 * Build the prior for `tree`. `fixed_age[v]` is the age of node v when it is
 * fixed (tips, and internal nodes pinned by a calibration), and `is_free[v]`
 * says which nodes the sampler moves. lex_logsd must be > 0.
 */
gd_tree_prior *gd_tree_prior_create(const gd_tree *tree, const bool *is_free,
                                    const double *fixed_age,
                                    double lex_logmean, double lex_logsd,
                                    gd_error *err);
void gd_tree_prior_destroy(gd_tree_prior *tp);

/*
 * The sensitivity alternative (docs/roadmap.md E1.3): a Yule (pure-birth,
 * completely sampled) prior conditioned on the root. Given the root age T,
 * each other internal age is independently Exponential(lambda) truncated to
 * (0, T) (Gernhard 2008). The rate is written as lambda = kappa / T, so the
 * prior has no calendar unit and lexical-time mode stays scale-invariant;
 * kappa -> 0 recovers the uniform prior. Defined only when every tip is
 * modern and no internal node is pinned (dated tips need the fossilized
 * birth-death prior, roadmap X3). Call after gd_tree_prior_create.
 */
int gd_tree_prior_set_yule(gd_tree_prior *tp, double kappa, gd_error *err);

/* log V(x), up to a constant; -HUGE_VAL at or below the floor. */
double gd_tree_prior_log_volume(const gd_tree_prior *tp, double root_age);

/*
 * The log prior of the ages and the clock, up to a constant that does not
 * depend on them:
 *   log N(log(rho * a_root); lex_logmean, lex_logsd) - log a_root - log V(a_root)
 * (the clock prior on log rho itself is the sampler's). With a fixed root the
 * volume term is constant and only the lexical-root factor remains.
 */
double gd_tree_prior_logpdf(const gd_tree_prior *tp, double root_age, double log_rate);

/* The full prior at the tree's current ages: gd_tree_prior_logpdf plus, for
 * the Yule prior, the terms of the non-root ages. */
double gd_tree_prior_logpdf_tree(const gd_tree_prior *tp, double log_rate);

/* The Yule kappa that matches a Yule tree's expected root age given its
 * rate: E[T_root] = (1/2 + 1/3 + ... + 1/n) / lambda for n tips. */
double gd_tree_prior_default_kappa(size_t n_tips);

#endif /* GD_TREE_PRIOR_H */

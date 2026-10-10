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

#endif /* GD_TREE_PRIOR_H */

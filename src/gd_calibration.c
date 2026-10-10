/*
 * gd_calibration.c — calibration distributions and chronology for Clepsydra (M5).
 *
 * Implements calibration type parsing, feasible-age propagation, prior
 * evaluation, and prior-only sampling. All times are years BP.
 */
#include "gd_calibration.h"
#include "gd_tree.h"
#include "gd_tsv.h"
#include "gd_memory.h"
#include "gd.h"
#include "gd_rng.h"
#include "gd_math.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * A soft uniform window (uniform_node with is_hard = 0) is flat on [p1, p2]
 * and falls off outside it in half-normal tails, each holding
 * GD_CAL_SOFT_TAIL of the mass (the MCMCTree convention of 2.5% per side).
 * With equal tail scales s on both sides the density is c on [p1, p2] and
 * c exp(-z^2/2) beyond, so c (p2 - p1) = 1 - 2 tail and c s sqrt(pi/2) = tail.
 */
#define GD_CAL_SOFT_TAIL 0.025

static double soft_uniform_scale(double p1, double p2)
{
    return GD_CAL_SOFT_TAIL * (p2 - p1) / ((1.0 - 2.0 * GD_CAL_SOFT_TAIL) * sqrt(0.5 * M_PI));
}

/*
 * Densities with unbounded support still need a finite window for the
 * sampler. It is cut where the tail beyond holds about 1e-9 of the mass
 * (6 standard deviations, on the log scale for the lognormals), which is
 * below anything an MCMC run can resolve.
 */
#define GD_CAL_SUPPORT_SD 6.0

gd_calibration_set *gd_calset_create(size_t ncals, gd_error *err)
{
    gd_calibration_set *cs = (gd_calibration_set *)gd_calloc_checked(1, sizeof(*cs), err);
    if (cs == NULL) {
        gd_error_set(err, GD_E_NOMEM, "calset: alloc");
        return NULL;
    }
    cs->ncals = ncals;
    cs->cals = (gd_calibration *)gd_calloc_checked(ncals, sizeof(gd_calibration), err);
    cs->node_lo = NULL;
    cs->node_hi = NULL;
    cs->nnodes = 0;
    cs->feasible = false;
    cs->lca_built = false;
    if (cs->cals == NULL && ncals > 0) {
        gd_error_set(err, GD_E_NOMEM, "calset: cals alloc");
        free(cs);
        return NULL;
    }
    return cs;
}

void gd_calset_destroy(gd_calibration_set *cs)
{
    if (cs == NULL) return;
    for (size_t i = 0; i < cs->ncals; i++) {
        free(cs->cals[i].calibration_id);
        free(cs->cals[i].target_id);
        free(cs->cals[i].source_key);
    }
    free(cs->cals);
    free(cs->node_lo);
    free(cs->node_hi);
    free(cs);
}

bool gd_cal_type_from_string(const char *s, size_t len, gd_cal_type *out)
{
    static const struct { const char *name; gd_cal_type type; } TYPES[] = {
        { "fixed_tip",             GD_CAL_FIXED_TIP             },
        { "uniform_node",          GD_CAL_UNIFORM_NODE          },
        { "normal_node",           GD_CAL_NORMAL_NODE           },
        { "lognormal_node",        GD_CAL_LOGNORMAL_NODE        },
        { "offset_lognormal_node", GD_CAL_OFFSET_LOGNORMAL_NODE },
        { "hard_minimum",          GD_CAL_HARD_MINIMUM          },
        { "hard_maximum",          GD_CAL_HARD_MAXIMUM          },
        { "soft_minimum",          GD_CAL_SOFT_MINIMUM          },
        { "soft_maximum",          GD_CAL_SOFT_MAXIMUM          },
    };
    if (out == NULL) return false;
    for (size_t i = 0; i < sizeof TYPES / sizeof TYPES[0]; i++) {
        if (gd_tsv_field_eq(s, len, TYPES[i].name)) {
            *out = TYPES[i].type;
            return true;
        }
    }
    return false;
}

int gd_calibrations_resolve(gd_calibration_set *cs, const gd_tree *tree,
                            gd_error *err)
{
    if (cs == NULL || tree == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "resolve: NULL");
        return GD_E_INVALID_ARG;
    }
    for (size_t i = 0; i < cs->ncals; i++) {
        size_t idx = gd_tree_find_label(tree, cs->cals[i].target_id);
        if (idx == GD_NODE_NONE) {
            gd_error_set(err, GD_E_CALIBRATION,
                         "calibration target not found in tree");
            cs->cals[i].resolved = false;
            return GD_E_CALIBRATION;
        }
        cs->cals[i].node_index = idx;
        cs->cals[i].resolved = true;
    }
    return 0;
}

/* Validate a single calibration's distribution parameters. Returns 0 if the
 * parameters are usable, GD_E_CALIBRATION otherwise. This rejects the inputs
 * that would make gd_cal_logprior divide by zero or return NaN (sd <= 0,
 * inverted uniform bounds, non-finite parameters) before they can enter the
 * feasible-bound computation or the MCMC posterior. */
static int validate_cal_params(const gd_calibration *cal, gd_error *err)
{
    double p1 = cal->p1, p2 = cal->p2, p3 = cal->p3;
    switch (cal->type) {
        case GD_CAL_FIXED_TIP:
            if (!isfinite(p1) || p1 < 0.0) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: fixed_tip age must be finite and >= 0");
                return GD_E_CALIBRATION;
            }
            break;
        case GD_CAL_UNIFORM_NODE:
            if (!isfinite(p1) || !isfinite(p2) || p1 >= p2) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: uniform_node needs finite p1 < p2");
                return GD_E_CALIBRATION;
            }
            break;
        case GD_CAL_NORMAL_NODE:
            if (!isfinite(p1) || !isfinite(p2) || p2 <= 0.0) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: normal_node needs finite mean and sd > 0");
                return GD_E_CALIBRATION;
            }
            break;
        case GD_CAL_LOGNORMAL_NODE:
            if (!isfinite(p1) || !isfinite(p2) || p2 <= 0.0) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: lognormal_node needs finite logmean and logsd > 0");
                return GD_E_CALIBRATION;
            }
            break;
        case GD_CAL_OFFSET_LOGNORMAL_NODE:
            if (!isfinite(p1) || !isfinite(p2) || !isfinite(p3) || p3 <= 0.0 || p1 < 0.0) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: offset_lognormal_node needs finite params, offset >= 0, logsd > 0");
                return GD_E_CALIBRATION;
            }
            break;
        case GD_CAL_HARD_MINIMUM:
        case GD_CAL_HARD_MAXIMUM:
            if (!isfinite(p1)) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: hard bound must be finite");
                return GD_E_CALIBRATION;
            }
            break;
        case GD_CAL_SOFT_MINIMUM:
        case GD_CAL_SOFT_MAXIMUM:
            if (!isfinite(p1) || !isfinite(p2) || p2 <= 0.0) {
                gd_error_set(err, GD_E_CALIBRATION,
                             "calibration: soft bound needs finite threshold and scale > 0");
                return GD_E_CALIBRATION;
            }
            break;
        default:
            gd_error_set(err, GD_E_CALIBRATION, "calibration: unknown type");
            return GD_E_CALIBRATION;
    }
    return 0;
}

int gd_calibrations_propagate(gd_calibration_set *cs, const gd_tree *tree,
                              gd_error *err)
{
    if (cs == NULL || tree == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "propagate: NULL");
        return GD_E_INVALID_ARG;
    }
    /* Validate every calibration's parameters before they touch the bounds or
     * the posterior. */
    for (size_t i = 0; i < cs->ncals; i++) {
        int vrc = validate_cal_params(&cs->cals[i], err);
        if (vrc != 0) return vrc;
    }
    size_t n = tree->nnodes;

    /* Allocate bounds if not already done. */
    if (cs->node_lo == NULL || cs->nnodes != n) {
        free(cs->node_lo);
        free(cs->node_hi);
        cs->node_lo = (double *)gd_malloc_array_checked(n, sizeof(double), err);
        cs->node_hi = (double *)gd_malloc_array_checked(n, sizeof(double), err);
        if (cs->node_lo == NULL || cs->node_hi == NULL) {
            gd_error_set(err, GD_E_NOMEM, "propagate: alloc");
            return GD_E_NOMEM;
        }
        cs->nnodes = n;
    }

    /* Initialize: [0, +inf) for all nodes. */
    for (size_t i = 0; i < n; i++) {
        cs->node_lo[i] = 0.0;
        cs->node_hi[i] = 1e15;  /* effectively +inf */
    }

    /* A tip no calibration mentions is an extant language: age 0. Leaving it
     * free let the sampler move "modern" languages into the past, where
     * they absorbed lexical change that belongs to the internal branches. A
     * tip with an uncertain age says so with its own calibration. */
    for (size_t i = 0; i < n; i++) {
        if (!tree->nodes[i].is_leaf) continue;
        bool targeted = false;
        for (size_t c = 0; c < cs->ncals && !targeted; c++)
            targeted = cs->cals[c].resolved && cs->cals[c].node_index == i;
        if (!targeted) cs->node_hi[i] = 0.0;
    }

    /* Apply calibrations. */
    for (size_t i = 0; i < cs->ncals; i++) {
        if (!cs->cals[i].resolved) continue;
        size_t idx = cs->cals[i].node_index;
        double p1 = cs->cals[i].p1;
        double p2 = cs->cals[i].p2;
        switch (cs->cals[i].type) {
            case GD_CAL_FIXED_TIP:
                cs->node_lo[idx] = p1;
                cs->node_hi[idx] = p1;
                break;
            case GD_CAL_UNIFORM_NODE:
                if (cs->cals[i].is_hard) {
                    if (p1 > cs->node_lo[idx]) cs->node_lo[idx] = p1;
                    if (p2 < cs->node_hi[idx]) cs->node_hi[idx] = p2;
                } else {
                    double sc = soft_uniform_scale(p1, p2);
                    double lo = p1 - GD_CAL_SUPPORT_SD * sc;
                    double hi = p2 + GD_CAL_SUPPORT_SD * sc;
                    if (lo > cs->node_lo[idx]) cs->node_lo[idx] = lo;
                    if (hi < cs->node_hi[idx]) cs->node_hi[idx] = hi;
                }
                break;
            case GD_CAL_NORMAL_NODE:
                /* The normal is in years, not log-years: its window is
                 * mean +/- k sd. (It used to be exp(mean + 5 sd), which for a
                 * 2000 +/- 200 calibration is exp(3000).) */
                {
                    double lo = p1 - GD_CAL_SUPPORT_SD * p2;
                    double hi = p1 + GD_CAL_SUPPORT_SD * p2;
                    if (lo > cs->node_lo[idx]) cs->node_lo[idx] = lo;
                    if (hi < cs->node_hi[idx]) cs->node_hi[idx] = hi;
                }
                break;
            case GD_CAL_LOGNORMAL_NODE:
            case GD_CAL_OFFSET_LOGNORMAL_NODE:
                {
                    double offset = (cs->cals[i].type == GD_CAL_OFFSET_LOGNORMAL_NODE) ? p1 : 0.0;
                    double lm = (cs->cals[i].type == GD_CAL_OFFSET_LOGNORMAL_NODE) ? p2 : p1;
                    double ls = (cs->cals[i].type == GD_CAL_OFFSET_LOGNORMAL_NODE) ? cs->cals[i].p3 : p2;
                    double lo = offset;
                    double hi = offset + exp(lm + GD_CAL_SUPPORT_SD * ls);
                    if (lo > cs->node_lo[idx]) cs->node_lo[idx] = lo;
                    if (hi < cs->node_hi[idx]) cs->node_hi[idx] = hi;
                }
                break;
            case GD_CAL_HARD_MINIMUM:
                if (p1 > cs->node_lo[idx]) cs->node_lo[idx] = p1;
                break;
            case GD_CAL_HARD_MAXIMUM:
                if (p1 < cs->node_hi[idx]) cs->node_hi[idx] = p1;
                break;
            case GD_CAL_SOFT_MINIMUM:
            case GD_CAL_SOFT_MAXIMUM:
                /* Soft bounds: don't change feasible range. */
                break;
        }
    }

    /* Propagate: parents must be older than children.
     * Upward pass: parent_lo = max(parent_lo, max(child_lo)).
     * Downward pass: child_hi = min(child_hi, parent_hi).
     * Use the postorder (children before parents) for upward, and
     * preorder (parents before children) for downward. */
    /* Upward: postorder ensures children processed before parents. */
    for (size_t i = 0; i < n; i++) {
        size_t node = tree->postorder[i];
        if (node == tree->root) continue;
        size_t parent = tree->nodes[node].parent;
        if (parent != GD_NODE_NONE) {
            if (cs->node_lo[node] > cs->node_lo[parent]) {
                cs->node_lo[parent] = cs->node_lo[node];
            }
        }
    }
    /* Downward: preorder ensures parents processed before children. */
    for (size_t i = 0; i < n; i++) {
        size_t node = tree->preorder[i];
        for (size_t c = 0; c < tree->nodes[node].nchildren; c++) {
            size_t child = tree->nodes[node].children[c];
            if (cs->node_hi[node] < cs->node_hi[child]) {
                cs->node_hi[child] = cs->node_hi[node];
            }
        }
    }

    /* Check feasibility: lo <= hi for all nodes, and parent_lo >= child_lo. */
    cs->feasible = true;
    for (size_t i = 0; i < n; i++) {
        if (cs->node_lo[i] > cs->node_hi[i]) {
            cs->feasible = false;
            gd_error_set(err, GD_E_CALIBRATION,
                         "infeasible: node bounds violated");
            return GD_E_CALIBRATION;
        }
        /* Parent must be older than child. Tested on the parent rather than
         * on "is this the root": gd_tree is public, so a caller can hand one
         * over, and this is where a second parentless node used to index at
         * GD_NODE_NONE. gd_tree_build_traversals rejects such a tree now; this
         * does not depend on that. */
        size_t parent = tree->nodes[i].parent;
        if (i != tree->root && parent != GD_NODE_NONE) {
            if (cs->node_lo[i] >= cs->node_hi[parent] && cs->node_hi[i] > cs->node_lo[parent]) {
                /* Child's minimum is at or above parent's maximum. */
                /* Allow equality for tips with fixed_tip (age = parent). */
                if (cs->node_lo[i] > cs->node_hi[parent]) {
                    cs->feasible = false;
                    gd_error_set(err, GD_E_CALIBRATION,
                                 "infeasible: child older than parent allows");
                    return GD_E_CALIBRATION;
                }
            }
        }
    }

    cs->lca_built = true;
    return 0;
}

double gd_cal_logprior(const gd_calibration *cal, double age)
{
    if (cal == NULL) return 0.0;
    double p1 = cal->p1, p2 = cal->p2, p3 = cal->p3;
    switch (cal->type) {
        case GD_CAL_FIXED_TIP:
            return (fabs(age - p1) < 1e-10) ? 0.0 : -HUGE_VAL;
        case GD_CAL_UNIFORM_NODE:
            /* Degenerate/inverted bounds carry no density. */
            if (!(p2 > p1)) return -HUGE_VAL;
            if (cal->is_hard) {
                if (age < p1 || age > p2) return -HUGE_VAL;
                return -log(p2 - p1);
            }
            {
                double base = log1p(-2.0 * GD_CAL_SOFT_TAIL) - log(p2 - p1);
                double sc = soft_uniform_scale(p1, p2);
                double z = (age < p1) ? (p1 - age) / sc
                         : (age > p2) ? (age - p2) / sc : 0.0;
                return base - 0.5 * z * z;
            }
        case GD_CAL_NORMAL_NODE:
            /* A non-positive sd would divide by zero / take log of <= 0. */
            if (!(p2 > 0.0)) return -HUGE_VAL;
            {
                double z = (age - p1) / p2;
                return -0.5 * z * z - log(p2 * sqrt(2.0 * M_PI));
            }
        case GD_CAL_LOGNORMAL_NODE:
            if (age <= 0) return -HUGE_VAL;
            if (!(p2 > 0.0)) return -HUGE_VAL;
            {
                double lr = log(age);
                double z = (lr - p1) / p2;
                return -0.5 * z * z - log(age * p2 * sqrt(2.0 * M_PI));
            }
        case GD_CAL_OFFSET_LOGNORMAL_NODE:
            if (age <= p1) return -HUGE_VAL;
            if (!(p3 > 0.0)) return -HUGE_VAL;
            {
                double lr = log(age - p1);
                double z = (lr - p2) / p3;
                return -0.5 * z * z - log((age - p1) * p3 * sqrt(2.0 * M_PI));
            }
        case GD_CAL_HARD_MINIMUM:
            return (age >= p1) ? 0.0 : -HUGE_VAL;
        case GD_CAL_HARD_MAXIMUM:
            return (age <= p1) ? 0.0 : -HUGE_VAL;
        case GD_CAL_SOFT_MINIMUM:
            {
                /* Soft penalty: if age < p1, penalize by Normal(0, p2). */
                if (age >= p1) return 0.0;
                if (!(p2 > 0.0)) return -HUGE_VAL;
                double z = (p1 - age) / p2;
                return -0.5 * z * z;
            }
        case GD_CAL_SOFT_MAXIMUM:
            {
                if (age <= p1) return 0.0;
                if (!(p2 > 0.0)) return -HUGE_VAL;
                double z = (age - p1) / p2;
                return -0.5 * z * z;
            }
    }
    return 0.0;
}

int gd_cal_sample_prior(gd_rng *rng, gd_calibration_set *cs,
                        const gd_tree *tree, double *ages, gd_error *err)
{
    if (rng == NULL || cs == NULL || tree == NULL || ages == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "sample_prior: NULL");
        return GD_E_INVALID_ARG;
    }
    if (!cs->feasible) {
        gd_error_set(err, GD_E_CALIBRATION, "sample_prior: infeasible");
        return GD_E_CALIBRATION;
    }
    /* Sample from the root down (preorder), respecting bounds.
     * For each node, sample from its calibration prior, truncated to
     * [node_lo, node_hi] and also >= parent's sampled age. */
    for (size_t i = 0; i < tree->nnodes; i++) {
        size_t node = tree->preorder[i];
        double lo = cs->node_lo[node];
        double hi = cs->node_hi[node];
        /* Must be older than any child's sampled age (we sample root-first,
         * so children aren't sampled yet). Actually for parent > child,
         * the parent must be OLDER (larger age). So we sample parent first,
         * then children must be YOUNGER (smaller age) than parent.
         * But we're going in preorder (root first), so parent is sampled
         * before children. Children's upper bound is min(child_hi, parent_age).
         * This is handled in the downward propagation. But we need to
         * enforce it during sampling too. */
        /* For the root, parent doesn't exist. For others, parent is already
         * sampled (preorder). */
        size_t parent = tree->nodes[node].parent;
        if (node != tree->root && parent != GD_NODE_NONE) {
            if (ages[parent] < hi) hi = ages[parent];
        }
        if (lo > hi) {
            /* This shouldn't happen if bounds are propagated correctly,
             * but floating-point may cause issues. */
            ages[node] = lo;
            continue;
        }

        /* Find the calibration for this node (if any). */
        const gd_calibration *cal = NULL;
        for (size_t c = 0; c < cs->ncals; c++) {
            if (cs->cals[c].resolved && cs->cals[c].node_index == node) {
                cal = &cs->cals[c];
                break;
            }
        }

        /* Sample via rejection within [lo, hi]. */
        double age;
        int attempts = 0;
        for (;;) {
            if (cal != NULL) {
                switch (cal->type) {
                    case GD_CAL_FIXED_TIP:
                        age = cal->p1;
                        break;
                    case GD_CAL_UNIFORM_NODE:
                        age = lo + (hi - lo) * gd_rng_uniform_closed_open(rng);
                        break;
                    case GD_CAL_NORMAL_NODE:
                        age = gd_rng_normal(rng, cal->p1, cal->p2);
                        break;
                    case GD_CAL_LOGNORMAL_NODE:
                        age = exp(gd_rng_normal(rng, cal->p1, cal->p2));
                        break;
                    case GD_CAL_OFFSET_LOGNORMAL_NODE:
                        age = cal->p1 + exp(gd_rng_normal(rng, cal->p2, cal->p3));
                        break;
                    case GD_CAL_HARD_MINIMUM:
                    case GD_CAL_HARD_MAXIMUM:
                    case GD_CAL_SOFT_MINIMUM:
                    case GD_CAL_SOFT_MAXIMUM:
                        /* Uniform within bounds. */
                        age = lo + (hi - lo) * gd_rng_uniform_closed_open(rng);
                        break;
                    default:
                        age = lo + (hi - lo) * gd_rng_uniform_closed_open(rng);
                        break;
                }
            } else {
                /* No calibration: uniform within bounds. */
                age = lo + (hi - lo) * gd_rng_uniform_closed_open(rng);
            }
            if (age >= lo && age <= hi) break;
            if (++attempts > 1000) {
                age = lo;  /* fallback */
                break;
            }
        }
        ages[node] = age;
    }
    return 0;
}

bool gd_cal_check_order(const gd_tree *tree, const double *ages)
{
    if (tree == NULL || ages == NULL) return false;
    for (size_t i = 0; i < tree->nnodes; i++) {
        if (i == tree->root) continue;
        size_t parent = tree->nodes[i].parent;
        if (parent != GD_NODE_NONE) {
            if (ages[parent] <= ages[i]) return false;
        }
    }
    return true;
}

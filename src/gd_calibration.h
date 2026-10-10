/*
 * gd_calibration.h — calibration distributions and chronology for Clepsydra (M5).
 *
 * All times are numeric years BP. The C engine never parses calendar strings.
 * See docs/calibrations.md for the full specification.
 */
#ifndef GD_CALIBRATION_H
#define GD_CALIBRATION_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"

#include <stdbool.h>
#include <stddef.h>

/* gd_cal_type is public (gd_types.h): a caller enumerating a calibration table
 * reads it back, so the vocabulary has to be the same on both sides. */

typedef struct gd_calibration {
    char  *calibration_id;
    char  *target_id;     /* label of the node/tip */
    gd_cal_type type;
    double p1, p2, p3;
    bool   is_hard;
    char  *source_key;
    /* Resolved node index (set by gd_calibrations_resolve). */
    size_t node_index;
    bool   resolved;
    /* Prior probability that this calibration is wrong (docs/roadmap.md
     * E2.6). With q > 0 its density is (1 - q) cal(a) + q wide(a), where
     * wide is the same calibration with its spread multiplied by
     * GD_CAL_OUTLIER_WIDEN (after OxCal's outlier model). 0 = trusted. */
    double outlier;
} gd_calibration;

#define GD_CAL_OUTLIER_WIDEN 10.0

typedef struct gd_calibration_set {
    gd_calibration *cals;
    size_t          ncals;
    /* Feasible bounds per node (after propagation). */
    double *node_lo;   /* lower bound */
    double *node_hi;   /* upper bound */
    size_t  nnodes;
    bool    feasible;  /* 1 if a feasible chronology exists */
    bool    lca_built; /* track if bounds computed */
} gd_calibration_set;

/* Lifecycle. */
gd_calibration_set *gd_calset_create(size_t ncals, gd_error *err);
void gd_calset_destroy(gd_calibration_set *cs);

/* Map a calibration-type name from a table onto the enum. A predicate: true
 * when the name is recognised and *out is set, false when it is not and *out
 * is untouched. There is no -1 sentinel to mistake for a type. */
bool gd_cal_type_from_string(const char *s, size_t len, gd_cal_type *out);

/* Resolve calibrations to node indices in the tree. Returns 0 on success;
 * GD_E_CALIBRATION if a target_id is not found in the tree. */
int gd_calibrations_resolve(gd_calibration_set *cs, const gd_tree *tree,
                            gd_error *err);

/*
 * Propagate feasible bounds:
 *  - hard bounds set node_lo / node_hi directly;
 *  - parent_lo = max(parent_lo, max(child_lo)) for each child;
 *  - parent_hi = min(parent_hi, parent's own hard max);
 *  - child_hi = min(child_hi, parent_hi);
 *  - detect contradictions (child_lo > parent_hi).
 *
 * Iterative (no deep recursion). Sets cs->feasible.
 */
int gd_calibrations_propagate(gd_calibration_set *cs, const gd_tree *tree,
                              gd_error *err);

/*
 * Evaluate the log-prior density of a node age under a calibration.
 * Returns log p(age | calibration). If no calibration targets this node,
 * returns 0.0 (flat prior). */
double gd_cal_logprior(const gd_calibration *cal, double age);

/*
 * Prior-only sampling: draw node ages from the calibration priors,
 * respecting parent > child. Fills ages[0..nnodes-1].
 * Uses rejection sampling within feasible bounds.
 */
int gd_cal_sample_prior(gd_rng *rng, gd_calibration_set *cs,
                        const gd_tree *tree, double *ages, gd_error *err);

/* Loosen a calibration's spread by factor f >= 1: the half-width of a
 * uniform window about its centre, the sd of a normal or soft bound, the
 * log-sd of a lognormal or offset lognormal. Fixed tips and hard bounds are
 * left as they are. */
void gd_cal_widen(gd_calibration *cal, double f);

/* Posterior probability, at age `age`, that an outlier-prone calibration is
 * the outlier: q wide(a) / ((1 - q) cal(a) + q wide(a)). 0 when q = 0. */
double gd_cal_outlier_posterior(const gd_calibration *cal, double age);

/* Do the ages satisfy parent > child on every edge? A predicate: true when
 * the ordering holds. */
bool gd_cal_check_order(const gd_tree *tree, const double *ages);

#endif /* GD_CALIBRATION_H */

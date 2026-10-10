/*
 * gd_sensitivity.c — scoring a chronology re-dated without one calibration.
 *
 * Moved out of cli/commands.c in D6. It was three static functions there: the
 * divergence ceiling, the shift measure, and the covered / unbounded / absent
 * decision. All three are scientific policy, none of them could be reached by
 * a test, and the browser -- which can already ask for the variant runs
 * through cal_drop -- had no way to score what came back.
 */
#include "gd_analysis.h"
#include "gd_tree.h"
#include "gd_math.h"

#include <math.h>

double gd_analysis_divergence_ceiling(const gd_analysis_tree_result *baseline,
                                      const gd_analysis_calibration *cals,
                                      size_t ncal)
{
    double scale = 0.0;
    if (baseline != NULL) {
        for (size_t i = 0; i < baseline->nages; i++) {
            if (baseline->ages[i].age.median > scale)
                scale = baseline->ages[i].age.median;
        }
    }
    if (cals != NULL) {
        for (size_t r = 0; r < ncal; r++) {
            if (cals[r].p1 > scale) scale = cals[r].p1;
            if (cals[r].p2 > scale) scale = cals[r].p2;
        }
    }
    return (scale > 0.0 ? scale : 1000.0) * 100.0;
}

/* Median age of a node, or NAN when that node's age was not sampled. */
static double node_median(const gd_analysis_tree_result *res, size_t node)
{
    const gd_analysis_node_age *a = gd_analysis_result_node(res, node);
    return (a != NULL) ? a->age.median : NAN;
}

/* How far the rest of the chronology moved. Nodes beyond the ceiling are
 * excluded: once an age is unidentified, the distance it wandered is not a
 * measurement of anything. */
static double max_other_shift(const gd_analysis_tree_result *base,
                              const gd_analysis_tree_result *var,
                              size_t skip_node, double ceiling)
{
    double worst = 0.0;
    for (size_t v = 0; v < base->n_nodes; v++) {
        if (v == skip_node) continue;
        double b = node_median(base, v), l = node_median(var, v);
        if (isnan(b) || isnan(l) || l > ceiling) continue;
        double d = fabs(l - b);
        if (d > worst) worst = d;
    }
    return worst;
}

int gd_analysis_sensitivity_score(const gd_analysis_tree_result *baseline,
                                  const gd_analysis_tree_result *variant,
                                  const gd_analysis_calibration *withheld,
                                  double ceiling,
                                  gd_sensitivity_row *out, gd_error *err)
{
    if (baseline == NULL || variant == NULL || withheld == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG,
                            "sensitivity: null argument");
    }
    if (baseline->n_nodes != variant->n_nodes) {
        return gd_error_set(err, GD_E_INVALID_ARG,
                            "sensitivity: the tree changed between runs");
    }

    gd_sensitivity_row row;
    row.outcome = GD_SENS_ABSENT;
    row.median = row.hpd_lo = row.hpd_hi = NAN;
    row.covered = false;
    row.abs_err = NAN;

    size_t t = GD_NODE_NONE;
    if (variant->tree != NULL) {
        gd_error ignore = GD_ERROR_INIT;   /* unresolvable: reported as absent */
        if (gd_tree_resolve_target(variant->tree, withheld->target, &t, &ignore) != 0)
            t = GD_NODE_NONE;
    }
    const gd_analysis_node_age *a =
        (t != GD_NODE_NONE) ? gd_analysis_result_node(variant, t) : NULL;

    row.max_other_shift = max_other_shift(baseline, variant, t, ceiling);

    if (a != NULL && gd_is_finite(a->age.median) && a->age.median <= ceiling) {
        row.outcome = GD_SENS_SCORED;
        row.median  = a->age.median;
        row.hpd_lo  = a->age.hpd_lo;
        row.hpd_hi  = a->age.hpd_hi;
        row.covered = (withheld->point >= a->age.hpd_lo &&
                       withheld->point <= a->age.hpd_hi);
        row.abs_err = fabs(a->age.median - withheld->point);
    } else if (a != NULL) {
        /* Load-bearing: nothing else constrains this node, so dropping its
         * calibration leaves the age unidentified. Said plainly rather than
         * reported as a number. */
        row.outcome = GD_SENS_UNBOUNDED;
    }

    *out = row;
    return GD_E_OK;
}

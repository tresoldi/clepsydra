/*
 * gd_grid.h — summarising a posterior held as a log-density on a 1-D grid.
 *
 * Three models produce a divergence-time posterior the same way: evaluate a
 * log-density on a grid of times, then reduce it to moments, quantiles and
 * intervals. Before D2 each of them carried its own copy of that reduction —
 * gd_pairwise.c (P0), gd_pairwise_p1.c (P1) and gd_decay.c — and the copies
 * had drifted: the inverse CDF was written twice line-for-line, and only P0
 * ever filled the HPD fields its shared result struct declares.
 *
 * This module owns the reduction. A model supplies the grid and its
 * log-density; everything below is computed here, once.
 *
 * The log-density need not be normalised: every quantity divides through by
 * the integrated mass, so an unnormalised log-posterior gives the same answer
 * as a normalised one.
 */
#ifndef GD_GRID_H
#define GD_GRID_H

#include "gd_types.h"
#include "gd_error.h"

#include <stddef.h>

/* Everything a 1-D posterior grid can say about itself.
 *
 * Equal-tailed intervals (et_*) cut equal probability from each tail; HPD
 * intervals (hpd_*) are the narrowest intervals holding the level's mass, so
 * on a skewed posterior they are narrower and not centred on the median. */
typedef struct gd_grid_summary {
    double mean, sd, median, mode;
    double et_lo_50, et_hi_50;
    double et_lo_80, et_hi_80;
    double et_lo_90, et_hi_90;
    double et_lo_95, et_hi_95;
    double hpd_lo_80, hpd_hi_80;
    double hpd_lo_90, hpd_hi_90;
    double hpd_lo_95, hpd_hi_95;
} gd_grid_summary;

/*
 * Reduce a log-density on a grid to a summary.
 *
 *   t        strictly increasing grid of n points; spacing need not be even
 *   logdens  log-density at each point, normalised or not; -HUGE_VAL is fine
 *   n        number of grid points; at least 2
 *
 * Returns GD_E_OK and fills *out, or a GD_E_* code. Allocates only scratch,
 * which it releases before returning: on failure *out is untouched and this
 * function owns nothing.
 */
int gd_grid_summarize(const double *t, const double *logdens, size_t n,
                      gd_grid_summary *out, gd_error *err);

#endif /* GD_GRID_H */

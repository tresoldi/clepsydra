/*
 * gd_math.h — internal log-space numerical helpers for Clepsydra.
 *
 * All likelihood computations use log-space to avoid underflow. The helpers
 * here are the required set from the brief §11.4 / model-spec §3 (numerical
 * stability): logsumexp, log1p, expm1, stable log(1-exp(-x)), logdiffexp,
 * explicit NaN/Inf checks, no silent underflow.
 */
#ifndef GD_MATH_H
#define GD_MATH_H

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

/* log(sum_i exp(x_i)) computed stably via the max-subtraction trick.
 * Returns -HUGE_VAL if n==0. */
double gd_logsumexp(const double *x, size_t n);

/* log(exp(a) + exp(b)) for two values, stable. */
double gd_logaddexp(double a, double b);

/* log(exp(a) - exp(b)) for a >= b. Returns -HUGE_VAL if a == b.
 * If a < b, returns NaN (caller must guard). */
double gd_logdiffexp(double a, double b);

/* log(1 - exp(x)) for x <= 0, stable. (i.e. log(1-exp(-|x|)) form.)
 * For x > 0 we use log(-expm1(x)) which requires x < 0 for a real result.
 * Convention: input is the *log* of the inner term, i.e. computes
 * log(1 - exp(logval)). For the likelihood, the common case is
 * log(1 - exp(-rate*t)) = log1mexp(-rate*t) = log1mexp(logval) with
 * logval = -rate*t < 0.
 *   For logval near 0^- (i.e. exp(logval) near 1): use log(-expm1(logval)).
 *   For logval << 0  (i.e. exp(logval) near 0): use log1p(-exp(logval)).
 * If logval >= 0, returns NaN (1-exp(>=0) <= 0 has no real log). */
double gd_log1mexp(double logval);

/* Stable log(1+x) — wraps log1p, handling NaN/Inf. */
double gd_safe_log1p(double x);

/* Stable expm1. */
double gd_safe_expm1(double x);

/* Is x finite and not NaN? (isfinite wrapper, C99). A predicate. */
bool   gd_is_finite(double x);

/* Clamp x to [lo, hi]. */
double gd_clamp(double x, double lo, double hi);

/* Sign function: -1, 0, or 1. */
int    gd_sign(double x);

#endif /* GD_MATH_H */

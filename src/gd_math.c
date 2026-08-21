/*
 * gd_math.c — log-space numerical helpers for Clepsydra.
 */
#include "gd_math.h"

#include <float.h>
#include <math.h>

double gd_logsumexp(const double *x, size_t n)
{
    if (n == 0 || x == NULL) {
        return -HUGE_VAL;
    }
    /* Find the max. */
    double m = x[0];
    for (size_t i = 1; i < n; i++) {
        if (x[i] > m) {
            m = x[i];
        }
    }
    /* If max is -inf, the sum is -inf (avoid -inf - (-inf) = NaN). */
    if (m == -HUGE_VAL) {
        return -HUGE_VAL;
    }
    /* If max is +inf, the sum is +inf. */
    if (m == HUGE_VAL) {
        return HUGE_VAL;
    }
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += exp(x[i] - m);
    }
    return m + log(sum);
}

double gd_logaddexp(double a, double b)
{
    if (a == -HUGE_VAL) {
        return b;
    }
    if (b == -HUGE_VAL) {
        return a;
    }
    double m = (a > b) ? a : b;
    return m + log(exp(a - m) + exp(b - m));
}

double gd_logdiffexp(double a, double b)
{
    /* log(exp(a) - exp(b)); requires a > b for a positive result. */
    if (a == b) {
        return -HUGE_VAL;
    }
    if (a < b) {
        return NAN;  /* caller must guard */
    }
    if (a == -HUGE_VAL) {
        return -HUGE_VAL;  /* both -inf */
    }
    /* exp(a) - exp(b) = exp(a) * (1 - exp(b-a)); a >= b so b-a <= 0. */
    return a + gd_log1mexp(b - a);
}

double gd_log1mexp(double logval)
{
    /* log(1 - exp(logval)); requires logval < 0 for a real, positive result.
     * If logval >= 0, exp(logval) >= 1, so 1-exp <= 0: return NaN. */
    if (isnan(logval)) {
        return NAN;
    }
    if (logval >= 0.0) {
        /* exp(logval) >= 1; 1 - exp <= 0 */
        return (logval == 0.0) ? -HUGE_VAL : NAN;
    }
    /* logval < 0. Two regimes (Mächler 2012):
     *   if logval < -log(2):  1 - exp(logval) ~ 1  -> use log1p(-exp(logval))
     *   if logval >= -log(2): exp(logval) ~ 1      -> use log(-expm1(logval))
     */
    const double log_half = -0.69314718055994530942;  /* log(0.5) */
    if (logval < log_half) {
        return log1p(-exp(logval));
    }
    return log(-expm1(logval));
}

double gd_safe_log1p(double x)
{
    if (isnan(x)) {
        return NAN;
    }
    if (x == -1.0) {
        return -HUGE_VAL;
    }
    return log1p(x);
}

double gd_safe_expm1(double x)
{
    if (isnan(x)) {
        return NAN;
    }
    return expm1(x);
}

bool gd_is_finite(double x)
{
    return isfinite(x) != 0;
}

double gd_clamp(double x, double lo, double hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

int gd_sign(double x)
{
    if (x > 0.0) {
        return 1;
    }
    if (x < 0.0) {
        return -1;
    }
    return 0;
}

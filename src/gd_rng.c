/*
 * gd_rng.c — PCG32 random number generator for Clepsydra.
 *
 * PCG32 reference (Melissa O'Neill, 2014: http://www.pcg-random.org).
 * 64-bit state, 64-bit increment (must be odd), 32-bit output via XSH-RR.
 * Period 2^64, streams 2^63.
 *
 * The public API in gd.h wraps these _impl functions.
 */
#include "gd.h"
#include "gd_rng.h"
#include "gd_math.h"
#include "gd_memory.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* PCG32 step. */
static uint32_t pcg32_next(gd_rng *r)
{
    uint64_t oldstate = r->state;
    /* LCG advance: state = oldstate * 6364136223846793005 + inc */
    r->state = oldstate * 6364136223846793005ULL + r->inc;
    /* XSH-RR output function. */
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-(int32_t)rot) & 31));
}

gd_rng *gd_rng_create(uint64_t seed, uint64_t stream, gd_error *err)
{
    gd_rng *r = (gd_rng *)gd_malloc_checked(sizeof(gd_rng), err);
    if (r == NULL) {
        return NULL;
    }
    /* inc must be odd; derive from stream. */
    r->inc = (stream << 1u) | 1u;   /* odd */
    r->state = 0U;
    /* The standard pcg32_srandom_r: state = 0; then step with inc; then
     * add seed; then step again. */
    pcg32_next(r);
    r->state += seed;
    pcg32_next(r);
    return r;
}

void gd_rng_destroy(gd_rng *rng)
{
    gd_free(rng);
}

uint32_t gd_rng_next32(gd_rng *rng)
{
    return pcg32_next(rng);
}

uint32_t gd_rng_bounded(gd_rng *rng, uint32_t bound)
{
    /* Lemire's method: no modulo bias, no division by zero.
     * Compute 2^64-range product, take high 32 bits, reject if below bound. */
    uint32_t range = bound;
    /* An empty range has no value to return; 0 is the documented answer, and
     * no draw is consumed. */
    if (range == 0) return 0;
    /* Lemire: threshold = 2^32 mod range, computed in unsigned arithmetic.
     * It used to negate through int32_t, which is undefined behaviour for
     * range = 0x80000000. */
    uint32_t threshold = (0u - range) % range;
    for (;;) {
        uint32_t x = pcg32_next(rng);
        /* Use 64-bit multiply to get high and low halves. */
        uint64_t m = (uint64_t)x * (uint64_t)range;
        uint32_t l = (uint32_t)m;
        if (l >= threshold) {
            return (uint32_t)(m >> 32);
        }
    }
}

double gd_rng_uniform_closed_open(gd_rng *rng)
{
    /* [0, 1): 24 or 53 bits of randomness. Use 53 bits for double precision:
     * combine two 32-bit draws, take top 53 bits, divide by 2^53.
     * This gives a uniform double in [0, 1). */
    uint32_t hi = pcg32_next(rng);
    uint32_t lo = pcg32_next(rng);
    uint64_t bits = ((uint64_t)hi << 21u) | ((uint64_t)lo >> 11u);
    return (double)bits * (1.0 / 9007199254740992.0);  /* 2^53 */
}

double gd_rng_uniform_open(gd_rng *rng)
{
    /* (0, 1): draw in [0,1) and retry if 0 (extremely rare). */
    double u;
    do {
        u = gd_rng_uniform_closed_open(rng);
    } while (u == 0.0);
    return u;
}

/* Ziggurat-free normal via Box-Muller (simpler, sufficient for our needs;
 * we are single-threaded and correctness matters more than peak speed). */
double gd_rng_normal(gd_rng *rng, double mean, double sd)
{
    if (sd <= 0.0) {
        return mean;
    }
    double u1 = gd_rng_uniform_open(rng);
    double u2 = gd_rng_uniform_closed_open(rng);
    double r = sqrt(-2.0 * log(u1));
    double z = r * cos(2.0 * 3.14159265358979323846 * u2);
    return mean + sd * z;
}

double gd_rng_exponential(gd_rng *rng, double rate)
{
    if (rate <= 0.0) {
        return HUGE_VAL;
    }
    double u = gd_rng_uniform_open(rng);
    return -log(u) / rate;
}

/* Marsaglia-Tsang gamma sampler (shape >= 1). For shape < 1, use the
 * boosting: gamma(shape) = gamma(shape+1) * U^(1/shape). */
double gd_rng_gamma(gd_rng *rng, double shape, double scale)
{
    if (shape <= 0.0 || scale <= 0.0) {
        return 0.0;
    }
    if (shape < 1.0) {
        /* Boost: g = gamma(shape+1); return g * U^(1/shape). */
        double g = gd_rng_gamma(rng, shape + 1.0, 1.0);
        double u = gd_rng_uniform_open(rng);
        return scale * g * pow(u, 1.0 / shape);
    }
    double d = shape - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do {
            x = gd_rng_normal(rng, 0.0, 1.0);
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v;
        double u = gd_rng_uniform_closed_open(rng);
        double x2 = x * x;
        if (u < 1.0 - 0.0331 * x2 * x2) {
            return scale * d * v;
        }
        if (log(u) < 0.5 * x2 + d * (1.0 - v + log(v))) {
            return scale * d * v;
        }
        /* else: reject and retry. */
    }
}

/* One Gamma(shape, 1) draw, returned both as a value and as its log. The
 * value is computed exactly as gd_rng_gamma computes it (same draws, same
 * arithmetic), so callers that only need the value are bit-identical; the
 * log survives the underflow of tiny shapes, where the value is 0. */
static double gamma_draw_log(gd_rng *rng, double shape, double *value)
{
    if (shape < 1.0) {
        double g = gd_rng_gamma(rng, shape + 1.0, 1.0);
        double u = gd_rng_uniform_open(rng);
        *value = 1.0 * g * pow(u, 1.0 / shape);
        return log(g) + log(u) / shape;
    }
    *value = gd_rng_gamma(rng, shape, 1.0);
    return log(*value);
}

double gd_rng_beta(gd_rng *rng, double a, double b)
{
    if (a <= 0.0 || b <= 0.0) {
        return NAN;
    }
    /* Beta(a,b) = X/(X+Y), X~Gamma(a,1), Y~Gamma(b,1). */
    double x, y;
    double lx = gamma_draw_log(rng, a, &x);
    double ly = gamma_draw_log(rng, b, &y);
    double s = x + y;
    if (s > 0.0 && isfinite(s)) {
        return x / s;
    }
    /* Both draws underflowed (very small shapes). The ratio is still defined
     * on the log scale: X/(X+Y) = 1/(1 + exp(lY - lX)). This used to return
     * 0.5, a point mass where a U-shaped Beta has almost none. */
    return 1.0 / (1.0 + exp(ly - lx));
}

size_t gd_rng_categorical(gd_rng *rng, const double *probs, size_t n)
{
    if (n == 0 || probs == NULL) {
        return 0;
    }
    /* Compute cumulative sum, then draw uniform and binary search.
     * We do a linear scan here for simplicity (n is small in practice). */
    double u = gd_rng_uniform_closed_open(rng);
    double cum = 0.0;
    for (size_t i = 0; i < n; i++) {
        cum += probs[i];
        if (u < cum) {
            return i;
        }
    }
    return n - 1;  /* floating-point roundup fallback */
}

int gd_rng_serialize(const gd_rng *rng, unsigned char *buf, size_t len)
{
    if (rng == NULL || buf == NULL) {
        return GD_E_INVALID_ARG;
    }
    if (len < GD_RNG_SERIAL_SIZE) {
        return GD_E_RANGE;
    }
    /* Little-endian serialization of state (8 bytes) + inc (8 bytes). */
    for (int i = 0; i < 8; i++) {
        buf[i]     = (unsigned char)(rng->state >> (8 * i));
        buf[i + 8] = (unsigned char)(rng->inc >> (8 * i));
    }
    return GD_E_OK;
}

int gd_rng_restore(gd_rng *rng, const unsigned char *buf, size_t len)
{
    if (rng == NULL || buf == NULL) {
        return GD_E_INVALID_ARG;
    }
    if (len < GD_RNG_SERIAL_SIZE) {
        return GD_E_RANGE;
    }
    uint64_t state = 0, inc = 0;
    for (int i = 0; i < 8; i++) {
        state |= ((uint64_t)buf[i]) << (8 * i);
        inc    |= ((uint64_t)buf[i + 8]) << (8 * i);
    }
    rng->state = state;
    rng->inc = inc;
    return GD_E_OK;
}

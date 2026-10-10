/*
 * gd_rng.h — internal RNG API (PCG32) for Clepsydra.
 *
 * PCG32: Permuted Conggruential Generator, 32-bit output from a 64-bit
 * state and 64-bit stream/increment. Reference implementation by Melissa
 * O'Neill (PCG, 2014). We use the standard pcg32 (LCG + XSH-RR output).
 *
 * Properties required by the brief §11.5:
 *  - independent streams (the stream/increment parameter);
 *  - published test vectors (tests/unit/test_rng.c);
 *  - deterministic serialization (serialize/restore);
 *  - uniform integers without modulo bias (rejection);
 *  - open-interval and closed-open uniform doubles;
 *  - normal, exponential, gamma, beta, categorical sampling.
 */
#ifndef GD_RNG_H
#define GD_RNG_H

#include "gd_types.h"
#include "gd_error.h"

#include <stddef.h>
#include <stdint.h>

struct gd_rng {
    uint64_t state;
    uint64_t inc;    /* stream; must be odd */
};

/* Create with a seed and a stream id. stream is used to derive the odd
 * increment. Two chains use the same seed with different stream ids, or
 * different seeds. */
gd_rng *gd_rng_create(uint64_t seed, uint64_t stream, gd_error *err);
void    gd_rng_destroy(gd_rng *rng);

/* Core 32-bit output. */
uint32_t gd_rng_next32(gd_rng *rng);

/* Unbiased bounded integer in [0, bound) via Lemire's method (no modulo
 * bias). bound = 0 is an empty range: it returns 0 and consumes no draw. */
uint32_t gd_rng_bounded(gd_rng *rng, uint32_t bound);

/* Uniform doubles. */
double gd_rng_uniform_closed_open(gd_rng *rng);  /* [0, 1) */
double gd_rng_uniform_open(gd_rng *rng);         /* (0, 1) */

/* Distributions. */
double gd_rng_normal(gd_rng *rng, double mean, double sd);
double gd_rng_exponential(gd_rng *rng, double rate);
double gd_rng_gamma(gd_rng *rng, double shape, double scale);
double gd_rng_beta(gd_rng *rng, double a, double b);
size_t gd_rng_categorical(gd_rng *rng, const double *probs, size_t n);

/* Serialized size of a generator: state (8 bytes) + inc (8), little-endian. */
#define GD_RNG_SERIAL_SIZE 16

/* Write / read the generator state. `buf` must have room for at least
 * GD_RNG_SERIAL_SIZE bytes; a shorter one is GD_E_RANGE. Returns GD_E_OK on
 * success. */
int gd_rng_serialize(const gd_rng *rng, unsigned char *buf, size_t len);
int gd_rng_restore(gd_rng *rng, const unsigned char *buf, size_t len);

#endif /* GD_RNG_H */

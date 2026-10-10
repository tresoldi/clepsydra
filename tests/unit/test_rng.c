/*
 * test_rng.c — tests for PCG32 RNG and distributions.
 *
 * PCG32 has published test vectors from Melissa O'Neill's reference
 * implementation. We verify the exact output sequence for a known
 * seed+stream, plus distribution moments and serialization round-trip.
 */
#include "test_runner.h"
#include "gd_rng.h"
#include "gd.h"

#include <stdint.h>

TEST(test_rng_next32_known_vector)
{
    /* PCG32 with seed=42, stream=54 (inc = 54<<1 | 1 = 109).
     * Reference output from O'Neill's pcg32 with the same seeding:
     *   pcg32_srandom_r(rng, 42, 54) then successive pcg32_next_r.
     * The first few 32-bit outputs are:
     *   0xa521ca5e, 0x86c26885, 0x5f28b921, 0x18f6f4eb, 0xe0285d9a, ...
     * (These match the canonical pcg32 demo sequence.)
     */
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(42u, 54u, &err);
    CHECK(r != NULL);

    /* The canonical pcg32 demo (seed=42, stream=54) outputs, verified
     * against an independent PCG32 implementation using the same seeding
     * (state=0; step; state+=seed; step; then first output): */
    uint32_t expected[] = {
        0xa15c02b7u, 0x7b47f409u, 0xba1d3330u, 0x83d2f293u,
        0xbfa4784bu, 0xcbed606eu, 0xbfc6a3adu, 0x812fff6du
    };
    for (int i = 0; i < 8; i++) {
        uint32_t v = gd_rng_next32(r);
        if (v != expected[i]) {
            printf("  FAIL: rng output[%d] = 0x%08x, expected 0x%08x\n",
                   i, v, expected[i]);
            g_test_failed++;
            gd_rng_destroy(r);
            return;
        }
    }
    gd_rng_destroy(r);
}

TEST(test_rng_bounded_unbiased)
{
    /* Bounded sampling: every value in [0, bound) should appear. */
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(123, 1, &err);
    CHECK(r != NULL);
    uint32_t bound = 7;
    int counts[7] = {0};
    for (int i = 0; i < 100000; i++) {
        uint32_t v = gd_rng_bounded(r, bound);
        CHECK(v < bound);
        counts[v]++;
    }
    /* Every value should appear at least once in 100k draws with bound=7. */
    for (int i = 0; i < 7; i++) {
        CHECK(counts[i] > 0);
    }
    /* Rough uniformity: each count within 10% of expected (14285). */
    for (int i = 0; i < 7; i++) {
        CHECK(counts[i] > 12000 && counts[i] < 16500);
    }
    gd_rng_destroy(r);
}

TEST(test_rng_uniform_range)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(999, 2, &err);
    CHECK(r != NULL);
    for (int i = 0; i < 100000; i++) {
        double u = gd_rng_uniform_closed_open(r);
        CHECK(u >= 0.0 && u < 1.0);
        double o = gd_rng_uniform_open(r);
        CHECK(o > 0.0 && o < 1.0);
    }
    gd_rng_destroy(r);
}

TEST(test_rng_uniform_moments)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(7, 3, &err);
    CHECK(r != NULL);
    double sum = 0.0, sum2 = 0.0;
    int N = 1000000;
    for (int i = 0; i < N; i++) {
        double u = gd_rng_uniform_closed_open(r);
        sum  += u;
        sum2 += u * u;
    }
    double mean = sum / N;
    double var  = sum2 / N - mean * mean;
    /* E[U[0,1)] = 0.5, Var = 1/12 ≈ 0.08333 */
    CHECK_APPROX(mean, 0.5, 0.005);
    CHECK_APPROX(var, 1.0 / 12.0, 0.003);
    gd_rng_destroy(r);
}

TEST(test_rng_normal_moments)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(11, 4, &err);
    CHECK(r != NULL);
    double sum = 0.0, sum2 = 0.0;
    int N = 200000;
    double mean0 = 5.0, sd = 2.0;
    for (int i = 0; i < N; i++) {
        double z = gd_rng_normal(r, mean0, sd);
        sum  += z;
        sum2 += (z - mean0) * (z - mean0);
    }
    double mean = sum / N;
    double var  = sum2 / N;
    CHECK_APPROX(mean, mean0, 0.05);
    CHECK_APPROX(sqrt(var), sd, 0.05);
    gd_rng_destroy(r);
}

TEST(test_rng_exponential_moments)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(13, 5, &err);
    CHECK(r != NULL);
    double sum = 0.0;
    int N = 200000;
    double rate = 2.5;
    for (int i = 0; i < N; i++) {
        sum += gd_rng_exponential(r, rate);
    }
    double mean = sum / N;
    /* E[Exp(rate)] = 1/rate = 0.4 */
    CHECK_APPROX(mean, 1.0 / rate, 0.01);
    gd_rng_destroy(r);
}

TEST(test_rng_gamma_moments)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(17, 6, &err);
    CHECK(r != NULL);
    double sum = 0.0;
    int N = 200000;
    double shape = 3.0, scale = 2.0;
    for (int i = 0; i < N; i++) {
        sum += gd_rng_gamma(r, shape, scale);
    }
    double mean = sum / N;
    /* E[Gamma(k,θ)] = k*θ = 6.0 */
    CHECK_APPROX(mean, shape * scale, 0.1);
    gd_rng_destroy(r);
}

TEST(test_rng_beta_moments)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(19, 7, &err);
    CHECK(r != NULL);
    double sum = 0.0;
    int N = 200000;
    double a = 2.0, b = 5.0;
    for (int i = 0; i < N; i++) {
        double v = gd_rng_beta(r, a, b);
        CHECK(v >= 0.0 && v <= 1.0);
        sum += v;
    }
    double mean = sum / N;
    /* E[Beta(2,5)] = 2/(2+5) = 0.2857 */
    CHECK_APPROX(mean, a / (a + b), 0.005);
    gd_rng_destroy(r);
}

TEST(test_rng_categorical)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r = gd_rng_create(23, 8, &err);
    CHECK(r != NULL);
    double probs[] = {0.1, 0.2, 0.7};
    int counts[3] = {0, 0, 0};
    int N = 100000;
    for (int i = 0; i < N; i++) {
        size_t idx = gd_rng_categorical(r, probs, 3);
        CHECK(idx < 3);
        counts[idx]++;
    }
    /* Expected: 10000, 20000, 70000; allow 5% tolerance. */
    CHECK(counts[0] > 9000 && counts[0] < 11000);
    CHECK(counts[1] > 18000 && counts[1] < 22000);
    CHECK(counts[2] > 66000 && counts[2] < 74000);
    gd_rng_destroy(r);
}

TEST(test_rng_serialize_roundtrip)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r1 = gd_rng_create(42, 54, &err);
    CHECK(r1 != NULL);
    /* Advance a few steps. */
    for (int i = 0; i < 5; i++) gd_rng_next32(r1);

    unsigned char buf[16];
    CHECK(gd_rng_serialize(r1, buf, GD_RNG_SERIAL_SIZE) == 0);

    gd_rng *r2 = gd_rng_create(1, 1, &err);
    CHECK(r2 != NULL);
    CHECK(gd_rng_restore(r2, buf, GD_RNG_SERIAL_SIZE) == 0);

    /* Both should now produce identical sequences. */
    for (int i = 0; i < 100; i++) {
        uint32_t a = gd_rng_next32(r1);
        uint32_t b = gd_rng_next32(r2);
        CHECK(a == b);
    }
    gd_rng_destroy(r1);
    gd_rng_destroy(r2);
}

TEST(test_rng_deterministic_same_seed)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r1 = gd_rng_create(42, 54, &err);
    gd_rng *r2 = gd_rng_create(42, 54, &err);
    CHECK(r1 != NULL);
    CHECK(r2 != NULL);
    for (int i = 0; i < 100; i++) {
        CHECK(gd_rng_next32(r1) == gd_rng_next32(r2));
    }
    gd_rng_destroy(r1);
    gd_rng_destroy(r2);
}

TEST(test_rng_different_streams_diverge)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *r1 = gd_rng_create(42, 1, &err);
    gd_rng *r2 = gd_rng_create(42, 2, &err);
    CHECK(r1 != NULL);
    CHECK(r2 != NULL);
    int differ = 0;
    for (int i = 0; i < 100; i++) {
        if (gd_rng_next32(r1) != gd_rng_next32(r2)) differ++;
    }
    CHECK(differ > 90);  /* almost all should differ */
    gd_rng_destroy(r1);
    gd_rng_destroy(r2);
}

/* E0.14 (codex C28): the bounded draw is defined at its edges, and a Beta
 * with tiny shapes is U-shaped instead of collapsing to a point at 0.5. */
TEST(test_rng_edge_cases)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(11u, 3u, &err);
    CHECK(rng != NULL);
    CHECK(gd_rng_bounded(rng, 0u) == 0u);
    CHECK(gd_rng_bounded(rng, 1u) == 0u);
    for (int i = 0; i < 1000; i++) {
        uint32_t v = gd_rng_bounded(rng, 0x80000000u);
        CHECK(v < 0x80000000u);
    }
    size_t n = 20000, near_half = 0, low = 0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double x = gd_rng_beta(rng, 1e-3, 1e-3);
        CHECK(x >= 0.0 && x <= 1.0);
        sum += x;
        if (fabs(x - 0.5) < 0.25) near_half++;
        if (x < 0.5) low++;
    }
    gd_rng_destroy(rng);
    /* Beta(0.001, 0.001) puts about 0.1% of its mass in (0.25, 0.75). */
    CHECK(near_half < n / 100);
    CHECK_APPROX(sum / (double)n, 0.5, 0.02);
    CHECK(low > n / 3 && low < 2 * n / 3);
}

static void run_tests(void)
{
    RUN(test_rng_edge_cases);
    RUN(test_rng_next32_known_vector);
    RUN(test_rng_bounded_unbiased);
    RUN(test_rng_uniform_range);
    RUN(test_rng_uniform_moments);
    RUN(test_rng_normal_moments);
    RUN(test_rng_exponential_moments);
    RUN(test_rng_gamma_moments);
    RUN(test_rng_beta_moments);
    RUN(test_rng_categorical);
    RUN(test_rng_serialize_roundtrip);
    RUN(test_rng_deterministic_same_seed);
    RUN(test_rng_different_streams_diverge);
}

TEST_MAIN()

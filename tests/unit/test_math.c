/*
 * test_math.c — tests for log-space numerical helpers.
 */
#include "test_runner.h"
#include "gd_math.h"

TEST(test_logsumexp_basic)
{
    double x[] = {0.0, 0.0, 0.0};
    /* log(exp(0)+exp(0)+exp(0)) = log(3) */
    CHECK_APPROX(gd_logsumexp(x, 3), log(3.0), 1e-12);
}

TEST(test_logsumexp_single)
{
    double x[] = {5.5};
    CHECK_APPROX(gd_logsumexp(x, 1), 5.5, 1e-12);
}

TEST(test_logsumexp_empty)
{
    CHECK(gd_logsumexp(NULL, 0) == -HUGE_VAL);
    double x[] = {1.0};
    CHECK(gd_logsumexp(x, 0) == -HUGE_VAL);
}

TEST(test_logsumexp_neginf)
{
    double x[] = {-HUGE_VAL, -HUGE_VAL};
    CHECK(gd_logsumexp(x, 2) == -HUGE_VAL);
}

TEST(test_logsumexp_large)
{
    double x[] = {1000.0, 1000.0};
    /* Should not overflow: result = 1000 + log(2) */
    CHECK_APPROX(gd_logsumexp(x, 2), 1000.0 + log(2.0), 1e-9);
}

TEST(test_logaddexp)
{
    CHECK_APPROX(gd_logaddexp(0.0, 0.0), log(2.0), 1e-12);
    CHECK_APPROX(gd_logaddexp(-HUGE_VAL, 3.0), 3.0, 1e-12);
    CHECK_APPROX(gd_logaddexp(3.0, -HUGE_VAL), 3.0, 1e-12);
}

TEST(test_logdiffexp)
{
    /* log(exp(3) - exp(2)) = log(20.0855 - 7.3890) = log(12.6965) */
    double expected = log(exp(3.0) - exp(2.0));
    CHECK_APPROX(gd_logdiffexp(3.0, 2.0), expected, 1e-9);
    /* Equal values -> -inf */
    CHECK(gd_logdiffexp(5.0, 5.0) == -HUGE_VAL);
    /* a < b -> NaN */
    CHECK(isnan(gd_logdiffexp(2.0, 3.0)));
}

TEST(test_log1mexp_near_zero)
{
    /* log(1 - exp(-0.001)) = log(1 - 0.99900) ≈ log(0.0009995) ≈ -6.908 */
    /* For logval = -0.001 (near 0): use log(-expm1(-0.001)) */
    double x = -0.001;
    double expected = log(1.0 - exp(x));
    CHECK_APPROX(gd_log1mexp(x), expected, 1e-12);
}

TEST(test_log1mexp_far)
{
    /* log(1 - exp(-5)) = log(1 - 0.006738) ≈ log(0.993262) ≈ -0.00676 */
    double x = -5.0;
    double expected = log(1.0 - exp(x));
    CHECK_APPROX(gd_log1mexp(x), expected, 1e-12);
}

TEST(test_log1mexp_nan)
{
    /* logval >= 0 should give NaN (or -inf for exactly 0) */
    CHECK(isnan(gd_log1mexp(1.0)));
    CHECK(gd_log1mexp(0.0) == -HUGE_VAL);
    CHECK(isnan(gd_log1mexp(NAN)));
}

TEST(test_safe_log1p)
{
    CHECK_APPROX(gd_safe_log1p(0.0), 0.0, 1e-12);
    CHECK_APPROX(gd_safe_log1p(1.0), log(2.0), 1e-12);
    CHECK(gd_safe_log1p(-1.0) == -HUGE_VAL);
    CHECK(isnan(gd_safe_log1p(NAN)));
}

TEST(test_clamp_sign)
{
    CHECK_APPROX(gd_clamp(5.0, 0.0, 10.0), 5.0, 0);
    CHECK_APPROX(gd_clamp(-5.0, 0.0, 10.0), 0.0, 0);
    CHECK_APPROX(gd_clamp(15.0, 0.0, 10.0), 10.0, 0);
    CHECK(gd_sign(3.0) == 1);
    CHECK(gd_sign(-2.0) == -1);
    CHECK(gd_sign(0.0) == 0);
}

TEST(test_is_finite)
{
    CHECK(gd_is_finite(1.0) == 1);
    CHECK(gd_is_finite(HUGE_VAL) == 0);
    CHECK(gd_is_finite(NAN) == 0);
}

static void run_tests(void)
{
    RUN(test_logsumexp_basic);
    RUN(test_logsumexp_single);
    RUN(test_logsumexp_empty);
    RUN(test_logsumexp_neginf);
    RUN(test_logsumexp_large);
    RUN(test_logaddexp);
    RUN(test_logdiffexp);
    RUN(test_log1mexp_near_zero);
    RUN(test_log1mexp_far);
    RUN(test_log1mexp_nan);
    RUN(test_safe_log1p);
    RUN(test_clamp_sign);
    RUN(test_is_finite);
}

TEST_MAIN()

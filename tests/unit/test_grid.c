/*
 * test_grid.c — tests for the posterior-grid summariser.
 *
 * Most cases use a grid whose answer is known analytically, so a wrong
 * quadrature shows up as a number rather than as "different from before".
 */
#include "test_runner.h"
#include "gd_grid.h"

#include <math.h>
#include <stdlib.h>

/* Fill t[] evenly over [lo, hi] and logdens[] from a normal log-density. */
static void make_normal(double *t, double *logdens, size_t n,
                        double lo, double hi, double mu, double sigma)
{
    double step = (hi - lo) / (double)(n - 1);
    for (size_t i = 0; i < n; i++) {
        t[i] = lo + (double)i * step;
        double z = (t[i] - mu) / sigma;
        logdens[i] = -0.5 * z * z;      /* unnormalised on purpose */
    }
}

TEST(test_grid_uniform_moments)
{
    /* A flat density on [0, 10]: mean 5, sd = 10/sqrt(12), median 5. */
    enum { N = 1001 };
    double t[N], ld[N];
    for (size_t i = 0; i < N; i++) {
        t[i] = 10.0 * (double)i / (double)(N - 1);
        ld[i] = 0.0;
    }
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    CHECK_APPROX(s.mean, 5.0, 1e-9);
    CHECK_APPROX(s.median, 5.0, 1e-9);
    CHECK_APPROX(s.sd, 10.0 / sqrt(12.0), 1e-3);
    CHECK_APPROX(s.et_lo_90, 0.5, 1e-6);
    CHECK_APPROX(s.et_hi_90, 9.5, 1e-6);
}

TEST(test_grid_normal_recovers_mean_and_sd)
{
    enum { N = 2001 };
    double t[N], ld[N];
    make_normal(t, ld, N, 0.0, 20.0, 10.0, 2.0);   /* well inside the grid */
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    CHECK_APPROX(s.mean, 10.0, 1e-6);
    CHECK_APPROX(s.median, 10.0, 1e-6);
    CHECK_APPROX(s.mode, 10.0, 0.02);
    CHECK_APPROX(s.sd, 2.0, 1e-4);
    /* 90% equal-tailed is mu +- 1.6449 sigma. */
    CHECK_APPROX(s.et_lo_90, 10.0 - 1.6449 * 2.0, 1e-2);
    CHECK_APPROX(s.et_hi_90, 10.0 + 1.6449 * 2.0, 1e-2);
}

TEST(test_grid_unnormalised_matches_normalised)
{
    /* Every quantity divides through by the mass, so shifting the
     * log-density by a constant must change nothing. */
    enum { N = 501 };
    double t[N], ld[N], ld_shifted[N];
    make_normal(t, ld, N, 0.0, 20.0, 8.0, 3.0);
    for (size_t i = 0; i < N; i++) ld_shifted[i] = ld[i] + 137.0;

    gd_grid_summary a, b;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &a, &err) == GD_E_OK);
    CHECK(gd_grid_summarize(t, ld_shifted, N, &b, &err) == GD_E_OK);
    CHECK_APPROX(a.mean, b.mean, 1e-12);
    CHECK_APPROX(a.sd, b.sd, 1e-12);
    CHECK_APPROX(a.median, b.median, 1e-12);
    CHECK_APPROX(a.hpd_lo_90, b.hpd_lo_90, 1e-12);
    CHECK_APPROX(a.hpd_hi_90, b.hpd_hi_90, 1e-12);
}

TEST(test_grid_quantiles_are_ordered)
{
    enum { N = 801 };
    double t[N], ld[N];
    make_normal(t, ld, N, 0.0, 20.0, 9.0, 2.5);
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    CHECK(s.et_lo_95 <= s.et_lo_90);
    CHECK(s.et_lo_90 <= s.et_lo_80);
    CHECK(s.et_lo_80 <= s.et_lo_50);
    CHECK(s.et_lo_50 <= s.median);
    CHECK(s.median   <= s.et_hi_50);
    CHECK(s.et_hi_50 <= s.et_hi_80);
    CHECK(s.et_hi_80 <= s.et_hi_90);
    CHECK(s.et_hi_90 <= s.et_hi_95);
}

TEST(test_grid_hpd_is_narrower_than_equal_tailed_when_skewed)
{
    /* A strongly right-skewed density: the HPD should be narrower than the
     * equal-tailed interval at the same level, and shifted toward the mode.
     * Nothing tested this before D2 — only P0 ever filled the HPD fields,
     * and no shipped model called P0. */
    enum { N = 2001 };
    double t[N], ld[N];
    double lo = 0.05, hi = 20.0;
    double step = (hi - lo) / (double)(N - 1);
    for (size_t i = 0; i < N; i++) {
        t[i] = lo + (double)i * step;
        /* log-normal(0, 0.75): heavy right tail. */
        double z = (log(t[i]) - 0.0) / 0.75;
        ld[i] = -0.5 * z * z - log(t[i]);
    }
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);

    double hpd_w = s.hpd_hi_90 - s.hpd_lo_90;
    double et_w  = s.et_hi_90  - s.et_lo_90;
    CHECK(hpd_w > 0.0);
    CHECK(hpd_w < et_w);
    /* Skewed right, so the HPD sits below the equal-tailed interval. */
    CHECK(s.hpd_lo_90 < s.et_lo_90);
    CHECK(s.hpd_hi_90 < s.et_hi_90);
    /* Wider levels contain narrower ones. */
    CHECK(s.hpd_lo_95 <= s.hpd_lo_90);
    CHECK(s.hpd_hi_95 >= s.hpd_hi_90);
    CHECK(s.hpd_lo_90 <= s.hpd_lo_80);
    CHECK(s.hpd_hi_90 >= s.hpd_hi_80);
}

TEST(test_grid_hpd_equals_equal_tailed_when_symmetric)
{
    /* On a symmetric density the two intervals agree to within a grid cell. */
    enum { N = 2001 };
    double t[N], ld[N];
    make_normal(t, ld, N, 0.0, 20.0, 10.0, 2.0);
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    double cell = 20.0 / (double)(N - 1);
    CHECK_APPROX(s.hpd_lo_90, s.et_lo_90, 5.0 * cell);
    CHECK_APPROX(s.hpd_hi_90, s.et_hi_90, 5.0 * cell);
}

TEST(test_grid_point_mass_collapses)
{
    /* All mass on one cell: every summary should land on it, and nothing
     * should divide by zero. */
    enum { N = 101 };
    double t[N], ld[N];
    for (size_t i = 0; i < N; i++) {
        t[i] = (double)i;
        ld[i] = -HUGE_VAL;
    }
    ld[40] = 0.0;
    ld[41] = 0.0;
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    CHECK(s.mean >= 40.0 && s.mean <= 41.0);
    CHECK(s.median >= 40.0 && s.median <= 41.0);
    CHECK(s.mode >= 40.0 && s.mode <= 41.0);
    CHECK(s.sd >= 0.0 && s.sd < 1.0);
}

TEST(test_grid_all_zero_density_does_not_divide_by_zero)
{
    enum { N = 64 };
    double t[N], ld[N];
    for (size_t i = 0; i < N; i++) { t[i] = (double)i; ld[i] = -HUGE_VAL; }
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    /* No mass anywhere: the summariser reports the grid rather than NaN. */
    CHECK(s.mean == 0.0);
    CHECK(isfinite(s.sd));
    CHECK(s.hpd_lo_90 == t[0]);
    CHECK(s.hpd_hi_90 == t[N - 1]);
}

TEST(test_grid_rejects_bad_arguments)
{
    double t[4] = {0, 1, 2, 3};
    double ld[4] = {0, 0, 0, 0};
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(NULL, ld, 4, &s, &err) == GD_E_INVALID_ARG);
    CHECK(gd_grid_summarize(t, NULL, 4, &s, &err) == GD_E_INVALID_ARG);
    CHECK(gd_grid_summarize(t, ld, 4, NULL, &err) == GD_E_INVALID_ARG);
    CHECK(gd_grid_summarize(t, ld, 1, &s, &err) == GD_E_INVALID_ARG);
    /* A NULL error is always allowed; the code still comes back. */
    CHECK(gd_grid_summarize(t, ld, 0, &s, NULL) == GD_E_INVALID_ARG);
}

TEST(test_grid_handles_uneven_spacing)
{
    /* Spacing need not be even: a flat density on a stretched grid still
     * has its median at the midpoint of the range by mass. */
    enum { N = 5 };
    double t[N]  = {0.0, 1.0, 2.0, 6.0, 10.0};
    double ld[N] = {0.0, 0.0, 0.0, 0.0, 0.0};
    gd_grid_summary s;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_grid_summarize(t, ld, N, &s, &err) == GD_E_OK);
    CHECK_APPROX(s.mean, 5.0, 1e-9);
    CHECK_APPROX(s.median, 5.0, 1e-9);
}

static void run_tests(void)
{
    RUN(test_grid_uniform_moments);
    RUN(test_grid_normal_recovers_mean_and_sd);
    RUN(test_grid_unnormalised_matches_normalised);
    RUN(test_grid_quantiles_are_ordered);
    RUN(test_grid_hpd_is_narrower_than_equal_tailed_when_skewed);
    RUN(test_grid_hpd_equals_equal_tailed_when_symmetric);
    RUN(test_grid_point_mass_collapses);
    RUN(test_grid_all_zero_density_does_not_divide_by_zero);
    RUN(test_grid_rejects_bad_arguments);
    RUN(test_grid_handles_uneven_spacing);
}

TEST_MAIN()

/*
 * test_loo.c — PSIS-LOO and stacking (src/gd_loo.c) against known answers.
 */
#include "test_runner.h"
#include "gd_loo.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>

static int dcmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* The GPD fit recovers the shape of exact GPD quantiles. */
TEST(test_gpd_fit_recovers_shape)
{
    const double ks[3] = { 0.2, 0.5, 0.9 };
    for (int t = 0; t < 3; t++) {
        size_t n = 4000;
        double *x = (double *)malloc(n * sizeof(double));
        for (size_t i = 0; i < n; i++) {
            double p = ((double)i + 0.5) / (double)n;
            x[i] = 2.0 * expm1(-ks[t] * log1p(-p)) / ks[t];
        }
        qsort(x, n, sizeof(double), dcmp);
        double sigma, k = gd_gpd_fit(x, n, &sigma);
        printf("    GPD k %.2f: fitted %.3f, sigma %.3f (want 2)\n", ks[t], k, sigma);
        CHECK_APPROX(k, ks[t], 0.05);
        CHECK_APPROX(sigma, 2.0, 0.15);
        free(x);
    }
}

/*
 * Conjugate normal model: y_i ~ N(mu, 1), mu ~ N(0, 10^2). Exact LOO
 * predictive for y_i is N(m_{-i}, 1 + v_{-i}) from the other n-1 points.
 * PSIS-LOO from exact posterior draws must match it.
 */
TEST(test_psis_loo_matches_exact_loo)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(11u, 1u, &err);
    const size_t N = 40, S = 4000;
    double y[40], sum = 0.0;
    for (size_t i = 0; i < N; i++) { y[i] = 1.0 + gd_rng_normal(rng, 0, 1); sum += y[i]; }
    y[0] = 5.0;                                    /* one surprising point */
    sum = 0.0; for (size_t i = 0; i < N; i++) sum += y[i];
    double prec = 1.0 / 100.0 + (double)N, v = 1.0 / prec, m = v * sum;
    double *ll = (double *)malloc(S * N * sizeof(double));
    for (size_t s = 0; s < S; s++) {
        double mu = gd_rng_normal(rng, m, sqrt(v));
        for (size_t i = 0; i < N; i++)
            ll[s * N + i] = -0.5 * log(2 * 3.14159265358979323846) - 0.5 * (y[i] - mu) * (y[i] - mu);
    }
    double elpd[40], k[40], tot = 0.0, exact = 0.0, worst = 0.0, kmax = 0.0;
    CHECK(gd_psis_loo(ll, S, N, elpd, k) == 0);
    for (size_t i = 0; i < N; i++) {
        double pv = 1.0 / (1.0 / 100.0 + (double)(N - 1)), pm = pv * (sum - y[i]);
        double sd2 = 1.0 + pv;
        double ex = -0.5 * log(2 * 3.14159265358979323846 * sd2) - 0.5 * (y[i] - pm) * (y[i] - pm) / sd2;
        tot += elpd[i]; exact += ex;
        if (fabs(elpd[i] - ex) > worst) worst = fabs(elpd[i] - ex);
        if (k[i] > kmax) kmax = k[i];
    }
    printf("    PSIS-LOO %.3f vs exact %.3f (worst point %.4f, max k %.2f)\n", tot, exact, worst, kmax);
    CHECK_APPROX(tot, exact, 0.05);
    CHECK(worst < 0.02);
    CHECK(kmax < 0.7);
    free(ll);
    gd_rng_destroy(rng);
}

/* Stacking: a model better on every concept takes all the weight; two
 * models each perfect on half the concepts split it evenly; the optimum of
 * a two-model, two-concept case matches its closed form. */
TEST(test_stacking_weights)
{
    double w[2];
    double dom[2 * 4] = { -1, -1, -1, -1,   -2, -2, -2, -2 };
    gd_stacking_weights(dom, 2, 4, w);
    CHECK(w[0] > 0.999);
    double half[2 * 4] = { 0, 0, -30, -30,   -30, -30, 0, 0 };
    gd_stacking_weights(half, 2, 4, w);
    CHECK_APPROX(w[0], 0.5, 1e-6);
    /* Concept 1: densities (a1, b1) = (0.8, 0.2); concept 2: (0.3, 0.6).
     * Maximise log(0.2 + 0.6 w) + log(0.6 - 0.3 w): w* = (0.6*0.6 - 0.3*0.2) /
     * (2 * 0.6 * 0.3) = 0.8333. */
    double cf[2 * 2] = { log(0.8), log(0.3),   log(0.2), log(0.6) };
    gd_stacking_weights(cf, 2, 2, w);
    printf("    stacking closed form: w %.5f (want 0.83333)\n", w[0]);
    CHECK_APPROX(w[0], 0.3 / 0.36, 1e-5);
    CHECK_APPROX(w[0] + w[1], 1.0, 1e-12);
}

/* A concept every draw scores the same (an excluded one scores 0) has an
 * exact LOO estimate: uniform weights, k = 0, elpd = the common value. */
TEST(test_psis_constant_ratios)
{
    double ll[200 * 2];
    for (size_t s = 0; s < 200; s++) { ll[s * 2] = 0.0; ll[s * 2 + 1] = -1.25; }
    double elpd[2], k[2];
    CHECK(gd_psis_loo(ll, 200, 2, elpd, k) == 0);
    CHECK(k[0] == 0.0 && k[1] == 0.0);
    CHECK_APPROX(elpd[0], 0.0, 1e-12);
    CHECK_APPROX(elpd[1], -1.25, 1e-12);
}

static void run_tests(void)
{
    RUN(test_psis_constant_ratios);
    RUN(test_gpd_fit_recovers_shape);
    RUN(test_psis_loo_matches_exact_loo);
    RUN(test_stacking_weights);
}

TEST_MAIN()

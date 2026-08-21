/*
 * test_pairwise.c — tests for the P0 pairwise engine (M2).
 *
 * Verifies:
 *  - agreement with the classical Swadesh/Lees formula;
 *  - posterior normalization (integral = 1);
 *  - quantiles behave correctly;
 *  - simulation round-trips through the likelihood;
 *  - uncertain cognacy and contamination.
 */
#include "test_runner.h"
#include "gd_pairwise.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <string.h>

TEST(test_p0_loglik_match)
{
    /* P(same) = exp(-r*d). For match, loglik = log(exp(-r*d)) = -r*d. */
    double ll = gd_p0_loglik_concept(0.1, 10.0, 1, 1.0, 0.0, 0.5, NULL);
    CHECK_APPROX(ll, -1.0, 1e-12);  /* -0.1 * 10 = -1.0 */
}

TEST(test_p0_loglik_mismatch)
{
    /* P(mismatch) = 1 - exp(-r*d). loglik = log(1 - exp(-r*d)). */
    double r = 0.1, d = 10.0;
    double expected = log(1.0 - exp(-r * d));
    double ll = gd_p0_loglik_concept(r, d, 0, 1.0, 0.0, 0.5, NULL);
    CHECK_APPROX(ll, expected, 1e-12);
}

TEST(test_p0_loglik_uncertain_cognacy)
{
    /* y = 0.5, match: L = 0.5 * p + 0.5 * (1-p) = 0.5 (for any p). */
    double ll = gd_p0_loglik_concept(0.1, 10.0, 1, 0.5, 0.0, 0.5, NULL);
    CHECK_APPROX(ll, log(0.5), 1e-12);
}

TEST(test_p0_loglik_contamination)
{
    /* eps=0.3, q=0.5: p_match = 0.7*exp(-rd) + 0.3*0.5. */
    double r = 0.1, d = 10.0, eps = 0.3, q = 0.5;
    double s = exp(-r * d);
    double p_match = (1.0 - eps) * s + eps * q;
    double expected = log(p_match);  /* match, y=1 */
    double ll = gd_p0_loglik_concept(r, d, 1, 1.0, eps, q, NULL);
    CHECK_APPROX(ll, expected, 1e-12);
}

TEST(test_p0_posterior_normalization)
{
    /* The posterior density must integrate to 1. */
    double rates[] = {0.15, 0.15, 0.15, 0.15, 0.15};
    int matches[] = {1, 1, 0, 1, 0};
    gd_error err = GD_ERROR_INIT;
    gd_pair_result result;
    int rc = gd_p0_posterior(rates, matches, NULL, 5, 0, 0, 100, 5000,
                             0.0, NULL, 256, &result, &err);
    CHECK(rc == 0);

    /* Integrate the density via trapezoid. */
    double total = 0.0;
    for (size_t i = 1; i < result.grid_n; i++) {
        double dT = result.grid_T[i] - result.grid_T[i - 1];
        double f0 = result.grid_density[i - 1];
        double f1 = result.grid_density[i];
        total += 0.5 * (f0 + f1) * dT;
    }
    CHECK_APPROX(total, 1.0, 1e-6);

    free((void *)result.grid_T);
    free((void *)result.grid_density);
}

TEST(test_p0_posterior_classical_formula)
{
    /* Classical glottochronology: for n concepts with m matches, the MLE
     * of T is T_hat = -m / (n * r * 2) ... actually the classical formula
     * gives the retention fraction: m/n ≈ exp(-2*r*T), so T ≈ -log(m/n)/(2r).
     *
     * With 100 concepts, 86 matches, r=0.1508228897/kyr (Swadesh-100):
     *   r = 0.1508228897/kyr = 0.0001508228897/year.
     *   T ≈ -log(0.86) / (2 * 0.0001508228897) ≈ 500 years.
     *
     * The posterior mode should be near 500. */
    size_t n = 100;
    double rates[100];
    int matches[100];
    double r_per_year = 0.1508228897 / 1000.0;
    for (size_t i = 0; i < n; i++) {
        rates[i] = r_per_year;
        matches[i] = (i < 86) ? 1 : 0;  /* 86 matches */
    }
    gd_error err = GD_ERROR_INIT;
    gd_pair_result result;
    int rc = gd_p0_posterior(rates, matches, NULL, n, 0, 0, 100, 5000,
                             0.0, NULL, 512, &result, &err);
    CHECK(rc == 0);

    double T_classical = -log(86.0 / 100.0) / (2.0 * r_per_year);
    /* The posterior mode should be close to the classical estimate. */
    CHECK_APPROX(result.mode, T_classical, 50.0);  /* within 50 years */

    free((void *)result.grid_T);
    free((void *)result.grid_density);
}

TEST(test_p0_posterior_quantiles_ordered)
{
    double rates[] = {0.15, 0.20, 0.10, 0.25, 0.18};
    int matches[] = {1, 0, 1, 0, 1};
    gd_error err = GD_ERROR_INIT;
    gd_pair_result result;
    int rc = gd_p0_posterior(rates, matches, NULL, 5, 0, 0, 50, 5000,
                             0.0, NULL, 256, &result, &err);
    CHECK(rc == 0);
    /* Quantiles must be ordered: wider intervals contain narrower ones.
     * et_lo_95 <= et_lo_90 <= et_lo_80 <= et_lo_50 <= median
     * median <= et_hi_50 <= et_hi_80 <= et_hi_90 <= et_hi_95 */
    CHECK(result.et_lo_95 <= result.et_lo_90);
    CHECK(result.et_lo_90 <= result.et_lo_80);
    CHECK(result.et_lo_80 <= result.et_lo_50);
    CHECK(result.et_lo_50 <= result.median);
    CHECK(result.median <= result.et_hi_50);
    CHECK(result.et_hi_50 <= result.et_hi_80);
    CHECK(result.et_hi_80 <= result.et_hi_90);
    CHECK(result.et_hi_90 <= result.et_hi_95);
    /* Mean and median should be positive and within range. */
    CHECK(result.mean > 0.0);
    CHECK(result.sd > 0.0);
    free((void *)result.grid_T);
    free((void *)result.grid_density);
}

TEST(test_p0_simulate_roundtrip)
{
    /* Simulate data with known T and r, then check the posterior recovers T. */
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(42, 1, &err);
    CHECK(rng != NULL);

    double r = 0.15 / 1000.0; /* per year (= 0.15/kyr) */
    double T_true = 1000.0;  /* years */
    double d = 2.0 * T_true;  /* 2000 */
    size_t n = 200;
    int matches[200];
    double rates[200];
    for (size_t i = 0; i < n; i++) rates[i] = r;

    int rc = gd_p0_simulate(rng, r, d, n, matches, &err);
    CHECK(rc == 0);

    /* Expected matches ≈ n * exp(-r*d) = 200 * exp(-0.00015*2000)
     * = 200 * exp(-0.3) ≈ 200 * 0.7408 ≈ 148. */
    /* The posterior should put mass near T=1000. */
    gd_pair_result result;
    rc = gd_p0_posterior(rates, matches, NULL, n, 0, 0, 100, 5000,
                         0.0, NULL, 512, &result, &err);
    CHECK(rc == 0);
    /* The posterior mean should be within ~30% of T_true. */
    CHECK(fabs(result.mean - T_true) < 0.3 * T_true);

    free((void *)result.grid_T);
    free((void *)result.grid_density);
    gd_rng_destroy(rng);
}

TEST(test_p0_zero_path_length)
{
    /* T = max(a1, a2): path length d = 0. P(match) = 1. */
    double ll = gd_p0_loglik_concept(0.1, 0.0, 1, 1.0, 0.0, 0.5, NULL);
    CHECK_APPROX(ll, 0.0, 1e-12);  /* log(1) = 0 */
    /* Mismatch at d=0: P = 0, log(0) = -inf. */
    ll = gd_p0_loglik_concept(0.1, 0.0, 0, 1.0, 0.0, 0.5, NULL);
    CHECK(ll == -HUGE_VAL);
}

TEST(test_p0_extreme_old_date)
{
    /* Very large T: d → ∞, exp(-r*d) → 0. Match → log(0) = -inf.
     * Mismatch → log(1) = 0. */
    double ll = gd_p0_loglik_concept(0.1, 1e15, 1, 1.0, 0.0, 0.5, NULL);
    CHECK(ll == -HUGE_VAL);
    ll = gd_p0_loglik_concept(0.1, 1e15, 0, 1.0, 0.0, 0.5, NULL);
    CHECK_APPROX(ll, 0.0, 1e-12);
}

TEST(test_p0_very_slow_concept)
{
    /* r → 0: nearly no change. P(match) ≈ 1 even for large d. */
    double ll = gd_p0_loglik_concept(1e-10, 1000.0, 1, 1.0, 0.0, 0.5, NULL);
    /* exp(-1e-10 * 1000) = exp(-1e-7) ≈ 1 - 1e-7. log ≈ -1e-7. */
    CHECK_APPROX(ll, -1e-7, 1e-9);
}

TEST(test_p0_pairwise_symmetry)
{
    /* Swapping a1 and a2 should give the same result (d = 2T - a1 - a2). */
    double rates[] = {0.1, 0.2, 0.15};
    int matches[] = {1, 0, 1};
    gd_error err = GD_ERROR_INIT;
    gd_pair_result r1, r2;
    int rc = gd_p0_posterior(rates, matches, NULL, 3, 100, 300, 400, 5000,
                             0.0, NULL, 256, &r1, &err);
    CHECK(rc == 0);
    rc = gd_p0_posterior(rates, matches, NULL, 3, 300, 100, 400, 5000,
                         0.0, NULL, 256, &r2, &err);
    CHECK(rc == 0);
    CHECK_APPROX(r1.mean, r2.mean, 1e-6);
    CHECK_APPROX(r1.median, r2.median, 1e-6);
    free((void *)r1.grid_T); free((void *)r1.grid_density);
    free((void *)r2.grid_T); free((void *)r2.grid_density);
}

TEST(test_p0_threshold_probability)
{
    /* P(T > threshold) should be computable from the CDF. We check
     * that the CDF at the last grid point is 1 (or very close). */
    double rates[] = {0.00015, 0.00020, 0.00010, 0.00018};
    int matches[] = {1, 0, 1, 0};
    gd_error err = GD_ERROR_INIT;
    gd_pair_result result;
    int rc = gd_p0_posterior(rates, matches, NULL, 4, 0, 0, 100, 5000,
                             0.0, NULL, 256, &result, &err);
    CHECK(rc == 0);
    /* Integrate density from T_lo to some threshold and check 0 < mass < 1. */
    double mass_below = 0.0, total = 0.0;
    for (size_t i = 1; i < result.grid_n; i++) {
        double dT = result.grid_T[i] - result.grid_T[i - 1];
        double f0 = result.grid_density[i - 1];
        double f1 = result.grid_density[i];
        double m = 0.5 * (f0 + f1) * dT;
        total += m;
        if (result.grid_T[i] < 2000.0) mass_below += m;
    }
    CHECK_APPROX(total, 1.0, 1e-6);
    CHECK(mass_below > 0.0 && mass_below < 1.0);
    free((void *)result.grid_T);
    free((void *)result.grid_density);
}

static void run_tests(void)
{
    RUN(test_p0_loglik_match);
    RUN(test_p0_loglik_mismatch);
    RUN(test_p0_loglik_uncertain_cognacy);
    RUN(test_p0_loglik_contamination);
    RUN(test_p0_posterior_normalization);
    RUN(test_p0_posterior_classical_formula);
    RUN(test_p0_posterior_quantiles_ordered);
    RUN(test_p0_simulate_roundtrip);
    RUN(test_p0_zero_path_length);
    RUN(test_p0_extreme_old_date);
    RUN(test_p0_very_slow_concept);
    RUN(test_p0_pairwise_symmetry);
    RUN(test_p0_threshold_probability);
}

TEST_MAIN()

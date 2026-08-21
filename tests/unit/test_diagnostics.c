/*
 * test_diagnostics.c — MCMC convergence diagnostics (M8).
 *
 * Checks the inverse-normal quantile, rank-normalized split-R-hat, and bulk
 * ESS on constructed cases with known answers: converged chains give R-hat ~ 1
 * and ESS ~ the total draw count; a stuck / offset chain is detected by a large
 * R-hat; strong autocorrelation collapses the ESS.
 */
#include "test_runner.h"
#include "gd_diagnostics.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>

TEST(test_inv_normal_cdf_known)
{
    CHECK_APPROX(gd_inv_normal_cdf(0.5), 0.0, 1e-12);
    CHECK_APPROX(gd_inv_normal_cdf(0.975), 1.959963984540054, 1e-9);
    CHECK_APPROX(gd_inv_normal_cdf(0.025), -1.959963984540054, 1e-9);
    CHECK_APPROX(gd_inv_normal_cdf(0.8413447460685429), 1.0, 1e-9);
    CHECK_APPROX(gd_inv_normal_cdf(0.1586552539314571), -1.0, 1e-9);
    CHECK(gd_inv_normal_cdf(0.0) == -HUGE_VAL);
    CHECK(gd_inv_normal_cdf(1.0) == HUGE_VAL);
    CHECK(isnan(gd_inv_normal_cdf(-0.1)));
    CHECK(isnan(gd_inv_normal_cdf(1.1)));
}

TEST(test_rhat_converged_vs_stuck)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(7u, 1u, &err);
    CHECK(rng != NULL);
    size_t nc = 4, n = 500;
    double *conv = (double *)malloc(nc * n * sizeof(double));
    double *stuck = (double *)malloc(nc * n * sizeof(double));
    for (size_t c = 0; c < nc; c++) {
        for (size_t i = 0; i < n; i++) {
            double x = gd_rng_normal(rng, 0.0, 1.0);
            conv[c * n + i] = x;
            /* One chain sits at a different location: an unconverged ensemble. */
            stuck[c * n + i] = (c == 0) ? x + 6.0 : x;
        }
    }
    double rc = gd_rhat_rank(conv, nc, n);
    double rs = gd_rhat_rank(stuck, nc, n);
    CHECK(rc < 1.05);        /* converged: near 1 */
    CHECK(rs > 1.2);         /* one offset chain is flagged */
    /* Classic (non-rank) split-R-hat behaves the same way here. */
    CHECK(gd_split_rhat(conv, nc, n) < 1.05);
    CHECK(gd_split_rhat(stuck, nc, n) > 1.2);
    /* Degenerate inputs. */
    CHECK(isnan(gd_rhat_rank(conv, 1, n)));       /* need >= 2 chains */
    CHECK(isnan(gd_rhat_rank(conv, nc, 3)));      /* need >= 4 draws */
    free(conv); free(stuck);
    gd_rng_destroy(rng);
}

TEST(test_ess_iid_vs_autocorrelated)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(11u, 2u, &err);
    CHECK(rng != NULL);
    size_t nc = 4, n = 1000;
    size_t total = nc * n;
    double *iid = (double *)malloc(total * sizeof(double));
    double *ar1 = (double *)malloc(total * sizeof(double));
    double phi = 0.9;
    for (size_t c = 0; c < nc; c++) {
        double prev = gd_rng_normal(rng, 0.0, 1.0);
        for (size_t i = 0; i < n; i++) {
            double e = gd_rng_normal(rng, 0.0, 1.0);
            iid[c * n + i] = e;
            prev = phi * prev + sqrt(1.0 - phi * phi) * e;
            ar1[c * n + i] = prev;
        }
    }
    double ess_iid = gd_ess_bulk(iid, nc, n);
    double ess_ar1 = gd_ess_bulk(ar1, nc, n);
    /* iid draws: ESS is a large fraction of the total. */
    CHECK(ess_iid > 0.6 * (double)total);
    CHECK(ess_iid < 1.4 * (double)total);
    /* AR(1) phi=0.9: theoretical ESS ~ total*(1-phi)/(1+phi) ~ 0.053*total. */
    CHECK(ess_ar1 < 0.25 * (double)total);
    CHECK(ess_ar1 > 0.0);
    free(iid); free(ar1);
    gd_rng_destroy(rng);
}

static void run_tests(void)
{
    RUN(test_inv_normal_cdf_known);
    RUN(test_rhat_converged_vs_stuck);
    RUN(test_ess_iid_vs_autocorrelated);
}

TEST_MAIN()

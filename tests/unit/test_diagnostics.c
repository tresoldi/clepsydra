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
#include <string.h>

#include "diag_golden.h"

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

/* E0.9: every estimator agrees with ArviZ to rounding, on cases chosen to
 * exercise the corners: an odd chain length (the middle draw is dropped on
 * split), strong autocorrelation, one offset chain, a scale mismatch that only
 * the folded R-hat sees, ties, and antithetic draws (ESS above N). */
TEST(test_diagnostics_match_arviz)
{
    size_t ncase = sizeof GOLDEN / sizeof GOLDEN[0];
    for (size_t k = 0; k < ncase; k++) {
        const golden_case *g = &GOLDEN[k];
        double r  = gd_rhat_rank(g->draws, GOLDEN_C, GOLDEN_N);
        double eb = gd_ess_bulk(g->draws, GOLDEN_C, GOLDEN_N);
        double et = gd_ess_tail(g->draws, GOLDEN_C, GOLDEN_N);
        double em = gd_ess_mean(g->draws, GOLDEN_C, GOLDEN_N);
        double mc = gd_mcse_mean(g->draws, GOLDEN_C, GOLDEN_N);
        CHECK_APPROX(r,  g->rhat,      1e-9 * g->rhat);
        CHECK_APPROX(eb, g->ess_bulk,  1e-7 * g->ess_bulk);
        CHECK_APPROX(et, g->ess_tail,  1e-7 * g->ess_tail);
        CHECK_APPROX(em, g->ess_mean,  1e-7 * g->ess_mean);
        CHECK_APPROX(mc, g->mcse_mean, 1e-7 * g->mcse_mean);
    }
}

/* The folded R-hat is what catches chains that agree in location and
 * disagree in spread; the bulk R-hat alone does not. */
TEST(test_rhat_tail_catches_scale_mismatch)
{
    const golden_case *g = NULL;
    for (size_t k = 0; k < sizeof GOLDEN / sizeof GOLDEN[0]; k++)
        if (strcmp(GOLDEN[k].name, "scale_mismatch") == 0) g = &GOLDEN[k];
    CHECK(g != NULL);
    CHECK(gd_rhat_bulk(g->draws, GOLDEN_C, GOLDEN_N) < 1.05);
    CHECK(gd_rhat_tail(g->draws, GOLDEN_C, GOLDEN_N) > 1.1);
}

static void run_tests(void)
{
    RUN(test_diagnostics_match_arviz);
    RUN(test_rhat_tail_catches_scale_mismatch);
    RUN(test_inv_normal_cdf_known);
    RUN(test_rhat_converged_vs_stuck);
    RUN(test_ess_iid_vs_autocorrelated);
}

TEST_MAIN()

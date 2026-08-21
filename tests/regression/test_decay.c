/*
 * test_decay.c — alternative pairwise decay laws (M12 extension).
 *
 * The Swadesh / Starostin / Vasilyev / two-component retention laws ported from
 * distglot, each wrapped in a Bayesian posterior. Checks the retention-curve
 * properties, round-trip inversion, posterior recovery, and — the point of
 * porting them — that the competing laws give *different* dates on the same
 * data, so they can be compared rather than assumed.
 */
#include "test_runner.h"
#include "gd_decay.h"
#include "gd_memory.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>

TEST(test_decay_curve_properties)
{
    gd_decay_params p; gd_decay_params_default(&p);
    gd_decay_model models[4] = { GD_DECAY_SWADESH, GD_DECAY_STAROSTIN,
                                 GD_DECAY_VASILYEV, GD_DECAY_TWO_COMPONENT };
    for (int k = 0; k < 4; k++) {
        /* s(0) = 1 for all laws. */
        CHECK_APPROX(gd_decay_similarity(models[k], 0.0, &p), 1.0, 1e-9);
        /* Monotone non-increasing and within [0,1]. */
        double prev = 1.0;
        for (double T = 0.0; T <= 15000.0; T += 250.0) {
            double s = gd_decay_similarity(models[k], T, &p);
            CHECK(s <= prev + 1e-12 && s >= -1e-12 && s <= 1.0 + 1e-12);
            prev = s;
        }
    }
    /* Vasilyev has an ~8% ultra-stable floor; Swadesh decays to ~0. */
    CHECK(gd_decay_similarity(GD_DECAY_VASILYEV, 1e6, &p) > 0.079);
    CHECK(gd_decay_similarity(GD_DECAY_SWADESH, 20000.0, &p) < 0.05);
    /* Two-component floor = w_s * stable_base. */
    CHECK_APPROX(gd_decay_similarity(GD_DECAY_TWO_COMPONENT, 1e7, &p),
                 p.tc_stable_weight * p.tc_stable_base, 1e-6);
}

TEST(test_decay_inversion_roundtrip)
{
    gd_decay_params p; gd_decay_params_default(&p);
    gd_decay_model models[3] = { GD_DECAY_SWADESH, GD_DECAY_STAROSTIN, GD_DECAY_VASILYEV };
    for (int k = 0; k < 3; k++) {
        for (double T = 200.0; T <= 6000.0; T += 400.0) {
            double s = gd_decay_similarity(models[k], T, &p);
            double Tinv = gd_decay_invert(models[k], s, &p, 20000.0);
            CHECK(fabs(Tinv - T) < 5.0);          /* years */
        }
    }
    /* A similarity below the Vasilyev floor is unreachable. */
    CHECK(isnan(gd_decay_invert(GD_DECAY_VASILYEV, 0.05, &p, 20000.0)));
    CHECK(isnan(gd_decay_invert(GD_DECAY_SWADESH, 1.5, &p, 20000.0)));
}

TEST(test_decay_models_disagree)
{
    /* The whole reason to port them: at the same observed similarity the three
     * laws imply materially different divergence times. */
    gd_decay_params p; gd_decay_params_default(&p);
    double sim = 0.55;
    double tsw = gd_decay_invert(GD_DECAY_SWADESH, sim, &p, 20000.0);
    double tst = gd_decay_invert(GD_DECAY_STAROSTIN, sim, &p, 20000.0);
    double tva = gd_decay_invert(GD_DECAY_VASILYEV, sim, &p, 20000.0);
    CHECK(tsw > 0 && tst > 0 && tva > 0);
    CHECK(fabs(tsw - tst) > 200.0);
    CHECK(fabs(tsw - tva) > 200.0);
    printf("    similarity %.2f -> swadesh=%.0f starostin=%.0f vasilyev=%.0f BP\n",
           sim, tsw, tst, tva);
}

TEST(test_decay_posterior_recovers_and_normalizes)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(2024u, 9u, &err);
    CHECK(rng != NULL);
    gd_decay_params p; gd_decay_params_default(&p);

    const size_t NCON = 200;
    const double T_true = 1500.0;
    double s_true = gd_decay_similarity(GD_DECAY_SWADESH, T_true, &p);
    int *matches = (int *)malloc(NCON * sizeof(int));
    for (size_t c = 0; c < NCON; c++)
        matches[c] = (gd_rng_uniform_closed_open(rng) < s_true) ? 1 : 0;

    gd_pair_result r;
    CHECK(gd_decay_posterior(GD_DECAY_SWADESH, &p, matches, NULL, NCON,
                             100.0, 6000.0, 1024, &r, &err) == 0);
    /* Recovery: posterior mean near the truth. */
    CHECK(fabs(r.mean - T_true) < 0.25 * T_true);
    /* Ordered intervals. */
    CHECK(r.et_lo_90 < r.median && r.median < r.et_hi_90);
    /* Normalized density integrates to ~1. */
    double integ = 0.0;
    for (size_t i = 1; i < r.grid_n; i++) {
        double dT = r.grid_T[i] - r.grid_T[i - 1];
        integ += 0.5 * (r.grid_density[i - 1] + r.grid_density[i]) * dT;
    }
    CHECK_APPROX(integ, 1.0, 1e-6);
    gd_free((void *)r.grid_T);
    gd_free((void *)r.grid_density);

    /* The same data under different laws gives different posterior means. */
    gd_pair_result r_sw, r_va;
    CHECK(gd_decay_posterior(GD_DECAY_SWADESH, &p, matches, NULL, NCON,
                             100.0, 8000.0, 1024, &r_sw, &err) == 0);
    CHECK(gd_decay_posterior(GD_DECAY_VASILYEV, &p, matches, NULL, NCON,
                             100.0, 8000.0, 1024, &r_va, &err) == 0);
    CHECK(fabs(r_sw.mean - r_va.mean) > 100.0);
    gd_free((void *)r_sw.grid_T); gd_free((void *)r_sw.grid_density);
    gd_free((void *)r_va.grid_T); gd_free((void *)r_va.grid_density);

    free(matches);
    gd_rng_destroy(rng);
}

TEST(test_decay_posterior_fills_hpd_like_p1)
{
    /* Before D2 the decay laws left gd_pair_result's HPD fields at zero:
     * only P0 filled them, and no pipeline called P0. Every model goes
     * through gd_grid now, so the fields mean the same thing for all. */
    enum { NCONC = 40 };
    int m40[NCONC];
    for (size_t i = 0; i < NCONC; i++) m40[i] = (i % 3 == 0) ? 1 : 0;

    gd_decay_model laws[4] = { GD_DECAY_SWADESH, GD_DECAY_STAROSTIN,
                               GD_DECAY_VASILYEV, GD_DECAY_TWO_COMPONENT };
    for (size_t k = 0; k < 4; k++) {
        gd_pair_result r;
        gd_error err = GD_ERROR_INIT;
        int rc = gd_decay_posterior(laws[k], NULL, m40, NULL, NCONC,
                                    100.0, 8000.0, 1024, &r, &err);
        /* Copy out and release before asserting: CHECK returns from the test
         * on failure, so anything still held here would leak and ASan's
         * report would bury the assertion that actually fired. */
        double lo90 = 0, hi90 = 0, lo95 = 0, hi95 = 0, elo90 = 0, ehi90 = 0;
        if (rc == GD_E_OK) {
            lo90 = r.hpd_lo_90; hi90 = r.hpd_hi_90;
            lo95 = r.hpd_lo_95; hi95 = r.hpd_hi_95;
            elo90 = r.et_lo_90; ehi90 = r.et_hi_90;
            gd_free((void *)r.grid_T);
            gd_free((void *)r.grid_density);
        }

        CHECK(rc == GD_E_OK);
        CHECK(lo90 < hi90);
        CHECK(lo90 >= 100.0 && hi90 <= 8000.0);
        CHECK(lo95 <= lo90);
        CHECK(hi95 >= hi90);
        /* The HPD is never wider than the equal-tailed interval. */
        CHECK((hi90 - lo90) <= (ehi90 - elo90) + 1e-9);
    }
}

static void run_tests(void)
{
    RUN(test_decay_curve_properties);
    RUN(test_decay_inversion_roundtrip);
    RUN(test_decay_models_disagree);
    RUN(test_decay_posterior_recovers_and_normalizes);
    RUN(test_decay_posterior_fills_hpd_like_p1);
}

TEST_MAIN()

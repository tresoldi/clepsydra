/*
 * test_calibration.c — tests for calibrations and chronology (M5).
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_calibration.h"
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_rng.h"
#include "gd_math.h"
#include "gd.h"

#include <string.h>
#include <math.h>

/* strdup is not C99; use a simple helper. */
static char *dup_str(const char *s) {
    size_t n = strlen(s);
    char *p = (char *)malloc(n + 1);
    if (p) { memcpy(p, s, n); p[n] = 0; }
    return p;
}

/* Helper: build a simple Romance-style tree. */
static gd_tree *make_romance_tree(gd_error *err)
{
    const char *nwk = "((italian,spanish)romance_west,romanian)romance_root;";
    gd_tree *t = NULL;
    gd_newick_parse(nwk, strlen(nwk), &t, err);
    return t;
}

TEST(test_cal_type_from_string)
{
    gd_cal_type t;
    CHECK(gd_cal_type_from_string("fixed_tip", 9, &t) && t == GD_CAL_FIXED_TIP);
    CHECK(gd_cal_type_from_string("uniform_node", 12, &t) && t == GD_CAL_UNIFORM_NODE);
    CHECK(gd_cal_type_from_string("normal_node", 11, &t) && t == GD_CAL_NORMAL_NODE);
    CHECK(gd_cal_type_from_string("lognormal_node", 14, &t) && t == GD_CAL_LOGNORMAL_NODE);
    CHECK(gd_cal_type_from_string("offset_lognormal_node", 21, &t) && t == GD_CAL_OFFSET_LOGNORMAL_NODE);
    CHECK(gd_cal_type_from_string("hard_minimum", 12, &t) && t == GD_CAL_HARD_MINIMUM);
    CHECK(gd_cal_type_from_string("hard_maximum", 12, &t) && t == GD_CAL_HARD_MAXIMUM);

    /* An unrecognised name is refused and leaves the output alone, so a caller
     * that ignores the result cannot end up with a plausible-looking type. */
    t = GD_CAL_SOFT_MAXIMUM;
    CHECK(!gd_cal_type_from_string("unknown_type", 12, &t));
    CHECK(t == GD_CAL_SOFT_MAXIMUM);
}

TEST(test_cal_resolve_targets)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(1, &err);
    CHECK(cs != NULL);
    cs->cals[0].target_id = dup_str("romance_root");
    cs->cals[0].type = GD_CAL_UNIFORM_NODE;
    cs->cals[0].p1 = 2000;
    cs->cals[0].p2 = 5000;
    cs->cals[0].is_hard = false;

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == 0);
    CHECK(cs->cals[0].resolved == 1);

    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_resolve_missing_target)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(1, &err);
    CHECK(cs != NULL);
    cs->cals[0].target_id = dup_str("nonexistent");
    cs->cals[0].type = GD_CAL_UNIFORM_NODE;
    cs->cals[0].p1 = 2000;
    cs->cals[0].p2 = 5000;

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == GD_E_CALIBRATION);

    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_propagate_feasible)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(2, &err);
    CHECK(cs != NULL);
    /* Root: uniform [2000, 5000]. */
    cs->cals[0].target_id = dup_str("romance_root");
    cs->cals[0].type = GD_CAL_UNIFORM_NODE;
    cs->cals[0].p1 = 2000;
    cs->cals[0].p2 = 5000;
    /* Italian tip: fixed at 0 (modern). */
    cs->cals[1].target_id = dup_str("italian");
    cs->cals[1].type = GD_CAL_FIXED_TIP;
    cs->cals[1].p1 = 0;

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == 0);
    rc = gd_calibrations_propagate(cs, t, &err);
    CHECK(rc == 0);
    CHECK(cs->feasible == 1);

    /* Root should be >= 2000. */
    size_t root = gd_tree_find_label(t, "romance_root");
    CHECK(cs->node_lo[root] >= 2000);

    /* Italian tip should be exactly 0. */
    size_t it = gd_tree_find_label(t, "italian");
    CHECK(cs->node_lo[it] == 0);
    CHECK(cs->node_hi[it] == 0);

    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_contradiction_detected)
{
    /* Hard minimum on a child that exceeds the parent's hard maximum. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(2, &err);
    CHECK(cs != NULL);
    /* Parent (root): hard max 1000. */
    cs->cals[0].target_id = dup_str("romance_root");
    cs->cals[0].type = GD_CAL_HARD_MAXIMUM;
    cs->cals[0].p1 = 1000;
    /* Child (italian): hard min 2000 → contradiction. */
    cs->cals[1].target_id = dup_str("italian");
    cs->cals[1].type = GD_CAL_HARD_MINIMUM;
    cs->cals[1].p1 = 2000;

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == 0);
    rc = gd_calibrations_propagate(cs, t, &err);
    CHECK(rc == GD_E_CALIBRATION);
    CHECK(cs->feasible == 0);

    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_reject_invalid_normal_sd)
{
    /* normal_node with sd = 0 would divide by zero in the log-prior. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(1, &err);
    CHECK(cs != NULL);
    cs->cals[0].target_id = dup_str("romance_root");
    cs->cals[0].type = GD_CAL_NORMAL_NODE;
    cs->cals[0].p1 = 3000;
    cs->cals[0].p2 = 0;   /* invalid sd */

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == 0);
    rc = gd_calibrations_propagate(cs, t, &err);
    CHECK(rc == GD_E_CALIBRATION);

    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_reject_inverted_uniform)
{
    /* uniform_node with p1 >= p2 has no valid density. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(1, &err);
    CHECK(cs != NULL);
    cs->cals[0].target_id = dup_str("romance_root");
    cs->cals[0].type = GD_CAL_UNIFORM_NODE;
    cs->cals[0].p1 = 5000;
    cs->cals[0].p2 = 2000;   /* inverted */

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == 0);
    rc = gd_calibrations_propagate(cs, t, &err);
    CHECK(rc == GD_E_CALIBRATION);

    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_sample_prior_order)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    gd_calibration_set *cs = gd_calset_create(1, &err);
    CHECK(cs != NULL);
    cs->cals[0].target_id = dup_str("romance_root");
    cs->cals[0].type = GD_CAL_UNIFORM_NODE;
    cs->cals[0].p1 = 2000;
    cs->cals[0].p2 = 5000;

    int rc = gd_calibrations_resolve(cs, t, &err);
    CHECK(rc == 0);
    rc = gd_calibrations_propagate(cs, t, &err);
    CHECK(rc == 0);

    gd_rng *rng = gd_rng_create(42, 1, &err);
    double ages[5];
    rc = gd_cal_sample_prior(rng, cs, t, ages, &err);
    CHECK(rc == 0);

    /* Parent > child on every edge. */
    CHECK(gd_cal_check_order(t, ages) == 1);
    /* Root age in [2000, 5000]. */
    CHECK(ages[t->root] >= 2000 && ages[t->root] <= 5000);

    gd_rng_destroy(rng);
    gd_calset_destroy(cs);
    gd_tree_destroy(t);
}

TEST(test_cal_logprior_normal)
{
    gd_calibration cal = {NULL, NULL, GD_CAL_NORMAL_NODE, 3000, 500, 0, 0, NULL, 0, 0};
    /* At the mean, logprior = -log(500 * sqrt(2pi)). */
    double expected = -log(500 * sqrt(2 * 3.14159265358979));
    CHECK_APPROX(gd_cal_logprior(&cal, 3000), expected, 1e-6);
    /* Outside support: still finite (normal has infinite support). */
    CHECK(gd_is_finite(gd_cal_logprior(&cal, 10000)));
}

TEST(test_cal_logprior_uniform)
{
    gd_calibration cal = {NULL, NULL, GD_CAL_UNIFORM_NODE, 2000, 5000, 0, 0, NULL, 0, 0};
    /* Inside: -log(3000). */
    CHECK_APPROX(gd_cal_logprior(&cal, 3000), -log(3000), 1e-12);
    /* Outside: -inf. */
    CHECK(gd_cal_logprior(&cal, 1000) == -HUGE_VAL);
    CHECK(gd_cal_logprior(&cal, 6000) == -HUGE_VAL);
}

TEST(test_cal_logprior_fixed_tip)
{
    gd_calibration cal = {NULL, NULL, GD_CAL_FIXED_TIP, 1900, 0, 0, 1, NULL, 0, 0};
    CHECK(gd_cal_logprior(&cal, 1900) == 0.0);
    CHECK(gd_cal_logprior(&cal, 1901) == -HUGE_VAL);
}

TEST(test_cal_mrca_extraction)
{
    /* For date-pair-context: the MRCA of two leaves gives the pair's
     * divergence time. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_romance_tree(&err);
    CHECK(t != NULL);
    size_t it = gd_tree_find_label(t, "italian");
    size_t sp = gd_tree_find_label(t, "spanish");
    size_t mrca = gd_tree_mrca(t, it, sp);
    size_t rw = gd_tree_find_label(t, "romance_west");
    CHECK(mrca == rw);
    /* The pair's divergence = age of romance_west. */
    gd_tree_destroy(t);
}

static void run_tests(void)
{
    RUN(test_cal_type_from_string);
    RUN(test_cal_resolve_targets);
    RUN(test_cal_resolve_missing_target);
    RUN(test_cal_propagate_feasible);
    RUN(test_cal_contradiction_detected);
    RUN(test_cal_reject_invalid_normal_sd);
    RUN(test_cal_reject_inverted_uniform);
    RUN(test_cal_sample_prior_order);
    RUN(test_cal_logprior_normal);
    RUN(test_cal_logprior_uniform);
    RUN(test_cal_logprior_fixed_tip);
    RUN(test_cal_mrca_extraction);
}

TEST_MAIN()

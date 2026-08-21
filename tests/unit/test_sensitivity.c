/*
 * test_sensitivity.c — scoring a chronology re-dated without one calibration.
 *
 * Until D6 this was three static functions in cli/commands.c: scientific
 * policy -- a 100x divergence ceiling, a shift measure, and a three-way
 * verdict -- that no test could reach. The only way to exercise it was to run
 * the binary and read a TSV.
 *
 * Every case destroys its tree BEFORE asserting. CHECK returns from the test
 * on failure, so a tree held across an assertion leaks and LeakSanitizer's
 * report buries the FAIL line that says what broke. That is not hypothetical:
 * it is exactly what the first draft did when the shift rule was broken on
 * purpose to confirm these tests can fail.
 */
#include "test_runner.h"
#include "gd_analysis.h"
#include "gd_tree.h"
#include "gd_newick.h"
#include "gd_tree_build.h"   /* gd_tree_destroy: tests build trees directly */

#include <math.h>
#include <string.h>

/* A four-node tree: (a,b)inner, with `inner` the calibrated target. */
static gd_tree *small_tree(void)
{
    static const char NWK[] = "(a:1,b:1)inner;";
    gd_tree *t = NULL;
    gd_error err = GD_ERROR_INIT;
    if (gd_newick_parse(NWK, sizeof NWK - 1, &t, &err) != GD_E_OK) return NULL;
    return t;
}

/* A result carrying the given ages, over the given tree. Nothing is owned:
 * the caller keeps the arrays alive. */
static void make_result(gd_analysis_tree_result *res, gd_tree *tree,
                        gd_analysis_node_age *ages, size_t nages)
{
    memset(res, 0, sizeof *res);
    res->tree = tree;
    res->ages = ages;
    res->nages = nages;
    res->n_nodes = tree->nnodes;
}

static void age(gd_analysis_node_age *a, size_t node, const char *label,
                double median, double lo, double hi)
{
    memset(a, 0, sizeof *a);
    a->node = node;
    a->label = label;
    a->age.median = median;
    a->age.mean = median;
    a->age.hpd_lo = lo;
    a->age.hpd_hi = hi;
}

static void calib(gd_analysis_calibration *c, const char *target, double point)
{
    memset(c, 0, sizeof *c);
    c->target = (char *)target;
    c->id = (char *)"cal1";
    c->type_name = (char *)"normal_node";
    c->point = point;
    c->p1 = point;
}

TEST(test_sensitivity_ceiling_is_a_hundred_times_the_scale)
{
    gd_analysis_node_age ages[1];
    age(&ages[0], 0, "inner", 2000.0, 1800.0, 2200.0);
    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    res.ages = ages; res.nages = 1;

    gd_analysis_calibration c;
    calib(&c, "inner", 1500.0);
    /* The largest of the baseline medians and the calibration bounds. */
    CHECK_APPROX(gd_analysis_divergence_ceiling(&res, &c, 1), 200000.0, 1e-9);

    /* A calibration older than anything sampled sets the scale instead. */
    c.p1 = 9000.0;
    CHECK_APPROX(gd_analysis_divergence_ceiling(&res, &c, 1), 900000.0, 1e-9);

    /* Nothing to go on: a stated fallback, not a zero ceiling that would call
     * every age unbounded. */
    CHECK_APPROX(gd_analysis_divergence_ceiling(NULL, NULL, 0), 100000.0, 1e-9);
}

TEST(test_sensitivity_scores_a_recovered_calibration)
{
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);
    size_t inner = gd_tree_find_label(tree, "inner");

    gd_analysis_node_age b[1], v[1];
    age(&b[0], inner, "inner", 2000.0, 1900.0, 2100.0);
    age(&v[0], inner, "inner", 1950.0, 1700.0, 2200.0);
    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 1);
    make_result(&var, tree, v, 1);

    gd_analysis_calibration c;
    calib(&c, "inner", 2050.0);        /* inside [1700, 2200] */

    gd_sensitivity_row row;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &row, &err);
    gd_tree_destroy(tree);

    CHECK(rc == GD_E_OK);
    CHECK(row.outcome == GD_SENS_SCORED);
    CHECK(row.covered);
    CHECK_APPROX(row.median, 1950.0, 1e-9);
    CHECK_APPROX(row.abs_err, 100.0, 1e-9);
}

TEST(test_sensitivity_reports_a_missed_calibration)
{
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);
    size_t inner = gd_tree_find_label(tree, "inner");

    gd_analysis_node_age b[1], v[1];
    age(&b[0], inner, "inner", 2000.0, 1900.0, 2100.0);
    age(&v[0], inner, "inner", 1200.0, 1100.0, 1300.0);
    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 1);
    make_result(&var, tree, v, 1);

    gd_analysis_calibration c;
    calib(&c, "inner", 2050.0);        /* outside [1100, 1300] */

    gd_sensitivity_row row;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &row, &err);
    gd_tree_destroy(tree);

    CHECK(rc == GD_E_OK);
    CHECK(row.outcome == GD_SENS_SCORED);
    CHECK(!row.covered);
    CHECK_APPROX(row.abs_err, 850.0, 1e-9);
}

TEST(test_sensitivity_ceiling_fires_and_does_not_fire)
{
    /* The load-bearing case: with its only calibration withheld the age is
     * unidentified and the sampler wanders. Past the ceiling that is reported
     * as unbounded rather than as a spurious date. */
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);
    size_t inner = gd_tree_find_label(tree, "inner");

    gd_analysis_node_age b[1], v[1];
    age(&b[0], inner, "inner", 2000.0, 1900.0, 2100.0);
    age(&v[0], inner, "inner", 500000.0, 1000.0, 900000.0);
    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 1);
    make_result(&var, tree, v, 1);

    gd_analysis_calibration c;
    calib(&c, "inner", 2050.0);
    gd_sensitivity_row above, below, at;
    gd_error err = GD_ERROR_INIT;
    int rc_a = gd_analysis_sensitivity_score(&base, &var, &c, 600000.0, &above, &err);
    int rc_b = gd_analysis_sensitivity_score(&base, &var, &c, 400000.0, &below, &err);
    int rc_c = gd_analysis_sensitivity_score(&base, &var, &c, 500000.0, &at, &err);
    gd_tree_destroy(tree);

    /* Ceiling above the wandered age: still a measurement. */
    CHECK(rc_a == GD_E_OK);
    CHECK(above.outcome == GD_SENS_SCORED);
    /* Ceiling below it: unbounded. */
    CHECK(rc_b == GD_E_OK);
    CHECK(below.outcome == GD_SENS_UNBOUNDED);
    CHECK(isnan(below.median));
    CHECK(!below.covered);
    /* Exactly at the ceiling is still scored: the bound is inclusive. */
    CHECK(rc_c == GD_E_OK);
    CHECK(at.outcome == GD_SENS_SCORED);
}

TEST(test_sensitivity_absent_target_is_not_scored)
{
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);

    gd_analysis_node_age b[1], v[1];
    age(&b[0], 0, "a", 0.0, 0.0, 0.0);
    age(&v[0], 0, "a", 0.0, 0.0, 0.0);
    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 1);
    make_result(&var, tree, v, 1);

    gd_analysis_calibration c;
    calib(&c, "no_such_node", 2050.0);

    gd_sensitivity_row row;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &row, &err);
    gd_tree_destroy(tree);

    CHECK(rc == GD_E_OK);
    CHECK(row.outcome == GD_SENS_ABSENT);
    CHECK(isnan(row.median));
}

TEST(test_sensitivity_shift_skips_the_target_and_the_unbounded)
{
    /* max_other_shift answers "how far did the REST of the chronology move",
     * so it excludes the withheld node itself, and excludes nodes that ran
     * past the ceiling -- the distance an unidentified age wandered is not a
     * measurement either. */
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);
    size_t inner = gd_tree_find_label(tree, "inner");
    size_t leaf_a = gd_tree_find_label(tree, "a");

    gd_analysis_node_age b[2], v[2];
    age(&b[0], inner, "inner", 2000.0, 1900.0, 2100.0);
    age(&b[1], leaf_a, "a",     500.0,  450.0,  550.0);
    /* The target moves a long way; the other node moves 30. */
    age(&v[0], inner, "inner", 9000.0, 1000.0, 12000.0);
    age(&v[1], leaf_a, "a",     530.0,  480.0,  580.0);

    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 2);
    make_result(&var, tree, v, 2);

    gd_analysis_calibration c;
    calib(&c, "inner", 2050.0);
    gd_sensitivity_row near, past;
    gd_error err = GD_ERROR_INIT;
    int rc_n = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &near, &err);
    /* Now push the other node past the ceiling: it stops counting, and the
     * shift falls back to zero rather than reporting where it wandered. */
    v[1].age.median = 900000.0;
    int rc_p = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &past, &err);
    gd_tree_destroy(tree);

    CHECK(rc_n == GD_E_OK);
    /* 30, not 7000: the withheld node is skipped. */
    CHECK_APPROX(near.max_other_shift, 30.0, 1e-9);
    CHECK(rc_p == GD_E_OK);
    CHECK_APPROX(past.max_other_shift, 0.0, 1e-9);
}

TEST(test_sensitivity_refuses_mismatched_results)
{
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);
    gd_analysis_node_age b[1];
    age(&b[0], 0, "a", 1.0, 0.0, 2.0);
    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 1);
    make_result(&var, tree, b, 1);
    var.n_nodes = base.n_nodes + 1;      /* a different topology */

    gd_analysis_calibration c;
    calib(&c, "a", 1.0);
    gd_sensitivity_row row;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_analysis_sensitivity_score(&base, &var, &c, 1000.0, &row, &err);
    gd_tree_destroy(tree);

    CHECK(rc == GD_E_INVALID_ARG);
    CHECK(strstr(err.message, "tree changed") != NULL);
}

TEST(test_sensitivity_rejects_null_arguments)
{
    gd_sensitivity_row row;
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    gd_analysis_calibration c;
    calib(&c, "x", 1.0);
    CHECK(gd_analysis_sensitivity_score(NULL, &res, &c, 1.0, &row, &err) == GD_E_INVALID_ARG);
    CHECK(gd_analysis_sensitivity_score(&res, NULL, &c, 1.0, &row, &err) == GD_E_INVALID_ARG);
    CHECK(gd_analysis_sensitivity_score(&res, &res, NULL, 1.0, &row, &err) == GD_E_INVALID_ARG);
    CHECK(gd_analysis_sensitivity_score(&res, &res, &c, 1.0, NULL, &err) == GD_E_INVALID_ARG);
    /* A NULL error is always allowed. */
    CHECK(gd_analysis_sensitivity_score(NULL, &res, &c, 1.0, &row, NULL) == GD_E_INVALID_ARG);
}

TEST(test_sensitivity_coverage_bounds_are_inclusive)
{
    /* A withheld age sitting exactly on an HPD endpoint counts as recovered.
     * Nothing tested this until flipping >= to > in the coverage test left
     * every case still passing. */
    gd_tree *tree = small_tree();
    CHECK(tree != NULL);
    size_t inner = gd_tree_find_label(tree, "inner");

    gd_analysis_node_age b[1], v[1];
    age(&b[0], inner, "inner", 2000.0, 1900.0, 2100.0);
    age(&v[0], inner, "inner", 1950.0, 1700.0, 2200.0);
    gd_analysis_tree_result base, var;
    make_result(&base, tree, b, 1);
    make_result(&var, tree, v, 1);

    gd_sensitivity_row lo, hi, just_out;
    gd_error err = GD_ERROR_INIT;
    gd_analysis_calibration c;

    calib(&c, "inner", 1700.0);        /* exactly hpd_lo */
    int rc_lo = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &lo, &err);
    calib(&c, "inner", 2200.0);        /* exactly hpd_hi */
    int rc_hi = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &hi, &err);
    calib(&c, "inner", 2200.5);        /* just outside */
    int rc_out = gd_analysis_sensitivity_score(&base, &var, &c, 200000.0, &just_out, &err);
    gd_tree_destroy(tree);

    CHECK(rc_lo == GD_E_OK);
    CHECK(lo.covered);
    CHECK(rc_hi == GD_E_OK);
    CHECK(hi.covered);
    CHECK(rc_out == GD_E_OK);
    CHECK(!just_out.covered);
}

static void run_tests(void)
{
    RUN(test_sensitivity_ceiling_is_a_hundred_times_the_scale);
    RUN(test_sensitivity_scores_a_recovered_calibration);
    RUN(test_sensitivity_reports_a_missed_calibration);
    RUN(test_sensitivity_coverage_bounds_are_inclusive);
    RUN(test_sensitivity_ceiling_fires_and_does_not_fire);
    RUN(test_sensitivity_absent_target_is_not_scored);
    RUN(test_sensitivity_shift_skips_the_target_and_the_unbounded);
    RUN(test_sensitivity_refuses_mismatched_results);
    RUN(test_sensitivity_rejects_null_arguments);
}

TEST_MAIN()

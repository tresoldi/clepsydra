/*
 * test_tree_prior.c — the uniform-conditional volume and the lexical root.
 *
 * The volume V(x) of the ordered region below a root at age x is checked
 * against closed forms (all tips modern) and against Monte Carlo (dated tips,
 * where V is a genuine polynomial with mixed-sign coefficients).
 */
#include "test_runner.h"
#include "gd_tree_prior.h"
#include "gd_newick.h"
#include "gd_tree_build.h"
#include "gd_tree.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static gd_tree *parse(const char *nwk)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    return t;
}

/* Free = internal; fixed = tips at the given ages (0 unless named). */
static void masks(const gd_tree *t, bool *is_free, double *age,
                  const char *dated, double dated_age)
{
    for (size_t v = 0; v < t->nnodes; v++) {
        is_free[v] = !t->nodes[v].is_leaf;
        age[v] = 0.0;
        if (dated != NULL && t->nodes[v].label != NULL &&
            strcmp(t->nodes[v].label, dated) == 0) age[v] = dated_age;
    }
}

/* All tips modern: V(x) = (#ranked histories / m!) x^m, m free non-root
 * nodes. The balanced quartet has 2 rankings of its two cherries: V = x^2. */
TEST(test_tree_prior_volume_modern_tips)
{
    gd_tree *t = parse("((a,b)r1,(c,d)r2)root;");
    CHECK(t != NULL);
    bool fr[7]; double ag[7];
    masks(t, fr, ag, NULL, 0.0);
    gd_error err = GD_ERROR_INIT;
    gd_tree_prior *tp = gd_tree_prior_create(t, fr, ag, 0.0, 1.0, &err);
    CHECK(tp != NULL);
    double d1 = gd_tree_prior_log_volume(tp, 3000.0) - gd_tree_prior_log_volume(tp, 1500.0);
    gd_tree_prior_destroy(tp);
    gd_tree_destroy(t);
    CHECK_APPROX(d1, 2.0 * log(2.0), 1e-12);

    /* Caterpillar (((a,b)x,c)y,d)root: x < y < root, volume x^2 / 2. */
    t = parse("(((a,b)x,c)y,d)root;");
    masks(t, fr, ag, NULL, 0.0);
    tp = gd_tree_prior_create(t, fr, ag, 0.0, 1.0, &err);
    CHECK(tp != NULL);
    double d2 = gd_tree_prior_log_volume(tp, 900.0) - gd_tree_prior_log_volume(tp, 300.0);
    gd_tree_prior_destroy(tp);
    gd_tree_destroy(t);
    CHECK_APPROX(d2, 2.0 * log(3.0), 1e-12);
}

/* A dated tip makes V a polynomial with a floor: checked by Monte Carlo. */
TEST(test_tree_prior_volume_dated_tip_vs_monte_carlo)
{
    gd_tree *t = parse("(((a,b)x,c)y,(d,e)w)root;");
    CHECK(t != NULL);
    size_t n = t->nnodes;
    bool *fr = (bool *)malloc(n * sizeof(bool));
    double *ag = (double *)malloc(n * sizeof(double));
    masks(t, fr, ag, "a", 600.0);           /* a dated tip under x */
    gd_error err = GD_ERROR_INIT;
    gd_tree_prior *tp = gd_tree_prior_create(t, fr, ag, 0.0, 1.0, &err);
    CHECK(tp != NULL);
    size_t x = gd_tree_find_label(t, "x"), y = gd_tree_find_label(t, "y");
    size_t w = gd_tree_find_label(t, "w");
    gd_rng *rng = gd_rng_create(9u, 1u, &err);
    const double R[2] = { 1000.0, 2500.0 };
    double vol[2];
    for (int k = 0; k < 2; k++) {
        size_t hit = 0, N = 2000000;
        for (size_t i = 0; i < N; i++) {
            double ax = R[k] * gd_rng_uniform_closed_open(rng);
            double ay = R[k] * gd_rng_uniform_closed_open(rng);
            double aw = R[k] * gd_rng_uniform_closed_open(rng);
            if (ax > 600.0 && ay > ax) hit++;
            (void)aw;
        }
        vol[k] = (double)hit / (double)N * R[k] * R[k] * R[k];
    }
    double exact = gd_tree_prior_log_volume(tp, R[1]) - gd_tree_prior_log_volume(tp, R[0]);
    double mc = log(vol[1]) - log(vol[0]);
    /* Closed form for this tree: (R - 600)^2 / 2 * R. */
    double closed = log(0.5 * (R[1] - 600) * (R[1] - 600) * R[1])
                  - log(0.5 * (R[0] - 600) * (R[0] - 600) * R[0]);
    CHECK(isinf(gd_tree_prior_log_volume(tp, 599.0)));
    gd_rng_destroy(rng);
    gd_tree_prior_destroy(tp);
    gd_tree_destroy(t);
    free(fr); free(ag);
    (void)x; (void)y; (void)w;
    CHECK_APPROX(exact, closed, 1e-10);
    CHECK_APPROX(exact, mc, 0.01);
}

/* The lexical-root factor: with no volume (root fixed), the log prior as a
 * function of log rate is the normal density of log(rho * a) on its log scale. */
TEST(test_tree_prior_lexical_root)
{
    gd_tree *t = parse("(a,b)root;");
    bool fr[3] = { false, false, false };
    double ag[3] = { 0, 0, 0 };
    gd_error err = GD_ERROR_INIT;
    gd_tree_prior *tp = gd_tree_prior_create(t, fr, ag, log(0.3), 1.5, &err);
    CHECK(tp != NULL);
    double a = 2000.0;
    double lr1 = log(0.3 / a), lr2 = log(0.3 / a) + 1.5;
    double d = gd_tree_prior_logpdf(tp, a, lr2) - gd_tree_prior_logpdf(tp, a, lr1);
    gd_tree_prior_destroy(tp);
    gd_tree_destroy(t);
    CHECK_APPROX(d, -0.5, 1e-12);
    CHECK(gd_tree_prior_create(t = parse("(a,b)r;"), fr, ag, 0.0, 0.0, &err) == NULL);
    gd_tree_destroy(t);
}

static void run_tests(void)
{
    RUN(test_tree_prior_volume_modern_tips);
    RUN(test_tree_prior_volume_dated_tip_vs_monte_carlo);
    RUN(test_tree_prior_lexical_root);
}

TEST_MAIN()

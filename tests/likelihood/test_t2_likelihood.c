/*
 * test_t2_likelihood.c — validation of the binary cognate-set likelihood (T2)
 * for Clepsydra (M7).
 *
 * The M7 acceptance gate (PLAN.md) is brute-force agreement, simulation
 * recovery, a documented within-meaning dependence warning, and no
 * cross-encoding comparison without a shared observable space.
 *
 * The brute force here is a full enumeration over all 2^n node-state
 * assignments — it shares nothing with Felsenstein pruning, so agreement
 * between the two is a genuine independent check (product brief §3.4,
 * docs/validation.md §3). It recomputes the CTMC transition matrix from the
 * closed form directly, independently of gd_t2_transition_logs.
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_model_binary.h"
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_rng.h"
#include "gd_math.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TOL 1e-9

/* ------------------------------------------------------------ brute force */

/* Closed-form 2x2 transition matrix, computed linearly and independently. */
static void bf_P(double lg, double ll, double t, double P[2][2])
{
    double mu = lg + ll;
    double e = exp(-mu * t);
    P[0][1] = (lg / mu) * (1.0 - e);   /* 0 -> 1 (gain) */
    P[0][0] = 1.0 - P[0][1];
    P[1][0] = (ll / mu) * (1.0 - e);   /* 1 -> 0 (loss) */
    P[1][1] = 1.0 - P[1][0];
}

static void bf_root_prior(const gd_t2_options *opt, double lg, double ll,
                          double *pi0, double *pi1)
{
    if (opt->root_prior == GD_T2_ROOT_FIXED) {
        *pi1 = opt->root_pi1;
        *pi0 = 1.0 - opt->root_pi1;
    } else {
        double mu = lg + ll;
        *pi0 = ll / mu;
        *pi1 = lg / mu;
    }
}

/* Raw (uncorrected) log-likelihood of an observation vector, by enumerating
 * every assignment of {0,1} to all nodes and keeping those consistent with the
 * data. obs[v]: 0 absent, >0 present, <0 missing (free). */
static double bf_raw(const gd_tree *t, double lg, double ll,
                     const int *obs, const gd_t2_options *opt)
{
    size_t n = t->nnodes;
    double pi0, pi1;
    bf_root_prior(opt, lg, ll, &pi0, &pi1);
    /* calloc, and checked: the brute force is the independent oracle the
     * pruning is validated against, so it must not read uninitialised state or
     * walk off a failed allocation. */
    int *g = (int *)calloc(n, sizeof(int));
    if (g == NULL) return NAN;
    double total = -HUGE_VAL;
    size_t npat = (size_t)1 << n;
    for (size_t mask = 0; mask < npat; mask++) {
        int ok = 1;
        for (size_t v = 0; v < n; v++) {
            g[v] = (int)((mask >> v) & 1u);
            if (obs[v] == 0 && g[v] != 0) { ok = 0; break; }
            if (obs[v] > 0 && g[v] != 1) { ok = 0; break; }
        }
        if (!ok) continue;
        double lp = log((g[t->root] == 1) ? pi1 : pi0);
        for (size_t v = 0; v < n; v++) {
            if (v == t->root) continue;
            size_t u = t->nodes[v].parent;
            double dt = t->nodes[u].age - t->nodes[v].age;
            if (dt < 0) dt = 0;
            double P[2][2];
            bf_P(lg, ll, dt, P);
            lp += log(P[g[u]][g[v]]);
        }
        total = gd_logaddexp(total, lp);
    }
    free(g);
    return total;
}

/* Ascertainment-corrected brute-force log-likelihood. */
static double bf_class(const gd_tree *t, double lg, double ll,
                       const int *obs, const gd_t2_options *opt)
{
    double raw = bf_raw(t, lg, ll, obs, opt);
    if (!opt->ascertainment) return raw;
    size_t n = t->nnodes;
    int *allabs = (int *)malloc(n * sizeof(int));
    for (size_t i = 0; i < n; i++) allabs[i] = (obs[i] < 0) ? -1 : 0;
    double la = bf_raw(t, lg, ll, allabs, opt);
    free(allabs);
    double denom = gd_log1mexp(la);
    if (!isfinite(denom)) return raw;   /* mirror the library's fallback */
    return raw - denom;
}

/* Library per-class log-likelihood. */
static double lib_class(gd_t2_ws *ws, double lg, double ll, const int *obs,
                        size_t n, const gd_t2_options *opt, int *rc_out)
{
    gd_binary_class cls = { obs, n };
    gd_t2_class_result res;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_t2_loglik_class(ws, lg, ll, &cls, opt, &res, &err);
    if (rc_out) *rc_out = rc;
    if (rc != 0) return -HUGE_VAL;
    return res.log_lik;
}

static int agrees(gd_tree *t, double lg, double ll, const int *obs,
                  const gd_t2_options *opt)
{
    gd_error err = GD_ERROR_INIT;
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    if (ws == NULL) return 0;
    int rc = 0;
    double lib = lib_class(ws, lg, ll, obs, t->nnodes, opt, &rc);
    gd_t2_ws_destroy(ws);
    if (rc != 0) return 0;
    double bf = bf_class(t, lg, ll, obs, opt);
    int ok = (fabs(lib - bf) <= TOL);
    if (!ok) {
        printf("    lg=%g ll=%g lib=%.12g bf=%.12g asc=%d\n",
               lg, ll, lib, bf, opt->ascertainment);
        printf("    obs =");
        for (size_t i = 0; i < t->nnodes; i++) printf(" %d", obs[i]);
        printf("\n");
    }
    return ok;
}

/* ------------------------------------------------------------- fixtures */

static void set_age(gd_tree *t, const char *label, double age)
{
    size_t i = gd_tree_find_label(t, label);
    if (i != GD_NODE_NONE) t->nodes[i].age = age;
}

/* ((a,b)r1,c)root — 5 nodes, binary. */
static gd_tree *tree_simple(void)
{
    const char *nwk = "((a,b)r1,c)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    set_age(t, "root", 3000);
    set_age(t, "r1", 2000);
    return t;
}

/* ((a,b)r1,(c,d)r2)root — 7 nodes, binary. */
static gd_tree *tree_quartet(void)
{
    const char *nwk = "((a,b)r1,(c,d)r2)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    set_age(t, "root", 4000);
    set_age(t, "r1", 2500);
    set_age(t, "r2", 1500);
    return t;
}

/* ((a,b,c)r1,d)root — a polytomy below the root, ancient tip b. */
static gd_tree *tree_polytomy(void)
{
    const char *nwk = "((a,b,c)r1,d)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    set_age(t, "root", 5000);
    set_age(t, "r1", 3000);
    set_age(t, "b", 800);
    return t;
}

/* Balanced binary tree of the given depth; node i has children 2i+1, 2i+2. */
static gd_tree *tree_balanced(size_t depth, gd_error *err)
{
    size_t nnodes = ((size_t)1 << (depth + 1)) - 1;
    gd_tree *t = gd_tree_create(nnodes, err);
    if (t == NULL) return NULL;
    t->root = 0;
    for (size_t i = 0; i < nnodes; i++) {
        size_t l = 2 * i + 1, r = 2 * i + 2;
        if (r < nnodes) {
            if (gd_tree_add_child(t, i, l, err) != 0) goto fail;
            if (gd_tree_add_child(t, i, r, err) != 0) goto fail;
        }
    }
    if (gd_tree_build_traversals(t, err) != 0) goto fail;
    for (size_t i = 0; i < nnodes; i++) {
        size_t d = 0, x = i;
        while (x != 0) { x = (x - 1) / 2; d++; }
        t->nodes[i].age = (double)(depth - d) * 1000.0;
    }
    return t;
fail:
    gd_tree_destroy(t);
    return NULL;
}

/* --------------------------------------------------- transition matrix */

TEST(test_t2_transition_matrix)
{
    gd_error err = GD_ERROR_INIT;
    double lg = 0.0007, ll = 0.0011;
    double logP[2][2];
    CHECK(gd_t2_transition_logs(lg, ll, 1200.0, logP, &err) == 0);
    /* Rows are probability distributions. */
    CHECK_APPROX(exp(logP[0][0]) + exp(logP[0][1]), 1.0, 1e-12);
    CHECK_APPROX(exp(logP[1][0]) + exp(logP[1][1]), 1.0, 1e-12);
    /* Against the independent closed form. */
    double P[2][2];
    bf_P(lg, ll, 1200.0, P);
    CHECK_APPROX(exp(logP[0][0]), P[0][0], 1e-12);
    CHECK_APPROX(exp(logP[0][1]), P[0][1], 1e-12);
    CHECK_APPROX(exp(logP[1][0]), P[1][0], 1e-12);
    CHECK_APPROX(exp(logP[1][1]), P[1][1], 1e-12);
    /* t = 0 is the identity. */
    CHECK(gd_t2_transition_logs(lg, ll, 0.0, logP, &err) == 0);
    CHECK_APPROX(exp(logP[0][0]), 1.0, 1e-12);
    CHECK(logP[0][1] == -HUGE_VAL);
    CHECK(logP[1][0] == -HUGE_VAL);
    CHECK_APPROX(exp(logP[1][1]), 1.0, 1e-12);
    /* t -> infinity approaches the stationary distribution on every row. */
    CHECK(gd_t2_transition_logs(lg, ll, 1e9, logP, &err) == 0);
    double mu = lg + ll;
    CHECK_APPROX(exp(logP[0][1]), lg / mu, 1e-9);
    CHECK_APPROX(exp(logP[1][1]), lg / mu, 1e-9);
    /* Bad arguments. */
    CHECK(gd_t2_transition_logs(0.0, ll, 10.0, logP, &err) == GD_E_RANGE);
    CHECK(gd_t2_transition_logs(lg, ll, -1.0, logP, &err) == GD_E_RANGE);
}

/* ------------------------------------------------- brute-force agreement */

TEST(test_t2_bf_binary_all_patterns)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    size_t tip[4];
    tip[0] = gd_tree_find_label(t, "a");
    tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c");
    tip[3] = gd_tree_find_label(t, "d");
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    /* Every presence/absence pattern of the four tips, both with and without
     * the ascertainment correction, at two rate settings. */
    for (unsigned mask = 0; mask < 16u; mask++) {
        int obs[7];
        for (size_t i = 0; i < 7; i++) obs[i] = -1;
        for (size_t i = 0; i < 4; i++) obs[tip[i]] = (int)((mask >> i) & 1u);
        opt.ascertainment = true;
        CHECK(agrees(t, 0.0006, 0.0009, obs, &opt));
        CHECK(agrees(t, 0.002, 0.0005, obs, &opt));
        opt.ascertainment = false;
        CHECK(agrees(t, 0.0006, 0.0009, obs, &opt));
    }
    gd_tree_destroy(t);
}

TEST(test_t2_bf_missing_and_internal)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    size_t a = gd_tree_find_label(t, "a"), b = gd_tree_find_label(t, "b");
    size_t c = gd_tree_find_label(t, "c"), d = gd_tree_find_label(t, "d");
    size_t r1 = gd_tree_find_label(t, "r1"), r2 = gd_tree_find_label(t, "r2");

    /* A missing tip is marginalized. */
    int obs[7];
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[a] = 1; obs[b] = 0; obs[c] = 1;      /* d missing */
    CHECK(agrees(t, 0.0007, 0.0008, obs, &opt));

    /* An observed internal node (an ancient / attested language). */
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[a] = 1; obs[b] = 1; obs[c] = 0; obs[d] = 0;
    obs[r1] = 1;                              /* the a-b ancestor carried it */
    CHECK(agrees(t, 0.0007, 0.0008, obs, &opt));
    obs[r2] = 0;                              /* the c-d ancestor did not */
    CHECK(agrees(t, 0.0007, 0.0008, obs, &opt));
    gd_tree_destroy(t);
}

TEST(test_t2_bf_root_prior_fixed)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    opt.root_prior = GD_T2_ROOT_FIXED;
    opt.root_pi1 = 0.3;
    size_t a = gd_tree_find_label(t, "a"), b = gd_tree_find_label(t, "b");
    size_t c = gd_tree_find_label(t, "c");
    int obs[5];
    for (unsigned mask = 0; mask < 8u; mask++) {
        for (size_t i = 0; i < 5; i++) obs[i] = -1;
        obs[a] = (int)(mask & 1u);
        obs[b] = (int)((mask >> 1) & 1u);
        obs[c] = (int)((mask >> 2) & 1u);
        CHECK(agrees(t, 0.0008, 0.0006, obs, &opt));
    }
    gd_tree_destroy(t);
}

TEST(test_t2_bf_polytomy_native)
{
    /* Felsenstein pruning combines any number of children directly: the
     * polytomy is handled without any ladder, unlike T1. Checked against the
     * brute force on the same unresolved polytomy. */
    gd_tree *t = tree_polytomy();
    CHECK(t != NULL);
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    size_t idx[4];
    idx[0] = gd_tree_find_label(t, "a");
    idx[1] = gd_tree_find_label(t, "b");
    idx[2] = gd_tree_find_label(t, "c");
    idx[3] = gd_tree_find_label(t, "d");
    for (unsigned mask = 0; mask < 16u; mask++) {
        int obs[6];
        for (size_t i = 0; i < 6; i++) obs[i] = -1;
        for (size_t i = 0; i < 4; i++) obs[idx[i]] = (int)((mask >> i) & 1u);
        CHECK(agrees(t, 0.0005, 0.0007, obs, &opt));
    }
    gd_tree_destroy(t);
}

TEST(test_t2_bf_randomized)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(20260731u, 13u, &err);
    CHECK(rng != NULL);
    const char *shapes[] = {
        "((a,b)r1,(c,d)r2)root;",
        "(((a,b)r1,c)r2,d)root;",
        "((a,b,c)r1,d)root;",
        "(a,b,(c,d)r1)root;",
        "((a,b)r1,c,d)root;"
    };
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    int checked = 0;
    for (size_t s = 0; s < sizeof shapes / sizeof shapes[0]; s++) {
        for (int rep = 0; rep < 12; rep++) {
            gd_tree *t = NULL;
            CHECK(gd_newick_parse(shapes[s], strlen(shapes[s]), &t, &err) == 0);
            for (size_t i = 0; i < t->nnodes; i++) {
                size_t d = 0, x = i;
                while (t->nodes[x].parent != GD_NODE_NONE) {
                    x = t->nodes[x].parent; d++;
                }
                double base = 4000.0 - 900.0 * (double)d;
                t->nodes[i].age = t->nodes[i].nchildren == 0
                    ? gd_rng_uniform_closed_open(rng) * 400.0
                    : base + gd_rng_uniform_closed_open(rng) * 200.0;
            }
            for (size_t i = 0; i < t->nnodes; i++) {
                size_t v = t->postorder[i];
                size_t p = t->nodes[v].parent;
                if (p != GD_NODE_NONE && t->nodes[p].age <= t->nodes[v].age) {
                    t->nodes[p].age = t->nodes[v].age + 100.0;
                }
            }
            int *obs = (int *)malloc(t->nnodes * sizeof(int));
            for (size_t i = 0; i < t->nnodes; i++) {
                double u = gd_rng_uniform_closed_open(rng);
                obs[i] = (u < 0.25) ? -1 : (u < 0.6 ? 0 : 1);
            }
            double lg = 0.0002 + gd_rng_uniform_closed_open(rng) * 0.002;
            double ll = 0.0002 + gd_rng_uniform_closed_open(rng) * 0.002;
            opt.ascertainment = (rep % 2);
            CHECK(agrees(t, lg, ll, obs, &opt));
            checked++;
            free(obs);
            gd_tree_destroy(t);
        }
    }
    CHECK(checked == 60);
    gd_rng_destroy(rng);
}

/* --------------------------------------------------------- normalization */

/* With the ascertainment correction the class likelihood is a distribution
 * over the *observable* tip patterns (every pattern except all-absent), so it
 * sums to one over them — an absolute check needing no reference. */
TEST(test_t2_normalization)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    size_t tip[4];
    tip[0] = gd_tree_find_label(t, "a");
    tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c");
    tip[3] = gd_tree_find_label(t, "d");
    gd_t2_options opt;
    gd_t2_options_default(&opt);   /* ascertainment on */

    double lgs[2] = {0.0006, 0.0018};
    double lls[2] = {0.0009, 0.0004};
    for (int k = 0; k < 2; k++) {
        double mass = 0.0;
        for (unsigned mask = 1; mask < 16u; mask++) {   /* skip all-absent */
            int obs[7];
            for (size_t i = 0; i < 7; i++) obs[i] = -1;
            for (size_t i = 0; i < 4; i++) obs[tip[i]] = (int)((mask >> i) & 1u);
            gd_binary_class cls = { obs, 7 };
            gd_t2_class_result res;
            CHECK(gd_t2_loglik_class(ws, lgs[k], lls[k], &cls, &opt, &res, &err) == 0);
            mass += exp(res.log_lik);
        }
        CHECK_APPROX(mass, 1.0, 1e-9);
    }

    /* Without ascertainment the raw likelihood sums to one over *all* patterns
     * including all-absent. */
    opt.ascertainment = false;
    double mass_all = 0.0;
    for (unsigned mask = 0; mask < 16u; mask++) {
        int obs[7];
        for (size_t i = 0; i < 7; i++) obs[i] = -1;
        for (size_t i = 0; i < 4; i++) obs[tip[i]] = (int)((mask >> i) & 1u);
        gd_binary_class cls = { obs, 7 };
        gd_t2_class_result res;
        CHECK(gd_t2_loglik_class(ws, 0.0006, 0.0009, &cls, &opt, &res, &err) == 0);
        mass_all += exp(res.log_lik);
    }
    CHECK_APPROX(mass_all, 1.0, 1e-9);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t2_ascertainment_effect)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    size_t a = gd_tree_find_label(t, "a"), b = gd_tree_find_label(t, "b");
    size_t c = gd_tree_find_label(t, "c");
    int obs[5];
    for (size_t i = 0; i < 5; i++) obs[i] = -1;
    obs[a] = 1; obs[b] = 0; obs[c] = 0;
    gd_binary_class cls = { obs, 5 };

    gd_t2_options opt;
    gd_t2_options_default(&opt);
    gd_t2_class_result res;
    CHECK(gd_t2_loglik_class(ws, 0.0008, 0.0012, &cls, &opt, &res, &err) == 0);
    /* Conditioning on observability raises the likelihood: dividing by a
     * denominator < 1 adds a positive amount in log space. */
    CHECK(res.ascertained == 1);
    CHECK(res.log_lik > res.log_lik_raw);
    CHECK(res.log_p_all_absent < 0.0);
    CHECK_APPROX(res.log_lik,
                 res.log_lik_raw - gd_log1mexp(res.log_p_all_absent), 1e-12);
    CHECK(res.all_absent == 0);

    /* A class observed present nowhere is degenerate and flagged. */
    for (size_t i = 0; i < 5; i++) obs[i] = 0;
    CHECK(gd_t2_loglik_class(ws, 0.0008, 0.0012, &cls, &opt, &res, &err) == 0);
    CHECK(res.all_absent == 1);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

/* --------------------------------------------------------------- dataset */

TEST(test_t2_dataset_sums_and_scaling)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    size_t tip[4];
    tip[0] = gd_tree_find_label(t, "a");
    tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c");
    tip[3] = gd_tree_find_label(t, "d");
    /* Three cognate classes over two meanings: meaning 0 has two classes
     * (synonymy), meaning 1 has one. */
    int c0[7], c1[7], c2[7];
    for (size_t i = 0; i < 7; i++) { c0[i] = -1; c1[i] = -1; c2[i] = -1; }
    const int p0[4] = {1, 1, 0, 0};
    const int p1[4] = {0, 0, 1, 1};
    const int p2[4] = {1, 0, 1, 0};
    for (size_t i = 0; i < 4; i++) {
        c0[tip[i]] = p0[i]; c1[tip[i]] = p1[i]; c2[tip[i]] = p2[i];
    }
    gd_binary_class classes[3] = { {c0, 7}, {c1, 7}, {c2, 7} };
    size_t concept_ids[3] = { 0, 0, 1 };
    double rates[2] = { 1.0, 1.0 };
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    gd_t2_class_result per[3];
    gd_t2_dataset_result sum;
    double total = 0.0;
    CHECK(gd_t2_loglik_dataset(ws, 0.0007, 0.0009, classes, concept_ids, 3,
                               rates, 2, &opt, per, &sum, &total, &err) == 0);
    CHECK_APPROX(total, per[0].log_lik + per[1].log_lik + per[2].log_lik, 1e-12);
    CHECK(sum.n_classes == 3);
    CHECK(sum.n_concepts == 2);
    CHECK(sum.max_classes_per_concept == 2);
    CHECK(sum.within_meaning_dependence_warning == 1);   /* meaning 0 has two */
    CHECK(sum.cross_encoding_incomparable == 1);

    /* A concept rate scales that meaning's rates: with all-ones vs a doubled
     * meaning-0 rate, only the two classes of meaning 0 change. */
    double rates2[2] = { 2.0, 1.0 };
    gd_t2_class_result per2[3];
    double total2 = 0.0;
    CHECK(gd_t2_loglik_dataset(ws, 0.0007, 0.0009, classes, concept_ids, 3,
                               rates2, 2, &opt, per2, NULL, &total2, &err) == 0);
    CHECK(fabs(per2[0].log_lik - per[0].log_lik) > 1e-9);
    CHECK(fabs(per2[1].log_lik - per[1].log_lik) > 1e-9);
    CHECK_APPROX(per2[2].log_lik, per[2].log_lik, 1e-12);   /* meaning 1 fixed */

    /* Scaling meaning 0 by rate 2 equals passing 2x base rates for those two
     * classes directly. */
    gd_binary_class m0a = { c0, 7 };
    gd_t2_class_result direct;
    CHECK(gd_t2_loglik_class(ws, 0.0014, 0.0018, &m0a, &opt, &direct, &err) == 0);
    CHECK_APPROX(per2[0].log_lik, direct.log_lik, 1e-12);

    /* No within-meaning dependence when every meaning has one class. */
    size_t ids_single[3] = { 0, 1, 2 };
    CHECK(gd_t2_loglik_dataset(ws, 0.0007, 0.0009, classes, ids_single, 3,
                               NULL, 3, &opt, NULL, &sum, &total, &err) == 0);
    CHECK(sum.max_classes_per_concept == 1);
    CHECK(sum.within_meaning_dependence_warning == 0);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t2_partition_weights)
{
    gd_error err = GD_ERROR_INIT;
    /* Meaning 0: 3 classes; meaning 1: 1 class; meaning 2: 0 classes. */
    size_t ids[4] = { 0, 0, 0, 1 };
    double w[3];
    /* Meaning weighting (King 2026 default): every populated meaning counts 1;
     * a meaning with no cognate class contributes nothing. */
    CHECK(gd_t2_partition_weights(ids, 4, 3, GD_T2_WEIGHT_MEANING, w, &err) == 0);
    CHECK_APPROX(w[0], 1.0, 1e-15);
    CHECK_APPROX(w[1], 1.0, 1e-15);
    CHECK_APPROX(w[2], 0.0, 1e-15);
    /* Cognate-count weighting (opt-in legacy): the count of classes. */
    CHECK(gd_t2_partition_weights(ids, 4, 3, GD_T2_WEIGHT_COGNATE, w, &err) == 0);
    CHECK_APPROX(w[0], 3.0, 1e-15);
    CHECK_APPROX(w[1], 1.0, 1e-15);
    CHECK_APPROX(w[2], 0.0, 1e-15);
    /* Out-of-range concept id is refused. */
    size_t bad[1] = { 5 };
    CHECK(gd_t2_partition_weights(bad, 1, 3, GD_T2_WEIGHT_MEANING, w, &err)
          == GD_E_INVALID_ARG);
}

TEST(test_t2_weight_from_string)
{
    gd_error err = GD_ERROR_INIT;
    gd_t2_partition_weight m;
    CHECK(gd_t2_weight_from_string("meaning", &m, &err) == 0);
    CHECK(m == GD_T2_WEIGHT_MEANING);
    CHECK(gd_t2_weight_from_string("cognate", &m, &err) == 0);
    CHECK(m == GD_T2_WEIGHT_COGNATE);
    CHECK(gd_t2_weight_from_string(NULL, &m, &err) == 0);
    CHECK(m == GD_T2_WEIGHT_MEANING);          /* the default */
    CHECK(gd_t2_weight_from_string("classes", &m, &err) == GD_E_INVALID_ARG);
}

/* ------------------------------------------------------------- simulation */

TEST(test_t2_simulate_matches_likelihood)
{
    /* Monte-Carlo frequencies of each observable tip pattern must match the
     * ascertainment-corrected likelihood of that pattern (all-absent draws are
     * discarded, which is exactly the ascertainment condition). */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    gd_rng *rng = gd_rng_create(31u, 5u, &err);
    CHECK(rng != NULL);
    size_t tip[3];
    tip[0] = gd_tree_find_label(t, "a");
    tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c");
    const double lg = 0.0009, ll = 0.0006;
    gd_t2_options opt;
    gd_t2_options_default(&opt);

    const size_t NREP = 120000;
    size_t counts[8];
    for (size_t i = 0; i < 8; i++) counts[i] = 0;
    size_t observable = 0;
    int sim[5];
    for (size_t rep = 0; rep < NREP; rep++) {
        CHECK(gd_t2_simulate_class(rng, t, lg, ll, &opt, sim, &err) == 0);
        unsigned pat = 0;
        for (size_t i = 0; i < 3; i++) pat |= (unsigned)(sim[tip[i]] << i);
        if (pat == 0) continue;   /* all-absent: unobservable */
        counts[pat]++;
        observable++;
    }
    for (unsigned pat = 1; pat < 8u; pat++) {
        int obs[5];
        for (size_t i = 0; i < 5; i++) obs[i] = -1;
        for (size_t i = 0; i < 3; i++) obs[tip[i]] = (int)((pat >> i) & 1u);
        gd_binary_class cls = { obs, 5 };
        gd_t2_class_result res;
        CHECK(gd_t2_loglik_class(ws, lg, ll, &cls, &opt, &res, &err) == 0);
        double p = exp(res.log_lik);
        double freq = (double)counts[pat] / (double)observable;
        double se = sqrt(p * (1.0 - p) / (double)observable);
        CHECK(fabs(freq - p) < 4.0 * se + 1e-12);
    }
    gd_rng_destroy(rng);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t2_simulation_recovers_rate)
{
    /* Simulate many observable classes at a known rate scale, refit the scale
     * by profile likelihood, check the peak is near the truth. Gain and loss
     * share the scale (their ratio is fixed), mirroring a lexical clock. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_balanced(3, &err);      /* 15 nodes, 8 tips */
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    gd_rng *rng = gd_rng_create(2026u, 17u, &err);
    CHECK(rng != NULL);

    const double lg = 0.0006, ll = 0.0009;
    const size_t NCLASS = 400;
    size_t n = t->nnodes;
    gd_t2_options opt;
    gd_t2_options_default(&opt);

    int *store = (int *)malloc(NCLASS * n * sizeof(int));
    gd_binary_class *classes =
        (gd_binary_class *)malloc(NCLASS * sizeof(gd_binary_class));
    int *sim = (int *)malloc(n * sizeof(int));
    size_t got = 0;
    while (got < NCLASS) {
        CHECK(gd_t2_simulate_class(rng, t, lg, ll, &opt, sim, &err) == 0);
        int any = 0;
        for (size_t i = 0; i < n; i++) {
            if (t->nodes[i].is_leaf && sim[i] == 1) { any = 1; break; }
        }
        if (!any) continue;                    /* ascertainment: skip all-absent */
        int *row = store + got * n;
        for (size_t i = 0; i < n; i++) {
            row[i] = t->nodes[i].is_leaf ? sim[i] : -1;   /* tips only */
        }
        classes[got].present = row;
        classes[got].nnodes = n;
        got++;
    }

    size_t ids_zero[1] = { 0 };   /* one meaning; rate scale applied per class */
    (void)ids_zero;
    double best_ll = -HUGE_VAL, best_s = 0.0;
    for (int k = 0; k < 41; k++) {
        double s = pow(10.0, -0.9 + 1.8 * (double)k / 40.0);  /* ~0.13 .. ~7.9 */
        double ll_tot = 0.0;
        for (size_t c = 0; c < NCLASS; c++) {
            gd_t2_class_result res;
            CHECK(gd_t2_loglik_class(ws, lg * s, ll * s, &classes[c], &opt,
                                     &res, &err) == 0);
            ll_tot += res.log_lik;
        }
        if (ll_tot > best_ll) { best_ll = ll_tot; best_s = s; }
    }
    CHECK(best_s > 0.6);
    CHECK(best_s < 1.7);

    free(sim); free(classes); free(store);
    gd_rng_destroy(rng);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

/* --------------------------------------------------------------- misc */

TEST(test_t2_bad_arguments)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[5] = { -1, -1, 1, 0, 1 };
    gd_binary_class cls = { obs, 5 };
    gd_binary_class wrong = { obs, 4 };
    gd_t2_class_result res;
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    CHECK(gd_t2_loglik_class(ws, 0.0, 0.001, &cls, &opt, &res, &err) == GD_E_RANGE);
    CHECK(gd_t2_loglik_class(ws, 0.001, -1.0, &cls, &opt, &res, &err) == GD_E_RANGE);
    CHECK(gd_t2_loglik_class(ws, 0.001, 0.001, &wrong, &opt, &res, &err)
          == GD_E_INVALID_ARG);
    CHECK(gd_t2_loglik_class(NULL, 0.001, 0.001, &cls, &opt, &res, &err)
          == GD_E_INVALID_ARG);
    /* Fixed root prior out of (0,1). */
    opt.root_prior = GD_T2_ROOT_FIXED;
    opt.root_pi1 = 1.0;
    CHECK(gd_t2_loglik_class(ws, 0.001, 0.001, &cls, &opt, &res, &err) == GD_E_RANGE);
    /* Chronology inversion is reported, not silently clamped. */
    gd_t2_options_default(&opt);
    t->nodes[gd_tree_find_label(t, "r1")].age = 9000;   /* older than the root */
    CHECK(gd_t2_loglik_class(ws, 0.001, 0.001, &cls, &opt, &res, &err) == GD_E_RANGE);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t2_ages_are_reread)
{
    /* The workspace must follow node ages as an MCMC chain moves them. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[5] = { -1, -1, 1, 0, 1 };
    gd_binary_class cls = { obs, 5 };
    gd_t2_options opt;
    gd_t2_options_default(&opt);
    gd_t2_class_result res;
    CHECK(gd_t2_loglik_class(ws, 0.001, 0.001, &cls, &opt, &res, &err) == 0);
    double before = res.log_lik;
    t->nodes[gd_tree_find_label(t, "root")].age = 8000;
    CHECK(gd_t2_loglik_class(ws, 0.001, 0.001, &cls, &opt, &res, &err) == 0);
    CHECK(fabs(before - res.log_lik) > 1e-6);
    CHECK(agrees(t, 0.001, 0.001, obs, &opt));   /* still matches brute force */
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t2_mcmc_adapter_matches_the_likelihood)
{
    /*
     * T2 as a gd_mcmc_model. This is the seam's second adapter, and the reason
     * the seam is real rather than hypothetical: T2 has been validated against
     * full enumeration since M7 and could reach no sampler until E3, because
     * gd_mcmc_state carried one scalar log-rate and this likelihood needs two.
     *
     * The adapter must report exactly what gd_t2_loglik_dataset reports for
     * the same rates -- it is a view onto the likelihood, not a second one.
     */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    CHECK(ws != NULL);

    size_t tip[4];
    tip[0] = gd_tree_find_label(t, "a"); tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c"); tip[3] = gd_tree_find_label(t, "d");
    int c0[7], c1[7];
    for (size_t i = 0; i < 7; i++) { c0[i] = -1; c1[i] = -1; }
    const int p0[4] = {1, 1, 0, 0}, p1[4] = {1, 0, 1, 0};
    for (size_t i = 0; i < 4; i++) { c0[tip[i]] = p0[i]; c1[tip[i]] = p1[i]; }
    gd_binary_class classes[2] = { {c0, 7}, {c1, 7} };
    size_t concept_ids[2] = { 0, 1 };

    gd_t2_options opt;
    gd_t2_options_default(&opt);

    const double lg = 0.0007, ll = 0.0009;
    double direct = 0.0;
    int rc_direct = gd_t2_loglik_dataset(ws, lg, ll, classes, concept_ids, 2,
                                         NULL, 2, &opt, NULL, NULL, &direct, &err);

    gd_mcmc_t2 *m = gd_mcmc_t2_create(ws, classes, concept_ids, 2, 2, &opt, &err);
    gd_mcmc_model model;
    memset(&model, 0, sizeof model);
    if (m != NULL) gd_mcmc_t2_model(m, &model);

    /* The sampler hands rates in log space. */
    double rates[2] = { log(lg), log(ll) };
    gd_mcmc_state st;
    memset(&st, 0, sizeof st);
    st.rates = rates;
    st.n_rates = 2;
    double via_adapter = 0.0;
    int rc_adapter = (model.loglik != NULL)
        ? model.loglik(model.ctx, t, &st, NULL, &via_adapter, &err) : GD_E_INTERNAL;

    /* Pointwise output must line up with the classes, one per class. */
    double pw[2] = { 0.0, 0.0 };
    double with_pw = 0.0;
    int rc_pw = (model.loglik != NULL)
        ? model.loglik(model.ctx, t, &st, pw, &with_pw, &err) : GD_E_INTERNAL;

    /* One rate is not enough for a two-rate likelihood, and that is refused
     * rather than read past the end of the vector. */
    st.n_rates = 1;
    gd_error err1 = GD_ERROR_INIT;
    double ignored = 0.0;
    int rc_one = (model.loglik != NULL)
        ? model.loglik(model.ctx, t, &st, NULL, &ignored, &err1) : GD_E_INTERNAL;

    size_t n_pointwise = model.n_pointwise;
    gd_mcmc_t2_destroy(m);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);

    CHECK(rc_direct == GD_E_OK);
    CHECK(m != NULL);
    CHECK(rc_adapter == GD_E_OK);
    CHECK_APPROX(via_adapter, direct, 1e-12);
    CHECK(rc_pw == GD_E_OK);
    CHECK_APPROX(with_pw, direct, 1e-12);
    CHECK_APPROX(pw[0] + pw[1], direct, 1e-9);
    CHECK(n_pointwise == 2);
    CHECK(rc_one == GD_E_INVALID_ARG);
}

TEST(test_t2_mcmc_adapter_rejects_bad_arguments)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_t2_ws *ws = gd_t2_ws_create(t, &err);
    int c0[7];
    for (size_t i = 0; i < 7; i++) c0[i] = -1;
    gd_binary_class classes[1] = { {c0, 7} };
    size_t ids[1] = { 0 };

    gd_mcmc_t2 *a = gd_mcmc_t2_create(NULL, classes, ids, 1, 1, NULL, &err);
    gd_mcmc_t2 *b = gd_mcmc_t2_create(ws, NULL, ids, 1, 1, NULL, &err);
    gd_mcmc_t2 *c = gd_mcmc_t2_create(ws, classes, ids, 0, 1, NULL, &err);
    gd_t2_ws_destroy(ws);
    gd_tree_destroy(t);
    gd_mcmc_t2_destroy(NULL);   /* allowed */

    CHECK(a == NULL);
    CHECK(b == NULL);
    CHECK(c == NULL);
}

static void run_tests(void)
{
    RUN(test_t2_transition_matrix);
    RUN(test_t2_bf_binary_all_patterns);
    RUN(test_t2_bf_missing_and_internal);
    RUN(test_t2_bf_root_prior_fixed);
    RUN(test_t2_bf_polytomy_native);
    RUN(test_t2_bf_randomized);
    RUN(test_t2_normalization);
    RUN(test_t2_ascertainment_effect);
    RUN(test_t2_dataset_sums_and_scaling);
    RUN(test_t2_partition_weights);
    RUN(test_t2_weight_from_string);
    RUN(test_t2_simulate_matches_likelihood);
    RUN(test_t2_simulation_recovers_rate);
    RUN(test_t2_bad_arguments);
    RUN(test_t2_ages_are_reread);
    RUN(test_t2_mcmc_adapter_matches_the_likelihood);
    RUN(test_t2_mcmc_adapter_rejects_bad_arguments);
}

TEST_MAIN()

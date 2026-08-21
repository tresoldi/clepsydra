/*
 * test_t1_likelihood.c — validation of the strict novel-state likelihood (T1)
 * for Clepsydra (M6).
 *
 * The acceptance gate for M6 (PLAN.md) is brute-force agreement on tiny
 * trees, where "brute force" must not reuse the p_v/q_v recursion: agreement
 * between two implementations of the same misunderstood recurrence proves
 * nothing (product brief §3.4, docs/validation.md §3).
 *
 * Two independent brute-force implementations are used here, and they are
 * cross-checked against each other as well as against the library:
 *
 *   bf_changes     — enumerates the 2^|E| replacement patterns of the
 *                    generative process directly (each replacement opens a
 *                    new class), and sums the probabilities of those whose
 *                    labelling matches the observed partition.
 *   bf_partitions  — enumerates *all* labellings of the nodes (as restricted
 *                    growth strings, i.e. all set partitions of V), keeps
 *                    those with the path property (Proposition 3) that match
 *                    the data, and multiplies the per-edge probabilities.
 *
 * Neither uses p_v/q_v, the data extension, LCA preprocessing, or any log
 * space helper other than logaddexp.
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_model_novel_state.h"
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_rng.h"
#include "gd_math.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Tolerance for brute-force agreement: the two routes differ only by
 * floating-point summation order, so this is generous. */
#define TOL 1e-9

/* ------------------------------------------------------------ brute force */

/* Independent ancestry helpers — plain parent walks, no LCA structure. */
static int is_ancestor(const gd_tree *t, size_t a, size_t d)
{
    for (size_t x = d; x != GD_NODE_NONE; x = t->nodes[x].parent) {
        if (x == a) return 1;
    }
    return 0;
}

static size_t lca_naive(const gd_tree *t, size_t a, size_t b)
{
    for (size_t x = a; x != GD_NODE_NONE; x = t->nodes[x].parent) {
        if (is_ancestor(t, x, b)) return x;
    }
    return GD_NODE_NONE;
}

/* Does the labelling g match the observed partition exactly?
 * (f(u) = f(w)  <=>  g(u) = g(w) for all observed u, w.) */
static int matches_data(const int *g, const gd_cognate_data *data)
{
    for (size_t i = 0; i < data->nnodes; i++) {
        if (data->labels[i] < 0) continue;
        for (size_t j = i + 1; j < data->nnodes; j++) {
            if (data->labels[j] < 0) continue;
            int data_same = (data->labels[i] == data->labels[j]);
            int g_same    = (g[i] == g[j]);
            if (data_same != g_same) return 0;
        }
    }
    return 1;
}

/* Path property (Proposition 3): whenever g(u) = g(w), every node on the
 * path u..w carries that class. */
static int has_path_property(const gd_tree *t, const int *g)
{
    for (size_t i = 0; i < t->nnodes; i++) {
        for (size_t j = i + 1; j < t->nnodes; j++) {
            if (g[i] != g[j]) continue;
            size_t l = lca_naive(t, i, j);
            for (size_t x = i; x != l; x = t->nodes[x].parent) {
                if (g[x] != g[i]) return 0;
            }
            for (size_t x = j; x != l; x = t->nodes[x].parent) {
                if (g[x] != g[i]) return 0;
            }
            if (g[l] != g[i]) return 0;
        }
    }
    return 1;
}

/* Product over edges of ret (child kept the parent class) or ch. */
static double edge_logprob(const gd_tree *t, double rate, const int *g)
{
    double lp = 0.0;
    for (size_t v = 0; v < t->nnodes; v++) {
        if (v == t->root) continue;
        size_t u = t->nodes[v].parent;
        double dt = t->nodes[u].age - t->nodes[v].age;
        if (dt < 0) dt = 0;
        double ret = exp(-rate * dt);
        lp += (g[u] == g[v]) ? log(ret) : log(1.0 - ret);
    }
    return lp;
}

/* Brute force 1: enumerate replacement patterns over the edges. */
static double bf_changes(const gd_tree *t, double rate,
                         const gd_cognate_data *data)
{
    size_t n = t->nnodes;
    int *g = (int *)malloc(n * sizeof(int));
    size_t *edge = (size_t *)malloc(n * sizeof(size_t));
    size_t nedges = 0;
    for (size_t i = 0; i < n; i++) {
        if (i != t->root) edge[nedges++] = i;
    }
    double total = -HUGE_VAL;
    size_t npat = (size_t)1 << nedges;
    for (size_t mask = 0; mask < npat; mask++) {
        g[t->root] = 0;
        int next_class = 1;
        double lp = 0.0;
        /* Edges are listed in index order, and a parent always has a smaller
         * index than its child in the trees built here; walk preorder to be
         * safe about the parent being labelled first. */
        for (size_t i = 0; i < n; i++) {
            size_t v = t->preorder[i];
            if (v == t->root) continue;
            size_t e = 0;
            while (edge[e] != v) e++;
            size_t u = t->nodes[v].parent;
            double dt = t->nodes[u].age - t->nodes[v].age;
            if (dt < 0) dt = 0;
            double ret = exp(-rate * dt);
            if ((mask >> e) & 1) {
                g[v] = next_class++;
                lp += log(1.0 - ret);
            } else {
                g[v] = g[u];
                lp += log(ret);
            }
        }
        if (matches_data(g, data)) total = gd_logaddexp(total, lp);
    }
    free(g);
    free(edge);
    return total;
}

/* Next restricted growth string: enumerates every set partition of n items. */
static int rgs_next(int *g, size_t n)
{
    if (n <= 1) return 0;
    for (size_t i = n; i-- > 1; ) {
        int mx = 0;
        for (size_t j = 0; j < i; j++) {
            if (g[j] > mx) mx = g[j];
        }
        if (g[i] < mx + 1) {
            g[i]++;
            for (size_t j = i + 1; j < n; j++) g[j] = 0;
            return 1;
        }
    }
    return 0;
}

/* Brute force 2: enumerate all labellings, filter by the path property. */
static double bf_partitions(const gd_tree *t, double rate,
                            const gd_cognate_data *data)
{
    size_t n = t->nnodes;
    int *g = (int *)calloc(n, sizeof(int));
    double total = -HUGE_VAL;
    do {
        if (!has_path_property(t, g)) continue;
        if (!matches_data(g, data)) continue;
        total = gd_logaddexp(total, edge_logprob(t, rate, g));
    } while (rgs_next(g, n));
    free(g);
    return total;
}

/* Relaxed-clock per-edge probability: the branch above child v uses the rate
 * scaled by the per-branch multiplier bm[v]. */
static double edge_logprob_bm(const gd_tree *t, double rate, const double *bm,
                              const int *g)
{
    double lp = 0.0;
    for (size_t v = 0; v < t->nnodes; v++) {
        if (v == t->root) continue;
        size_t u = t->nodes[v].parent;
        double dt = t->nodes[u].age - t->nodes[v].age;
        if (dt < 0) dt = 0;
        double ret = exp(-(rate * bm[v]) * dt);
        lp += (g[u] == g[v]) ? log(ret) : log(1.0 - ret);
    }
    return lp;
}

/* Relaxed brute force: enumerate all labellings, filter, per-branch rates. */
static double bf_partitions_bm(const gd_tree *t, double rate, const double *bm,
                               const gd_cognate_data *data)
{
    size_t n = t->nnodes;
    int *g = (int *)calloc(n, sizeof(int));
    double total = -HUGE_VAL;
    do {
        if (!has_path_property(t, g)) continue;
        if (!matches_data(g, data)) continue;
        total = gd_logaddexp(total, edge_logprob_bm(t, rate, bm, g));
    } while (rgs_next(g, n));
    free(g);
    return total;
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

/* (a,b,c)root — a hard polytomy. */
static gd_tree *tree_polytomy3(void)
{
    const char *nwk = "(a,b,c)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    set_age(t, "root", 3000);
    return t;
}

/* ((a,b,c)r1,d)root — a polytomy below the root, with different tip ages. */
static gd_tree *tree_polytomy_nested(void)
{
    const char *nwk = "((a,b,c)r1,d)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    set_age(t, "root", 5000);
    set_age(t, "r1", 3000);
    set_age(t, "b", 800);          /* an ancient (dated) tip */
    return t;
}

/* (a,b,c,d)root — a 4-way polytomy (two inserted ladder nodes). */
static gd_tree *tree_polytomy4(void)
{
    const char *nwk = "(a,b,c,d)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    set_age(t, "root", 2500);
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
    if (gd_tree_build_lca(t, err) != 0) goto fail;
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

/* Tree-component log-likelihood, contamination switched off. The exclude
 * policy is what lets a non-extendable concept come back as log 0 instead of
 * a diagnostic (a zero-weight contamination component is a dead end, and the
 * library says so rather than returning -inf in silence). */
static double t1_tree_loglik(gd_t1_ws *ws, double rate,
                             const gd_cognate_data *data, int *rc_out)
{
    gd_t1_options opt;
    gd_t1_concept_result res;
    gd_error err = GD_ERROR_INIT;
    gd_t1_options_default(&opt);
    opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    int rc = gd_t1_loglik_concept(ws, rate, data, &opt, &res, &err);
    if (rc_out) *rc_out = rc;
    if (rc != 0) return -HUGE_VAL;
    return res.log_lik_tree;
}

/* Compare the library against both brute forces on one dataset. */
static int agrees(gd_tree *t, double rate, const int *labels)
{
    gd_error err = GD_ERROR_INIT;
    gd_cognate_data data = { labels, t->nnodes };
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    if (ws == NULL) return 0;
    int rc = 0;
    double lib = t1_tree_loglik(ws, rate, &data, &rc);
    gd_t1_ws_destroy(ws);
    if (rc != 0) return 0;

    double b1 = bf_changes(t, rate, &data);
    double b2 = bf_partitions(t, rate, &data);
    int ok;
    if (fabs(b1 - b2) > TOL) ok = 0;                      /* the two agree */
    else if (b1 == -HUGE_VAL) ok = (lib == -HUGE_VAL);
    else ok = (fabs(lib - b1) <= TOL);
    if (!ok) {
        printf("    rate=%g lib=%.12g bf_changes=%.12g bf_partitions=%.12g\n"
               "    labels =", rate, lib, b1, b2);
        for (size_t i = 0; i < t->nnodes; i++) printf(" %d", labels[i]);
        printf("\n");
    }
    return ok;
}

/* ------------------------------------------------- brute-force agreement */

TEST(test_t1_bf_binary_one_class)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    int obs[5] = {-1, -1, 0, 0, 0};      /* root, r1 missing; a=b=c */
    CHECK(agrees(t, 0.001, obs));
    CHECK(agrees(t, 0.0001, obs));
    CHECK(agrees(t, 0.01, obs));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_binary_several_classes)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    int two[5]  = {-1, -1, 0, 0, 1};
    int all[5]  = {-1, -1, 0, 1, 2};
    int mix[5]  = {-1, -1, 0, 1, 0};
    CHECK(agrees(t, 0.001, two));
    CHECK(agrees(t, 0.0005, all));
    CHECK(agrees(t, 0.0005, mix));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_missing_leaf)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    int obs[5] = {-1, -1, 0, -1, 1};     /* b unobserved */
    CHECK(agrees(t, 0.001, obs));
    int none[5] = {-1, -1, -1, -1, -1};  /* nothing observed: likelihood 1 */
    CHECK(agrees(t, 0.001, none));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_observed_internal_node)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    size_t r1 = gd_tree_find_label(t, "r1");
    size_t root = gd_tree_find_label(t, "root");
    int obs[5];
    for (size_t i = 0; i < 5; i++) obs[i] = -1;
    obs[gd_tree_find_label(t, "a")] = 0;
    obs[gd_tree_find_label(t, "b")] = 1;
    obs[gd_tree_find_label(t, "c")] = 0;
    obs[r1] = 0;                          /* an attested ancestor */
    CHECK(agrees(t, 0.001, obs));
    obs[root] = 0;                        /* the root itself observed */
    CHECK(agrees(t, 0.001, obs));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_unobserved_ancestral_state)
{
    /* Only the internal node carries class 1; no tip does. The tree
     * likelihood must still account for it. */
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    int obs[7];
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[gd_tree_find_label(t, "r1")] = 1;
    obs[gd_tree_find_label(t, "a")] = 0;
    obs[gd_tree_find_label(t, "b")] = 2;
    obs[gd_tree_find_label(t, "c")] = 0;
    CHECK(agrees(t, 0.0004, obs));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_quartet_patterns)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    size_t a = gd_tree_find_label(t, "a"), b = gd_tree_find_label(t, "b");
    size_t c = gd_tree_find_label(t, "c"), d = gd_tree_find_label(t, "d");
    int obs[7];
    /* Every partition of the four tips, at two rates. */
    int part[4];
    for (size_t i = 0; i < 4; i++) part[i] = 0;
    do {
        for (size_t i = 0; i < 7; i++) obs[i] = -1;
        obs[a] = part[0]; obs[b] = part[1];
        obs[c] = part[2]; obs[d] = part[3];
        CHECK(agrees(t, 0.0003, obs));
        CHECK(agrees(t, 0.002, obs));
    } while (rgs_next(part, 4));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_polytomy)
{
    /* The ladder resolution is checked against brute-force enumeration on the
     * *unresolved* polytomy (docs/model-spec.md §3.3, §12 item 10). */
    gd_tree *t = tree_polytomy3();
    CHECK(t != NULL);
    size_t a = gd_tree_find_label(t, "a"), b = gd_tree_find_label(t, "b");
    size_t c = gd_tree_find_label(t, "c");
    int obs[4];
    int part[3];
    for (size_t i = 0; i < 3; i++) part[i] = 0;
    do {
        for (size_t i = 0; i < 4; i++) obs[i] = -1;
        obs[a] = part[0]; obs[b] = part[1]; obs[c] = part[2];
        CHECK(agrees(t, 0.0002, obs));
        CHECK(agrees(t, 0.0015, obs));
    } while (rgs_next(part, 3));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_polytomy_four_way)
{
    gd_tree *t = tree_polytomy4();
    CHECK(t != NULL);
    size_t idx[4];
    idx[0] = gd_tree_find_label(t, "a");
    idx[1] = gd_tree_find_label(t, "b");
    idx[2] = gd_tree_find_label(t, "c");
    idx[3] = gd_tree_find_label(t, "d");
    int obs[5];
    int part[4];
    for (size_t i = 0; i < 4; i++) part[i] = 0;
    do {
        for (size_t i = 0; i < 5; i++) obs[i] = -1;
        for (size_t i = 0; i < 4; i++) obs[idx[i]] = part[i];
        CHECK(agrees(t, 0.0004, obs));
    } while (rgs_next(part, 4));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_polytomy_nested_with_ancient_tip)
{
    gd_tree *t = tree_polytomy_nested();
    CHECK(t != NULL);
    size_t idx[4];
    idx[0] = gd_tree_find_label(t, "a");
    idx[1] = gd_tree_find_label(t, "b");
    idx[2] = gd_tree_find_label(t, "c");
    idx[3] = gd_tree_find_label(t, "d");
    int obs[6];
    int part[4];
    for (size_t i = 0; i < 4; i++) part[i] = 0;
    do {
        for (size_t i = 0; i < 6; i++) obs[i] = -1;
        for (size_t i = 0; i < 4; i++) obs[idx[i]] = part[i];
        CHECK(agrees(t, 0.00025, obs));
    } while (rgs_next(part, 4));
    /* Also with the polytomy's own parent observed. */
    for (size_t i = 0; i < 6; i++) obs[i] = -1;
    obs[gd_tree_find_label(t, "r1")] = 0;
    obs[idx[0]] = 0; obs[idx[1]] = 1; obs[idx[2]] = 0; obs[idx[3]] = 2;
    CHECK(agrees(t, 0.00025, obs));
    gd_tree_destroy(t);
}

TEST(test_t1_bf_randomized)
{
    /* Random ages and random partial data on binary and polytomous trees,
     * including many non-extendable patterns. */
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(20260731u, 7u, &err);
    CHECK(rng != NULL);
    const char *shapes[] = {
        "((a,b)r1,(c,d)r2)root;",
        "(((a,b)r1,c)r2,d)root;",
        "((a,b,c)r1,d)root;",
        "(a,b,(c,d)r1)root;",
        "((a,b)r1,c,d)root;"
    };
    int checked = 0;
    for (size_t s = 0; s < sizeof shapes / sizeof shapes[0]; s++) {
        for (int rep = 0; rep < 12; rep++) {
            gd_tree *t = NULL;
            CHECK(gd_newick_parse(shapes[s], strlen(shapes[s]), &t, &err) == 0);
            /* Ages: strictly decreasing with depth, jittered. */
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
            /* Keep parents strictly older than children. */
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
                obs[i] = (u < 0.25) ? -1 : (int)(gd_rng_uniform_closed_open(rng) * 3.0);
            }
            double rate = 0.0001 + gd_rng_uniform_closed_open(rng) * 0.002;
            CHECK(agrees(t, rate, obs));
            checked++;
            free(obs);
            gd_tree_destroy(t);
        }
    }
    CHECK(checked == 60);
    gd_rng_destroy(rng);
}

/* --------------------------------------------------------- normalization */

/* The tree likelihood is a probability distribution over the partitions of
 * the observed nodes, so it must sum to one over them. This is an absolute
 * check that needs no reference implementation at all. */
static double partition_mass(gd_tree *t, double rate, const size_t *obs_nodes,
                             size_t nobs, double epsilon)
{
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    if (ws == NULL) return -1.0;
    gd_t1_options opt;
    gd_t1_options_default(&opt);
    opt.epsilon = epsilon;
    /* With no contamination weight, sum the tree component alone (the
     * partitions that are not extendable contribute exactly zero). */
    if (epsilon <= 0.0) opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    int *labels = (int *)malloc(t->nnodes * sizeof(int));
    int *part = (int *)calloc(nobs, sizeof(int));
    double mass = 0.0;
    do {
        for (size_t i = 0; i < t->nnodes; i++) labels[i] = -1;
        for (size_t i = 0; i < nobs; i++) labels[obs_nodes[i]] = part[i];
        gd_cognate_data data = { labels, t->nnodes };
        gd_t1_concept_result res;
        if (gd_t1_loglik_concept(ws, rate, &data, &opt, &res, &err) != 0) {
            mass = -1.0;
            break;
        }
        mass += exp((epsilon <= 0.0) ? res.log_lik_tree : res.log_lik);
    } while (rgs_next(part, nobs));
    free(labels);
    free(part);
    gd_t1_ws_destroy(ws);
    return mass;
}

TEST(test_t1_normalization_binary)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    size_t obs[4];
    obs[0] = gd_tree_find_label(t, "a");
    obs[1] = gd_tree_find_label(t, "b");
    obs[2] = gd_tree_find_label(t, "c");
    obs[3] = gd_tree_find_label(t, "d");
    CHECK_APPROX(partition_mass(t, 0.0006, obs, 4, 0.0), 1.0, 1e-12);
    CHECK_APPROX(partition_mass(t, 0.004, obs, 4, 0.0), 1.0, 1e-12);
    /* Including an observed internal node in the observed set. */
    size_t obs5[5];
    obs5[0] = obs[0]; obs5[1] = obs[1]; obs5[2] = obs[2]; obs5[3] = obs[3];
    obs5[4] = gd_tree_find_label(t, "r1");
    CHECK_APPROX(partition_mass(t, 0.0006, obs5, 5, 0.0), 1.0, 1e-12);
    gd_tree_destroy(t);
}

TEST(test_t1_normalization_polytomy)
{
    gd_tree *t = tree_polytomy_nested();
    CHECK(t != NULL);
    size_t obs[4];
    obs[0] = gd_tree_find_label(t, "a");
    obs[1] = gd_tree_find_label(t, "b");
    obs[2] = gd_tree_find_label(t, "c");
    obs[3] = gd_tree_find_label(t, "d");
    CHECK_APPROX(partition_mass(t, 0.0005, obs, 4, 0.0), 1.0, 1e-12);
    gd_tree_destroy(t);
}

TEST(test_t1_noise_model_normalized)
{
    /* The CRP noise component lives on the same observable space, so it too
     * sums to one over the partitions of the observed nodes. */
    for (size_t nobs = 1; nobs <= 5; nobs++) {
        for (int ai = 0; ai < 3; ai++) {
            double alpha = (ai == 0) ? 0.4 : (ai == 1) ? 1.0 : 3.5;
            int *labels = (int *)calloc(nobs, sizeof(int));
            double mass = 0.0;
            do {
                gd_cognate_data data = { labels, nobs };
                mass += exp(gd_t1_noise_loglik(&data, alpha));
            } while (rgs_next(labels, nobs));
            free(labels);
            CHECK_APPROX(mass, 1.0, 1e-12);
        }
    }
}

TEST(test_t1_mixture_normalized)
{
    /* Both components are normalized over the partitions of the observed
     * nodes, so the contamination mixture is too — including the partitions
     * that are not extendable at all. */
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    size_t obs[4];
    obs[0] = gd_tree_find_label(t, "a");
    obs[1] = gd_tree_find_label(t, "b");
    obs[2] = gd_tree_find_label(t, "c");
    obs[3] = gd_tree_find_label(t, "d");
    CHECK_APPROX(partition_mass(t, 0.0006, obs, 4, 0.05), 1.0, 1e-12);
    CHECK_APPROX(partition_mass(t, 0.0006, obs, 4, 0.5), 1.0, 1e-12);
    gd_tree_destroy(t);
}

/* ---------------------------------------------- extension and diagnostics */

TEST(test_t1_extension_determined_nodes)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    size_t a = gd_tree_find_label(t, "a"), c = gd_tree_find_label(t, "c");
    size_t b = gd_tree_find_label(t, "b"), d = gd_tree_find_label(t, "d");
    size_t r1 = gd_tree_find_label(t, "r1"), r2 = gd_tree_find_label(t, "r2");
    size_t root = gd_tree_find_label(t, "root");

    /* a = c = 0 determines the whole path a-r1-root-r2-c. b and d differ,
     * and a class seen once determines only its own node. */
    int obs[7];
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[a] = 0; obs[c] = 0; obs[b] = 1; obs[d] = 2;
    gd_cognate_data data = { obs, 7 };
    gd_t1_extension_info info;
    CHECK(gd_t1_extension(ws, &data, &info, &err) == 0);
    CHECK(info.extendable == 1);
    CHECK(info.determined[a] == 0);
    CHECK(info.determined[r1] == 0);
    CHECK(info.determined[root] == 0);
    CHECK(info.determined[r2] == 0);
    CHECK(info.determined[c] == 0);
    CHECK(info.determined[b] == 1);
    CHECK(info.determined[d] == 2);
    CHECK(info.n_determined == 7);
    CHECK(info.conflict_node == GD_NODE_NONE);

    /* A singleton class determines only itself. */
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[a] = 5;
    CHECK(gd_t1_extension(ws, &data, &info, &err) == 0);
    CHECK(info.extendable == 1);
    CHECK(info.n_determined == 1);
    CHECK(info.determined[a] == 5);
    CHECK(info.determined[r1] == -1);

    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_incompatible_detected)
{
    /* a=c=0 forces the root to 0; b=d=1 forces it to 1. */
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[7];
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[gd_tree_find_label(t, "a")] = 0;
    obs[gd_tree_find_label(t, "b")] = 1;
    obs[gd_tree_find_label(t, "c")] = 0;
    obs[gd_tree_find_label(t, "d")] = 1;
    gd_cognate_data data = { obs, 7 };

    CHECK(gd_t1_is_compatible(ws, &data) == 0);
    gd_t1_extension_info info;
    CHECK(gd_t1_extension(ws, &data, &info, &err) == 0);
    CHECK(info.extendable == 0);
    CHECK(info.conflict_node != GD_NODE_NONE);
    CHECK((info.conflict_class_a == 0 && info.conflict_class_b == 1) ||
          (info.conflict_class_a == 1 && info.conflict_class_b == 0));

    /* The brute force agrees that no history produces this pattern. */
    CHECK(bf_changes(t, 0.001, &data) == -HUGE_VAL);
    CHECK(bf_partitions(t, 0.001, &data) == -HUGE_VAL);

    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_incompatible_policies)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[7];
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[gd_tree_find_label(t, "a")] = 0;
    obs[gd_tree_find_label(t, "b")] = 1;
    obs[gd_tree_find_label(t, "c")] = 0;
    obs[gd_tree_find_label(t, "d")] = 1;
    gd_cognate_data data = { obs, 7 };

    gd_t1_options opt;
    gd_t1_concept_result res;

    /* error: aborts with a diagnostic naming the node and both classes. */
    gd_t1_options_default(&opt);
    opt.incompatible = GD_INCOMPATIBLE_ERROR;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err)
          == GD_E_INCOMPATIBLE);
    CHECK(err.code == GD_E_INCOMPATIBLE);
    CHECK(err.message != NULL);
    CHECK(strstr(err.message, "not extendable") != NULL);
    CHECK(strstr(err.message, "class") != NULL);

    /* exclude: dropped, recorded, contributes nothing. */
    gd_t1_options_default(&opt);
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == 0);
    CHECK(res.excluded == 1);
    CHECK(res.compatible == 0);
    CHECK(res.log_lik == 0.0);

    /* contamination (default): the noise component carries the concept. */
    gd_t1_options_default(&opt);
    CHECK(opt.incompatible == GD_INCOMPATIBLE_CONTAMINATION);
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == 0);
    CHECK(res.compatible == 0);
    CHECK(res.excluded == 0);
    CHECK(res.log_lik_tree == -HUGE_VAL);
    CHECK(gd_is_finite(res.log_lik));
    CHECK_APPROX(res.log_lik, log(opt.epsilon) + res.log_lik_noise, 1e-12);
    CHECK_APPROX(res.p_contaminated, 1.0, 1e-12);

    /* contamination with zero weight is a dead end, and says so. */
    opt.epsilon = 0.0;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err)
          == GD_E_INCOMPATIBLE);
    CHECK(strstr(err.message, "epsilon") != NULL);

    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_contamination_mixture)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[5] = {-1, -1, 0, 0, 1};
    gd_cognate_data data = { obs, 5 };
    gd_t1_options opt;
    gd_t1_concept_result res;
    gd_t1_options_default(&opt);
    opt.epsilon = 0.1;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == 0);
    CHECK(res.compatible == 1);
    double expect = gd_logaddexp(log(0.9) + res.log_lik_tree,
                                 log(0.1) + res.log_lik_noise);
    CHECK_APPROX(res.log_lik, expect, 1e-12);
    CHECK_APPROX(res.p_contaminated,
                 exp(log(0.1) + res.log_lik_noise - res.log_lik), 1e-12);

    /* A concept the tree explains well is attributed to the tree; one it
     * explains badly is pushed towards the noise component. */
    int same[5] = {-1, -1, 0, 0, 0};
    int allx[5] = {-1, -1, 0, 1, 2};
    gd_cognate_data d_same = { same, 5 };
    gd_cognate_data d_allx = { allx, 5 };
    gd_t1_concept_result r_same, r_allx;
    CHECK(gd_t1_loglik_concept(ws, 0.0001, &d_same, &opt, &r_same, &err) == 0);
    CHECK(gd_t1_loglik_concept(ws, 0.0001, &d_allx, &opt, &r_allx, &err) == 0);
    CHECK(r_same.p_contaminated < 0.2);
    CHECK(r_allx.p_contaminated > r_same.p_contaminated);

    /* With epsilon = 0 the mixture disappears. */
    opt.epsilon = 0.0;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == 0);
    CHECK_APPROX(res.log_lik, res.log_lik_tree, 1e-15);
    CHECK(res.p_contaminated == 0.0);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_dataset_sums_concepts)
{
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    size_t tip[4];
    tip[0] = gd_tree_find_label(t, "a");
    tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c");
    tip[3] = gd_tree_find_label(t, "d");
    int c0[7], c1[7], c2[7];
    for (size_t i = 0; i < 7; i++) { c0[i] = -1; c1[i] = -1; c2[i] = -1; }
    /* a=b | c=d ; all four distinct ; a=c | b=d (not extendable). */
    const int p0[4] = {0, 0, 1, 1};
    const int p1[4] = {0, 1, 2, 3};
    const int p2[4] = {0, 1, 0, 1};
    for (size_t i = 0; i < 4; i++) {
        c0[tip[i]] = p0[i];
        c1[tip[i]] = p1[i];
        c2[tip[i]] = p2[i];
    }
    gd_cognate_data concepts[3] = { {c0, 7}, {c1, 7}, {c2, 7} };
    double rates[3] = { 0.0004, 0.0011, 0.0007 };
    gd_t1_concept_result per[3];
    gd_t1_options opt;
    gd_t1_options_default(&opt);
    double total = 0.0;
    CHECK(gd_t1_loglik_dataset(ws, rates, concepts, 3, &opt, per, &total, &err) == 0);
    CHECK_APPROX(total, per[0].log_lik + per[1].log_lik + per[2].log_lik, 1e-12);
    CHECK(per[0].compatible == 1);
    CHECK(per[1].compatible == 1);
    CHECK(per[2].compatible == 0);
    CHECK(per[2].p_contaminated == 1.0);
    /* The incompatible concept is reported, never silently dropped. */
    CHECK(gd_is_finite(total));

    /* Under exclude, the excluded concept contributes exactly nothing. */
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    double total_excl = 0.0;
    CHECK(gd_t1_loglik_dataset(ws, rates, concepts, 3, &opt, per, &total_excl, &err) == 0);
    CHECK(per[2].excluded == 1);
    CHECK_APPROX(total_excl, per[0].log_lik + per[1].log_lik, 1e-12);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

/* ------------------------------------------------------------- behaviour */

TEST(test_t1_rate_limits)
{
    /* rate -> 0: every node keeps the root class, so only the all-one-class
     * partition has any mass. rate -> infinity: every branch replaces, so
     * only the all-distinct partition does. */
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    int same[5] = {-1, -1, 0, 0, 0};
    int diff[5] = {-1, -1, 0, 1, 2};
    gd_cognate_data d_same = { same, 5 };
    gd_cognate_data d_diff = { diff, 5 };
    int rc = 0;
    CHECK_APPROX(exp(t1_tree_loglik(ws, 0.0, &d_same, &rc)), 1.0, 1e-12);
    CHECK(rc == 0);
    CHECK(t1_tree_loglik(ws, 0.0, &d_diff, &rc) == -HUGE_VAL);
    CHECK_APPROX(exp(t1_tree_loglik(ws, 5.0, &d_diff, &rc)), 1.0, 1e-9);
    CHECK(exp(t1_tree_loglik(ws, 5.0, &d_same, &rc)) < 1e-9);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_bad_arguments)
{
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[5] = {-1, -1, 0, 0, 1};
    gd_cognate_data data = { obs, 5 };
    gd_cognate_data wrong = { obs, 4 };
    gd_t1_concept_result res;
    gd_t1_options opt;
    gd_t1_options_default(&opt);
    CHECK(gd_t1_loglik_concept(ws, -1.0, &data, &opt, &res, &err) == GD_E_RANGE);
    CHECK(gd_t1_loglik_concept(ws, 0.001, &wrong, &opt, &res, &err)
          == GD_E_INVALID_ARG);
    CHECK(gd_t1_loglik_concept(NULL, 0.001, &data, &opt, &res, &err)
          == GD_E_INVALID_ARG);
    opt.epsilon = 1.0;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == GD_E_RANGE);
    gd_t1_options_default(&opt);
    opt.noise_alpha = 0.0;
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == GD_E_RANGE);
    /* An inverted chronology is reported, not silently clamped. */
    gd_t1_options_default(&opt);
    t->nodes[gd_tree_find_label(t, "r1")].age = 5000;   /* older than the root */
    CHECK(gd_t1_loglik_concept(ws, 0.001, &data, &opt, &res, &err) == GD_E_RANGE);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_ages_are_reread)
{
    /* The workspace must follow node ages as an MCMC chain moves them. */
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    int obs[5] = {-1, -1, 0, 0, 1};
    gd_cognate_data data = { obs, 5 };
    int rc = 0;
    double before = t1_tree_loglik(ws, 0.001, &data, &rc);
    t->nodes[gd_tree_find_label(t, "root")].age = 6000;
    double after = t1_tree_loglik(ws, 0.001, &data, &rc);
    CHECK(rc == 0);
    CHECK(fabs(before - after) > 1e-6);
    CHECK(agrees(t, 0.001, obs));           /* still matches brute force */
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

/* ------------------------------------------------------------- simulation */

TEST(test_t1_simulate_matches_likelihood)
{
    /* Monte-Carlo frequencies of each tip partition must match the
     * likelihood of that partition. */
    gd_tree *t = tree_simple();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    gd_rng *rng = gd_rng_create(11u, 3u, &err);
    CHECK(rng != NULL);
    size_t tip[3];
    tip[0] = gd_tree_find_label(t, "a");
    tip[1] = gd_tree_find_label(t, "b");
    tip[2] = gd_tree_find_label(t, "c");

    const size_t NREP = 40000;
    const double rate = 0.0004;
    /* Bell(3) = 5 partitions of the three tips, canonical RGS order. */
    size_t counts[5];
    for (size_t i = 0; i < 5; i++) counts[i] = 0;
    int sim[5];
    for (size_t rep = 0; rep < NREP; rep++) {
        CHECK(gd_t1_simulate(rng, t, rate, sim, &err) == 0);
        /* Canonicalize the tip labels into an RGS. */
        int canon[3], next = 0, seen[3];
        for (size_t i = 0; i < 3; i++) canon[i] = -1;
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < i; j++) {
                if (sim[tip[j]] == sim[tip[i]]) { canon[i] = canon[j]; break; }
            }
            if (canon[i] < 0) canon[i] = next++;
        }
        (void)seen;
        int part[3] = {0, 0, 0};
        size_t idx = 0, found = 5;
        do {
            if (part[0] == canon[0] && part[1] == canon[1] && part[2] == canon[2]) {
                found = idx;
                break;
            }
            idx++;
        } while (rgs_next(part, 3));
        CHECK(found < 5);
        counts[found]++;
    }

    int part[3] = {0, 0, 0};
    size_t idx = 0;
    int labels[5];
    do {
        for (size_t i = 0; i < 5; i++) labels[i] = -1;
        for (size_t i = 0; i < 3; i++) labels[tip[i]] = part[i];
        gd_cognate_data data = { labels, 5 };
        int rc = 0;
        double p = exp(t1_tree_loglik(ws, rate, &data, &rc));
        CHECK(rc == 0);
        double freq = (double)counts[idx] / (double)NREP;
        double se = sqrt(p * (1.0 - p) / (double)NREP);
        CHECK(fabs(freq - p) < 4.0 * se + 1e-12);
        idx++;
    } while (rgs_next(part, 3));
    CHECK(idx == 5);

    gd_rng_destroy(rng);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_simulation_recovers_rate)
{
    /* Simulate many concepts at a known rate, observe only the tips, and
     * check that the profile likelihood peaks near the truth. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_balanced(3, &err);      /* 15 nodes, 8 tips */
    CHECK(t != NULL);
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    gd_rng *rng = gd_rng_create(2026u, 11u, &err);
    CHECK(rng != NULL);

    const size_t NCONCEPT = 300;
    const double true_rate = 0.0005;
    size_t n = t->nnodes;
    int *sim = (int *)malloc(n * sizeof(int));
    int *store = (int *)malloc(NCONCEPT * n * sizeof(int));
    gd_cognate_data *concepts =
        (gd_cognate_data *)malloc(NCONCEPT * sizeof(gd_cognate_data));
    for (size_t c = 0; c < NCONCEPT; c++) {
        CHECK(gd_t1_simulate(rng, t, true_rate, sim, &err) == 0);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) {
            row[i] = t->nodes[i].is_leaf ? sim[i] : -1;   /* tips only */
        }
        concepts[c].labels = row;
        concepts[c].nnodes = n;
    }

    gd_t1_options opt;
    gd_t1_options_default(&opt);
    opt.epsilon = 0.0;
    double *rates = (double *)malloc(NCONCEPT * sizeof(double));
    double best_ll = -HUGE_VAL, best_rate = 0.0;
    for (int k = 0; k < 41; k++) {
        double r = 0.0001 * pow(10.0, 1.5 * (double)k / 40.0); /* 1e-4 .. ~3e-3 */
        for (size_t c = 0; c < NCONCEPT; c++) rates[c] = r;
        double ll = 0.0;
        CHECK(gd_t1_loglik_dataset(ws, rates, concepts, NCONCEPT, &opt,
                                   NULL, &ll, &err) == 0);
        if (ll > best_ll) { best_ll = ll; best_rate = r; }
    }
    CHECK(best_rate > 0.7 * true_rate);
    CHECK(best_rate < 1.4 * true_rate);

    free(rates);
    free(concepts);
    free(store);
    free(sim);
    gd_rng_destroy(rng);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

/* ---------------------------------------------------------------- scaling */

TEST(test_t1_scaling_is_linear)
{
    /* Acceptance gate: ~linear scaling in nodes and in concepts for
     * compatible data. The recursion visits each working node exactly once;
     * the extension is bounded by a small multiple of the node count even in
     * its worst case (every tip in one class, so the spanning subtree is the
     * whole tree). */
    gd_error err = GD_ERROR_INIT;
    for (size_t depth = 2; depth <= 8; depth++) {
        gd_tree *t = tree_balanced(depth, &err);
        CHECK(t != NULL);
        gd_t1_ws *ws = gd_t1_ws_create(t, &err);
        CHECK(ws != NULL);
        size_t n = t->nnodes;
        CHECK(gd_t1_ws_working_nodes(ws) == n);   /* binary: nothing inserted */
        int *obs = (int *)malloc(n * sizeof(int));
        for (size_t i = 0; i < n; i++) obs[i] = t->nodes[i].is_leaf ? 0 : -1;
        gd_cognate_data data = { obs, n };
        gd_t1_concept_result res;
        gd_t1_options opt;
        gd_t1_options_default(&opt);
        CHECK(gd_t1_loglik_concept(ws, 0.0005, &data, &opt, &res, &err) == 0);
        CHECK(res.compatible == 1);
        CHECK(res.stability_warning == 0);
        CHECK(res.ops_recursion == n);
        CHECK(res.ops_extension <= 6 * n);
        free(obs);
        gd_t1_ws_destroy(ws);
        gd_tree_destroy(t);
    }
}

TEST(test_t1_scaling_in_concepts)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_balanced(5, &err);
    CHECK(t != NULL);
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    size_t n = t->nnodes;
    const size_t NC = 64;
    int *store = (int *)malloc(NC * n * sizeof(int));
    gd_cognate_data *concepts =
        (gd_cognate_data *)malloc(NC * sizeof(gd_cognate_data));
    double *rates = (double *)malloc(NC * sizeof(double));
    gd_t1_concept_result *per =
        (gd_t1_concept_result *)malloc(NC * sizeof(gd_t1_concept_result));
    for (size_t c = 0; c < NC; c++) {
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) {
            row[i] = t->nodes[i].is_leaf ? (int)(i % 3) : -1;
        }
        concepts[c].labels = row;
        concepts[c].nnodes = n;
        rates[c] = 0.0004;
    }
    gd_t1_options opt;
    gd_t1_options_default(&opt);
    double ll = 0.0;
    CHECK(gd_t1_loglik_dataset(ws, rates, concepts, NC, &opt, per, &ll, &err) == 0);
    /* Identical concepts: the total is exactly NC times one of them, and the
     * per-concept work is identical — the cost is linear in the concept
     * count by construction. */
    CHECK_APPROX(ll, (double)NC * per[0].log_lik, 1e-9);
    for (size_t c = 1; c < NC; c++) {
        CHECK(per[c].ops_recursion == per[0].ops_recursion);
        CHECK(per[c].ops_extension == per[0].ops_extension);
    }
    free(per); free(rates); free(concepts); free(store);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_policy_from_string)
{
    gd_error err = GD_ERROR_INIT;
    gd_incompatible_policy p;
    CHECK(gd_t1_policy_from_string("error", &p, &err) == 0);
    CHECK(p == GD_INCOMPATIBLE_ERROR);
    CHECK(gd_t1_policy_from_string("exclude", &p, &err) == 0);
    CHECK(p == GD_INCOMPATIBLE_EXCLUDE);
    CHECK(gd_t1_policy_from_string("contamination", &p, &err) == 0);
    CHECK(p == GD_INCOMPATIBLE_CONTAMINATION);
    CHECK(gd_t1_policy_from_string(NULL, &p, &err) == 0);
    CHECK(p == GD_INCOMPATIBLE_CONTAMINATION);          /* the default */
    /* An unknown setting is refused rather than guessed. */
    CHECK(gd_t1_policy_from_string("drop", &p, &err) == GD_E_INVALID_ARG);
    CHECK(gd_t1_policy_from_string("Exclude", &p, &err) == GD_E_INVALID_ARG);
}

TEST(test_t1_polytomy_inserts_ladder_nodes)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_polytomy4();
    CHECK(t != NULL);
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    /* 5 original nodes, a 4-way polytomy: two nodes inserted. */
    CHECK(gd_t1_ws_working_nodes(ws) == t->nnodes + 2);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

/* ------------------------------------------------- relaxed branch clock (M9) */

/* Relaxed tree-component loglik (contamination off). */
static double t1_tree_loglik_bm(gd_t1_ws *ws, double rate, const double *bm,
                                const gd_cognate_data *data, int *rc_out)
{
    gd_t1_options opt;
    gd_t1_concept_result res;
    gd_error err = GD_ERROR_INIT;
    gd_t1_options_default(&opt);
    opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    int rc = gd_t1_loglik_concept_relaxed(ws, rate, bm, data, &opt, &res, &err);
    if (rc_out) *rc_out = rc;
    if (rc != 0) return -HUGE_VAL;
    return res.log_lik_tree;
}

static int agrees_bm(gd_tree *t, double rate, const double *bm, const int *labels)
{
    gd_error err = GD_ERROR_INIT;
    gd_cognate_data data = { labels, t->nnodes };
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    if (ws == NULL) return 0;
    int rc = 0;
    double lib = t1_tree_loglik_bm(ws, rate, bm, &data, &rc);
    gd_t1_ws_destroy(ws);
    if (rc != 0) return 0;
    double bf = bf_partitions_bm(t, rate, bm, &data);
    int ok;
    if (bf == -HUGE_VAL) ok = (lib == -HUGE_VAL);
    else ok = (fabs(lib - bf) <= TOL);
    if (!ok) printf("    bm: rate=%g lib=%.12g bf=%.12g\n", rate, lib, bf);
    return ok;
}

TEST(test_t1_relaxed_strict_equivalence)
{
    /* All-ones multipliers reproduce the strict clock exactly. */
    gd_tree *t = tree_quartet();
    CHECK(t != NULL);
    gd_error err = GD_ERROR_INIT;
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    CHECK(ws != NULL);
    double bm[7];
    for (size_t i = 0; i < 7; i++) bm[i] = 1.0;
    int obs[7];
    for (size_t i = 0; i < 7; i++) obs[i] = -1;
    obs[gd_tree_find_label(t, "a")] = 0;
    obs[gd_tree_find_label(t, "b")] = 0;
    obs[gd_tree_find_label(t, "c")] = 1;
    obs[gd_tree_find_label(t, "d")] = 1;
    gd_cognate_data data = { obs, 7 };
    int rc1 = 0, rc2 = 0;
    double strict = t1_tree_loglik(ws, 0.0006, &data, &rc1);
    double relaxed = t1_tree_loglik_bm(ws, 0.0006, bm, &data, &rc2);
    CHECK(rc1 == 0 && rc2 == 0);
    CHECK_APPROX(strict, relaxed, 1e-12);
    gd_t1_ws_destroy(ws);
    gd_tree_destroy(t);
}

TEST(test_t1_relaxed_matches_bruteforce)
{
    gd_error err = GD_ERROR_INIT;
    gd_rng *rng = gd_rng_create(20260801u, 4u, &err);
    CHECK(rng != NULL);
    const char *shapes[] = {
        "((a,b)r1,(c,d)r2)root;",
        "(((a,b)r1,c)r2,d)root;",
        "((a,b,c)r1,d)root;"          /* polytomy: zero-length ladder edges */
    };
    int checked = 0;
    for (size_t s = 0; s < sizeof shapes / sizeof shapes[0]; s++) {
        for (int rep = 0; rep < 12; rep++) {
            gd_tree *t = NULL;
            CHECK(gd_newick_parse(shapes[s], strlen(shapes[s]), &t, &err) == 0);
            for (size_t i = 0; i < t->nnodes; i++) {
                size_t d = 0, x = i;
                while (t->nodes[x].parent != GD_NODE_NONE) { x = t->nodes[x].parent; d++; }
                t->nodes[i].age = t->nodes[i].nchildren == 0
                    ? gd_rng_uniform_closed_open(rng) * 300.0
                    : 3600.0 - 800.0 * (double)d + gd_rng_uniform_closed_open(rng) * 150.0;
            }
            for (size_t i = 0; i < t->nnodes; i++) {
                size_t v = t->postorder[i];
                size_t p = t->nodes[v].parent;
                if (p != GD_NODE_NONE && t->nodes[p].age <= t->nodes[v].age)
                    t->nodes[p].age = t->nodes[v].age + 100.0;
            }
            /* Per-branch multipliers ~ exp(N(0, 0.5)); wide, so branches differ. */
            double *bm = (double *)malloc(t->nnodes * sizeof(double));
            for (size_t i = 0; i < t->nnodes; i++)
                bm[i] = exp(gd_rng_normal(rng, 0.0, 0.5));
            int *obs = (int *)malloc(t->nnodes * sizeof(int));
            for (size_t i = 0; i < t->nnodes; i++) {
                double u = gd_rng_uniform_closed_open(rng);
                obs[i] = (u < 0.25) ? -1 : (int)(gd_rng_uniform_closed_open(rng) * 3.0);
            }
            double rate = 0.0002 + gd_rng_uniform_closed_open(rng) * 0.0015;
            CHECK(agrees_bm(t, rate, bm, obs));
            checked++;
            free(bm); free(obs);
            gd_tree_destroy(t);
        }
    }
    CHECK(checked == 36);
    gd_rng_destroy(rng);
}

static void run_tests(void)
{
    RUN(test_t1_bf_binary_one_class);
    RUN(test_t1_bf_binary_several_classes);
    RUN(test_t1_bf_missing_leaf);
    RUN(test_t1_bf_observed_internal_node);
    RUN(test_t1_bf_unobserved_ancestral_state);
    RUN(test_t1_bf_quartet_patterns);
    RUN(test_t1_bf_polytomy);
    RUN(test_t1_bf_polytomy_four_way);
    RUN(test_t1_bf_polytomy_nested_with_ancient_tip);
    RUN(test_t1_bf_randomized);
    RUN(test_t1_normalization_binary);
    RUN(test_t1_normalization_polytomy);
    RUN(test_t1_noise_model_normalized);
    RUN(test_t1_mixture_normalized);
    RUN(test_t1_extension_determined_nodes);
    RUN(test_t1_incompatible_detected);
    RUN(test_t1_incompatible_policies);
    RUN(test_t1_contamination_mixture);
    RUN(test_t1_dataset_sums_concepts);
    RUN(test_t1_rate_limits);
    RUN(test_t1_bad_arguments);
    RUN(test_t1_ages_are_reread);
    RUN(test_t1_simulate_matches_likelihood);
    RUN(test_t1_simulation_recovers_rate);
    RUN(test_t1_scaling_is_linear);
    RUN(test_t1_scaling_in_concepts);
    RUN(test_t1_policy_from_string);
    RUN(test_t1_polytomy_inserts_ladder_nodes);
    RUN(test_t1_relaxed_strict_equivalence);
    RUN(test_t1_relaxed_matches_bruteforce);
}

TEST_MAIN()

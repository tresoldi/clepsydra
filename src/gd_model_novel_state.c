/*
 * gd_model_novel_state.c — strict novel-state multistate likelihood (T1).
 *
 * Independent re-implementation of the Rönchen et al. 2024 algorithm:
 *   1. data extension (Proposition 7): mark every determined node, detect
 *      non-extendable data;
 *   2. the p_v/q_v double recursion (Theorem 10, cases i/iia/iib/iic);
 *   3. the robust contamination mixture of docs/model-spec.md §6.
 *
 * All arithmetic is in log space (gd_math helpers); every traversal is
 * iterative (browser stack safety); the likelihood allocates nothing (all
 * scratch lives in gd_t1_ws). See docs/model-spec.md §10 for the pseudocode
 * this follows and §11 for the independent brute-force check that validates
 * it.
 *
 * Polytomies: a node with k > 2 children is resolved, once, into a ladder of
 * k-2 inserted nodes carrying the polytomy's own age, so every inserted
 * branch has duration zero. The generative process of §3.3 is unchanged by
 * that rewriting — each child still independently retains the parent class
 * with probability exp(-r t_child) or opens a new class — so the binary
 * recursion applies to the ladder. This is the documented resolution of
 * §12 item 10; the exact nonbinary recursion (Thm A.8/A.9) still awaits the
 * appendix supplement, and tests/likelihood/test_t1_likelihood.c checks the
 * ladder against brute-force enumeration on the *unresolved* polytomy.
 */
#include "gd_model_novel_state.h"
#include "gd_memory.h"
#include "gd_tree.h"
#include "gd_tree_build.h"
#include "gd_math.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Largest tolerated negative value of a p - q*ret cancellation (log scale)
 * before it is reported as a precision problem rather than clamped in
 * silence. Mathematically these differences are non-negative. */
#define GD_T1_CANCEL_TOL 1e-9

/* Ages must be non-increasing from parent to child; anything below this is a
 * genuine chronology inversion rather than rounding. */
#define GD_T1_AGE_TOL 1e-9

typedef struct t1_key {
    int    label;
    size_t node;
} t1_key;

struct gd_t1_ws {
    gd_tree *tree;
    size_t   n;              /* nodes of the original tree */
    size_t   m;              /* nodes of the binary working tree (n + inserted) */
    size_t   root;

    /* Working tree (indices < n are the original nodes, unchanged). */
    size_t  *parent;         /* m */
    size_t  *left;           /* m; GD_NODE_NONE if none */
    size_t  *right;          /* m; GD_NODE_NONE if none (leaf or unary node) */
    size_t  *orig;           /* m -> original node; inserted nodes map to the
                              * polytomy they resolve */
    unsigned char *inserted; /* m; 1 for ladder nodes */
    size_t  *post;           /* m; postorder, children before parents */
    size_t  *stack;          /* m; traversal scratch */
    double  *age;            /* m; refreshed from the tree on every call */

    /* Per-call scratch. */
    int     *det;            /* m; determined class or -1 */
    double  *logp;           /* m */
    double  *logq;           /* m */
    t1_key  *keys;           /* n; observed (label, node) sorted by label */
    size_t   nkeys;          /* number of observed nodes in `keys` */
    int     *det_pub;        /* n; extension mapped back to original nodes */
};

/* ------------------------------------------------------------------ helpers */

void gd_t1_options_default(gd_t1_options *opt)
{
    if (opt == NULL) return;
    opt->incompatible = GD_INCOMPATIBLE_CONTAMINATION;
    opt->epsilon      = 0.05;   /* mean of the default Beta(1,19) prior (§6) */
    opt->noise_alpha  = 1.0;
}

int gd_t1_policy_from_string(const char *s, gd_incompatible_policy *out,
                             gd_error *err)
{
    if (out == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1: NULL policy output");
        return GD_E_INVALID_ARG;
    }
    if (s == NULL || s[0] == '\0') {
        *out = GD_INCOMPATIBLE_CONTAMINATION;
        return 0;
    }
    if (strcmp(s, "error") == 0)         { *out = GD_INCOMPATIBLE_ERROR; return 0; }
    if (strcmp(s, "exclude") == 0)       { *out = GD_INCOMPATIBLE_EXCLUDE; return 0; }
    if (strcmp(s, "contamination") == 0) { *out = GD_INCOMPATIBLE_CONTAMINATION; return 0; }
    gd_error_set(err, GD_E_INVALID_ARG,
                 "t1: incompatible must be error, exclude or contamination");
    return GD_E_INVALID_ARG;
}

static void result_init(gd_t1_concept_result *res)
{
    res->log_lik           = -HUGE_VAL;
    res->log_lik_tree      = -HUGE_VAL;
    res->log_lik_noise     = -HUGE_VAL;
    res->p_contaminated    = 0.0;
    res->compatible        = true;
    res->excluded          = false;
    res->ops_extension     = 0;
    res->ops_recursion     = 0;
    res->stability_warning = false;
}

/* log(exp(a) - exp(b)) where the difference is non-negative by construction.
 * Clamps a vanishing difference to log 0 and flags a measurable violation. */
static double t1_logdiff(double a, double b, bool *stability)
{
    if (b == -HUGE_VAL) return a;
    if (a <= b) {
        if (a > -HUGE_VAL && (b - a) > GD_T1_CANCEL_TOL) *stability = true;
        return -HUGE_VAL;
    }
    return gd_logdiffexp(a, b);
}

static int key_cmp(const void *va, const void *vb)
{
    const t1_key *a = (const t1_key *)va;
    const t1_key *b = (const t1_key *)vb;
    if (a->label != b->label) return (a->label < b->label) ? -1 : 1;
    if (a->node  != b->node)  return (a->node  < b->node)  ? -1 : 1;
    return 0;
}

/* --------------------------------------------------------------- workspace */

static void ws_free_arrays(gd_t1_ws *ws)
{
    free(ws->parent);   free(ws->left);   free(ws->right);
    free(ws->orig);     free(ws->inserted);
    free(ws->post);     free(ws->stack);  free(ws->age);
    free(ws->det);      free(ws->logp);   free(ws->logq);
    free(ws->keys);     free(ws->det_pub);
}

gd_t1_ws *gd_t1_ws_create(gd_tree *tree, gd_error *err)
{
    if (tree == NULL || tree->nnodes == 0 || tree->nodes == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1: NULL or empty tree");
        return NULL;
    }
    if (!tree->lca_built && gd_tree_build_lca(tree, err) != 0) return NULL;

    size_t n = tree->nnodes;
    size_t extra = 0;
    for (size_t i = 0; i < n; i++) {
        size_t k = tree->nodes[i].nchildren;
        if (k > 2) extra += k - 2;
    }
    size_t m = n + extra;

    gd_t1_ws *ws = (gd_t1_ws *)gd_calloc_checked(1, sizeof(*ws), err);
    if (ws == NULL) {
        gd_error_set(err, GD_E_NOMEM, "t1: workspace alloc");
        return NULL;
    }
    ws->tree = tree;
    ws->n = n;
    ws->m = m;
    ws->root = tree->root;

    ws->parent   = (size_t *)gd_malloc_array_checked(m, sizeof(size_t), err);
    ws->left     = (size_t *)gd_malloc_array_checked(m, sizeof(size_t), err);
    ws->right    = (size_t *)gd_malloc_array_checked(m, sizeof(size_t), err);
    ws->orig     = (size_t *)gd_malloc_array_checked(m, sizeof(size_t), err);
    ws->inserted = (unsigned char *)gd_malloc_array_checked(m, sizeof(unsigned char), err);
    ws->post     = (size_t *)gd_malloc_array_checked(m, sizeof(size_t), err);
    ws->stack    = (size_t *)gd_malloc_array_checked(m, sizeof(size_t), err);
    ws->age      = (double *)gd_malloc_array_checked(m, sizeof(double), err);
    ws->det      = (int *)gd_malloc_array_checked(m, sizeof(int), err);
    ws->logp     = (double *)gd_malloc_array_checked(m, sizeof(double), err);
    ws->logq     = (double *)gd_malloc_array_checked(m, sizeof(double), err);
    ws->keys     = (t1_key *)gd_malloc_array_checked(n, sizeof(t1_key), err);
    ws->det_pub  = (int *)gd_malloc_array_checked(n, sizeof(int), err);
    if (ws->parent == NULL || ws->left == NULL || ws->right == NULL ||
        ws->orig == NULL || ws->inserted == NULL || ws->post == NULL ||
        ws->stack == NULL || ws->age == NULL || ws->det == NULL ||
        ws->logp == NULL || ws->logq == NULL || ws->keys == NULL ||
        ws->det_pub == NULL) {
        ws_free_arrays(ws);
        free(ws);
        gd_error_set(err, GD_E_NOMEM, "t1: workspace arrays");
        return NULL;
    }

    for (size_t i = 0; i < m; i++) {
        ws->parent[i]   = GD_NODE_NONE;
        ws->left[i]     = GD_NODE_NONE;
        ws->right[i]    = GD_NODE_NONE;
        ws->orig[i]     = (i < n) ? i : GD_NODE_NONE;
        ws->inserted[i] = (i < n) ? (unsigned char)0 : (unsigned char)1;
    }

    /* Wire the working tree; resolve polytomies into zero-length ladders. */
    size_t next = n;
    for (size_t u = 0; u < n; u++) {
        size_t k = tree->nodes[u].nchildren;
        const size_t *ch = tree->nodes[u].children;
        if (k == 0) continue;
        if (k == 1) {
            ws->left[u] = ch[0];
            ws->parent[ch[0]] = u;
            continue;
        }
        if (k == 2) {
            ws->left[u]  = ch[0];
            ws->right[u] = ch[1];
            ws->parent[ch[0]] = u;
            ws->parent[ch[1]] = u;
            continue;
        }
        size_t cur = u;
        for (size_t i = 0; i + 2 < k; i++) {
            size_t h = next++;
            ws->orig[h] = u;              /* the polytomy it stands for */
            ws->left[cur]   = ch[i];
            ws->parent[ch[i]] = cur;
            ws->right[cur]  = h;
            ws->parent[h]   = cur;
            cur = h;
        }
        ws->left[cur]  = ch[k - 2];
        ws->right[cur] = ch[k - 1];
        ws->parent[ch[k - 2]] = cur;
        ws->parent[ch[k - 1]] = cur;
    }

    /* Postorder: reverse of the (root, right, left) preorder. Iterative. */
    size_t sp = 0, idx = m;
    ws->stack[sp++] = ws->root;
    while (sp > 0) {
        size_t u = ws->stack[--sp];
        ws->post[--idx] = u;
        if (ws->left[u]  != GD_NODE_NONE) ws->stack[sp++] = ws->left[u];
        if (ws->right[u] != GD_NODE_NONE) ws->stack[sp++] = ws->right[u];
    }
    if (idx != 0) {   /* the tree was not connected / root was wrong */
        ws_free_arrays(ws);
        free(ws);
        gd_error_set(err, GD_E_INVALID_ARG, "t1: tree is not a single rooted tree");
        return NULL;
    }
    return ws;
}

void gd_t1_ws_destroy(gd_t1_ws *ws)
{
    if (ws == NULL) return;
    ws_free_arrays(ws);
    free(ws);
}

size_t gd_t1_ws_working_nodes(const gd_t1_ws *ws)
{
    return (ws != NULL) ? ws->m : 0;
}

/* ------------------------------------------------------- data extension */

/*
 * Mark every determined node (Proposition 7). For each cognate class, the set
 * of nodes determined to it is the minimal subtree spanning the observed
 * nodes carrying it, i.e. the union of the paths from each of them up to
 * their common ancestor — so one LCA per class plus one upward walk per
 * observed node suffices, and each walk stops as soon as it reaches ground
 * already covered for that class. Observed nodes are marked by their own walk
 * (a class observed exactly once determines just its own node), which is what
 * dom(f~) ⊇ dom(f_m) requires.
 *
 * Returns 1 if extendable, 0 otherwise (with the conflict reported).
 */
/* Predicate, not a status: true when the data extends onto the tree. The two
 * false returns are the two ways it does not, and *conflict_node names where. */
static bool t1_extend(gd_t1_ws *ws, const gd_cognate_data *data, size_t *ops,
                      size_t *conflict_node, int *class_a, int *class_b)
{
    size_t n = ws->n;
    size_t work = 0;

    for (size_t i = 0; i < ws->m; i++) ws->det[i] = -1;
    *conflict_node = GD_NODE_NONE;
    *class_a = -1;
    *class_b = -1;

    size_t nk = 0;
    for (size_t i = 0; i < n; i++) {
        if (data->labels[i] >= 0) {
            ws->keys[nk].label = data->labels[i];
            ws->keys[nk].node  = i;
            nk++;
        }
    }
    ws->nkeys = nk;
    if (nk > 1) qsort(ws->keys, nk, sizeof(t1_key), key_cmp);

    size_t s = 0;
    while (s < nk) {
        size_t e = s + 1;
        while (e < nk && ws->keys[e].label == ws->keys[s].label) e++;
        int label = ws->keys[s].label;

        /* Common ancestor of the whole class (O(1) per pair via the tree's
         * Euler-tour RMQ; observed nodes are always original nodes). */
        size_t anc = ws->keys[s].node;
        for (size_t j = s + 1; j < e; j++) {
            anc = gd_tree_mrca(ws->tree, anc, ws->keys[j].node);
            work++;
            if (anc == GD_NODE_NONE) return false;  /* cannot happen on a tree */
        }

        for (size_t j = s; j < e; j++) {
            size_t cur = ws->keys[j].node;
            for (;;) {
                work++;
                if (ws->det[cur] == label) break;      /* already covered */
                if (ws->det[cur] >= 0) {               /* two classes: no f~ */
                    *conflict_node = ws->inserted[cur] ? ws->orig[cur] : cur;
                    *class_a = ws->det[cur];
                    *class_b = label;
                    *ops += work;
                    return false;
                }
                ws->det[cur] = label;
                if (cur == anc) break;
                cur = ws->parent[cur];
                if (cur == GD_NODE_NONE) break;
            }
        }
        s = e;
    }
    *ops += work;
    return true;
}

int gd_t1_extension(gd_t1_ws *ws, const gd_cognate_data *data,
                    gd_t1_extension_info *info, gd_error *err)
{
    if (ws == NULL || data == NULL || info == NULL || data->labels == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1_extension: NULL argument");
        return GD_E_INVALID_ARG;
    }
    if (data->nnodes != ws->n) {
        gd_error_set(err, GD_E_INVALID_ARG,
                     "t1_extension: data length does not match the tree");
        return GD_E_INVALID_ARG;
    }
    size_t ops = 0, conflict = GD_NODE_NONE;
    int ca = -1, cb = -1;
    bool ok = t1_extend(ws, data, &ops, &conflict, &ca, &cb);

    size_t ndet = 0;
    for (size_t i = 0; i < ws->n; i++) {
        ws->det_pub[i] = ws->det[i];
        if (ws->det[i] >= 0) ndet++;
    }
    info->determined      = ws->det_pub;
    info->nnodes          = ws->n;
    info->n_determined    = ndet;
    info->extendable      = ok;
    info->conflict_node   = conflict;
    info->conflict_class_a = ca;
    info->conflict_class_b = cb;
    return 0;
}

bool gd_t1_is_compatible(gd_t1_ws *ws, const gd_cognate_data *data)
{
    gd_t1_extension_info info;
    gd_error err = GD_ERROR_INIT;
    if (gd_t1_extension(ws, data, &info, &err) != GD_E_OK) return false;
    return info.extendable != 0;
}

/* ----------------------------------------------------------- the recursion */

/* Branch quantities for the edge above working node w. */
/* Predicate: true when the branch above node w has a usable duration. False
 * means a chronology inversion, which the caller turns into a real error. */
static bool branch_logs(const gd_t1_ws *ws, size_t w, double rate,
                       double *log_ret, double *log_ch)
{
    double dt = ws->age[ws->parent[w]] - ws->age[w];
    if (dt < 0.0) {
        if (dt < -GD_T1_AGE_TOL) return false;   /* chronology inversion */
        dt = 0.0;
    }
    double x = -rate * dt;                   /* = log ret, exactly */
    *log_ret = x;
    *log_ch  = (x < 0.0) ? gd_log1mexp(x) : -HUGE_VAL;
    return true;
}

/*
 * Bottom-up p_v/q_v recursion over the binary working tree (Theorem 10).
 *   p_v = P(the subtree of v is compatible with f~)
 *   q_v = P(that, and v still carries a class that some determined node of
 *          its subtree carries)
 * A missing second child is treated as a child with p = 1, q = 0 and no
 * determination, which reproduces the correct unary-node formulas.
 */
static int t1_recursion(gd_t1_ws *ws, double rate, const double *bm,
                        double *out_logp,
                        bool *stability, size_t *ops, gd_error *err)
{
    /* Ages are re-read every call: an MCMC chain moves them under us. */
    for (size_t i = 0; i < ws->m; i++) {
        ws->age[i] = ws->tree->nodes[ws->orig[i]].age;
    }

    for (size_t i = 0; i < ws->m; i++) {
        size_t u = ws->post[i];
        size_t l = ws->left[u];
        size_t r = ws->right[u];
        int det_u = ws->det[u];

        if (l == GD_NODE_NONE) {                       /* leaf */
            ws->logp[u] = 0.0;
            ws->logq[u] = (det_u >= 0) ? 0.0 : -HUGE_VAL;
            continue;
        }

        /* Relaxed clock: the rate on the edge above child w is scaled by the
         * per-branch multiplier bm[orig(w)] (bm == NULL => strict clock).
         * Inserted ladder nodes carry zero-duration edges, so their multiplier
         * never matters. */
        double p1 = ws->logp[l], q1 = ws->logq[l];
        int det1 = ws->det[l];
        double lr1, lc1;
        double effr1 = rate * (bm != NULL ? bm[ws->orig[l]] : 1.0);
        if (!branch_logs(ws, l, effr1, &lr1, &lc1)) {
            gd_error_set(err, GD_E_RANGE,
                         "t1: child is older than its parent (chronology inversion)");
            return GD_E_RANGE;
        }

        double p2 = 0.0, q2 = -HUGE_VAL, lr2 = 0.0, lc2 = -HUGE_VAL;
        int det2 = -1;
        if (r != GD_NODE_NONE) {
            p2 = ws->logp[r];
            q2 = ws->logq[r];
            det2 = ws->det[r];
            double effr2 = rate * (bm != NULL ? bm[ws->orig[r]] : 1.0);
            if (!branch_logs(ws, r, effr2, &lr2, &lc2)) {
                gd_error_set(err, GD_E_RANGE,
                             "t1: child is older than its parent (chronology inversion)");
                return GD_E_RANGE;
            }
        }

        double pu, qu;
        if (det_u < 0) {
            /* Case i: u is not determined. */
            double pp  = p1 + p2;
            double sub = q1 + q2 + lr1 + lr2;
            pu = t1_logdiff(pp, sub, stability);

            double ta = q1 + lr1;
            double tb = q2 + lr2;
            ta = (ta == -HUGE_VAL) ? -HUGE_VAL
                                   : ta + t1_logdiff(p2, q2 + lr2, stability);
            tb = (tb == -HUGE_VAL) ? -HUGE_VAL
                                   : tb + t1_logdiff(p1, q1 + lr1, stability);
            qu = gd_logaddexp(ta, tb);
        } else {
            /* Case ii: u is determined, hence q_u = p_u. */
            double prob1 = (det1 == det_u) ? lr1 : lc1;
            double prob2 = (det2 == det_u) ? lr2 : lc2;
            if (det1 >= 0 && det2 >= 0) {                       /* iia */
                pu = p1 + prob1 + p2 + prob2;
            } else if (det1 < 0 && det2 < 0) {                  /* iib */
                pu = t1_logdiff(p1, q1 + lr1, stability)
                   + t1_logdiff(p2, q2 + lr2, stability);
            } else if (det1 >= 0) {                             /* iic */
                pu = p1 + prob1 + t1_logdiff(p2, q2 + lr2, stability);
            } else {                                            /* iic, mirrored */
                pu = p2 + prob2 + t1_logdiff(p1, q1 + lr1, stability);
            }
            qu = pu;
        }
        ws->logp[u] = pu;
        ws->logq[u] = qu;
    }

    *ops += ws->m;
    *out_logp = ws->logp[ws->root];
    return 0;
}

/* -------------------------------------------------------------- noise model */

/*
 * The noise component is a Chinese-restaurant process of concentration alpha
 * over the partition of the observed nodes:
 *   log P = B log a + sum_b log Gamma(n_b) - [log Gamma(a + n) - log Gamma(a)]
 * with B blocks of sizes n_b over n observed nodes. This is the K -> infinity
 * limit of the symmetric Dirichlet-multinomial (docs/model-spec.md §6) and is
 * normalized over exactly the observable space of the tree component.
 *
 * Hot path: block sizes read off the label-sorted keys of the extension.
 */
static double crp_logp_from_keys(const t1_key *keys, size_t nk, double alpha)
{
    if (nk == 0) return 0.0;
    double lp = 0.0;
    size_t nblocks = 0, s = 0;
    while (s < nk) {
        size_t e = s + 1;
        while (e < nk && keys[e].label == keys[s].label) e++;
        lp += lgamma((double)(e - s));
        nblocks++;
        s = e;
    }
    lp += (double)nblocks * log(alpha);
    lp -= lgamma(alpha + (double)nk) - lgamma(alpha);
    return lp;
}

double gd_t1_noise_loglik(const gd_cognate_data *data, double alpha)
{
    if (data == NULL || data->labels == NULL || !(alpha > 0.0)) return -HUGE_VAL;

    /* Standalone entry point (tests, diagnostics): counts blocks in place,
     * quadratically in the number of observed nodes. The likelihood itself
     * uses crp_logp_from_keys on the already-sorted keys. */
    size_t nobs = 0, nblocks = 0;
    double lp = 0.0;
    for (size_t i = 0; i < data->nnodes; i++) {
        if (data->labels[i] < 0) continue;
        nobs++;
        int first = 1;
        size_t count = 0;
        for (size_t j = 0; j < data->nnodes; j++) {
            if (data->labels[j] != data->labels[i]) continue;
            if (j < i) { first = 0; break; }
            count++;
        }
        if (first) {
            nblocks++;
            lp += lgamma((double)count);
        }
    }
    if (nobs == 0) return 0.0;
    lp += (double)nblocks * log(alpha);
    lp -= lgamma(alpha + (double)nobs) - lgamma(alpha);
    return lp;
}

/* ------------------------------------------------------------- likelihood */

static int check_options(const gd_t1_options *opt, gd_error *err)
{
    if (opt->epsilon < 0.0 || opt->epsilon >= 1.0 || !gd_is_finite(opt->epsilon)) {
        gd_error_set(err, GD_E_RANGE, "t1: contamination epsilon must be in [0,1)");
        return GD_E_RANGE;
    }
    if (!(opt->noise_alpha > 0.0) || !gd_is_finite(opt->noise_alpha)) {
        gd_error_set(err, GD_E_RANGE, "t1: noise alpha must be > 0");
        return GD_E_RANGE;
    }
    if (opt->incompatible != GD_INCOMPATIBLE_ERROR &&
        opt->incompatible != GD_INCOMPATIBLE_EXCLUDE &&
        opt->incompatible != GD_INCOMPATIBLE_CONTAMINATION) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1: unknown incompatible policy");
        return GD_E_INVALID_ARG;
    }
    return 0;
}

static int t1_concept_impl(gd_t1_ws *ws, double rate, const double *bm,
                           const gd_cognate_data *data,
                           const gd_t1_options *opt,
                           gd_t1_concept_result *res, gd_error *err)
{
    gd_t1_options defaults;
    gd_t1_concept_result local;

    if (ws == NULL || data == NULL || data->labels == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1: NULL argument");
        return GD_E_INVALID_ARG;
    }
    if (data->nnodes != ws->n) {
        gd_error_set(err, GD_E_INVALID_ARG,
                     "t1: data length does not match the tree");
        return GD_E_INVALID_ARG;
    }
    if (!gd_is_finite(rate) || rate < 0.0) {
        gd_error_set(err, GD_E_RANGE, "t1: rate must be finite and >= 0");
        return GD_E_RANGE;
    }
    if (opt == NULL) {
        gd_t1_options_default(&defaults);
        opt = &defaults;
    }
    int orc = check_options(opt, err);
    if (orc != 0) return orc;
    if (res == NULL) res = &local;
    result_init(res);

    size_t conflict = GD_NODE_NONE;
    int ca = -1, cb = -1;
    res->compatible = t1_extend(ws, data, &res->ops_extension,
                                &conflict, &ca, &cb);

    /* The extension label-sorts the observed nodes before it can fail, so the
     * noise term is available whether or not the data is extendable. */
    res->log_lik_noise = crp_logp_from_keys(ws->keys, ws->nkeys,
                                            opt->noise_alpha);

    if (!res->compatible) {
        const char *lbl = NULL;
        if (conflict != GD_NODE_NONE && conflict < ws->n) {
            lbl = ws->tree->nodes[conflict].label;
        }
        if (opt->incompatible == GD_INCOMPATIBLE_ERROR) {
            return gd_error_setf(err, GD_E_INCOMPATIBLE,
                     "t1: concept data is not extendable — node %s%s%s (index %zu) "
                     "is determined to both cognate class %d and class %d",
                     lbl ? "'" : "", lbl ? lbl : "(unlabelled)", lbl ? "'" : "",
                     conflict, ca, cb);
        }
        if (opt->incompatible == GD_INCOMPATIBLE_EXCLUDE) {
            res->excluded = true;
            res->log_lik = 0.0;          /* contributes nothing to the total */
            return 0;                    /* dropped, not attributed to noise */
        }
        /* Contamination: the tree component is zero; the noise part carries
         * the concept — but only if it has any weight. */
        if (opt->epsilon <= 0.0) {
            return gd_error_setf(err, GD_E_INCOMPATIBLE,
                     "t1: concept data is not extendable (node index %zu determined "
                     "to classes %d and %d) and contamination epsilon is 0; set "
                     "epsilon > 0 or use incompatible=exclude",
                     conflict, ca, cb);
        }
        res->log_lik = log(opt->epsilon) + res->log_lik_noise;
        res->p_contaminated = 1.0;
        return 0;
    }

    int rc = t1_recursion(ws, rate, bm, &res->log_lik_tree,
                          &res->stability_warning, &res->ops_recursion, err);
    if (rc != 0) return rc;

    if (opt->epsilon <= 0.0) {
        res->log_lik = res->log_lik_tree;
        res->p_contaminated = 0.0;
        return 0;
    }
    double a = log1p(-opt->epsilon) + res->log_lik_tree;
    double b = log(opt->epsilon)    + res->log_lik_noise;
    res->log_lik = gd_logaddexp(a, b);
    res->p_contaminated = (res->log_lik == -HUGE_VAL) ? 0.0
                                                      : exp(b - res->log_lik);
    return 0;
}

int gd_t1_loglik_concept(gd_t1_ws *ws, double rate,
                         const gd_cognate_data *data,
                         const gd_t1_options *opt,
                         gd_t1_concept_result *res, gd_error *err)
{
    return t1_concept_impl(ws, rate, NULL, data, opt, res, err);
}

int gd_t1_loglik_concept_relaxed(gd_t1_ws *ws, double rate,
                                 const double *branch_mult,
                                 const gd_cognate_data *data,
                                 const gd_t1_options *opt,
                                 gd_t1_concept_result *res, gd_error *err)
{
    return t1_concept_impl(ws, rate, branch_mult, data, opt, res, err);
}

int gd_t1_loglik_dataset(gd_t1_ws *ws, const double *rates,
                         const gd_cognate_data *concepts, size_t nconcepts,
                         const gd_t1_options *opt,
                         gd_t1_concept_result *per_concept,
                         double *log_lik_total, gd_error *err)
{
    return gd_t1_loglik_dataset_relaxed(ws, rates, NULL, concepts, nconcepts,
                                        opt, per_concept, log_lik_total, err);
}

int gd_t1_loglik_dataset_relaxed(gd_t1_ws *ws, const double *rates,
                                 const double *branch_mult,
                                 const gd_cognate_data *concepts, size_t nconcepts,
                                 const gd_t1_options *opt,
                                 gd_t1_concept_result *per_concept,
                                 double *log_lik_total, gd_error *err)
{
    if (ws == NULL || rates == NULL || concepts == NULL || log_lik_total == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1_dataset: NULL argument");
        return GD_E_INVALID_ARG;
    }
    double total = 0.0;
    for (size_t c = 0; c < nconcepts; c++) {
        gd_t1_concept_result scratch;
        gd_t1_concept_result *res = (per_concept != NULL) ? &per_concept[c]
                                                          : &scratch;
        int rc = t1_concept_impl(ws, rates[c], branch_mult, &concepts[c],
                                 opt, res, err);
        if (rc != 0) return rc;
        total += res->log_lik;
    }
    *log_lik_total = total;
    return 0;
}

/* ------------------------------------------------------------- simulation */

int gd_t1_simulate(gd_rng *rng, const gd_tree *tree, double rate,
                   int *labels, gd_error *err)
{
    if (rng == NULL || tree == NULL || labels == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t1_simulate: NULL argument");
        return GD_E_INVALID_ARG;
    }
    if (!gd_is_finite(rate) || rate < 0.0) {
        gd_error_set(err, GD_E_RANGE, "t1_simulate: rate must be finite and >= 0");
        return GD_E_RANGE;
    }
    int next_class = 1;
    for (size_t i = 0; i < tree->nnodes; i++) {
        size_t node = tree->preorder[i];
        if (node == tree->root) {
            labels[node] = 0;
            continue;
        }
        size_t parent = tree->nodes[node].parent;
        if (parent == GD_NODE_NONE) {
            return gd_error_set(err, GD_E_INVALID_ARG,
                                "t1: non-root node without a parent");
        }
        double dt = tree->nodes[parent].age - tree->nodes[node].age;
        if (dt < 0.0) {
            if (dt < -GD_T1_AGE_TOL) {
                gd_error_set(err, GD_E_RANGE,
                             "t1_simulate: child is older than its parent");
                return GD_E_RANGE;
            }
            dt = 0.0;
        }
        double ret = exp(-rate * dt);
        labels[node] = (gd_rng_uniform_closed_open(rng) < ret)
                     ? labels[parent]
                     : next_class++;
    }
    return 0;
}

/* --------------------------------------------------- MCMC adapter ----------
 *
 * T1 seen as a gd_mcmc_model. This lived inside gd_mcmc.c until E3, which made
 * the sampler the only implementation of its own interface and forced
 * gd_mcmc.h to include this header. A likelihood owns its adapter now, and the
 * sampler knows nothing about which likelihoods exist.
 */

struct gd_mcmc_t1 {
    gd_t1_ws *ws;
    const gd_cognate_data *concepts;
    size_t nconcepts;
    size_t nnodes;
    gd_t1_options opt;
    double *rates;
    double *bm;
    gd_t1_concept_result *per;
};

static int t1_model_loglik(void *vctx, gd_tree *tree, const gd_mcmc_state *st,
                           double *pw, double *out, gd_error *err)
{
    gd_mcmc_t1 *m = (gd_mcmc_t1 *)vctx;
    (void)tree;
    for (size_t c = 0; c < m->nconcepts; c++) {
        double dev = (st->concept_logdev != NULL) ? st->concept_logdev[c] : 0.0;
        m->rates[c] = exp(st->rates[0] + dev);
    }
    const double *bm = NULL;
    if (st->branch_logmult != NULL) {
        for (size_t i = 0; i < m->nnodes; i++) m->bm[i] = exp(st->branch_logmult[i]);
        bm = m->bm;
    }
    gd_t1_options opt = m->opt;
    opt.epsilon = st->epsilon;
    gd_t1_concept_result *per = (pw != NULL) ? m->per : NULL;
    double tot;
    int rc = gd_t1_loglik_dataset_relaxed(m->ws, m->rates, bm, m->concepts,
                                          m->nconcepts, &opt, per, &tot, err);
    if (rc != 0) return rc;
    if (pw != NULL) for (size_t c = 0; c < m->nconcepts; c++) pw[c] = m->per[c].log_lik;
    *out = tot;
    return 0;
}

gd_mcmc_t1 *gd_mcmc_t1_create(gd_t1_ws *ws, const gd_cognate_data *concepts,
                              size_t nconcepts, const gd_t1_options *opt,
                              gd_error *err)
{
    if (ws == NULL || concepts == NULL || nconcepts == 0) {
        gd_error_set(err, GD_E_INVALID_ARG, "mcmc_t1: NULL/zero"); return NULL;
    }
    gd_mcmc_t1 *m = (gd_mcmc_t1 *)gd_calloc_checked(1, sizeof(*m), err);
    if (m == NULL) { gd_error_set(err, GD_E_NOMEM, "mcmc_t1"); return NULL; }
    m->ws = ws;
    m->concepts = concepts;
    m->nconcepts = nconcepts;
    m->nnodes = concepts[0].nnodes;
    if (opt != NULL) m->opt = *opt; else gd_t1_options_default(&m->opt);
    m->rates = (double *)gd_malloc_array_checked(nconcepts, sizeof(double), err);
    m->bm = (double *)gd_malloc_array_checked(m->nnodes ? m->nnodes : 1,
                                              sizeof(double), err);
    m->per = (gd_t1_concept_result *)gd_malloc_array_checked(nconcepts, sizeof(gd_t1_concept_result), err);
    if (m->rates == NULL || m->bm == NULL || m->per == NULL) {
        gd_error_set(err, GD_E_NOMEM, "mcmc_t1: scratch");
        free(m->rates); free(m->bm); free(m->per); free(m);
        return NULL;
    }
    return m;
}

void gd_mcmc_t1_destroy(gd_mcmc_t1 *m)
{
    if (m == NULL) return;
    free(m->rates); free(m->bm); free(m->per); free(m);
}

void gd_mcmc_t1_model(gd_mcmc_t1 *m, gd_mcmc_model *out)
{
    if (m == NULL || out == NULL) return;
    out->ctx = m;
    out->n_pointwise = m->nconcepts;
    out->loglik = t1_model_loglik;
}

/*
 * gd_model_binary.c — binary cognate-set likelihood (T2).
 *
 * Each cognate class is a binary present/absent trait evolving under the
 * asymmetric two-state CTMC of docs/model-spec.md §3.4:
 *
 *     Q = [ -lambda_gain   lambda_gain ]     states 0 = absent, 1 = present.
 *         [  lambda_loss  -lambda_loss ]
 *
 * The per-class likelihood is Felsenstein's pruning algorithm in log space,
 * iterative postorder (browser stack safety), computing nothing on the heap in
 * the inner loop. Felsenstein pruning combines a node's children directly, so
 * a polytomy needs no special handling — there is deliberately no ladder here,
 * unlike T1 (gd_model_novel_state.c).
 *
 * Classes present in no sampled language are unobservable, so the likelihood
 * is conditioned on observability (the "all-absent" ascertainment correction).
 * Partition rates default to meaning-weighting, not cognate-count weighting,
 * per King 2026 (docs/paper-evaluation.md). Cognate classes within a meaning
 * are modelled as independent binary characters, which they are not; the
 * dataset result flags this whenever it is actually relied upon.
 */
#include "gd_model_binary.h"
#include "gd_memory.h"
#include "gd_tree.h"
#include "gd_tree_build.h"
#include "gd_math.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Ages must be non-increasing parent-to-child; anything below this is a
 * genuine chronology inversion rather than rounding. */
#define GD_T2_AGE_TOL 1e-9

struct gd_t2_ws {
    gd_tree *tree;
    size_t   n;        /* nodes of the tree */
    size_t   root;
    double  *logL0;    /* n; partial log-likelihood, state 0 (absent) */
    double  *logL1;    /* n; partial log-likelihood, state 1 (present) */
    int     *obsbuf;   /* n; scratch observation vector (all-absent term) */
};

/* ------------------------------------------------------------------ options */

void gd_t2_options_default(gd_t2_options *opt)
{
    if (opt == NULL) return;
    opt->ascertainment = true;
    opt->root_prior    = GD_T2_ROOT_STATIONARY;
    opt->root_pi1      = 0.5;
    opt->weight        = GD_T2_WEIGHT_MEANING;
}

int gd_t2_weight_from_string(const char *s, gd_t2_partition_weight *out,
                             gd_error *err)
{
    if (out == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2: NULL weight output");
        return GD_E_INVALID_ARG;
    }
    if (s == NULL || s[0] == '\0') { *out = GD_T2_WEIGHT_MEANING; return 0; }
    if (strcmp(s, "meaning") == 0)  { *out = GD_T2_WEIGHT_MEANING; return 0; }
    if (strcmp(s, "cognate") == 0)  { *out = GD_T2_WEIGHT_COGNATE; return 0; }
    gd_error_set(err, GD_E_INVALID_ARG,
                 "t2: partition weight must be meaning or cognate");
    return GD_E_INVALID_ARG;
}

/* --------------------------------------------------------- transition matrix */

int gd_t2_transition_logs(double lambda_gain, double lambda_loss, double t,
                          double logP[2][2], gd_error *err)
{
    if (!gd_is_finite(lambda_gain) || !gd_is_finite(lambda_loss) ||
        lambda_gain <= 0.0 || lambda_loss <= 0.0) {
        gd_error_set(err, GD_E_RANGE, "t2: gain and loss rates must be > 0");
        return GD_E_RANGE;
    }
    if (!gd_is_finite(t) || t < 0.0) {
        gd_error_set(err, GD_E_RANGE, "t2: branch duration must be finite and >= 0");
        return GD_E_RANGE;
    }
    double mu = lambda_gain + lambda_loss;
    double log_mu = log(mu);
    /* log(1 - exp(-mu t)) is the shared factor of both off-diagonal entries. */
    double log1me = gd_log1mexp(-mu * t);   /* -inf when t == 0 (identity) */
    /* P(0->1) = (lambda_gain/mu)(1 - e^{-mu t}); P(0->0) = 1 - P(0->1). */
    double logP01 = log(lambda_gain) - log_mu + log1me;
    double logP10 = log(lambda_loss) - log_mu + log1me;
    logP[0][1] = logP01;
    logP[0][0] = gd_log1mexp(logP01);       /* log(1 - P(0->1)) */
    logP[1][0] = logP10;
    logP[1][1] = gd_log1mexp(logP10);       /* log(1 - P(1->0)) */
    return 0;
}

/* --------------------------------------------------------------- workspace */

gd_t2_ws *gd_t2_ws_create(gd_tree *tree, gd_error *err)
{
    if (tree == NULL || tree->nnodes == 0 || tree->nodes == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2: NULL or empty tree");
        return NULL;
    }
    if (tree->postorder == NULL &&
        gd_tree_build_traversals(tree, err) != 0) {
        return NULL;
    }
    gd_t2_ws *ws = (gd_t2_ws *)gd_calloc_checked(1, sizeof(*ws), err);
    if (ws == NULL) {
        gd_error_set(err, GD_E_NOMEM, "t2: workspace alloc");
        return NULL;
    }
    ws->tree = tree;
    ws->n    = tree->nnodes;
    ws->root = tree->root;
    ws->logL0  = (double *)gd_malloc_array_checked(ws->n, sizeof(double), err);
    ws->logL1  = (double *)gd_malloc_array_checked(ws->n, sizeof(double), err);
    ws->obsbuf = (int *)gd_malloc_array_checked(ws->n, sizeof(int), err);
    if (ws->logL0 == NULL || ws->logL1 == NULL || ws->obsbuf == NULL) {
        free(ws->logL0); free(ws->logL1); free(ws->obsbuf); free(ws);
        gd_error_set(err, GD_E_NOMEM, "t2: workspace arrays");
        return NULL;
    }
    return ws;
}

void gd_t2_ws_destroy(gd_t2_ws *ws)
{
    if (ws == NULL) return;
    free(ws->logL0); free(ws->logL1); free(ws->obsbuf);
    free(ws);
}

/* ---------------------------------------------------------------- pruning */

/* Root prior log-probabilities of states 0 and 1. */
static void root_logprior(const gd_t2_options *opt, double lg, double ll,
                          double *logpi0, double *logpi1)
{
    if (opt->root_prior == GD_T2_ROOT_FIXED) {
        *logpi1 = log(opt->root_pi1);
        *logpi0 = log1p(-opt->root_pi1);
    } else {
        double mu = lg + ll;
        *logpi0 = log(ll) - log(mu);   /* stationary P(absent)  = ll/mu */
        *logpi1 = log(lg) - log(mu);   /* stationary P(present) = lg/mu */
    }
}

/*
 * Felsenstein pruning of one observation vector. obs[v] in {0 absent, >=1
 * present, <0 missing}. Fills the log-likelihood at the root into *out_ll.
 */
static int t2_prune(gd_t2_ws *ws, double lg, double ll, const int *obs,
                    const gd_t2_options *opt, double *out_ll, size_t *ops,
                    gd_error *err)
{
    const gd_tree *t = ws->tree;
    for (size_t i = 0; i < ws->n; i++) {
        size_t v = t->postorder[i];
        const gd_node *node = &t->nodes[v];
        double m0 = 0.0, m1 = 0.0;     /* log 1: empty product over children */
        for (size_t j = 0; j < node->nchildren; j++) {
            size_t c = node->children[j];
            double dt = node->age - t->nodes[c].age;
            if (dt < 0.0) {
                if (dt < -GD_T2_AGE_TOL) {
                    gd_error_set(err, GD_E_RANGE,
                                 "t2: child is older than its parent "
                                 "(chronology inversion)");
                    return GD_E_RANGE;
                }
                dt = 0.0;
            }
            double logP[2][2];
            if (gd_t2_transition_logs(lg, ll, dt, logP, err) != 0) return GD_E_RANGE;
            m0 += gd_logaddexp(logP[0][0] + ws->logL0[c],
                               logP[0][1] + ws->logL1[c]);
            m1 += gd_logaddexp(logP[1][0] + ws->logL0[c],
                               logP[1][1] + ws->logL1[c]);
        }
        int o = obs[v];
        if (o == 0) {                  /* observed absent */
            ws->logL0[v] = m0;
            ws->logL1[v] = -HUGE_VAL;
        } else if (o > 0) {            /* observed present */
            ws->logL0[v] = -HUGE_VAL;
            ws->logL1[v] = m1;
        } else {                       /* missing: marginalized */
            ws->logL0[v] = m0;
            ws->logL1[v] = m1;
        }
    }
    double logpi0, logpi1;
    root_logprior(opt, lg, ll, &logpi0, &logpi1);
    *out_ll = gd_logaddexp(logpi0 + ws->logL0[ws->root],
                           logpi1 + ws->logL1[ws->root]);
    *ops += ws->n;
    return 0;
}

/* ------------------------------------------------------------- likelihood */

int gd_t2_loglik_class(gd_t2_ws *ws, double lambda_gain, double lambda_loss,
                       const gd_binary_class *cls, const gd_t2_options *opt,
                       gd_t2_class_result *res, gd_error *err)
{
    gd_t2_options defaults;
    gd_t2_class_result local;

    if (ws == NULL || cls == NULL || cls->present == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2: NULL argument");
        return GD_E_INVALID_ARG;
    }
    if (cls->nnodes != ws->n) {
        gd_error_set(err, GD_E_INVALID_ARG,
                     "t2: class length does not match the tree");
        return GD_E_INVALID_ARG;
    }
    if (!gd_is_finite(lambda_gain) || !gd_is_finite(lambda_loss) ||
        lambda_gain <= 0.0 || lambda_loss <= 0.0) {
        gd_error_set(err, GD_E_RANGE, "t2: gain and loss rates must be > 0");
        return GD_E_RANGE;
    }
    if (opt == NULL) { gd_t2_options_default(&defaults); opt = &defaults; }
    if (opt->root_prior == GD_T2_ROOT_FIXED &&
        (!(opt->root_pi1 > 0.0) || !(opt->root_pi1 < 1.0))) {
        gd_error_set(err, GD_E_RANGE, "t2: fixed root_pi1 must be in (0,1)");
        return GD_E_RANGE;
    }
    if (res == NULL) res = &local;
    res->log_lik = res->log_lik_raw = res->log_p_all_absent = -HUGE_VAL;
    res->ascertained = false;
    res->all_absent = true;
    res->ops = 0;
    res->stability_warning = false;

    /* Is the class observed present anywhere? (degenerate if not.) */
    for (size_t i = 0; i < ws->n; i++) {
        if (cls->present[i] > 0) { res->all_absent = false; break; }
    }

    int rc = t2_prune(ws, lambda_gain, lambda_loss, cls->present, opt,
                      &res->log_lik_raw, &res->ops, err);
    if (rc != 0) return rc;

    /* All-absent term: force every observed node to absent, keep missing free.
     * P(no observed node carries the class) — the unobservable pattern. */
    for (size_t i = 0; i < ws->n; i++) {
        ws->obsbuf[i] = (cls->present[i] < 0) ? -1 : 0;
    }
    rc = t2_prune(ws, lambda_gain, lambda_loss, ws->obsbuf, opt,
                  &res->log_p_all_absent, &res->ops, err);
    if (rc != 0) return rc;

    if (!opt->ascertainment) {
        res->log_lik = res->log_lik_raw;
        return 0;
    }
    /* Condition on observability: divide by 1 - P(all-absent). */
    res->ascertained = true;
    double log_denom = gd_log1mexp(res->log_p_all_absent);
    if (!gd_is_finite(log_denom)) {
        /* P(all-absent) ~ 1: essentially no observable mass. The correction is
         * numerically unreliable; fall back to the raw likelihood and flag. */
        res->stability_warning = true;
        res->log_lik = res->log_lik_raw;
        return 0;
    }
    res->log_lik = res->log_lik_raw - log_denom;
    return 0;
}

int gd_t2_partition_weights(const size_t *concept_ids, size_t nclasses,
                            size_t nconcepts, gd_t2_partition_weight mode,
                            double *weights, gd_error *err)
{
    if (weights == NULL || (nclasses > 0 && concept_ids == NULL)) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2_weights: NULL argument");
        return GD_E_INVALID_ARG;
    }
    for (size_t c = 0; c < nconcepts; c++) weights[c] = 0.0;
    for (size_t k = 0; k < nclasses; k++) {
        size_t c = concept_ids[k];
        if (c >= nconcepts) {
            gd_error_set(err, GD_E_INVALID_ARG,
                         "t2_weights: concept id out of range");
            return GD_E_INVALID_ARG;
        }
        weights[c] += 1.0;   /* class count */
    }
    if (mode == GD_T2_WEIGHT_MEANING) {
        /* Every meaning that has at least one class counts exactly once. */
        for (size_t c = 0; c < nconcepts; c++) {
            weights[c] = (weights[c] > 0.0) ? 1.0 : 0.0;
        }
    }
    return 0;
}

int gd_t2_loglik_dataset(gd_t2_ws *ws,
                         double lambda_gain, double lambda_loss,
                         const gd_binary_class *classes,
                         const size_t *concept_ids, size_t nclasses,
                         const double *concept_rates, size_t nconcepts,
                         const gd_t2_options *opt,
                         gd_t2_class_result *per_class,
                         gd_t2_dataset_result *summary,
                         double *log_lik_total, gd_error *err)
{
    if (ws == NULL || classes == NULL || log_lik_total == NULL ||
        (nclasses > 0 && concept_ids == NULL)) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2_dataset: NULL argument");
        return GD_E_INVALID_ARG;
    }
    size_t *counts = NULL;
    if (nconcepts > 0) {
        counts = (size_t *)gd_calloc_checked(nconcepts, sizeof(size_t), err);
        if (counts == NULL) {
            gd_error_set(err, GD_E_NOMEM, "t2_dataset: counts");
            return GD_E_NOMEM;
        }
    }
    double total = 0.0;
    int rc = 0;
    for (size_t k = 0; k < nclasses; k++) {
        size_t c = concept_ids[k];
        if (c >= nconcepts) {
            gd_error_set(err, GD_E_INVALID_ARG,
                         "t2_dataset: concept id out of range");
            rc = GD_E_INVALID_ARG;
            break;
        }
        counts[c]++;
        double scale = (concept_rates != NULL) ? concept_rates[c] : 1.0;
        if (!gd_is_finite(scale) || scale <= 0.0) {
            gd_error_set(err, GD_E_RANGE, "t2_dataset: concept rate must be > 0");
            rc = GD_E_RANGE;
            break;
        }
        gd_t2_class_result scratch;
        gd_t2_class_result *cr = (per_class != NULL) ? &per_class[k] : &scratch;
        rc = gd_t2_loglik_class(ws, lambda_gain * scale, lambda_loss * scale,
                                &classes[k], opt, cr, err);
        if (rc != 0) break;
        total += cr->log_lik;
    }
    if (summary != NULL) {
        summary->n_classes = nclasses;
        summary->n_concepts = nconcepts;
        size_t maxc = 0;
        for (size_t c = 0; c < nconcepts; c++) {
            if (counts != NULL && counts[c] > maxc) maxc = counts[c];
        }
        summary->max_classes_per_concept = maxc;
        summary->within_meaning_dependence_warning = (maxc > 1);
        summary->cross_encoding_incomparable = true;
    }
    free(counts);
    if (rc != 0) return rc;
    *log_lik_total = total;
    return 0;
}

/* ------------------------------------------------------------- simulation */

int gd_t2_simulate_class(gd_rng *rng, const gd_tree *tree,
                         double lambda_gain, double lambda_loss,
                         const gd_t2_options *opt, int *present,
                         gd_error *err)
{
    gd_t2_options defaults;
    if (rng == NULL || tree == NULL || present == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2_simulate: NULL argument");
        return GD_E_INVALID_ARG;
    }
    if (!gd_is_finite(lambda_gain) || !gd_is_finite(lambda_loss) ||
        lambda_gain <= 0.0 || lambda_loss <= 0.0) {
        gd_error_set(err, GD_E_RANGE, "t2_simulate: rates must be > 0");
        return GD_E_RANGE;
    }
    if (opt == NULL) { gd_t2_options_default(&defaults); opt = &defaults; }
    if (tree->preorder == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "t2_simulate: tree has no traversal");
        return GD_E_INVALID_ARG;
    }
    double mu = lambda_gain + lambda_loss;
    double pi1 = (opt->root_prior == GD_T2_ROOT_FIXED) ? opt->root_pi1
                                                       : lambda_gain / mu;
    for (size_t i = 0; i < tree->nnodes; i++) {
        size_t v = tree->preorder[i];
        if (v == tree->root) {
            present[v] = (gd_rng_uniform_closed_open(rng) < pi1) ? 1 : 0;
            continue;
        }
        size_t p = tree->nodes[v].parent;
        if (p == GD_NODE_NONE) {
            return gd_error_set(err, GD_E_INVALID_ARG,
                                "t2: non-root node without a parent");
        }
        double dt = tree->nodes[p].age - tree->nodes[v].age;
        if (dt < 0.0) {
            if (dt < -GD_T2_AGE_TOL) {
                gd_error_set(err, GD_E_RANGE,
                             "t2_simulate: child is older than its parent");
                return GD_E_RANGE;
            }
            dt = 0.0;
        }
        double e = exp(-mu * dt);
        /* P(child present | parent state): P(0->1)=pg(1-e), P(1->1)=1-pl(1-e). */
        double pg = lambda_gain / mu;
        double pl = lambda_loss / mu;
        double p_child_present = (present[p] == 1)
            ? (1.0 - pl * (1.0 - e))
            : (pg * (1.0 - e));
        present[v] = (gd_rng_uniform_closed_open(rng) < p_child_present) ? 1 : 0;
    }
    return 0;
}


/* --------------------------------------------------- MCMC adapter ---------- */

struct gd_mcmc_t2 {
    gd_t2_ws               *ws;
    const gd_binary_class  *classes;
    const size_t           *concept_ids;
    size_t                  nclasses, nconcepts;
    gd_t2_options           opt;
    gd_t2_class_result     *per;
};

static int t2_model_loglik(void *vctx, gd_tree *tree, const gd_mcmc_state *st,
                           double *pw, double *out, gd_error *err)
{
    gd_mcmc_t2 *m = (gd_mcmc_t2 *)vctx;
    (void)tree;   /* ages are read off the tree inside the workspace */
    if (st->n_rates < 2) {
        return gd_error_set(err, GD_E_INVALID_ARG,
                            "mcmc_t2: needs a gain and a loss rate");
    }
    double lambda_gain = exp(st->rates[0]);
    double lambda_loss = exp(st->rates[1]);
    gd_t2_class_result *per = (pw != NULL) ? m->per : NULL;
    double tot;
    int rc = gd_t2_loglik_dataset(m->ws, lambda_gain, lambda_loss,
                                  m->classes, m->concept_ids, m->nclasses,
                                  NULL, m->nconcepts, &m->opt,
                                  per, NULL, &tot, err);
    if (rc != GD_E_OK) return rc;
    if (pw != NULL)
        for (size_t c = 0; c < m->nclasses; c++) pw[c] = m->per[c].log_lik;
    *out = tot;
    return GD_E_OK;
}

gd_mcmc_t2 *gd_mcmc_t2_create(gd_t2_ws *ws,
                              const gd_binary_class *classes,
                              const size_t *concept_ids, size_t nclasses,
                              size_t nconcepts,
                              const gd_t2_options *opt, gd_error *err)
{
    if (ws == NULL || classes == NULL || concept_ids == NULL ||
        nclasses == 0 || nconcepts == 0) {
        gd_error_set(err, GD_E_INVALID_ARG, "mcmc_t2: NULL/zero");
        return NULL;
    }
    gd_mcmc_t2 *m = (gd_mcmc_t2 *)gd_calloc_checked(1, sizeof *m, err);
    if (m == NULL) return NULL;
    m->ws = ws;
    m->classes = classes;
    m->concept_ids = concept_ids;
    m->nclasses = nclasses;
    m->nconcepts = nconcepts;
    if (opt != NULL) m->opt = *opt; else gd_t2_options_default(&m->opt);
    m->per = (gd_t2_class_result *)gd_malloc_array_checked(
                    nclasses, sizeof(gd_t2_class_result), err);
    if (m->per == NULL) { gd_free(m); return NULL; }
    return m;
}

void gd_mcmc_t2_destroy(gd_mcmc_t2 *m)
{
    if (m == NULL) return;
    gd_free(m->per);
    gd_free(m);
}

void gd_mcmc_t2_model(gd_mcmc_t2 *m, gd_mcmc_model *out)
{
    if (m == NULL || out == NULL) return;
    out->ctx = m;
    out->n_pointwise = m->nclasses;
    out->loglik = t2_model_loglik;
}

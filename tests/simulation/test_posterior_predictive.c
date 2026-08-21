/*
 * test_posterior_predictive.c — adequacy checks and model comparison (M11).
 *
 * SBC (tests/simulation/test_mcmc_*.c) establishes *implementation validity*:
 * the sampler targets the right posterior. It does NOT establish *empirical
 * adequacy*: a validly-implemented model can still be wrong for the data
 * (King 2026). This suite separates the two.
 *
 *   1. Posterior-predictive calibration — for a well-specified fit, an observed
 *      discrepancy statistic is typical of data simulated from the posterior
 *      (a moderate posterior-predictive p-value).
 *   2. A documented failure regime — a homogeneous-rate model fit to strongly
 *      concept-rate-heterogeneous data is validly implemented but inadequate;
 *      the same posterior-predictive check flags it (an extreme p-value).
 *   3. Model comparison in a shared observable space — on strictly-clocked
 *      data, the more flexible relaxed clock does NOT beat the strict clock on
 *      held-out cell states. "Most flexible" is not "best".
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_mcmc.h"
#include "gd_model_novel_state.h"
#include "gd_tree.h"
#include "gd_calibration.h"
#include "gd_rng.h"
#include "gd_math.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- fixtures */

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

static gd_calibration_set *bounds_all_fixed(const gd_tree *t, gd_error *err)
{
    gd_calibration_set *cs = gd_calset_create(0, err);
    if (cs == NULL) return NULL;
    size_t n = t->nnodes;
    cs->nnodes = n;
    cs->node_lo = (double *)malloc(n * sizeof(double));
    cs->node_hi = (double *)malloc(n * sizeof(double));
    if (!cs->node_lo || !cs->node_hi) { gd_calset_destroy(cs); return NULL; }
    for (size_t i = 0; i < n; i++) { cs->node_lo[i] = cs->node_hi[i] = t->nodes[i].age; }
    cs->feasible = true;
    return cs;
}

/* ------------------------------------------------- discrepancy statistics */

/* Number of distinct cognate classes among the observed tips of a concept. */
static int n_classes(const int *labels, const size_t *tips, size_t ntip)
{
    int seen[64]; size_t ns = 0;
    for (size_t i = 0; i < ntip; i++) {
        int lab = labels[tips[i]];
        if (lab < 0) continue;
        int dup = 0;
        for (size_t j = 0; j < ns; j++) if (seen[j] == lab) { dup = 1; break; }
        if (!dup && ns < 64) seen[ns++] = lab;
    }
    return (int)ns;
}

/* Variance across concepts of the per-concept distinct-class count. A model
 * with a single rate predicts a narrow distribution of class counts; genuine
 * concept-rate heterogeneity inflates this variance. */
static double stat_class_count_variance(const int *store, size_t ncon, size_t n,
                                        const size_t *tips, size_t ntip)
{
    double mean = 0.0;
    double *cc = (double *)malloc(ncon * sizeof(double));
    for (size_t c = 0; c < ncon; c++) {
        cc[c] = (double)n_classes(store + c * n, tips, ntip);
        mean += cc[c];
    }
    mean /= (double)ncon;
    double var = 0.0;
    for (size_t c = 0; c < ncon; c++) { double d = cc[c] - mean; var += d * d; }
    var /= (double)ncon;
    free(cc);
    return var;
}

/* ------------------------------------------------- posterior-predictive p */

/* Fit a single global rate (ages fixed) and return the posterior-predictive
 * p-value P(T(y_rep) >= T(y_obs)) for the class-count-variance statistic. */
static double pp_pvalue(gd_tree *t, gd_calibration_set *cs,
                        const int *store, size_t ncon, const size_t *tips,
                        size_t ntip, double rate_logmean, uint64_t seed)
{
    gd_error err = GD_ERROR_INIT;
    size_t n = t->nnodes;
    gd_cognate_data *concepts = (gd_cognate_data *)malloc(ncon * sizeof(gd_cognate_data));
    for (size_t c = 0; c < ncon; c++) { concepts[c].labels = store + c * n; concepts[c].nnodes = n; }
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt); opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, ncon, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);

    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.seed = seed; cfg.n_warmup = 400; cfg.n_samples = 200; cfg.thin = 1;
    cfg.rate_logmean = rate_logmean; cfg.rate_logsd = 1.0;
    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);

    double t_obs = stat_class_count_variance(store, ncon, n, tips, ntip);

    size_t total = res->n_chains * res->n_draws, np = res->n_param;
    size_t p_rate = 0;   /* single sampled param: the log-rate */
    gd_rng *rng = gd_rng_create(seed ^ 0x9e3779b9u, 42u, &err);
    int *rep = (int *)malloc(ncon * n * sizeof(int));
    int *sim = (int *)malloc(n * sizeof(int));
    size_t ge = 0;
    for (size_t s = 0; s < total; s++) {
        double r = exp(res->samples[s * np + p_rate]);
        for (size_t c = 0; c < ncon; c++) {
            gd_t1_simulate(rng, t, r, sim, &err);
            int *row = rep + c * n;
            for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
        }
        double t_rep = stat_class_count_variance(rep, ncon, n, tips, ntip);
        if (t_rep >= t_obs) ge++;
    }
    double p = (double)ge / (double)total;

    free(rep); free(sim); gd_rng_destroy(rng);
    gd_mcmc_result_destroy(res);
    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws); free(concepts);
    return p;
}

/* Collect tip node indices. */
static size_t collect_tips(const gd_tree *t, size_t *tips)
{
    size_t ntip = 0;
    for (size_t i = 0; i < t->nnodes; i++) if (t->nodes[i].is_leaf) tips[ntip++] = i;
    return ntip;
}

TEST(test_pp_wellspecified_calibrated)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_balanced(3, &err);        /* 8 tips */
    CHECK(t != NULL);
    size_t n = t->nnodes, tips[8]; size_t ntip = collect_tips(t, tips);
    gd_calibration_set *cs = bounds_all_fixed(t, &err);
    gd_rng *rng = gd_rng_create(12345u, 1u, &err);

    const size_t NCON = 90;
    const double rate = 0.0005;
    int *store = (int *)malloc(NCON * n * sizeof(int));
    int *sim = (int *)malloc(n * sizeof(int));
    for (size_t c = 0; c < NCON; c++) {         /* homogeneous rate: well-specified */
        gd_t1_simulate(rng, t, rate, sim, &err);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
    }
    free(sim);
    double p = pp_pvalue(t, cs, store, NCON, tips, ntip, log(rate), 700u);
    printf("    well-specified PP p-value = %.3f\n", p);
    CHECK(p > 0.05 && p < 0.95);                /* not extreme: adequate */

    free(store); gd_calset_destroy(cs); gd_tree_destroy(t); gd_rng_destroy(rng);
}

TEST(test_pp_detects_rate_heterogeneity)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_balanced(3, &err);
    CHECK(t != NULL);
    size_t n = t->nnodes, tips[8]; size_t ntip = collect_tips(t, tips);
    gd_calibration_set *cs = bounds_all_fixed(t, &err);
    gd_rng *rng = gd_rng_create(24680u, 2u, &err);

    const size_t NCON = 90;
    const double base = 0.0005;
    int *store = (int *)malloc(NCON * n * sizeof(int));
    int *sim = (int *)malloc(n * sizeof(int));
    for (size_t c = 0; c < NCON; c++) {
        /* Strong concept-rate heterogeneity: half very slow, half very fast.
         * A single-rate model cannot reproduce the resulting class-count spread. */
        double r = (c % 2 == 0) ? base * 0.15 : base * 6.0;
        gd_t1_simulate(rng, t, r, sim, &err);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
    }
    free(sim);
    /* Fit a single global rate; the PP check should flag the misspecification. */
    double p = pp_pvalue(t, cs, store, NCON, tips, ntip, log(base), 900u);
    printf("    heterogeneous-data PP p-value = %.3f (misspecified)\n", p);
    CHECK(p < 0.05);                            /* extreme: inadequacy detected */

    free(store); gd_calset_destroy(cs); gd_tree_destroy(t); gd_rng_destroy(rng);
}

/* --------------------------------- model comparison in shared observable space */

/* Held-out predictive log-score sum over test concepts, given posterior draws
 * of the rate and (optionally) per-branch log-multipliers. */
static double heldout_lppd(gd_t1_ws *ws, gd_tree *t,
                           const gd_cognate_data *test, size_t ntest,
                           const gd_mcmc_result *res, int relaxed)
{
    gd_error err = GD_ERROR_INIT;
    size_t total = res->n_chains * res->n_draws, np = res->n_param, n = t->nnodes;
    size_t p_rate = (size_t)-1;
    for (size_t p = 0; p < np; p++) if (res->param_type[p] == GD_MCMC_P_LOGRATE) p_rate = p;
    double *bm = (double *)malloc(n * sizeof(double));
    gd_t1_options opt; gd_t1_options_default(&opt); opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;

    double total_lppd = 0.0;
    double *perdraw = (double *)malloc(total * sizeof(double));
    for (size_t k = 0; k < ntest; k++) {
        for (size_t s = 0; s < total; s++) {
            double r = exp(res->samples[s * np + p_rate]);
            gd_t1_concept_result cr;
            if (relaxed) {
                for (size_t i = 0; i < n; i++) bm[i] = 1.0;
                for (size_t p = 0; p < np; p++)
                    if (res->param_type[p] == GD_MCMC_P_BRANCH_LOGMULT)
                        bm[res->param_node[p]] = exp(res->samples[s * np + p]);
                gd_t1_loglik_concept_relaxed(ws, r, bm, &test[k], &opt, &cr, &err);
            } else {
                gd_t1_loglik_concept(ws, r, &test[k], &opt, &cr, &err);
            }
            perdraw[s] = cr.log_lik_tree;
        }
        /* lppd_k = log( mean_s exp(loglik_s) ) = logsumexp - log(S). */
        total_lppd += gd_logsumexp(perdraw, total) - log((double)total);
    }
    free(perdraw); free(bm);
    return total_lppd;
}

TEST(test_model_comparison_shared_observable)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = tree_balanced(3, &err);
    CHECK(t != NULL);
    size_t n = t->nnodes;
    gd_calibration_set *cs = bounds_all_fixed(t, &err);
    gd_rng *rng = gd_rng_create(55555u, 3u, &err);

    /* Strictly-clocked truth. Split concepts into train / test. */
    const size_t NTRAIN = 120, NTEST = 60, NCON = NTRAIN + NTEST;
    const double rate = 0.0005;
    int *store = (int *)malloc(NCON * n * sizeof(int));
    int *sim = (int *)malloc(n * sizeof(int));
    gd_cognate_data *concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
    for (size_t c = 0; c < NCON; c++) {
        gd_t1_simulate(rng, t, rate, sim, &err);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
        concepts[c].labels = row; concepts[c].nnodes = n;
    }
    free(sim);
    const gd_cognate_data *train = concepts, *test = concepts + NTRAIN;

    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt); opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, train, NTRAIN, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);

    /* Fit the strict clock. */
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.seed = 11u; cfg.n_warmup = 500; cfg.n_samples = 300; cfg.thin = 1;
    cfg.rate_logmean = log(rate); cfg.rate_logsd = 1.0;
    gd_mcmc_result *res_strict = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);
    CHECK(res_strict != NULL);

    /* Fit the relaxed clock (branch multipliers + sigma_b). */
    cfg.sample_branch = true;
    gd_mcmc_result *res_relaxed = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);
    CHECK(res_relaxed != NULL);

    double lppd_strict = heldout_lppd(ws, t, test, NTEST, res_strict, 0);
    double lppd_relaxed = heldout_lppd(ws, t, test, NTEST, res_relaxed, 1);
    printf("    held-out lppd: strict=%.2f  relaxed=%.2f  (delta=%.2f)\n",
           lppd_strict, lppd_relaxed, lppd_relaxed - lppd_strict);
    /* On strictly-clocked data the extra flexibility does not help out-of-sample:
     * the relaxed clock must not beat the strict clock by any real margin. */
    CHECK(lppd_relaxed <= lppd_strict + 1.0);

    gd_mcmc_result_destroy(res_strict);
    gd_mcmc_result_destroy(res_relaxed);
    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    free(store); free(concepts);
    gd_calset_destroy(cs); gd_tree_destroy(t); gd_rng_destroy(rng);
}

static void run_tests(void)
{
    RUN(test_pp_wellspecified_calibrated);
    RUN(test_pp_detects_rate_heterogeneity);
    RUN(test_model_comparison_shared_observable);
}

TEST_MAIN()

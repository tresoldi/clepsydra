/*
 * test_mcmc_relaxed.c — relaxed clock + robust observations (M9).
 *
 * The M9 acceptance gate (PLAN.md): simulated branch-rate recovery where
 * identifiable; a warning when branch rates and node ages are confounded;
 * robustness tests with injected loans / miscoding.
 *
 * Setup: the quartet ((a,b)r1,(c,d)r2)root under T1. Branch multipliers m_b
 * scale each edge's rate (log m_b ~ N(-sigma_b^2/2, sigma_b^2)); the sampler
 * recovers them when the node ages are fixed (identifiable) and flags the
 * confounding when the ages are free. Contamination epsilon absorbs injected
 * loan/miscoded concepts.
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_mcmc.h"
#include "gd_model_novel_state.h"
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_calibration.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- fixtures */

static gd_tree *make_quartet(void)
{
    const char *nwk = "((a,b)r1,(c,d)r2)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    for (size_t i = 0; i < t->nnodes; i++) t->nodes[i].age = 0.0;
    t->nodes[gd_tree_find_label(t, "root")].age = 3000.0;
    t->nodes[gd_tree_find_label(t, "r1")].age = 2000.0;
    t->nodes[gd_tree_find_label(t, "r2")].age = 1500.0;
    return t;
}

/* Bounds fixing every node at its current age (a fully fixed chronology). */
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

/* Bounds: root fixed at R, internal non-root free in (0,R), tips at 0. */
static gd_calibration_set *bounds_free_internal(const gd_tree *t, double R, gd_error *err)
{
    gd_calibration_set *cs = gd_calset_create(0, err);
    if (cs == NULL) return NULL;
    size_t n = t->nnodes;
    cs->nnodes = n;
    cs->node_lo = (double *)malloc(n * sizeof(double));
    cs->node_hi = (double *)malloc(n * sizeof(double));
    if (!cs->node_lo || !cs->node_hi) { gd_calset_destroy(cs); return NULL; }
    for (size_t i = 0; i < n; i++) {
        if (t->nodes[i].is_leaf) { cs->node_lo[i] = cs->node_hi[i] = 0.0; }
        else if (i == t->root)   { cs->node_lo[i] = cs->node_hi[i] = R; }
        else                     { cs->node_lo[i] = 0.0; cs->node_hi[i] = R; }
    }
    cs->feasible = true;
    return cs;
}

/* Forward simulation with per-branch multipliers (bm indexed by child node). */
static void simulate_relaxed(gd_rng *rng, const gd_tree *t, double rate,
                             const double *bm, int *labels)
{
    int next = 1;
    for (size_t i = 0; i < t->nnodes; i++) {
        size_t v = t->preorder[i];
        if (v == t->root) { labels[v] = 0; continue; }
        size_t par = t->nodes[v].parent;
        double dt = t->nodes[par].age - t->nodes[v].age;
        if (dt < 0) dt = 0;
        double er = rate * (bm ? bm[v] : 1.0);
        double ret = exp(-er * dt);
        labels[v] = (gd_rng_uniform_closed_open(rng) < ret) ? labels[par] : next++;
    }
}

/* --------------------------------------------------- branch-rate recovery */

TEST(test_m9_branch_rate_recovery)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_quartet();
    CHECK(t != NULL);
    size_t n = t->nnodes;
    gd_rng *rng = gd_rng_create(9001u, 3u, &err);
    CHECK(rng != NULL);

    /* True per-branch multipliers (log scale, moderate spread). */
    double *true_logm = (double *)calloc(n, sizeof(double));
    double *bm = (double *)malloc(n * sizeof(double));
    for (size_t v = 0; v < n; v++) {
        true_logm[v] = (v == t->root) ? 0.0 : gd_rng_normal(rng, 0.0, 0.45);
        bm[v] = exp(true_logm[v]);
    }
    const double rate = 0.0005;
    const size_t NCON = 300;
    int *store = (int *)malloc(NCON * n * sizeof(int));
    gd_cognate_data *concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
    int *sim = (int *)malloc(n * sizeof(int));
    for (size_t c = 0; c < NCON; c++) {
        simulate_relaxed(rng, t, rate, bm, sim);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
        concepts[c].labels = row; concepts[c].nnodes = n;
    }
    free(sim);

    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt); opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, NCON, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);
    gd_calibration_set *cs = bounds_all_fixed(t, &err);   /* ages fixed => identifiable */
    CHECK(cs != NULL);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 314u;
    cfg.n_warmup = 800; cfg.n_samples = 800; cfg.thin = 1;
    cfg.rate_logmean = log(rate); cfg.rate_logsd = 1.0;
    cfg.sample_branch = true;

    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);
    CHECK(res != NULL);
    CHECK(res->rhat != NULL);

    /* Posterior-mean log m_b per branch, correlated against the truth. */
    size_t total = res->n_chains * res->n_draws, np = res->n_param;
    double post[16], truth[16]; size_t nb = 0;
    for (size_t p = 0; p < np; p++) {
        if (res->param_type[p] != GD_MCMC_P_BRANCH_LOGMULT) continue;
        double s = 0.0;
        for (size_t i = 0; i < total; i++) s += res->samples[i * np + p];
        post[nb] = s / (double)total;
        truth[nb] = true_logm[res->param_node[p]];
        nb++;
    }
    CHECK(nb == n - 1);   /* one multiplier per non-root branch */
    /* Pearson correlation of recovered vs. true branch rates. */
    double mp = 0, mt = 0;
    for (size_t i = 0; i < nb; i++) { mp += post[i]; mt += truth[i]; }
    mp /= (double)nb; mt /= (double)nb;
    double spp = 0, stt = 0, spt = 0;
    for (size_t i = 0; i < nb; i++) {
        double dp = post[i] - mp, dt = truth[i] - mt;
        spp += dp * dp; stt += dt * dt; spt += dp * dt;
    }
    double corr = spt / sqrt(spp * stt);
    printf("    branch-rate recovery corr = %.3f (nb=%zu)\n", corr, nb);
    CHECK(corr > 0.5);
    /* Ages fixed => no age variance => not flagged confounded. */
    CHECK(res->branch_rate_confounded == 0);

    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);
    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    free(true_logm); free(bm); free(store); free(concepts);
    gd_tree_destroy(t); gd_rng_destroy(rng);
}

/* ------------------------------------------------------ confounding warning */

TEST(test_m9_confounding_when_ages_free)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_quartet();
    CHECK(t != NULL);
    size_t n = t->nnodes;
    const double R = 3000.0;
    gd_rng *rng = gd_rng_create(9002u, 7u, &err);
    CHECK(rng != NULL);

    /* Simulate strict data; the confounding is a property of the posterior
     * geometry when both branch rates and ages are free, not of the data. */
    const double rate = 0.0006;
    const size_t NCON = 250;   /* pins the rate*duration product => sharp ridge */
    int *store = (int *)malloc(NCON * n * sizeof(int));
    gd_cognate_data *concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
    int *sim = (int *)malloc(n * sizeof(int));
    for (size_t c = 0; c < NCON; c++) {
        simulate_relaxed(rng, t, rate, NULL, sim);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
        concepts[c].labels = row; concepts[c].nnodes = n;
    }
    free(sim);

    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt); opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, NCON, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 271u;
    cfg.n_warmup = 900; cfg.n_samples = 900; cfg.thin = 1;
    /* Hold the global rate fixed so the branch-rate/date trade-off concentrates
     * on the multipliers (with the rate free, it spreads across rate+ages too).
     * Widen the sigma_b prior so the multipliers can travel the ridge. */
    cfg.sample_rate = false; cfg.fixed_rate = rate;
    cfg.sample_branch = true;
    cfg.sigma_b_prior_sd = 0.60;

    /* Free internal ages: branch rate x duration trade off => confounded. */
    gd_calibration_set *cs = bounds_free_internal(t, R, &err);
    CHECK(cs != NULL);
    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);
    CHECK(res != NULL);
    printf("    confounding: min branch-duration corr = %.3f, flag = %d\n",
           res->min_branch_duration_corr, res->branch_rate_confounded);
    CHECK(!isnan(res->min_branch_duration_corr));
    CHECK(res->min_branch_duration_corr < -0.6);   /* strong negative trade-off */
    CHECK(res->branch_rate_confounded == 1);
    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);

    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    free(store); free(concepts);
    gd_tree_destroy(t); gd_rng_destroy(rng);
}

/* ------------------------------------------------- sigma_b shrinkage (strict) */

TEST(test_m9_sigma_b_shrinks_to_strict)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_quartet();
    CHECK(t != NULL);
    size_t n = t->nnodes;
    gd_rng *rng = gd_rng_create(9003u, 2u, &err);
    CHECK(rng != NULL);

    /* Strictly-clocked data (all multipliers = 1). */
    const double rate = 0.0005;
    const size_t NCON = 250;
    int *store = (int *)malloc(NCON * n * sizeof(int));
    gd_cognate_data *concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
    int *sim = (int *)malloc(n * sizeof(int));
    for (size_t c = 0; c < NCON; c++) {
        simulate_relaxed(rng, t, rate, NULL, sim);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
        concepts[c].labels = row; concepts[c].nnodes = n;
    }
    free(sim);

    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt); opt.epsilon = 0.0;
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, NCON, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);
    gd_calibration_set *cs = bounds_all_fixed(t, &err);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 555u;
    cfg.n_warmup = 800; cfg.n_samples = 800; cfg.thin = 1;
    cfg.rate_logmean = log(rate); cfg.rate_logsd = 1.0;
    cfg.sample_branch = true;
    cfg.sigma_b_prior_sd = 0.30;

    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);
    CHECK(res != NULL);
    /* Posterior-mean sigma_b (from the LOG_SIGMA_B column, exp'd). */
    size_t total = res->n_chains * res->n_draws, np = res->n_param;
    size_t p_sig = (size_t)-1;
    for (size_t p = 0; p < np; p++)
        if (res->param_type[p] == GD_MCMC_P_LOG_SIGMA_B) p_sig = p;
    CHECK(p_sig != (size_t)-1);
    double mean_sigma = 0.0;
    for (size_t i = 0; i < total; i++) mean_sigma += exp(res->samples[i * np + p_sig]);
    mean_sigma /= (double)total;
    printf("    strict-data sigma_b posterior mean = %.3f\n", mean_sigma);
    /* With strictly-clocked data the shrinkage prior pulls sigma_b toward 0. */
    CHECK(mean_sigma < 0.30);
    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);

    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    free(store); free(concepts);
    gd_tree_destroy(t); gd_rng_destroy(rng);
}

/* ------------------------------------------------------ robust contamination */

/* Fit a rate + contamination model; return the posterior mean epsilon and the
 * posterior mean global log-rate. Data = clean concepts plus `nloan` injected
 * miscoded concepts (random incompatible tip patterns). */
static void fit_contamination(uint64_t seed, size_t nclean, size_t nloan,
                              double *eps_mean_out, double *logr_mean_out)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_quartet();
    size_t n = t->nnodes;
    gd_rng *rng = gd_rng_create(seed, 4u, &err);
    const double rate = 0.0005;
    size_t NCON = nclean + nloan;
    int *store = (int *)malloc(NCON * n * sizeof(int));
    gd_cognate_data *concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
    int *sim = (int *)malloc(n * sizeof(int));
    size_t a = gd_tree_find_label(t, "a"), b = gd_tree_find_label(t, "b");
    size_t c_ = gd_tree_find_label(t, "c"), d = gd_tree_find_label(t, "d");
    for (size_t c = 0; c < nclean; c++) {
        simulate_relaxed(rng, t, rate, NULL, sim);
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = t->nodes[i].is_leaf ? sim[i] : -1;
        concepts[c].labels = row; concepts[c].nnodes = n;
    }
    /* Injected loans: the incompatible pattern a=c, b=d (no tree history makes
     * it; a borrowing / miscoding signature). */
    for (size_t c = nclean; c < NCON; c++) {
        int *row = store + c * n;
        for (size_t i = 0; i < n; i++) row[i] = -1;
        row[a] = 0; row[c_] = 0; row[b] = 1; row[d] = 1;
        concepts[c].labels = row; concepts[c].nnodes = n;
    }
    free(sim);

    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt);
    opt.incompatible = GD_INCOMPATIBLE_CONTAMINATION;   /* absorb loans */
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, NCON, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);
    gd_calibration_set *cs = bounds_all_fixed(t, &err);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = seed + 1u;
    cfg.n_warmup = 700; cfg.n_samples = 700; cfg.thin = 1;
    cfg.rate_logmean = log(rate); cfg.rate_logsd = 1.0;
    cfg.sample_contamination = true;

    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 2, &err);
    size_t total = res->n_chains * res->n_draws, np = res->n_param;
    size_t p_eps = (size_t)-1, p_rate = (size_t)-1;
    for (size_t p = 0; p < np; p++) {
        if (res->param_type[p] == GD_MCMC_P_LOGIT_EPS) p_eps = p;
        if (res->param_type[p] == GD_MCMC_P_LOGRATE) p_rate = p;
    }
    double me = 0.0, mr = 0.0;
    for (size_t i = 0; i < total; i++) {
        double z = res->samples[i * np + p_eps];
        me += 1.0 / (1.0 + exp(-z));
        mr += res->samples[i * np + p_rate];
    }
    *eps_mean_out = me / (double)total;
    *logr_mean_out = mr / (double)total;

    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);
    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    free(store); free(concepts);
    gd_tree_destroy(t); gd_rng_destroy(rng);
}

TEST(test_m9_robust_to_injected_loans)
{
    double eps_clean, logr_clean, eps_loans, logr_loans;
    fit_contamination(7000u, 150, 0, &eps_clean, &logr_clean);
    fit_contamination(7000u, 150, 40, &eps_loans, &logr_loans);
    printf("    epsilon: clean=%.3f  with-loans=%.3f | logr: clean=%.3f loans=%.3f\n",
           eps_clean, eps_loans, logr_clean, logr_loans);
    /* Injected loans push the contamination weight up. */
    CHECK(eps_loans > eps_clean);
    CHECK(eps_loans > 0.10);
    /* And the rate estimate stays robust (near the truth log(0.0005)), because
     * the contamination component absorbs the loans instead of biasing it. */
    CHECK(fabs(logr_loans - log(0.0005)) < 0.6);
}

static void run_tests(void)
{
    RUN(test_m9_branch_rate_recovery);
    RUN(test_m9_confounding_when_ages_free);
    RUN(test_m9_sigma_b_shrinks_to_strict);
    RUN(test_m9_robust_to_injected_loans);
}

TEST_MAIN()

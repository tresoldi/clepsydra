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

/*
 * Prior-only sampling of the relaxed clock must return its hierarchical prior,
 * sigma_b ~ HalfNormal(0.3) and b_i | sigma_b ~ N(-sigma_b^2/2, sigma_b), and
 * mix while doing so. That is Neal's funnel, where single-site moves on
 * sigma_b stall near zero: the non-centred joint move is what makes this pass.
 */
TEST(test_m9_prior_only_funnel)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_quartet();
    CHECK(t != NULL);
    size_t n = t->nnodes;
    int *store = (int *)malloc(4 * n * sizeof(int));
    gd_cognate_data concepts[4];
    for (size_t c = 0; c < 4; c++) {
        for (size_t i = 0; i < n; i++) store[c * n + i] = t->nodes[i].is_leaf ? (int)(i % 2) : -1;
        concepts[c].labels = store + c * n; concepts[c].nnodes = n;
    }
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt);
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, 4, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);
    gd_calibration_set *cs = bounds_all_fixed(t, &err);
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.prior_only = true; cfg.sample_branch = true; cfg.sigma_b_prior_sd = 0.30;
    cfg.rate_logmean = log(5e-4); cfg.rate_logsd = 0.5;
    cfg.n_warmup = 2000; cfg.n_samples = 10000; cfg.seed = 21u;
    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 4, &err);
    CHECK(res != NULL);
    size_t total = res->n_chains * res->n_draws, np = res->n_param, p_sig = (size_t)-1, p_b = (size_t)-1;
    size_t p_r = (size_t)-1;
    for (size_t p = 0; p < np; p++) {
        if (res->param_type[p] == GD_MCMC_P_LOGRATE) p_r = p;
        if (res->param_type[p] == GD_MCMC_P_LOG_SIGMA_B) p_sig = p;
        if (res->param_type[p] == GD_MCMC_P_BRANCH_LOGMULT && p_b == (size_t)-1) p_b = p;
    }
    CHECK(p_sig != (size_t)-1 && p_b != (size_t)-1);
    /* HalfNormal(0.3) quartiles: 0.3 * Phi^-1(0.625, 0.75, 0.875). */
    const double q[3] = { 0.3 * 0.318639, 0.3 * 0.674490, 0.3 * 1.150349 };
    double below[3] = {0, 0, 0}, btail = 0.0;
    for (size_t i = 0; i < total; i++) {
        double s = exp(res->samples[i * np + p_sig]);
        for (int k = 0; k < 3; k++) if (s < q[k]) below[k] += 1.0;
        if (res->samples[i * np + p_b] < -0.3) btail += 1.0;
    }
    /* Reference P(b < -0.3) by direct simulation of the prior. */
    gd_rng *rng = gd_rng_create(77u, 1u, &err);
    double ref = 0.0; const size_t NREF = 400000;
    for (size_t i = 0; i < NREF; i++) {
        double s = fabs(gd_rng_normal(rng, 0.0, 0.3));
        if (gd_rng_normal(rng, -0.5 * s * s, s) < -0.3) ref += 1.0;
    }
    ref /= (double)NREF;
    printf("    sigma_b quartile fractions %.3f %.3f %.3f (want 0.25 0.50 0.75); "
           "P(b<-0.3) %.4f (want %.4f); R-hat %.3f, bulk ESS %.0f\n",
           below[0] / (double)total, below[1] / (double)total, below[2] / (double)total,
           btail / (double)total, ref, res->rhat[p_sig], res->ess[p_sig]);
    for (int k = 0; k < 3; k++) CHECK_APPROX(below[k] / (double)total, 0.25 * (double)(k + 1), 0.03);
    /* The rate-shift move moves the log rate and every multiplier together.
     * It must leave the target alone: the log rate's marginal (which the
     * lexical-root prior tilts away from N(log 5e-4, 0.5)) is the same with
     * the move off. */
    CHECK(p_r != (size_t)-1);
    cfg.shift_move = false;
    gd_mcmc_result *off = gd_mcmc_run(&model, t, cs, &cfg, 4, &err);
    CHECK(off != NULL);
    double mom[2][2];
    for (int r = 0; r < 2; r++) {
        const gd_mcmc_result *x = r == 0 ? res : off;
        size_t tot = x->n_chains * x->n_draws;
        double mu = 0.0, v = 0.0;
        for (size_t i = 0; i < tot; i++) mu += x->samples[i * np + p_r];
        mu /= (double)tot;
        for (size_t i = 0; i < tot; i++) {
            double d = x->samples[i * np + p_r] - mu; v += d * d;
        }
        mom[r][0] = mu; mom[r][1] = sqrt(v / (double)(tot - 1));
    }
    printf("    log rate with shift move %.3f (sd %.3f, ESS %.0f), without %.3f (sd %.3f, ESS %.0f)\n",
           mom[0][0], mom[0][1], res->ess[p_r], mom[1][0], mom[1][1], off->ess[p_r]);
    CHECK_APPROX(mom[0][0], mom[1][0], 0.03);
    CHECK_APPROX(mom[0][1], mom[1][1], 0.02);
    gd_mcmc_result_destroy(off);
    CHECK_APPROX(btail / (double)total, ref, 0.01);
    CHECK(res->rhat[p_sig] < 1.01);
    CHECK(res->ess[p_sig] > 2000.0);
    gd_rng_destroy(rng);
    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);
    gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    free(store);
    gd_tree_destroy(t);
}

/* Two relaxed-clock moves reset branch multipliers along with ages. The clock
 * move scales heights above the oldest tip below each node and resets every
 * branch whose length does not scale with them (here the branch to tip a,
 * dated at 1000). The rate-compensated node move resets the branches meeting
 * the node it moves. Prior only, so the likelihood cannot hide a wrong
 * Jacobian or prior term, and the marginals are known exactly. Given the
 * root T, the uniform-conditional prior makes the other ages uniform on their
 * region and divides by its volume, so T is log-normal with log-mean
 * lex_logmean - rate_logmean and log-variance rate_logsd^2 + lex_logsd^2,
 * truncated to the root window [2000, 5000]; log rho given T is normal. By
 * quadrature: E[T] = 3126.9 (sd 832.0), E[log rho] = -7.7622 (sd 0.4751).
 * b(a) is N(-s^2/2, s) with s ~ HalfNormal(0.3), checked by simulation. */
TEST(test_m9_rate_compensated_moves)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = make_quartet();
    CHECK(t != NULL);
    size_t n = t->nnodes, ia = gd_tree_find_label(t, "a");
    t->nodes[ia].age = 1000.0;
    int *store = (int *)malloc(4 * n * sizeof(int));
    gd_cognate_data concepts[4];
    for (size_t c = 0; c < 4; c++) {
        for (size_t i = 0; i < n; i++) store[c * n + i] = t->nodes[i].is_leaf ? (int)(i % 2) : -1;
        concepts[c].labels = store + c * n; concepts[c].nnodes = n;
    }
    gd_t1_ws *ws = gd_t1_ws_create(t, &err);
    gd_t1_options opt; gd_t1_options_default(&opt);
    opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, concepts, 4, &opt, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(m, &model);
    gd_calibration_set *cs = bounds_free_internal(t, 3000.0, &err);
    CHECK(cs != NULL);
    size_t ir = t->root, i1 = gd_tree_find_label(t, "r1");
    cs->node_lo[ia] = cs->node_hi[ia] = 1000.0;
    cs->node_lo[ir] = 2000.0; cs->node_hi[ir] = 5000.0;
    cs->node_hi[i1] = 5000.0; cs->node_hi[gd_tree_find_label(t, "r2")] = 5000.0;
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.prior_only = true; cfg.sample_branch = true; cfg.sigma_b_prior_sd = 0.30;
    cfg.rate_logmean = log(5e-4); cfg.rate_logsd = 0.5;
    cfg.n_warmup = 2000; cfg.n_samples = 20000; cfg.seed = 33u;
    CHECK(cfg.root_lex_logmean == -1.2 && cfg.root_lex_logsd == 1.5);
    gd_mcmc_result *res = gd_mcmc_run(&model, t, cs, &cfg, 4, &err);
    CHECK(res != NULL);
    size_t np = res->n_param, idx[3] = { (size_t)-1, (size_t)-1, (size_t)-1 };
    for (size_t p = 0; p < np; p++) {
        if (res->param_type[p] == GD_MCMC_P_LOGRATE) idx[0] = p;
        if (res->param_type[p] == GD_MCMC_P_AGE && res->param_node[p] == ir) idx[1] = p;
        if (res->param_type[p] == GD_MCMC_P_BRANCH_LOGMULT && res->param_node[p] == ia) idx[2] = p;
    }
    CHECK(idx[0] != (size_t)-1 && idx[1] != (size_t)-1 && idx[2] != (size_t)-1);
    /* Exact moments; b(a) by direct simulation of its prior. */
    double want[3][2] = { { -7.7622, 0.4751 }, { 3126.9, 832.0 }, { 0.0, 0.0 } };
    gd_rng *rng = gd_rng_create(78u, 1u, &err);
    const size_t NREF = 400000;
    double s1 = 0.0, s2 = 0.0;
    for (size_t i = 0; i < NREF; i++) {
        double s = fabs(gd_rng_normal(rng, 0.0, 0.3));
        double b = gd_rng_normal(rng, -0.5 * s * s, s);
        s1 += b; s2 += b * b;
    }
    want[2][0] = s1 / (double)NREF;
    want[2][1] = sqrt(s2 / (double)NREF - want[2][0] * want[2][0]);
    const char *name[3] = { "log rate", "root age", "b(a)" };
    size_t tot = res->n_chains * res->n_draws;
    for (int k = 0; k < 3; k++) {
        double mu = 0.0, v = 0.0;
        for (size_t i = 0; i < tot; i++) mu += res->samples[i * np + idx[k]];
        mu /= (double)tot;
        for (size_t i = 0; i < tot; i++) {
            double d = res->samples[i * np + idx[k]] - mu; v += d * d;
        }
        double sd = sqrt(v / (double)(tot - 1));
        /* The clock move's runs are lumpier than ESS says (rare long
         * excursions along the ridge), hence 5 standard errors. */
        double se = want[k][1] / sqrt(res->ess[idx[k]]);
        printf("    %-8s %.4g (sd %.4g, ESS %.0f); exact %.4g (sd %.4g)\n",
               name[k], mu, sd, res->ess[idx[k]], want[k][0], want[k][1]);
        CHECK(fabs(mu - want[k][0]) < 5.0 * se);
        CHECK_APPROX(sd / want[k][1], 1.0, 0.04);
    }
    gd_rng_destroy(rng);
    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs); gd_mcmc_t1_destroy(m); gd_t1_ws_destroy(ws);
    gd_tree_destroy(t); free(store);
}

static void run_tests(void)
{
    RUN(test_m9_branch_rate_recovery);
    RUN(test_m9_confounding_when_ages_free);
    RUN(test_m9_sigma_b_shrinks_to_strict);
    RUN(test_m9_robust_to_injected_loans);
    RUN(test_m9_prior_only_funnel);
    RUN(test_m9_rate_compensated_moves);
}

TEST_MAIN()

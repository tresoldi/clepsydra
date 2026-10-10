/*
 * test_mcmc_tree.c — fixed-tree posterior sampler (M8).
 *
 * The M8 acceptance gate (PLAN.md): SBC for a small strict-clock tree;
 * reproducible seeded chains; no adaptation after warm-up; diagnostics detect
 * deliberately stuck chains (the last is covered in test_diagnostics.c).
 *
 * The model: the quartet ((a,b)r1,(c,d)r2)root with tips at age 0 and the root
 * fixed at R. The free parameters are the two internal ages r1, r2 (each with a
 * Uniform(0,R) prior via flat bounds) and the global strict-clock log-rate
 * (Normal prior). Data is simulated under T1 and refit; SBC checks the rank of
 * each true value among the posterior draws is uniform.
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_mcmc.h"
#include "gd_model_novel_state.h"
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_calibration.h"
#include "gd_diagnostics.h"
#include "gd_rng.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- fixtures */

static gd_tree *make_quartet(double R)
{
    const char *nwk = "((a,b)r1,(c,d)r2)root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    if (gd_newick_parse(nwk, strlen(nwk), &t, &err) != 0) return NULL;
    for (size_t i = 0; i < t->nnodes; i++) {
        t->nodes[i].age = t->nodes[i].is_leaf ? 0.0 : R / 2.0;
    }
    t->nodes[t->root].age = R;
    return t;
}

/* Flat bounds: tips and root fixed, the two internal non-root nodes free in
 * (0, R) — a Uniform(0,R) prior each (no calibration => flat over bounds). */
static gd_calibration_set *make_bounds(const gd_tree *t, double R, gd_error *err)
{
    gd_calibration_set *cs = gd_calset_create(0, err);
    if (cs == NULL) return NULL;
    size_t n = t->nnodes;
    cs->nnodes = n;
    cs->node_lo = (double *)malloc(n * sizeof(double));
    cs->node_hi = (double *)malloc(n * sizeof(double));
    if (cs->node_lo == NULL || cs->node_hi == NULL) {
        gd_calset_destroy(cs);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        if (t->nodes[i].is_leaf) { cs->node_lo[i] = 0.0; cs->node_hi[i] = 0.0; }
        else if (i == t->root)   { cs->node_lo[i] = R;   cs->node_hi[i] = R; }
        else                     { cs->node_lo[i] = 0.0; cs->node_hi[i] = R; }
    }
    cs->feasible = true;
    return cs;
}

/* One likelihood recompute helper: build the strict-clock T1 model over data. */
typedef struct {
    gd_tree *t;
    gd_t1_ws *ws;
    int *store;
    gd_cognate_data *concepts;
    size_t ncon;
    gd_mcmc_t1 *m;
    gd_t1_options opt;
} sbc_model;

static void sbc_model_free(sbc_model *sm)
{
    if (sm->m) gd_mcmc_t1_destroy(sm->m);
    if (sm->ws) gd_t1_ws_destroy(sm->ws);
    free(sm->store);
    free(sm->concepts);
    if (sm->t) gd_tree_destroy(sm->t);
}

/* --------------------------------------------------------------------- SBC */

TEST(test_mcmc_sbc_strict_clock_tree)
{
    const double R = 3000.0;
    const double logrmean = log(0.0004);
    const double logrsd = 0.4;
    const size_t NREP = 120;
    const size_t NCON = 50;
    const size_t L = 39;           /* retained draws per replicate: L + 1 = 40 ranks, 5 per bin */
    const size_t B = 8;            /* uniformity bins */
    size_t bins[8];
    for (size_t i = 0; i < B; i++) bins[i] = 0;
    size_t counted = 0;

    for (size_t rep = 0; rep < NREP; rep++) {
        gd_error err = GD_ERROR_INIT;
        gd_rng *rng = gd_rng_create(1000u + rep, 99u, &err);
        CHECK(rng != NULL);
        sbc_model sm;
        memset(&sm, 0, sizeof sm);
        sm.t = make_quartet(R);
        CHECK(sm.t != NULL);
        size_t n = sm.t->nnodes;
        size_t r1 = gd_tree_find_label(sm.t, "r1");
        size_t r2 = gd_tree_find_label(sm.t, "r2");

        /* Draw truth from the priors. */
        double a1_true = R * gd_rng_uniform_closed_open(rng);
        double a2_true = R * gd_rng_uniform_closed_open(rng);
        double logr_true = logrmean + gd_rng_normal(rng, 0.0, logrsd);
        sm.t->nodes[r1].age = a1_true;
        sm.t->nodes[r2].age = a2_true;
        sm.t->nodes[sm.t->root].age = R;

        /* Simulate NCON concepts; observe tips only. */
        sm.ncon = NCON;
        sm.store = (int *)malloc(NCON * n * sizeof(int));
        sm.concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
        CHECK(sm.store != NULL && sm.concepts != NULL);
        int *sim = (int *)malloc(n * sizeof(int));
        for (size_t c = 0; c < NCON; c++) {
            CHECK(gd_t1_simulate(rng, sm.t, exp(logr_true), sim, &err) == 0);
            int *row = sm.store + c * n;
            for (size_t i = 0; i < n; i++)
                row[i] = sm.t->nodes[i].is_leaf ? sim[i] : -1;
            sm.concepts[c].labels = row;
            sm.concepts[c].nnodes = n;
        }
        free(sim);

        /* Strict-clock T1 model, pure tree likelihood (no contamination). */
        sm.ws = gd_t1_ws_create(sm.t, &err);
        CHECK(sm.ws != NULL);
        gd_t1_options_default(&sm.opt);
        sm.opt.epsilon = 0.0;
        sm.opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
        sm.m = gd_mcmc_t1_create(sm.ws, sm.concepts, NCON, &sm.opt, &err);
        CHECK(sm.m != NULL);
        gd_mcmc_model model;
        gd_mcmc_t1_model(sm.m, &model);

        gd_calibration_set *cs = make_bounds(sm.t, R, &err);
        CHECK(cs != NULL);

        gd_mcmc_config cfg;
        gd_mcmc_config_default(&cfg);
        cfg.seed = 500u + rep;
        cfg.n_warmup = 350;
        cfg.n_samples = L;
        cfg.thin = 3;
        cfg.sample_rate = true;
        cfg.rate_logmean = logrmean;
        cfg.rate_logsd = logrsd;
        /* Truth is drawn flat below a fixed root, so the sampler is given that
         * same (bounded, hence proper) measure. */
        cfg.tree_prior = GD_TREE_PRIOR_NONE;
        cfg.clock_move = false;

        gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0, &err);
        CHECK(ch != NULL);
        CHECK(gd_mcmc_chain_run(ch, &err) == 0);
        CHECK(gd_mcmc_chain_nsamples(ch) == L);

        /* Identify the parameter columns. */
        size_t p_rate = (size_t)-1, p_a1 = (size_t)-1, p_a2 = (size_t)-1;
        for (size_t p = 0; p < gd_mcmc_chain_param_count(ch); p++) {
            size_t node = gd_mcmc_chain_param_node(ch, p);
            if (node == GD_NODE_NONE) p_rate = p;
            else if (node == r1) p_a1 = p;
            else if (node == r2) p_a2 = p;
        }
        CHECK(p_rate != (size_t)-1 && p_a1 != (size_t)-1 && p_a2 != (size_t)-1);

        /* Rank of each true value among the draws (in [0, L]). */
        struct { size_t p; double truth; } params[3] = {
            { p_rate, logr_true }, { p_a1, a1_true }, { p_a2, a2_true }
        };
        for (int k = 0; k < 3; k++) {
            size_t rank = 0;
            for (size_t d = 0; d < L; d++)
                if (gd_mcmc_chain_sample(ch, d, params[k].p) < params[k].truth) rank++;
            size_t b = rank * B / (L + 1);
            if (b >= B) b = B - 1;
            bins[b]++;
            counted++;
        }

        gd_mcmc_chain_destroy(ch);
        gd_calset_destroy(cs);
        sbc_model_free(&sm);
        gd_rng_destroy(rng);
    }

    /* Rank uniformity over B = 8 bins (dof = 7): the 0.99 quantile is ~18.5 and
     * the 0.999 is ~24. With 120 replicates × 3 parameters = 360 counts and
     * thin = 3 (to cut autocorrelation), a correctly-calibrated sampler sits
     * near the mean of 7; a biased sampler (wrong Jacobian, bad acceptance,
     * warmup adaptation leaking into sampling) spikes past the gate. The gate
     * is tightened from the old lenient 22 (dof 5) to 20 with the extra bins,
     * replicates, and thinning giving the power to catch subtler miscalibration. */
    CHECK(counted == 3 * NREP);
    double expected = (double)counted / (double)B;
    double chisq = 0.0;
    for (size_t i = 0; i < B; i++) {
        double d = (double)bins[i] - expected;
        chisq += d * d / expected;
    }
    printf("    SBC bins:");
    for (size_t i = 0; i < B; i++) printf(" %zu", bins[i]);
    printf("  chisq=%.2f (expected/bin=%.1f)\n", chisq, expected);
    CHECK(chisq < 20.0);
}

/* -------------------------------------------------- recovery + diagnostics */

/* Build a fixed simulated dataset once and return an owned model + bounds. */
static int build_dataset(sbc_model *sm, gd_calibration_set **cs_out,
                         double R, double logr_true, double a1, double a2,
                         size_t ncon, uint64_t seed, gd_error *err)
{
    memset(sm, 0, sizeof *sm);
    sm->t = make_quartet(R);
    if (sm->t == NULL) return -1;
    size_t n = sm->t->nnodes;
    size_t r1 = gd_tree_find_label(sm->t, "r1");
    size_t r2 = gd_tree_find_label(sm->t, "r2");
    sm->t->nodes[r1].age = a1;
    sm->t->nodes[r2].age = a2;
    sm->t->nodes[sm->t->root].age = R;
    gd_rng *rng = gd_rng_create(seed, 5u, err);
    if (rng == NULL) return -1;
    sm->ncon = ncon;
    sm->store = (int *)malloc(ncon * n * sizeof(int));
    sm->concepts = (gd_cognate_data *)malloc(ncon * sizeof(gd_cognate_data));
    int *sim = (int *)malloc(n * sizeof(int));
    for (size_t c = 0; c < ncon; c++) {
        gd_t1_simulate(rng, sm->t, exp(logr_true), sim, err);
        int *row = sm->store + c * n;
        for (size_t i = 0; i < n; i++)
            row[i] = sm->t->nodes[i].is_leaf ? sim[i] : -1;
        sm->concepts[c].labels = row;
        sm->concepts[c].nnodes = n;
    }
    free(sim);
    gd_rng_destroy(rng);
    sm->ws = gd_t1_ws_create(sm->t, err);
    gd_t1_options_default(&sm->opt);
    sm->opt.epsilon = 0.0;
    sm->opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    sm->m = gd_mcmc_t1_create(sm->ws, sm->concepts, ncon, &sm->opt, err);
    *cs_out = make_bounds(sm->t, R, err);
    return (sm->ws && sm->m && *cs_out) ? 0 : -1;
}

TEST(test_mcmc_recovery_and_rhat)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0, logr_true = log(0.0005);
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, logr_true, 1800.0, 1200.0, 150, 2024u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 4242u;
    cfg.n_warmup = 600;
    cfg.n_samples = 500;
    cfg.thin = 1;
    cfg.rate_logmean = log(0.0005);
    cfg.rate_logsd = 1.0;

    gd_mcmc_result *res = gd_mcmc_run(&model, sm.t, cs, &cfg, 4, &err);
    CHECK(res != NULL);
    CHECK(res->n_param == 3);
    CHECK(res->n_draws == 500 && res->n_chains == 4);

    /* Find the rate column and its posterior mean over all chains/draws. */
    size_t total = res->n_chains * res->n_draws;
    double mean_logr = 0.0;
    for (size_t i = 0; i < total; i++) mean_logr += res->samples[i * res->n_param + 0];
    mean_logr /= (double)total;
    /* Column 0 is the rate (sample_rate=1). Posterior mean near the truth. */
    CHECK(fabs(mean_logr - logr_true) < 0.5);      /* within ~0.5 in log space */

    /* Converged: R-hat near 1 and a healthy ESS on every parameter. */
    for (size_t p = 0; p < res->n_param; p++) {
        CHECK(res->rhat[p] < 1.1);
        CHECK(res->ess[p] > 50.0);
        CHECK(res->acc_rate[p] > 0.1 && res->acc_rate[p] < 0.9);
    }

    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

/* ------------------------------------------------------------ reproducible */

TEST(test_mcmc_reproducible)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1700.0, 1300.0, 40, 77u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 555u;
    cfg.n_warmup = 100;
    cfg.n_samples = 60;

    gd_mcmc_chain *a = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 3u, &err);
    CHECK(a != NULL && gd_mcmc_chain_run(a, &err) == 0);
    /* A second chain with the same seed and stream must reproduce it exactly. */
    gd_mcmc_chain *b = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 3u, &err);
    CHECK(b != NULL && gd_mcmc_chain_run(b, &err) == 0);
    size_t np = gd_mcmc_chain_param_count(a);
    CHECK(gd_mcmc_chain_nsamples(a) == gd_mcmc_chain_nsamples(b));
    for (size_t d = 0; d < gd_mcmc_chain_nsamples(a); d++)
        for (size_t p = 0; p < np; p++)
            CHECK(gd_mcmc_chain_sample(a, d, p) == gd_mcmc_chain_sample(b, d, p));
    gd_mcmc_chain_destroy(a);
    gd_mcmc_chain_destroy(b);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

/* ---------------------------------------------- no post-warm-up adaptation */

TEST(test_mcmc_no_post_warmup_adaptation)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1600.0, 1400.0, 40, 88u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 909u;
    cfg.n_warmup = 400;
    cfg.n_samples = 300;
    cfg.thin = 1;

    gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0, &err);
    CHECK(ch != NULL);
    size_t np = gd_mcmc_chain_param_count(ch);

    /* Step exactly through warm-up; record the (adapted) step sizes. */
    CHECK(gd_mcmc_chain_step(ch, cfg.n_warmup, &err) == 0);
    double *frozen = (double *)malloc(np * sizeof(double));
    for (size_t p = 0; p < np; p++) frozen[p] = gd_mcmc_chain_step_size(ch, p);

    /* Step through the sampling phase; the steps must not move. */
    CHECK(gd_mcmc_chain_step(ch, cfg.n_samples, &err) == 0);
    for (size_t p = 0; p < np; p++)
        CHECK(gd_mcmc_chain_step_size(ch, p) == frozen[p]);
    CHECK(gd_mcmc_chain_nsamples(ch) == cfg.n_samples);

    free(frozen);
    gd_mcmc_chain_destroy(ch);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

/* ------------------------------------------------------ checkpoint/resume */

TEST(test_mcmc_checkpoint_roundtrip)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1750.0, 1250.0, 40, 33u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 1234u;
    cfg.n_warmup = 200;
    cfg.n_samples = 80;

    /* Uninterrupted reference run. */
    gd_mcmc_chain *ref = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    CHECK(ref != NULL && gd_mcmc_chain_run(ref, &err) == 0);
    size_t np = gd_mcmc_chain_param_count(ref);
    size_t last = gd_mcmc_chain_nsamples(ref) - 1;
    double *ref_last = (double *)malloc(np * sizeof(double));
    for (size_t p = 0; p < np; p++) ref_last[p] = gd_mcmc_chain_sample(ref, last, p);
    double ref_ll = gd_mcmc_chain_sample_loglik(ref, last);
    gd_mcmc_chain_destroy(ref);

    /* Interrupted run: step halfway through warm-up, checkpoint, restore into a
     * fresh chain (same seed/stream/config), continue to the end. */
    gd_mcmc_chain *b = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    CHECK(b != NULL);
    CHECK(gd_mcmc_chain_step(b, cfg.n_warmup / 2, &err) == 0);
    unsigned char *buf = NULL; size_t len = 0;
    CHECK(gd_mcmc_chain_serialize(b, &buf, &len, &err) == 0);
    CHECK(buf != NULL && len > 0);
    gd_mcmc_chain_destroy(b);

    gd_mcmc_chain *c = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    CHECK(c != NULL);
    CHECK(gd_mcmc_chain_restore(c, buf, len, &err) == 0);
    /* Remaining iterations: total = warmup + samples*thin, already did warmup/2. */
    size_t total = cfg.n_warmup + cfg.n_samples * cfg.thin;
    CHECK(gd_mcmc_chain_step(c, total - cfg.n_warmup / 2, &err) == 0);
    CHECK(gd_mcmc_chain_nsamples(c) == cfg.n_samples);
    size_t clast = gd_mcmc_chain_nsamples(c) - 1;
    for (size_t p = 0; p < np; p++)
        CHECK(gd_mcmc_chain_sample(c, clast, p) == ref_last[p]);
    CHECK(gd_mcmc_chain_sample_loglik(c, clast) == ref_ll);

    /* E3.3: a checkpoint taken during sampling carries the draws already
     * stored, so the resumed chain's every draw matches the reference. */
    cfg.store_pointwise = true;
    gd_mcmc_chain *r2 = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    CHECK(r2 != NULL && gd_mcmc_chain_run(r2, &err) == 0);
    gd_mcmc_chain *e = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    CHECK(e != NULL);
    CHECK(gd_mcmc_chain_step(e, cfg.n_warmup + cfg.n_samples / 2, &err) == 0);
    unsigned char *buf2 = NULL; size_t len2 = 0;
    CHECK(gd_mcmc_chain_serialize(e, &buf2, &len2, &err) == 0);
    gd_mcmc_chain_destroy(e);
    gd_mcmc_chain *f = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    CHECK(gd_mcmc_chain_restore(f, buf2, len2, &err) == 0);
    CHECK(gd_mcmc_chain_nsamples(f) == cfg.n_samples / 2);
    CHECK(gd_mcmc_chain_step(f, total - cfg.n_warmup - cfg.n_samples / 2, &err) == 0);
    CHECK(gd_mcmc_chain_nsamples(f) == cfg.n_samples);
    size_t ndiff = 0;
    for (size_t k = 0; k < cfg.n_samples; k++) {
        for (size_t p = 0; p < np; p++)
            if (gd_mcmc_chain_sample(f, k, p) != gd_mcmc_chain_sample(r2, k, p)) ndiff++;
        if (gd_mcmc_chain_sample_loglik(f, k) != gd_mcmc_chain_sample_loglik(r2, k)) ndiff++;
        const double *pa = gd_mcmc_chain_sample_pointwise(f, k);
        const double *pb = gd_mcmc_chain_sample_pointwise(r2, k);
        CHECK(pa != NULL && pb != NULL);
        for (size_t i = 0; i < model.n_pointwise; i++) if (pa[i] != pb[i]) ndiff++;
    }
    CHECK(ndiff == 0);
    gd_mcmc_chain_destroy(f); gd_mcmc_chain_destroy(r2);

    /* The checkpoint is bound to its inputs: other data, another seed or
     * another chain stream are refused even though the shapes agree. */
    gd_mcmc_chain *g = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 2u, &err);
    CHECK(gd_mcmc_chain_restore(g, buf2, len2, &err) == GD_E_CHECKPOINT);
    gd_mcmc_chain_destroy(g);
    gd_mcmc_config cfg2 = cfg; cfg2.seed = cfg.seed + 1;
    gd_mcmc_chain *h = gd_mcmc_chain_create(&model, sm.t, cs, &cfg2, 1u, &err);
    CHECK(gd_mcmc_chain_restore(h, buf2, len2, &err) == GD_E_CHECKPOINT);
    gd_mcmc_chain_destroy(h);
    sbc_model sm2; gd_calibration_set *cs2 = NULL;
    CHECK(build_dataset(&sm2, &cs2, R, log(0.0005), 1750.0, 1250.0, 40, 34u, &err) == 0);
    gd_mcmc_model model2; gd_mcmc_t1_model(sm2.m, &model2);
    gd_mcmc_chain *k2 = gd_mcmc_chain_create(&model2, sm2.t, cs2, &cfg, 1u, &err);
    CHECK(k2 != NULL);
    CHECK(gd_mcmc_chain_restore(k2, buf2, len2, &err) == GD_E_CHECKPOINT);
    gd_mcmc_chain_destroy(k2);
    gd_calset_destroy(cs2);
    sbc_model_free(&sm2);
    free(buf2);

    /* A corrupted checksum is rejected. */
    cfg.store_pointwise = false;
    gd_mcmc_chain *d = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 1u, &err);
    buf[len - 1] ^= 0xFF;
    CHECK(gd_mcmc_chain_restore(d, buf, len, &err) == GD_E_CHECKPOINT);
    gd_mcmc_chain_destroy(d);

    free(buf);
    free(ref_last);
    gd_mcmc_chain_destroy(c);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

/* --------------------------------------------------- pointwise + fixed rate */

TEST(test_mcmc_pointwise_and_fixed_rate)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1700.0, 1300.0, 30, 44u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 2468u;
    cfg.n_warmup = 100;
    cfg.n_samples = 50;
    cfg.store_pointwise = true;
    /* Hold the rate fixed: sample only the two ages. */
    cfg.sample_rate = false;
    cfg.fixed_rate = 0.0005;

    gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0, &err);
    CHECK(ch != NULL);
    CHECK(gd_mcmc_chain_param_count(ch) == 2);   /* two ages, no rate */
    CHECK(gd_mcmc_chain_param_node(ch, 0) != GD_NODE_NONE);
    CHECK(gd_mcmc_chain_run(ch, &err) == 0);

    for (size_t d = 0; d < gd_mcmc_chain_nsamples(ch); d++) {
        const double *pw = gd_mcmc_chain_sample_pointwise(ch, d);
        CHECK(pw != NULL);
        double s = 0.0;
        for (size_t c = 0; c < sm.ncon; c++) s += pw[c];
        CHECK_APPROX(s, gd_mcmc_chain_sample_loglik(ch, d), 1e-9);
    }
    gd_mcmc_chain_destroy(ch);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

TEST(test_mcmc_block_scaling_same_target)
{
    /* The block age-scaling move must leave the posterior unchanged (correct
     * Jacobian / detailed balance): with it on or off, the posterior means of
     * the node ages agree within Monte-Carlo error. */
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1800.0, 1200.0, 150, 606u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.seed = 9u; cfg.n_warmup = 1500; cfg.n_samples = 1500; cfg.thin = 1;
    cfg.rate_logmean = log(0.0005); cfg.rate_logsd = 1.0;

    cfg.block_age_scaling = false;
    gd_mcmc_result *off = gd_mcmc_run(&model, sm.t, cs, &cfg, 3, &err);
    CHECK(off != NULL);
    cfg.block_age_scaling = true;
    gd_mcmc_result *on = gd_mcmc_run(&model, sm.t, cs, &cfg, 3, &err);
    CHECK(on != NULL);

    size_t np = off->n_param, tot_off = off->n_chains * off->n_draws;
    size_t tot_on = on->n_chains * on->n_draws;
    for (size_t p = 0; p < np; p++) {
        if (on->param_type[p] != GD_MCMC_P_AGE) continue;
        double mo = 0.0, mn = 0.0;
        for (size_t i = 0; i < tot_off; i++) mo += off->samples[i * np + p];
        for (size_t i = 0; i < tot_on; i++) mn += on->samples[i * np + p];
        mo /= (double)tot_off; mn /= (double)tot_on;
        /* Same target: means agree within ~4% of the age scale. */
        CHECK(fabs(mo - mn) < 0.04 * R);
        /* Both must have converged. */
        CHECK(on->rhat[p] < 1.1 && off->rhat[p] < 1.1);
    }
    gd_mcmc_result_destroy(off);
    gd_mcmc_result_destroy(on);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

/* ------------------------------------------------- progress and cancel (M10) */

/* What the callback saw: enough to check monotonicity, completion and which
 * chains reported, without depending on the reporting interval. */
typedef struct {
    size_t calls;
    size_t last_iter[8];      /* per chain */
    size_t total_seen;
    size_t chains_seen;
    int    monotone;
    size_t cancel_after;      /* 0 = never cancel */
} progress_log;

static int progress_cb(void *user, size_t chain, size_t iter, size_t total)
{
    progress_log *pl = (progress_log *)user;
    pl->calls++;
    pl->total_seen = total;
    if (chain < 8) {
        if (iter <= pl->last_iter[chain] && pl->last_iter[chain] != 0)
            pl->monotone = 0;
        if (pl->last_iter[chain] == 0) pl->chains_seen++;
        pl->last_iter[chain] = iter;
    }
    if (pl->cancel_after != 0 && pl->calls >= pl->cancel_after) return 1;
    return 0;
}

TEST(test_mcmc_progress_reports_and_completes)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1700.0, 1300.0, 30, 61u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    progress_log pl;
    memset(&pl, 0, sizeof pl);
    pl.monotone = 1;

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 4242u;
    cfg.n_warmup = 200;
    cfg.n_samples = 200;
    cfg.on_progress = progress_cb;
    cfg.progress_user = &pl;
    cfg.progress_every = 50;

    gd_mcmc_result *res = gd_mcmc_run(&model, sm.t, cs, &cfg, 2, &err);
    CHECK(res != NULL);

    /* 400 iterations per chain, reported every 50, over two chains. */
    CHECK(pl.calls >= 16);
    CHECK(pl.monotone == 1);
    CHECK(pl.chains_seen == 2);
    CHECK(pl.total_seen == 400);
    /* Each chain must report its final iteration, so a bar can reach 100%. */
    CHECK(pl.last_iter[0] == 400);
    CHECK(pl.last_iter[1] == 400);

    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

TEST(test_mcmc_progress_does_not_change_draws)
{
    /* The point of the callback's default being NULL: watching a run must not
     * perturb it. Same seed, same streams, with and without reporting. */
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1700.0, 1300.0, 30, 62u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 31337u;
    cfg.n_warmup = 150;
    cfg.n_samples = 150;

    gd_mcmc_result *quiet = gd_mcmc_run(&model, sm.t, cs, &cfg, 2, &err);
    CHECK(quiet != NULL);

    progress_log pl;
    memset(&pl, 0, sizeof pl);
    pl.monotone = 1;
    cfg.on_progress = progress_cb;
    cfg.progress_user = &pl;
    cfg.progress_every = 1;          /* the noisiest possible reporting */
    gd_mcmc_result *watched = gd_mcmc_run(&model, sm.t, cs, &cfg, 2, &err);
    CHECK(watched != NULL);

    CHECK(pl.calls == 600);          /* every iteration of both chains */
    CHECK(quiet->n_param == watched->n_param);
    CHECK(quiet->n_draws == watched->n_draws);
    size_t n = quiet->n_chains * quiet->n_draws * quiet->n_param;
    for (size_t i = 0; i < n; i++)
        CHECK(quiet->samples[i] == watched->samples[i]);   /* bit-identical */

    gd_mcmc_result_destroy(quiet);
    gd_mcmc_result_destroy(watched);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

TEST(test_mcmc_progress_can_cancel)
{
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1700.0, 1300.0, 30, 63u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    progress_log pl;
    memset(&pl, 0, sizeof pl);
    pl.monotone = 1;
    pl.cancel_after = 3;

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 909u;
    cfg.n_warmup = 500;
    cfg.n_samples = 500;
    cfg.on_progress = progress_cb;
    cfg.progress_user = &pl;
    cfg.progress_every = 10;

    /* A cancelled run yields no result at all: a half-finished chain is not a
     * posterior, and handing one back would invite it to be used as one. */
    gd_mcmc_result *res = gd_mcmc_run(&model, sm.t, cs, &cfg, 2, &err);
    CHECK(res == NULL);
    CHECK(err.code == GD_E_CANCELLED);
    CHECK(pl.calls == 3);
    CHECK(pl.last_iter[0] == 30);    /* stopped at the third report, not later */

    /* The same cancel through the chunked API. */
    err.code = 0;
    memset(&pl, 0, sizeof pl);
    pl.monotone = 1;
    pl.cancel_after = 1;
    gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0u, &err);
    CHECK(ch != NULL);
    CHECK(gd_mcmc_chain_step(ch, 100, &err) == GD_E_CANCELLED);
    CHECK(err.code == GD_E_CANCELLED);
    gd_mcmc_chain_destroy(ch);

    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

TEST(test_mcmc_result_accessors_agree_with_the_raw_layout)
{
    /* gd_mcmc_result.samples is flat, and its stride used to be spelled out at
     * seven call sites -- two of them inside gd_mcmc.c, where the sampler
     * decoded its own output. These accessors are that arithmetic written
     * once, so they have to agree with it exactly, and the finder has to agree
     * with the param_type/param_node scan it replaced. */
    gd_error err = GD_ERROR_INIT;
    const double R = 3000.0;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, R, log(0.0005), 1700.0, 1300.0, 30, 77u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

    gd_mcmc_config cfg;
    gd_mcmc_config_default(&cfg);
    cfg.seed = 909u;
    cfg.n_warmup = 60;
    cfg.n_samples = 60;
    gd_mcmc_result *res = gd_mcmc_run(&model, sm.t, cs, &cfg, 2, &err);
    CHECK(res != NULL);

    size_t np = res->n_param, total = res->n_chains * res->n_draws;
    int at_wrong = 0, col_wrong = 0, find_wrong = 0;

    for (size_t p = 0; p < np; p++) {
        for (size_t i = 0; i < total; i++) {
            if (gd_mcmc_result_at(res, p, i) != res->samples[i * np + p]) at_wrong++;
        }
    }

    double *col = (double *)malloc(total * sizeof(double));
    CHECK(col != NULL);
    if (col != NULL) {
        for (size_t p = 0; p < np; p++) {
            if (gd_mcmc_result_column(res, p, col, total) != GD_E_OK) { col_wrong++; continue; }
            for (size_t i = 0; i < total; i++)
                if (col[i] != res->samples[i * np + p]) col_wrong++;
        }
    }

    /* The finder must return the same column the old hand-built scan did. */
    for (size_t p = 0; p < np; p++) {
        if (res->param_type[p] != GD_MCMC_P_AGE) continue;
        if (gd_mcmc_result_find(res, GD_MCMC_P_AGE, res->param_node[p]) != p) find_wrong++;
    }
    size_t missing = gd_mcmc_result_find(res, GD_MCMC_P_AGE, (size_t)1u << 30);
    size_t lograte = gd_mcmc_result_find(res, GD_MCMC_P_LOGRATE, 0);
    /* Out of range is NAN and a rejected count, not a read past the end. */
    double past_p = gd_mcmc_result_at(res, np, 0);
    double past_i = gd_mcmc_result_at(res, 0, total);
    int bad_n = gd_mcmc_result_column(res, 0, col, total + 1);

    free(col);
    gd_mcmc_result_destroy(res);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);

    CHECK(at_wrong == 0);
    CHECK(col_wrong == 0);
    CHECK(find_wrong == 0);
    CHECK(missing == GD_MCMC_NO_PARAM);
    /* A global log-rate has no node, so the node argument is ignored. */
    CHECK(lograte != GD_MCMC_NO_PARAM);
    CHECK(isnan(past_p));
    CHECK(isnan(past_i));
    CHECK(bad_n == GD_E_INVALID_ARG);
}

static int dead_loglik(void *ctx, gd_tree *tree, const gd_mcmc_state *st,
                       double *pointwise, double *out, gd_error *err)
{
    (void)ctx; (void)tree; (void)st; (void)pointwise; (void)err;
    *out = -HUGE_VAL;
    return 0;
}

/* E0.12: a chain never marks the caller's calibration set feasible, and it
 * refuses to start from a state of zero density. */
TEST(test_mcmc_refuses_infeasible_and_zero_density_starts)
{
    gd_error err = GD_ERROR_INIT;
    sbc_model sm; gd_calibration_set *cs = NULL;
    CHECK(build_dataset(&sm, &cs, 3000.0, log(0.0005), 1700.0, 1300.0, 40, 78u, &err) == 0);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.n_warmup = 10; cfg.n_samples = 10;

    cs->feasible = false;
    gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0u, &err);
    bool still_infeasible = !cs->feasible;
    gd_mcmc_chain_destroy(ch);
    cs->feasible = true;
    CHECK(ch == NULL);
    CHECK(still_infeasible);
    CHECK(err.code == GD_E_CALIBRATION);

    /* A likelihood that is -inf at every state. */
    gd_mcmc_model dead = model;
    dead.loglik = dead_loglik;
    gd_error_clear(&err);
    ch = gd_mcmc_chain_create(&dead, sm.t, cs, &cfg, 0u, &err);
    int refused = (ch == NULL);
    gd_mcmc_chain_destroy(ch);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
    CHECK(refused);
}

/* ------------------------------------------------- E1: the proper prior */

/* Tips fixed at 0, every internal node free and unbounded above. */
static gd_calibration_set *make_free_bounds(const gd_tree *t, gd_error *err)
{
    gd_calibration_set *cs = gd_calset_create(0, err);
    if (cs == NULL) return NULL;
    size_t n = t->nnodes;
    cs->nnodes = n;
    cs->node_lo = (double *)malloc(n * sizeof(double));
    cs->node_hi = (double *)malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        cs->node_lo[i] = 0.0;
        cs->node_hi[i] = t->nodes[i].is_leaf ? 0.0 : HUGE_VAL;
    }
    cs->feasible = true;
    return cs;
}

static double phi_cdf(double z) { return 0.5 * erfc(-z / sqrt(2.0)); }

/*
 * Prior-only sampling reproduces the prior's analytic marginals:
 *   log a_root = log L - log rho ~ N(muL - mu_rho, sdL^2 + sd_rho^2)
 *   r1 / a_root ~ Uniform(0, 1)   (balanced quartet, uniform below the root)
 * This is the test that "posterior = prior" can be detected against.
 */
TEST(test_mcmc_prior_only_matches_analytic_prior)
{
    gd_error err = GD_ERROR_INIT;
    sbc_model sm; gd_calibration_set *cs0 = NULL;
    CHECK(build_dataset(&sm, &cs0, 3000.0, log(0.0005), 1700.0, 1300.0, 20, 11u, &err) == 0);
    gd_calset_destroy(cs0);
    gd_calibration_set *cs = make_free_bounds(sm.t, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.prior_only = true;
    cfg.rate_logmean = log(2e-4); cfg.rate_logsd = 0.5;
    cfg.root_lex_logmean = log(0.4); cfg.root_lex_logsd = 0.8;
    cfg.n_warmup = 2000; cfg.n_samples = 20000; cfg.seed = 3u;
    gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0u, &err);
    CHECK(ch != NULL);
    CHECK(gd_mcmc_chain_run(ch, &err) == 0);
    size_t root = sm.t->root, r1 = gd_tree_find_label(sm.t, "r1");
    size_t p_root = (size_t)-1, p_r1 = (size_t)-1;
    for (size_t p = 0; p < gd_mcmc_chain_param_count(ch); p++) {
        if (gd_mcmc_chain_param_node(ch, p) == root) p_root = p;
        if (gd_mcmc_chain_param_node(ch, p) == r1) p_r1 = p;
    }
    CHECK(p_root != (size_t)-1 && p_r1 != (size_t)-1);
    double mu = log(0.4) - log(2e-4), sd = sqrt(0.8 * 0.8 + 0.5 * 0.5);
    size_t N = gd_mcmc_chain_nsamples(ch);
    double m = 0, v = 0, below[3] = {0, 0, 0}, ratio_mean = 0, ratio_lo = 0;
    const double z[3] = { -1.0, 0.0, 1.0 };
    for (size_t d = 0; d < N; d++) {
        double lr = log(gd_mcmc_chain_sample(ch, d, p_root));
        m += lr; v += lr * lr;
        for (int k = 0; k < 3; k++) if (lr < mu + z[k] * sd) below[k] += 1.0;
        double q = gd_mcmc_chain_sample(ch, d, p_r1) / gd_mcmc_chain_sample(ch, d, p_root);
        ratio_mean += q;
        if (q < 0.25) ratio_lo += 1.0;
    }
    m /= (double)N; v = v / (double)N - m * m;
    gd_mcmc_chain_destroy(ch);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
    printf("    prior-only log root: mean %.3f (want %.3f), sd %.3f (want %.3f)\n",
           m, mu, sqrt(v), sd);
    CHECK_APPROX(m, mu, 0.08);
    CHECK_APPROX(sqrt(v), sd, 0.06);
    for (int k = 0; k < 3; k++) CHECK_APPROX(below[k] / (double)N, phi_cdf(z[k]), 0.04);
    CHECK_APPROX(ratio_mean / (double)N, 0.5, 0.03);
    CHECK_APPROX(ratio_lo / (double)N, 0.25, 0.03);
}

/*
 * SBC under the proper prior (Gate E1), over the scenarios the gate names:
 *   FREE      — tips modern, no calibration: the lexical-time case;
 *   NORMAL    — plus a normal_node calibration on the root;
 *   LOGNORMAL — plus an offset_lognormal_node calibration on the root;
 *   ANCIENT   — tip a dated 300 BP;
 *   YULE      — the Yule alternative prior, kappa 2.
 * Truth is drawn from each scenario's prior by exact rejection, without the
 * sampler: rho and L from their priors, a_root = L / rho, accepted with
 * probability cal(a_root) / max cal (and a_root above the oldest tip); the
 * internal ages from their conditional prior. Ranks are binned per parameter,
 * and the joint log-likelihood is an extra test quantity (Modrak et al.
 * 2025): it is what fails when the posterior is just the prior.
 */
typedef enum { SC_FREE, SC_NORMAL, SC_LOGNORMAL, SC_ANCIENT, SC_YULE, SC_MIXRATES } sbc_scenario;

static void run_sbc_scenario(sbc_scenario sc, size_t nrep, double *chisq_out)
{
    const double lrm = log(3e-4), lrs = 0.3, lxm = log(0.5), lxs = 0.5;
    const size_t NCON = 60, Ld = 39, B = 8;
    size_t bins[5][8];
    memset(bins, 0, sizeof bins);
    const double TIP = 300.0, CAL_MEAN = 1800.0, CAL_SD = 250.0;
    const double OFF = 900.0, OLM = log(800.0), OLS = 0.35;
    for (size_t rep = 0; rep < nrep; rep++) {
        gd_error err = GD_ERROR_INIT;
        gd_rng *rng = gd_rng_create(30000u + 1000u * (uint64_t)sc + rep, 17u, &err);
        sbc_model sm;
        memset(&sm, 0, sizeof sm);
        sm.t = make_quartet(1000.0);
        size_t n = sm.t->nnodes, root = sm.t->root;
        size_t r1 = gd_tree_find_label(sm.t, "r1"), r2 = gd_tree_find_label(sm.t, "r2");
        size_t ta = gd_tree_find_label(sm.t, "a");
        double floor_r1 = (sc == SC_ANCIENT) ? TIP : 0.0;

        /* Truth by rejection from the scenario's prior. */
        double logr, R;
        for (;;) {
            logr = lrm + gd_rng_normal(rng, 0.0, lrs);
            R = exp(lxm + gd_rng_normal(rng, 0.0, lxs)) / exp(logr);
            double acc = 1.0;
            if (sc == SC_NORMAL) {
                double z = (R - CAL_MEAN) / CAL_SD;
                acc = exp(-0.5 * z * z);
            } else if (sc == SC_LOGNORMAL) {
                /* density / its maximum, the max taken at the mode */
                if (R <= OFF) acc = 0.0;
                else {
                    double x = R - OFF, mode = exp(OLM - OLS * OLS);
                    double lz = (log(x) - OLM) / OLS, lzm = (log(mode) - OLM) / OLS;
                    acc = exp(-0.5 * lz * lz - log(x) + 0.5 * lzm * lzm + log(mode));
                }
            } else if (sc == SC_ANCIENT) {
                acc = (R > TIP) ? 1.0 : 0.0;
            }
            if (gd_rng_uniform_open(rng) < acc) break;
        }
        double a1, a2;
        if (sc == SC_YULE) {
            /* Truncated Exponential(kappa / R) on (0, R), by inversion. */
            double k = 2.0;
            double u1 = gd_rng_uniform_open(rng), u2 = gd_rng_uniform_open(rng);
            a1 = -R / k * log(1.0 - u1 * (1.0 - exp(-k)));
            a2 = -R / k * log(1.0 - u2 * (1.0 - exp(-k)));
        } else {
            a1 = floor_r1 + (R - floor_r1) * gd_rng_uniform_open(rng);
            a2 = R * gd_rng_uniform_open(rng);
        }
        sm.t->nodes[root].age = R;
        sm.t->nodes[r1].age = a1;
        sm.t->nodes[r2].age = a2;
        sm.t->nodes[ta].age = (sc == SC_ANCIENT) ? TIP : 0.0;

        sm.ncon = NCON;
        sm.store = (int *)malloc(NCON * n * sizeof(int));
        sm.concepts = (gd_cognate_data *)malloc(NCON * sizeof(gd_cognate_data));
        int *sim = (int *)malloc(n * sizeof(int));
        /* SC_MIXRATES: sigma_c ~ HalfNormal(1), each concept's rate times a
         * log-normal multiplier with mean 1, as the mixture model assumes. */
        double sc_true = (sc == SC_MIXRATES) ? fabs(gd_rng_normal(rng, 0.0, 1.0)) : 0.0;
        for (size_t c = 0; c < NCON; c++) {
            double lam = (sc == SC_MIXRATES)
                       ? exp(sc_true * gd_rng_normal(rng, 0.0, 1.0) - 0.5 * sc_true * sc_true) : 1.0;
            CHECK(gd_t1_simulate(rng, sm.t, exp(logr) * lam, sim, &err) == 0);
            int *row = sm.store + c * n;
            for (size_t i = 0; i < n; i++) row[i] = sm.t->nodes[i].is_leaf ? sim[i] : -1;
            sm.concepts[c].labels = row;
            sm.concepts[c].nnodes = n;
        }
        free(sim);
        sm.ws = gd_t1_ws_create(sm.t, &err);
        gd_t1_options_default(&sm.opt);
        sm.opt.epsilon = 0.0;
        sm.opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
        sm.m = gd_mcmc_t1_create(sm.ws, sm.concepts, NCON, &sm.opt, &err);
        gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);

        /* Log-likelihood at the truth, for the extra test quantity. */
        double ll_true;
        {
            gd_mcmc_state st;
            memset(&st, 0, sizeof st);
            st.rates = &logr; st.n_rates = 1; st.epsilon = 0.0;
            st.rate_mix_sd = sc_true;
            CHECK(model.loglik(model.ctx, sm.t, &st, NULL, &ll_true, &err) == 0);
        }

        /* The calibration set the sampler sees. */
        gd_calibration_set *cs = NULL;
        size_t ncal = (sc == SC_NORMAL || sc == SC_LOGNORMAL) ? 1 : 0;
        cs = gd_calset_create(ncal + (sc == SC_ANCIENT ? 1 : 0), &err);
        size_t k = 0;
        if (ncal) {
            gd_calibration *cl = &cs->cals[k++];
            cl->target_id = (char *)malloc(5); memcpy(cl->target_id, "root", 5);
            if (sc == SC_NORMAL) { cl->type = GD_CAL_NORMAL_NODE; cl->p1 = CAL_MEAN; cl->p2 = CAL_SD; }
            else { cl->type = GD_CAL_OFFSET_LOGNORMAL_NODE; cl->p1 = OFF; cl->p2 = OLM; cl->p3 = OLS; }
        }
        if (sc == SC_ANCIENT) {
            gd_calibration *cl = &cs->cals[k++];
            cl->target_id = (char *)malloc(2); memcpy(cl->target_id, "a", 2);
            cl->type = GD_CAL_FIXED_TIP; cl->p1 = TIP; cl->is_hard = true;
        }
        CHECK(gd_calibrations_resolve(cs, sm.t, &err) == 0);
        CHECK(gd_calibrations_propagate(cs, sm.t, &err) == 0);

        gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
        cfg.seed = 900u + 1000u * (uint64_t)sc + rep;
        cfg.n_warmup = 600; cfg.n_samples = Ld; cfg.thin = 6;
        cfg.rate_logmean = lrm; cfg.rate_logsd = lrs;
        cfg.root_lex_logmean = lxm; cfg.root_lex_logsd = lxs;
        if (sc == SC_YULE) { cfg.tree_prior = GD_TREE_PRIOR_YULE; cfg.yule_kappa = 2.0; }
        if (sc == SC_MIXRATES) { cfg.sample_rate_mix = true; cfg.rate_mix_prior_sd = 1.0; }
        gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0, &err);
        CHECK(ch != NULL);
        CHECK(gd_mcmc_chain_run(ch, &err) == 0);
        size_t pc[4] = { (size_t)-1, (size_t)-1, (size_t)-1, (size_t)-1 };
        for (size_t p = 0; p < gd_mcmc_chain_param_count(ch); p++) {
            size_t node = gd_mcmc_chain_param_node(ch, p);
            if (gd_mcmc_chain_param_type(ch, p) == GD_MCMC_P_LOGRATE) pc[0] = p;
            else if (node == GD_NODE_NONE) continue;
            else if (node == root) pc[1] = p;
            else if (node == r1) pc[2] = p;
            else if (node == r2) pc[3] = p;
        }
        double truth[4] = { logr, R, a1, a2 };
        for (int q = 0; q < 5; q++) {
            size_t rank = 0;
            for (size_t d = 0; d < Ld; d++) {
                double v = (q < 4) ? gd_mcmc_chain_sample(ch, d, pc[q])
                                   : gd_mcmc_chain_sample_loglik(ch, d);
                double t = (q < 4) ? truth[q] : ll_true;
                if (v < t) rank++;
            }
            size_t b = rank * B / (Ld + 1);
            bins[q][b < B ? b : B - 1]++;
        }
        gd_mcmc_chain_destroy(ch);
        gd_calset_destroy(cs);
        sbc_model_free(&sm);
        gd_rng_destroy(rng);
    }
    static const char *const names[5] = { "log_rate", "root", "r1", "r2", "loglik" };
    static const char *const scn[6] = { "free", "normal", "lognormal", "ancient", "yule", "mixrates" };
    for (int q = 0; q < 5; q++) {
        double e = (double)nrep / (double)B, c2 = 0.0;
        printf("    SBC %-9s %-8s:", scn[sc], names[q]);
        for (size_t i = 0; i < B; i++) {
            double d = (double)bins[q][i] - e;
            c2 += d * d / e;
            printf(" %zu", bins[q][i]);
        }
        printf("  chisq=%.2f\n", c2);
        chisq_out[q] = c2;
    }
}

TEST(test_mcmc_sbc_proper_prior_scenarios)
{
    const sbc_scenario all[6] = { SC_FREE, SC_NORMAL, SC_LOGNORMAL, SC_ANCIENT, SC_YULE, SC_MIXRATES };
    for (int s = 0; s < 6; s++) {
        double c2[5];
        run_sbc_scenario(all[s], 160, c2);
        /* dof 7: 0.999 quantile 24.3; 30 tests, ~3% family-wise. */
        for (int q = 0; q < 5; q++) CHECK(c2[q] < 24.3);
    }
}

/* E1.3: the Yule alternative. Under the prior alone each non-root internal
 * age, as a fraction s of the root's, is Exponential(kappa) truncated to
 * (0, 1): E[s] = 1/kappa - 1/(e^kappa - 1). The root is the same lexical-root
 * prior as before. A dated tip is refused. */
TEST(test_mcmc_yule_prior_ratios)
{
    gd_error err = GD_ERROR_INIT;
    sbc_model sm; gd_calibration_set *cs0 = NULL;
    CHECK(build_dataset(&sm, &cs0, 3000.0, log(0.0005), 1700.0, 1300.0, 20, 12u, &err) == 0);
    gd_calset_destroy(cs0);
    gd_calibration_set *cs = make_free_bounds(sm.t, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.prior_only = true;
    cfg.tree_prior = GD_TREE_PRIOR_YULE;
    cfg.yule_kappa = 3.0;
    cfg.n_warmup = 2000; cfg.n_samples = 20000; cfg.seed = 5u;
    gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0u, &err);
    CHECK(ch != NULL);
    CHECK(gd_mcmc_chain_run(ch, &err) == 0);
    size_t root = sm.t->root, r1 = gd_tree_find_label(sm.t, "r1");
    size_t p_root = (size_t)-1, p_r1 = (size_t)-1;
    for (size_t p = 0; p < gd_mcmc_chain_param_count(ch); p++) {
        if (gd_mcmc_chain_param_node(ch, p) == root) p_root = p;
        if (gd_mcmc_chain_param_node(ch, p) == r1) p_r1 = p;
    }
    size_t N = gd_mcmc_chain_nsamples(ch);
    double mean = 0.0;
    for (size_t d = 0; d < N; d++)
        mean += gd_mcmc_chain_sample(ch, d, p_r1) / gd_mcmc_chain_sample(ch, d, p_root);
    mean /= (double)N;
    gd_mcmc_chain_destroy(ch);
    double want = 1.0 / 3.0 - 1.0 / (exp(3.0) - 1.0);
    printf("    Yule prior-only E[s] = %.4f (want %.4f)\n", mean, want);
    CHECK_APPROX(mean, want, 0.02);

    /* A dated tip: refused, with a pointer to the uniform prior. */
    cs->node_lo[gd_tree_find_label(sm.t, "a")] = 100.0;
    cs->node_hi[gd_tree_find_label(sm.t, "a")] = 100.0;
    gd_error_clear(&err);
    ch = gd_mcmc_chain_create(&model, sm.t, cs, &cfg, 0u, &err);
    int refused = (ch == NULL);
    gd_mcmc_chain_destroy(ch);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
    CHECK(refused);
    CHECK(strstr(err.message, "uniform") != NULL);
}

/* The subtree-scaling move keeps the prior: on the caterpillar
 * (((a,b)x,c)y,d)root under the uniform prior, x/root ~ Beta(1,2) (mean 1/3)
 * and y/root ~ Beta(2,1) (mean 2/3). The mover y scales y and x together, so
 * a wrong Jacobian shows up here. Run with and without the move. */
TEST(test_mcmc_subtree_move_keeps_prior)
{
    for (int with = 0; with < 2; with++) {
        gd_error err = GD_ERROR_INIT;
        const char *nwk = "(((a,b)x,c)y,d)root;";
        gd_tree *t = NULL;
        CHECK(gd_newick_parse(nwk, strlen(nwk), &t, &err) == 0);
        size_t n = t->nnodes;
        /* A dummy concept: prior-only never evaluates it. */
        int *lab = (int *)malloc(n * sizeof(int));
        for (size_t i = 0; i < n; i++) lab[i] = t->nodes[i].is_leaf ? 0 : -1;
        gd_cognate_data d = { lab, n };
        gd_t1_ws *ws = gd_t1_ws_create(t, &err);
        gd_t1_options opt; gd_t1_options_default(&opt);
        gd_mcmc_t1 *m = gd_mcmc_t1_create(ws, &d, 1, &opt, &err);
        gd_mcmc_model model; gd_mcmc_t1_model(m, &model);
        gd_calibration_set *cs = make_free_bounds(t, &err);
        gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
        cfg.prior_only = true;
        cfg.subtree_move = with ? true : false;
        cfg.n_warmup = 2000; cfg.n_samples = 30000; cfg.seed = 21u;
        gd_mcmc_chain *ch = gd_mcmc_chain_create(&model, t, cs, &cfg, 0u, &err);
        CHECK(ch != NULL);
        CHECK(gd_mcmc_chain_run(ch, &err) == 0);
        size_t px = (size_t)-1, py = (size_t)-1, pr = (size_t)-1;
        size_t x = gd_tree_find_label(t, "x"), y = gd_tree_find_label(t, "y");
        for (size_t p = 0; p < gd_mcmc_chain_param_count(ch); p++) {
            size_t node = gd_mcmc_chain_param_node(ch, p);
            if (node == x) px = p; else if (node == y) py = p; else if (node == t->root) pr = p;
        }
        double mx = 0, my = 0;
        size_t N = gd_mcmc_chain_nsamples(ch);
        for (size_t k = 0; k < N; k++) {
            double R = gd_mcmc_chain_sample(ch, k, pr);
            mx += gd_mcmc_chain_sample(ch, k, px) / R;
            my += gd_mcmc_chain_sample(ch, k, py) / R;
        }
        mx /= (double)N; my /= (double)N;
        gd_mcmc_chain_destroy(ch);
        gd_calset_destroy(cs);
        gd_mcmc_t1_destroy(m);
        gd_t1_ws_destroy(ws);
        free(lab);
        gd_tree_destroy(t);
        printf("    subtree move %s: E[x/R] = %.4f (1/3), E[y/R] = %.4f (2/3)\n",
               with ? "on " : "off", mx, my);
        CHECK_APPROX(mx, 1.0 / 3.0, 0.015);
        CHECK_APPROX(my, 2.0 / 3.0, 0.015);
    }
}

/*
 * Parallel chains (GD_THREADS) give the same draws as sequential ones, bit for
 * bit: each chain's random stream depends only on (seed, chain), and each runs
 * on its own tree copy and model clone. Relaxed clock and concept rates on, so
 * every move type is exercised. Without threads both runs are sequential.
 */
TEST(test_mcmc_parallel_chains_match_sequential)
{
    gd_error err = GD_ERROR_INIT;
    sbc_model sm; gd_calibration_set *cs0 = NULL;
    CHECK(build_dataset(&sm, &cs0, 3000.0, log(0.0005), 1700.0, 1300.0, 60, 17u, &err) == 0);
    gd_calset_destroy(cs0);
    gd_calibration_set *cs = make_free_bounds(sm.t, &err);
    gd_mcmc_model model; gd_mcmc_t1_model(sm.m, &model);
    gd_mcmc_config cfg; gd_mcmc_config_default(&cfg);
    cfg.rate_logmean = log(5e-4); cfg.rate_logsd = 0.5;
    cfg.sample_branch = true; cfg.sample_rate_mix = true;
    cfg.n_warmup = 300; cfg.n_samples = 300; cfg.seed = 99u;
    cfg.store_pointwise = true;
    cfg.threads = 1;
    gd_mcmc_result *a = gd_mcmc_run(&model, sm.t, cs, &cfg, 4, &err);
    CHECK(a != NULL);
    double end_a = sm.t->nodes[sm.t->root].age;
    cfg.threads = 0;
    gd_mcmc_result *b = gd_mcmc_run(&model, sm.t, cs, &cfg, 4, &err);
    CHECK(b != NULL);
    double end_b = sm.t->nodes[sm.t->root].age;
    CHECK(a->n_param == b->n_param && a->n_draws == b->n_draws);
    size_t n = a->n_chains * a->n_draws * a->n_param, diff = 0;
    for (size_t i = 0; i < n; i++) if (a->samples[i] != b->samples[i]) diff++;
    for (size_t i = 0; i < a->n_chains * a->n_draws; i++) if (a->loglik[i] != b->loglik[i]) diff++;
    for (size_t i = 0; i < a->n_chains * a->n_draws * a->n_concepts; i++)
        if (a->pointwise[i] != b->pointwise[i]) diff++;
    printf("    parallel vs sequential: %zu differing values of %zu; final root age %s\n",
           diff, n, end_a == end_b ? "equal" : "differs");
    CHECK(diff == 0);
    CHECK(end_a == end_b);
    gd_mcmc_result_destroy(a); gd_mcmc_result_destroy(b);
    gd_calset_destroy(cs);
    sbc_model_free(&sm);
}

static void run_tests(void)
{
    RUN(test_mcmc_reproducible);
    RUN(test_mcmc_refuses_infeasible_and_zero_density_starts);
    RUN(test_mcmc_block_scaling_same_target);
    RUN(test_mcmc_no_post_warmup_adaptation);
    RUN(test_mcmc_checkpoint_roundtrip);
    RUN(test_mcmc_pointwise_and_fixed_rate);
    RUN(test_mcmc_recovery_and_rhat);
    RUN(test_mcmc_progress_reports_and_completes);
    RUN(test_mcmc_result_accessors_agree_with_the_raw_layout);
    RUN(test_mcmc_progress_does_not_change_draws);
    RUN(test_mcmc_progress_can_cancel);
    RUN(test_mcmc_sbc_strict_clock_tree);
    RUN(test_mcmc_prior_only_matches_analytic_prior);
    RUN(test_mcmc_parallel_chains_match_sequential);
    RUN(test_mcmc_yule_prior_ratios);
    RUN(test_mcmc_subtree_move_keeps_prior);
    RUN(test_mcmc_sbc_proper_prior_scenarios);
}

TEST_MAIN()

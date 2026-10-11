/*
 * gd_mcmc.c — fixed-tree Metropolis-within-Gibbs sampler for Clepsydra
 * (M8 core + M9 relaxed clock / robust observations).
 *
 * A typed-parameter engine: each active block contributes one or more
 * parameters (log-rate, node ages, log sigma_b, logit epsilon, branch
 * log-multipliers, concept deviations). Each sweep does a single-site
 * random-walk Metropolis update of every parameter on its natural sampling
 * scale, recomputing the likelihood through the model callback (except the
 * hyperparameter sigma_b, which enters only the prior). Adaptation is
 * Robbins-Monro during warm-up only; per-chain PCG32 streams make runs
 * bit-reproducible; the live state serializes (versioned + checksummed).
 *
 * See gd_mcmc.h and docs/inference.md §2.5, docs/model-spec.md §4/§7.
 */
#include "gd_mcmc.h"
#include "gd_memory.h"
#include "gd_math.h"
#include "gd_tree.h"
#include "gd_tree_build.h"
#include "gd_calibration.h"
#include "gd_diagnostics.h"
#include "gd_rng.h"
#include "gd_tree_prior.h"
#include "gd.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef GD_THREADS
#include <pthread.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GD_MCMC_TARGET_ACCEPT 0.44
#define GD_MCMC_ADAPT_T0      10.0
#define GD_MCMC_FIX_TOL       1e-9
#define GD_MCMC_CKPT_MAGIC    0x394d4447u    /* "GDM9" */
#define GD_MCMC_CKPT_VERSION  8u
#define GD_MCMC_CONFOUND_THRESH 0.6          /* |corr| beyond this => confounded */

/* ------------------------------------------------------------- parameters */

typedef struct {
    gd_mcmc_ptype type;
    size_t node;        /* AGE / BRANCH_LOGMULT: node index; else GD_NODE_NONE */
    size_t concept;     /* CONCEPT_DEV: concept index; else 0 */
    double value;       /* current value on the sampling scale (AGE mirrors age) */
    double step, adapt_count, acc, prop;
} gd_param;

struct gd_mcmc_chain {
    const gd_mcmc_model *model;
    gd_tree *tree;
    gd_calibration_set *cal;
    gd_mcmc_config cfg;
    gd_rng *rng;
    uint64_t stream;          /* chain index, for progress reports */
    uint64_t binding;         /* hash of everything a checkpoint must match */

    size_t nparam;
    gd_param *params;

    /* Assembled state, kept in sync with the parameter values. */
    double logr;
    double sigma_b;
    double sigma_c;           /* among-concept spread; 0 when not sampled */
    double *tip_hazard;       /* nnodes, or NULL: -log(1 - e_v) per tip */
    double epsilon;
    double *branch_logmult;   /* nnodes, or NULL */
    double *concept_dev;      /* nconcepts, or NULL */
    int have_branch, have_concept;
    size_t nnodes, nconcepts;

    double cur_loglik, cur_logprior;
    uint64_t iter;

    /* Block age-scaling move state. */
    int     have_block;
    size_t  n_age;            /* free AGE params of internal nodes */
    size_t *age_params;       /* n_age indices into params[] */
    bool   *scaled;           /* nnodes: does the scale move move this node? */
    double *age_scratch;      /* n_age; proposed ages */
    double  block_step, block_adapt_count, block_acc, block_prop;

    /* Clock-preserving move state. */
    int     have_clock;
    double  clock_step, clock_adapt_count, clock_acc, clock_prop;

    /* Subtree-scaling move: sub_start[k]..sub_start[k+1] index sub_nodes,
     * the free internal nodes under the k-th mover sub_root[k]. */
    size_t  n_sub;
    size_t *sub_root, *sub_start, *sub_nodes;
    double *sub_scratch;
    double  sub_step, sub_adapt_count, sub_acc, sub_prop;

    /* Non-centred sigma_b move (relaxed clock): sigma_b and every branch
     * multiplier move together with the standardised deviations fixed. */
    size_t  funnel_param;     /* index of LOG_SIGMA_B; nparam if absent */
    double *funnel_scratch;   /* nparam: old branch values */
    double  funnel_step, funnel_adapt_count, funnel_acc, funnel_prop;

    /* Rate-shift move (relaxed clock): log rate up by delta, every branch
     * multiplier down by delta, so each branch's effective rate is fixed. */
    size_t  shift_param;      /* index of LOGRATE; nparam if the move is off */
    double  shift_step, shift_adapt_count, shift_acc, shift_prop;

    /* Proper prior on the ages (NULL under GD_TREE_PRIOR_NONE), and its
     * current value, kept as part of cur_logprior. */
    gd_tree_prior *tp;
    double  cur_tp;

    size_t cap, n_stored;
    double *samples, *sample_ll, *pointwise, *pw_scratch;
};

/* ------------------------------------------------------------- densities */

static double norm_logpdf(double x, double mean, double sd)
{
    double z = (x - mean) / sd;
    return -0.5 * z * z - log(sd) - 0.5 * log(2.0 * M_PI);
}

static double halfnormal_logpdf(double x, double sd)
{
    if (x < 0.0) return -HUGE_VAL;
    double z = x / sd;
    return log(2.0) - 0.5 * z * z - log(sd) - 0.5 * log(2.0 * M_PI);
}

static double beta_logpdf(double x, double a, double b)
{
    if (!(x > 0.0) || !(x < 1.0)) return -HUGE_VAL;
    double lbeta = gd_lgamma(a) + gd_lgamma(b) - gd_lgamma(a + b);
    return (a - 1.0) * log(x) + (b - 1.0) * log1p(-x) - lbeta;
}

static double sigmoid(double z) { return 1.0 / (1.0 + exp(-z)); }
static double logit(double p)   { return log(p) - log1p(-p); }

/* ----------------------------------------------------------------- config */

void gd_mcmc_config_default(gd_mcmc_config *cfg)
{
    if (cfg == NULL) return;
    memset(cfg, 0, sizeof *cfg);
    cfg->seed = 0;
    cfg->n_warmup = 1000;
    cfg->n_samples = 1000;
    cfg->thin = 1;
    cfg->adapt = true;
    cfg->store_pointwise = false;
    cfg->threads = 0;
    cfg->sample_rate = true;
    cfg->fixed_rate = 1e-4;
    cfg->rate_logmean = -9.0;
    cfg->rate_logsd = 1.0;
    cfg->rate_step0 = 0.3;
    cfg->tree_prior = GD_TREE_PRIOR_UNIFORM;
    cfg->root_lex_logmean = -1.2;   /* exp(-1.2) = 0.30 replacements per concept */
    cfg->root_lex_logsd = 1.5;
    cfg->prior_only = false;
    cfg->age_step0_frac = 0.3;
    cfg->clock_move = true;
    cfg->shift_move = true;
    cfg->clock_step0 = 0.10;
    cfg->subtree_move = true;
    cfg->subtree_step0 = 0.10;
    cfg->block_age_scaling = false;
    cfg->block_step0 = 0.10;
    cfg->sample_branch = false;
    cfg->sigma_b_prior_sd = 0.30;
    cfg->sigma_b_init = 0.10;
    cfg->sigma_b_step0 = 0.30;
    cfg->branch_step0 = 0.20;
    cfg->sample_contamination = false;
    cfg->fixed_epsilon = 0.0;
    cfg->eps_beta_a = 1.0;
    cfg->eps_beta_b = 19.0;
    cfg->eps_init = 0.05;
    cfg->eps_step0 = 0.40;
    cfg->sample_rate_mix = false;
    cfg->rate_mix_prior_sd = 1.0;
    cfg->rate_mix_prior_logmean = 0.0;
    cfg->rate_mix_prior_logsd = 0.0;
    cfg->rate_mix_init = 0.5;
    cfg->rate_mix_step0 = 0.30;
    cfg->sample_concept_rates = false;
    cfg->concept_rate_sd = 0.50;
    cfg->concept_step0 = 0.30;
}

/* -------------------------------------------------------- state assembly */

/*
 * The statics below assert their handles rather than checking them. That is
 * the deliberate split (gd_error.h): a chain that exists has been through
 * gd_mcmc_chain_create, so a NULL model or parameter array here is a bug in
 * this file, not something a caller can provoke through the API. Bad data
 * from a user still comes back as a GD_E_* code, never as an assertion.
 */
static void fill_state(const gd_mcmc_chain *ch, gd_mcmc_state *st)
{
    assert(ch != NULL && st != NULL);
    /* Points straight at the chain's own storage: same double, no copy, and
     * the value the sampler is about to propose against. */
    st->rates = &ch->logr;
    st->n_rates = 1;
    st->concept_logdev = ch->have_concept ? ch->concept_dev : NULL;
    st->branch_logmult = ch->have_branch ? ch->branch_logmult : NULL;
    st->epsilon = ch->epsilon;
    st->rate_mix_sd = ch->sigma_c;
    st->rate_mix_family = ch->cfg.rate_mix_family;
    st->tip_hazard = ch->tip_hazard;
}

static int chain_eval_loglik(const gd_mcmc_chain *ch, double *pw, double *out,
                             gd_error *err)
{
    gd_mcmc_state st;
    assert(ch != NULL && ch->model != NULL && ch->model->loglik != NULL);
    if (ch->cfg.prior_only) {
        if (pw != NULL)
            for (size_t c = 0; c < ch->nconcepts; c++) pw[c] = 0.0;
        *out = 0.0;
        return 0;
    }
    fill_state(ch, &st);
    return ch->model->loglik(ch->model->ctx, ch->tree, &st, pw, out, err);
}

/* Rebuild the assembled buffers from the parameter values (after init/restore). */
static void rebuild_assembled(gd_mcmc_chain *ch)
{
    for (size_t p = 0; p < ch->nparam; p++) {
        const gd_param *pp = &ch->params[p];
        switch (pp->type) {
            case GD_MCMC_P_LOGRATE:      ch->logr = pp->value; break;
            case GD_MCMC_P_LOG_SIGMA_B:  ch->sigma_b = exp(pp->value); break;
            case GD_MCMC_P_LOGIT_EPS:    ch->epsilon = sigmoid(pp->value); break;
            case GD_MCMC_P_LOG_SIGMA_C:  ch->sigma_c = exp(pp->value); break;
            case GD_MCMC_P_LOGIT_TIPERR:
                ch->tip_hazard[pp->node] = -log1p(-sigmoid(pp->value)); break;
            case GD_MCMC_P_BRANCH_LOGMULT:
                ch->branch_logmult[pp->node] = pp->value; break;
            case GD_MCMC_P_CONCEPT_DEV:
                ch->concept_dev[pp->concept] = pp->value; break;
            case GD_MCMC_P_AGE:
                ch->tree->nodes[pp->node].age = pp->value; break;
        }
    }
}

/* ------------------------------------------------------------- log-prior */

static double node_logprior(const gd_calibration_set *cal, size_t v, double age)
{
    /* Sum the log-priors of every calibration bound to this node. A node may
     * carry more than one calibration (e.g. a soft lower bound and a normal
     * prior); each is an independent factor in the posterior. Returning only
     * the first would silently drop the others from the density while
     * propagation still applied all of them to the feasible bounds. Hard
     * constraints contribute -HUGE_VAL when violated, restricting support. */
    double s = 0.0;
    for (size_t c = 0; c < cal->ncals; c++) {
        if (cal->cals[c].resolved && cal->cals[c].node_index == v) {
            double lp = gd_cal_logprior(&cal->cals[c], age);
            if (lp == -HUGE_VAL) return -HUGE_VAL;
            s += lp;
        }
    }
    return s;
}

/* Beta(m k, (1 - m) k) on a tip's error probability, on the logit scale. */
static double tiperr_logprior(const gd_mcmc_chain *ch, size_t v, double z)
{
    double m = ch->cfg.tip_err_mean[v], k = ch->cfg.tip_err_conc[v];
    double e = sigmoid(z);
    return beta_logpdf(e, m * k, (1.0 - m) * k) + log(e) + log1p(-e);
}

/* Prior on log sigma_c (the sampled value): half-normal on sigma_c with its
 * Jacobian, or a normal on log sigma_c. */
static double sigc_logprior(const gd_mcmc_chain *ch, double logs)
{
    if (ch->cfg.rate_mix_prior_logsd > 0.0)
        return norm_logpdf(logs, ch->cfg.rate_mix_prior_logmean, ch->cfg.rate_mix_prior_logsd);
    return halfnormal_logpdf(exp(logs), ch->cfg.rate_mix_prior_sd) + logs;
}

/* The proper age prior at the chain's current root age and rate. */
static double tree_lp(const gd_mcmc_chain *ch)
{
    if (ch->tp == NULL) return 0.0;
    return gd_tree_prior_logpdf_tree(ch->tp, ch->logr);
}

/* Sum of branch-multiplier log-priors at a given sigma_b. */
static double branchmult_prior_sum(const gd_mcmc_chain *ch, double sigma_b)
{
    double s = 0.0;
    double mean = -0.5 * sigma_b * sigma_b;
    for (size_t p = 0; p < ch->nparam; p++)
        if (ch->params[p].type == GD_MCMC_P_BRANCH_LOGMULT)
            s += norm_logpdf(ch->params[p].value, mean, sigma_b);
    return s;
}

/* Log-prior of one parameter at its current value (with change-of-variable
 * Jacobians for the transformed scales). */
static double param_logprior(const gd_mcmc_chain *ch, const gd_param *pp)
{
    assert(ch != NULL && pp != NULL);
    switch (pp->type) {
        case GD_MCMC_P_LOGRATE:
            return norm_logpdf(pp->value, ch->cfg.rate_logmean, ch->cfg.rate_logsd);
        case GD_MCMC_P_AGE:
            return node_logprior(ch->cal, pp->node, ch->tree->nodes[pp->node].age);
        case GD_MCMC_P_LOG_SIGMA_B: {
            double s = exp(pp->value);
            return halfnormal_logpdf(s, ch->cfg.sigma_b_prior_sd) + pp->value;
        }
        case GD_MCMC_P_LOGIT_EPS: {
            double e = sigmoid(pp->value);
            return beta_logpdf(e, ch->cfg.eps_beta_a, ch->cfg.eps_beta_b)
                 + log(e) + log1p(-e);
        }
        case GD_MCMC_P_BRANCH_LOGMULT:
            return norm_logpdf(pp->value, -0.5 * ch->sigma_b * ch->sigma_b,
                               ch->sigma_b);
        case GD_MCMC_P_CONCEPT_DEV:
            return norm_logpdf(pp->value, 0.0, ch->cfg.concept_rate_sd);
        case GD_MCMC_P_LOG_SIGMA_C:
            return sigc_logprior(ch, pp->value);
        case GD_MCMC_P_LOGIT_TIPERR:
            return tiperr_logprior(ch, pp->node, pp->value);
    }
    return 0.0;
}

static double chain_total_logprior(const gd_mcmc_chain *ch)
{
    double lp = 0.0;
    for (size_t p = 0; p < ch->nparam; p++)
        lp += param_logprior(ch, &ch->params[p]);
    return lp + tree_lp(ch);
}

/* -------------------------------------------------------------- lifecycle */

void gd_mcmc_chain_destroy(gd_mcmc_chain *ch)
{
    if (ch == NULL) return;
    if (ch->rng != NULL) gd_rng_destroy(ch->rng);
    gd_tree_prior_destroy(ch->tp);
    free(ch->params);
    free(ch->age_params);
    free(ch->scaled);
    free(ch->sub_root); free(ch->sub_start); free(ch->sub_nodes); free(ch->sub_scratch);
    free(ch->tip_hazard);
    free(ch->age_scratch);
    free(ch->funnel_scratch);
    free(ch->branch_logmult);
    free(ch->concept_dev);
    free(ch->samples); free(ch->sample_ll); free(ch->pointwise);
    free(ch->pw_scratch);
    free(ch);
}

/* Count the parameters implied by the config and tree. */
static size_t count_params(const gd_mcmc_chain *ch, size_t *nfree_out)
{
    size_t nfree = 0;
    for (size_t i = 0; i < ch->tree->nnodes; i++)
        if (ch->cal->node_hi[i] - ch->cal->node_lo[i] > GD_MCMC_FIX_TOL) nfree++;
    if (nfree_out) *nfree_out = nfree;
    size_t np = nfree;
    if (ch->cfg.sample_rate) np += 1;
    if (ch->cfg.sample_branch) np += (ch->tree->nnodes - 1) + 1; /* branches + sigma_b */
    if (ch->cfg.sample_contamination) np += 1;
    if (ch->cfg.sample_concept_rates) np += ch->nconcepts;
    if (ch->cfg.sample_rate_mix) np += 1;
    if (ch->cfg.tip_err_mean != NULL)
        for (size_t v = 0; v < ch->tree->nnodes; v++)
            if (ch->cfg.tip_err_mean[v] > 0.0) np += 1;
    return np;
}

static uint64_t chain_binding(const gd_mcmc_chain *ch);

gd_mcmc_chain *gd_mcmc_chain_create(const gd_mcmc_model *model, gd_tree *tree,
                                    gd_calibration_set *cal,
                                    const gd_mcmc_config *cfg, uint64_t stream,
                                    gd_error *err)
{
    if (model == NULL || model->loglik == NULL || tree == NULL || cal == NULL ||
        cfg == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "mcmc: NULL argument");
        return NULL;
    }
    if (cal->node_lo == NULL || cal->node_hi == NULL || cal->nnodes != tree->nnodes) {
        gd_error_set(err, GD_E_INVALID_ARG,
                     "mcmc: calibration bounds not propagated for this tree");
        return NULL;
    }
    if (cfg->thin == 0) { gd_error_set(err, GD_E_INVALID_ARG, "mcmc: thin >= 1"); return NULL; }
    if (cfg->sample_rate && !(cfg->rate_logsd > 0.0)) {
        gd_error_set(err, GD_E_RANGE, "mcmc: rate_logsd > 0 required"); return NULL;
    }
    if (!cfg->sample_rate && !(cfg->fixed_rate > 0.0)) {
        gd_error_set(err, GD_E_RANGE, "mcmc: fixed_rate > 0 required"); return NULL;
    }
    if (cfg->sample_rate_mix && !(cfg->rate_mix_prior_sd > 0.0 && cfg->rate_mix_init > 0.0)) {
        gd_error_set(err, GD_E_RANGE, "mcmc: rate_mix_prior_sd and rate_mix_init > 0 required"); return NULL;
    }
    if (cfg->sample_branch && !(cfg->sigma_b_prior_sd > 0.0)) {
        gd_error_set(err, GD_E_RANGE, "mcmc: sigma_b_prior_sd > 0 required"); return NULL;
    }

    gd_mcmc_chain *ch = (gd_mcmc_chain *)gd_calloc_checked(1, sizeof(*ch), err);
    if (ch == NULL) { gd_error_set(err, GD_E_NOMEM, "mcmc: chain alloc"); return NULL; }
    ch->model = model;
    ch->tree = tree;
    ch->cal = cal;
    ch->cfg = *cfg;
    ch->nnodes = tree->nnodes;
    ch->nconcepts = model->n_pointwise;
    ch->have_branch = cfg->sample_branch ? 1 : 0;
    ch->have_concept = cfg->sample_concept_rates ? 1 : 0;
    ch->cap = cfg->n_samples;
    ch->sigma_b = cfg->sample_branch ? cfg->sigma_b_init : 0.0;
    ch->epsilon = cfg->sample_contamination ? cfg->eps_init : cfg->fixed_epsilon;
    ch->logr = cfg->sample_rate ? cfg->rate_logmean : log(cfg->fixed_rate);

    size_t nfree = 0;
    ch->nparam = count_params(ch, &nfree);

    ch->stream = stream;
    ch->rng = gd_rng_create(cfg->seed, stream, err);
    ch->params = (gd_param *)gd_calloc_checked(ch->nparam ? ch->nparam : 1, sizeof(gd_param), err);
    /* Retained draws x parameters: the first factor is a run setting, the
     * second comes from the tree, so the product is checked before it becomes
     * a size. */
    size_t cap1 = ch->cap ? ch->cap : 1;
    size_t nsamp;
    if (!gd_size_mul(cap1, ch->nparam ? ch->nparam : 1, &nsamp)) {
        gd_error_set(err, GD_E_NOMEM, "chain: draw storage size overflows");
        gd_mcmc_chain_destroy(ch); return NULL;
    }
    ch->samples = (double *)gd_malloc_array_checked(nsamp, sizeof(double), err);
    ch->sample_ll = (double *)gd_malloc_array_checked(cap1, sizeof(double), err);
    if (ch->have_branch)
        ch->branch_logmult = (double *)gd_calloc_checked(ch->nnodes, sizeof(double), err);
    if (ch->have_concept && ch->nconcepts > 0)
        ch->concept_dev = (double *)gd_calloc_checked(ch->nconcepts, sizeof(double), err);
    if (ch->nconcepts > 0) {
        ch->pw_scratch = (double *)gd_malloc_array_checked(ch->nconcepts, sizeof(double), err);
        if (cfg->store_pointwise) {
            size_t npw;
            if (!gd_size_mul(cap1, ch->nconcepts, &npw)) {
                gd_error_set(err, GD_E_NOMEM, "chain: pointwise size overflows");
                gd_mcmc_chain_destroy(ch); return NULL;
            }
            ch->pointwise = (double *)gd_malloc_array_checked(npw, sizeof(double), err);
        }
    }
    if (ch->rng == NULL || ch->params == NULL || ch->samples == NULL ||
        ch->sample_ll == NULL ||
        (ch->have_branch && ch->branch_logmult == NULL) ||
        (ch->have_concept && ch->nconcepts > 0 && ch->concept_dev == NULL) ||
        (ch->nconcepts > 0 && ch->pw_scratch == NULL) ||
        (cfg->store_pointwise && ch->nconcepts > 0 && ch->pointwise == NULL)) {
        gd_error_set(err, GD_E_NOMEM, "mcmc: chain arrays");
        gd_mcmc_chain_destroy(ch);
        return NULL;
    }

    /* The calibration set is the caller's: a chain reads it and never marks
     * it feasible on the caller's behalf (it used to, which let an
     * unpropagated or contradictory set through to sampling). */
    if (!cal->feasible) {
        gd_error_set(err, GD_E_CALIBRATION,
                     "mcmc: calibrations are not marked feasible; propagate them first");
        gd_mcmc_chain_destroy(ch); return NULL;
    }

    /* Initial rate, drawn first: the initial root is drawn through it. */
    if (cfg->sample_rate)
        ch->logr = cfg->rate_logmean + gd_rng_normal(ch->rng, 0.0, cfg->rate_logsd);

    /* The age prior. */
    if (cfg->tree_prior == GD_TREE_PRIOR_UNIFORM || cfg->tree_prior == GD_TREE_PRIOR_YULE) {
        bool *is_free = (bool *)gd_malloc_array_checked(tree->nnodes, sizeof(bool), err);
        double *fixed = (double *)gd_malloc_array_checked(tree->nnodes, sizeof(double), err);
        if (is_free == NULL || fixed == NULL) {
            free(is_free); free(fixed);
            gd_error_set(err, GD_E_NOMEM, "mcmc: tree prior"); gd_mcmc_chain_destroy(ch); return NULL;
        }
        for (size_t v = 0; v < tree->nnodes; v++) {
            is_free[v] = (cal->node_hi[v] - cal->node_lo[v] > GD_MCMC_FIX_TOL);
            fixed[v] = cal->node_lo[v];
        }
        ch->tp = gd_tree_prior_create(tree, is_free, fixed, cfg->root_lex_logmean,
                                      cfg->root_lex_logsd, err);
        free(is_free); free(fixed);
        if (ch->tp == NULL) { gd_mcmc_chain_destroy(ch); return NULL; }
        if (cfg->tree_prior == GD_TREE_PRIOR_YULE) {
            size_t ntip = 0;
            for (size_t v = 0; v < tree->nnodes; v++) ntip += tree->nodes[v].is_leaf ? 1u : 0u;
            double kappa = (cfg->yule_kappa > 0.0) ? cfg->yule_kappa
                                                   : gd_tree_prior_default_kappa(ntip);
            if (gd_tree_prior_set_yule(ch->tp, kappa, err) != 0) {
                gd_mcmc_chain_destroy(ch); return NULL;
            }
        }
    } else {
        /* The flat measure is proper only inside bounded calibration windows. */
        for (size_t v = 0; v < tree->nnodes; v++) {
            if (cal->node_hi[v] - cal->node_lo[v] > GD_MCMC_FIX_TOL && !isfinite(cal->node_hi[v])) {
                gd_error_set(err, GD_E_RANGE,
                             "mcmc: tree prior 'none' needs every free age bounded above; "
                             "the posterior would be improper");
                gd_mcmc_chain_destroy(ch); return NULL;
            }
        }
    }

    /* Initial ages, root first. Each free age is drawn inside the window its
     * parent and calibrations leave it; the root, under the proper prior, as
     * L / rho with L from the lexical-root prior. Different streams start
     * different chains in different places. */
    for (size_t i = 0; i < tree->nnodes; i++) {
        size_t v = tree->preorder[i];
        double lo = cal->node_lo[v], hi = cal->node_hi[v];
        if (!(hi - lo > GD_MCMC_FIX_TOL)) { tree->nodes[v].age = lo; continue; }
        size_t par = tree->nodes[v].parent;
        if (v != tree->root && par != GD_NODE_NONE && tree->nodes[par].age < hi)
            hi = tree->nodes[par].age;
        double a;
        if (v == tree->root && ch->tp != NULL) {
            double L = exp(cfg->root_lex_logmean
                           + gd_rng_normal(ch->rng, 0.0, cfg->root_lex_logsd));
            a = L / exp(ch->logr);
            /* Keep it strictly inside the window; the likelihood decides the
             * rest during warm-up. */
            if (!(a > lo)) a = lo + (isfinite(hi) ? 0.5 * (hi - lo) : (lo > 0.0 ? lo : 1.0));
            if (!(a < hi)) a = lo + (hi - lo) * (0.25 + 0.5 * gd_rng_uniform_closed_open(ch->rng));
        } else {
            a = lo + (hi - lo) * (0.05 + 0.9 * gd_rng_uniform_closed_open(ch->rng));
        }
        tree->nodes[v].age = a;
    }

    /* Lay out the parameter vector and initialise values + steps. */
    size_t p = 0;
    if (cfg->sample_rate) {
        ch->params[p].type = GD_MCMC_P_LOGRATE;
        ch->params[p].node = GD_NODE_NONE;
        ch->params[p].value = ch->logr;
        ch->params[p].step = cfg->rate_step0;
        p++;
    }
    if (cfg->sample_branch) {
        /* sigma_b first, then one multiplier per non-root node. */
        double s = fabs(gd_rng_normal(ch->rng, 0.0, cfg->sigma_b_prior_sd));
        if (s < 1e-3) s = 1e-3;
        ch->sigma_b = s;
        ch->params[p].type = GD_MCMC_P_LOG_SIGMA_B;
        ch->params[p].node = GD_NODE_NONE;
        ch->params[p].value = log(s);
        ch->params[p].step = cfg->sigma_b_step0;
        p++;
        for (size_t v = 0; v < tree->nnodes; v++) {
            if (v == tree->root) continue;
            double lm = gd_rng_normal(ch->rng, -0.5 * s * s, s);
            ch->branch_logmult[v] = lm;
            ch->params[p].type = GD_MCMC_P_BRANCH_LOGMULT;
            ch->params[p].node = v;
            ch->params[p].value = lm;
            ch->params[p].step = cfg->branch_step0;
            p++;
        }
    }
    if (cfg->tip_err_mean != NULL) {
        ch->tip_hazard = (double *)gd_calloc_checked(tree->nnodes, sizeof(double), err);
        if (ch->tip_hazard == NULL) { gd_error_set(err, GD_E_NOMEM, "mcmc: tip error"); gd_mcmc_chain_destroy(ch); return NULL; }
        for (size_t v = 0; v < tree->nnodes; v++) {
            double m = cfg->tip_err_mean[v];
            if (!(m > 0.0)) continue;
            double k = (cfg->tip_err_conc != NULL) ? cfg->tip_err_conc[v] : 0.0;
            if (!(m < 1.0) || !(k > 0.0) || !tree->nodes[v].is_leaf) {
                gd_error_set(err, GD_E_RANGE,
                             "mcmc: a tip error prior needs a tip, a mean in (0, 1) and a concentration > 0");
                gd_mcmc_chain_destroy(ch); return NULL;
            }
            ch->params[p].type = GD_MCMC_P_LOGIT_TIPERR;
            ch->params[p].node = v;
            ch->params[p].value = logit(m);
            ch->params[p].step = 0.5;
            ch->tip_hazard[v] = -log1p(-m);
            p++;
        }
    }
    if (cfg->sample_rate_mix) {
        ch->sigma_c = cfg->rate_mix_init;
        ch->params[p].type = GD_MCMC_P_LOG_SIGMA_C;
        ch->params[p].node = GD_NODE_NONE;
        ch->params[p].value = log(ch->sigma_c);
        ch->params[p].step = cfg->rate_mix_step0;
        p++;
    }
    if (cfg->sample_contamination) {
        ch->epsilon = cfg->eps_init;
        ch->params[p].type = GD_MCMC_P_LOGIT_EPS;
        ch->params[p].node = GD_NODE_NONE;
        ch->params[p].value = logit(cfg->eps_init);
        ch->params[p].step = cfg->eps_step0;
        p++;
    }
    for (size_t i = 0; i < tree->nnodes; i++) {
        if (cal->node_hi[i] - cal->node_lo[i] > GD_MCMC_FIX_TOL) {
            ch->params[p].type = GD_MCMC_P_AGE;
            ch->params[p].node = i;
            ch->params[p].value = tree->nodes[i].age;
            ch->params[p].step = cfg->age_step0_frac;   /* dimensionless; see chain_sweep */
            p++;
        }
    }
    if (cfg->sample_concept_rates) {
        for (size_t c = 0; c < ch->nconcepts; c++) {
            double dev = gd_rng_normal(ch->rng, 0.0, cfg->concept_rate_sd);
            ch->concept_dev[c] = dev;
            ch->params[p].type = GD_MCMC_P_CONCEPT_DEV;
            ch->params[p].node = GD_NODE_NONE;
            ch->params[p].concept = c;
            ch->params[p].value = dev;
            ch->params[p].step = cfg->concept_step0;
            p++;
        }
    }

    /* Collect the ages the scale moves act on: free internal nodes only. A
     * free tip (an old list dated to a narrow window) would leave its window
     * under any useful scale factor and get every scale move rejected; it is
     * left to its own single-site updates. */
    ch->n_age = 0;
    ch->scaled = (bool *)gd_calloc_checked(tree->nnodes ? tree->nnodes : 1, sizeof(bool), err);
    if (ch->scaled == NULL) { gd_error_set(err, GD_E_NOMEM, "mcmc: scale set"); gd_mcmc_chain_destroy(ch); return NULL; }
    for (size_t q = 0; q < ch->nparam; q++)
        if (ch->params[q].type == GD_MCMC_P_AGE && !tree->nodes[ch->params[q].node].is_leaf) {
            ch->n_age++;
            ch->scaled[ch->params[q].node] = true;
        }
    ch->have_block = (cfg->block_age_scaling && ch->n_age >= 2) ? 1 : 0;
    ch->block_step = cfg->block_step0;
    ch->have_clock = (cfg->clock_move && cfg->sample_rate && ch->n_age >= 1) ? 1 : 0;
    ch->clock_step = cfg->clock_step0;

    /* Subtree movers: free internal non-root nodes with at least one free
     * internal descendant (a lone node gains nothing over its own move). */
    ch->sub_step = cfg->subtree_step0;
    if (cfg->subtree_move) {
        size_t n = tree->nnodes;
        ch->sub_root = (size_t *)gd_malloc_array_checked(n, sizeof(size_t), err);
        ch->sub_start = (size_t *)gd_malloc_array_checked(n + 1, sizeof(size_t), err);
        ch->sub_nodes = (size_t *)gd_malloc_array_checked(n * n, sizeof(size_t), err);
        ch->sub_scratch = (double *)gd_malloc_array_checked(n, sizeof(double), err);
        if (ch->sub_root == NULL || ch->sub_start == NULL || ch->sub_nodes == NULL ||
            ch->sub_scratch == NULL) {
            gd_error_set(err, GD_E_NOMEM, "mcmc: subtree move"); gd_mcmc_chain_destroy(ch); return NULL;
        }
        size_t used = 0;
        ch->sub_start[0] = 0;
        for (size_t v = 0; v < n; v++) {
            if (v == tree->root || !ch->scaled[v]) continue;
            size_t begin = used;
            /* Descendants of v (v included) that the scale moves act on. */
            for (size_t u = 0; u < n; u++) {
                if (!ch->scaled[u]) continue;
                size_t w = u;
                while (w != GD_NODE_NONE && w != v) w = tree->nodes[w].parent;
                if (w == v) ch->sub_nodes[used++] = u;
            }
            if (used - begin >= 2) {
                ch->sub_root[ch->n_sub] = v;
                ch->sub_start[++ch->n_sub] = used;
            } else {
                used = begin;
            }
        }
    }
    if (ch->n_age > 0) {
        ch->age_params = (size_t *)gd_malloc_array_checked(ch->n_age, sizeof(size_t), err);
        ch->age_scratch = (double *)gd_malloc_array_checked(ch->n_age, sizeof(double), err);
        if (ch->age_params == NULL || ch->age_scratch == NULL) {
            gd_error_set(err, GD_E_NOMEM, "mcmc: block-move scratch");
            gd_mcmc_chain_destroy(ch); return NULL;
        }
        size_t a = 0;
        for (size_t q = 0; q < ch->nparam; q++)
            if (ch->params[q].type == GD_MCMC_P_AGE && ch->scaled[ch->params[q].node])
                ch->age_params[a++] = q;
    }

    ch->funnel_param = ch->nparam;
    ch->funnel_step = 0.30;
    for (size_t q = 0; q < ch->nparam; q++)
        if (ch->params[q].type == GD_MCMC_P_LOG_SIGMA_B) ch->funnel_param = q;
    if (ch->funnel_param < ch->nparam) {
        ch->funnel_scratch = (double *)gd_malloc_array_checked(ch->nparam, sizeof(double), err);
        if (ch->funnel_scratch == NULL) { gd_mcmc_chain_destroy(ch); return NULL; }
    }
    ch->shift_param = ch->nparam;
    ch->shift_step = 0.10;
    if (cfg->shift_move && ch->funnel_param < ch->nparam)
        for (size_t q = 0; q < ch->nparam; q++)
            if (ch->params[q].type == GD_MCMC_P_LOGRATE) ch->shift_param = q;

    if (chain_eval_loglik(ch, NULL, &ch->cur_loglik, err) != 0) {
        gd_mcmc_chain_destroy(ch); return NULL;
    }
    ch->cur_tp = tree_lp(ch);
    ch->cur_logprior = chain_total_logprior(ch);
    /* A chain that starts at zero density cannot be trusted to find its way:
     * every proposal compares against -inf (or NaN), and the first accepted
     * move is arbitrary. Refuse rather than sample from there. */
    if (!isfinite(ch->cur_loglik) || !isfinite(ch->cur_logprior)) {
        gd_error_setf(err, GD_E_RANGE,
                      "mcmc: initial state has non-finite %s (log-likelihood %g, log-prior %g)",
                      !isfinite(ch->cur_loglik) ? "likelihood" : "prior",
                      ch->cur_loglik, ch->cur_logprior);
        gd_mcmc_chain_destroy(ch); return NULL;
    }
    ch->binding = chain_binding(ch);
    return ch;
}

/* The inputs a checkpoint is bound to (E3.3): the parameter layout, the tree,
 * the calibration windows, the model's data, every setting that shapes the
 * trajectory, and the chain's stream. Not the thread count or progress
 * settings, which change nothing a chain draws. */
static uint64_t chain_binding(const gd_mcmc_chain *ch)
{
    uint64_t h = 1469598103934665603ULL;
#define BIND(ptr, n) do { const unsigned char *b_ = (const unsigned char *)(ptr); \
        for (size_t i_ = 0; i_ < (n); i_++) { h ^= b_[i_]; h *= 1099511628211ULL; } } while (0)
#define BIND_D(x) do { double d_ = (double)(x); BIND(&d_, sizeof d_); } while (0)
    const gd_mcmc_config *c = &ch->cfg;
    BIND_D(ch->stream); BIND_D(c->seed);
    BIND_D(c->n_warmup); BIND_D(c->n_samples); BIND_D(c->thin);
    BIND_D(c->adapt); BIND_D(c->store_pointwise);
    BIND_D(c->sample_rate); BIND_D(c->fixed_rate);
    BIND_D(c->rate_logmean); BIND_D(c->rate_logsd); BIND_D(c->rate_step0);
    BIND_D(c->tree_prior); BIND_D(c->root_lex_logmean); BIND_D(c->root_lex_logsd);
    BIND_D(c->yule_kappa); BIND_D(c->prior_only); BIND_D(c->age_step0_frac);
    BIND_D(c->clock_move); BIND_D(c->clock_step0); BIND_D(c->shift_move);
    BIND_D(c->subtree_move); BIND_D(c->subtree_step0);
    BIND_D(c->block_age_scaling); BIND_D(c->block_step0);
    BIND_D(c->sample_branch); BIND_D(c->sigma_b_prior_sd); BIND_D(c->sigma_b_init);
    BIND_D(c->sigma_b_step0); BIND_D(c->branch_step0);
    BIND_D(c->sample_contamination); BIND_D(c->fixed_epsilon);
    BIND_D(c->eps_beta_a); BIND_D(c->eps_beta_b); BIND_D(c->eps_init); BIND_D(c->eps_step0);
    BIND_D(c->sample_rate_mix); BIND_D(c->rate_mix_family); BIND_D(c->rate_mix_prior_sd);
    BIND_D(c->rate_mix_prior_logmean); BIND_D(c->rate_mix_prior_logsd);
    BIND_D(c->rate_mix_init); BIND_D(c->rate_mix_step0);
    BIND_D(c->sample_concept_rates); BIND_D(c->concept_rate_sd); BIND_D(c->concept_step0);
    for (size_t v = 0; v < ch->tree->nnodes; v++) {
        if (c->tip_err_mean != NULL) BIND_D(c->tip_err_mean[v]);
        if (c->tip_err_conc != NULL) BIND_D(c->tip_err_conc[v]);
    }
    const gd_tree *t = ch->tree;
    BIND_D(t->nnodes); BIND_D(t->root);
    for (size_t v = 0; v < t->nnodes; v++) {
        BIND_D(t->nodes[v].parent); BIND_D(t->nodes[v].is_leaf);
        if (t->nodes[v].label != NULL) BIND(t->nodes[v].label, strlen(t->nodes[v].label) + 1);
        BIND_D(ch->cal->node_lo[v]); BIND_D(ch->cal->node_hi[v]);
    }
    BIND_D(ch->nparam);
    for (size_t p = 0; p < ch->nparam; p++) {
        BIND_D(ch->params[p].type); BIND_D(ch->params[p].node);
    }
    if (ch->model->fingerprint != NULL) {
        uint64_t f = ch->model->fingerprint(ch->model->ctx);
        BIND(&f, sizeof f);
    }
#undef BIND_D
#undef BIND
    return h;
}

/* --------------------------------------------------------------- sweep */

static void adapt_step(gd_param *pp, int accepted)
{
    pp->adapt_count += 1.0;
    double gamma = 1.0 / pow(pp->adapt_count + GD_MCMC_ADAPT_T0, 0.6);
    double logs = log(pp->step) + gamma * ((double)accepted - GD_MCMC_TARGET_ACCEPT);
    if (logs < -30.0) logs = -30.0;
    if (logs > 30.0) logs = 30.0;
    pp->step = exp(logs);
}

/* Metropolis accept test with a fresh uniform. */
static int accept_move(gd_mcmc_chain *ch, double log_alpha)
{
    return log(gd_rng_uniform_open(ch->rng)) < log_alpha;
}

/*
 * Scale every free age by s = exp(delta). With `clock` set the log-rate moves
 * by -delta at the same time, so rho * t is unchanged on every branch: a move
 * along the rate-age ridge that the likelihood barely sees. Rejected if any
 * order or calibration window is violated. The Jacobian of scaling m ages is
 * s^m (the log-rate shift has unit Jacobian), so the ratio carries m * delta.
 */
static int scale_move(gd_mcmc_chain *ch, int adapting, int clock, gd_error *err)
{
    double *step = clock ? &ch->clock_step : &ch->block_step;
    double *adapt_count = clock ? &ch->clock_adapt_count : &ch->block_adapt_count;
    double delta = gd_rng_normal(ch->rng, 0.0, *step);
    double s = exp(delta);
    bool feasible = true;
    for (size_t a = 0; a < ch->n_age; a++)
        ch->age_scratch[a] = ch->tree->nodes[ch->params[ch->age_params[a]].node].age;
    for (size_t a = 0; a < ch->n_age && feasible; a++) {
        size_t v = ch->params[ch->age_params[a]].node;
        gd_node *node = &ch->tree->nodes[v];
        double na = ch->age_scratch[a] * s;
        double lo = ch->cal->node_lo[v], hi = ch->cal->node_hi[v];
        for (size_t k = 0; k < node->nchildren; k++) {
            size_t c = node->children[k];
            double cage = ch->scaled[c] ? ch->tree->nodes[c].age * s : ch->tree->nodes[c].age;
            if (cage > lo) lo = cage;
        }
        if (v != ch->tree->root) {
            size_t par = node->parent;
            double page = ch->scaled[par] ? ch->tree->nodes[par].age * s : ch->tree->nodes[par].age;
            if (page < hi) hi = page;
        }
        if (!(na > lo && na < hi)) feasible = false;
    }
    int accepted = 0;
    if (feasible) {
        double dprior = 0.0;
        for (size_t a = 0; a < ch->n_age; a++) {
            size_t v = ch->params[ch->age_params[a]].node;
            double np = ch->age_scratch[a] * s;
            dprior += node_logprior(ch->cal, v, np)
                    - node_logprior(ch->cal, v, ch->age_scratch[a]);
            ch->tree->nodes[v].age = np;
        }
        double old_logr = ch->logr;
        if (clock) {
            ch->logr = old_logr - delta;
            dprior += norm_logpdf(ch->logr, ch->cfg.rate_logmean, ch->cfg.rate_logsd)
                    - norm_logpdf(old_logr, ch->cfg.rate_logmean, ch->cfg.rate_logsd);
        }
        double new_tp = tree_lp(ch);
        dprior += new_tp - ch->cur_tp;
        double new_ll;
        if (chain_eval_loglik(ch, NULL, &new_ll, err) != 0) {
            ch->logr = old_logr;
            for (size_t a = 0; a < ch->n_age; a++)
                ch->tree->nodes[ch->params[ch->age_params[a]].node].age = ch->age_scratch[a];
            return GD_E_INTERNAL;
        }
        double log_alpha = (new_ll - ch->cur_loglik) + dprior
                         + (double)ch->n_age * delta;   /* Jacobian s^m */
        if (accept_move(ch, log_alpha)) {
            for (size_t a = 0; a < ch->n_age; a++)
                ch->params[ch->age_params[a]].value =
                    ch->tree->nodes[ch->params[ch->age_params[a]].node].age;
            if (clock) {
                for (size_t p = 0; p < ch->nparam; p++)
                    if (ch->params[p].type == GD_MCMC_P_LOGRATE) ch->params[p].value = ch->logr;
            }
            ch->cur_loglik = new_ll;
            ch->cur_logprior += dprior;
            ch->cur_tp = new_tp;
            accepted = 1;
        } else {
            ch->logr = old_logr;
            for (size_t a = 0; a < ch->n_age; a++)   /* exact restore */
                ch->tree->nodes[ch->params[ch->age_params[a]].node].age =
                    ch->age_scratch[a];
        }
    }
    if (clock) { ch->clock_prop += 1.0; ch->clock_acc += (double)accepted; }
    else       { ch->block_prop += 1.0; ch->block_acc += (double)accepted; }
    if (adapting) {
        *adapt_count += 1.0;
        double gamma = 1.0 / pow(*adapt_count + GD_MCMC_ADAPT_T0, 0.6);
        double logs = log(*step) + gamma * ((double)accepted - GD_MCMC_TARGET_ACCEPT);
        if (logs < -30.0) logs = -30.0;
        if (logs > 30.0) logs = 30.0;
        *step = exp(logs);
    }
    return 0;
}

/*
 * The relaxed clock's funnel: branch log-multipliers b_i ~ N(-s^2/2, s), so as
 * s shrinks they are pinned together and single-site moves on s and on each
 * b_i stall (Neal's funnel). This move changes log s by delta and maps every
 * b_i through its standardised deviation z_i = (b_i + s^2/2) / s, which it
 * keeps: b_i' = -s'^2/2 + s' z_i. For fixed delta the map is triangular with
 * diagonal (1, e^delta, ..., e^delta), so the Jacobian is n_b * delta.
 */
static int funnel_move(gd_mcmc_chain *ch, int adapting, gd_error *err)
{
    gd_param *ps = &ch->params[ch->funnel_param];
    double delta = gd_rng_normal(ch->rng, 0.0, ch->funnel_step);
    double old_ls = ps->value, new_ls = old_ls + delta;
    double s0 = exp(old_ls), s1 = exp(new_ls);
    double m0 = -0.5 * s0 * s0, m1 = -0.5 * s1 * s1;
    double dprior = (halfnormal_logpdf(s1, ch->cfg.sigma_b_prior_sd) + new_ls)
                  - (halfnormal_logpdf(s0, ch->cfg.sigma_b_prior_sd) + old_ls);
    size_t nb = 0;
    for (size_t p = 0; p < ch->nparam; p++) {
        gd_param *pp = &ch->params[p];
        if (pp->type != GD_MCMC_P_BRANCH_LOGMULT) continue;
        double b0 = pp->value, b1 = m1 + s1 * (b0 - m0) / s0;
        ch->funnel_scratch[p] = b0;
        ch->branch_logmult[pp->node] = b1;
        dprior += norm_logpdf(b1, m1, s1) - norm_logpdf(b0, m0, s0);
        nb++;
    }
    double new_ll;
    int accepted = 0;
    if (chain_eval_loglik(ch, NULL, &new_ll, err) != 0) {
        for (size_t p = 0; p < ch->nparam; p++)
            if (ch->params[p].type == GD_MCMC_P_BRANCH_LOGMULT)
                ch->branch_logmult[ch->params[p].node] = ch->funnel_scratch[p];
        return GD_E_INTERNAL;
    }
    if (accept_move(ch, (new_ll - ch->cur_loglik) + dprior + (double)nb * delta)) {
        ps->value = new_ls; ch->sigma_b = s1;
        for (size_t p = 0; p < ch->nparam; p++)
            if (ch->params[p].type == GD_MCMC_P_BRANCH_LOGMULT)
                ch->params[p].value = ch->branch_logmult[ch->params[p].node];
        ch->cur_loglik = new_ll;
        ch->cur_logprior += dprior;
        accepted = 1;
    } else {
        for (size_t p = 0; p < ch->nparam; p++)
            if (ch->params[p].type == GD_MCMC_P_BRANCH_LOGMULT)
                ch->branch_logmult[ch->params[p].node] = ch->funnel_scratch[p];
    }
    ch->funnel_prop += 1.0; ch->funnel_acc += (double)accepted;
    if (adapting) {
        ch->funnel_adapt_count += 1.0;
        double gamma = 1.0 / pow(ch->funnel_adapt_count + GD_MCMC_ADAPT_T0, 0.6);
        double logs = log(ch->funnel_step) + gamma * ((double)accepted - GD_MCMC_TARGET_ACCEPT);
        if (logs < -30.0) logs = -30.0;
        if (logs > 30.0) logs = 30.0;
        ch->funnel_step = exp(logs);
    }
    return 0;
}

/*
 * The relaxed clock's other ridge: a branch's rate is exp(logr + b_i), so
 * logr + delta with every b_i - delta leaves the likelihood unchanged, and
 * only the priors (the clock, the branch multipliers, and the lexical-root
 * age prior through rho) see the move. Single-site updates cross this ridge
 * one coordinate at a time. A translation: unit Jacobian. The likelihood is
 * still evaluated, so the move stays exact if that invariance ever breaks.
 */
static int shift_move(gd_mcmc_chain *ch, int adapting, gd_error *err)
{
    gd_param *pr = &ch->params[ch->shift_param];
    double delta = gd_rng_normal(ch->rng, 0.0, ch->shift_step);
    double old_logr = ch->logr, s = ch->sigma_b, m = -0.5 * s * s;
    double dprior = norm_logpdf(old_logr + delta, ch->cfg.rate_logmean, ch->cfg.rate_logsd)
                  - norm_logpdf(old_logr, ch->cfg.rate_logmean, ch->cfg.rate_logsd);
    ch->logr = old_logr + delta;
    for (size_t p = 0; p < ch->nparam; p++) {
        gd_param *pp = &ch->params[p];
        if (pp->type != GD_MCMC_P_BRANCH_LOGMULT) continue;
        double b0 = pp->value, b1 = b0 - delta;
        ch->funnel_scratch[p] = b0;
        ch->branch_logmult[pp->node] = b1;
        dprior += norm_logpdf(b1, m, s) - norm_logpdf(b0, m, s);
    }
    double new_tp = tree_lp(ch);
    dprior += new_tp - ch->cur_tp;
    double new_ll;
    int accepted = 0, rc = chain_eval_loglik(ch, NULL, &new_ll, err);
    if (rc == 0 && accept_move(ch, (new_ll - ch->cur_loglik) + dprior)) {
        pr->value = ch->logr;
        for (size_t p = 0; p < ch->nparam; p++)
            if (ch->params[p].type == GD_MCMC_P_BRANCH_LOGMULT)
                ch->params[p].value = ch->branch_logmult[ch->params[p].node];
        ch->cur_loglik = new_ll;
        ch->cur_logprior += dprior;
        ch->cur_tp = new_tp;
        accepted = 1;
    } else {
        ch->logr = old_logr;
        for (size_t p = 0; p < ch->nparam; p++)
            if (ch->params[p].type == GD_MCMC_P_BRANCH_LOGMULT)
                ch->branch_logmult[ch->params[p].node] = ch->funnel_scratch[p];
        if (rc != 0) return GD_E_INTERNAL;
    }
    ch->shift_prop += 1.0; ch->shift_acc += (double)accepted;
    if (adapting) {
        ch->shift_adapt_count += 1.0;
        double gamma = 1.0 / pow(ch->shift_adapt_count + GD_MCMC_ADAPT_T0, 0.6);
        double logs = log(ch->shift_step) + gamma * ((double)accepted - GD_MCMC_TARGET_ACCEPT);
        if (logs < -30.0) logs = -30.0;
        if (logs > 30.0) logs = 30.0;
        ch->shift_step = exp(logs);
    }
    return 0;
}

/* Scale the free internal ages under one mover by exp(delta); see
 * gd_mcmc_config.subtree_move. Jacobian s^k for k scaled ages. */
static int subtree_move(gd_mcmc_chain *ch, size_t k, int adapting, gd_error *err)
{
    size_t a0 = ch->sub_start[k], a1 = ch->sub_start[k + 1];
    double delta = gd_rng_normal(ch->rng, 0.0, ch->sub_step);
    double s = exp(delta);
    bool feasible = true;
    for (size_t i = a0; i < a1; i++) ch->sub_scratch[i - a0] = ch->tree->nodes[ch->sub_nodes[i]].age;
    /* Membership test for "is this node scaled in this move". */
    for (size_t i = a0; i < a1 && feasible; i++) {
        size_t v = ch->sub_nodes[i];
        gd_node *node = &ch->tree->nodes[v];
        double na = node->age * s;
        double lo = ch->cal->node_lo[v], hi = ch->cal->node_hi[v];
        for (size_t c = 0; c < node->nchildren; c++) {
            size_t ch_i = node->children[c];
            bool in = false;
            for (size_t j = a0; j < a1; j++) if (ch->sub_nodes[j] == ch_i) { in = true; break; }
            double cage = in ? ch->tree->nodes[ch_i].age * s : ch->tree->nodes[ch_i].age;
            if (cage > lo) lo = cage;
        }
        if (v == ch->sub_root[k]) {
            double pa = ch->tree->nodes[node->parent].age;
            if (pa < hi) hi = pa;
        }
        if (!(na > lo && na < hi)) feasible = false;
    }
    int accepted = 0;
    if (feasible) {
        double dprior = 0.0;
        for (size_t i = a0; i < a1; i++) {
            size_t v = ch->sub_nodes[i];
            double old = ch->sub_scratch[i - a0];
            dprior += node_logprior(ch->cal, v, old * s) - node_logprior(ch->cal, v, old);
            ch->tree->nodes[v].age = old * s;
        }
        double new_tp = tree_lp(ch);
        dprior += new_tp - ch->cur_tp;
        double new_ll;
        if (chain_eval_loglik(ch, NULL, &new_ll, err) != 0) {
            for (size_t i = a0; i < a1; i++) ch->tree->nodes[ch->sub_nodes[i]].age = ch->sub_scratch[i - a0];
            return GD_E_INTERNAL;
        }
        double log_alpha = (new_ll - ch->cur_loglik) + dprior + (double)(a1 - a0) * delta;
        if (accept_move(ch, log_alpha)) {
            for (size_t p = 0; p < ch->nparam; p++)
                if (ch->params[p].type == GD_MCMC_P_AGE)
                    ch->params[p].value = ch->tree->nodes[ch->params[p].node].age;
            ch->cur_loglik = new_ll;
            ch->cur_logprior += dprior;
            ch->cur_tp = new_tp;
            accepted = 1;
        } else {
            for (size_t i = a0; i < a1; i++) ch->tree->nodes[ch->sub_nodes[i]].age = ch->sub_scratch[i - a0];
        }
    }
    ch->sub_prop += 1.0;
    ch->sub_acc += (double)accepted;
    if (adapting) {
        ch->sub_adapt_count += 1.0;
        double gamma = 1.0 / pow(ch->sub_adapt_count + GD_MCMC_ADAPT_T0, 0.6);
        double logs = log(ch->sub_step) + gamma * ((double)accepted - GD_MCMC_TARGET_ACCEPT);
        if (logs < -30.0) logs = -30.0;
        if (logs > 30.0) logs = 30.0;
        ch->sub_step = exp(logs);
    }
    return 0;
}

static int chain_sweep(gd_mcmc_chain *ch, int adapting, gd_error *err)
{
    for (size_t p = 0; p < ch->nparam; p++) {
        gd_param *pp = &ch->params[p];
        double old_v = pp->value;
        double prop = old_v + gd_rng_normal(ch->rng, 0.0, pp->step);
        int accepted = 0;

        if (pp->type == GD_MCMC_P_AGE) {
            size_t v = pp->node;
            gd_node *node = &ch->tree->nodes[v];
            double lo = ch->cal->node_lo[v], hi = ch->cal->node_hi[v];
            for (size_t k = 0; k < node->nchildren; k++) {
                double ca = ch->tree->nodes[node->children[k]].age;
                if (ca > lo) lo = ca;
            }
            bool is_root = (v == ch->tree->root);
            if (!is_root) {
                double pa = ch->tree->nodes[node->parent].age;
                if (pa < hi) hi = pa;
            }
            /* The window [lo, hi] does not depend on this node's own age, so
             * both proposals below are valid Metropolis-Hastings moves:
             *   non-root: a' = a + N(0, step * (hi - lo)), symmetric;
             *   root:     a' - lo = (a - lo) exp(N(0, step)), with Hastings
             *             factor (a' - lo) / (a - lo).
             * The steps are dimensionless, so they adapt the same way on a
             * 500-year family and a 9000-year one, and nothing depends on a
             * numerical upper bound. */
            double old_age = node->age;
            double log_hastings = 0.0;
            if (is_root) {
                double h = old_age - lo;
                double hn = h * exp(gd_rng_normal(ch->rng, 0.0, pp->step));
                prop = lo + hn;
                log_hastings = log(hn) - log(h);
            } else {
                prop = old_age + gd_rng_normal(ch->rng, 0.0, pp->step * (hi - lo));
            }
            if (prop > lo && prop < hi) {
                double old_pr = node_logprior(ch->cal, v, old_age);
                node->age = prop;
                double new_ll;
                if (chain_eval_loglik(ch, NULL, &new_ll, err) != 0) { node->age = old_age; return GD_E_INTERNAL; }
                double new_tp = tree_lp(ch);
                double dprior = node_logprior(ch->cal, v, prop) - old_pr
                              + (new_tp - ch->cur_tp);
                if (accept_move(ch, (new_ll - ch->cur_loglik) + dprior + log_hastings)) {
                    pp->value = prop; ch->cur_loglik = new_ll;
                    ch->cur_logprior += dprior; ch->cur_tp = new_tp; accepted = 1;
                } else {
                    node->age = old_age;
                }
            }
        } else if (pp->type == GD_MCMC_P_LOG_SIGMA_B) {
            /* Hyperparameter: affects only the prior (branch-mult terms). */
            double sb_old = ch->sigma_b, sb_new = exp(prop);
            double own_new = halfnormal_logpdf(sb_new, ch->cfg.sigma_b_prior_sd) + prop;
            double own_old = halfnormal_logpdf(sb_old, ch->cfg.sigma_b_prior_sd) + old_v;
            double dprior = (own_new - own_old)
                          + (branchmult_prior_sum(ch, sb_new)
                             - branchmult_prior_sum(ch, sb_old));
            if (accept_move(ch, dprior)) {
                pp->value = prop; ch->sigma_b = sb_new;
                ch->cur_logprior += dprior; accepted = 1;
            }
        } else {
            /* Likelihood-affecting scalar/vector blocks. Tentatively apply the
             * proposal to the assembled buffer, evaluate, restore on reject. */
            double new_ll;
            double dprior = 0.0;
            double new_tp = ch->cur_tp;
            switch (pp->type) {
                case GD_MCMC_P_LOGRATE:
                    ch->logr = prop;
                    /* The lexical root rho * a_root moves with the rate. */
                    new_tp = tree_lp(ch);
                    dprior = norm_logpdf(prop, ch->cfg.rate_logmean, ch->cfg.rate_logsd)
                           - norm_logpdf(old_v, ch->cfg.rate_logmean, ch->cfg.rate_logsd)
                           + (new_tp - ch->cur_tp);
                    break;
                case GD_MCMC_P_LOGIT_EPS:
                    ch->epsilon = sigmoid(prop);
                    dprior = (beta_logpdf(sigmoid(prop), ch->cfg.eps_beta_a, ch->cfg.eps_beta_b)
                              + log(sigmoid(prop)) + log1p(-sigmoid(prop)))
                           - (beta_logpdf(sigmoid(old_v), ch->cfg.eps_beta_a, ch->cfg.eps_beta_b)
                              + log(sigmoid(old_v)) + log1p(-sigmoid(old_v)));
                    break;
                case GD_MCMC_P_BRANCH_LOGMULT: {
                    ch->branch_logmult[pp->node] = prop;
                    double mean = -0.5 * ch->sigma_b * ch->sigma_b;
                    dprior = norm_logpdf(prop, mean, ch->sigma_b)
                           - norm_logpdf(old_v, mean, ch->sigma_b);
                    break;
                }
                case GD_MCMC_P_CONCEPT_DEV:
                    ch->concept_dev[pp->concept] = prop;
                    dprior = norm_logpdf(prop, 0.0, ch->cfg.concept_rate_sd)
                           - norm_logpdf(old_v, 0.0, ch->cfg.concept_rate_sd);
                    break;
                case GD_MCMC_P_LOGIT_TIPERR:
                    ch->tip_hazard[pp->node] = -log1p(-sigmoid(prop));
                    dprior = tiperr_logprior(ch, pp->node, prop) - tiperr_logprior(ch, pp->node, old_v);
                    break;
                case GD_MCMC_P_LOG_SIGMA_C:
                    ch->sigma_c = exp(prop);
                    dprior = sigc_logprior(ch, prop) - sigc_logprior(ch, old_v);
                    break;
                default: break;
            }
            if (chain_eval_loglik(ch, NULL, &new_ll, err) != 0) return GD_E_INTERNAL;
            if (accept_move(ch, (new_ll - ch->cur_loglik) + dprior)) {
                pp->value = prop; ch->cur_loglik = new_ll;
                ch->cur_logprior += dprior; ch->cur_tp = new_tp; accepted = 1;
            } else {
                /* Restore the assembled buffer. */
                switch (pp->type) {
                    case GD_MCMC_P_LOGRATE:        ch->logr = old_v; break;
                    case GD_MCMC_P_LOGIT_EPS:      ch->epsilon = sigmoid(old_v); break;
                    case GD_MCMC_P_BRANCH_LOGMULT: ch->branch_logmult[pp->node] = old_v; break;
                    case GD_MCMC_P_CONCEPT_DEV:    ch->concept_dev[pp->concept] = old_v; break;
                    case GD_MCMC_P_LOG_SIGMA_C:    ch->sigma_c = exp(old_v); break;
                    case GD_MCMC_P_LOGIT_TIPERR:
                        ch->tip_hazard[pp->node] = -log1p(-sigmoid(old_v)); break;
                    default: break;
                }
            }
        }

        pp->prop += 1.0;
        pp->acc += (double)accepted;
        if (adapting) adapt_step(pp, accepted);
    }

    /* --- block age-scaling move --- */
    if (ch->have_block) {
        int rc = scale_move(ch, adapting, 0, err);
        if (rc != 0) return rc;
    }
    /* --- clock-preserving move --- */
    if (ch->have_clock) {
        int rc = scale_move(ch, adapting, 1, err);
        if (rc != 0) return rc;
    }
    /* --- non-centred relaxed-clock spread --- */
    if (ch->funnel_param < ch->nparam) {
        int rc = funnel_move(ch, adapting, err);
        if (rc != 0) return rc;
    }
    /* --- relaxed-clock rate shift --- */
    if (ch->shift_param < ch->nparam) {
        int rc = shift_move(ch, adapting, err);
        if (rc != 0) return rc;
    }
    /* --- subtree-scaling moves --- */
    for (size_t k = 0; k < ch->n_sub; k++) {
        int rc = subtree_move(ch, k, adapting, err);
        if (rc != 0) return rc;
    }
    return 0;
}

static void store_draw(gd_mcmc_chain *ch, gd_error *err)
{
    if (ch->n_stored >= ch->cap) return;
    double *row = ch->samples + ch->n_stored * ch->nparam;
    for (size_t p = 0; p < ch->nparam; p++) row[p] = ch->params[p].value;
    ch->sample_ll[ch->n_stored] = ch->cur_loglik;
    if (ch->pointwise != NULL) {
        double tot;
        if (chain_eval_loglik(ch, ch->pw_scratch, &tot, err) == 0)
            memcpy(ch->pointwise + ch->n_stored * ch->nconcepts,
                   ch->pw_scratch, ch->nconcepts * sizeof(double));
    }
    ch->n_stored++;
}

int gd_mcmc_chain_step(gd_mcmc_chain *ch, size_t n_iterations, gd_error *err)
{
    if (ch == NULL) { gd_error_set(err, GD_E_INVALID_ARG, "mcmc: NULL chain"); return GD_E_INVALID_ARG; }

    /* Progress bookkeeping only; with no callback none of this is reached and
     * the loop below is the one M8 shipped. */
    size_t total = ch->cfg.n_warmup + ch->cfg.n_samples * ch->cfg.thin;
    size_t every = ch->cfg.progress_every;
    if (every == 0) every = (total >= 100) ? total / 100 : 1;

    for (size_t it = 0; it < n_iterations; it++) {
        uint64_t g = ch->iter;
        int adapting = (ch->cfg.adapt && g < ch->cfg.n_warmup) ? 1 : 0;
        int rc = chain_sweep(ch, adapting, err);
        if (rc != 0) return rc;
        if (g >= ch->cfg.n_warmup) {
            uint64_t post = g - ch->cfg.n_warmup;
            if (post % (uint64_t)ch->cfg.thin == 0) store_draw(ch, err);
        }
        ch->iter++;

        if (ch->cfg.on_progress != NULL &&
            (ch->iter % (uint64_t)every == 0 || it + 1 == n_iterations)) {
            if (ch->cfg.on_progress(ch->cfg.progress_user, (size_t)ch->stream,
                                    (size_t)ch->iter, total) != 0) {
                gd_error_set(err, GD_E_CANCELLED, "mcmc: cancelled by caller");
                return GD_E_CANCELLED;
            }
        }
    }
    return 0;
}

int gd_mcmc_chain_run(gd_mcmc_chain *ch, gd_error *err)
{
    if (ch == NULL) { gd_error_set(err, GD_E_INVALID_ARG, "mcmc: NULL chain"); return GD_E_INVALID_ARG; }
    size_t total = ch->cfg.n_warmup + ch->cfg.n_samples * ch->cfg.thin;
    return gd_mcmc_chain_step(ch, total, err);
}

/* --------------------------------------------------------------- accessors */

size_t gd_mcmc_chain_param_count(const gd_mcmc_chain *ch) { return ch ? ch->nparam : 0; }
size_t gd_mcmc_chain_nsamples(const gd_mcmc_chain *ch) { return ch ? ch->n_stored : 0; }

double gd_mcmc_chain_sample(const gd_mcmc_chain *ch, size_t draw, size_t p)
{
    if (ch == NULL || draw >= ch->n_stored || p >= ch->nparam) return NAN;
    return ch->samples[draw * ch->nparam + p];
}

double gd_mcmc_chain_sample_loglik(const gd_mcmc_chain *ch, size_t draw)
{
    if (ch == NULL || draw >= ch->n_stored) return NAN;
    return ch->sample_ll[draw];
}

const double *gd_mcmc_chain_sample_pointwise(const gd_mcmc_chain *ch, size_t draw)
{
    if (ch == NULL || ch->pointwise == NULL || draw >= ch->n_stored) return NULL;
    return ch->pointwise + draw * ch->nconcepts;
}

double gd_mcmc_chain_acceptance(const gd_mcmc_chain *ch, size_t p)
{
    if (ch == NULL || p >= ch->nparam || ch->params[p].prop == 0.0) return NAN;
    return ch->params[p].acc / ch->params[p].prop;
}

double gd_mcmc_chain_step_size(const gd_mcmc_chain *ch, size_t p)
{
    if (ch == NULL || p >= ch->nparam) return NAN;
    return ch->params[p].step;
}

size_t gd_mcmc_chain_param_node(const gd_mcmc_chain *ch, size_t p)
{
    if (ch == NULL || p >= ch->nparam) return GD_NODE_NONE;
    return ch->params[p].node;
}

gd_mcmc_ptype gd_mcmc_chain_param_type(const gd_mcmc_chain *ch, size_t p)
{
    if (ch == NULL || p >= ch->nparam) return GD_MCMC_P_LOGRATE;
    return ch->params[p].type;
}

/* ------------------------------------------------------- checkpoint/resume */

static uint64_t fnv1a(const unsigned char *p, size_t n)
{
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) { h ^= (uint64_t)p[i]; h *= 1099511628211ULL; }
    return h;
}

/* Draws stored so far: parameter rows, their log-likelihoods, and the
 * pointwise rows when the chain keeps them. Part of a checkpoint so a resume
 * during sampling loses nothing. */
static size_t ckpt_draws_size(const gd_mcmc_chain *ch, size_t n_stored)
{
    size_t per = ch->nparam + 1 + (ch->pointwise != NULL ? ch->nconcepts : 0);
    return n_stored * per * 8;
}

static size_t ckpt_body_size(const gd_mcmc_chain *ch, size_t n_stored)
{
    return 4 + 4 + 4 + 8 + 16                 /* magic, ver, nparam, iter, rng */
         + 8 + 8                              /* binding, n_stored */
         + ckpt_draws_size(ch, n_stored)
         + ch->nparam * (8 + 8 + 8 + 8)       /* value, step, adapt, acc */
         + 8 + 8                              /* block_step, block_adapt_count */
         + 8 + 8                              /* clock_step, clock_adapt_count */
         + 8 + 8                              /* sub_step, sub_adapt_count */
         + 8 + 8                              /* funnel_step, funnel_adapt_count */
         + 8 + 8;                             /* shift_step, shift_adapt_count */
}

int gd_mcmc_chain_serialize(const gd_mcmc_chain *ch, unsigned char **buf,
                            size_t *len, gd_error *err)
{
    if (ch == NULL || buf == NULL || len == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "mcmc_serialize: NULL"); return GD_E_INVALID_ARG;
    }
    size_t body = ckpt_body_size(ch, ch->n_stored), total = body + 8;
    unsigned char *out = *buf;
    if (out == NULL) {
        out = (unsigned char *)gd_malloc_checked(total, err);
        if (out == NULL) { return GD_E_NOMEM; }
    } else if (*len < total) {
        gd_error_set(err, GD_E_INVALID_ARG, "mcmc_serialize: buffer too small"); return GD_E_INVALID_ARG;
    }
    unsigned char *q = out;
    uint32_t u32; uint64_t u64; double d;
    u32 = GD_MCMC_CKPT_MAGIC;   memcpy(q, &u32, 4); q += 4;
    u32 = GD_MCMC_CKPT_VERSION; memcpy(q, &u32, 4); q += 4;
    u32 = (uint32_t)ch->nparam; memcpy(q, &u32, 4); q += 4;
    u64 = ch->iter;             memcpy(q, &u64, 8); q += 8;
    u64 = ch->binding;          memcpy(q, &u64, 8); q += 8;
    u64 = ch->n_stored;         memcpy(q, &u64, 8); q += 8;
    if (gd_rng_serialize(ch->rng, q, GD_RNG_SERIAL_SIZE) != 0) {
        if (*buf == NULL) free(out);
        gd_error_set(err, GD_E_CHECKPOINT, "mcmc_serialize: rng"); return GD_E_CHECKPOINT;
    }
    q += 16;
    for (size_t p = 0; p < ch->nparam; p++) {
        d = ch->params[p].value;       memcpy(q, &d, 8); q += 8;
        d = ch->params[p].step;        memcpy(q, &d, 8); q += 8;
        d = ch->params[p].adapt_count; memcpy(q, &d, 8); q += 8;
        d = ch->params[p].acc;         memcpy(q, &d, 8); q += 8;
        /* `prop` (a running acceptance-stat denominator) is deliberately not
         * serialized: it never feeds the sampling trajectory, only diagnostics,
         * so a resumed run reproduces the draws regardless. */
    }
    d = ch->block_step;        memcpy(q, &d, 8); q += 8;
    d = ch->block_adapt_count; memcpy(q, &d, 8); q += 8;
    d = ch->clock_step;        memcpy(q, &d, 8); q += 8;
    d = ch->clock_adapt_count; memcpy(q, &d, 8); q += 8;
    d = ch->sub_step;          memcpy(q, &d, 8); q += 8;
    d = ch->sub_adapt_count;   memcpy(q, &d, 8); q += 8;
    d = ch->funnel_step;       memcpy(q, &d, 8); q += 8;
    d = ch->funnel_adapt_count; memcpy(q, &d, 8); q += 8;
    d = ch->shift_step;        memcpy(q, &d, 8); q += 8;
    d = ch->shift_adapt_count; memcpy(q, &d, 8); q += 8;
    for (size_t k = 0; k < ch->n_stored; k++) {
        memcpy(q, ch->samples + k * ch->nparam, ch->nparam * 8); q += ch->nparam * 8;
        memcpy(q, &ch->sample_ll[k], 8); q += 8;
        if (ch->pointwise != NULL) {
            memcpy(q, ch->pointwise + k * ch->nconcepts, ch->nconcepts * 8);
            q += ch->nconcepts * 8;
        }
    }
    u64 = fnv1a(out, body); memcpy(q, &u64, 8);
    *buf = out; *len = total;
    return 0;
}

int gd_mcmc_chain_restore(gd_mcmc_chain *ch, const unsigned char *buf,
                          size_t len, gd_error *err)
{
    if (ch == NULL || buf == NULL) { gd_error_set(err, GD_E_INVALID_ARG, "mcmc_restore: NULL"); return GD_E_INVALID_ARG; }
    const size_t head = 4 + 4 + 4 + 8 + 8 + 8;
    if (len < head) { gd_error_set(err, GD_E_CHECKPOINT, "mcmc_restore: short"); return GD_E_CHECKPOINT; }
    const unsigned char *q = buf;
    uint32_t magic, ver, nparam; uint64_t iter, binding, n_stored; double d;
    memcpy(&magic, q, 4); q += 4;
    memcpy(&ver, q, 4); q += 4;
    memcpy(&nparam, q, 4); q += 4;
    if (magic != GD_MCMC_CKPT_MAGIC || ver != GD_MCMC_CKPT_VERSION ||
        nparam != (uint32_t)ch->nparam) {
        gd_error_set(err, GD_E_CHECKPOINT, "mcmc_restore: version/shape mismatch");
        return GD_E_CHECKPOINT;
    }
    memcpy(&iter, q, 8); q += 8;
    memcpy(&binding, q, 8); q += 8;
    memcpy(&n_stored, q, 8); q += 8;
    if (binding != ch->binding) {
        gd_error_set(err, GD_E_CHECKPOINT,
                     "mcmc_restore: checkpoint was made with different data, tree, "
                     "calibrations, settings or chain");
        return GD_E_CHECKPOINT;
    }
    if (n_stored > ch->cap) {
        gd_error_set(err, GD_E_CHECKPOINT, "mcmc_restore: more draws than the chain holds");
        return GD_E_CHECKPOINT;
    }
    size_t body = ckpt_body_size(ch, (size_t)n_stored);
    if (len < body + 8) { gd_error_set(err, GD_E_CHECKPOINT, "mcmc_restore: short"); return GD_E_CHECKPOINT; }
    uint64_t sum_expect; memcpy(&sum_expect, buf + body, 8);
    if (fnv1a(buf, body) != sum_expect) {
        gd_error_set(err, GD_E_CHECKPOINT, "mcmc_restore: checksum mismatch"); return GD_E_CHECKPOINT;
    }
    ch->iter = iter;
    if (gd_rng_restore(ch->rng, q, GD_RNG_SERIAL_SIZE) != 0) {
        gd_error_set(err, GD_E_CHECKPOINT, "mcmc_restore: rng"); return GD_E_CHECKPOINT;
    }
    q += 16;
    for (size_t p = 0; p < ch->nparam; p++) {
        memcpy(&d, q, 8); q += 8; ch->params[p].value = d;
        memcpy(&d, q, 8); q += 8; ch->params[p].step = d;
        memcpy(&d, q, 8); q += 8; ch->params[p].adapt_count = d;
        memcpy(&d, q, 8); q += 8; ch->params[p].acc = d;
    }
    memcpy(&d, q, 8); q += 8; ch->block_step = d;
    memcpy(&d, q, 8); q += 8; ch->block_adapt_count = d;
    memcpy(&d, q, 8); q += 8; ch->clock_step = d;
    memcpy(&d, q, 8); q += 8; ch->clock_adapt_count = d;
    memcpy(&d, q, 8); q += 8; ch->sub_step = d;
    memcpy(&d, q, 8); q += 8; ch->sub_adapt_count = d;
    memcpy(&d, q, 8); q += 8; ch->funnel_step = d;
    memcpy(&d, q, 8); q += 8; ch->funnel_adapt_count = d;
    memcpy(&d, q, 8); q += 8; ch->shift_step = d;
    memcpy(&d, q, 8); q += 8; ch->shift_adapt_count = d;
    for (size_t k = 0; k < n_stored; k++) {
        memcpy(ch->samples + k * ch->nparam, q, ch->nparam * 8); q += ch->nparam * 8;
        memcpy(&ch->sample_ll[k], q, 8); q += 8;
        if (ch->pointwise != NULL) {
            memcpy(ch->pointwise + k * ch->nconcepts, q, ch->nconcepts * 8);
            q += ch->nconcepts * 8;
        }
    }
    ch->n_stored = (size_t)n_stored;
    /* Programmer invariant, and the reason the last q += 8 is not dead: the
     * cursor must land exactly on the end of the body that ckpt_body_size
     * predicted. A serializer and a deserializer that disagree about the
     * layout would otherwise restore a chain full of shifted values and carry
     * on sampling from it. */
    assert((size_t)(q - buf) == body);
    rebuild_assembled(ch);
    if (chain_eval_loglik(ch, NULL, &ch->cur_loglik, err) != 0) return GD_E_INTERNAL;
    ch->cur_tp = tree_lp(ch);
    ch->cur_logprior = chain_total_logprior(ch);
    return 0;
}

/* ----------------------------------------------------- multi-chain driver */

void gd_mcmc_result_destroy(gd_mcmc_result *res)
{
    if (res == NULL) return;
    free(res->samples); free(res->pointwise); free(res->loglik);
    free(res->rhat); free(res->ess); free(res->ess_tail); free(res->mcse);
    free(res->acc_rate);
    free(res->param_type); free(res->param_node);
    free(res);
}

/* Pearson correlation of two length-n vectors; NAN if either has no variance. */
static double pearson(const double *x, const double *y, size_t n)
{
    double mx = 0.0, my = 0.0;
    for (size_t i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
    mx /= (double)n; my /= (double)n;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dx = x[i] - mx, dy = y[i] - my;
        sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
    }
    if (!(sxx > 0.0) || !(syy > 0.0)) return NAN;
    return sxy / sqrt(sxx * syy);
}

/* Confounding diagnostic: for each sampled branch multiplier, the correlation
 * between its log-value and its log-duration across all draws. The most
 * negative such correlation is reported; a strong negative one means the branch
 * rate and the node ages are trading off (docs/model-spec.md §7.2). */
size_t gd_mcmc_result_find(const gd_mcmc_result *r, gd_mcmc_ptype type,
                           size_t node)
{
    if (r == NULL) return GD_MCMC_NO_PARAM;
    for (size_t p = 0; p < r->n_param; p++) {
        if (r->param_type[p] != type) continue;
        if (type == GD_MCMC_P_AGE || type == GD_MCMC_P_BRANCH_LOGMULT ||
            type == GD_MCMC_P_CONCEPT_DEV) {
            if (r->param_node[p] != node) continue;
        }
        return p;
    }
    return GD_MCMC_NO_PARAM;
}

double gd_mcmc_result_at(const gd_mcmc_result *r, size_t p, size_t i)
{
    if (r == NULL || p >= r->n_param) return NAN;
    if (i >= r->n_chains * r->n_draws) return NAN;
    return r->samples[i * r->n_param + p];
}

int gd_mcmc_result_column(const gd_mcmc_result *r, size_t p,
                          double *out, size_t n)
{
    if (r == NULL || out == NULL || p >= r->n_param) return GD_E_INVALID_ARG;
    size_t total = r->n_chains * r->n_draws;
    if (n != total) return GD_E_INVALID_ARG;
    size_t np = r->n_param;
    for (size_t i = 0; i < total; i++) out[i] = r->samples[i * np + p];
    return GD_E_OK;
}

static void compute_confounding(gd_mcmc_result *res, gd_tree *tree,
                                gd_calibration_set *cal)
{
    res->min_branch_duration_corr = NAN;
    res->branch_rate_confounded = false;
    size_t np = res->n_param, total = res->n_chains * res->n_draws;
    if (total < 4) return;

    /* Map each node to its age column (SIZE_MAX if the age is fixed). */
    /* NULL error: this routine reports a diagnostic on `res`, not a status, and
     * a run whose confounding check could not allocate still has its chronology.
     * The checked helper is here for the overflow guard. */
    size_t *age_col = (size_t *)gd_malloc_array_checked(tree->nnodes, sizeof(size_t), NULL);
    if (age_col == NULL) return;
    for (size_t i = 0; i < tree->nnodes; i++)
        age_col[i] = gd_mcmc_result_find(res, GD_MCMC_P_AGE, i);

    double *logm = (double *)gd_malloc_array_checked(total, sizeof(double), NULL);
    double *logt = (double *)gd_malloc_array_checked(total, sizeof(double), NULL);
    if (logm == NULL || logt == NULL) { free(age_col); free(logm); free(logt); return; }

    double min_corr = 2.0;
    bool any = false;
    for (size_t p = 0; p < np; p++) {
        if (res->param_type[p] != GD_MCMC_P_BRANCH_LOGMULT) continue;
        size_t v = res->param_node[p];
        if (v >= tree->nnodes) continue;
        size_t par = tree->nodes[v].parent;
        if (par == GD_NODE_NONE) continue;   /* no branch above the root */
        size_t k = 0;
        bool ok = true;
        for (size_t i = 0; i < total; i++) {
            double av = (age_col[v] != GD_MCMC_NO_PARAM)
                        ? gd_mcmc_result_at(res, age_col[v], i) : cal->node_lo[v];
            double ap = (par == GD_NODE_NONE) ? cal->node_lo[v]  /* unused */
                      : (age_col[par] != GD_MCMC_NO_PARAM)
                        ? gd_mcmc_result_at(res, age_col[par], i) : cal->node_lo[par];
            double t = ap - av;
            if (!(t > 0.0)) { ok = false; break; }
            logm[k] = gd_mcmc_result_at(res, p, i);
            logt[k] = log(t);
            k++;
        }
        if (!ok) continue;
        double c = pearson(logm, logt, k);
        if (!isnan(c)) { any = true; if (c < min_corr) min_corr = c; }
    }
    free(age_col); free(logm); free(logt);
    if (any) {
        res->min_branch_duration_corr = min_corr;
        res->branch_rate_confounded = (min_corr < -GD_MCMC_CONFOUND_THRESH);
    }
}

/* One chain of a parallel run: its own tree copy, model clone and chain. */
typedef struct {
    const gd_mcmc_model *model;
    gd_tree *src;
    gd_calibration_set *cal;
    const gd_mcmc_config *cfg;
    size_t c;
    gd_tree *tree;
    void *ctx;
    gd_mcmc_model local;
    gd_mcmc_chain *ch;
    gd_error err;
    int rc;
} chain_job;

static void chain_jobs_free(chain_job *jobs, size_t n, const gd_mcmc_model *model)
{
    if (jobs == NULL) return;
    for (size_t c = 0; c < n; c++) {
        gd_mcmc_chain_destroy(jobs[c].ch);
        if (jobs[c].ctx != NULL) model->clone_free(jobs[c].ctx);
        gd_tree_destroy(jobs[c].tree);
    }
    free(jobs);
}

#ifdef GD_THREADS
static void *chain_job_main(void *arg)
{
    chain_job *j = (chain_job *)arg;
    j->rc = GD_E_INTERNAL;
    j->tree = gd_tree_copy(j->src, &j->err);
    if (j->tree == NULL) return NULL;
    j->ctx = j->model->clone(j->model->ctx, j->tree, &j->err);
    if (j->ctx == NULL) return NULL;
    j->local = *j->model;
    j->local.ctx = j->ctx;
    j->ch = gd_mcmc_chain_create(&j->local, j->tree, j->cal, j->cfg, j->c, &j->err);
    if (j->ch == NULL) return NULL;
    j->rc = gd_mcmc_chain_run(j->ch, &j->err);
    return NULL;
}

/* Run the chains in batches of `nthreads`; report the first failure. */
static int run_chains_parallel(chain_job *jobs, const gd_mcmc_model *model, gd_tree *tree,
                               gd_calibration_set *cal, const gd_mcmc_config *cfg,
                               size_t n_chains, size_t nthreads, gd_error *err)
{
    pthread_t *th = (pthread_t *)gd_malloc_array_checked(nthreads, sizeof(pthread_t), err);
    if (th == NULL) return GD_E_NOMEM;
    for (size_t c = 0; c < n_chains; c++) {
        jobs[c].model = model; jobs[c].src = tree; jobs[c].cal = cal; jobs[c].cfg = cfg;
        jobs[c].c = c;
        gd_error e0 = GD_ERROR_INIT; jobs[c].err = e0;
    }
    int rc = 0;
    for (size_t b = 0; b < n_chains && rc == 0; b += nthreads) {
        size_t e = (b + nthreads < n_chains) ? b + nthreads : n_chains;
        size_t started = b;
        for (size_t c = b; c < e; c++) {
            if (pthread_create(&th[c - b], NULL, chain_job_main, &jobs[c]) != 0) {
                rc = gd_error_set(err, GD_E_INTERNAL, "mcmc_run: cannot start a chain thread");
                break;
            }
            started = c + 1;
        }
        for (size_t c = b; c < started; c++) pthread_join(th[c - b], NULL);
        for (size_t c = b; c < started && rc == 0; c++) {
            if (jobs[c].rc != 0) {
                *err = jobs[c].err;
                rc = jobs[c].rc;
            }
        }
    }
    free(th);
    return rc;
}
#else
static int run_chains_parallel(chain_job *jobs, const gd_mcmc_model *model, gd_tree *tree,
                               gd_calibration_set *cal, const gd_mcmc_config *cfg,
                               size_t n_chains, size_t nthreads, gd_error *err)
{
    (void)jobs; (void)model; (void)tree; (void)cal; (void)cfg; (void)n_chains; (void)nthreads;
    return gd_error_set(err, GD_E_INTERNAL, "mcmc_run: built without threads");
}
#endif

gd_mcmc_result *gd_mcmc_run(const gd_mcmc_model *model, gd_tree *tree,
                            gd_calibration_set *cal, const gd_mcmc_config *cfg,
                            size_t n_chains, gd_error *err)
{
    if (model == NULL || tree == NULL || cal == NULL || cfg == NULL || n_chains == 0) {
        gd_error_set(err, GD_E_INVALID_ARG, "mcmc_run: NULL/zero argument"); return NULL;
    }
    gd_mcmc_result *res = (gd_mcmc_result *)gd_calloc_checked(1, sizeof(*res), err);
    if (res == NULL) { gd_error_set(err, GD_E_NOMEM, "mcmc_run"); return NULL; }
    res->n_chains = n_chains;
    res->n_draws = cfg->n_samples;
    res->n_concepts = cfg->store_pointwise ? model->n_pointwise : 0;
    res->min_branch_duration_corr = NAN;

    gd_mcmc_chain *probe = gd_mcmc_chain_create(model, tree, cal, cfg, 0, err);
    if (probe == NULL) { free(res); return NULL; }
    res->n_param = probe->nparam;
    size_t np = res->n_param, nd = res->n_draws;
    res->param_type = (gd_mcmc_ptype *)gd_malloc_array_checked(np, sizeof(gd_mcmc_ptype), err);
    res->param_node = (size_t *)gd_malloc_array_checked(np, sizeof(size_t), err);
    if (res->param_type && res->param_node)
        for (size_t p = 0; p < np; p++) {
            res->param_type[p] = probe->params[p].type;
            res->param_node[p] = probe->params[p].node;
        }
    gd_mcmc_chain_destroy(probe);

    /* chains x draws is the base of both the sample matrix and the pointwise
     * matrix, and every factor is a run setting the caller chose. Multiplying
     * them out under a check keeps an absurd request a clean GD_E_NOMEM rather
     * than a wrapped size and a buffer the loops below would run past. */
    size_t ndraws_all, nsamples, npointwise = 0;
    if (!gd_size_mul(n_chains, nd, &ndraws_all) ||
        !gd_size_mul(ndraws_all, np, &nsamples) ||
        (res->n_concepts > 0 &&
         !gd_size_mul(ndraws_all, res->n_concepts, &npointwise))) {
        gd_error_set(err, GD_E_NOMEM, "mcmc_run: sample matrix size overflows");
        gd_mcmc_result_destroy(res); return NULL;
    }
    res->samples = (double *)gd_malloc_array_checked(nsamples, sizeof(double), err);
    res->loglik = (double *)gd_malloc_array_checked(ndraws_all ? ndraws_all : 1, sizeof(double), err);
    res->rhat = (double *)gd_malloc_array_checked(np, sizeof(double), err);
    res->ess = (double *)gd_malloc_array_checked(np, sizeof(double), err);
    res->ess_tail = (double *)gd_malloc_array_checked(np, sizeof(double), err);
    res->mcse = (double *)gd_malloc_array_checked(np, sizeof(double), err);
    res->acc_rate = (double *)gd_calloc_checked(np, sizeof(double), err);
    if (res->n_concepts > 0)
        res->pointwise = (double *)gd_malloc_array_checked(npointwise, sizeof(double), err);
    if (res->samples == NULL || res->loglik == NULL || res->rhat == NULL || res->ess == NULL ||
        res->ess_tail == NULL || res->mcse == NULL ||
        res->acc_rate == NULL || res->param_type == NULL || res->param_node == NULL ||
        (res->n_concepts > 0 && res->pointwise == NULL)) {
        gd_error_set(err, GD_E_NOMEM, "mcmc_run: arrays");
        gd_mcmc_result_destroy(res); return NULL;
    }

    size_t nthreads = cfg->threads == 0 ? n_chains : cfg->threads;
    bool parallel = false;
#ifdef GD_THREADS
    /* A progress callback is the caller's code, written for one chain at a
     * time (and its cancel stops "the" run): with one, stay sequential. */
    parallel = nthreads > 1 && n_chains > 1 && model->clone != NULL &&
               model->clone_free != NULL && cfg->on_progress == NULL;
#endif
    chain_job *jobs = NULL;
    if (parallel) {
        jobs = (chain_job *)gd_calloc_checked(n_chains, sizeof(chain_job), err);
        if (jobs == NULL) { gd_mcmc_result_destroy(res); return NULL; }
        int rc = run_chains_parallel(jobs, model, tree, cal, cfg, n_chains, nthreads, err);
        if (rc != 0) {
            chain_jobs_free(jobs, n_chains, model);
            gd_mcmc_result_destroy(res); return NULL;
        }
        /* Leave the caller's tree as the sequential path does: holding the
         * last chain's final ages. */
        for (size_t v = 0; v < tree->nnodes; v++)
            tree->nodes[v].age = jobs[n_chains - 1].tree->nodes[v].age;
    }

    for (size_t c = 0; c < n_chains; c++) {
        gd_mcmc_chain *ch;
        if (parallel) {
            ch = jobs[c].ch;
        } else {
            ch = gd_mcmc_chain_create(model, tree, cal, cfg, c, err);
            if (ch == NULL) { gd_mcmc_result_destroy(res); return NULL; }
            if (gd_mcmc_chain_run(ch, err) != 0) {
                gd_mcmc_chain_destroy(ch); gd_mcmc_result_destroy(res); return NULL;
            }
        }
        for (size_t draw = 0; draw < nd; draw++) {
            double *dst = res->samples + (c * nd + draw) * np;
            for (size_t p = 0; p < np; p++) dst[p] = gd_mcmc_chain_sample(ch, draw, p);
            res->loglik[c * nd + draw] = gd_mcmc_chain_sample_loglik(ch, draw);
            if (res->n_concepts > 0) {
                const double *pw = gd_mcmc_chain_sample_pointwise(ch, draw);
                if (pw != NULL)
                    memcpy(res->pointwise + (c * nd + draw) * res->n_concepts,
                           pw, res->n_concepts * sizeof(double));
            }
        }
        for (size_t p = 0; p < np; p++) {
            double a = gd_mcmc_chain_acceptance(ch, p);
            res->acc_rate[p] += (isnan(a) ? 0.0 : a) / (double)n_chains;
        }
        if (!parallel) gd_mcmc_chain_destroy(ch);
    }
    if (parallel) chain_jobs_free(jobs, n_chains, model);

    /* ndraws_all == n_chains * nd, already checked above. */
    double *col = (double *)gd_malloc_array_checked(ndraws_all, sizeof(double), err);
    if (col == NULL) { gd_error_set(err, GD_E_NOMEM, "mcmc_run: col"); gd_mcmc_result_destroy(res); return NULL; }
    for (size_t p = 0; p < np; p++) {
        (void)gd_mcmc_result_column(res, p, col, ndraws_all);
        res->rhat[p] = gd_rhat_rank(col, n_chains, nd);
        res->ess[p] = gd_ess_bulk(col, n_chains, nd);
        res->ess_tail[p] = gd_ess_tail(col, n_chains, nd);
        res->mcse[p] = gd_mcse_mean(col, n_chains, nd);
    }
    free(col);

    compute_confounding(res, tree, cal);
    return res;
}


/*
 * gd_mcmc.h — fixed-tree posterior sampler for Clepsydra (M8 + M9).
 *
 * Metropolis-within-Gibbs over a fixed rooted topology (docs/inference.md §2,
 * docs/model-spec.md §9). The sampler is a typed-parameter engine; which blocks
 * are active is chosen by the config:
 *
 *   - internal node ages a(v)                         (always, when free)
 *   - global strict-clock log-rate                    (cfg.sample_rate)
 *   - relaxed branch multipliers log m_b + sigma_b     (cfg.sample_branch, M9)
 *   - global contamination epsilon                    (cfg.sample_contamination, M9)
 *   - per-concept log-rate deviations mu_c            (cfg.sample_concept_rates, M9)
 *
 * The relaxed clock is `log m_b ~ N(-sigma_b^2/2, sigma_b^2)` (so E[m_b]=1) with
 * a HalfNormal shrinkage prior on sigma_b (docs/model-spec.md §4); the effective
 * per-branch rate for a concept is `exp(logr + mu_c) * m_b`. Branch multipliers
 * and node ages are confounded (docs/model-spec.md §7.2); gd_mcmc_run reports a
 * confounding diagnostic.
 *
 * Properties: single-site random-walk Metropolis (ages on the natural scale
 * with parent/child order by rejection, rates/scales on the log scale, epsilon
 * on the logit scale); adaptation ONLY during warm-up, frozen before retained
 * draws; per-chain PCG32 streams (bit-reproducible); versioned + checksummed
 * checkpoint/resume; optional per-concept pointwise loglik. The library never
 * calls exit()/abort().
 */
#ifndef GD_MCMC_H
#define GD_MCMC_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"
#include "gd_calibration.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* gd_mcmc_ptype is public (gd_types.h): it labels parameters in a reported
 * result, so a caller reading one needs the same vocabulary the sampler uses. */

#include "gd_mcmc_model.h"

typedef enum {
    GD_TREE_PRIOR_UNIFORM = 0,
    GD_TREE_PRIOR_NONE    = 1,
    GD_TREE_PRIOR_YULE    = 2    /* sensitivity alternative; tips at 0 only */
} gd_tree_prior_kind;


typedef struct gd_mcmc_config {
    uint64_t seed;
    size_t   n_warmup;
    size_t   n_samples;
    size_t   thin;              /* >= 1 */
    bool     adapt;
    bool     store_pointwise;

    /* Global strict-clock rate (log r ~ N(rate_logmean, rate_logsd^2)). */
    bool     sample_rate;
    double   fixed_rate;
    double   rate_logmean, rate_logsd, rate_step0;

    /*
     * Prior on node ages (docs/roadmap.md E1). GD_TREE_PRIOR_UNIFORM, the
     * default, is proper: a log-normal "lexical root" L = rho * a_root and
     * uniform ages below the root (gd_tree_prior.h). GD_TREE_PRIOR_NONE is
     * the pre-0.3 flat measure on the ordered region; it is refused unless
     * every free age is bounded by calibrations, since otherwise the
     * posterior is improper.
     */
    gd_tree_prior_kind tree_prior;
    double   root_lex_logmean;  /* default -1.2: 0.30 replacements per concept */
    double   root_lex_logsd;    /* default 1.5 */
    double   yule_kappa;        /* GD_TREE_PRIOR_YULE: lambda * root; 0 = H(n) - 1 */

    /* Likelihood switched off: the sampler then draws from the prior the
     * posterior is built on, exactly (E1.5). */
    bool     prior_only;

    /* Node-age proposals: a non-root age moves by N(0, step * window), where
     * the window is the interval its parent, children and calibrations leave
     * it; the root moves on the log scale of its height above its oldest
     * child. age_step0_frac is the initial step for both. */
    double   age_step0_frac;

    /* Clock-preserving move: scale every free age by exp(delta) and the rate
     * by exp(-delta), which leaves lexical branch lengths unchanged. It is
     * what moves a chain along the rate-age ridge. Needs sample_rate. */
    bool     clock_move;
    double   clock_step0;       /* default 0.10 (log scale) */

    /* Subtree-scaling move: for each free internal non-root node, scale the
     * free internal ages of its subtree by one factor. Nested nodes with
     * close ages, which single-site updates move only one at a time, then
     * move together. */
    bool     subtree_move;
    double   subtree_step0;     /* default 0.10 (log scale) */

    /* Block age-scaling move (docs/inference.md §2.2): once per sweep, scale
     * every free node age by a shared factor exp(N(0, block_step0)), which
     * decorrelates nested ages that single-site moves mix through slowly. Off
     * by default; the CLI enables it for deep trees. Requires >= 2 free ages. */
    bool     block_age_scaling;
    double   block_step0;   /* default 0.10 (log scale) */

    /* Relaxed branch clock (M9): log m_b ~ N(-sigma_b^2/2, sigma_b^2),
     * sigma_b ~ HalfNormal(0, sigma_b_prior_sd). */
    bool     sample_branch;
    double   sigma_b_prior_sd;  /* default 0.30 */
    double   sigma_b_init;      /* default 0.10 */
    double   sigma_b_step0;     /* default 0.30 (log scale) */
    double   branch_step0;      /* default 0.20 (log scale) */

    /* Contamination (M9): epsilon ~ Beta(eps_beta_a, eps_beta_b). */
    bool     sample_contamination;
    double   fixed_epsilon;     /* used when !sample_contamination (default 0) */
    double   eps_beta_a, eps_beta_b;  /* default 1, 19 */
    double   eps_init;          /* default 0.05 */
    double   eps_step0;         /* default 0.40 (logit scale) */

    /* Among-concept rate heterogeneity, integrated out (E2.1): each concept's
     * rate is rho * lambda_c, log lambda_c ~ N(-sigma_c^2/2, sigma_c^2), summed
     * over an 8-point quadrature in the likelihood; sigma_c ~
     * HalfNormal(rate_mix_prior_sd), sampled on the log scale. */
    bool     sample_rate_mix;
    int      rate_mix_family;   /* GD_T1_MIX_LOGNORMAL (0) or _GAMMA (1) */
    double   rate_mix_prior_sd; /* default 1.0 */
    double   rate_mix_init;     /* default 0.5 */
    double   rate_mix_step0;    /* default 0.30 (log scale) */

    /* Observation error on tips (E2.9): for a tip v with tip_err_mean[v] > 0,
     * e_v ~ Beta(m k, (1 - m) k) with m = tip_err_mean[v], k = tip_err_conc[v],
     * sampled on the logit scale. Arrays of nnodes, borrowed; NULL = none. */
    const double *tip_err_mean;
    const double *tip_err_conc;

    /* Per-concept rate deviations (M9): mu_c ~ N(0, concept_rate_sd^2). */
    bool     sample_concept_rates;
    double   concept_rate_sd;   /* default 0.50 */
    double   concept_step0;     /* default 0.30 (log scale) */

    /*
     * Optional progress reporting (M10). `on_progress` is called during
     * sampling with the chain's stream index, the iterations that chain has
     * completed, and the number it will run in total; returning nonzero
     * cancels, and the run reports GD_E_CANCELLED rather than handing back a
     * partial chain. `progress_every` sets the reporting interval in
     * iterations (0 = roughly every 1% of a chain, plus the last iteration of
     * each call).
     *
     * NULL — the default — means no reporting and not one line of the
     * sampler's arithmetic changes; the draws are bit-identical either way.
     */
    int    (*on_progress)(void *user, size_t chain, size_t iter, size_t total);
    void    *progress_user;
    size_t   progress_every;
} gd_mcmc_config;

void gd_mcmc_config_default(gd_mcmc_config *cfg);

/* ---- Single chain ---- */

typedef struct gd_mcmc_chain gd_mcmc_chain;

gd_mcmc_chain *gd_mcmc_chain_create(const gd_mcmc_model *model, gd_tree *tree,
                                    gd_calibration_set *cal,
                                    const gd_mcmc_config *cfg, uint64_t stream,
                                    gd_error *err);
void gd_mcmc_chain_destroy(gd_mcmc_chain *chain);

int gd_mcmc_chain_step(gd_mcmc_chain *chain, size_t n_iterations, gd_error *err);
int gd_mcmc_chain_run(gd_mcmc_chain *chain, gd_error *err);

size_t gd_mcmc_chain_param_count(const gd_mcmc_chain *chain);
size_t gd_mcmc_chain_nsamples(const gd_mcmc_chain *chain);
double gd_mcmc_chain_sample(const gd_mcmc_chain *chain, size_t draw, size_t p);
double gd_mcmc_chain_sample_loglik(const gd_mcmc_chain *chain, size_t draw);
const double *gd_mcmc_chain_sample_pointwise(const gd_mcmc_chain *chain, size_t draw);
double gd_mcmc_chain_acceptance(const gd_mcmc_chain *chain, size_t p);
double gd_mcmc_chain_step_size(const gd_mcmc_chain *chain, size_t p);

/* Parameter metadata. `param_node` is the node index for AGE / BRANCH_LOGMULT
 * (GD_NODE_NONE otherwise); `param_type` is the kind. */
size_t        gd_mcmc_chain_param_node(const gd_mcmc_chain *chain, size_t p);
gd_mcmc_ptype gd_mcmc_chain_param_type(const gd_mcmc_chain *chain, size_t p);

/* ---- Checkpoint / resume ---- */

int gd_mcmc_chain_serialize(const gd_mcmc_chain *chain, unsigned char **buf,
                            size_t *len, gd_error *err);
int gd_mcmc_chain_restore(gd_mcmc_chain *chain, const unsigned char *buf,
                          size_t len, gd_error *err);

/* ---- Multi-chain driver + diagnostics ---- */

typedef struct gd_mcmc_result {
    size_t  n_chains;
    size_t  n_draws;
    size_t  n_param;
    size_t  n_concepts;
    double *samples;     /* [chain][draw][param] */
    double *pointwise;   /* [chain][draw][concept] or NULL */
    double *loglik;      /* [chain][draw]: total log-likelihood of each draw */
    double *rhat;        /* n_param: max of bulk and folded (tail) R-hat */
    double *ess;         /* n_param: bulk ESS */
    double *ess_tail;    /* n_param: tail ESS (5% and 95% quantiles) */
    double *mcse;        /* n_param: Monte Carlo SE of the mean */
    double *acc_rate;    /* n_param */
    gd_mcmc_ptype *param_type; /* n_param */
    size_t *param_node;        /* n_param (GD_NODE_NONE where N/A) */

    /* Confounding diagnostic (docs/model-spec.md §7.2). When branch multipliers
     * are sampled, the most negative posterior correlation between a branch's
     * log-multiplier and its log-duration across draws; a strong negative value
     * means the branch rate and the node ages are trading off (confounded).
     * branch_rate_confounded is 1 when that correlation crosses the threshold. */
    bool   branch_rate_confounded;
    double min_branch_duration_corr;  /* NAN if branches not sampled */
} gd_mcmc_result;

/*
 * Reading a result back.
 *
 * `samples` is a flat array and its stride -- [(chain * n_draws + draw) *
 * n_param + p] -- used to be spelled out at seven call sites, two of them
 * inside this file, where the sampler decoded its own output. These three
 * functions are that arithmetic, written once.
 *
 * gd_mcmc_result_column is not a convenience over gd_mcmc_result_at: it is the
 * one to use for a whole parameter's draws, because it is a single call rather
 * than one per draw, and both of the bulk sites run np * n_chains * n_draws
 * times.
 */
#define GD_MCMC_NO_PARAM ((size_t)-1)

/* The parameter column of this kind at this node, or GD_MCMC_NO_PARAM.
 * `node` is ignored for kinds that do not have one. */
size_t gd_mcmc_result_find(const gd_mcmc_result *r, gd_mcmc_ptype type,
                           size_t node);

/* One stored draw. `i` is a flat draw index in [0, n_chains * n_draws);
 * chain c's draw d sits at c * n_draws + d. NAN when out of range. */
double gd_mcmc_result_at(const gd_mcmc_result *r, size_t p, size_t i);

/* Copy parameter p's draws into out[0 .. n), in flat draw order, so chains
 * stay adjacent and gd_rhat_rank can read them directly. `n` must be
 * n_chains * n_draws. Returns GD_E_OK, or GD_E_INVALID_ARG. */
int gd_mcmc_result_column(const gd_mcmc_result *r, size_t p,
                          double *out, size_t n);

gd_mcmc_result *gd_mcmc_run(const gd_mcmc_model *model, gd_tree *tree,
                            gd_calibration_set *cal, const gd_mcmc_config *cfg,
                            size_t n_chains, gd_error *err);
void gd_mcmc_result_destroy(gd_mcmc_result *res);



#endif /* GD_MCMC_H */

/*
 * gd_analysis.h — end-to-end analysis pipelines over in-memory inputs.
 *
 * This is libclepsydra's public API: the shared runner behind both the native
 * CLI (`cli/commands.c`) and the WebAssembly shim (`web/clepsydra_wasm.c`),
 * and the surface any other caller should build against. See PLAN.md M10.
 *
 * Everything here works from buffers, never from paths: the browser has no
 * filesystem, so the caller does its own I/O and its own formatting and this
 * module does the science. That is why the page and the CLI cannot compute
 * different numbers.
 *
 * Two pipelines, mirroring the two CLI commands:
 *   gd_analysis_tree_run — fixed-topology T1 chronology (MCMC over node ages)
 *   gd_analysis_pair_run — P1 pairwise divergence, or an alternative decay law
 *
 * Like the rest of the library, nothing here calls exit(), printf() or the
 * filesystem; every failure comes back through gd_error.
 */
#ifndef GD_ANALYSIS_H
#define GD_ANALYSIS_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/*
 * "This buffer is NUL-terminated; measure it." Use it instead of a length when
 * the buffer came from something that terminates, such as a C string literal
 * or JavaScript's stringToNewUTF8.
 *
 * This used to be spelled 0, which made a genuinely empty buffer impossible to
 * express and made a zeroed input struct silently mean "call strlen on it".
 * SIZE_MAX cannot collide with a real length, so 0 now means what it says.
 */
#define GD_LEN_NUL_TERMINATED ((size_t)-1)

/* A posterior summary of one scalar quantity. `hpd_lo`/`hpd_hi` bound the
 * narrowest interval holding `hpd_level` of the draws (config, default 0.90). */
typedef struct gd_analysis_summary {
    double mean, median, sd, hpd_lo, hpd_hi;
} gd_analysis_summary;

/* ------------------------------------------------- calibration variation */

/*
 * One calibration as read from the table, before it is resolved onto a tree.
 * Reading them out separately lets a caller enumerate the calibrations —
 * to drop each in turn, or to report which one it withheld — without parsing
 * the TSV a second time and risking a different reading of it.
 */
typedef struct gd_analysis_calibration {
    char  *id;
    char  *target;
    char       *type_name; /* as written in the table */
    gd_cal_type type;
    double p1, p2, p3;
    bool   is_hard;
    double outlier;        /* outlier_prob column; 0 when absent or '.' */
    /* The table's source_key column, or "." when it has none. Provenance for
     * the calibration, carried through so the chronology intake does not have
     * to decode this table a second time to reach it. */
    char  *source_key;
    /* The age this calibration pins for its target: a uniform window's centre,
     * otherwise p1. This is the value to score against when the calibration is
     * withheld and its target re-dated. */
    double point;
    /* 1 for a fixed tip at age 0 — a living language, whose "recovery" is not
     * an interesting test. */
    bool   is_extant_tip;
} gd_analysis_calibration;

/* `len` follows the same rule as the input buffers below: a byte count, or
 * GD_LEN_NUL_TERMINATED to have it measured. On success *out is an owned array
 * of *n calibrations, released with gd_analysis_calibrations_free. */
int  gd_analysis_calibrations_read(const char *calibrations_tsv, size_t len,
                                   gd_analysis_calibration **out, size_t *n,
                                   gd_error *err);
void gd_analysis_calibrations_free(gd_analysis_calibration *cals, size_t n);

/*
 * A systematic variation of the calibration set, for sensitivity analysis.
 *
 *   drop_index — omit this row (leave-one-calibration-out); GD_ANALYSIS_KEEP_ALL
 *                keeps every calibration.
 *   widen      — multiply the width of every *soft* window by this factor
 *                (a uniform window about its centre, a normal/lognormal sd).
 *                Fixed tips and hard bounds have no soft window and are left
 *                alone. 1.0 changes nothing.
 *
 * The two compose: a run can drop one calibration and loosen the rest.
 */
#define GD_ANALYSIS_KEEP_ALL ((size_t)-1)

typedef struct gd_analysis_cal_variant {
    size_t drop_index;
    double widen;
} gd_analysis_cal_variant;

/* ------------------------------------------------------------- date-tree */

/*
 * Input buffers.
 *
 * Each `*_len` is the length of its buffer in bytes, or GD_LEN_NUL_TERMINATED
 * to have the library measure it. A length of 0 means an empty buffer, which
 * is a real input and will be reported as such rather than measured.
 *
 * The `*_name` fields are labels for diagnostics ("forms.tsv: ...") and may be
 * NULL; they are never opened — this module reads only what it is handed.
 */
typedef struct gd_analysis_tree_input {
    const char *newick;           size_t newick_len;
    const char *forms_tsv;        size_t forms_len;
    const char *calibrations_tsv; size_t cal_len;
    const char *newick_name, *forms_name, *cal_name;
    /* Which coding to read when the forms hold several (analysis_id column);
     * NULL is an error in that case. */
    const char *analysis_id;
} gd_analysis_tree_input;

typedef struct gd_analysis_tree_config {
    size_t   n_warmup;          /* 2000 */
    size_t   n_samples;         /* 2000 */
    size_t   n_chains;          /* 4 */
    uint64_t seed;              /* 42 */
    double   rate_logmean;      /* -8.80 = Swadesh/year */
    double   rate_logsd;        /* 0.80 */
    double   epsilon;           /* 0.05: contamination weight when not sampled */
    /* Sample epsilon ~ Beta(1, 19) instead of fixing it (E2.2). Default 1. */
    bool     sample_epsilon;
    bool     relaxed;           /* relaxed branch clock (M9); 0 = strict */
    bool     block_age_scaling; /* 1: decorrelate nested ages (deep trees) */
    /* Prior on the lexical root L = rho * a_root, in expected replacements
     * per concept along a root-to-tip path (gd_tree_prior.h). */
    double   root_lex_logmean;  /* -1.2, i.e. 0.30 */
    double   root_lex_logsd;    /* 1.5 */
    bool     prior_only;        /* 0: sample with the likelihood switched off */
    /* Sensitivity alternative to the uniform age prior: Yule (pure birth,
     * complete sampling) conditioned on the root, lambda = kappa / root.
     * Tips must all be modern. kappa 0 means H(n) - 1 for n tips. */
    bool     yule;              /* 0 */
    /* Among-concept rate heterogeneity, integrated out (E2.1): log-normal
     * multipliers with mean 1 and a sampled spread sigma_c ~ HalfNormal(sd).
     * On by default; 0 is the homogeneous model. */
    bool     concept_rates;     /* 1 */
    /* Forms rows at or above this loan probability are read as missing
     * (E2.4); above 1 keeps them all. */
    double   loan_threshold;    /* 0.5 */
    double   concept_rate_prior_sd; /* 1.0 */
    double   yule_kappa;        /* 0 */
    bool     compare_prior;     /* 1: also run the prior alone and report it */
    double   hpd_level;         /* 0.90 */

    /*
     * Progress reporting, passed straight to the sampler: called with the
     * chain index, that chain's completed iterations and its total; a nonzero
     * return cancels the run, which then fails with GD_E_CANCELLED instead of
     * returning a half-finished chronology. `progress_every` is the interval
     * in iterations (0 = about every 1%). NULL leaves the draws bit-identical
     * to a run without reporting — the browser needs a progress bar, and a
     * progress bar must not change the answer.
     */
    int    (*on_progress)(void *user, size_t chain, size_t iter, size_t total);
    void    *progress_user;
    size_t   progress_every;

    /* Keep at most this many evenly spaced draws per chain per parameter, so a
     * caller can draw a trace. 0 (the default) keeps none: a chronology summary
     * does not need them, and a long run's full sample is large. */
    size_t   trace_draws;

    /* Sensitivity: which calibrations this run is given. Defaults to all of
     * them, unwidened. */
    gd_analysis_cal_variant calibrations;
} gd_analysis_tree_config;

void gd_analysis_tree_config_default(gd_analysis_tree_config *cfg);

/* One node's posterior age, in years BP. `label` is borrowed from the result's
 * tree and is NULL for an unlabelled node. */
typedef struct gd_analysis_node_age {
    size_t              node;
    const char         *label;
    gd_analysis_summary age;
    double              rhat, ess, ess_tail, mcse;

    /*
     * The three kinds of age (docs/roadmap.md, "The goal"):
     *   age      — calendar years BP (above);
     *   lexical  — rho * age: expected replacements per concept since the
     *              node, which the cognates measure directly;
     *   relative — age / root age, which needs no calendar scale at all.
     * Under a relaxed clock rho is the mean rate, so lexical and relative
     * ages partly reflect the prior on branch multipliers.
     */
    gd_analysis_summary lexical, relative;

    /* The same node under the prior alone (a likelihood-free run), and how
     * far the data moved it: `prior_overlap` is the overlap coefficient of the
     * prior and posterior densities of log age (1 = indistinguishable), and
     * `contraction` is 1 - var_post / var_prior of log age. NAN when no
     * prior run was made. */
    gd_analysis_summary prior;
    double              prior_overlap, contraction;
} gd_analysis_node_age;

/* One sampled parameter, in the sampler's own scale. For an age or a branch
 * multiplier, `node`/`label` name the node; otherwise `node` is GD_NODE_NONE
 * and `label` is NULL. `type_name` is a stable short string ("log_rate",
 * "age", "log_sigma_b", "logit_eps", "branch_logmult", "concept_dev"). */
typedef struct gd_analysis_param {
    gd_mcmc_ptype       type;
    const char         *type_name;
    size_t              node;
    const char         *label;
    gd_analysis_summary stat;
    /* rhat is max(bulk, folded) rank R-hat; ess the bulk ESS; ess_tail the
     * ESS of the 5% and 95% quantiles; mcse the Monte Carlo SE of the mean
     * (Vehtari et al. 2021, matched to ArviZ). */
    double              rhat, ess, ess_tail, mcse, acceptance;
} gd_analysis_param;

enum { GD_SENS_CLOCK = 0, GD_SENS_LEXROOT = 1, GD_SENS_CALIBRATIONS = 2, GD_SENS_LIKELIHOOD = 3 };
#define GD_SENS_THRESHOLD 0.2

#define GD_CONV_RHAT_MAX      1.01
#define GD_CONV_ESS_PER_CHAIN 100.0

/* One concept after the run (E2.2): averaged over thinned posterior draws.
 * `p_contaminated` is the posterior probability that the concept came from
 * the noise (contamination) component rather than the tree. `rate_mult` is
 * its posterior rate relative to the mean rate under the concept-rate
 * mixture (1 when the mixture is off). `n_observed` counts the tips with a
 * form. `name` is owned by the result. */
typedef struct gd_analysis_concept {
    char   *name;
    size_t  n_observed;
    bool    incompatible;   /* no T1 history on this topology fits; excluded */
    double  p_contaminated, rate_mult, mean_loglik;
} gd_analysis_concept;

/* One posterior predictive check (E2.10): a statistic of the observed data,
 * its distribution over data replicated from posterior draws, and the
 * posterior predictive p-value P(replicated >= observed). A p-value near 0
 * or 1 means the model does not reproduce that feature of the data. */
typedef struct gd_analysis_check {
    char   name[96];
    double observed, rep_mean, rep_lo90, rep_hi90, ppp;
} gd_analysis_check;

typedef struct gd_analysis_tree_result {
    gd_tree *tree;                    /* owned by the result */

    gd_analysis_node_age *ages;       /* one per sampled node age */
    size_t                nages;
    gd_analysis_param    *params;     /* every sampled parameter */
    size_t                nparams;

    size_t n_languages, n_nodes, n_concepts;
    size_t n_chains, n_warmup, n_samples, n_draws;
    size_t n_skipped_rows;            /* form rows for languages not in the tree */
    /* Form rows whose primary_weight was present but unreadable, and so was
     * taken as 1.0. A typo in one cell of a large forms file is otherwise
     * silent: the run succeeds and the wrong synonym may have been chosen. */
    size_t n_unreadable_weights;
    /* Form rows that carried no observation ('.', '?', empty, or status
     * missing/excluded/absent). Skipped, not read as a cognate class. */
    size_t n_missing_cognates;
    size_t n_loans;                   /* form rows dropped as loans (E2.4) */
    bool   relaxed;

    /* Convergence verdict over every sampled parameter (Vehtari et al. 2021):
     * converged when max R-hat < GD_CONV_RHAT_MAX and both the bulk and the
     * tail ESS are at least GD_CONV_ESS_PER_CHAIN per chain. A run that fails
     * is still returned -- its draws are evidence of the failure -- but the
     * CLI exits non-zero unless told otherwise. */
    /*
     * Where the calendar scale comes from. `scale_source` names the anchors
     * present: "calibration" (an internal node is calibrated), "tip" (a tip
     * has a non-zero date), both joined by "+", or "clock_prior" when there
     * is neither and calendar ages are the lexical ones divided by the clock
     * prior. Static storage.
     *
     * `rate_contraction` is 1 - var_post / var_prior of log rho. Near 0 means
     * the data and calibrations left the clock where its prior put it, so
     * every calendar age carries the prior's width. `scale_from_prior` flags
     * that (contraction < 0.5).
     */
    gd_analysis_concept *concepts;    /* n_concepts entries, or NULL */
    gd_analysis_check   *checks;      /* posterior predictive checks */
    size_t               n_checks;

    /* One row per calibration applied to a sampled node (E2.6): the target's
     * posterior median and HPD, and for an outlier-prone calibration its
     * prior and posterior outlier probability. Strings borrowed from the
     * tree or owned by the result (id). */
    struct gd_analysis_cal_check {
        char       *id;
        const char *target;
        double      outlier_prior, outlier_post;
        gd_analysis_summary target_age;
    } *cal_checks;
    size_t n_cal_checks;
    size_t n_incompatible;            /* concepts excluded as incompatible */

    const char *scale_source;
    double rate_contraction;
    bool   scale_from_prior;
    gd_analysis_summary lexical_root; /* rho * root age */

    /*
     * Where the calendar root's width comes from. log a_root = log L - log rho,
     * so Var(log a_root) = Var(log L) + Var(log rho) - 2 Cov. `var_share_data`
     * is Var(log L) / Var(log a_root), the part the cognates set;
     * `var_share_clock` is Var(log rho) / Var(log a_root); the remainder is
     * their covariance. NAN when the rate is not sampled.
     */
    double var_share_data, var_share_clock;

    /*
     * Power-scaling sensitivity (after Kallioinen et al. 2024, first order):
     * raising one component of the posterior to a power alpha moves the
     * posterior of a quantity g by d/d(log alpha) at alpha = 1, which is a
     * covariance with that component's log density. Reported as the larger of
     * the shifts in mean and sd, in posterior sd, for g = log root age
     * (`sens_root`) and g = log lexical root (`sens_lexroot`). Index with
     * GD_SENS_*. Above 0.2 a component is influential. A quantity that is
     * prior-sensitive but not likelihood-sensitive is prior-dominated; one
     * sensitive to both shows prior-data conflict.
     */
    double sens_root[4], sens_lexroot[4];

    bool   converged;
    size_t n_unconverged;             /* parameters failing any criterion */
    double max_rhat, min_ess_bulk, min_ess_tail;

    /* Branch-rate × node-age confounding (docs/model-spec.md §7.2).
     * `min_branch_duration_corr` is NAN when branch multipliers are not
     * sampled. */
    bool   branch_rate_confounded;
    double min_branch_duration_corr;

    /* Thinned draws, or NULL when cfg.trace_draws was 0. Indexed
     * [(param * n_chains + chain) * trace_draws + k]; the draws are evenly
     * spaced through each chain, so a plot of them shows the whole run rather
     * than its first moments. Chains are kept apart on purpose: two chains that
     * disagree are what R̂ is measuring, and a merged trace hides it. */
    double *trace;
    size_t  trace_draws;
} gd_analysis_tree_result;

/*
 * A prepared chronology run: the data is loaded, the concepts are built and the
 * calibrations are resolved, but the sampler has not started. The split exists
 * because both callers want the dataset's size before committing to a run that
 * takes a minute — the CLI prints it, the page shows it beside the progress bar.
 *
 * A job runs once. On a successful run the tree moves into the result, and
 * everything that pointed at it inside the job is torn down at that moment, so
 * the job and the result can be destroyed in either order.
 */
typedef struct gd_analysis_tree_job gd_analysis_tree_job;

int  gd_analysis_tree_prepare(const gd_analysis_tree_input *in,
                              const gd_analysis_tree_config *cfg,
                              gd_analysis_tree_job **out, gd_error *err);
void gd_analysis_tree_job_destroy(gd_analysis_tree_job *job);

/* Known after prepare, before the sampler runs. */
size_t gd_analysis_tree_job_languages(const gd_analysis_tree_job *job);
size_t gd_analysis_tree_job_nodes(const gd_analysis_tree_job *job);
size_t gd_analysis_tree_job_concepts(const gd_analysis_tree_job *job);
/* Form rows naming a language that is not a tip of the tree; they are skipped,
 * and a caller that cares should say so rather than let them vanish. */
size_t gd_analysis_tree_job_skipped_rows(const gd_analysis_tree_job *job);

/* Convenience for reading a result back by tree node rather than by parameter
 * order. Returns NULL when that node's age was not sampled (a fixed tip, or a
 * node the sampler held constant). */
const gd_analysis_node_age *gd_analysis_result_node(
        const gd_analysis_tree_result *res, size_t node);

/*
 * One thinned draw: parameter `param`, chain `chain`, draw `k`. NAN when any
 * index is out of range or the run kept no trace.
 *
 * `trace` and `trace_draws` stay public and readable; this is the stride
 * written once instead of at each call site. The formula documented on the
 * field is still correct -- it is just no longer something a caller has to
 * reproduce to plot a chain.
 */
double gd_analysis_trace_draw(const gd_analysis_tree_result *res,
                              size_t param, size_t chain, size_t k);

int  gd_analysis_tree_job_run(gd_analysis_tree_job *job,
                              gd_analysis_tree_result **out, gd_error *err);

/* Convenience: prepare, run, destroy the job. */
int  gd_analysis_tree_run(const gd_analysis_tree_input *in,
                          const gd_analysis_tree_config *cfg,
                          gd_analysis_tree_result **out, gd_error *err);

void gd_analysis_tree_result_destroy(gd_analysis_tree_result *res);

/* ------------------------------------------------------------- date-pair */

typedef struct gd_analysis_pair_input {
    const char *forms_tsv; size_t forms_len; const char *forms_name;
    const char *lang_a, *lang_b;
    const char *analysis_id;      /* as in gd_analysis_tree_input */
} gd_analysis_pair_input;

/*
 * `model` selects the pairwise model; GD_PAIR_P1 (the default, and the zero
 * value) is the hierarchical constant-rate posterior, and the GD_PAIR_DECAY_*
 * members select an alternative retention law to compare against it.
 *
 * The rate prior is the log replacement rate PER YEAR, the same unit as
 * date-tree (it was per kiloyear before 0.3; subtract log(1000) = 6.9078 from
 * an old value). `t_lo`/`t_hi`, the tip ages and every reported time are in
 * years BP.
 *
 * `age_a`/`age_b` are the sampling dates of the two lists (0 = modern). With
 * dated lists the divergence time T is at least the older of the two, and
 * the lexical evidence measures the lineage separation 2T - age_a - age_b,
 * which the result also reports. For a list that samples the other's own
 * ancestor, that separation is the elapsed time between the two lists.
 * The alternative decay laws assume contemporaneous lists and refuse tip
 * ages.
 */
typedef struct gd_analysis_pair_config {
    double        t_lo;          /* 100 */
    double        t_hi;          /* 8000 */
    double        rate_logmean;  /* -8.9 per year (about -2.0 per kiloyear) */
    double        rate_logsd;    /* 0.6 */
    double        age_a;         /* 0, years BP */
    double        age_b;         /* 0, years BP */
    double        loan_threshold; /* 0.5, as for date-tree */
    gd_pair_model model;         /* GD_PAIR_P1 */
} gd_analysis_pair_config;

/*
 * Map a user-facing model name onto the enum:
 *   NULL or ""       -> GD_PAIR_P1
 *   "p1"             -> GD_PAIR_P1
 *   "swadesh" | "starostin" | "vasilyev" | "two-component"
 *                    -> the matching GD_PAIR_DECAY_*
 *
 * An unrecognised name is GD_E_INVALID_ARG and *out is untouched — a typo in
 * --decay must never quietly select a different model than the one asked for.
 * This lives in the library rather than in each caller so the CLI and the
 * browser cannot drift into accepting different spellings.
 */
int gd_pair_model_from_string(const char *name, gd_pair_model *out,
                              gd_error *err);

/* The stable name of a model, for a manifest or a report. Borrowed static
 * storage; never NULL. */
const char *gd_pair_model_name(gd_pair_model m);

void gd_analysis_pair_config_default(gd_analysis_pair_config *cfg);

typedef struct gd_analysis_pair_result {
    /* The model that produced this result, as gd_pair_model_name spells it.
     * Borrowed static storage: it outlives the result and must not be freed.
     * A copy would buy nothing — the names are a closed set. */
    const char   *model;
    gd_pair_model model_id;
    size_t n_shared, n_cognate;
    double retention;
    size_t n_missing_cognates;   /* rows skipped as unobserved, whole file */
    size_t n_synonym_concepts;   /* shared concepts with more than one form on a side */

    /* Years BP. */
    double mean, median, sd, et_lo_90, et_hi_90;

    /* Lineage separation 2T - age_a - age_b, in years: the total time along
     * the path between the two lists, which is what the cognates measure. */
    double age_a, age_b;
    double sep_median, sep_et_lo_90, sep_et_hi_90;

    bool   absolute_anchored;
    double prior_dominance;

    /* The posterior over divergence time, years BP against density: what the
     * summaries above are computed from, and what a density plot draws. */
    double *grid_t, *grid_density;
    size_t  ngrid;
} gd_analysis_pair_result;

/* ------------------------------------------- calibration sensitivity scoring */

/*
 * Asking what a chronology owes to its calibrations means re-dating it with
 * one withheld and comparing. The comparison is not arithmetic: dropping the
 * only calibration that bounds a node leaves that node's age unidentifiable,
 * and the sampler wanders to values that are not a measurement of anything.
 * Scoring such a run as though it produced a date would be worse than saying
 * nothing.
 *
 * The rule lived in three static functions in cli/commands.c until D6, where
 * it was scientific policy no test could reach and the browser could not use,
 * even though it can already ask for the variant runs.
 */

/* Where a re-dated age stops being a measurement: 100x the largest age the
 * baseline or the calibrations name. Falls back to 100,000 years when the
 * baseline says nothing. */
double gd_analysis_divergence_ceiling(const gd_analysis_tree_result *baseline,
                                      const gd_analysis_calibration *cals,
                                      size_t ncal);

typedef enum gd_sensitivity_outcome {
    GD_SENS_SCORED = 0,   /* the withheld target was re-dated; the row scores */
    GD_SENS_UNBOUNDED,    /* load-bearing: the age ran past the ceiling */
    GD_SENS_ABSENT        /* the target is not a node of this result */
} gd_sensitivity_outcome;

typedef struct gd_sensitivity_row {
    gd_sensitivity_outcome outcome;
    /* Meaningful only when outcome is GD_SENS_SCORED. */
    double median, hpd_lo, hpd_hi;
    bool   covered;       /* the withheld age fell inside the re-dated HPD */
    double abs_err;       /* |re-dated median - withheld age| */
    /* Always meaningful: how far the REST of the chronology moved. Nodes past
     * the ceiling are excluded, because the distance an unidentified age
     * wandered is not a measurement either. */
    double max_other_shift;
} gd_sensitivity_row;

/*
 * Score one withheld calibration. `baseline` and `variant` must be results for
 * the same topology; a mismatch is GD_E_INVALID_ARG rather than a comparison
 * of unrelated nodes.
 */
int gd_analysis_sensitivity_score(const gd_analysis_tree_result *baseline,
                                  const gd_analysis_tree_result *variant,
                                  const gd_analysis_calibration *withheld,
                                  double ceiling,
                                  gd_sensitivity_row *out, gd_error *err);

int  gd_analysis_pair_run(const gd_analysis_pair_input *in,
                          const gd_analysis_pair_config *cfg,
                          gd_analysis_pair_result **out, gd_error *err);
void gd_analysis_pair_result_destroy(gd_analysis_pair_result *res);

#endif /* GD_ANALYSIS_H */

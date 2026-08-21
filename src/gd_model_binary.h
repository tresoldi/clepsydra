/*
 * gd_model_binary.h — binary cognate-set likelihood (T2).
 *
 * The alternative to the strict novel-state model (T1). Each observed cognate
 * class is a binary trait (present / absent) at each language, evolving down
 * the tree under an asymmetric two-state continuous-time Markov chain
 *
 *     0  --lambda_gain-->  1        (an absent class is gained)
 *     1  --lambda_loss-->  0        (a present class is lost)
 *
 * Because a language may carry several cognate classes for one meaning at once,
 * T2 handles synonymy / polymorphism natively — where T1 assumes a single
 * primary cognate per language-meaning. See docs/model-spec.md §3.4.
 *
 * Contrast with T1: the pruning recursion of Felsenstein handles a node with
 * any number of children directly, so T2 needs *no* polytomy resolution —
 * there is no ladder here (cf. gd_model_novel_state.c). Ages are re-read from
 * the tree on every call, so a workspace survives an MCMC chain moving ages.
 *
 * Three model facts are enforced through the API, following Benedict King 2026
 * (doi:10.1371/journal.pcbi.1014312) and docs/paper-evaluation.md:
 *
 *   1. Ascertainment. A cognate class present in *no* sampled language is
 *      never in the input, so the per-class likelihood is conditioned on the
 *      class being observed somewhere: L / (1 - P(absent everywhere observed)).
 *      On by default.
 *   2. Partition rates are weighted by the number of *meanings*, not the
 *      number of cognate classes (GD_T2_WEIGHT_MEANING, the default). The
 *      cognate-count weighting that miscalibrates under ascertainment is
 *      available only as an explicit opt-in for reproducing old analyses.
 *   3. Cognate classes within one meaning are *not* independent (they compete
 *      for the same slot). T2 models them as independent binary characters
 *      anyway; the dataset result raises within_meaning_dependence_warning
 *      whenever any meaning contributes more than one class, and raw T2
 *      marginal likelihoods are never comparable across encodings (T1 vs T2)
 *      without a shared observable space.
 */
#ifndef GD_MODEL_BINARY_H
#define GD_MODEL_BINARY_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"
#include "gd_mcmc_model.h"

#include <stdbool.h>
#include <stddef.h>

/*
 * One binary cognate class: its presence at each node of the tree.
 *   present[i] == 1  present at node i;
 *   present[i] == 0  absent at node i;
 *   present[i] <  0  missing (unknown) — marginalized, contributes factor 1.
 * Internal nodes may be observed (ancient / attested languages).
 */
typedef struct gd_binary_class {
    const int *present;    /* length nnodes; <0 = missing */
    size_t     nnodes;
} gd_binary_class;

/* Relative-rate weighting of meaning partitions (King 2026). */
typedef enum gd_t2_partition_weight {
    GD_T2_WEIGHT_MEANING = 0,  /* default: every meaning counts once */
    GD_T2_WEIGHT_COGNATE = 1   /* opt-in legacy: weight by cognate-class count */
} gd_t2_partition_weight;

/* Prior over the root state. */
typedef enum gd_t2_root_prior {
    GD_T2_ROOT_STATIONARY = 0, /* default: CTMC stationary dist (l_loss,l_gain)/mu */
    GD_T2_ROOT_FIXED      = 1   /* use root_pi1 (P(present) at the root) */
} gd_t2_root_prior;

typedef struct gd_t2_options {
    bool   ascertainment;          /* 1 = condition on observability (default 1) */
    gd_t2_root_prior root_prior;   /* default GD_T2_ROOT_STATIONARY */
    double root_pi1;               /* P(present) at root when root_prior==FIXED */
    gd_t2_partition_weight weight; /* default GD_T2_WEIGHT_MEANING */
} gd_t2_options;

void gd_t2_options_default(gd_t2_options *opt);

/*
 * Map the user-facing partition-weight string
 * ("meaning" | "cognate") onto the enum. NULL/empty selects the default
 * (meaning). An unrecognised string is an error — never guessed.
 */
int gd_t2_weight_from_string(const char *s, gd_t2_partition_weight *out,
                             gd_error *err);

/* Per-class likelihood result and diagnostics. */
typedef struct gd_t2_class_result {
    double log_lik;          /* ascertainment-corrected class log-likelihood */
    double log_lik_raw;      /* before the ascertainment correction */
    double log_p_all_absent; /* log P(absent at every observed node) */
    bool   ascertained;      /* 1 if the ascertainment correction was applied */
    bool   all_absent;       /* 1 if the class is observed present nowhere */
    size_t ops;              /* work counter: nodes visited by the recursion */
    bool   stability_warning;/* 1 if the ascertainment denominator underflowed */
} gd_t2_class_result;

/* Whole-dataset diagnostics (independent of the numeric total). */
typedef struct gd_t2_dataset_result {
    size_t n_classes;
    size_t n_concepts;
    size_t max_classes_per_concept;
    /* 1 if any meaning contributes more than one cognate class, so the
     * (false) within-meaning independence assumption is actually relied on. */
    bool   within_meaning_dependence_warning;
    /* Always 1: a standing reminder that raw T2 likelihoods are not
     * comparable across encodings without a shared observable space. */
    bool   cross_encoding_incomparable;
} gd_t2_dataset_result;

/* ---- Workspace ---- */

typedef struct gd_t2_ws gd_t2_ws;

/*
 * Build a workspace for a fixed topology. Traversals are built here if absent
 * (hence non-const tree). No polytomy resolution is needed. Node ages are
 * re-read from the tree on every likelihood call.
 */
gd_t2_ws *gd_t2_ws_create(gd_tree *tree, gd_error *err);
void      gd_t2_ws_destroy(gd_t2_ws *ws);

/* ---- Transition matrix ---- */

/*
 * Log of the 2x2 transition probability matrix P(t) = exp(Q t) of the
 * asymmetric two-state CTMC, in closed form:
 *   logP[i][j] = log P(state j at the child | state i at the parent),
 *   states 0 = absent, 1 = present, over a branch of duration t >= 0.
 * Requires lambda_gain > 0, lambda_loss > 0, t >= 0.
 */
int gd_t2_transition_logs(double lambda_gain, double lambda_loss, double t,
                          double logP[2][2], gd_error *err);

/* ---- Likelihood ---- */

/*
 * Log-likelihood of one binary cognate class under T2.
 *   ws     — workspace over the tree (ages re-read from the tree).
 *   lambda_gain, lambda_loss — effective gain/loss rates (1/year) already
 *            scaled by any concept/family/branch multiplier by the caller.
 *   cls    — presence/absence, one entry per node of the tree.
 *   opt    — options; NULL means gd_t2_options_default().
 *   res    — per-class result and diagnostics; may be NULL.
 *
 * Returns 0 on success, GD_E_RANGE on a chronology inversion or bad rate,
 * GD_E_INVALID_ARG on a shape mismatch.
 */
int gd_t2_loglik_class(gd_t2_ws *ws, double lambda_gain, double lambda_loss,
                       const gd_binary_class *cls, const gd_t2_options *opt,
                       gd_t2_class_result *res, gd_error *err);

/*
 * Log-likelihood of a whole dataset: the sum over cognate classes, with each
 * class tagged by a meaning (concept) id in [0, nconcepts). Both rates for a
 * class are scaled by concept_rates[concept_ids[k]] (the hierarchy of §4 is
 * applied by the caller; pass all-ones for a strict clock). concept_rates may
 * be NULL, meaning all rates equal 1.
 *
 * `per_class` (may be NULL) receives one gd_t2_class_result per class;
 * `summary` (may be NULL) receives the dataset diagnostics, including the
 * within-meaning dependence warning.
 */
int gd_t2_loglik_dataset(gd_t2_ws *ws,
                         double lambda_gain, double lambda_loss,
                         const gd_binary_class *classes,
                         const size_t *concept_ids, size_t nclasses,
                         const double *concept_rates, size_t nconcepts,
                         const gd_t2_options *opt,
                         gd_t2_class_result *per_class,
                         gd_t2_dataset_result *summary,
                         double *log_lik_total, gd_error *err);

/*
 * Relative-rate normalization weights per meaning (King 2026). Fills
 * weights[c] for c in [0, nconcepts): 1.0 under GD_T2_WEIGHT_MEANING (every
 * meaning counts once, regardless of how many cognate classes it has), or the
 * class count of meaning c under GD_T2_WEIGHT_COGNATE. A mean-one relative-rate
 * constraint (M8/M9) uses these weights; supplying them here keeps the
 * King-2026 default in one place.
 */
int gd_t2_partition_weights(const size_t *concept_ids, size_t nclasses,
                            size_t nconcepts, gd_t2_partition_weight mode,
                            double *weights, gd_error *err);

/* ---- Simulation ---- */

/*
 * Simulate one binary cognate class forward under the CTMC of §3.4: the root
 * is drawn from the prior (stationary by default), then each child's state is
 * drawn from P(t) given its parent's. Fills present[0..nnodes-1] with 0/1.
 * The draw may be all-absent; observed-data simulation rejects those (the
 * ascertainment condition). Used for SBC and posterior-predictive checks.
 */
int gd_t2_simulate_class(gd_rng *rng, const gd_tree *tree,
                         double lambda_gain, double lambda_loss,
                         const gd_t2_options *opt, int *present,
                         gd_error *err);


/* --------------------------------------------------- MCMC adapter ----------
 *
 * T2 seen as a gd_mcmc_model. T2 has been validated against full enumeration
 * since M7 and could reach no sampler until E3, because gd_mcmc_state carried
 * a single scalar log-rate and this likelihood needs two: a gain rate and a
 * loss rate. It reads st->rates[0] and st->rates[1].
 *
 * This is the seam's second adapter, and the reason the seam is real rather
 * than hypothetical.
 */
typedef struct gd_mcmc_t2 gd_mcmc_t2;

/*
 * Wrap a T2 workspace and its classes as a likelihood model. `classes`,
 * `concept_ids` and `opt` must outlive the adapter. `concept_ids[i]` names the
 * meaning class i belongs to, in [0, nconcepts).
 */
gd_mcmc_t2 *gd_mcmc_t2_create(gd_t2_ws *ws,
                              const gd_binary_class *classes,
                              const size_t *concept_ids, size_t nclasses,
                              size_t nconcepts,
                              const gd_t2_options *opt, gd_error *err);
void gd_mcmc_t2_destroy(gd_mcmc_t2 *m);
void gd_mcmc_t2_model(gd_mcmc_t2 *m, gd_mcmc_model *out);

#endif /* GD_MODEL_BINARY_H */

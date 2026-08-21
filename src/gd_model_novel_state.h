/*
 * gd_model_novel_state.h — strict novel-state multistate likelihood (T1).
 *
 * Independent re-implementation of the Rönchen et al. 2024 likelihood
 * algorithm (Theorem 10 + Proposition 7). See docs/model-spec.md §3.3
 * (model), §10 (pseudocode) and §11 (brute-force validation).
 *
 * Generative process, per concept, on a rooted tree with node ages:
 *   - the root carries an arbitrary cognate class;
 *   - along an edge of duration t with concept rate r, the child keeps the
 *     parent class with probability ret = exp(-r t), and otherwise takes a
 *     class that never occurred anywhere in the tree (prob ch = 1 - ret).
 *
 * The observable for a concept is therefore the *partition* of the observed
 * nodes into classes; the concrete label integers are arbitrary. The tree
 * likelihood is a proper probability distribution over those partitions.
 *
 * Structure of this module:
 *   gd_t1_ws          — per-tree workspace; all scratch memory is allocated
 *                       once so the likelihood itself allocates nothing.
 *   gd_t1_extension   — data extension / extendability (debug mode).
 *   gd_t1_loglik_*    — per-concept and whole-dataset likelihoods, including
 *                       the robust contamination mixture (§6).
 *   gd_t1_simulate    — forward simulation for SBC / posterior-predictive.
 *
 * Polytomies are handled by deterministically resolving each polytomy into a
 * ladder of zero-length branches inside the workspace (docs/model-spec.md
 * §3.3, §12 item 10). The exact nonbinary recursion (Thm A.8/A.9) awaits the
 * appendix supplement; see docs/paper-evaluation.md.
 */
#ifndef GD_MODEL_NOVEL_STATE_H
#define GD_MODEL_NOVEL_STATE_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"
#include "gd_mcmc_model.h"

#include <stdbool.h>
#include <stddef.h>

/* Per-concept cognate data: assigns each node of the tree a cognate class.
 * Labels are arbitrary non-negative integers; only their equality structure
 * matters. A label < 0 marks an unobserved (missing) node. Internal nodes may
 * be observed (ancient languages). */
typedef struct gd_cognate_data {
    const int *labels;     /* length nnodes; < 0 = missing */
    size_t     nnodes;
} gd_cognate_data;

/* Policy for concepts whose data cannot be extended (docs/model-spec.md
 * §5.5). No policy silently removes a concept: `exclude` records it. */
typedef enum gd_incompatible_policy {
    GD_INCOMPATIBLE_ERROR         = 0, /* abort with a diagnostic */
    GD_INCOMPATIBLE_EXCLUDE       = 1, /* drop the concept, record a warning */
    GD_INCOMPATIBLE_CONTAMINATION = 2  /* default: L_tree = 0, noise carries it */
} gd_incompatible_policy;

/*
 * Options for the T1 likelihood.
 *
 * The contamination mixture (docs/model-spec.md §6) is
 *   L_c = (1 - eps) * L_tree,c + eps * L_noise,c.
 *
 * The noise component is a Chinese-restaurant process with concentration
 * `noise_alpha` over the partition of the *observed* nodes — the K -> inf
 * limit of the symmetric Dirichlet-multinomial the spec allows, and the only
 * form that is normalized over exactly the observable space of the tree
 * component (a distribution over partitions, not over labelings). It is an
 * error model, not a claim that the real event was random.
 */
typedef struct gd_t1_options {
    gd_incompatible_policy incompatible; /* default CONTAMINATION */
    double epsilon;                      /* contamination weight; default 0.05
                                          * (mean of the Beta(1,19) prior) */
    double noise_alpha;                  /* CRP concentration; default 1.0 */
} gd_t1_options;

void gd_t1_options_default(gd_t1_options *opt);

/*
 * Map the user-facing incompatible-concept setting
 * ("error" | "exclude" | "contamination") onto the policy. NULL or an empty
 * string selects the default. An unrecognised string is an error — the
 * setting is never guessed.
 */
int gd_t1_policy_from_string(const char *s, gd_incompatible_policy *out,
                             gd_error *err);

/* Per-concept likelihood result and diagnostics. */
typedef struct gd_t1_concept_result {
    double log_lik;          /* concept contribution after the mixture */
    double log_lik_tree;     /* tree component; -HUGE_VAL if incompatible */
    double log_lik_noise;    /* noise component */
    double p_contaminated;   /* posterior P(concept came from the noise part) */
    bool   compatible;       /* 1 if the data is extendable */
    bool   excluded;         /* 1 if dropped under GD_INCOMPATIBLE_EXCLUDE */
    /* Work counters, for the scaling acceptance test. `ops_recursion` is
     * exactly the number of working nodes; `ops_extension` counts the node
     * visits of the data extension. */
    size_t ops_extension;
    size_t ops_recursion;
    /* 1 if a p - q*ret cancellation went measurably negative (should not
     * happen; indicates loss of precision rather than a wrong answer). */
    bool   stability_warning;
} gd_t1_concept_result;

/* Extended-state constraints, for the debug mode (docs/model-spec.md §10).
 * `determined` has one entry per node of the *original* tree: the class the
 * node is determined to, or -1 if undetermined. Storage is owned by the
 * workspace and is valid until the next call on that workspace. */
typedef struct gd_t1_extension_info {
    const int *determined;
    size_t     nnodes;
    size_t     n_determined;
    bool       extendable;
    size_t     conflict_node;     /* GD_NODE_NONE if extendable */
    int        conflict_class_a;  /* the two classes in conflict; -1 if none */
    int        conflict_class_b;
} gd_t1_extension_info;

/* ---- Workspace ---- */

typedef struct gd_t1_ws gd_t1_ws;

/*
 * Build a workspace for a fixed tree topology. Polytomies are resolved into
 * zero-length ladders here, once. Node *ages* are re-read from the tree on
 * every likelihood call, so the workspace stays valid while an MCMC chain
 * moves node ages; it must be rebuilt if the topology changes.
 *
 * The tree needs LCA support; it is built here if absent (hence non-const).
 */
gd_t1_ws *gd_t1_ws_create(gd_tree *tree, gd_error *err);
void      gd_t1_ws_destroy(gd_t1_ws *ws);

/* Number of nodes in the internal binary working tree (>= tree->nnodes; the
 * difference is the number of nodes inserted to resolve polytomies). */
size_t gd_t1_ws_working_nodes(const gd_t1_ws *ws);

/* ---- Data extension (Proposition 7) ---- */

/*
 * Compute the extension f~ of the concept data: every node determined by a
 * pair of equally-labelled observed nodes gets that class. Fills `info`
 * (debug mode: this is the extended-state constraint set).
 *
 * Returns 0 on success — including when the data is *not* extendable, which
 * is reported through info->extendable rather than as an error.
 */
int gd_t1_extension(gd_t1_ws *ws, const gd_cognate_data *data,
                    gd_t1_extension_info *info, gd_error *err);

/* Convenience: 1 if the concept data is extendable (compatible), else 0. */
bool gd_t1_is_compatible(gd_t1_ws *ws, const gd_cognate_data *data);

/* ---- Likelihood ---- */

/*
 * Log-likelihood of one concept.
 *
 *   ws    — workspace over the tree (ages are re-read from the tree).
 *   rate  — the concept's replacement rate on every branch (1/year). Branch
 *           and family rate multipliers are applied by the caller (M9).
 *   data  — cognate labels, one per node of the original tree.
 *   opt   — options; NULL means gd_t1_options_default().
 *   res   — per-concept result and diagnostics; may be NULL.
 *
 * Returns 0 on success, GD_E_INCOMPATIBLE if the data is not extendable and
 * the policy is GD_INCOMPATIBLE_ERROR (the message names the conflicting
 * node and the two classes), or another GD_E_* on invalid input.
 */
int gd_t1_loglik_concept(gd_t1_ws *ws, double rate,
                         const gd_cognate_data *data,
                         const gd_t1_options *opt,
                         gd_t1_concept_result *res, gd_error *err);

/*
 * Relaxed-clock variant: the rate on the edge above child node w is scaled by
 * branch_mult[w] (a per-branch multiplier m_b, indexed by the child node of
 * the branch; the root entry is unused). branch_mult == NULL is the strict
 * clock, identical to gd_t1_loglik_concept. Under docs/model-spec.md §4 the
 * effective per-branch rate for the concept is `rate * m_b`; the multipliers
 * are shared across concepts and supplied by the sampler (M9).
 */
int gd_t1_loglik_concept_relaxed(gd_t1_ws *ws, double rate,
                                 const double *branch_mult,
                                 const gd_cognate_data *data,
                                 const gd_t1_options *opt,
                                 gd_t1_concept_result *res, gd_error *err);

/*
 * Log-likelihood of a whole dataset: the sum over concepts, with per-concept
 * rates. Concepts excluded under GD_INCOMPATIBLE_EXCLUDE contribute 0 and are
 * flagged in `per_concept`. `per_concept` may be NULL; it must otherwise have
 * `nconcepts` entries.
 */
int gd_t1_loglik_dataset(gd_t1_ws *ws, const double *rates,
                         const gd_cognate_data *concepts, size_t nconcepts,
                         const gd_t1_options *opt,
                         gd_t1_concept_result *per_concept,
                         double *log_lik_total, gd_error *err);

/*
 * Relaxed-clock dataset variant: per-concept `rates` with a shared per-branch
 * multiplier array `branch_mult` (indexed by child node; NULL = strict clock).
 */
int gd_t1_loglik_dataset_relaxed(gd_t1_ws *ws, const double *rates,
                                 const double *branch_mult,
                                 const gd_cognate_data *concepts, size_t nconcepts,
                                 const gd_t1_options *opt,
                                 gd_t1_concept_result *per_concept,
                                 double *log_lik_total, gd_error *err);

/*
 * Noise component: log CRP(alpha) probability of the partition induced by the
 * observed labels. Normalized over all partitions of the observed nodes, i.e.
 * over exactly the observable space of the tree component.
 */
double gd_t1_noise_loglik(const gd_cognate_data *data, double alpha);

/* ---- Simulation ---- */

/*
 * Simulate one concept forward under the generative process of §3.3.
 * Fills labels[0..nnodes-1] with class IDs >= 0 (root = 0; every replacement
 * opens a new class). Used for SBC and posterior-predictive checks.
 */
int gd_t1_simulate(gd_rng *rng, const gd_tree *tree, double rate,
                   int *labels, gd_error *err);

/* ---- T1 as a gd_mcmc_model ----
 *
 * Declared here, with the likelihood, rather than in gd_mcmc.h. Until E3 the
 * sampler's header included THIS one so it could declare this adapter, which
 * is the compile-time dependency on T1 that docs/c99-alignment-plan.md
 * §Deferred names. The sampler no longer knows which likelihoods exist.
 */

typedef struct gd_mcmc_t1 gd_mcmc_t1;

/*
 * Wrap a T1 workspace and concept dataset as a likelihood model. Concept rates
 * are exp(logr + mu_c); branch multipliers (when sampled) scale each branch;
 * the contamination epsilon overrides opt->epsilon per evaluation. `concepts`
 * and `opt` must outlive the adapter. Set opt->incompatible = CONTAMINATION when
 * sampling contamination so incompatible/loan concepts are absorbed, not errors.
 */
gd_mcmc_t1 *gd_mcmc_t1_create(gd_t1_ws *ws, const gd_cognate_data *concepts,
                              size_t nconcepts, const gd_t1_options *opt,
                              gd_error *err);
void gd_mcmc_t1_destroy(gd_mcmc_t1 *m);
void gd_mcmc_t1_model(gd_mcmc_t1 *m, gd_mcmc_model *out);

#endif /* GD_MODEL_NOVEL_STATE_H */

/*
 * gd_types.h — types shared across the Clepsydra public API.
 *
 * Strict C99. Depends only on <stddef.h> and <stdint.h>.
 *
 * This header holds the small vocabulary that more than one public header
 * needs. Types belonging to a single module are declared with that module:
 * results and configuration in gd_analysis.h, the topology in gd_tree.h,
 * error codes in gd_error.h.
 */
#ifndef GD_TYPES_H
#define GD_TYPES_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Longest error message the library will record, including the terminator.
 * Anything longer is truncated. */
#define GD_ERROR_MSG_MAX 256

/*
 * Error reporting. The library never calls exit() or abort(); every operation
 * that can fail takes a gd_error* and reports through it. The full convention
 * — what an int return means, when a predicate returns bool, what assert() is
 * for — is documented at the top of gd_error.h.
 *
 * The error OWNS its message: gd_error_set copies, and gd_error_setf formats
 * straight into this storage. So a message may name a runtime value without
 * the caller supplying a buffer to hold it, and it stays readable for exactly
 * as long as the gd_error does. `message` is never NULL and is always
 * NUL-terminated, so it is always safe to print; it is "" when unset.
 *
 * Initialise with GD_ERROR_INIT (gd_error.h), never a positional brace list.
 */
typedef struct gd_error {
    int  code;                        /* GD_E_* from gd_error.h; 0 == OK */
    char message[GD_ERROR_MSG_MAX];   /* owned; "" when unset */
} gd_error;

/*
 * Kinds of sampled parameter. Public because a reported result labels each
 * parameter with one: a caller reading gd_analysis_param needs the same
 * vocabulary the sampler writes.
 */
typedef enum gd_mcmc_ptype {
    GD_MCMC_P_LOGRATE = 0,       /* global log clock-rate */
    GD_MCMC_P_AGE,               /* a node age */
    GD_MCMC_P_LOG_SIGMA_B,       /* log of the branch-clock SD */
    GD_MCMC_P_LOGIT_EPS,         /* logit of the contamination weight */
    GD_MCMC_P_BRANCH_LOGMULT,    /* log of a branch multiplier m_b */
    GD_MCMC_P_CONCEPT_DEV,       /* a concept log-rate deviation mu_c */
    GD_MCMC_P_LOG_SIGMA_C,       /* log of the among-concept rate spread */
    GD_MCMC_P_LOGIT_TIPERR       /* logit of a tip's observation-error probability */
} gd_mcmc_ptype;

/*
 * How a calibration constrains its target's age (docs/calibrations.md §1).
 * Public because a caller enumerating a calibration table with
 * gd_analysis_calibrations_read reads it back.
 */
typedef enum gd_cal_type {
    GD_CAL_FIXED_TIP = 0,
    GD_CAL_UNIFORM_NODE,
    GD_CAL_NORMAL_NODE,
    GD_CAL_LOGNORMAL_NODE,
    GD_CAL_OFFSET_LOGNORMAL_NODE,
    GD_CAL_HARD_MINIMUM,
    GD_CAL_HARD_MAXIMUM,
    GD_CAL_SOFT_MINIMUM,
    GD_CAL_SOFT_MAXIMUM,
    /* A dated tip that samples its sister lineage's own ancestor: the tip
     * is fixed at p1 and its parent is pinned to the same age, so the tip
     * hangs on a zero-length branch (roadmap E2.7). */
    GD_CAL_ANCESTOR_TIP,
    /* Not an age: a prior on a tip's observation-error probability (roadmap
     * E2.9), Beta with mean p1 and concentration p2. The tip's age is
     * whatever the other rows (or the default, 0) make it. */
    GD_CAL_TIP_ERROR
} gd_cal_type;

/*
 * Which model dates a language pair.
 *
 * GD_PAIR_P1 is the default and is deliberately 0, so a zeroed configuration
 * asks for the model the project actually recommends. The four decay laws are
 * alternatives to compare against it, not improvements on it — see
 * docs/model-cards.md on why "more flexible" is not "better".
 */
typedef enum gd_pair_model {
    GD_PAIR_P1 = 0,              /* hierarchical constant-rate (default) */
    GD_PAIR_DECAY_SWADESH,       /* constant rate, no concept hierarchy */
    GD_PAIR_DECAY_STAROSTIN,     /* replacement appears to accelerate */
    GD_PAIR_DECAY_VASILYEV,      /* single curve with a retention floor */
    GD_PAIR_DECAY_TWO_COMPONENT  /* stable core + regular vocabulary */
} gd_pair_model;

/*
 * Deterministic PCG32 generator with independent streams. Opaque: the state
 * layout is an implementation detail, and a serialized stream is the supported
 * way to carry one across a run boundary.
 */
struct gd_rng;
typedef struct gd_rng gd_rng;

#ifdef __cplusplus
}
#endif
#endif /* GD_TYPES_H */

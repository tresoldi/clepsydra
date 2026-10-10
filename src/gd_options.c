/*
 * gd_options.c — the option tables, and the type check that makes them safe.
 *
 * The interesting part is GD_OPT_ROW. Each row pairs offsetof(T, field) with
 * a gd_opt_type, and the two must agree or the parser writes the wrong number
 * of bytes into the config. C99 has no _Static_assert, so the check is the
 * classic negative-array-size trick: sizeof a char array whose length is 1
 * when the field's size matches the type's and -1 when it does not. A
 * mismatch is a compile error in an ordinary `make`, not a corruption found
 * later by a browser test.
 */
#include "gd_options.h"
#include "gd_analysis.h"
#include "gd_tsv.h"
#include "gd_math.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Evaluates to 0, and fails to compile when `cond` is false. */
#define GD_SIZE_CHECK(cond) \
    (0 * (int)sizeof(char[(cond) ? 1 : -1]))

/* `name_` is the canonical key -- what the WebAssembly shim accepts in a
 * key=value line and what the page sends. It is stated rather than derived
 * from the field name: the shim's key set is a published interface (see
 * docs/web.md) and predates this table, so `warmup` must not silently become
 * `n_warmup`. */
#define GD_OPT_ROW(T, field, gdtype, ctype, name_, flag_, lo_, hi_, loopen_, hiopen_, help_) \
    GD_OPT_ROW_NOTE(T, field, gdtype, ctype, name_, flag_, lo_, hi_, loopen_, hiopen_, NULL, help_)

#define GD_OPT_ROW_NOTE(T, field, gdtype, ctype, name_, flag_, lo_, hi_, loopen_, hiopen_, note_, help_) \
    { name_, flag_, gdtype,                                                 \
      offsetof(T, field) + (size_t)GD_SIZE_CHECK(                            \
          sizeof(((T *)0)->field) == sizeof(ctype)),                         \
      lo_, hi_, loopen_, hiopen_, note_, help_ }

#define GD_INF HUGE_VAL

/* ------------------------------------------------------------ tree options */

static const gd_option TREE_OPTIONS[] = {
    GD_OPT_ROW(gd_analysis_tree_config, n_warmup, GD_OPT_SIZE, size_t, "warmup",
        "--warmup", 0.0, GD_INF, false, false,
        "warm-up iterations per chain, discarded"),
    GD_OPT_ROW(gd_analysis_tree_config, n_samples, GD_OPT_SIZE, size_t, "samples",
        "--samples", 0.0, GD_INF, true, false,
        "retained iterations per chain"),
    GD_OPT_ROW(gd_analysis_tree_config, n_chains, GD_OPT_SIZE, size_t, "chains",
        "--chains", 0.0, GD_INF, true, false,
        "independent chains; R-hat needs at least 2"),
    GD_OPT_ROW(gd_analysis_tree_config, threads, GD_OPT_SIZE, size_t, "threads",
        "--threads", 0.0, GD_INF, false, false,
        "chains run at once (0 = all); draws do not depend on it"),
    GD_OPT_ROW(gd_analysis_tree_config, seed, GD_OPT_U64, uint64_t, "seed",
        "--seed", 0.0, GD_INF, false, false,
        "RNG seed; a seeded run is bit-reproducible"),
    GD_OPT_ROW(gd_analysis_tree_config, rate_logmean, GD_OPT_DOUBLE, double, "rate_logmean",
        "--rate-logmean", -GD_INF, GD_INF, false, false,
        "prior mean of the log clock rate (-8.80 = Swadesh/yr)"),
    GD_OPT_ROW(gd_analysis_tree_config, rate_logsd, GD_OPT_DOUBLE, double, "rate_logsd",
        "--rate-logsd", 0.0, GD_INF, true, false,
        "prior sd of the log clock rate; must be positive"),
    GD_OPT_ROW(gd_analysis_tree_config, epsilon, GD_OPT_DOUBLE, double, "epsilon",
        "--epsilon", 0.0, 1.0, false, true,
        "contamination weight in [0, 1), used with --sample-epsilon 0"),
    GD_OPT_ROW(gd_analysis_tree_config, sample_epsilon, GD_OPT_BOOL, bool, "sample_epsilon",
        "--sample-epsilon", 0.0, 1.0, false, false,
        "sample the contamination weight, Beta(1, 19) prior"),
    GD_OPT_ROW(gd_analysis_tree_config, relaxed, GD_OPT_BOOL, bool, "relaxed",
        "--relaxed", 0.0, 1.0, false, false,
        "relaxed branch clock (M9) instead of a strict one"),
    GD_OPT_ROW(gd_analysis_tree_config, block_age_scaling, GD_OPT_BOOL, bool, "block_age_scaling",
        "--block-age-scaling", 0.0, 1.0, false, false,
        "decorrelate nested ages; helps on deep trees"),
    GD_OPT_ROW(gd_analysis_tree_config, root_lex_logmean, GD_OPT_DOUBLE, double, "root_lex_logmean",
        "--root-lex-logmean", -GD_INF, GD_INF, false, false,
        "prior mean of log(rate x root age): replacements per concept, root to tip"),
    GD_OPT_ROW(gd_analysis_tree_config, root_lex_logsd, GD_OPT_DOUBLE, double, "root_lex_logsd",
        "--root-lex-logsd", 0.0, GD_INF, true, false,
        "prior sd of that log lexical root; must be positive"),
    GD_OPT_ROW(gd_analysis_tree_config, prior_only, GD_OPT_BOOL, bool, "prior_only",
        "--prior-only", 0.0, 1.0, false, false,
        "sample the prior only (likelihood off), to compare against"),
    GD_OPT_ROW(gd_analysis_tree_config, concept_rates, GD_OPT_BOOL, bool, "concept_rates",
        "--concept-rates", 0.0, 1.0, false, false,
        "integrate out among-concept rate variation (0 = one rate for all)"),
    GD_OPT_ROW(gd_analysis_tree_config, concept_rate_prior_logmean, GD_OPT_DOUBLE, double, "concept_rate_prior_logmean",
        "--concept-rate-prior-logmean", -GD_INF, GD_INF, false, false,
        "prior mean of log sigma_c, the concept-rate spread (empirical: -0.1)"),
    GD_OPT_ROW(gd_analysis_tree_config, concept_rate_prior_logsd, GD_OPT_DOUBLE, double, "concept_rate_prior_logsd",
        "--concept-rate-prior-logsd", 0.0, GD_INF, true, false,
        "prior sd of log sigma_c (empirical: 0.35)"),
    GD_OPT_ROW(gd_analysis_tree_config, concept_rate_gamma, GD_OPT_BOOL, bool, "concept_rate_gamma",
        "--concept-rate-gamma", 0.0, 1.0, false, false,
        "gamma instead of log-normal concept rates (sensitivity; heavier slow tail)"),
    GD_OPT_ROW(gd_analysis_tree_config, loan_threshold, GD_OPT_DOUBLE, double, "loan_threshold",
        "--loan-threshold", 0.0, GD_INF, false, false,
        "forms at or above this loan_probability are read as missing; >1 keeps all"),
    GD_OPT_ROW(gd_analysis_tree_config, yule, GD_OPT_BOOL, bool, "yule",
        "--yule", 0.0, 1.0, false, false,
        "Yule age prior instead of uniform (sensitivity; modern tips only)"),
    GD_OPT_ROW(gd_analysis_tree_config, yule_kappa, GD_OPT_DOUBLE, double, "yule_kappa",
        "--yule-kappa", 0.0, GD_INF, false, false,
        "Yule birth rate x root age; 0 = the expected value for the tip count"),
    GD_OPT_ROW(gd_analysis_tree_config, keep_pointwise, GD_OPT_BOOL, bool, "keep_pointwise",
        "--write-pointwise", 0.0, 1.0, false, false,
        "write per-draw, per-concept log-likelihoods (pointwise.tsv) for LOO/WAIC"),
    GD_OPT_ROW(gd_analysis_tree_config, compare_prior, GD_OPT_BOOL, bool, "compare_prior",
        "--compare-prior", 0.0, 1.0, false, false,
        "also sample the prior alone and report it beside each node"),
    GD_OPT_ROW(gd_analysis_tree_config, hpd_level, GD_OPT_DOUBLE, double, "hpd_level",
        "--hpd-level", 0.0, 1.0, true, true,
        "credible level for the reported interval"),
    GD_OPT_ROW(gd_analysis_tree_config, progress_every, GD_OPT_SIZE, size_t, "progress_every",
        "--progress-every", 0.0, GD_INF, false, false,
        "progress interval in iterations; 0 is about every 1%"),
    GD_OPT_ROW(gd_analysis_tree_config, trace_draws, GD_OPT_SIZE, size_t, "trace_draws",
        "--trace-draws", 0.0, GD_INF, false, false,
        "evenly spaced draws kept per chain per parameter"),
    GD_OPT_ROW_NOTE(gd_analysis_tree_config, calibrations.drop_index, GD_OPT_SIZE, size_t, "cal_drop",
        "--cal-drop", 0.0, GD_INF, false, false, "keep all",
        "withhold this calibration row, by index (leave-one-out)"),
    GD_OPT_ROW(gd_analysis_tree_config, calibrations.widen, GD_OPT_DOUBLE, double, "cal_widen",
        "--cal-widen", 0.0, GD_INF, true, false,
        "multiply every soft calibration window by this"),
};

/* ------------------------------------------------------------ pair options */

static const gd_option PAIR_OPTIONS[] = {
    GD_OPT_ROW(gd_analysis_pair_config, t_lo, GD_OPT_DOUBLE, double, "t_lo",
        "--t-lo", 0.0, GD_INF, false, false,
        "lower bound of the divergence-time grid, years BP"),
    GD_OPT_ROW(gd_analysis_pair_config, t_hi, GD_OPT_DOUBLE, double, "t_hi",
        "--t-hi", 0.0, GD_INF, true, false,
        "upper bound of the divergence-time grid, years BP"),
    GD_OPT_ROW(gd_analysis_pair_config, rate_logmean, GD_OPT_DOUBLE, double, "rate_logmean",
        "--rate-logmean", -GD_INF, GD_INF, false, false,
        "prior mean of the log rate, per year (as date-tree)"),
    GD_OPT_ROW(gd_analysis_pair_config, rate_logsd, GD_OPT_DOUBLE, double, "rate_logsd",
        "--rate-logsd", 0.0, GD_INF, true, false,
        "prior sd of the log rate; must be positive"),
    GD_OPT_ROW(gd_analysis_pair_config, loan_threshold, GD_OPT_DOUBLE, double, "loan_threshold",
        "--loan-threshold", 0.0, GD_INF, false, false,
        "forms at or above this loan_probability are read as missing; >1 keeps all"),
    GD_OPT_ROW(gd_analysis_pair_config, age_a, GD_OPT_DOUBLE, double, "age_a",
        "--age-a", 0.0, GD_INF, false, false,
        "sampling date of language A's list, years BP (0 = modern)"),
    GD_OPT_ROW(gd_analysis_pair_config, age_b, GD_OPT_DOUBLE, double, "age_b",
        "--age-b", 0.0, GD_INF, false, false,
        "sampling date of language B's list, years BP (0 = modern)"),
};

const gd_option *gd_options_tree(size_t *n)
{
    if (n != NULL) *n = sizeof TREE_OPTIONS / sizeof TREE_OPTIONS[0];
    return TREE_OPTIONS;
}

const gd_option *gd_options_pair(size_t *n)
{
    if (n != NULL) *n = sizeof PAIR_OPTIONS / sizeof PAIR_OPTIONS[0];
    return PAIR_OPTIONS;
}

const gd_option *gd_option_find(const gd_option *opts, size_t n,
                                const char *name, size_t name_len)
{
    if (opts == NULL || name == NULL) return NULL;
    if (name_len == GD_OPT_MEASURE) name_len = strlen(name);
    for (size_t i = 0; i < n; i++) {
        if (strlen(opts[i].name) == name_len &&
            memcmp(opts[i].name, name, name_len) == 0) {
            return &opts[i];
        }
    }
    return NULL;
}

const gd_option *gd_option_find_flag(const gd_option *opts, size_t n,
                                     const char *flag)
{
    if (opts == NULL || flag == NULL) return NULL;
    for (size_t i = 0; i < n; i++) {
        if (strcmp(opts[i].flag, flag) == 0) return &opts[i];
    }
    return NULL;
}

/* Where this option's field sits inside a given config. */
static void *field_of(const gd_option *opt, void *config)
{
    return (void *)((char *)config + opt->offset);
}

static const void *field_of_const(const gd_option *opt, const void *config)
{
    return (const void *)((const char *)config + opt->offset);
}

int gd_option_set(const gd_option *opt, void *config,
                  const char *value, size_t value_len, gd_error *err)
{
    if (opt == NULL || config == NULL || value == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "option: null argument");
    }
    if (value_len == GD_OPT_MEASURE) value_len = strlen(value);

    /*
     * Integer-typed options are parsed as integers, not through a double.
     * The table this replaced sent every value through strtod, which silently
     * mangles anything past 2^53 -- a seed, or the size_t sentinel that means
     * "keep every calibration", neither of which survives the trip.
     */
    bool is_int = (opt->type == GD_OPT_SIZE || opt->type == GD_OPT_U64);
    unsigned long long uv = 0;
    double v = 0.0;

    if (is_int) {
        char buf[32];
        if (value_len == 0 || value_len >= sizeof buf) {
            return gd_error_setf(err, GD_E_BAD_NUMERIC,
                                 "option '%s' is not a whole number: %.*s",
                                 opt->name, (int)value_len, value);
        }
        memcpy(buf, value, value_len);
        buf[value_len] = '\0';
        const char *p = buf;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '-' || *p == '+') {
            return gd_error_setf(err, GD_E_BAD_NUMERIC,
                                 "option '%s' is not a whole number: %s",
                                 opt->name, buf);
        }
        char *end = NULL;
        errno = 0;
        uv = strtoull(p, &end, 10);
        if (end == p || end == NULL || *end != '\0' || errno == ERANGE) {
            return gd_error_setf(err, GD_E_BAD_NUMERIC,
                                 "option '%s' is not a whole number: %s",
                                 opt->name, buf);
        }
        v = (double)uv;
    } else {
        if (gd_tsv_parse_double(value, value_len, &v) != GD_E_OK || !gd_is_finite(v)) {
            return gd_error_setf(err, GD_E_BAD_NUMERIC,
                                 "option '%s' is not a number: %.*s",
                                 opt->name, (int)value_len, value);
        }
    }

    bool below = opt->lo_open ? !(v > opt->lo) : !(v >= opt->lo);
    bool above = opt->hi_open ? !(v < opt->hi) : !(v <= opt->hi);
    if (below || above) {
        if (opt->hi == GD_INF) {
            return gd_error_setf(err, GD_E_INVALID_ARG,
                                 "option '%s' must be %s %g",
                                 opt->name, opt->lo_open ? "greater than" : "at least",
                                 opt->lo);
        }
        return gd_error_setf(err, GD_E_INVALID_ARG,
                             "option '%s' must be in %s%g, %g%s",
                             opt->name, opt->lo_open ? "(" : "[", opt->lo, opt->hi,
                             opt->hi_open ? ")" : "]");
    }

    void *slot = field_of(opt, config);
    switch (opt->type) {
    case GD_OPT_SIZE:   *(size_t *)slot   = (size_t)uv;   break;
    case GD_OPT_U64:    *(uint64_t *)slot = (uint64_t)uv; break;
    case GD_OPT_DOUBLE: *(double *)slot   = v;           break;
    case GD_OPT_BOOL:   *(bool *)slot     = (v != 0.0);  break;
    }
    return GD_E_OK;
}

int gd_option_format(const gd_option *opt, const void *config,
                     char *buf, size_t bufsz, gd_error *err)
{
    if (opt == NULL || config == NULL || buf == NULL || bufsz == 0) {
        return gd_error_set(err, GD_E_INVALID_ARG, "option: null argument");
    }
    const void *slot = field_of_const(opt, config);
    int wrote = -1;
    switch (opt->type) {
    case GD_OPT_SIZE:
        wrote = snprintf(buf, bufsz, "%zu", *(const size_t *)slot);
        break;
    case GD_OPT_U64:
        wrote = snprintf(buf, bufsz, "%llu",
                         (unsigned long long)*(const uint64_t *)slot);
        break;
    case GD_OPT_DOUBLE:
        wrote = snprintf(buf, bufsz, "%g", *(const double *)slot);
        break;
    case GD_OPT_BOOL:
        wrote = snprintf(buf, bufsz, "%d", *(const bool *)slot ? 1 : 0);
        break;
    }
    if (wrote < 0 || (size_t)wrote >= bufsz) {
        return gd_error_setf(err, GD_E_RANGE,
                             "option '%s': value does not fit", opt->name);
    }
    return GD_E_OK;
}

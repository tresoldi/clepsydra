/*
 * test_options.c — the run-option registry.
 *
 * The registry exists so the CLI, the WebAssembly shim and the page cannot
 * drift apart on what an option is called, what it defaults to, or what
 * values it accepts. These cases test the properties that would let them
 * drift again.
 */
#include "test_runner.h"
#include "gd_options.h"
#include "gd_analysis.h"

#include <math.h>
#include <string.h>

TEST(test_options_every_row_round_trips)
{
    /* set(format(x)) == x for every option, which is what lets a surface show
     * a default by reading it back instead of restating the number. */
    gd_analysis_tree_config a, b;
    gd_analysis_tree_config_default(&a);
    gd_analysis_tree_config_default(&b);
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    CHECK(n > 0);

    for (size_t i = 0; i < n; i++) {
        char buf[64];
        gd_error err = GD_ERROR_INIT;
        CHECK(gd_option_format(&opts[i], &a, buf, sizeof buf, &err) == GD_E_OK);
        CHECK(gd_option_set(&opts[i], &b, buf, GD_OPT_MEASURE, &err) == GD_E_OK);
    }
    /* Nothing moved: b was defaulted, then rewritten with a's own values. */
    CHECK(memcmp(&a, &b, sizeof a) == 0);
}

TEST(test_options_pair_rows_round_trip)
{
    gd_analysis_pair_config a, b;
    gd_analysis_pair_config_default(&a);
    gd_analysis_pair_config_default(&b);
    size_t n;
    const gd_option *opts = gd_options_pair(&n);
    CHECK(n > 0);
    for (size_t i = 0; i < n; i++) {
        char buf[64];
        gd_error err = GD_ERROR_INIT;
        CHECK(gd_option_format(&opts[i], &a, buf, sizeof buf, &err) == GD_E_OK);
        CHECK(gd_option_set(&opts[i], &b, buf, GD_OPT_MEASURE, &err) == GD_E_OK);
    }
    CHECK(memcmp(&a, &b, sizeof a) == 0);
}

TEST(test_options_names_and_flags_are_unique)
{
    /* Two rows sharing a name would make one of them unreachable, and the
     * lookup would silently pick the first. */
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            CHECK(strcmp(opts[i].name, opts[j].name) != 0);
            CHECK(strcmp(opts[i].flag, opts[j].flag) != 0);
        }
        /* Every flag is a long option; every name is bare. */
        CHECK(opts[i].flag[0] == '-' && opts[i].flag[1] == '-');
        CHECK(opts[i].name[0] != '-');
        CHECK(opts[i].help != NULL && opts[i].help[0] != '\0');
    }
}

TEST(test_options_offsets_land_inside_their_config)
{
    /* A row whose offset ran past its struct would corrupt the stack. The
     * sizeof-assert in the table catches a type mismatch at compile time;
     * this catches a row attached to the wrong struct entirely. */
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    for (size_t i = 0; i < n; i++) {
        CHECK(opts[i].offset < sizeof(gd_analysis_tree_config));
    }
    opts = gd_options_pair(&n);
    for (size_t i = 0; i < n; i++) {
        CHECK(opts[i].offset < sizeof(gd_analysis_pair_config));
    }
}

TEST(test_options_lookup_by_name_and_flag)
{
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    const gd_option *chains = gd_option_find(opts, n, "chains", GD_OPT_MEASURE);
    CHECK(chains != NULL);
    CHECK_STREQ(chains->flag, "--chains");
    /* A prefix must not match: "chain" is not "chains". */
    CHECK(gd_option_find(opts, n, "chain", GD_OPT_MEASURE) == NULL);
    /* An explicit length is honoured, so a caller can pass a slice of a
     * key=value line without copying it. */
    CHECK(gd_option_find(opts, n, "chains=4", 6) == chains);
    CHECK(gd_option_find_flag(opts, n, "--chains") == chains);
    CHECK(gd_option_find_flag(opts, n, "--nonesuch") == NULL);
}

TEST(test_options_reject_out_of_range)
{
    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    gd_error err = GD_ERROR_INIT;

    /* chains must be positive: the lower bound is exclusive. */
    const gd_option *chains = gd_option_find(opts, n, "chains", GD_OPT_MEASURE);
    CHECK(gd_option_set(chains, &cfg, "0", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
    CHECK(strstr(err.message, "chains") != NULL);
    CHECK(gd_option_set(chains, &cfg, "1", GD_OPT_MEASURE, &err) == GD_E_OK);

    /* epsilon is [0, 1): 0 is allowed, 1 is not. This is the bound the
     * browser did not have before the registry — the CLI refused an epsilon
     * of 1 and the page passed it straight to the sampler. */
    const gd_option *eps = gd_option_find(opts, n, "epsilon", GD_OPT_MEASURE);
    CHECK(gd_option_set(eps, &cfg, "0", GD_OPT_MEASURE, &err) == GD_E_OK);
    CHECK(gd_option_set(eps, &cfg, "1", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_set(eps, &cfg, "-0.1", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);

    /* rate_logsd must be positive, and unbounded above. */
    const gd_option *sd = gd_option_find(opts, n, "rate_logsd", GD_OPT_MEASURE);
    CHECK(gd_option_set(sd, &cfg, "0", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_set(sd, &cfg, "1e6", GD_OPT_MEASURE, &err) == GD_E_OK);

    /* rate_logmean is unbounded on both sides, negative included. */
    const gd_option *mu = gd_option_find(opts, n, "rate_logmean", GD_OPT_MEASURE);
    CHECK(gd_option_set(mu, &cfg, "-8.8", GD_OPT_MEASURE, &err) == GD_E_OK);
}

TEST(test_options_reject_non_numbers)
{
    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    const gd_option *seed = gd_option_find(opts, n, "seed", GD_OPT_MEASURE);
    gd_error err = GD_ERROR_INIT;

    CHECK(gd_option_set(seed, &cfg, "abc", GD_OPT_MEASURE, &err) == GD_E_BAD_NUMERIC);
    CHECK(strstr(err.message, "seed") != NULL);
    CHECK(gd_option_set(seed, &cfg, "12x", GD_OPT_MEASURE, &err) == GD_E_BAD_NUMERIC);
    CHECK(gd_option_set(seed, &cfg, "", GD_OPT_MEASURE, &err) == GD_E_BAD_NUMERIC);
    /* Not-a-number and infinity are refused rather than stored. */
    CHECK(gd_option_set(seed, &cfg, "nan", GD_OPT_MEASURE, &err) == GD_E_BAD_NUMERIC);
    CHECK(gd_option_set(seed, &cfg, "inf", GD_OPT_MEASURE, &err) == GD_E_BAD_NUMERIC);
}

TEST(test_options_write_the_right_field)
{
    /* The bug the registry replaced wrote the wrong number of bytes through a
     * void *, corrupting the field after the one it meant to set. Setting
     * each option must move its own field and leave every other one alone. */
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    gd_error err = GD_ERROR_INIT;

    for (size_t i = 0; i < n; i++) {
        gd_analysis_tree_config before, after;
        gd_analysis_tree_config_default(&before);
        gd_analysis_tree_config_default(&after);

        /* A value in range for this row. "1" suits most, but hpd_level is
         * (0, 1) exclusive, so fall back to something strictly inside. */
        int rc = gd_option_set(&opts[i], &after, "1", GD_OPT_MEASURE, &err);
        if (rc != GD_E_OK) {
            rc = gd_option_set(&opts[i], &after, "0.5", GD_OPT_MEASURE, &err);
        }
        CHECK(rc == GD_E_OK);

        /* Everything outside this option's field is byte-identical. */
        const unsigned char *b = (const unsigned char *)&before;
        const unsigned char *a = (const unsigned char *)&after;
        size_t width = 0;
        switch (opts[i].type) {
        case GD_OPT_SIZE:   width = sizeof(size_t);   break;
        case GD_OPT_U64:    width = sizeof(uint64_t); break;
        case GD_OPT_DOUBLE: width = sizeof(double);   break;
        case GD_OPT_BOOL:   width = sizeof(bool);     break;
        }
        for (size_t k = 0; k < sizeof before; k++) {
            if (k >= opts[i].offset && k < opts[i].offset + width) continue;
            CHECK(a[k] == b[k]);
        }
    }
}

TEST(test_options_bool_takes_only_zero_or_one)
{
    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    const gd_option *relaxed = gd_option_find(opts, n, "relaxed", GD_OPT_MEASURE);
    gd_error err = GD_ERROR_INIT;
    CHECK(relaxed->type == GD_OPT_BOOL);
    CHECK(gd_option_set(relaxed, &cfg, "1", GD_OPT_MEASURE, &err) == GD_E_OK);
    CHECK(cfg.relaxed == true);
    CHECK(gd_option_set(relaxed, &cfg, "0", GD_OPT_MEASURE, &err) == GD_E_OK);
    CHECK(cfg.relaxed == false);
    CHECK(gd_option_set(relaxed, &cfg, "2", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
}

TEST(test_options_cover_the_configurable_fields)
{
    /* A field added to the config without a row is reachable from nothing --
     * which is how the CLI ended up unable to set six options the shim could.
     * This is a count, so adding a field forces a decision about it rather
     * than letting it be forgotten. */
    size_t n;
    (void)gd_options_tree(&n);
    CHECK(n == 26);
    (void)gd_options_pair(&n);
    CHECK(n == 7);
}

TEST(test_options_reject_null_arguments)
{
    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    char buf[8];
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_option_set(NULL, &cfg, "1", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_set(&opts[0], NULL, "1", GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_set(&opts[0], &cfg, NULL, GD_OPT_MEASURE, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_format(NULL, &cfg, buf, sizeof buf, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_format(&opts[0], &cfg, buf, 0, &err) == GD_E_INVALID_ARG);
    CHECK(gd_option_find(NULL, 0, "x", GD_OPT_MEASURE) == NULL);
    /* A NULL error is always allowed. */
    CHECK(gd_option_set(&opts[0], &cfg, "zzz", GD_OPT_MEASURE, NULL) == GD_E_BAD_NUMERIC);
}

TEST(test_options_format_reports_a_short_buffer)
{
    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    size_t n;
    const gd_option *opts = gd_options_tree(&n);
    const gd_option *seed = gd_option_find(opts, n, "seed", GD_OPT_MEASURE);
    char tiny[2];
    gd_error err = GD_ERROR_INIT;
    /* "42" needs three bytes with its terminator. */
    CHECK(gd_option_format(seed, &cfg, tiny, sizeof tiny, &err) == GD_E_RANGE);
}

static void run_tests(void)
{
    RUN(test_options_every_row_round_trips);
    RUN(test_options_pair_rows_round_trip);
    RUN(test_options_names_and_flags_are_unique);
    RUN(test_options_offsets_land_inside_their_config);
    RUN(test_options_lookup_by_name_and_flag);
    RUN(test_options_reject_out_of_range);
    RUN(test_options_reject_non_numbers);
    RUN(test_options_write_the_right_field);
    RUN(test_options_bool_takes_only_zero_or_one);
    RUN(test_options_cover_the_configurable_fields);
    RUN(test_options_reject_null_arguments);
    RUN(test_options_format_reports_a_short_buffer);
}

TEST_MAIN()

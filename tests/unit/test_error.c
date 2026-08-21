/*
 * test_error.c — the error model's own contract.
 *
 * These tests exist because the convention documented at the top of
 * gd_error.h used to be documented and not true. Four different meanings of
 * an int return coexisted in the library, and a message that needed to name a
 * runtime value had nowhere to live unless the caller happened to pass a
 * scratch buffer. Both are things a comment cannot enforce, so they are
 * enforced here instead.
 */
#include "test_runner.h"
#include "gd.h"
#include "gd_tsv.h"

#include <string.h>

/* ---------------------------------------------------- the status convention */

TEST(test_set_returns_the_code_it_recorded)
{
    /* The whole point of the rule: what a function returns and what it puts in
     * the error are the same value, so a caller may read either. */
    gd_error err = GD_ERROR_INIT;
    int codes[] = { GD_E_NOMEM, GD_E_INVALID_ARG, GD_E_PARSE, GD_E_RANGE,
                    GD_E_CALIBRATION, GD_E_CANCELLED, GD_E_INTERNAL };
    size_t i;

    for (i = 0; i < sizeof codes / sizeof codes[0]; i++) {
        int returned = gd_error_set(&err, codes[i], "something went wrong");
        CHECK(returned == codes[i]);
        CHECK(err.code == codes[i]);
    }
}

TEST(test_setf_returns_the_code_it_recorded)
{
    gd_error err = GD_ERROR_INIT;
    int returned = gd_error_setf(&err, GD_E_BAD_NUMERIC, "bad value '%s'", "xyz");
    CHECK(returned == GD_E_BAD_NUMERIC);
    CHECK(err.code == GD_E_BAD_NUMERIC);
}

TEST(test_a_null_error_still_returns_the_code)
{
    /* Passing NULL drops the detail, never the classification — otherwise a
     * caller that does not want a message would be punished with a status it
     * cannot act on. */
    CHECK(gd_error_set(NULL, GD_E_PARSE, "ignored") == GD_E_PARSE);
    CHECK(gd_error_setf(NULL, GD_E_NOMEM, "ignored %d", 1) == GD_E_NOMEM);
}

TEST(test_ok_is_zero_so_a_plain_truth_test_works)
{
    /* Callers are allowed to write `if (f(...) != 0)` and be right. */
    CHECK(GD_E_OK == 0);
    CHECK(gd_error_set(NULL, GD_E_OK, NULL) == 0);
}

/* ------------------------------------------------------- owned message storage */

TEST(test_message_outlives_the_buffer_it_was_built_from)
{
    /* The regression the old design could not pass. The source string is a
     * local that goes out of scope; a borrowed pointer would dangle here. */
    gd_error err = GD_ERROR_INIT;

    {
        char scratch[64];
        snprintf(scratch, sizeof scratch, "concept '%s' is unusable", "hand");
        gd_error_set(&err, GD_E_INCOMPATIBLE, scratch);
        memset(scratch, 'X', sizeof scratch);   /* the buffer is reused */
    }

    CHECK_STREQ(err.message, "concept 'hand' is unusable");
}

TEST(test_setf_names_a_runtime_value_with_no_caller_buffer)
{
    /* This is what gd_analysis_diag existed to make possible. No caller
     * cooperation is required any more. */
    gd_error err = GD_ERROR_INIT;
    gd_error_setf(&err, GD_E_PARSE, "unknown calibration type '%s'", "not_a_type");
    CHECK(strstr(err.message, "not_a_type") != NULL);
}

TEST(test_message_is_never_null_and_is_empty_when_unset)
{
    /* Always safe to print, which is why no caller needs a NULL check. */
    gd_error err = GD_ERROR_INIT;
    CHECK(err.message[0] == '\0');
    CHECK_STREQ(gd_error_message(&err), "");
    CHECK_STREQ(gd_error_message(NULL), "");

    gd_error_set(&err, GD_E_PARSE, NULL);
    CHECK(err.message[0] == '\0');
}

TEST(test_an_overlong_message_truncates_rather_than_overflowing)
{
    gd_error err = GD_ERROR_INIT;
    char big[GD_ERROR_MSG_MAX * 2];

    memset(big, 'a', sizeof big - 1);
    big[sizeof big - 1] = '\0';
    gd_error_set(&err, GD_E_PARSE, big);

    CHECK(strlen(err.message) == GD_ERROR_MSG_MAX - 1);
    CHECK(err.code == GD_E_PARSE);

    gd_error_setf(&err, GD_E_PARSE, "%s", big);
    CHECK(strlen(err.message) == GD_ERROR_MSG_MAX - 1);
}

TEST(test_clear_resets_both_halves)
{
    gd_error err = GD_ERROR_INIT;
    gd_error_set(&err, GD_E_NOMEM, "out of memory");
    gd_error_clear(&err);
    CHECK(err.code == GD_E_OK);
    CHECK(err.message[0] == '\0');
    CHECK(!gd_error_failed(&err));
}

/* ------------------------------------------------------------- predicates */

TEST(test_predicates_answer_questions_not_statuses)
{
    /* gd_error_failed is a predicate: true means "there was a failure", the
     * opposite polarity to a status return, which is exactly why it must not
     * be an int that a reader might mistake for a code. */
    gd_error ok = GD_ERROR_INIT;
    gd_error bad = GD_ERROR_INIT;

    gd_error_set(&bad, GD_E_PARSE, "nope");
    CHECK(gd_error_failed(&bad));
    CHECK(!gd_error_failed(&ok));
    CHECK(!gd_error_failed(NULL));

    CHECK(gd_tsv_is_missing(".", 1));
    CHECK(!gd_tsv_is_missing("x", 1));
    CHECK(gd_tsv_field_eq("abc", 3, "abc"));
    CHECK(!gd_tsv_field_eq("abc", 3, "abd"));
}

/* ------------------------------------------ a real failure path end to end */

TEST(test_a_parser_reports_the_code_it_returns)
{
    /* The convention checked through a real function rather than through
     * gd_error_set alone: gd_tsv_parse's return value and err.code agree. */
    gd_error err = GD_ERROR_INIT;
    const char *no_header = "";
    int rc = gd_tsv_parse(no_header, 0, NULL, NULL, NULL, &err);

    if (rc != GD_E_OK) {
        CHECK(rc == err.code);
        CHECK(err.message[0] != '\0');
    }

    /* A malformed number is GD_E_BAD_NUMERIC, and not some other nonzero. */
    double v;
    CHECK(gd_tsv_parse_double("not a number", 12, &v) == GD_E_BAD_NUMERIC);
    /* A missing marker is distinguishable from a malformed value. */
    CHECK(gd_tsv_parse_double(".", 1, &v) == GD_E_MISSING_FIELD);
    /* And a good one succeeds. */
    CHECK(gd_tsv_parse_double("2.5", 3, &v) == GD_E_OK);
    CHECK_APPROX(v, 2.5, 1e-12);
}

static void run_tests(void)
{
    RUN(test_set_returns_the_code_it_recorded);
    RUN(test_setf_returns_the_code_it_recorded);
    RUN(test_a_null_error_still_returns_the_code);
    RUN(test_ok_is_zero_so_a_plain_truth_test_works);
    RUN(test_message_outlives_the_buffer_it_was_built_from);
    RUN(test_setf_names_a_runtime_value_with_no_caller_buffer);
    RUN(test_message_is_never_null_and_is_empty_when_unset);
    RUN(test_an_overlong_message_truncates_rather_than_overflowing);
    RUN(test_clear_resets_both_halves);
    RUN(test_predicates_answer_questions_not_statuses);
    RUN(test_a_parser_reports_the_code_it_returns);
}

TEST_MAIN()

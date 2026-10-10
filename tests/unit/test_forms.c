/*
 * test_forms.c — one reading of the forms table.
 *
 * The rule under test is "highest primary_weight wins, first row wins a tie",
 * plus the two orderings the rest of the engine depends on: concepts and
 * languages in first-seen row order.
 *
 * Every case snapshots what it needs and destroys the gd_forms BEFORE
 * asserting. CHECK returns from the test on failure, so a live allocation
 * held across an assertion leaks, and LeakSanitizer's report then buries the
 * FAIL line that says which assertion fired. That is not hypothetical: it is
 * how the first draft of this file behaved when the rule was broken on
 * purpose to confirm these tests can fail.
 */
#include "test_runner.h"
#include "gd_forms.h"

#include <string.h>

static const char *FORMS =
    "language_id\tconcept_id\tcognate_id\tprimary_weight\n"
    "latin\thand\tH1\t1.0\n"
    "latin\tfoot\tF1\t1.0\n"
    "spanish\thand\tH1\t0.4\n"
    "spanish\thand\tH9\t0.9\n"      /* higher weight wins */
    "spanish\tfoot\tF2\t1.0\n"
    "ghost\twater\tW1\t1.0\n"       /* a language no tree will mention */
    "ghost\twater\tW2\t1.0\n";      /* equal weight: the first row wins */

/* Snapshot of a cell or a name, so the gd_forms can die before we assert. */
typedef struct { char s[64]; int present; } snap;

static void take(snap *o, const char *s)
{
    if (s == NULL) { o->present = 0; o->s[0] = '\0'; return; }
    o->present = 1;
    strncpy(o->s, s, sizeof o->s - 1);
    o->s[sizeof o->s - 1] = '\0';
}

TEST(test_forms_first_seen_order)
{
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);

    size_t nconcept = 0, nlang = 0;
    snap c0, c1, c2, l0, l1, l2;
    memset(&c0, 0, sizeof c0); memset(&c1, 0, sizeof c1); memset(&c2, 0, sizeof c2);
    memset(&l0, 0, sizeof l0); memset(&l1, 0, sizeof l1); memset(&l2, 0, sizeof l2);
    if (rc == GD_E_OK) {
        nconcept = f->nconcept; nlang = f->nlang;
        if (nconcept > 2) { take(&c0, f->concepts[0]); take(&c1, f->concepts[1]); take(&c2, f->concepts[2]); }
        if (nlang > 2) { take(&l0, f->languages[0]); take(&l1, f->languages[1]); take(&l2, f->languages[2]); }
        gd_forms_destroy(f);
    }

    CHECK(rc == GD_E_OK);
    /* Concepts in the order the rows present them, over every row. */
    CHECK(nconcept == 3);
    CHECK_STREQ(c0.s, "hand");
    CHECK_STREQ(c1.s, "foot");
    CHECK_STREQ(c2.s, "water");
    /* Languages likewise, including one no tree would match. */
    CHECK(nlang == 3);
    CHECK_STREQ(l0.s, "latin");
    CHECK_STREQ(l1.s, "spanish");
    CHECK_STREQ(l2.s, "ghost");
}

TEST(test_forms_highest_weight_wins)
{
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);
    snap got; memset(&got, 0, sizeof got);
    size_t sp = SIZE_MAX;
    if (rc == GD_E_OK) {
        sp = gd_forms_language_index(f, "spanish");
        if (sp != SIZE_MAX) take(&got, gd_forms_primary(f, 0, sp));
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    CHECK(sp != SIZE_MAX);
    /* H9 at 0.9 beats H1 at 0.4, even though H1 came first. */
    CHECK_STREQ(got.s, "H9");
}

TEST(test_forms_first_row_wins_a_tie)
{
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);
    snap got; memset(&got, 0, sizeof got);
    if (rc == GD_E_OK) {
        size_t gh = gd_forms_language_index(f, "ghost");
        if (gh != SIZE_MAX) take(&got, gd_forms_primary(f, 2 /* water */, gh));
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    CHECK_STREQ(got.s, "W1");
}

TEST(test_forms_absent_cell_is_null)
{
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);
    int latin_water_present = 1, etruscan_known = 1;
    if (rc == GD_E_OK) {
        size_t latin = gd_forms_language_index(f, "latin");
        latin_water_present = (gd_forms_primary(f, 2 /* water */, latin) != NULL);
        etruscan_known = (gd_forms_language_index(f, "etruscan") != SIZE_MAX);
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    /* latin has no form for 'water'. */
    CHECK(!latin_water_present);
    /* And a language the file never mentions has no index at all. */
    CHECK(!etruscan_known);
}

TEST(test_forms_rows_per_language_counts_rows)
{
    /* The chronology intake reports skipped ROWS, not skipped languages, so
     * this count has to be per row. */
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);
    size_t rl = 0, rs = 0, rg = 0;
    if (rc == GD_E_OK) {
        rl = f->rows_per_lang[gd_forms_language_index(f, "latin")];
        rs = f->rows_per_lang[gd_forms_language_index(f, "spanish")];
        rg = f->rows_per_lang[gd_forms_language_index(f, "ghost")];
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    CHECK(rl == 2);
    CHECK(rs == 3);
    CHECK(rg == 2);
}

TEST(test_forms_assignments_are_in_file_order_and_include_overrides)
{
    /* A caller numbering cognate classes itself replays this sequence, so it
     * must contain the row that was later overridden (spanish/hand H1) and
     * must be in file order. */
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);
    size_t na = 0;
    snap a0, a2, a3, a5;
    int same_cell = 0;
    memset(&a0, 0, sizeof a0); memset(&a2, 0, sizeof a2);
    memset(&a3, 0, sizeof a3); memset(&a5, 0, sizeof a5);
    if (rc == GD_E_OK) {
        na = f->nassignments;
        if (na >= 6) {
            take(&a0, f->assignments[0].cognate);
            take(&a2, f->assignments[2].cognate);
            take(&a3, f->assignments[3].cognate);
            take(&a5, f->assignments[5].cognate);
            same_cell = (f->assignments[2].concept == f->assignments[3].concept)
                     && (f->assignments[2].lang    == f->assignments[3].lang);
        }
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    /* latin/hand, latin/foot, spanish/hand(H1), spanish/hand(H9),
     * spanish/foot, ghost/water(W1). W2 ties and never takes the cell. */
    CHECK(na == 6);
    CHECK_STREQ(a0.s, "H1");
    CHECK_STREQ(a2.s, "H1");   /* the override victim, still in the sequence */
    CHECK_STREQ(a3.s, "H9");
    CHECK(same_cell);
    CHECK_STREQ(a5.s, "W1");
}

TEST(test_forms_unreadable_weight_falls_back_and_is_counted)
{
    /* Settled: an unreadable weight is taken as 1.0 everywhere, and counted.
     * The pairwise pipeline used to refuse the whole run instead; see
     * gd_forms.h on why that changed and what it gives up. */
    static const char *bad =
        "language_id\tconcept_id\tcognate_id\tprimary_weight\n"
        "latin\thand\tH1\tnot-a-number\n"
        "latin\thand\tH2\t0.5\n";

    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(bad, strlen(bad), &f, &err);
    snap got; memset(&got, 0, sizeof got);
    size_t nbad = 0;
    if (rc == GD_E_OK) {
        take(&got, gd_forms_primary(f, 0, 0));
        nbad = f->n_unreadable_weights;
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    /* Taken as 1.0, so it keeps the cell against a readable 0.5 -- it is not
     * dropped, which is what reading it as zero would have done. */
    CHECK_STREQ(got.s, "H1");
    /* And the file's one bad cell is still reportable. */
    CHECK(nbad == 1);
}

TEST(test_forms_clean_file_counts_no_unreadable_weights)
{
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(FORMS, strlen(FORMS), &f, &err);
    size_t nbad = 1;
    if (rc == GD_E_OK) { nbad = f->n_unreadable_weights; gd_forms_destroy(f); }
    CHECK(rc == GD_E_OK);
    CHECK(nbad == 0);
}

TEST(test_forms_missing_weight_column_defaults_to_one)
{
    static const char *nw =
        "language_id\tconcept_id\tcognate_id\n"
        "latin\thand\tH1\n"
        "latin\thand\tH2\n";
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(nw, strlen(nw), &f, &err);
    snap got; memset(&got, 0, sizeof got);
    if (rc == GD_E_OK) { take(&got, gd_forms_primary(f, 0, 0)); gd_forms_destroy(f); }
    CHECK(rc == GD_E_OK);
    /* All weights equal, so the first row keeps the cell. */
    CHECK_STREQ(got.s, "H1");
}

TEST(test_forms_missing_required_column_is_refused)
{
    static const char *nolang =
        "concept_id\tcognate_id\n"
        "hand\tH1\n";
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(nolang, strlen(nolang), &f, &err);
    gd_forms_destroy(f);
    CHECK(rc == GD_E_UNKNOWN_COLUMN);
    CHECK(f == NULL);
}

TEST(test_forms_rejects_null_and_reports_empty)
{
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_forms_read(NULL, 0, &f, &err) == GD_E_INVALID_ARG);

    /* A header with no data rows is a real, empty table — not an error here.
     * The pipelines decide whether empty is acceptable. */
    static const char *hdr = "language_id\tconcept_id\tcognate_id\n";
    gd_error err2 = GD_ERROR_INIT;
    int rc = gd_forms_read(hdr, strlen(hdr), &f, &err2);
    size_t nc = 1, nl = 1, na = 1;
    if (rc == GD_E_OK) { nc = f->nconcept; nl = f->nlang; na = f->nassignments; gd_forms_destroy(f); }
    CHECK(rc == GD_E_OK);
    CHECK(nc == 0);
    CHECK(nl == 0);
    CHECK(na == 0);
}

TEST(test_forms_field_with_embedded_nul_does_not_overread)
{
    /* Regression: the first version of field_is() compared with strncmp and
     * then read s[len]. A field carrying an embedded NUL makes strncmp stop
     * early and report equality, after which s[len] runs past the end of a
     * shorter stored string. libFuzzer found it through fuzz_analysis; this
     * is the same shape as a plain unit test.
     *
     * The reading itself truncates at the NUL, which is what the string-table
     * loader gd_forms replaced did via gd_strdup_n plus strcmp. */
    /* Order matters: "latin" must be interned FIRST, so the later field is
     * longer than the stored string it is compared against. That is the
     * direction in which the old s[len] read ran off the end. */
    static const char raw[] =
        "language_id\tconcept_id\tcognate_id\n"
        "latin\thand\tH1\n"
        "latin\0extra\thand\tH2\n";
    size_t len = sizeof raw - 1;      /* keeps the embedded NUL */

    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(raw, len, &f, &err);
    size_t nlang = 0;
    snap l0, got;
    memset(&l0, 0, sizeof l0); memset(&got, 0, sizeof got);
    if (rc == GD_E_OK) {
        nlang = f->nlang;
        if (nlang > 0) take(&l0, f->languages[0]);
        take(&got, gd_forms_primary(f, 0, 0));
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    /* "latin\0extra" and "latin" are the same language: both truncate. */
    CHECK(nlang == 1);
    CHECK_STREQ(l0.s, "latin");
    /* Equal weights, so the first row keeps the cell. */
    CHECK_STREQ(got.s, "H1");
}

/* Codex finding 13: weights that are numbers but not weights, and rows with
 * no language or concept, are refused with their line number. */
TEST(test_forms_rejects_non_weights_and_empty_ids)
{
    static const char *const BAD[] = {
        "language_id\tconcept_id\tcognate_id\tprimary_weight\nA\tc\tk\tnan\n",
        "language_id\tconcept_id\tcognate_id\tprimary_weight\nA\tc\tk\tinf\n",
        "language_id\tconcept_id\tcognate_id\tprimary_weight\nA\tc\tk\t-0.5\n",
        "language_id\tconcept_id\tcognate_id\nA\tc\tk\n\tc\tk\n",
        "language_id\tconcept_id\tcognate_id\nA\tc\tk\nA\n",
    };
    for (size_t i = 0; i < sizeof BAD / sizeof BAD[0]; i++) {
        gd_forms *f = NULL;
        gd_error err = GD_ERROR_INIT;
        int rc = gd_forms_read(BAD[i], strlen(BAD[i]), &f, &err);
        gd_forms_destroy(f);
        CHECK(rc != GD_E_OK);
        CHECK(strstr(err.message, "line") != NULL);
    }
}

/* E0.1 at the reader: unobserved rows are counted and give no cell. */
TEST(test_forms_skips_unobserved_rows)
{
    static const char *F =
        "language_id\tconcept_id\tcognate_id\tstatus\n"
        "A\tc1\tk1\tobserved\n"
        "A\tc2\t.\tobserved\n"
        "A\tc3\tk3\tabsent\n";
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(F, strlen(F), &f, &err);
    size_t miss = 0, nconcept = 0, rows = 0;
    int c2_null = 0;
    if (rc == GD_E_OK) {
        miss = f->n_missing_cognates; nconcept = f->nconcept;
        rows = f->rows_per_lang[0];
        c2_null = (gd_forms_primary(f, 1, 0) == NULL);
        gd_forms_destroy(f);
    }
    CHECK(rc == GD_E_OK);
    CHECK(miss == 2);
    CHECK(nconcept == 3);   /* concepts keep the shape the file declares */
    CHECK(rows == 1);
    CHECK(c2_null);
}

/* E2.4: two codings in one file must not be mixed; a selected one is read
 * alone; loans become missing. */
TEST(test_forms_analysis_and_loans)
{
    static const char *F =
        "analysis_id\tlanguage_id\tconcept_id\tcognate_id\tloan_probability\n"
        "strict\tA\tc1\tk1\t0.0\n"
        "strict\tB\tc1\tk1\t0.9\n"
        "loose\tA\tc1\tk7\t.\n"
        "loose\tB\tc1\tk7\t.\n";
    gd_forms *f = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_forms_read(F, strlen(F), &f, &err);
    gd_forms_destroy(f); f = NULL;
    CHECK(rc != GD_E_OK);
    CHECK(strstr(err.message, "strict") != NULL && strstr(err.message, "loose") != NULL);

    gd_forms_options o; gd_forms_options_default(&o);
    o.analysis_id = "strict";
    gd_error_clear(&err);
    rc = gd_forms_read_opts(F, strlen(F), &o, &f, &err);
    size_t loans = 0, other = 0; int b_missing = 0, a_k1 = 0;
    if (rc == GD_E_OK) {
        loans = f->n_loans; other = f->n_other_analysis;
        size_t la = gd_forms_language_index(f, "A"), lb = gd_forms_language_index(f, "B");
        a_k1 = (la != SIZE_MAX && gd_forms_primary(f, 0, la) != NULL &&
                strcmp(gd_forms_primary(f, 0, la), "k1") == 0);
        b_missing = (lb != SIZE_MAX && gd_forms_primary(f, 0, lb) == NULL);
        gd_forms_destroy(f); f = NULL;
    }
    CHECK(rc == GD_E_OK);
    CHECK(loans == 1);
    CHECK(other == 2);
    CHECK(a_k1);
    CHECK(b_missing);

    /* A threshold above 1 keeps the loan. */
    o.loan_threshold = 1.5;
    rc = gd_forms_read_opts(F, strlen(F), &o, &f, &err);
    loans = (rc == GD_E_OK) ? f->n_loans : 99;
    gd_forms_destroy(f);
    CHECK(loans == 0);
}

static void run_tests(void)
{
    RUN(test_forms_analysis_and_loans);
    RUN(test_forms_rejects_non_weights_and_empty_ids);
    RUN(test_forms_skips_unobserved_rows);
    RUN(test_forms_first_seen_order);
    RUN(test_forms_highest_weight_wins);
    RUN(test_forms_first_row_wins_a_tie);
    RUN(test_forms_absent_cell_is_null);
    RUN(test_forms_rows_per_language_counts_rows);
    RUN(test_forms_assignments_are_in_file_order_and_include_overrides);
    RUN(test_forms_unreadable_weight_falls_back_and_is_counted);
    RUN(test_forms_clean_file_counts_no_unreadable_weights);
    RUN(test_forms_missing_weight_column_defaults_to_one);
    RUN(test_forms_missing_required_column_is_refused);
    RUN(test_forms_rejects_null_and_reports_empty);
    RUN(test_forms_field_with_embedded_nul_does_not_overread);
}

TEST_MAIN()

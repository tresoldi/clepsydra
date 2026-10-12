/*
 * test_report.c — the reported tables and their two renderers.
 *
 * Before D5 none of this was reachable from a test: every column name and
 * every format lived inside a command body in cli/commands.c, which has no C
 * test at all. The formats here are the ones the committed example outputs
 * were produced with, so a moved decimal is a change in the evidence.
 */
#include "test_runner.h"
#include "gd_report.h"

#include <stdlib.h>
#include <string.h>

/* The test's adapter at the sink seam: bytes into a fixed buffer. */
typedef struct {
    char   buf[4096];
    size_t len;
    int    fail_after;   /* -1 never; otherwise refuse once len exceeds it */
} membuf;

static int mem_write(void *ctx, const char *b, size_t n)
{
    membuf *m = (membuf *)ctx;
    if (m->fail_after >= 0 && m->len + n > (size_t)m->fail_after) return GD_E_IO;
    if (m->len + n >= sizeof m->buf) return GD_E_IO;
    memcpy(m->buf + m->len, b, n);
    m->len += n;
    m->buf[m->len] = '\0';
    return GD_E_OK;
}

static void mem_init(membuf *m) { m->len = 0; m->buf[0] = '\0'; m->fail_after = -1; }

static const gd_column COLS[] = {
    { "name",  "name",  0, GD_NUM_FIXED },
    { "fixed", "fixed", 2, GD_NUM_FIXED },
    { "gen",   "gen",   4, GD_NUM_GENERAL },
    { NULL,    "only_json", 0, GD_NUM_FIXED },
    { "only_tsv", NULL,  0, GD_NUM_FIXED },
};

static gd_table *sample(gd_error *err)
{
    gd_table *t = NULL;
    if (gd_table_create(COLS, 5, &t, err) != GD_E_OK) return NULL;
    if (gd_table_add_row(t, err) != GD_E_OK) { gd_table_destroy(t); return NULL; }
    gd_table_set_text(t, 0, "alpha", err);
    gd_table_set_num(t, 1, 2670.0, err);
    gd_table_set_num(t, 2, -9.002166, err);
    gd_table_set_uint(t, 3, 7, err);
    gd_table_set_text(t, 4, "tsv-only", err);
    return t;
}

TEST(test_report_tsv_uses_the_column_formats)
{
    gd_error err = GD_ERROR_INIT;
    gd_table *t = sample(&err);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    int rc = (t != NULL) ? gd_table_write_tsv(t, &sink, &err) : GD_E_INTERNAL;
    gd_table_destroy(t);

    CHECK(rc == GD_E_OK);
    /* json_name-only columns are absent; fixed keeps decimals, general keeps
     * significant digits. */
    CHECK_STREQ(m.buf, "name\tfixed\tgen\tonly_tsv\n"
                       "alpha\t2670.00\t-9.002\ttsv-only\n");
}

TEST(test_report_json_uses_the_other_names)
{
    gd_error err = GD_ERROR_INIT;
    gd_table *t = sample(&err);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    int rc = (t != NULL) ? gd_table_write_json(t, &sink, &err) : GD_E_INTERNAL;
    gd_table_destroy(t);

    CHECK(rc == GD_E_OK);
    /* tsv_name-only columns are absent, and numbers keep %.10g regardless of
     * the TSV's precision -- the two schemas are deliberately different. */
    CHECK_STREQ(m.buf,
        "[{\"name\":\"alpha\",\"fixed\":2670,\"gen\":-9.002166,\"only_json\":7}]");
}

TEST(test_report_json_escapes_and_nulls)
{
    static const gd_column c[] = { { "t", "t", 0, GD_NUM_FIXED },
                                   { "n", "n", 1, GD_NUM_FIXED } };
    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_table_create(c, 2, &t, &err) == GD_E_OK);
    CHECK(gd_table_add_row(t, &err) == GD_E_OK);
    gd_table_set_text(t, 0, "a\"b\\c\nd", &err);
    gd_table_set_num(t, 1, HUGE_VAL, &err);      /* not finite */

    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    int rc = gd_table_write_json(t, &sink, &err);
    gd_table_destroy(t);

    CHECK(rc == GD_E_OK);
    /* A non-finite number becomes null so JSON.parse cannot throw on a
     * degenerate run -- the shim's emitter has always done this. */
    CHECK_STREQ(m.buf, "[{\"t\":\"a\\\"b\\\\c\\nd\",\"n\":null}]");
}

TEST(test_report_cells_read_back_without_parsing)
{
    /* The table is the test surface: a case asserts on values, not on text it
     * had to parse back out of a rendering. */
    gd_error err = GD_ERROR_INIT;
    gd_table *t = sample(&err);
    CHECK(t != NULL);
    CHECK(gd_table_nrows(t) == 1);
    CHECK(gd_table_ncols(t) == 5);
    CHECK_STREQ(gd_table_text_at(t, 0, 0), "alpha");
    CHECK_APPROX(gd_table_num_at(t, 0, 1), 2670.0, 1e-12);
    CHECK_APPROX(gd_table_num_at(t, 0, 3), 7.0, 1e-12);
    /* Out of range reads are NULL / NaN rather than a crash. */
    CHECK(gd_table_text_at(t, 5, 0) == NULL);
    CHECK(gd_table_text_at(t, 0, 9) == NULL);
    gd_table_destroy(t);
}

TEST(test_report_a_failing_sink_stops_the_render)
{
    /* A sink that refuses must surface its code, not be ignored: on the
     * command line that sink is a file, and a silent short write would leave
     * a truncated table looking complete. */
    gd_error err = GD_ERROR_INIT;
    gd_table *t = sample(&err);
    membuf m; mem_init(&m);
    m.fail_after = 4;                 /* accepts the first few bytes only */
    gd_sink sink = { &m, mem_write };
    int rc = (t != NULL) ? gd_table_write_tsv(t, &sink, &err) : GD_E_INTERNAL;
    gd_table_destroy(t);
    CHECK(rc == GD_E_IO);
}

TEST(test_report_empty_table_still_writes_its_header)
{
    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_table_create(COLS, 5, &t, &err) == GD_E_OK);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    CHECK(gd_table_write_tsv(t, &sink, &err) == GD_E_OK);
    CHECK_STREQ(m.buf, "name\tfixed\tgen\tonly_tsv\n");
    mem_init(&m);
    CHECK(gd_table_write_json(t, &sink, &err) == GD_E_OK);
    CHECK_STREQ(m.buf, "[]");
    gd_table_destroy(t);
}

TEST(test_report_posterior_nodes_columns_are_the_shipped_ones)
{
    /* The header of posterior_nodes.tsv is a published schema:
     * tools/validation/calibration_benchmark.py reads it back by name. If
     * this changes, that script breaks silently. */
    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_report_nodes(&res, &t, &err) == GD_E_OK);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    CHECK(gd_table_write_tsv(t, &sink, &err) == GD_E_OK);
    gd_table_destroy(t);
    CHECK_STREQ(m.buf,
        "node\tmean_bp\tmedian_bp\tsd\thpd90_lo_bp\thpd90_hi_bp\trhat\tess_bulk\tess_tail\tlexical_median\tlexical_hpd90_lo\tlexical_hpd90_hi\trelative_median\trelative_hpd90_lo\trelative_hpd90_hi\tprior_median_bp\tprior_hpd90_lo_bp\tprior_hpd90_hi_bp\tprior_overlap\tcalendar_basis\n");
}

TEST(test_report_diagnostics_columns_are_the_shipped_ones)
{
    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_report_diagnostics(&res, &t, &err) == GD_E_OK);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    CHECK(gd_table_write_tsv(t, &sink, &err) == GD_E_OK);
    gd_table_destroy(t);
    CHECK_STREQ(m.buf, "param\ttype\tmean\tsd\trhat\tess_bulk\tess_tail\tmcse_mean\tacceptance\n");
}

TEST(test_report_pair_is_a_key_value_table)
{
    gd_analysis_pair_result res;
    memset(&res, 0, sizeof res);
    res.model = "p1";
    res.n_shared = 40;
    res.n_cognate = 13;
    res.retention = 0.325;
    res.median = 2656.4;
    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_report_pair(&res, "Latin", "Spanish", &t, &err) == GD_E_OK);

    /* Every row is key then value, and the numbers carry the precision the
     * shipped file has always used. */
    CHECK(gd_table_nrows(t) == 19);
    CHECK_STREQ(gd_table_text_at(t, 0, 0), "model");
    CHECK_STREQ(gd_table_text_at(t, 0, 1), "p1");
    CHECK_STREQ(gd_table_text_at(t, 1, 1), "Latin");
    CHECK_STREQ(gd_table_text_at(t, 2, 1), "Spanish");
    CHECK_STREQ(gd_table_text_at(t, 5, 0), "retention_fraction");
    CHECK_STREQ(gd_table_text_at(t, 5, 1), "0.3250");
    CHECK_STREQ(gd_table_text_at(t, 7, 0), "divergence_median_bp");
    CHECK_STREQ(gd_table_text_at(t, 7, 1), "2656.4");
    gd_table_destroy(t);
}

TEST(test_report_model_name_matches_the_relaxed_flag)
{
    /* This string used to be two independent literals, in cli/commands.c and
     * web/clepsydra_wasm.c, with nothing comparing them. */
    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    CHECK_STREQ(gd_report_tree_model_name(&res), "T1_novel_state");
    res.relaxed = true;
    CHECK_STREQ(gd_report_tree_model_name(&res), "T1_novel_state_relaxed");
}

TEST(test_report_rejects_bad_arguments)
{
    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_table_create(NULL, 2, &t, &err) == GD_E_INVALID_ARG);
    CHECK(gd_table_create(COLS, 0, &t, &err) == GD_E_INVALID_ARG);
    CHECK(gd_report_nodes(NULL, &t, &err) == GD_E_INVALID_ARG);
    CHECK(gd_report_diagnostics(NULL, &t, &err) == GD_E_INVALID_ARG);
    CHECK(gd_report_pair(NULL, "a", "b", &t, &err) == GD_E_INVALID_ARG);
    /* Setting a cell before any row exists is refused, not written past. */
    CHECK(gd_table_create(COLS, 5, &t, &err) == GD_E_OK);
    CHECK(gd_table_set_text(t, 0, "x", &err) == GD_E_INVALID_ARG);
    CHECK(gd_table_set_num(t, 99, 1.0, &err) == GD_E_INVALID_ARG);
    gd_table_destroy(t);
    /* Destroying NULL is allowed. */
    gd_table_destroy(NULL);
}

TEST(test_report_diagnostics_row_is_byte_for_byte_the_shipped_format)
{
    /*
     * The header alone does not pin the file: the number formats do too, and
     * they are not uniform. mean and sd are %.4g -- significant digits, not
     * decimals -- because they summarise parameters whose magnitudes differ by
     * orders: a log-rate near -9 and an age near 2670 share this column.
     * Rendering them as %.4f gives "-9.0022" and "2670.0000", which is what a
     * first draft of gd_report did; only the example diff caught it. This
     * case catches it in the unit suite.
     */
    gd_analysis_param params[2];
    memset(params, 0, sizeof params);
    params[0].type = GD_MCMC_P_LOGRATE;
    params[0].type_name = "log_rate";
    params[0].label = NULL;
    params[0].stat.mean = -9.002166;
    params[0].stat.sd = 0.1188447;
    params[0].rhat = 1.0009;
    params[0].ess = 834.4;
    params[0].ess_tail = 701.2;
    params[0].mcse = 0.004114;
    params[0].acceptance = 0.4441;

    params[1].type = GD_MCMC_P_AGE;
    params[1].type_name = "age";
    params[1].label = "root";
    params[1].stat.mean = 2670.31;
    params[1].stat.sd = 104.72;
    params[1].rhat = 1.0011;
    params[1].ess = 2295.4;
    params[1].ess_tail = 1999.6;
    params[1].mcse = 2.18571;
    params[1].acceptance = 0.4432;

    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    res.params = params;
    res.nparams = 2;

    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_report_diagnostics(&res, &t, &err) == GD_E_OK);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    CHECK(gd_table_write_tsv(t, &sink, &err) == GD_E_OK);
    gd_table_destroy(t);

    CHECK_STREQ(m.buf,
        "param\ttype\tmean\tsd\trhat\tess_bulk\tess_tail\tmcse_mean\tacceptance\n"
        "log_rate\tlog_rate\t-9.002\t0.1188\t1.001\t834\t701\t0.004114\t0.444\n"
        "age:root\tage\t2670\t104.7\t1.001\t2295\t2000\t2.186\t0.443\n");
}

TEST(test_report_node_row_is_byte_for_byte_the_shipped_format)
{
    gd_analysis_node_age ages[1];
    memset(ages, 0, sizeof ages);
    ages[0].node = 4;
    ages[0].label = "proto_romance";
    ages[0].age.mean = 1610.44;
    ages[0].age.median = 1601.27;
    ages[0].age.sd = 143.81;
    ages[0].age.hpd_lo = 1370.06;
    ages[0].age.hpd_hi = 1849.92;
    ages[0].rhat = 1.00449;
    ages[0].ess = 740.6;
    ages[0].ess_tail = 512.3;
    ages[0].lexical.median = 0.24133;
    ages[0].lexical.hpd_lo = 0.19021;
    ages[0].lexical.hpd_hi = 0.30477;
    ages[0].relative.median = 0.5871;
    ages[0].relative.hpd_lo = 0.5102;
    ages[0].relative.hpd_hi = 0.6648;
    ages[0].prior.median = 900.04;
    ages[0].prior.hpd_lo = 120.55;
    ages[0].prior.hpd_hi = 2500.97;
    ages[0].prior_overlap = 0.2149;

    gd_analysis_tree_result res;
    memset(&res, 0, sizeof res);
    res.ages = ages;
    res.nages = 1;

    gd_error err = GD_ERROR_INIT;
    gd_table *t = NULL;
    CHECK(gd_report_nodes(&res, &t, &err) == GD_E_OK);
    membuf m; mem_init(&m);
    gd_sink sink = { &m, mem_write };
    CHECK(gd_table_write_tsv(t, &sink, &err) == GD_E_OK);

    CHECK_STREQ(m.buf,
        "node\tmean_bp\tmedian_bp\tsd\thpd90_lo_bp\thpd90_hi_bp\trhat\tess_bulk\tess_tail\tlexical_median\tlexical_hpd90_lo\tlexical_hpd90_hi\trelative_median\trelative_hpd90_lo\trelative_hpd90_hi\tprior_median_bp\tprior_hpd90_lo_bp\tprior_hpd90_hi_bp\tprior_overlap\tcalendar_basis\n"
        "proto_romance\t1610.4\t1601.3\t143.8\t1370.1\t1849.9\t1.004\t741\t512"
        "\t0.2413\t0.1902\t0.3048\t0.587\t0.510\t0.665\t900.0\t120.5\t2501.0\t0.215\tanchored\n");

    /* The same row as the browser sees it: different key names, different
     * precision, and the node index the TSV has no column for. */
    mem_init(&m);
    CHECK(gd_table_write_json(t, &sink, &err) == GD_E_OK);
    gd_table_destroy(t);
    CHECK_STREQ(m.buf,
        "[{\"label\":\"proto_romance\",\"node\":4,\"mean\":1610.44,"
        "\"median\":1601.27,\"sd\":143.81,\"hpd_lo\":1370.06,"
        "\"hpd_hi\":1849.92,\"rhat\":1.00449,\"ess\":740.6,\"ess_tail\":512.3,"
        "\"lexical_median\":0.24133,\"lexical_hpd_lo\":0.19021,\"lexical_hpd_hi\":0.30477,"
        "\"relative_median\":0.5871,\"relative_hpd_lo\":0.5102,\"relative_hpd_hi\":0.6648,"
        "\"prior_median\":900.04,\"prior_hpd_lo\":120.55,\"prior_hpd_hi\":2500.97,"
        "\"prior_overlap\":0.2149,\"calendar_basis\":\"anchored\"}]");
}

static void run_tests(void)
{
    RUN(test_report_tsv_uses_the_column_formats);
    RUN(test_report_json_uses_the_other_names);
    RUN(test_report_json_escapes_and_nulls);
    RUN(test_report_cells_read_back_without_parsing);
    RUN(test_report_a_failing_sink_stops_the_render);
    RUN(test_report_empty_table_still_writes_its_header);
    RUN(test_report_posterior_nodes_columns_are_the_shipped_ones);
    RUN(test_report_diagnostics_columns_are_the_shipped_ones);
    RUN(test_report_diagnostics_row_is_byte_for_byte_the_shipped_format);
    RUN(test_report_node_row_is_byte_for_byte_the_shipped_format);
    RUN(test_report_pair_is_a_key_value_table);
    RUN(test_report_model_name_matches_the_relaxed_flag);
    RUN(test_report_rejects_bad_arguments);
}

TEST_MAIN()

/*
 * test_tsv.c — tests for the strict TSV parser.
 */
#include "test_runner.h"
#include "gd_tsv.h"

#include <string.h>

/* A simple collector for testing. */
typedef struct {
    int header_rows;
    int data_rows;
    int last_nfields;
    char last_field0[128];
    int last_line;
} collector;

static int header_cb(const gd_tsv_row *row, void *user, gd_error *err)
{
    (void)err;
    collector *c = (collector *)user;
    c->header_rows++;
    c->last_nfields = (int)row->nfields;
    return 0;
}

static int row_cb(const gd_tsv_row *row, void *user, gd_error *err)
{
    (void)err;
    collector *c = (collector *)user;
    c->data_rows++;
    c->last_nfields = (int)row->nfields;
    c->last_line = (int)row->line_no;
    if (row->nfields > 0 && row->lengths[0] < sizeof(c->last_field0)) {
        memcpy(c->last_field0, row->fields[0], row->lengths[0]);
        c->last_field0[row->lengths[0]] = '\0';
    }
    return 0;
}

TEST(test_tsv_basic_parse)
{
    const char *data = "a\tb\tc\n1\t2\t3\n4\t5\t6\n";
    gd_error err = GD_ERROR_INIT;
    collector c = {0};
    int rc = gd_tsv_parse(data, strlen(data), header_cb, row_cb, &c, &err);
    CHECK(rc == 0);
    CHECK(c.header_rows == 1);
    CHECK(c.data_rows == 2);
    CHECK(c.last_nfields == 3);
}

TEST(test_tsv_comments_and_blanks)
{
    const char *data = "# comment\n\na\tb\n# another\n\n1\t2\n";
    gd_error err = GD_ERROR_INIT;
    collector c = {0};
    int rc = gd_tsv_parse(data, strlen(data), header_cb, row_cb, &c, &err);
    CHECK(rc == 0);
    CHECK(c.header_rows == 1);
    CHECK(c.data_rows == 1);
}

TEST(test_tsv_missing_header)
{
    const char *data = "\n\n\n";  /* only blank lines */
    gd_error err = GD_ERROR_INIT;
    collector c = {0};
    int rc = gd_tsv_parse(data, strlen(data), header_cb, row_cb, &c, &err);
    CHECK(rc == GD_E_MISSING_HEADER);
}

TEST(test_tsv_missing_dot)
{
    CHECK(gd_tsv_is_missing(".", 1) == 1);
    CHECK(gd_tsv_is_missing("x", 1) == 0);
    CHECK(gd_tsv_is_missing("..", 2) == 0);
}

TEST(test_tsv_field_eq)
{
    CHECK(gd_tsv_field_eq("hello", 5, "hello") == 1);
    CHECK(gd_tsv_field_eq("hello", 5, "world") == 0);
    CHECK(gd_tsv_field_eq("hell", 4, "hello") == 0);
}

TEST(test_tsv_parse_double)
{
    double out;
    CHECK(gd_tsv_parse_double("3.14", 4, &out) == 0);
    CHECK_APPROX(out, 3.14, 1e-10);
    CHECK(gd_tsv_parse_double("-0.5", 4, &out) == 0);
    CHECK_APPROX(out, -0.5, 1e-10);
    CHECK(gd_tsv_parse_double("1e5", 3, &out) == 0);
    CHECK_APPROX(out, 100000.0, 1e-5);
    /* Missing */
    CHECK(gd_tsv_parse_double(".", 1, &out) == GD_E_MISSING_FIELD);
    /* Bad numeric */
    CHECK(gd_tsv_parse_double("abc", 3, &out) == GD_E_BAD_NUMERIC);
    CHECK(gd_tsv_parse_double("12x", 3, &out) == GD_E_BAD_NUMERIC);
}

TEST(test_tsv_parse_int)
{
    long out;
    CHECK(gd_tsv_parse_int("42", 2, &out) == 0);
    CHECK(out == 42);
    CHECK(gd_tsv_parse_int("-7", 2, &out) == 0);
    CHECK(out == -7);
    CHECK(gd_tsv_parse_int(".", 1, &out) == GD_E_MISSING_FIELD);
    CHECK(gd_tsv_parse_int("xyz", 3, &out) == GD_E_BAD_NUMERIC);
}

TEST(test_tsv_many_fields)
{
    /* 10 fields. */
    const char *data = "a\tb\tc\td\te\tf\tg\th\ti\tj\n1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n";
    gd_error err = GD_ERROR_INIT;
    collector c = {0};
    int rc = gd_tsv_parse(data, strlen(data), header_cb, row_cb, &c, &err);
    CHECK(rc == 0);
    CHECK(c.last_nfields == 10);
}

TEST(test_tsv_crlf)
{
    const char *data = "a\tb\r\n1\t2\r\n";
    gd_error err = GD_ERROR_INIT;
    collector c = {0};
    int rc = gd_tsv_parse(data, strlen(data), header_cb, row_cb, &c, &err);
    CHECK(rc == 0);
    CHECK(c.data_rows == 1);
    CHECK(c.last_nfields == 2);
}

TEST(test_tsv_empty_fields)
{
    /* Empty fields (two consecutive tabs). */
    const char *data = "a\tb\tc\n1\t\t3\n";
    gd_error err = GD_ERROR_INIT;
    collector c = {0};
    int rc = gd_tsv_parse(data, strlen(data), header_cb, row_cb, &c, &err);
    CHECK(rc == 0);
    CHECK(c.data_rows == 1);
    CHECK(c.last_nfields == 3);
}

static void run_tests(void)
{
    RUN(test_tsv_basic_parse);
    RUN(test_tsv_comments_and_blanks);
    RUN(test_tsv_missing_header);
    RUN(test_tsv_missing_dot);
    RUN(test_tsv_field_eq);
    RUN(test_tsv_parse_double);
    RUN(test_tsv_parse_int);
    RUN(test_tsv_many_fields);
    RUN(test_tsv_crlf);
    RUN(test_tsv_empty_fields);
}

TEST_MAIN()

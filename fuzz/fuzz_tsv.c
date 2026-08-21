/*
 * fuzz_tsv.c — the TSV parser against arbitrary bytes.
 *
 * Every dataset the project reads is TSV: forms, calibrations, priors,
 * concepts. The parser is callback-driven, so this harness supplies callbacks
 * that touch every field of every row — a parser can only be shown safe
 * against the way its output is actually used, and "hand back a field pointer
 * and a length that do not agree" is precisely the bug a status-only check
 * would miss.
 *
 * `path == NULL` keeps the filesystem out of it. Rows that fail are as
 * interesting as rows that parse; only a crash, a leak or undefined behaviour
 * is a failure here.
 */
#include "gd.h"
#include "gd_tsv.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    size_t rows;
    size_t cols_seen;
} counters;

/* Read every field the parser reports, through its own length rather than by
 * assuming termination — the header promises fields are NOT NUL-terminated in
 * place, and this is where that promise gets tested. */
static int touch_row(const gd_tsv_row *row, void *user, gd_error *err)
{
    counters *c = (counters *)user;
    (void)err;

    c->rows++;
    for (size_t i = 0; i < row->nfields; i++) {
        const char *f = row->fields[i];
        size_t n = row->lengths[i];
        volatile unsigned char sink = 0;

        for (size_t k = 0; k < n; k++) sink = (unsigned char)f[k];
        (void)sink;

        /* The field helpers are part of the parser's surface and take the same
         * (pointer, length) pair, so they are fuzzed alongside it. */
        (void)gd_tsv_is_missing(f, n);
        (void)gd_tsv_field_eq(f, n, "language_id");
        double d;
        long l;
        (void)gd_tsv_parse_double(f, n, &d);
        (void)gd_tsv_parse_int(f, n, &l);
        c->cols_seen++;
    }
    return GD_E_OK;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    gd_error err = GD_ERROR_INIT;
    counters c = { 0, 0 };

    (void)gd_tsv_parse((const char *)data, size, touch_row, touch_row,
                       &c, &err);
    return 0;
}

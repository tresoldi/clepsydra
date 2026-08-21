/*
 * gd_tsv.c — strict TSV parser for Clepsydra.
 *
 * Reads the buffer into a working copy, then scans line by line,
 * splitting on tab characters. No quoting, no multiline fields. Lines
 * starting with '#' are comments; blank lines are ignored. The header row
 * is mandatory.
 */
#include "gd_tsv.h"
#include "gd_memory.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

int gd_tsv_parse(const char *in_buf, size_t in_len,
                 gd_tsv_row_cb header_cb,
                 gd_tsv_row_cb row_cb,
                 void *user, gd_error *err)
{
    char *buf = NULL;
    size_t buf_len = 0;

    if (in_buf == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "tsv: no input");
    }
    /* Work on a copy so we can NUL-terminate safely. */
    if (in_len > SIZE_MAX - 1) {
        return gd_error_set(err, GD_E_NOMEM, "tsv: buffer too large");
    }
    buf = (char *)gd_malloc_checked(in_len + 1, err);
    if (buf == NULL) {
        return GD_E_NOMEM;
    }
    memcpy(buf, in_buf, in_len);
    buf[in_len] = '\0';
    buf_len = in_len;

    /* Determine max fields by counting tabs in the longest line. */
    /* We'll allocate field arrays dynamically per row; cap at 256 fields. */
    size_t max_fields = 256;

    const char **fields = (const char **)gd_malloc_checked(
        max_fields * sizeof(char *), err);
    size_t *lengths = (size_t *)gd_malloc_checked(
        max_fields * sizeof(size_t), err);
    if (fields == NULL || lengths == NULL) {
        gd_free(fields);
        gd_free(lengths);
        gd_free(buf);
        return GD_E_NOMEM;
    }

    int result = 0;
    size_t line_no = 0;
    char *p = buf;
    char *line_start = buf;
    int header_done = 0;
    gd_tsv_row row;
    row.fields = fields;
    row.lengths = lengths;

    while (p <= buf + buf_len) {
        bool is_end = (p == buf + buf_len);
        char c = is_end ? '\n' : *p;

        if (c == '\n' || is_end) {
            /* End of line. Process line_start .. p-1 (excluding the newline).
             * Handle \r\n: if the last char before p is \r, trim it. */
            size_t line_len = (size_t)(p - line_start);
            if (line_len > 0 && line_start[line_len - 1] == '\r') {
                line_len--;
            }
            line_no++;

            /* Skip blank lines and comment lines. */
            bool is_blank = true;
            for (size_t i = 0; i < line_len; i++) {
                if (!isspace((unsigned char)line_start[i])) {
                    is_blank = false;
                    break;
                }
            }
            bool is_comment = (line_len > 0 && line_start[0] == '#');

            if (!is_blank && !is_comment) {
                /* Split on tabs. */
                size_t nf = 0;
                char *fstart = line_start;
                for (size_t i = 0; i <= line_len; i++) {
                    if (i == line_len || line_start[i] == '\t') {
                        if (nf >= max_fields) {
                            gd_error_set(err, GD_E_PARSE,
                                         "tsv: too many fields");
                            result = GD_E_PARSE;
                            goto done;
                        }
                        fields[nf] = fstart;
                        lengths[nf] = (size_t)(&line_start[i] - fstart);
                        nf++;
                        if (i < line_len) {
                            fstart = &line_start[i + 1];
                        }
                    }
                }
                row.nfields = nf;
                row.line_no = line_no;

                if (!header_done) {
                    if (header_cb != NULL) {
                        result = header_cb(&row, user, err);
                        if (result != 0) {
                            goto done;
                        }
                    }
                    header_done = 1;
                } else {
                    if (row_cb != NULL) {
                        result = row_cb(&row, user, err);
                        if (result != 0) {
                            goto done;
                        }
                    }
                }
            }

            line_start = p + 1;
            if (is_end) {
                break;
            }
        }
        p++;
    }

    if (!header_done) {
        gd_error_set(err, GD_E_MISSING_HEADER, "tsv: no header row found");
        result = GD_E_MISSING_HEADER;
    }

done:
    gd_free(fields);
    gd_free(lengths);
    gd_free(buf);
    return result;
}

size_t gd_tsv_colmap_find(const gd_tsv_colmap *map, size_t n,
                          const char *name)
{
    if (map == NULL || name == NULL) {
        return SIZE_MAX;
    }
    for (size_t i = 0; i < n; i++) {
        if (strcmp(map[i].name, name) == 0) {
            return map[i].index;
        }
    }
    return SIZE_MAX;
}

bool gd_tsv_is_missing(const char *field, size_t len)
{
    return (len == 1 && field[0] == '.');
}

bool gd_tsv_field_eq(const char *field, size_t len, const char *str)
{
    size_t slen = strlen(str);
    if (slen != len) {
        return false;
    }
    return memcmp(field, str, len) == 0;
}

/* Parse a double from a field that is not NUL-terminated. We make a temp
 * copy. This is called infrequently (not in the inner loop). */
int gd_tsv_parse_double(const char *field, size_t len, double *out)
{
    if (gd_tsv_is_missing(field, len)) {
        return GD_E_MISSING_FIELD;
    }
    /* Copy to a small stack buffer. */
    if (len > 63) {
        return GD_E_BAD_NUMERIC;
    }
    char tmp[64];
    memcpy(tmp, field, len);
    tmp[len] = '\0';

    /* Validate: allow optional leading +/-, digits, one '.', optional
     * exponent. strtoend checks via strtod endptr.
     * Note: strtod is locale-dependent, so under a comma-decimal locale a
     * value like "3.14" is *rejected* (endptr stops at '.', → GD_E_BAD_NUMERIC),
     * never silently misparsed — there is no data-corruption path, only a hard
     * error the caller sees. Input files are ASCII '.'-decimal by spec. */
    char *endptr = NULL;
    double val = strtod(tmp, &endptr);
    if (endptr == tmp || *endptr != '\0') {
        return GD_E_BAD_NUMERIC;
    }
    *out = val;
    return 0;
}

int gd_tsv_parse_int(const char *field, size_t len, long *out)
{
    if (gd_tsv_is_missing(field, len)) {
        return GD_E_MISSING_FIELD;
    }
    if (len > 31) {
        return GD_E_BAD_NUMERIC;
    }
    char tmp[32];
    memcpy(tmp, field, len);
    tmp[len] = '\0';
    char *endptr = NULL;
    long val = strtol(tmp, &endptr, 10);
    if (endptr == tmp || *endptr != '\0') {
        return GD_E_BAD_NUMERIC;
    }
    *out = val;
    return 0;
}

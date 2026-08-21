/*
 * gd_tsv.h — internal strict TSV parser for Clepsydra.
 *
 * Brief §1.6 requirements:
 *  - input from an in-memory buffer (the caller does its own file I/O);
 *  - headers mandatory;
 *  - fields separated by literal tab;
 *  - UTF-8 passed through (no normalization);
 *  - '.' = missing;
 *  - blank lines ignored;
 *  - lines beginning with '#' are comments;
 *  - malformed numeric fields produce useful errors;
 *  - duplicate IDs detected (by caller via callbacks);
 *  - no quoted/multiline fields;
 *  - never exit()/abort().
 *
 * The parser is a zero-allocation streaming scanner: it calls back per row
 * with field pointers into the working copy. Callers own validation logic.
 */
#ifndef GD_TSV_H
#define GD_TSV_H

#include "gd_types.h"
#include "gd_error.h"

#include <stdbool.h>
#include <stddef.h>

/* A parsed row: array of field pointers (into the parser's buffer) and
 * their lengths. Fields are NOT NUL-terminated in place; use the length. */
typedef struct gd_tsv_row {
    const char **fields;   /* pointers into the buffer */
    size_t      *lengths;  /* length of each field */
    size_t       nfields;
    size_t       line_no;  /* 1-based line number in the source */
} gd_tsv_row;

/* Callback type: called for each data row. Return non-zero to stop parsing
 * (treated as an error). */
typedef int (*gd_tsv_row_cb)(const gd_tsv_row *row, void *user, gd_error *err);

/* Column map: maps a header name to its index. Built from the header row. */
typedef struct gd_tsv_colmap {
    const char *name;
    size_t      index;
} gd_tsv_colmap;

/* Parse an in-memory TSV buffer. `len` is a byte count: 0 means an empty
 * buffer, which is a real input, and no length is ever measured here — a
 * caller holding the public GD_LEN_NUL_TERMINATED sentinel resolves it before
 * calling (gd_analysis.c does). Calls header_cb (may be NULL) for the header
 * row, then row_cb for each data row. Returns GD_E_OK on success, a GD_E_*
 * code on failure.
 *
 * The parser never opens a file: nothing in src/ does I/O, so the same engine
 * runs natively, in a fuzzer, and in a browser with no filesystem. */
int gd_tsv_parse(const char *buf, size_t len,
                 gd_tsv_row_cb header_cb,
                 gd_tsv_row_cb row_cb,
                 void *user, gd_error *err);

/* Look up a column index by name in a colmap of n entries. Returns SIZE_MAX
 * if not found. */
size_t gd_tsv_colmap_find(const gd_tsv_colmap *map, size_t n,
                          const char *name);

/* Parse a numeric field. On success returns 0 and sets *out. On failure
 * returns GD_E_BAD_NUMERIC. '.' is treated as missing: returns GD_E_MISSING_FIELD
 * and does not set *out (caller decides how to handle). */
int gd_tsv_parse_double(const char *field, size_t len, double *out);
int gd_tsv_parse_int(const char *field, size_t len, long *out);

/* Is this field the "missing" marker '.'? A predicate. */
bool gd_tsv_is_missing(const char *field, size_t len);

/* Does this field equal the given C string? A predicate. */
bool gd_tsv_field_eq(const char *field, size_t len, const char *str);

#endif /* GD_TSV_H */

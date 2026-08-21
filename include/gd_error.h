/*
 * gd_error.h — the error model for Clepsydra (libclepsydra).
 *
 * ===========================================================================
 * The convention. One rule, no exceptions:
 *
 *   A function returning `int` returns a GD_E_* code. GD_E_OK (0) is success;
 *   anything else is that specific failure.
 *
 * Not "nonzero means it went wrong". The value IS the classification, and it
 * matches `err->code` whenever an err was supplied. A caller may branch on
 * `!= GD_E_OK` and ignore the detail, but a caller that switches on the code
 * gets the truth.
 *
 * Corollaries:
 *
 *   - Never return a bare 1 or -1 to mean "failed". 1 is GD_E_NOMEM and -1 is
 *     not a code at all; both lie to a caller reading the contract.
 *   - A function that answers a question rather than performing an operation
 *     returns `bool`, so it is visibly not a status function. gd_error_failed
 *     and gd_tsv_is_missing are predicates; gd_tsv_parse is not.
 *   - A function returning a pointer reports failure with NULL and puts the
 *     classification in `err`.
 *   - assert() is for programmer invariants only — conditions a caller cannot
 *     produce through the API. Bad input from a file, a user or a browser is
 *     never an assertion; it is a returned code.
 * ===========================================================================
 *
 * The library never calls exit() or abort(). Every operation that can fail
 * takes a gd_error* and reports through it. Passing NULL for the gd_error* is
 * always allowed: the code is still returned, only the detail is dropped.
 */
#ifndef GD_ERROR_H
#define GD_ERROR_H

#include "gd_types.h"

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GD_E_OK              0
#define GD_E_NOMEM           1   /* allocation failure / checked-mult overflow */
#define GD_E_INVALID_ARG     2
#define GD_E_PARSE           3   /* TSV / Newick parse error */
#define GD_E_MISSING_HEADER  4
#define GD_E_UNKNOWN_COLUMN  5   /* unknown required column */
#define GD_E_DUPLICATE_ID    6
#define GD_E_MISSING_FIELD   7   /* required field is '.' or empty */
#define GD_E_BAD_NUMERIC     8   /* malformed numeric field */
#define GD_E_INCOMPATIBLE    9   /* incompatible concept under incompatible=error */
#define GD_E_CALIBRATION     10  /* infeasible / conflicting calibrations */
#define GD_E_NOT_ANCHORED    11  /* absolute time not identifiable (informational) */
#define GD_E_RANGE           12  /* value out of feasible range */
#define GD_E_NOT_IMPLEMENTED 13  /* reserved for unimplemented (e.g. nonbinary T1 exact) */
#define GD_E_CHECKPOINT      14  /* version/checksum mismatch on restore */
#define GD_E_INTERNAL        15
#define GD_E_UNSUPPORTED     16
#define GD_E_CANCELLED       17  /* a progress callback asked the run to stop */
/*
 * A caller-supplied sink failed to accept bytes (gd_report.h). The library
 * itself never raises this -- it does no I/O -- but a renderer has to classify
 * a sink that said no, and "internal" would be a lie about whose fault it is.
 */
#define GD_E_IO              18

/*
 * Initialiser for a gd_error on the stack:
 *
 *     gd_error err = GD_ERROR_INIT;
 *
 * Use this rather than a positional brace list, which silently goes wrong when
 * the struct gains or loses a field.
 */
#define GD_ERROR_INIT { GD_E_OK, { '\0' } }

/*
 * Format-attribute shim. gcc and clang can check a printf-style call against
 * its format string; nothing else needs to. Isolated here rather than spread
 * through the code, so the library stays portable C99 with no extension in
 * sight at any call site.
 */
#if defined(__GNUC__)
#  define GD_PRINTF_LIKE(fmt_idx, first_arg) \
       __attribute__((format(printf, fmt_idx, first_arg)))
#else
#  define GD_PRINTF_LIKE(fmt_idx, first_arg)
#endif

/*
 * Record a failure. `message` is COPIED into the error's own storage, so it may
 * point at a stack buffer, and the caller may read err->message for as long as
 * the gd_error itself lives. A message longer than GD_ERROR_MSG_MAX-1 bytes is
 * truncated rather than rejected — a truncated explanation beats none.
 *
 * `message` may be NULL, which leaves the message empty. `err` may be NULL, in
 * which case nothing is recorded. Returns `code`, so a failing branch can be
 * written as:
 *
 *     return gd_error_set(err, GD_E_PARSE, "unterminated quote");
 */
int  gd_error_set(gd_error *err, int code, const char *message);

/* As gd_error_set, but builds the message from a printf-style format. This is
 * how an error names a value ("bad calibration type 'xyz'") — the formatted
 * text lands in the error's own storage, so no caller has to supply a buffer
 * for it. Returns `code`. */
int  gd_error_setf(gd_error *err, int code, const char *fmt, ...)
     GD_PRINTF_LIKE(3, 4);

/* Reset to success with an empty message. NULL is accepted. */
void gd_error_clear(gd_error *err);

/* Did this error record a failure? A predicate, not a status. */
bool gd_error_failed(const gd_error *err);

/* The recorded message, or "" when none was set. Never NULL, so it is always
 * safe to print. Borrowed from `err` and valid as long as it is. */
const char *gd_error_message(const gd_error *err);

#ifdef __cplusplus
}
#endif
#endif /* GD_ERROR_H */

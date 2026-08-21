/*
 * gd_error.c — error reporting for Clepsydra.
 *
 * The library never calls exit()/abort(); all errors flow through gd_error.
 * The convention every int-returning function follows is documented at the top
 * of gd_error.h.
 *
 * Messages are copied into the error's own storage rather than borrowed. That
 * costs a bounded copy on a path that has already failed, and buys the thing
 * the borrowed-pointer design could not give: a message may name a runtime
 * value ("bad calibration type 'xyz'") without the caller having to supply
 * somewhere for it to live.
 */
#include "gd_error.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/* Copy at most GD_ERROR_MSG_MAX-1 bytes and always terminate. strncpy is not
 * used: it does not terminate on truncation, which is the one case that
 * matters here. */
static void msg_copy(char *dst, const char *src)
{
    size_t n;

    /* Programmer invariant, not input validation: dst is always err->message,
     * which is an array inside a struct the caller already proved non-NULL.
     * A NULL here would mean this file called itself wrongly. */
    assert(dst != NULL);

    if (src == NULL) {
        dst[0] = '\0';
        return;
    }
    n = strlen(src);
    if (n > GD_ERROR_MSG_MAX - 1) {
        n = GD_ERROR_MSG_MAX - 1;
    }
    memcpy(dst, src, n);
    dst[n] = '\0';
}

int gd_error_set(gd_error *err, int code, const char *message)
{
    if (err != NULL) {
        err->code = code;
        msg_copy(err->message, message);
    }
    return code;
}

int gd_error_setf(gd_error *err, int code, const char *fmt, ...)
{
    va_list ap;

    if (err == NULL) {
        return code;
    }
    err->code = code;
    if (fmt == NULL) {
        err->message[0] = '\0';
        return code;
    }
    va_start(ap, fmt);
    /* vsnprintf always terminates within the given size and truncates rather
     * than overflowing; a negative return means an encoding error, and an
     * empty message is the honest answer to that. */
    if (vsnprintf(err->message, sizeof err->message, fmt, ap) < 0) {
        err->message[0] = '\0';
    }
    va_end(ap);
    return code;
}

void gd_error_clear(gd_error *err)
{
    if (err == NULL) {
        return;
    }
    err->code = GD_E_OK;
    err->message[0] = '\0';
}

bool gd_error_failed(const gd_error *err)
{
    return (err != NULL) && (err->code != GD_E_OK);
}

const char *gd_error_message(const gd_error *err)
{
    return (err != NULL) ? err->message : "";
}

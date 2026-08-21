/*
 * gd_memory.c — checked allocation for Clepsydra.
 *
 * Design (brief §14.2):
 *  - checked multiplication before allocation;
 *  - no VLAs depending on user input;
 *  - no compiler extensions.
 *
 * There was an arena here — a bump allocator over a block list, for
 * analysis-lifetime memory. It was removed in the C3 alignment milestone
 * because nothing ever allocated from it: every module owns its allocations
 * through the object that outlives them, which is a simpler lifetime model
 * than the arena would have imposed. An allocator abstraction has to be
 * earned by a real lifetime problem, and this one never was. It is in the
 * git history if a future milestone finds that problem.
 */
#include "gd_memory.h"
#include "gd_error.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------ heap allocation */

void *gd_malloc_checked(size_t nbytes, gd_error *err)
{
    void *p = malloc(nbytes);
    if (p == NULL && nbytes > 0) {
        gd_error_set(err, GD_E_NOMEM, "malloc failed");
    }
    return p;
}

void *gd_calloc_checked(size_t count, size_t elem_size, gd_error *err)
{
    if (count != 0 && elem_size > SIZE_MAX / count) {
        gd_error_set(err, GD_E_NOMEM, "calloc: size overflow");
        return NULL;
    }
    void *p = calloc(count, elem_size);
    if (p == NULL && count > 0 && elem_size > 0) {
        gd_error_set(err, GD_E_NOMEM, "calloc failed");
    }
    return p;
}

bool gd_size_mul(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a) {
        return false;
    }
    *out = a * b;
    return true;
}

void *gd_malloc_array_checked(size_t count, size_t elem_size, gd_error *err)
{
    if (count != 0 && elem_size > SIZE_MAX / count) {
        gd_error_set(err, GD_E_NOMEM, "malloc: array size overflow");
        return NULL;
    }
    void *p = malloc(count * elem_size);
    if (p == NULL && count > 0 && elem_size > 0) {
        gd_error_set(err, GD_E_NOMEM, "malloc failed");
    }
    return p;
}

void *gd_realloc_array_checked(void *ptr, size_t count, size_t elem_size,
                               gd_error *err)
{
    if (count != 0 && elem_size > SIZE_MAX / count) {
        gd_error_set(err, GD_E_NOMEM, "realloc: array size overflow");
        return NULL;
    }
    void *p = realloc(ptr, count * elem_size);
    if (p == NULL && count > 0 && elem_size > 0) {
        gd_error_set(err, GD_E_NOMEM, "realloc failed");
    }
    return p;
}

void gd_free(void *ptr)
{
    free(ptr);
}

/* ---------------------------------------------------- string duplication */

char *gd_strdup_n(const char *src, size_t len, gd_error *err)
{
    if (src == NULL) {
        return NULL;
    }
    /* Check for overflow: len + 1. */
    if (len > SIZE_MAX - 1) {
        gd_error_set(err, GD_E_NOMEM, "strdup: length overflow");
        return NULL;
    }
    char *dst = (char *)malloc(len + 1);
    if (dst == NULL) {
        gd_error_set(err, GD_E_NOMEM, "strdup: malloc failed");
        return NULL;
    }
    memcpy(dst, src, len);
    dst[len] = '\0';
    return dst;
}

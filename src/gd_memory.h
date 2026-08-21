/*
 * gd_memory.h — internal memory API for Clepsydra.
 *
 * Checked heap allocation. Internal (not part of the public API): it lives in
 * src/ alongside the implementation and is included by other src/ files.
 *
 * The strategy is deliberately plain — direct heap ownership, with each
 * allocation owned by the object that outlives it and released in that
 * object's _destroy. No arena, no allocator callbacks: the lifetimes here are
 * simple enough that either would be abstraction without a problem to solve.
 */
#ifndef GD_MEMORY_H
#define GD_MEMORY_H

#include "gd_types.h"
#include "gd_error.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/*
 * Heap helpers (caller frees with gd_free).
 *
 * The array forms exist because `malloc(n * sizeof(T))` is the shape that
 * silently wraps: every `n` here comes from parsed input — concepts, nodes,
 * draws — so a hostile or absurd dataset could otherwise produce a buffer far
 * smaller than the loop that fills it believes. Taking count and element size
 * separately means the multiplication is checked where it happens.
 *
 * All of them accept a NULL `err`: the NULL return still says it failed.
 */
/*
 * *out = a * b, or false on overflow with *out untouched. For the handful of
 * allocations whose element count is itself a product (chains x draws x
 * parameters), where checking only the final multiply would leave the first
 * one able to wrap.
 */
bool      gd_size_mul(size_t a, size_t b, size_t *out);

void     *gd_malloc_checked(size_t nbytes, gd_error *err);

/* count * elem_size bytes, uninitialised. NULL on overflow or failure. */
void     *gd_malloc_array_checked(size_t count, size_t elem_size, gd_error *err);

/* count * elem_size bytes, zeroed. NULL on overflow or failure. */
void     *gd_calloc_checked(size_t count, size_t elem_size, gd_error *err);

/*
 * Grow (or shrink) an array to count * elem_size bytes. On failure returns
 * NULL and leaves `ptr` valid and owned by the caller — as realloc does, so
 * the usual "assign to a new variable, and on NULL the old one is still mine"
 * pattern is unchanged.
 */
void     *gd_realloc_array_checked(void *ptr, size_t count, size_t elem_size,
                                   gd_error *err);

void      gd_free(void *ptr);

/* Duplicate at most len bytes; NUL-terminated. Caller frees with gd_free. */
char     *gd_strdup_n(const char *src, size_t len, gd_error *err);

#endif /* GD_MEMORY_H */

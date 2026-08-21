/*
 * gd_newick.h — Newick tree parser for Clepsydra.
 *
 * Brief §10.8: accept rooted Newick with labelled leaves and labelled
 * internal nodes. The parser supports:
 *  - leaf and internal labels;
 *  - branch lengths (after ':');
 *  - polytomies (nonbinary);
 *  - comments in [...] only if simple to support safely (rejected if nested
 *    or malformed — never silently misparsed).
 * Branch lengths may be ignored / used as starting values / interpreted as
 * durations per an explicit configuration (never guessed).
 *
 * Output: a gd_tree structure (see gd_tree.h) with nodes stored in a flat
 * array; traversal uses iterative postorder/preorder (no deep recursion).
 */
#ifndef GD_NEWICK_H
#define GD_NEWICK_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_tree.h"

#include <stddef.h>

/* Parse Newick from an in-memory buffer. `len` is a byte count; no length is
 * measured here (see gd_tsv.h). On success, *out points to a heap-allocated
 * gd_tree, which the caller frees with gd_tree_destroy. Returns GD_E_OK on
 * success.
 *
 * The parser never opens a file — see gd_tsv.h. */
int gd_newick_parse(const char *buf, size_t len,
                    gd_tree **out, gd_error *err);

#endif /* GD_NEWICK_H */

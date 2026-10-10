/*
 * gd_tree_build.h — constructing a gd_tree. Internal.
 *
 * The counterpart to the public `include/gd_tree.h`, which exposes only the
 * read surface. Everything that allocates, mutates or precomputes lives here,
 * because only the Newick parser, the calibration resolver, the likelihoods
 * and the tests ever build a tree — a library caller receives one already
 * built, owned by the result that handed it over.
 *
 * Ownership: gd_tree_create returns an owned tree; gd_tree_destroy frees it
 * along with every node label and child array it holds. A tree that has been
 * moved into a gd_analysis result is destroyed by that result, not here.
 */
#ifndef GD_TREE_BUILD_H
#define GD_TREE_BUILD_H

#include "gd_tree.h"

/* Allocate a tree with `nnodes` nodes, each initialised to a childless,
 * unlabelled, parentless leaf. Returns NULL and sets `err` on failure. */
gd_tree *gd_tree_create(size_t nnodes, gd_error *err);

/* Free a tree and everything it owns. NULL is accepted. */
void     gd_tree_destroy(gd_tree *t);

/* A deep copy with the same node indices, labels, ages and branch lengths;
 * traversals are rebuilt, and LCA support too when the source had it. Used to
 * give each parallel MCMC chain a tree of its own to move ages in. */
gd_tree *gd_tree_copy(const gd_tree *src, gd_error *err);

/* Add a child to a parent node, growing the children array as needed. */
int gd_tree_add_child(gd_tree *t, size_t parent_idx, size_t child_idx,
                      gd_error *err);

/* Build traversal orders (postorder, preorder) iteratively. */
int gd_tree_build_traversals(gd_tree *t, gd_error *err);

/* Build LCA support (depth, Euler tour, RMQ). Must be called after
 * traversals; gd_tree_mrca and gd_tree_depth need it. */
int gd_tree_build_lca(gd_tree *t, gd_error *err);

#endif /* GD_TREE_BUILD_H */

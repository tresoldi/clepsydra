/*
 * gd_tree.h — the phylogeny a result is reported against.
 *
 * Public because a chronology is only meaningful against its topology: a
 * caller reading `gd_analysis_tree_result` needs to walk the tree to know
 * which node an age belongs to, and both shipped consumers do exactly that
 * (the CLI resolves labels, the browser page draws the topology).
 *
 * The struct is exposed rather than made opaque on purpose. It is a flat
 * array of nodes with parent/child indices — value-like data a caller reads
 * directly — and wrapping every field in an accessor would buy nothing.
 *
 * What is NOT public is construction. Building a tree, growing its child
 * arrays and computing its traversal and LCA tables are the parser's and the
 * sampler's business; those declarations live in `src/gd_tree_build.h`. A
 * caller receives a tree that is already built and owned by the result that
 * handed it over, and destroys it through that result, never directly.
 *
 * Nodes are stored flat (cache-friendly) and traversal orders are precomputed
 * as index arrays so every walk is iterative — no deep recursion, because the
 * browser stack is small.
 */
#ifndef GD_TREE_H
#define GD_TREE_H

#include "gd_types.h"
#include "gd_error.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* "No such node": the parent of the root, and the result of a failed lookup. */
#define GD_NODE_NONE  SIZE_MAX

typedef struct gd_node {
    size_t  parent;        /* index of parent, or GD_NODE_NONE for root */
    size_t *children;      /* array of child indices */
    size_t  nchildren;     /* number of children */
    size_t  child_cap;     /* allocated capacity */
    char   *label;         /* node label (may be NULL); owned by the tree */
    double  branch_length; /* from parent to this node; 0 for root; -1 if absent */
    bool    is_leaf;
    /* Age (years BP): set by calibration/MCMC, not by the parser. */
    double  age;
} gd_node;

typedef struct gd_tree {
    gd_node *nodes;
    size_t   nnodes;
    size_t   root;         /* index of root node */
    /* Precomputed traversal orders (indices into nodes[]). */
    size_t  *postorder;    /* length nnodes */
    size_t  *preorder;     /* length nnodes */
    /* Depth/LCA support. */
    size_t  *depth;        /* depth of each node (root = 0) */
    /* Euler tour + RMQ for LCA (built on demand). */
    size_t  *euler;        /* Euler tour: node indices, length 2*nnodes-1 */
    size_t  *euler_depth;  /* depth at each Euler position */
    size_t  *first_visit;  /* first Euler position of each node */
    size_t  *rmq_table;    /* RMQ sparse table for LCA; rows * cols */
    size_t   rmq_rows;
    size_t   rmq_cols;
    bool     lca_built;
} gd_tree;

/* ------------------------------------------------------------- reading */

/* Find a node by label. Returns its index, or GD_NODE_NONE if no node carries
 * that label. `label` is compared exactly; NULL labels never match. */
size_t gd_tree_find_label(const gd_tree *t, const char *label);

/* Compute the MRCA (LCA) of two nodes. Requires the tree's LCA support to have
 * been built, which every tree from the parser has. */
size_t gd_tree_mrca(const gd_tree *t, size_t a, size_t b);

/* Depth of a node, root = 0. Returns GD_NODE_NONE for an out-of-range index or
 * a tree without depth information. */
static inline size_t gd_tree_depth(const gd_tree *t, size_t idx)
{
    return (t != NULL && t->depth != NULL && idx < t->nnodes)
        ? t->depth[idx] : GD_NODE_NONE;
}

/* Is node a leaf? Out-of-range indices answer 0 rather than trapping. */
static inline int gd_tree_is_leaf(const gd_tree *t, size_t idx)
{
    return (t != NULL && idx < t->nnodes) ? t->nodes[idx].is_leaf : 0;
}

#endif /* GD_TREE_H */

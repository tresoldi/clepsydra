/*
 * fuzz_newick.c — the Newick parser against arbitrary bytes.
 *
 * Why this one first: the browser page lets a visitor paste a tree into a
 * textarea, and that text reaches gd_newick_parse unaltered. It is the
 * shortest path from a stranger's keyboard to this C, so it is the one that
 * has to survive anything.
 *
 * The harness calls the production parser directly — no reimplementation, no
 * pre-validation — because a harness that cleans its input first only proves
 * the cleaner works. It passes `path == NULL` so nothing touches the
 * filesystem, which is exactly the property gd_newick_parse's buffer form
 * already had; that is why this file is twenty lines rather than a rewrite.
 *
 * Success and failure are both fine. The only unacceptable outcomes are a
 * crash, a leak, or undefined behaviour, which the sanitizers this is built
 * with will report.
 */
#include "gd.h"
#include "gd_newick.h"
#include "gd_tree_build.h"

#include <stddef.h>
#include <stdint.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *tree = NULL;

    if (gd_newick_parse((const char *)data, size, &tree, &err) == GD_E_OK) {
        /* On success the tree is owned here. Reading it back is part of the
         * test: a parser that returns a malformed tree passes a check that
         * only looks at the status. */
        if (tree != NULL) {
            for (size_t i = 0; i < tree->nnodes; i++) {
                const gd_node *n = &tree->nodes[i];
                (void)n->parent;
                (void)n->branch_length;
                (void)gd_tree_is_leaf(tree, i);
                (void)gd_tree_depth(tree, i);
                for (size_t c = 0; c < n->nchildren; c++) (void)n->children[c];
            }
            (void)gd_tree_find_label(tree, "x");
        }
    }
    gd_tree_destroy(tree);
    return 0;
}

/*
 * gd_tree.c — tree structure, iterative traversal, and LCA for Clepsydra.
 *
 * All traversals are iterative (explicit stack) to avoid deep recursion
 * that would overflow a browser stack. Traversal orders are precomputed
 * into flat index arrays so the inner likelihood loop can simply iterate.
 */
#include "gd_tree.h"
#include "gd_tree_build.h"
#include "gd_memory.h"

#include <stdlib.h>
#include <string.h>

gd_tree *gd_tree_create(size_t nnodes, gd_error *err)
{
    if (nnodes == 0) {
        gd_error_set(err, GD_E_INVALID_ARG, "tree: nnodes must be > 0");
        return NULL;
    }
    gd_tree *t = (gd_tree *)gd_calloc_checked(1, sizeof(gd_tree), err);
    if (t == NULL) {
        return NULL;
    }
    t->nodes = (gd_node *)gd_calloc_checked(nnodes, sizeof(gd_node), err);
    if (t->nodes == NULL) {
        gd_free(t);
        return NULL;
    }
    t->nnodes = nnodes;
    t->root = 0;
    for (size_t i = 0; i < nnodes; i++) {
        t->nodes[i].parent = GD_NODE_NONE;
        t->nodes[i].children = NULL;
        t->nodes[i].nchildren = 0;
        t->nodes[i].child_cap = 0;
        t->nodes[i].label = NULL;
        t->nodes[i].branch_length = -1.0;
        t->nodes[i].is_leaf = false;
        t->nodes[i].age = 0.0;
    }
    t->postorder = NULL;
    t->preorder = NULL;
    t->depth = NULL;
    t->euler = NULL;
    t->euler_depth = NULL;
    t->first_visit = NULL;
    t->rmq_table = NULL;
    t->rmq_rows = 0;
    t->rmq_cols = 0;
    t->lca_built = false;
    return t;
}

void gd_tree_destroy(gd_tree *t)
{
    if (t == NULL) {
        return;
    }
    for (size_t i = 0; i < t->nnodes; i++) {
        gd_free(t->nodes[i].children);
        gd_free(t->nodes[i].label);
    }
    gd_free(t->nodes);
    gd_free(t->postorder);
    gd_free(t->preorder);
    gd_free(t->depth);
    gd_free(t->euler);
    gd_free(t->euler_depth);
    gd_free(t->first_visit);
    gd_free(t->rmq_table);
    gd_free(t);
}

/* Add a child to a parent node, growing the children array as needed. */
int gd_tree_add_child(gd_tree *t, size_t parent_idx, size_t child_idx,
                      gd_error *err)
{
    if (t == NULL || parent_idx >= t->nnodes || child_idx >= t->nnodes) {
        gd_error_set(err, GD_E_INVALID_ARG, "tree: add_child bad index");
        return GD_E_INVALID_ARG;
    }
    gd_node *p = &t->nodes[parent_idx];
    if (p->nchildren >= p->child_cap) {
        size_t new_cap = (p->child_cap == 0) ? 4 : p->child_cap * 2;
        /* Check overflow. */
        if (new_cap > SIZE_MAX / sizeof(size_t)) {
            gd_error_set(err, GD_E_NOMEM, "tree: child cap overflow");
            return GD_E_NOMEM;
        }
        size_t *new_children = (size_t *)gd_realloc_array_checked(p->children,
            new_cap, sizeof(size_t), err);
        if (new_children == NULL) {
            gd_error_set(err, GD_E_NOMEM, "tree: realloc children failed");
            return GD_E_NOMEM;
        }
        p->children = new_children;
        p->child_cap = new_cap;
    }
    p->children[p->nchildren++] = child_idx;
    t->nodes[child_idx].parent = parent_idx;
    t->nodes[child_idx].is_leaf = false;  /* will be set after all children added */
    return 0;
}

/* Finalize: mark leaves (nodes with 0 children). Call after building. */
static void gd_tree_finalize_leaves(gd_tree *t)
{
    for (size_t i = 0; i < t->nnodes; i++) {
        t->nodes[i].is_leaf = (t->nodes[i].nchildren == 0);
    }
}

/*
 * Iterative postorder traversal using an explicit stack.
 * postorder: children left-to-right, then the node itself.
 */
int gd_tree_build_traversals(gd_tree *t, gd_error *err)
{
    if (t == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "tree: NULL");
        return GD_E_INVALID_ARG;
    }
    gd_tree_finalize_leaves(t);

    /*
     * Connectivity, checked once here rather than defended against everywhere.
     *
     * A tree has exactly one node without a parent, and it is the root. Code
     * downstream reads a node's parent whenever the node is not the root —
     * five places do — and a second parentless node would send every one of
     * them indexing at GD_NODE_NONE. Found by fuzzing: the Newick parser used
     * to accept "a,(b);", a string with two top-level nodes and no enclosing
     * parentheses, and produced exactly that shape.
     *
     * This is the boundary for it: every consumer calls build_traversals after
     * construction, whether the tree came from the parser or was assembled by
     * hand in a test.
     */
    if (t->root >= t->nnodes) {
        return gd_error_set(err, GD_E_PARSE, "tree: root index out of range");
    }
    for (size_t i = 0; i < t->nnodes; i++) {
        if (i == t->root) {
            if (t->nodes[i].parent != GD_NODE_NONE) {
                return gd_error_set(err, GD_E_PARSE, "tree: the root has a parent");
            }
        } else if (t->nodes[i].parent == GD_NODE_NONE) {
            return gd_error_set(err, GD_E_PARSE,
                                "tree: disconnected node (more than one root)");
        } else if (t->nodes[i].parent >= t->nnodes) {
            return gd_error_set(err, GD_E_PARSE, "tree: parent index out of range");
        }
    }

    /* Allocate postorder and preorder arrays. */
    t->postorder = (size_t *)gd_calloc_checked(t->nnodes, sizeof(size_t), err);
    t->preorder  = (size_t *)gd_calloc_checked(t->nnodes, sizeof(size_t), err);
    if (t->postorder == NULL || t->preorder == NULL) {
        return GD_E_NOMEM;
    }

    /* Explicit stack for DFS. Stack entries: (node_index, child_position).
     * We use two parallel arrays for the stack. */
    size_t stack_cap = t->nnodes + 1;
    size_t *stack_node = (size_t *)gd_calloc_checked(stack_cap,
        sizeof(size_t), err);
    size_t *stack_pos  = (size_t *)gd_calloc_checked(stack_cap,
        sizeof(size_t), err);
    if (stack_node == NULL || stack_pos == NULL) {
        gd_free(stack_node);
        gd_free(stack_pos);
        return GD_E_NOMEM;
    }

    /* ---- Preorder: root, then children left-to-right ---- */
    size_t sp = 0;
    size_t pre_idx = 0;
    stack_node[sp] = t->root;
    stack_pos[sp] = 0;
    sp++;
    while (sp > 0) {
        size_t node = stack_node[sp - 1];
        size_t pos = stack_pos[sp - 1];
        if (pos == 0) {
            /* First visit: emit in preorder. */
            t->preorder[pre_idx++] = node;
        }
        if (pos < t->nodes[node].nchildren) {
            stack_pos[sp - 1]++;
            stack_node[sp] = t->nodes[node].children[pos];
            stack_pos[sp] = 0;
            sp++;
        } else {
            sp--;
        }
    }

    /* ---- Postorder: children first, then node ----
     * Reuse the same stack approach but emit on "last visit". */
    sp = 0;
    size_t post_idx = 0;
    stack_node[sp] = t->root;
    stack_pos[sp] = 0;
    sp++;
    while (sp > 0) {
        size_t node = stack_node[sp - 1];
        size_t pos = stack_pos[sp - 1];
        if (pos < t->nodes[node].nchildren) {
            stack_pos[sp - 1]++;
            stack_node[sp] = t->nodes[node].children[pos];
            stack_pos[sp] = 0;
            sp++;
        } else {
            /* All children processed: emit this node. */
            t->postorder[post_idx++] = node;
            sp--;
        }
    }

    gd_free(stack_node);
    gd_free(stack_pos);
    return 0;
}

/* ---- LCA via Euler tour + RMQ ---- */

int gd_tree_build_lca(gd_tree *t, gd_error *err)
{
    if (t == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "tree: NULL");
        return GD_E_INVALID_ARG;
    }
    if (t->postorder == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "tree: build traversals first");
        return GD_E_INVALID_ARG;
    }

    size_t n = t->nnodes;
    if (n == 0) {
        return 0;
    }

    t->depth = (size_t *)gd_calloc_checked(n, sizeof(size_t), err);
    if (t->depth == NULL) {
        return GD_E_NOMEM;
    }

    /* Depth from preorder (parent always before children in preorder). */
    t->depth[t->root] = 0;
    for (size_t i = 0; i < n; i++) {
        size_t node = t->preorder[i];
        for (size_t c = 0; c < t->nodes[node].nchildren; c++) {
            size_t child = t->nodes[node].children[c];
            t->depth[child] = t->depth[node] + 1;
        }
    }

    /* Euler tour: length 2*n - 1. */
    size_t euler_len = 2 * n - 1;
    if (n == 1) {
        euler_len = 1;
    }
    t->euler = (size_t *)gd_calloc_checked(euler_len, sizeof(size_t), err);
    t->euler_depth = (size_t *)gd_calloc_checked(euler_len, sizeof(size_t), err);
    t->first_visit = (size_t *)gd_calloc_checked(n, sizeof(size_t), err);
    if (t->euler == NULL || t->euler_depth == NULL || t->first_visit == NULL) {
        return GD_E_NOMEM;
    }
    for (size_t i = 0; i < n; i++) {
        t->first_visit[i] = GD_NODE_NONE;
    }

    /* Iterative Euler tour (DFS). */
    size_t stack_cap = n + 1;
    size_t *stack_node = (size_t *)gd_calloc_checked(stack_cap,
        sizeof(size_t), err);
    size_t *stack_pos  = (size_t *)gd_calloc_checked(stack_cap,
        sizeof(size_t), err);
    if (stack_node == NULL || stack_pos == NULL) {
        gd_free(stack_node);
        gd_free(stack_pos);
        return GD_E_NOMEM;
    }

    size_t sp = 0;
    size_t ei = 0;  /* euler index */
    stack_node[sp] = t->root;
    stack_pos[sp] = 0;
    sp++;
    t->euler[ei] = t->root;
    t->euler_depth[ei] = 0;
    if (t->first_visit[t->root] == GD_NODE_NONE) {
        t->first_visit[t->root] = ei;
    }
    ei++;

    while (sp > 0) {
        size_t node = stack_node[sp - 1];
        size_t pos = stack_pos[sp - 1];
        if (pos < t->nodes[node].nchildren) {
            stack_pos[sp - 1]++;
            size_t child = t->nodes[node].children[pos];
            stack_node[sp] = child;
            stack_pos[sp] = 0;
            sp++;
            t->euler[ei] = child;
            t->euler_depth[ei] = t->depth[child];
            if (t->first_visit[child] == GD_NODE_NONE) {
                t->first_visit[child] = ei;
            }
            ei++;
        } else {
            sp--;
            if (sp > 0) {
                /* Returning to parent: emit parent again. */
                size_t parent = stack_node[sp - 1];
                t->euler[ei] = parent;
                t->euler_depth[ei] = t->depth[parent];
                ei++;
            }
        }
    }

    gd_free(stack_node);
    gd_free(stack_pos);

    /* Build RMQ sparse table for euler_depth (minimum query). */
    size_t L = euler_len;
    /* Number of columns = floor(log2(L)) + 1. */
    size_t logL = 0;
    size_t tmp = L;
    while (tmp > 1) {
        tmp >>= 1;
        logL++;
    }
    size_t cols = logL + 1;
    t->rmq_cols = cols;
    t->rmq_rows = L;

    /* Check overflow: L * cols. */
    if (L > 0 && cols > SIZE_MAX / L) {
        gd_error_set(err, GD_E_NOMEM, "tree: RMQ table overflow");
        return GD_E_NOMEM;
    }
    t->rmq_table = (size_t *)gd_calloc_checked(L * cols, sizeof(size_t), err);
    if (t->rmq_table == NULL) {
        return GD_E_NOMEM;
    }

    /* RMQ stores the index (into euler[]) of the minimum euler_depth. */
    for (size_t i = 0; i < L; i++) {
        t->rmq_table[i * cols + 0] = i;
    }
    for (size_t j = 1; j < cols; j++) {
        size_t step = (size_t)1 << (j - 1);
        for (size_t i = 0; i + ((size_t)1 << j) <= L; i++) {
            size_t left = t->rmq_table[i * cols + (j - 1)];
            size_t right_idx = i + step;
            size_t right = t->rmq_table[right_idx * cols + (j - 1)];
            t->rmq_table[i * cols + j] =
                (t->euler_depth[left] <= t->euler_depth[right]) ? left : right;
        }
    }

    t->lca_built = true;
    return 0;
}

size_t gd_tree_mrca(const gd_tree *t, size_t a, size_t b)
{
    if (t == NULL || !t->lca_built || a >= t->nnodes || b >= t->nnodes) {
        return GD_NODE_NONE;
    }
    size_t fa = t->first_visit[a];
    size_t fb = t->first_visit[b];
    if (fa > fb) {
        size_t tmp = fa;
        fa = fb;
        fb = tmp;
    }
    /* RMQ on euler_depth[fa..fb]. */
    size_t len = fb - fa + 1;
    size_t k = 0;
    size_t tmp2 = len;
    while (tmp2 > 1) {
        tmp2 >>= 1;
        k++;
    }
    size_t step = (size_t)1 << k;
    size_t left = t->rmq_table[fa * t->rmq_cols + k];
    size_t right = t->rmq_table[(fb - step + 1) * t->rmq_cols + k];
    size_t min_idx = (t->euler_depth[left] <= t->euler_depth[right])
        ? left : right;
    return t->euler[min_idx];
}

int gd_tree_resolve_target(const gd_tree *t, const char *target, size_t *out,
                           gd_error *err)
{
    if (t == NULL || target == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "resolve_target: NULL");
    }
    if (strncmp(target, "mrca:", 5) != 0) {
        size_t idx = gd_tree_find_label(t, target);
        if (idx == GD_NODE_NONE) {
            return gd_error_setf(err, GD_E_CALIBRATION,
                                 "target '%s' is not a node label in the tree", target);
        }
        *out = idx;
        return 0;
    }
    size_t acc = GD_NODE_NONE;
    const char *p = target + 5;
    while (*p != '\0') {
        while (*p == ' ' || *p == '\t') p++;
        const char *e = p;
        while (*e != '\0' && *e != ',') e++;
        const char *ee = e;
        while (ee > p && (ee[-1] == ' ' || ee[-1] == '\t')) ee--;
        size_t len = (size_t)(ee - p);
        if (len == 0) {
            return gd_error_setf(err, GD_E_CALIBRATION,
                                 "target '%s': empty tip name in the mrca list", target);
        }
        char name[256];
        if (len >= sizeof name) {
            return gd_error_setf(err, GD_E_CALIBRATION,
                                 "target '%s': a tip name is too long", target);
        }
        memcpy(name, p, len);
        name[len] = '\0';
        size_t v = gd_tree_find_label(t, name);
        if (v == GD_NODE_NONE || !t->nodes[v].is_leaf) {
            return gd_error_setf(err, GD_E_CALIBRATION,
                                 "target '%s': '%s' is not a tip of the tree", target, name);
        }
        acc = (acc == GD_NODE_NONE) ? v : gd_tree_mrca(t, acc, v);
        if (acc == GD_NODE_NONE) {
            return gd_error_set(err, GD_E_CALIBRATION, "target: tree has no LCA support");
        }
        p = (*e == ',') ? e + 1 : e;
    }
    if (acc == GD_NODE_NONE) {
        return gd_error_setf(err, GD_E_CALIBRATION, "target '%s' lists no tips", target);
    }
    *out = acc;
    return 0;
}

int gd_tree_name_internal_nodes(gd_tree *t, gd_error *err)
{
    if (t == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "name_internal_nodes: NULL");
    const char **minlab = (const char **)gd_calloc_checked(t->nnodes ? t->nnodes : 1,
                                                           sizeof(char *), err);
    if (minlab == NULL) return GD_E_NOMEM;
    /* Smallest tip label in each clade, children before parents. */
    for (size_t i = 0; i < t->nnodes; i++) {
        size_t v = t->postorder[i];
        if (t->nodes[v].is_leaf) { minlab[v] = t->nodes[v].label; continue; }
        for (size_t c = 0; c < t->nodes[v].nchildren; c++) {
            const char *m = minlab[t->nodes[v].children[c]];
            if (m != NULL && (minlab[v] == NULL || strcmp(m, minlab[v]) < 0)) minlab[v] = m;
        }
    }
    int rc = 0;
    for (size_t v = 0; v < t->nnodes && rc == 0; v++) {
        gd_node *nd = &t->nodes[v];
        if (nd->is_leaf || nd->label != NULL || nd->nchildren < 2) continue;
        const char *a = minlab[nd->children[0]], *b = minlab[nd->children[1]];
        if (a == NULL || b == NULL) continue;   /* unlabelled tips: nothing stable to say */
        if (strcmp(a, b) > 0) { const char *tmp = a; a = b; b = tmp; }
        size_t need = 5 + strlen(a) + 1 + strlen(b) + 1;
        char *lab = (char *)gd_malloc_checked(need, err);
        if (lab == NULL) { rc = GD_E_NOMEM; break; }
        memcpy(lab, "mrca:", 5);
        memcpy(lab + 5, a, strlen(a));
        lab[5 + strlen(a)] = ',';
        memcpy(lab + 6 + strlen(a), b, strlen(b) + 1);
        nd->label = lab;
    }
    gd_free((void *)minlab);
    return rc;
}

size_t gd_tree_find_label(const gd_tree *t, const char *label)
{
    if (t == NULL || label == NULL) {
        return GD_NODE_NONE;
    }
    for (size_t i = 0; i < t->nnodes; i++) {
        if (t->nodes[i].label != NULL &&
            strcmp(t->nodes[i].label, label) == 0) {
            return i;
        }
    }
    return GD_NODE_NONE;
}

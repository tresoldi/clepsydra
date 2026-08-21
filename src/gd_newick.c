/*
 * gd_newick.c — Newick tree parser for Clepsydra.
 *
 * Parses rooted Newick with:
 *  - leaf and internal labels (before ':');
 *  - branch lengths (after ':');
 *  - polytomies (nonbinary, e.g. (a,b,c)d);
 *  - comments in [...] (simple, non-nested; stripped);
 *  - terminal ';'.
 *
 * The parser is iterative (no deep recursion) and builds the tree in a
 * single pass using a node stack. It does not guess the meaning of branch
 * lengths (that is configured by the caller).
 */
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_tree_build.h"
#include "gd_memory.h"
#include "gd_tsv.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Skip whitespace and comments. Returns the new position. */
static const char *skip_ws_comments(const char *p, const char *end)
{
    while (p < end) {
        if (isspace((unsigned char)*p)) {
            p++;
        } else if (*p == '[') {
            /* Skip comment; reject nested. */
            int depth = 0;
            while (p < end) {
                if (*p == '[') {
                    depth++;
                    if (depth > 1) {
                        return NULL;  /* nested comment: reject */
                    }
                    p++;
                } else if (*p == ']') {
                    depth--;
                    p++;
                    if (depth == 0) {
                        break;
                    }
                } else {
                    p++;
                }
            }
            if (depth != 0) {
                return NULL;  /* unterminated comment */
            }
        } else {
            break;
        }
    }
    return p;
}

/* Parse a label (until we hit ':', ',', ')', or ';'). Returns the label
 * start and sets *label_len. */
static const char *parse_label(const char *p, const char *end,
                               const char **label, size_t *label_len)
{
    *label = p;
    *label_len = 0;
    while (p < end) {
        char c = *p;
        if (c == ':' || c == ',' || c == ')' || c == ';' || c == '(' ||
            isspace((unsigned char)c) || c == '[') {
            break;
        }
        p++;
    }
    *label_len = (size_t)(p - *label);
    return p;
}

/* Parse a branch length (after ':'). Sets *blen. */
static const char *parse_branch_length(const char *p, const char *end,
                                       double *blen, gd_error *err)
{
    *blen = -1.0;
    p = skip_ws_comments(p, end);
    if (p == NULL) {
        gd_error_set(err, GD_E_PARSE, "newick: bad comment in branch length");
        return NULL;
    }
    /* Expect a number. */
    const char *start = p;
    while (p < end && *p != ',' && *p != ')' && *p != ';' &&
           !isspace((unsigned char)*p) && *p != '[') {
        p++;
    }
    size_t len = (size_t)(p - start);
    if (len == 0) {
        gd_error_set(err, GD_E_PARSE, "newick: empty branch length");
        return NULL;
    }
    int rc = gd_tsv_parse_double(start, len, blen);
    if (rc != 0) {
        gd_error_set(err, GD_E_BAD_NUMERIC, "newick: bad branch length");
        return NULL;
    }
    /* Branch lengths must be finite and non-negative. strtod (used by
     * gd_tsv_parse_double) accepts "nan"/"inf" and negative values, which are
     * meaningless as durations and would silently poison downstream models. */
    if (!isfinite(*blen) || *blen < 0.0) {
        gd_error_set(err, GD_E_BAD_NUMERIC,
                     "newick: branch length must be finite and non-negative");
        return NULL;
    }
    return p;
}

/* Reject duplicate node labels. Labels are identifiers in the dating workflow
 * (calibrations and forms bind to nodes by label), so a repeated tip or
 * internal label would silently bind to the wrong node. Unlabeled internal
 * nodes (label == NULL) are permitted and not compared. */
static int check_unique_labels(const gd_tree *t, gd_error *err)
{
    for (size_t i = 0; i < t->nnodes; i++) {
        const char *li = t->nodes[i].label;
        if (li == NULL) continue;
        for (size_t j = i + 1; j < t->nnodes; j++) {
            const char *lj = t->nodes[j].label;
            if (lj != NULL && strcmp(li, lj) == 0) {
                gd_error_set(err, GD_E_DUPLICATE_ID,
                             "newick: duplicate node label");
                return GD_E_DUPLICATE_ID;
            }
        }
    }
    return 0;
}

int gd_newick_parse(const char *in_buf, size_t in_len,
                    gd_tree **out, gd_error *err)
{
    char *buf = NULL;
    size_t buf_len = 0;

    if (in_buf == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "newick: no input");
    }
    /* Work on a copy so we can NUL-terminate safely. */
    if (in_len > SIZE_MAX - 1) {
        return gd_error_set(err, GD_E_NOMEM, "newick: buffer too large");
    }
    buf = (char *)gd_malloc_checked(in_len + 1, err);
    if (buf == NULL) {
        return GD_E_NOMEM;
    }
    memcpy(buf, in_buf, in_len);
    buf[in_len] = '\0';
    buf_len = in_len;

    /* First pass: count nodes to allocate. We count '(' and labels. */
    /* Actually, we'll build dynamically: start with a reasonable capacity
     * and grow. */
    size_t cap = 16;
    gd_tree *t = gd_tree_create(cap, err);
    if (t == NULL) {
        gd_free(buf);
        return GD_E_NOMEM;
    }
    size_t node_count = 0;

    /* Allocate a node-index stack for the parser. */
    size_t stack_cap = 64;
    size_t *stack = (size_t *)gd_malloc_checked(stack_cap * sizeof(size_t), err);
    if (stack == NULL) {
        gd_tree_destroy(t);
        gd_free(buf);
        return GD_E_NOMEM;
    }
    size_t sp = 0;

    /* Helper: allocate a new node, growing the nodes array if needed. */
    /* We pre-allocate cap nodes; if we exceed, we realloc. */
    #define ENSURE_CAP() do { \
        if (node_count >= cap) { \
            size_t new_cap = cap * 2; \
            gd_node *nn = (gd_node *)gd_realloc_array_checked(t->nodes, \
                new_cap, sizeof(gd_node), err); \
            if (nn == NULL) { \
                gd_error_set(err, GD_E_NOMEM, "newick: node alloc"); \
                rc = GD_E_NOMEM; goto done; \
            } \
            for (size_t _i = cap; _i < new_cap; _i++) { \
                nn[_i].parent = GD_NODE_NONE; \
                nn[_i].children = NULL; \
                nn[_i].nchildren = 0; \
                nn[_i].child_cap = 0; \
                nn[_i].label = NULL; \
                nn[_i].branch_length = -1.0; \
                nn[_i].is_leaf = false; \
                nn[_i].age = 0.0; \
            } \
            t->nodes = nn; \
            t->nnodes = new_cap; \
            cap = new_cap; \
        } \
    } while (0)

    int rc = 0;
    const char *p = buf;
    const char *end = buf + buf_len;
    int found_semicolon = 0;
    size_t root_idx = GD_NODE_NONE;

    /* Grammar state: 1 when the parser expects an operand (a subtree '(' or a
     * leaf label), 0 when it expects a separator (',' or ')') or the terminal
     * ';'. This rejects empty children ("(a,,b)", "(,a)", "()"), missing
     * separators ("(a b)"), and stray operands, none of which the plain token
     * loop caught. */
    int expect_operand = 1;

    p = skip_ws_comments(p, end);
    if (p == NULL) {
        gd_error_set(err, GD_E_PARSE, "newick: bad comment at start");
        rc = GD_E_PARSE;
        goto done;
    }

    while (p < end) {
        p = skip_ws_comments(p, end);
        if (p == NULL) {
            gd_error_set(err, GD_E_PARSE, "newick: bad comment");
            rc = GD_E_PARSE;
            goto done;
        }
        if (p >= end) {
            break;
        }
        char c = *p;

        if (c == '(') {
            if (!expect_operand) {
                gd_error_set(err, GD_E_PARSE,
                             "newick: missing separator before '('");
                rc = GD_E_PARSE;
                goto done;
            }
            /* Start a new internal node. Push a placeholder on the stack;
             * we'll assign its index when we encounter its first child. */
            if (sp >= stack_cap) {
                size_t new_sc = stack_cap * 2;
                size_t *ns = (size_t *)gd_realloc_array_checked(stack, new_sc, sizeof(size_t), err);
                if (ns == NULL) {
                    gd_error_set(err, GD_E_NOMEM, "newick: stack overflow");
                    rc = GD_E_NOMEM;
                    goto done;
                }
                stack = ns;
                stack_cap = new_sc;
            }
            /* Create the internal node now. */
            ENSURE_CAP();
            size_t idx = node_count++;
            stack[sp++] = idx;
            p++;
        } else if (c == ',') {
            if (expect_operand) {
                gd_error_set(err, GD_E_PARSE,
                             "newick: empty child (misplaced ',')");
                rc = GD_E_PARSE;
                goto done;
            }
            expect_operand = 1;
            p++;
        } else if (c == ')') {
            /* Close current internal node. */
            if (sp == 0) {
                gd_error_set(err, GD_E_PARSE, "newick: unbalanced ')'");
                rc = GD_E_PARSE;
                goto done;
            }
            if (expect_operand) {
                gd_error_set(err, GD_E_PARSE,
                             "newick: empty node (no child before ')')");
                rc = GD_E_PARSE;
                goto done;
            }
            size_t idx = stack[--sp];
            p++;
            /* Parse the internal node's label and branch length. */
            p = skip_ws_comments(p, end);
            if (p == NULL) { rc = GD_E_PARSE; goto done; }
            if (p < end && *p != ',' && *p != ')' && *p != ';') {
                const char *label;
                size_t label_len;
                p = parse_label(p, end, &label, &label_len);
                if (label_len > 0) {
                    t->nodes[idx].label = gd_strdup_n(label, label_len, err);
                    if (t->nodes[idx].label == NULL) {
                        rc = GD_E_NOMEM;
                        goto done;
                    }
                }
            }
            p = skip_ws_comments(p, end);
            if (p == NULL) { rc = GD_E_PARSE; goto done; }
            if (p < end && *p == ':') {
                p++;
                double blen;
                p = parse_branch_length(p, end, &blen, err);
                if (p == NULL) { rc = GD_E_PARSE; goto done; }
                t->nodes[idx].branch_length = blen;
            }
            if (sp == 0) {
                root_idx = idx;
            }
            /* Link to parent. */
            if (sp > 0) {
                size_t parent = stack[sp - 1];
                rc = gd_tree_add_child(t, parent, idx, err);
                if (rc != 0) { goto done; }
            }
            /* A closed subtree is a completed operand. */
            expect_operand = 0;
        } else if (c == ';') {
            if (expect_operand) {
                gd_error_set(err, GD_E_PARSE,
                             "newick: missing operand before ';'");
                rc = GD_E_PARSE;
                goto done;
            }
            found_semicolon = 1;
            p++;
            break;
        } else {
            /* A leaf label. */
            if (!expect_operand) {
                gd_error_set(err, GD_E_PARSE,
                             "newick: missing separator before label");
                rc = GD_E_PARSE;
                goto done;
            }
            ENSURE_CAP();
            size_t idx = node_count++;
            const char *label;
            size_t label_len;
            p = parse_label(p, end, &label, &label_len);
            if (label_len > 0) {
                t->nodes[idx].label = gd_strdup_n(label, label_len, err);
                if (t->nodes[idx].label == NULL) {
                    rc = GD_E_NOMEM;
                    goto done;
                }
            }
            /* Branch length. */
            p = skip_ws_comments(p, end);
            if (p == NULL) { rc = GD_E_PARSE; goto done; }
            if (p < end && *p == ':') {
                p++;
                double blen;
                p = parse_branch_length(p, end, &blen, err);
                if (p == NULL) { rc = GD_E_PARSE; goto done; }
                t->nodes[idx].branch_length = blen;
            }
            /* Link to parent on the stack. */
            if (sp > 0) {
                size_t parent = stack[sp - 1];
                rc = gd_tree_add_child(t, parent, idx, err);
                if (rc != 0) { goto done; }
            } else {
                /* Single leaf tree. */
                root_idx = idx;
            }
            /* A leaf is a completed operand. */
            expect_operand = 0;
        }
    }

    /* After the terminal ';', only whitespace and comments may remain; any
     * further tree data (e.g. "(a,b)r;junk" or a second tree) is an error. */
    if (found_semicolon) {
        p = skip_ws_comments(p, end);
        if (p == NULL) {
            gd_error_set(err, GD_E_PARSE, "newick: bad comment after ';'");
            rc = GD_E_PARSE;
            goto done;
        }
        if (p < end) {
            gd_error_set(err, GD_E_PARSE, "newick: trailing data after ';'");
            rc = GD_E_PARSE;
            goto done;
        }
    }

    if (!found_semicolon) {
        gd_error_set(err, GD_E_PARSE, "newick: missing terminal ';'");
        rc = GD_E_PARSE;
        goto done;
    }
    if (root_idx == GD_NODE_NONE) {
        gd_error_set(err, GD_E_PARSE, "newick: no root node");
        rc = GD_E_PARSE;
        goto done;
    }
    if (sp != 0) {
        gd_error_set(err, GD_E_PARSE, "newick: unbalanced '('");
        rc = GD_E_PARSE;
        goto done;
    }

    t->root = root_idx;
    t->nnodes = node_count;

    /* Reject duplicate labels before any label-keyed lookups downstream. */
    rc = check_unique_labels(t, err);
    if (rc != 0) { goto done; }

    /* Build traversals. */
    rc = gd_tree_build_traversals(t, err);
    if (rc != 0) { goto done; }

    /* Build LCA. */
    rc = gd_tree_build_lca(t, err);
    if (rc != 0) { goto done; }

    *out = t;

done:
    if (rc != 0) {
        /* Ensure all allocated labels/children are freed: set nnodes to the
         * actual number of nodes created so far. */
        if (t != NULL) {
            t->nnodes = node_count;
        }
        gd_tree_destroy(t);
    }
    gd_free(stack);
    gd_free(buf);
    return rc;
}

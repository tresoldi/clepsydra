/*
 * test_newick.c — tests for the Newick parser and tree structure.
 */
#include "test_runner.h"
#include "gd_tree_build.h"
#include "gd_newick.h"
#include "gd_tree.h"

#include "gd.h"

#include <string.h>

TEST(test_newick_simple_binary)
{
    const char *nwk = "((italian,spanish)romance_west,romanian)romance_root;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    CHECK(t != NULL);
    /* 5 nodes: 3 leaves + 2 internal. */
    CHECK(t->nnodes == 5);

    /* Find the root. */
    size_t root = gd_tree_find_label(t, "romance_root");
    CHECK(root != GD_NODE_NONE);
    CHECK(t->root == root);

    /* Root has 2 children. */
    CHECK(t->nodes[root].nchildren == 2);

    /* Find romance_west; it has 2 children. */
    size_t rw = gd_tree_find_label(t, "romance_west");
    CHECK(rw != GD_NODE_NONE);
    CHECK(t->nodes[rw].nchildren == 2);

    /* Leaves. */
    size_t it = gd_tree_find_label(t, "italian");
    size_t sp = gd_tree_find_label(t, "spanish");
    size_t ro = gd_tree_find_label(t, "romanian");
    CHECK(it != GD_NODE_NONE);
    CHECK(sp != GD_NODE_NONE);
    CHECK(ro != GD_NODE_NONE);
    CHECK(gd_tree_is_leaf(t, it) == 1);
    CHECK(gd_tree_is_leaf(t, sp) == 1);
    CHECK(gd_tree_is_leaf(t, ro) == 1);

    /* MRCA of italian and spanish is romance_west. */
    size_t mrca = gd_tree_mrca(t, it, sp);
    CHECK(mrca == rw);

    /* MRCA of italian and romanian is the root. */
    mrca = gd_tree_mrca(t, it, ro);
    CHECK(mrca == root);

    gd_tree_destroy(t);
}

TEST(test_newick_branch_lengths)
{
    const char *nwk = "(a:1.5,b:2.5)r;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    size_t a = gd_tree_find_label(t, "a");
    size_t b = gd_tree_find_label(t, "b");
    CHECK_APPROX(t->nodes[a].branch_length, 1.5, 1e-10);
    CHECK_APPROX(t->nodes[b].branch_length, 2.5, 1e-10);
    /* Root has no branch length. */
    CHECK(t->nodes[t->root].branch_length == -1.0);
    gd_tree_destroy(t);
}

TEST(test_newick_polytomy)
{
    const char *nwk = "(a,b,c)d;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    size_t d = gd_tree_find_label(t, "d");
    CHECK(d != GD_NODE_NONE);
    CHECK(t->nodes[d].nchildren == 3);
    CHECK(gd_tree_is_leaf(t, d) == 0);
    gd_tree_destroy(t);
}

TEST(test_newick_missing_semicolon)
{
    const char *nwk = "(a,b)r";  /* no ; */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == GD_E_PARSE);
    CHECK(t == NULL);
}

TEST(test_newick_unbalanced)
{
    const char *nwk = "(a,b);";  /* missing root label is OK, but... */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    /* This should parse: (a,b) with the root unlabeled. */
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    CHECK(t->nnodes == 3);
    gd_tree_destroy(t);

    /* Unbalanced parens. */
    const char *bad = "((a,b);";
    rc = gd_newick_parse(bad, strlen(bad), &t, &err);
    CHECK(rc == GD_E_PARSE);
}

TEST(test_newick_single_leaf)
{
    const char *nwk = "latin;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    CHECK(t->nnodes == 1);
    CHECK(gd_tree_is_leaf(t, 0) == 1);
    gd_tree_destroy(t);
}

TEST(test_newick_traversal_order)
{
    /*        r
     *       / \
     *      x   d
     *     / \
     *    a   b
     */
    const char *nwk = "((a,b)x,d)r;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);

    /* Postorder: a, b, x, d, r (children before parent). */
    CHECK(t->postorder[0] == gd_tree_find_label(t, "a"));
    CHECK(t->postorder[1] == gd_tree_find_label(t, "b"));
    CHECK(t->postorder[2] == gd_tree_find_label(t, "x"));
    CHECK(t->postorder[3] == gd_tree_find_label(t, "d"));
    CHECK(t->postorder[4] == t->root);

    /* Preorder: r, x, a, b, d (parent before children). */
    CHECK(t->preorder[0] == t->root);
    CHECK(t->preorder[1] == gd_tree_find_label(t, "x"));
    CHECK(t->preorder[2] == gd_tree_find_label(t, "a"));
    CHECK(t->preorder[3] == gd_tree_find_label(t, "b"));
    CHECK(t->preorder[4] == gd_tree_find_label(t, "d"));

    gd_tree_destroy(t);
}

TEST(test_newick_depth)
{
    const char *nwk = "((a,b)x,d)r;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    CHECK(gd_tree_depth(t, t->root) == 0);
    CHECK(gd_tree_depth(t, gd_tree_find_label(t, "x")) == 1);
    CHECK(gd_tree_depth(t, gd_tree_find_label(t, "d")) == 1);
    CHECK(gd_tree_depth(t, gd_tree_find_label(t, "a")) == 2);
    gd_tree_destroy(t);
}

TEST(test_newick_comments_stripped)
{
    const char *nwk = "(a [leaf a],b)r;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    CHECK(rc == 0);
    size_t a = gd_tree_find_label(t, "a");
    CHECK(a != GD_NODE_NONE);
    gd_tree_destroy(t);
}

TEST(test_newick_nested_comment_rejected)
{
    const char *nwk = "(a[b[c]],b)r;";
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    /* Nested comments are rejected. */
    CHECK(rc == GD_E_PARSE);
}

TEST(test_newick_deep_tree_no_stack_overflow)
{
    /* A ladder of 10000 leaves: (((((...a)a1)a2...)aN;
     * Tests iterative traversal doesn't overflow the stack. */
    /* Build a long ladder string. */
    char buf[80002];
    int pos = 0;
    int depth = 5000;
    for (int i = 0; i < depth; i++) {
        buf[pos++] = '(';
    }
    buf[pos++] = 'a';
    for (int i = 0; i < depth; i++) {
        pos += (int)snprintf(buf + pos, sizeof(buf) - (size_t)pos, ")n%d", i);
    }
    buf[pos++] = ';';
    buf[pos] = '\0';

    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(buf, (size_t)pos, &t, &err);
    CHECK(rc == 0);
    CHECK(t->nnodes == (size_t)(depth + 1));
    /* Postorder should have a as first. */
    size_t a = gd_tree_find_label(t, "a");
    CHECK(t->postorder[0] == a);
    gd_tree_destroy(t);
}

/* --- strict-parsing regression tests (adversarial review) --- */

static int parse_fails(const char *nwk)
{
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse(nwk, strlen(nwk), &t, &err);
    if (rc == 0) { gd_tree_destroy(t); return 0; }
    return 1;
}

TEST(test_newick_reject_trailing_data)
{
    CHECK(parse_fails("((a,b)r);junk"));
    CHECK(parse_fails("(a,b)r;(c,d)s;"));
    /* Trailing whitespace/comments are still fine. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse("(a,b)r;  \n", 10, &t, &err);
    CHECK(rc == 0);
    gd_tree_destroy(t);
}

TEST(test_newick_reject_missing_separator)
{
    CHECK(parse_fails("(a b)r;"));
    CHECK(parse_fails("(a,b)r c;"));
}

TEST(test_newick_reject_empty_child)
{
    CHECK(parse_fails("(a,,b)r;"));
    CHECK(parse_fails("(,a)r;"));
    CHECK(parse_fails("(a,)r;"));
    CHECK(parse_fails("()r;"));
}

TEST(test_newick_reject_bad_branch_length)
{
    CHECK(parse_fails("(a,b)r:NaN;"));
    CHECK(parse_fails("(a:nan,b)r;"));
    CHECK(parse_fails("(a,b)r:-1;"));
    CHECK(parse_fails("(a:inf,b)r;"));
}

TEST(test_newick_reject_duplicate_labels)
{
    CHECK(parse_fails("(a,a)r;"));
    CHECK(parse_fails("((a,b)x,c)x;"));
}

TEST(test_newick_single_child_ok)
{
    /* Single-child (ladder) internal nodes remain valid. */
    gd_error err = GD_ERROR_INIT;
    gd_tree *t = NULL;
    int rc = gd_newick_parse("(a)b;", 5, &t, &err);
    CHECK(rc == 0);
    CHECK(t->nnodes == 2);
    gd_tree_destroy(t);
}

/*
 * Found by fuzzing (docs/fuzzing.md), 2026-08-21.
 *
 * "a,(b);" is not a tree: it has two top-level nodes and no enclosing
 * parentheses. The parser used to accept it and hand back a structure with a
 * second node whose parent was GD_NODE_NONE. Five places downstream read a
 * node's parent whenever the node is not the root, so each of them indexed an
 * array at SIZE_MAX — undefined behaviour, caught by UBSan inside
 * gd_calibrations_propagate.
 *
 * The check lives in gd_tree_build_traversals, which every consumer calls
 * after construction, so a hand-built tree is covered as well as a parsed one.
 */
TEST(test_newick_rejects_a_forest)
{
    static const char *const forests[] = {
        "a,(b);",           /* the reproducer libFuzzer found */
        "a,b;",
        "(a),(b);",
        "(a,b),c;",
    };
    for (size_t i = 0; i < sizeof forests / sizeof forests[0]; i++) {
        gd_error err = GD_ERROR_INIT;
        gd_tree *t = NULL;
        int rc = gd_newick_parse(forests[i], strlen(forests[i]), &t, &err);
        CHECK(rc != GD_E_OK);
        CHECK(t == NULL);
        gd_tree_destroy(t);
    }

    /* The same shapes with a root are still accepted, so the check rejects
     * forests rather than commas. */
    static const char *const trees[] = { "(a,b);", "((a),(b));", "((a,b),c);" };
    for (size_t i = 0; i < sizeof trees / sizeof trees[0]; i++) {
        gd_error err = GD_ERROR_INIT;
        gd_tree *t = NULL;
        CHECK(gd_newick_parse(trees[i], strlen(trees[i]), &t, &err) == GD_E_OK);
        CHECK(t != NULL);
        gd_tree_destroy(t);
    }
}

static void run_tests(void)
{
    RUN(test_newick_rejects_a_forest);
    RUN(test_newick_reject_trailing_data);
    RUN(test_newick_reject_missing_separator);
    RUN(test_newick_reject_empty_child);
    RUN(test_newick_reject_bad_branch_length);
    RUN(test_newick_reject_duplicate_labels);
    RUN(test_newick_single_child_ok);
    RUN(test_newick_simple_binary);
    RUN(test_newick_branch_lengths);
    RUN(test_newick_polytomy);
    RUN(test_newick_missing_semicolon);
    RUN(test_newick_unbalanced);
    RUN(test_newick_single_leaf);
    RUN(test_newick_traversal_order);
    RUN(test_newick_depth);
    RUN(test_newick_comments_stripped);
    RUN(test_newick_nested_comment_rejected);
    RUN(test_newick_deep_tree_no_stack_overflow);
}

TEST_MAIN()

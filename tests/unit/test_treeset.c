/*
 * test_treeset.c — reading tree samples (roadmap E3.4): NEXUS with a
 * Translate table and BEAST-style [&...] comments, and plain Newick.
 */
#include "test_runner.h"
#include "gd_treeset.h"

#include <stdlib.h>
#include <string.h>

TEST(test_treeset_nexus_translate_and_comments)
{
    static const char *NEX =
        "#NEXUS\n\nBegin taxa;\n  Dimensions ntax=3;\nEnd;\nBegin trees;\n"
        "  Translate\n     1 Alpha,\n     2 Beta,\n     3 'Gamma'\n;\n"
        "tree STATE_0 = [&R] ((1[&rate=1.0]:1.5,2:1.5)[&posterior=1]:2,3:3.5);\n"
        "tree STATE_10 = ((1:1,3:1):1,2:2);\n"
        "End;\n";
    gd_treeset *ts = NULL;
    gd_error err = GD_ERROR_INIT;
    int rc = gd_treeset_read(NEX, strlen(NEX), &ts, &err);
    size_t n = 0; char t0[128] = {0}, t1[128] = {0};
    if (rc == 0) {
        n = ts->ntrees;
        strncpy(t0, ts->newick[0], sizeof t0 - 1);
        strncpy(t1, ts->newick[1], sizeof t1 - 1);
        gd_treeset_destroy(ts);
    }
    CHECK(rc == 0);
    CHECK(n == 2);
    CHECK_STREQ(t0, "((Alpha:1.5,Beta:1.5):2,Gamma:3.5);");
    CHECK_STREQ(t1, "((Alpha:1,Gamma:1):1,Beta:2);");

    char *tips = NULL;
    CHECK(gd_newick_tip_list("((Alpha:1,Gamma:1):1,Beta:2);", &tips, &err) == 0);
    CHECK_STREQ(tips, "Alpha,Gamma,Beta");
    free(tips);
}

TEST(test_treeset_plain_newick)
{
    static const char *NWK = "((a,b),c);\n((a,c),b);\n\n";
    gd_treeset *ts = NULL;
    gd_error err = GD_ERROR_INIT;
    CHECK(gd_treeset_read(NWK, strlen(NWK), &ts, &err) == 0);
    size_t n = ts->ntrees;
    int same = (strcmp(ts->newick[1], "((a,c),b);") == 0);
    gd_treeset_destroy(ts);
    CHECK(n == 2);
    CHECK(same);
    ts = NULL;
    CHECK(gd_treeset_read("   \n", 4, &ts, &err) != 0);
    CHECK(ts == NULL);
}

static void run_tests(void)
{
    RUN(test_treeset_nexus_translate_and_comments);
    RUN(test_treeset_plain_newick);
}

TEST_MAIN()

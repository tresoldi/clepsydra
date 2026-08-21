/*
 * test_runner.h — minimal test framework for Clepsydra unit tests.
 */
#ifndef TEST_RUNNER_H
#define TEST_RUNNER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int g_test_total = 0;
static int g_test_failed = 0;
static const char *g_test_current = NULL;

/*
 * CHECK returns from the test on failure, so anything the test was holding --
 * a tree, a result, a gd_forms -- is leaked at that point. Under ASan the
 * leak report then runs to dozens of stack traces and buries the one line
 * that says which assertion fired. That cost three separate diagnoses during
 * the D-series before it was recognised for what it is.
 *
 * The leak is a consequence of the failure, not an independent finding. So
 * leak detection stays fully on for a passing run -- where a leak IS a
 * finding, and several real ones have been caught that way -- and is turned
 * off only once a test has already failed and has something more useful to
 * report. LeakSanitizer queries this hook when it runs its check at exit.
 */
#if defined(__has_feature)
# if __has_feature(address_sanitizer)
#  define TEST_UNDER_ASAN 1
# endif
#elif defined(__SANITIZE_ADDRESS__)
# define TEST_UNDER_ASAN 1
#endif

#ifdef TEST_UNDER_ASAN
int __lsan_is_turned_off(void);
int __lsan_is_turned_off(void) { return (g_test_failed > 0) ? 1 : 0; }
#endif

#define TEST(name) \
    static void name(void); \
    static void name(void)

#define RUN(name) do { \
    g_test_current = #name; \
    g_test_total++; \
    name(); \
} while (0)

#define CHECK(cond) do { \
    if (!(cond)) { \
        printf("  FAIL %s:%d: %s\n", g_test_current ? g_test_current : "?", \
               __LINE__, #cond); \
        g_test_failed++; \
        return; \
    } \
} while (0)

#define CHECK_APPROX(a, b, tol) do { \
    if (fabs((double)(a) - (double)(b)) > (tol)) { \
        printf("  FAIL %s:%d: %s != %s (got %g, tol %g)\n", \
               g_test_current ? g_test_current : "?", __LINE__, #a, #b, \
               (double)(a), (double)(tol)); \
        g_test_failed++; \
        return; \
    } \
} while (0)

#define CHECK_STREQ(a, b) do { \
    if (strcmp((a), (b)) != 0) { \
        printf("  FAIL %s:%d: \"%s\" != \"%s\"\n", \
               g_test_current ? g_test_current : "?", __LINE__, (a), (b)); \
        g_test_failed++; \
        return; \
    } \
} while (0)

#define TEST_MAIN() \
    int main(void) { \
        run_tests(); \
        printf("%d/%d tests passed\n", g_test_total - g_test_failed, g_test_total); \
        return (g_test_failed > 0) ? 1 : 0; \
    }

#endif /* TEST_RUNNER_H */

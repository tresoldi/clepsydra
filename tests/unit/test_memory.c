/*
 * test_memory.c — allocation-size arithmetic.
 *
 * The library allocates arrays whose element counts come from parsed input:
 * concepts, tree nodes, chains, draws. `malloc(n * sizeof(T))` on such an `n`
 * is the shape that wraps silently and hands back a buffer far smaller than
 * the loop that fills it believes — and unlike an ordinary out-of-memory, the
 * program does not stop, it corrupts. These tests pin the checks that make
 * that outcome impossible.
 *
 * They are cheap: every case here is rejected by arithmetic before any
 * allocator is asked for anything, so nothing tries to reserve an absurd
 * amount of memory.
 */
#include "test_runner.h"
#include "gd.h"
#include "gd_memory.h"

#include <stdint.h>

/* ------------------------------------------------------------ gd_size_mul */

TEST(test_size_mul_computes_the_product)
{
    size_t out = 12345;
    CHECK(gd_size_mul(0, 0, &out));       CHECK(out == 0);
    CHECK(gd_size_mul(0, 999, &out));     CHECK(out == 0);
    CHECK(gd_size_mul(999, 0, &out));     CHECK(out == 0);
    CHECK(gd_size_mul(1, 1, &out));       CHECK(out == 1);
    CHECK(gd_size_mul(6, 7, &out));       CHECK(out == 42);
}

TEST(test_size_mul_rejects_a_wrapping_product)
{
    size_t out = 12345;

    /* The exact boundary: SIZE_MAX/2 * 2 fits (or is one short of fitting),
     * one more element does not. */
    CHECK(gd_size_mul(SIZE_MAX, 1, &out));
    CHECK(out == SIZE_MAX);

    CHECK(!gd_size_mul(SIZE_MAX, 2, &out));
    CHECK(!gd_size_mul(2, SIZE_MAX, &out));
    CHECK(!gd_size_mul(SIZE_MAX / 2 + 1, 3, &out));

    /* A failed multiply leaves the output alone rather than writing a
     * half-computed size a careless caller might still use. */
    out = 4242;
    CHECK(!gd_size_mul(SIZE_MAX, 4, &out));
    CHECK(out == 4242);
}

/* ------------------------------------------------- checked array allocation */

TEST(test_array_alloc_rejects_a_wrapping_size)
{
    /* This is the case the raw form got wrong: n * sizeof(double) wraps to a
     * small number, malloc succeeds, and the caller writes past the end.
     * Here it is a clean NULL with GD_E_NOMEM and nothing is allocated. */
    gd_error err = GD_ERROR_INIT;
    /* SIZE_MAX/2 + 1 wraps against any element wider than one byte, so every
     * case below is rejected by the arithmetic and no allocator is ever asked
     * for anything. Picking a count that merely *is* enormous would instead
     * reach malloc, which is a different failure with a different owner. */
    size_t absurd = SIZE_MAX / 2 + 1;

    CHECK(gd_malloc_array_checked(absurd, sizeof(double), &err) == NULL);
    CHECK(err.code == GD_E_NOMEM);
    CHECK(err.message[0] != '\0');

    gd_error_clear(&err);
    CHECK(gd_calloc_checked(absurd, sizeof(double), &err) == NULL);
    CHECK(err.code == GD_E_NOMEM);

    gd_error_clear(&err);
    CHECK(gd_realloc_array_checked(NULL, absurd, sizeof(double), &err) == NULL);
    CHECK(err.code == GD_E_NOMEM);
}

TEST(test_array_alloc_reports_through_a_null_error)
{
    /* Several callers report a sentinel rather than a status and pass NULL.
     * The NULL return must still be the answer. */
    CHECK(gd_malloc_array_checked(SIZE_MAX / 2 + 1, sizeof(double), NULL) == NULL);
    CHECK(gd_calloc_checked(SIZE_MAX / 2 + 1, sizeof(double), NULL) == NULL);
}

TEST(test_array_alloc_succeeds_for_a_sane_request)
{
    gd_error err = GD_ERROR_INIT;
    double *p = (double *)gd_malloc_array_checked(16, sizeof(double), &err);

    CHECK(p != NULL);
    CHECK(err.code == GD_E_OK);
    for (size_t i = 0; i < 16; i++) p[i] = (double)i;   /* ASan checks the extent */
    CHECK_APPROX(p[15], 15.0, 1e-12);
    gd_free(p);
}

TEST(test_calloc_actually_zeroes)
{
    gd_error err = GD_ERROR_INIT;
    int *p = (int *)gd_calloc_checked(32, sizeof(int), &err);

    CHECK(p != NULL);
    for (size_t i = 0; i < 32; i++) CHECK(p[i] == 0);
    gd_free(p);
}

TEST(test_realloc_keeps_the_old_block_when_it_fails)
{
    /* realloc's contract, which the checked wrapper must not break: on
     * failure the original is still valid and still the caller's to free.
     * Getting this wrong leaks on every growth failure. */
    gd_error err = GD_ERROR_INIT;
    int *p = (int *)gd_malloc_array_checked(4, sizeof(int), &err);
    CHECK(p != NULL);
    p[0] = 7; p[3] = 9;

    int *grown = (int *)gd_realloc_array_checked(p, SIZE_MAX / 2 + 1, sizeof(int), &err);
    CHECK(grown == NULL);
    /* p is untouched and still readable. */
    CHECK(p[0] == 7);
    CHECK(p[3] == 9);
    gd_free(p);
}

TEST(test_realloc_grows_and_preserves_contents)
{
    gd_error err = GD_ERROR_INIT;
    int *p = (int *)gd_malloc_array_checked(4, sizeof(int), &err);
    CHECK(p != NULL);
    for (int i = 0; i < 4; i++) p[i] = i * 11;

    int *grown = (int *)gd_realloc_array_checked(p, 64, sizeof(int), &err);
    CHECK(grown != NULL);
    for (int i = 0; i < 4; i++) CHECK(grown[i] == i * 11);
    grown[63] = 1;                    /* ASan checks the new extent */
    gd_free(grown);
}

TEST(test_zero_sized_requests_are_not_failures)
{
    /* An empty dataset is a legitimate input, not an allocation error, so a
     * zero count must not be reported as one whichever way the allocator
     * chooses to answer it. */
    gd_error err = GD_ERROR_INIT;
    void *p = gd_malloc_array_checked(0, sizeof(double), &err);
    CHECK(err.code == GD_E_OK);
    gd_free(p);

    p = gd_calloc_checked(0, sizeof(double), &err);
    CHECK(err.code == GD_E_OK);
    gd_free(p);
}

/* -------------------------------------------------------------- strdup_n */

TEST(test_strdup_n_copies_exactly_len_bytes)
{
    gd_error err = GD_ERROR_INIT;
    char *s = gd_strdup_n("abcdef", 3, &err);

    CHECK(s != NULL);
    CHECK_STREQ(s, "abc");
    gd_free(s);

    CHECK(gd_strdup_n(NULL, 3, &err) == NULL);
}

static void run_tests(void)
{
    RUN(test_size_mul_computes_the_product);
    RUN(test_size_mul_rejects_a_wrapping_product);
    RUN(test_array_alloc_rejects_a_wrapping_size);
    RUN(test_array_alloc_reports_through_a_null_error);
    RUN(test_array_alloc_succeeds_for_a_sane_request);
    RUN(test_calloc_actually_zeroes);
    RUN(test_realloc_keeps_the_old_block_when_it_fails);
    RUN(test_realloc_grows_and_preserves_contents);
    RUN(test_zero_sized_requests_are_not_failures);
    RUN(test_strdup_n_copies_exactly_len_bytes);
}

TEST_MAIN()

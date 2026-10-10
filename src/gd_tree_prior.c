/*
 * gd_tree_prior.c — the uniform-conditional node-age prior and the lexical
 * root. See gd_tree_prior.h for the model.
 *
 * The volume recursion. For a free internal node v, let F_v(x) be the volume
 * of the free ages strictly below v given that v sits at age x. For a child c:
 *   c fixed (a tip, or pinned)  ->  contributes the factor 1 (its own subtree
 *                                   volume is a constant), and its age is a
 *                                   floor for v;
 *   c free                      ->  contributes G_c(x) = int_{floor_c}^x F_c(y) dy.
 * F_v = prod_c G_c, and floor_v = the largest floor among the children. On
 * x >= floor_v every factor is a polynomial in x, so F_v is one polynomial
 * there, of degree equal to the number of free internal nodes below v. The
 * root's F is V.
 *
 * Polynomials are kept in the scaled variable z = x / scale, which keeps the
 * powers in range on deep trees; constant factors are dropped throughout,
 * since only ratios V(x') / V(x) enter the sampler.
 */
#include "gd_tree_prior.h"
#include "gd_memory.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double *c; size_t n; } poly;   /* n coefficients, degree n-1 */

static double poly_eval(const poly *p, double z)
{
    double s = 0.0;
    for (size_t k = p->n; k-- > 0; ) s = s * z + p->c[k];
    return s;
}

/* out = a * b. */
static int poly_mul(const poly *a, const poly *b, poly *out)
{
    size_t n = a->n + b->n - 1;
    double *c = (double *)gd_calloc_checked(n, sizeof(double), NULL);
    if (c == NULL) return GD_E_NOMEM;
    for (size_t i = 0; i < a->n; i++)
        for (size_t j = 0; j < b->n; j++) c[i + j] += a->c[i] * b->c[j];
    out->c = c; out->n = n;
    return 0;
}

/* out(z) = int_{z0}^{z} p(t) dt. */
static int poly_integral_from(const poly *p, double z0, poly *out)
{
    size_t n = p->n + 1;
    double *c = (double *)gd_calloc_checked(n, sizeof(double), NULL);
    if (c == NULL) return GD_E_NOMEM;
    for (size_t k = 0; k < p->n; k++) c[k + 1] = p->c[k] / (double)(k + 1);
    poly q = { c, n };
    c[0] = -poly_eval(&q, z0);
    out->c = c; out->n = n;
    return 0;
}

gd_tree_prior *gd_tree_prior_create(const gd_tree *tree, const bool *is_free,
                                    const double *fixed_age,
                                    double lex_logmean, double lex_logsd,
                                    gd_error *err)
{
    if (tree == NULL || is_free == NULL || fixed_age == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "tree_prior: NULL argument");
        return NULL;
    }
    if (!(lex_logsd > 0.0) || !isfinite(lex_logsd) || !isfinite(lex_logmean)) {
        gd_error_set(err, GD_E_RANGE, "tree_prior: lexical root log-sd must be finite and > 0");
        return NULL;
    }
    size_t n = tree->nnodes;
    gd_tree_prior *tp = (gd_tree_prior *)gd_calloc_checked(1, sizeof *tp, err);
    poly *F = (poly *)gd_calloc_checked(n ? n : 1, sizeof(poly), err);
    double *fl = (double *)gd_calloc_checked(n ? n : 1, sizeof(double), err);
    if (tp == NULL || F == NULL || fl == NULL) {
        gd_error_set(err, GD_E_NOMEM, "tree_prior: alloc");
        free(tp); free(F); free(fl);
        return NULL;
    }
    tp->root = tree->root;
    tp->root_free = is_free[tree->root];
    tp->lex_logmean = lex_logmean;
    tp->lex_logsd = lex_logsd;

    /* Scale: the oldest fixed age, or 1000 years when every tip is modern. */
    double maxfixed = 0.0;
    for (size_t v = 0; v < n; v++)
        if (!is_free[v] && fixed_age[v] > maxfixed) maxfixed = fixed_age[v];
    tp->scale = (maxfixed > 1000.0) ? maxfixed : 1000.0;

    int rc = 0;
    for (size_t i = 0; i < n && rc == 0; i++) {
        size_t v = tree->postorder[i];
        if (!is_free[v]) continue;
        double one = 1.0;
        poly acc = { NULL, 0 };
        double *c0 = (double *)gd_malloc_checked(sizeof(double), NULL);
        if (c0 == NULL) { rc = GD_E_NOMEM; break; }
        c0[0] = one; acc.c = c0; acc.n = 1;
        double floor_v = 0.0;
        for (size_t k = 0; k < tree->nodes[v].nchildren && rc == 0; k++) {
            size_t c = tree->nodes[v].children[k];
            if (!is_free[c]) {
                if (fixed_age[c] > floor_v) floor_v = fixed_age[c];
                continue;
            }
            if (fl[c] > floor_v) floor_v = fl[c];
            poly G, prod;
            if (poly_integral_from(&F[c], fl[c] / tp->scale, &G) != 0) { rc = GD_E_NOMEM; break; }
            if (poly_mul(&acc, &G, &prod) != 0) { free(G.c); rc = GD_E_NOMEM; break; }
            free(G.c); free(acc.c);
            acc = prod;
        }
        if (rc != 0) { free(acc.c); break; }
        F[v] = acc;
        fl[v] = floor_v;
    }
    if (rc == 0 && tp->root_free) {
        tp->coef = F[tree->root].c;
        tp->ncoef = F[tree->root].n;
        F[tree->root].c = NULL;
        tp->floor = fl[tree->root];
    }
    for (size_t v = 0; v < n; v++) free(F[v].c);
    free(F); free(fl);
    if (rc != 0) {
        gd_error_set(err, GD_E_NOMEM, "tree_prior: polynomial alloc");
        gd_tree_prior_destroy(tp);
        return NULL;
    }
    return tp;
}

void gd_tree_prior_destroy(gd_tree_prior *tp)
{
    if (tp == NULL) return;
    free(tp->coef);
    free(tp);
}

double gd_tree_prior_log_volume(const gd_tree_prior *tp, double root_age)
{
    if (tp == NULL || !tp->root_free) return 0.0;
    if (!(root_age > tp->floor) || !isfinite(root_age)) return -HUGE_VAL;
    poly p = { tp->coef, tp->ncoef };
    double v = poly_eval(&p, root_age / tp->scale);
    /* Cancellation can leave a tiny non-positive value just above the floor,
     * where the true volume is a positive power of (x - floor). */
    if (!(v > 0.0)) return -HUGE_VAL;
    return log(v);
}

double gd_tree_prior_logpdf(const gd_tree_prior *tp, double root_age, double log_rate)
{
    if (tp == NULL) return 0.0;
    if (!(root_age > 0.0) || !isfinite(root_age)) return -HUGE_VAL;
    double ll = log_rate + log(root_age);
    double z = (ll - tp->lex_logmean) / tp->lex_logsd;
    double lp = -0.5 * z * z - log(tp->lex_logsd) - 0.5 * log(2.0 * M_PI) - log(root_age);
    if (tp->root_free) lp -= gd_tree_prior_log_volume(tp, root_age);
    return lp;
}

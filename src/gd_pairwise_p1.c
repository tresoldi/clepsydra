/*
 * gd_pairwise_p1.c — hierarchical pairwise (P1) engine for Clepsydra (M4).
 *
 * Extends P0 with the log-rate hierarchy: log r_{f,c} = μ_c + β_f + u_{f,c}.
 * For a single-family pair, β_f is a scalar and u_{f,c} is absorbed into the
 * per-concept reference prior SD. The posterior over T is computed by grid
 * integration, with per-concept rates marginalized via Gauss-Legendre
 * quadrature in log-rate space.
 *
 * Sensitivity presets (tight/standard/wide) scale the hyperprior SDs.
 */
#include "gd_pairwise.h"
#include "gd_grid.h"
#include "gd_memory.h"
#include "gd.h"
#include "gd_math.h"
#include "gd_rng.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Gauss-Hermite nodes used to marginalize the shared family clock offset β. */
#define GD_P1_BETA_NODES 7

/* ---- Gauss-Legendre quadrature nodes/weights on [-1,1] ---- */
/* We use a fixed 20-point rule (sufficient for smooth lognormal integrands). */

static void gauss_legendre_nodes(double *nodes, double *weights, size_t n)
{
    /* Simple Newton-Raphson computation of GL nodes (for small n). */
    for (size_t i = 0; i < n; i++) {
        /* Initial guess: Chebyshev point. */
        double x = cos(M_PI * ((double)(i + 1) - 0.25) / ((double)n + 0.5));
        double dp;
        for (int iter = 0; iter < 100; iter++) {
            /* Legendre polynomial P_n(x) and derivative P'_n(x). */
            double p0 = 1.0, p1 = x;
            for (size_t k = 1; k < n; k++) {
                double p2 = ((2.0 * (double)k + 1.0) * x * p1 - (double)k * p0) / ((double)k + 1.0);
                p0 = p1;
                p1 = p2;
            }
            dp = (double)n * (x * p1 - p0) / (x * x - 1.0);
            double dx = p1 / dp;
            x -= dx;
            if (fabs(dx) < 1e-15) break;
        }
        nodes[i] = x;
        weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
}

/* ---- Gauss-Hermite quadrature (for integrating against a Gaussian) ---- */
/* Nodes/weights for  ∫ e^{-x^2} f(x) dx ≈ Σ w_i f(x_i)  (Σ w_i = √π).
 * Newton iteration on the normalized Hermite recurrence (Numerical Recipes
 * gauher). Used to marginalize the shared family clock offset β analytically
 * rather than conditioning on a point value. Requires n >= 2. */
static void gauss_hermite_nodes(double *x, double *w, size_t n)
{
    const double EPS = 1e-14;
    const double PIM4 = 0.7511255444649425;   /* π^{-1/4} */
    const int MAXIT = 100;
    size_t m = (n + 1) / 2;
    double z = 0.0, pp = 0.0;
    for (size_t i = 0; i < m; i++) {
        double nn = (double)n;
        if (i == 0)
            z = sqrt(2.0 * nn + 1.0) - 1.85575 * pow(2.0 * nn + 1.0, -0.16667);
        else if (i == 1)
            z -= 1.14 * pow(nn, 0.426) / z;
        else if (i == 2)
            z = 1.86 * z - 0.86 * x[0];
        else if (i == 3)
            z = 1.91 * z - 0.91 * x[1];
        else
            z = 2.0 * z - x[i - 2];
        double p1 = PIM4;
        for (int its = 0; its < MAXIT; its++) {
            p1 = PIM4;
            double p2 = 0.0, p3;
            for (size_t j = 1; j <= n; j++) {
                double jj = (double)j;
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0 / jj) * p2 - sqrt((jj - 1.0) / jj) * p3;
            }
            pp = sqrt(2.0 * nn) * p2;
            double z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) <= EPS) break;
        }
        x[i] = z;
        x[n - 1 - i] = -z;
        w[i] = 2.0 / (pp * pp);
        w[n - 1 - i] = w[i];
    }
}

/* Transform [-1,1] node to [lo,hi]. */
static double transform_node(double node, double lo, double hi)
{
    return 0.5 * (lo + hi) + 0.5 * (hi - lo) * node;
}

static double transform_weight(double weight, double lo, double hi)
{
    return 0.5 * (hi - lo) * weight;
}

double gd_p1_loglik_concept_marginal(double d, int match, double y_c,
                                     double log_rate_mean, double log_rate_sd,
                                     double eps, double q_c,
                                     const double *quad_nodes,
                                     const double *quad_weights,
                                     size_t n_nodes)
{
    if (d < 0.0) return -HUGE_VAL;
    if (d == 0.0) {
        /* P(match) = 1 for any rate. */
        double p = (match) ? 1.0 : 0.0;
        if (eps > 0.0) {
            double p_match = (1.0 - eps) * p + eps * q_c;
            double p_obs = (match) ? p_match : (1.0 - p_match);
            double L = (y_c >= 1.0) ? p_obs : ((y_c <= 0.0) ? (1.0 - p_obs) : (y_c * p_obs + (1.0 - y_c) * (1.0 - p_obs)));
            return (L > 0) ? log(L) : -HUGE_VAL;
        }
        double p_obs = (match) ? 1.0 : 0.0;
        double L = (y_c >= 1.0) ? p_obs : ((y_c <= 0.0) ? (1.0 - p_obs) : (y_c * p_obs + (1.0 - y_c) * (1.0 - p_obs)));
        return (L > 0) ? log(L) : -HUGE_VAL;
    }

    /* Marginalize over log r ~ N(log_rate_mean, log_rate_sd^2).
     * The integrand is L_c(r) * phi(log r; mean, sd).
     * We integrate in log-rate space over [mean - 5*sd, mean + 5*sd]. */
    double lo = log_rate_mean - 5.0 * log_rate_sd;
    double hi = log_rate_mean + 5.0 * log_rate_sd;

    double total = 0.0;  /* linear sum (not log) for the integral */
    for (size_t i = 0; i < n_nodes; i++) {
        double lr = transform_node(quad_nodes[i], lo, hi);
        double w = transform_weight(quad_weights[i], lo, hi);
        double r = exp(lr);
        /* L_c(r) — the P0 per-concept likelihood (not log). */
        double s = exp(-r * d);
        double p_match = (eps > 0.0) ? ((1.0 - eps) * s + eps * q_c) : s;
        double p_obs = (match) ? p_match : (1.0 - p_match);
        double L;
        if (y_c >= 1.0) L = p_obs;
        else if (y_c <= 0.0) L = 1.0 - p_obs;
        else L = y_c * p_obs + (1.0 - y_c) * (1.0 - p_obs);
        /* Normal density in log-rate space. */
        double z = (lr - log_rate_mean) / log_rate_sd;
        double phi = exp(-0.5 * z * z) / (log_rate_sd * sqrt(2.0 * M_PI));
        total += L * phi * w;
    }
    return (total > 0.0) ? log(total) : -HUGE_VAL;
}

/* ---- P1 posterior ---- */

int gd_p1_posterior(const gd_concept_ref *refs, size_t n_concepts,
                    const int *matches, const double *y_c,
                    double a1, double a2,
                    double T_lo, double T_hi,
                    double beta_mean, double beta_sd,
                    double eps, double eps_a, double eps_b,
                    gd_sensitivity sens,
                    size_t grid_T_n, size_t grid_rate_n,
                    gd_pair_result *result, gd_error *err)
{
    (void)eps_a; (void)eps_b;  /* future: integrate over eps (ε fixed for now) */
    if (refs == NULL || matches == NULL || result == NULL || n_concepts == 0) {
        gd_error_set(err, GD_E_INVALID_ARG, "p1_posterior: bad args");
        return GD_E_INVALID_ARG;
    }
    if (T_lo < a1 || T_lo < a2) T_lo = (a1 > a2) ? a1 : a2;
    if (T_hi <= T_lo) {
        gd_error_set(err, GD_E_RANGE, "p1_posterior: T_hi <= T_lo");
        return GD_E_RANGE;
    }
    if (grid_T_n < 16) grid_T_n = 256;
    if (grid_rate_n < 10) grid_rate_n = 20;

    /* Sensitivity presets: scale the hyperprior SDs. */
    double sd_scale;
    switch (sens) {
        case GD_SENSITIVITY_TIGHT:    sd_scale = 0.5; break;
        case GD_SENSITIVITY_STANDARD: sd_scale = 1.0; break;
        case GD_SENSITIVITY_WIDE:     sd_scale = 2.0; break;
        default:                      sd_scale = 1.0; break;
    }

    /* Gauss-Legendre quadrature nodes/weights for rate marginalization. */
    double *gl_nodes = (double *)gd_malloc_array_checked(grid_rate_n, sizeof(double), err);
    double *gl_weights = (double *)gd_malloc_array_checked(grid_rate_n, sizeof(double), err);
    double *log_post = (double *)gd_malloc_array_checked(grid_T_n, sizeof(double), err);
    double *grid_T = (double *)gd_malloc_array_checked(grid_T_n, sizeof(double), err);
    if (gl_nodes == NULL || gl_weights == NULL || log_post == NULL || grid_T == NULL) {
        free(gl_nodes); free(gl_weights); free(log_post); free(grid_T);
        gd_error_set(err, GD_E_NOMEM, "p1: alloc");
        return GD_E_NOMEM;
    }
    gauss_legendre_nodes(gl_nodes, gl_weights, grid_rate_n);

    /* Grid in T. */
    double dT = (T_hi - T_lo) / (double)(grid_T_n - 1);
    for (size_t i = 0; i < grid_T_n; i++) {
        grid_T[i] = T_lo + (double)i * dT;
    }

    /* Family clock offset β ~ N(beta_mean, beta_sd²), shared across all
     * concepts. We integrate it out with Gauss-Hermite quadrature rather than
     * conditioning on a point value: b_nodes[k] are the β samples and
     * b_logw[k] = log(w_k/√π) are the log-weights of the Gaussian (Σ exp = 1).
     * When beta_sd == 0 this collapses to a single node at beta_mean, exactly
     * reproducing the point-β behaviour. */
    size_t nb;
    double b_nodes[GD_P1_BETA_NODES];
    double b_logw[GD_P1_BETA_NODES];
    if (beta_sd > 0.0) {
        double hx[GD_P1_BETA_NODES], hw[GD_P1_BETA_NODES];
        nb = GD_P1_BETA_NODES;
        gauss_hermite_nodes(hx, hw, nb);
        double log_sqrt_pi = 0.5 * log(M_PI);
        for (size_t k = 0; k < nb; k++) {
            b_nodes[k] = beta_mean + sqrt(2.0) * beta_sd * hx[k];
            b_logw[k] = log(hw[k]) - log_sqrt_pi;
        }
    } else {
        nb = 1;
        b_nodes[0] = beta_mean;
        b_logw[0] = 0.0;
    }

    /* For each T, the marginal log-likelihood integrates β out (outer
     * logsumexp over the β nodes) of the per-concept marginals (inner sum). */
    for (size_t i = 0; i < grid_T_n; i++) {
        double T = grid_T[i];
        double d = 2.0 * T - a1 - a2;
        if (d < 0.0) { log_post[i] = -HUGE_VAL; continue; }
        double acc = -HUGE_VAL;
        for (size_t k = 0; k < nb; k++) {
            double ll = b_logw[k];
            for (size_t c = 0; c < n_concepts; c++) {
                double lr_mean = refs[c].log_rate_mean + b_nodes[k];
                double lr_sd = refs[c].log_rate_sd * sd_scale;
                double y = (y_c != NULL) ? y_c[c] : 1.0;
                double lc = gd_p1_loglik_concept_marginal(
                    d, matches[c], y, lr_mean, lr_sd, eps, 0.5,
                    gl_nodes, gl_weights, grid_rate_n);
                if (lc == -HUGE_VAL) { ll = -HUGE_VAL; break; }
                ll += lc;
            }
            acc = gd_logaddexp(acc, ll);
        }
        log_post[i] = acc;
    }

    /* Normalize as a density (trapezoidal integral = 1). */
    double m = -HUGE_VAL;
    for (size_t i = 0; i < grid_T_n; i++)
        if (log_post[i] > m) m = log_post[i];
    if (m == -HUGE_VAL) {
        free(gl_nodes); free(gl_weights); free(log_post); free(grid_T);
        gd_error_set(err, GD_E_RANGE, "p1: all-zero likelihood");
        return GD_E_RANGE;
    }
    for (size_t i = 0; i < grid_T_n; i++) log_post[i] -= m;
    double integral = 0.0;
    for (size_t i = 1; i < grid_T_n; i++) {
        integral += 0.5 * (exp(log_post[i - 1]) + exp(log_post[i])) * dT;
    }
    if (integral <= 0.0) {
        free(gl_nodes); free(gl_weights); free(log_post); free(grid_T);
        gd_error_set(err, GD_E_RANGE, "p1: zero integral");
        return GD_E_RANGE;
    }
    double log_int = log(integral);
    for (size_t i = 0; i < grid_T_n; i++) log_post[i] -= log_int;

    /* Summaries: the reduction P0 and P1 always shared, now in one module. */
    gd_grid_summary gs;
    {
        int grc = gd_grid_summarize(grid_T, log_post, grid_T_n, &gs, err);
        if (grc != GD_E_OK) {
            free(gl_nodes); free(gl_weights); free(log_post); free(grid_T);
            return grc;
        }
    }

    memset(result, 0, sizeof(*result));
    result->grid_T = grid_T;
    result->grid_n = grid_T_n;

    result->mean   = gs.mean;
    result->sd     = gs.sd;
    result->median = gs.median;
    result->mode   = gs.mode;
    result->et_lo_50 = gs.et_lo_50; result->et_hi_50 = gs.et_hi_50;
    result->et_lo_80 = gs.et_lo_80; result->et_hi_80 = gs.et_hi_80;
    result->et_lo_90 = gs.et_lo_90; result->et_hi_90 = gs.et_hi_90;
    result->et_lo_95 = gs.et_lo_95; result->et_hi_95 = gs.et_hi_95;
    /* P1 fills these now. Before D2 only P0 did, so the HPD fields this
     * result declares were zero for every model that actually shipped. */
    result->hpd_lo_80 = gs.hpd_lo_80; result->hpd_hi_80 = gs.hpd_hi_80;
    result->hpd_lo_90 = gs.hpd_lo_90; result->hpd_hi_90 = gs.hpd_hi_90;
    result->hpd_lo_95 = gs.hpd_lo_95; result->hpd_hi_95 = gs.hpd_hi_95;

    /* Density output. */
    double *density = (double *)gd_malloc_array_checked(grid_T_n, sizeof(double), err);
    if (density != NULL) {
        for (size_t i = 0; i < grid_T_n; i++) density[i] = exp(log_post[i]);
    }
    result->grid_density = density;

    /* Prior-to-posterior comparison: prior mean for T (uniform) is (T_lo+T_hi)/2. */
    result->prior_mean = (T_lo + T_hi) / 2.0;

    /* Prior dominance: if posterior SD ≈ prior SD, the data is uninformative.
     * Prior SD (uniform) = (T_hi - T_lo) / sqrt(12). */
    double prior_sd = (T_hi - T_lo) / sqrt(12.0);
    result->prior_dominance = (prior_sd > 0) ? (result->sd / prior_sd) : 1.0;

    /* Absolute anchoring. */
    result->absolute_anchored = (a1 > 0.0 || a2 > 0.0);

    free(gl_nodes); free(gl_weights); free(log_post);
    /* grid_T and density owned by result. */
    return 0;
}

/* ---- SBC ---- */

int gd_p1_sbc(gd_rng *rng, size_t n_sbc, size_t n_concepts,
              double T_lo, double T_hi,
              double log_rate_mean, double log_rate_sd,
              double beta_sd, double eps,
              gd_sensitivity sens,
              uint32_t *ranks, gd_error *err)
{
    if (rng == NULL || ranks == NULL || n_sbc == 0) {
        gd_error_set(err, GD_E_INVALID_ARG, "sbc: bad args");
        return GD_E_INVALID_ARG;
    }
    /* For each replicate:
     * 1. Draw T_true ~ Uniform(T_lo, T_hi).
     * 2. Draw beta ~ N(0, beta_sd^2)  (the unknown, to-be-marginalized offset).
     * 3. For each concept, draw log r_c ~ N(log_rate_mean + beta, log_rate_sd^2).
     * 4. Simulate matches with rate r_c and d = 2*T_true.
     * 5. Fit P1 posterior, marginalizing beta over its prior N(0, beta_sd^2).
     * 6. Compute rank of T_true in the posterior (how many posterior samples
     *    are below T_true, divided by total → rank in [0, n_grid).
     */
    gd_concept_ref *refs = (gd_concept_ref *)gd_malloc_array_checked(n_concepts, sizeof(gd_concept_ref), err);
    int *matches = (int *)gd_malloc_array_checked(n_concepts, sizeof(int), err);
    if (refs == NULL || matches == NULL) {
        free(refs); free(matches);
        gd_error_set(err, GD_E_NOMEM, "sbc: alloc");
        return GD_E_NOMEM;
    }

    for (size_t rep = 0; rep < n_sbc; rep++) {
        /* Draw true parameters. */
        double T_true = T_lo + (T_hi - T_lo) * gd_rng_uniform_closed_open(rng);
        double beta = gd_rng_normal(rng, 0.0, beta_sd);

        /* Draw per-concept rates and simulate. */
        double d = 2.0 * T_true;
        for (size_t c = 0; c < n_concepts; c++) {
            double log_r = log_rate_mean + beta + gd_rng_normal(rng, 0.0, log_rate_sd);
            double r = exp(log_r);
            double s = exp(-r * d);
            double u = gd_rng_uniform_closed_open(rng);
            matches[c] = (u < s) ? 1 : 0;
            /* The reference prior for fitting: same as the generating prior. */
            refs[c].log_rate_mean = log_rate_mean;
            refs[c].log_rate_sd = log_rate_sd;
        }

        /* Fit P1 posterior. Inference is given only the β *prior*
         * (beta_mean = 0, beta_sd), NOT the drawn β — the family clock offset
         * is unknown at inference time and must be marginalized. Passing the
         * true β here would make SBC test a posterior conditional on the truth
         * and hide any miscalibration of the β-marginal. */
        gd_pair_result result;
        int rc = gd_p1_posterior(refs, n_concepts, matches, NULL, 0, 0,
                                 T_lo, T_hi, 0.0, beta_sd, eps, 1.0, 19.0,
                                 sens, 160, 14, &result, err);
        if (rc != 0) {
            free(refs); free(matches);
            return rc;
        }

        /* Compute rank: fraction of posterior mass below T_true. */
        double mass_below = 0.0;
        for (size_t i = 1; i < result.grid_n; i++) {
            double dTi = result.grid_T[i] - result.grid_T[i - 1];
            double f0 = result.grid_density[i - 1];
            double f1 = result.grid_density[i];
            double m = 0.5 * (f0 + f1) * dTi;
            if (result.grid_T[i] < T_true) mass_below += m;
        }
        /* Rank on a fixed [0, 256) resolution, independent of the internal
         * T-grid size, so callers can bin without knowing that grid. */
        double r256 = mass_below * 256.0;
        if (r256 > 255.0) r256 = 255.0;
        if (r256 < 0.0) r256 = 0.0;
        ranks[rep] = (uint32_t)r256;

        free((void *)result.grid_T);
        free((void *)result.grid_density);
    }

    free(refs);
    free(matches);
    return 0;
}

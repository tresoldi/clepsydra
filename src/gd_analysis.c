/*
 * gd_analysis.c — the shared runner behind the CLI and the WebAssembly shim.
 *
 * This is the pipeline that used to live inside cli/commands.c, moved here
 * unchanged in its arithmetic (PLAN.md M10): the acceptance gate for the move
 * is that every example reproduces its previous output byte for byte, so the
 * order of interning, the tie-breaking on primary_weight and the numeric
 * parsing are all deliberately preserved rather than tidied.
 *
 * The caller supplies buffers and formats the results; nothing here touches the
 * filesystem, stdio or exit().
 */
#include "gd_analysis.h"
#include "gd_memory.h"
#include "gd_tree_build.h"

#include "gd_tsv.h"
#include "gd_forms.h"
#include "gd_newick.h"
#include "gd_calibration.h"
#include "gd_model_novel_state.h"
#include "gd_mcmc.h"
#include "gd_rng.h"
#include "gd_pairwise.h"
#include "gd_decay.h"

#include <math.h>
#include <stdio.h>      /* snprintf, to compose a labelled error message */
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------- utilities */

/* strdup is POSIX, not C99. */
/*
 * The status a failed path reports. Every `goto done` in this file is taken
 * after something recorded the reason in `err`, so the code to return is
 * whatever the error says; GD_E_INTERNAL is the answer only if nothing did,
 * which would itself be a bug. This keeps the failure branches readable
 * without any of them having to restate the code they just set.
 */
static int fail_code(const gd_error *err)
{
    return (err != NULL && err->code != GD_E_OK) ? err->code : GD_E_INTERNAL;
}

static char *dupstr(const char *s)
{
    return gd_strdup_n(s, strlen(s), NULL);
}

/* Prefixes an error with the input's name ("forms.tsv: ..."), so a parse
 * failure still says which of the three inputs it came from even though the
 * parsers never see a filename. Unconditional now that the error carries its
 * own message: there is no caller buffer to be missing. */
static void label_error(gd_error *err, const char *name)
{
    if (name == NULL || err == NULL) return;
    /* Composed in a scratch buffer, not straight into err->message: the
     * message is one of the inputs, and vsnprintf may not read a buffer it is
     * writing. The explicit field widths keep the whole thing inside one
     * message and stop the name crowding out the reason. */
    char tmp[GD_ERROR_MSG_MAX];
    snprintf(tmp, sizeof tmp, "%.64s: %.180s", name,
             err->message[0] != '\0' ? err->message : "parse failed");
    gd_error_set(err, err->code, tmp);
}

/* -------------------------------------------------------------- TSV table */

typedef struct {
    char  **header;
    size_t  ncols;
    char ***rows;
    size_t  nrows, cap;
} tsv_table;

static char *field_dup(const gd_tsv_row *row, size_t i)
{
    size_t len = (i < row->nfields) ? row->lengths[i] : 0;
    const char *src = (i < row->nfields) ? row->fields[i] : "";
    char *s = gd_strdup_n(src, len, NULL);
    if (s == NULL) return NULL;
    return s;
}

static int tbl_header_cb(const gd_tsv_row *row, void *user, gd_error *err)
{
    tsv_table *t = (tsv_table *)user;
    t->ncols = row->nfields;
    t->header = (char **)gd_calloc_checked(row->nfields ? row->nfields : 1, sizeof(char *), err);
    if (t->header == NULL) { return gd_error_set(err, GD_E_NOMEM, "out of memory"); }
    for (size_t i = 0; i < row->nfields; i++) {
        t->header[i] = field_dup(row, i);
        if (t->header[i] == NULL) {
            return gd_error_set(err, GD_E_NOMEM, "out of memory");
        }
    }
    return 0;
}

static int tbl_row_cb(const gd_tsv_row *row, void *user, gd_error *err)
{
    tsv_table *t = (tsv_table *)user;
    if (t->nrows == t->cap) {
        size_t nc = t->cap ? t->cap * 2 : 64;
        char ***nr = (char ***)gd_realloc_array_checked(t->rows, nc, sizeof(char **), err);
        if (nr == NULL) { return gd_error_set(err, GD_E_NOMEM, "out of memory"); }
        t->rows = nr; t->cap = nc;
    }
    char **cells = (char **)gd_calloc_checked(t->ncols ? t->ncols : 1, sizeof(char *), err);
    if (cells == NULL) { return gd_error_set(err, GD_E_NOMEM, "out of memory"); }
    for (size_t i = 0; i < t->ncols; i++) {
        cells[i] = field_dup(row, i);
        if (cells[i] == NULL) {
            for (size_t k = 0; k < i; k++) free(cells[k]);
            free(cells);
            return gd_error_set(err, GD_E_NOMEM, "out of memory");
        }
    }
    t->rows[t->nrows++] = cells;
    return 0;
}

/*
 * Loads a whole TSV buffer into memory. `len == 0` means NUL-terminated.
 *
 * The path argument is deliberately NULL: gd_tsv_parse reads from `path` when
 * it is given one and ignores the buffer entirely, which would quietly turn
 * this module back into a filesystem client — exactly what the browser cannot
 * be. The input's name is only a label; label_error() folds it into the message.
 */
static int table_load(const char *buf, size_t len, tsv_table *t, gd_error *err)
{
    memset(t, 0, sizeof *t);
    if (buf == NULL) { return gd_error_set(err, GD_E_INVALID_ARG, "missing input buffer"); }
    if (len == GD_LEN_NUL_TERMINATED) len = strlen(buf);
    return gd_tsv_parse(buf, len, tbl_header_cb, tbl_row_cb, t, err);
}

static void table_free(tsv_table *t)
{
    for (size_t r = 0; r < t->nrows; r++) {
        if (t->rows[r] == NULL) continue;
        for (size_t c = 0; c < t->ncols; c++) free(t->rows[r][c]);
        free(t->rows[r]);
    }
    free(t->rows);
    if (t->header != NULL)
        for (size_t c = 0; c < t->ncols; c++) free(t->header[c]);
    free(t->header);
    memset(t, 0, sizeof *t);
}

static size_t table_col(const tsv_table *t, const char *name)
{
    for (size_t i = 0; i < t->ncols; i++)
        if (t->header[i] && strcmp(t->header[i], name) == 0) return i;
    return SIZE_MAX;
}

/*
 * A numeric cell, or `def` when the column is absent or the cell is the missing
 * marker. Returns -1 on text that is not a finite number: a malformed
 * calibration parameter silently becoming 0.0 is how a typo turns into a
 * confidently wrong date.
 */
static int cell_double(const tsv_table *t, size_t col, size_t r,
                       double def, double *out)
{
    if (col == SIZE_MAX) { *out = def; return GD_E_OK; }
    const char *s = t->rows[r][col];
    if (s == NULL || strcmp(s, ".") == 0) { *out = def; return GD_E_OK; }
    double v;
    if (gd_tsv_parse_double(s, strlen(s), &v) != 0 || !isfinite(v)) {
        return GD_E_BAD_NUMERIC;
    }
    *out = v;
    return GD_E_OK;
}

/* -------------------------------------------------------------- summaries */

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Posterior summary with the narrowest `level`-mass interval. Identical to the
 * summarize() that lived in cli/commands.c. */
static gd_analysis_summary summarize(const double *x, size_t n, double level)
{
    gd_analysis_summary s = {0, 0, 0, 0, 0};
    if (n == 0) return s;
    double m = 0.0;
    for (size_t i = 0; i < n; i++) m += x[i];
    m /= (double)n;
    double v = 0.0;
    for (size_t i = 0; i < n; i++) { double d = x[i] - m; v += d * d; }
    s.mean = m;
    s.sd = (n > 1) ? sqrt(v / (double)(n - 1)) : 0.0;
    /* No err here on purpose: summarize() reports a degenerate summary, not
     * a status, so the NULL is the whole answer. */
    double *sorted = (double *)gd_malloc_array_checked(n, sizeof(double), NULL);
    if (sorted == NULL) return s;
    memcpy(sorted, x, n * sizeof(double));
    qsort(sorted, n, sizeof(double), cmp_double);
    s.median = (n % 2) ? sorted[n / 2]
                       : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
    size_t w = (size_t)ceil(level * (double)n);
    if (w < 1) w = 1;
    if (w > n) w = n;
    double best = HUGE_VAL; size_t bi = 0;
    for (size_t i = 0; i + w - 1 < n; i++) {
        double width = sorted[i + w - 1] - sorted[i];
        if (width < best) { best = width; bi = i; }
    }
    s.hpd_lo = sorted[bi];
    s.hpd_hi = sorted[bi + w - 1];
    free(sorted);
    return s;
}

/* Overlap coefficient of two samples' densities, on a shared histogram of
 * log values: sum over bins of min(p, q). 1 when indistinguishable, 0 when
 * disjoint. */
static double overlap_log(const double *a, size_t na, const double *b, size_t nb)
{
    enum { NB = 48 };
    double lo = HUGE_VAL, hi = -HUGE_VAL;
    for (size_t i = 0; i < na; i++) if (a[i] > 0) { double l = log(a[i]); if (l < lo) lo = l; if (l > hi) hi = l; }
    for (size_t i = 0; i < nb; i++) if (b[i] > 0) { double l = log(b[i]); if (l < lo) lo = l; if (l > hi) hi = l; }
    if (!(hi > lo)) return 1.0;
    double ha[NB] = {0}, hb[NB] = {0}, w = (hi - lo) / NB;
    for (size_t i = 0; i < na; i++) if (a[i] > 0) { size_t k = (size_t)((log(a[i]) - lo) / w); ha[k < NB ? k : NB - 1] += 1.0 / (double)na; }
    for (size_t i = 0; i < nb; i++) if (b[i] > 0) { size_t k = (size_t)((log(b[i]) - lo) / w); hb[k < NB ? k : NB - 1] += 1.0 / (double)nb; }
    double s = 0.0;
    for (int k = 0; k < NB; k++) s += (ha[k] < hb[k]) ? ha[k] : hb[k];
    return s;
}

/* First-order power-scaling sensitivity of g to a component with per-draw
 * log density lp: max(|Cov(g, lp)|, |Cov((g - mean)^2, lp)| / (2 sd)) / sd. */
static double power_sens(const double *g, const double *lp, size_t n)
{
    double mg = 0, ml = 0;
    for (size_t i = 0; i < n; i++) { mg += g[i]; ml += lp[i]; }
    mg /= (double)n; ml /= (double)n;
    double vg = 0, c1 = 0, c2 = 0, m2 = 0;
    for (size_t i = 0; i < n; i++) { double d = g[i] - mg; vg += d * d; m2 += d * d; }
    vg /= (double)(n - 1); m2 /= (double)n;
    if (!(vg > 0.0)) return NAN;
    for (size_t i = 0; i < n; i++) {
        double d = g[i] - mg, e = lp[i] - ml;
        c1 += d * e;
        c2 += (d * d - m2) * e;
    }
    c1 /= (double)(n - 1); c2 /= (double)(n - 1);
    double sd = sqrt(vg);
    double s_mean = fabs(c1) / sd, s_sd = fabs(c2) / (2.0 * sd) / sd;
    return (s_mean > s_sd) ? s_mean : s_sd;
}

static double var_log(const double *x, size_t n)
{
    double m = 0.0, v = 0.0;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) if (x[i] > 0) { m += log(x[i]); k++; }
    if (k < 2) return NAN;
    m /= (double)k;
    for (size_t i = 0; i < n; i++) if (x[i] > 0) { double d = log(x[i]) - m; v += d * d; }
    return v / (double)(k - 1);
}

/* --------------------------------------------- posterior predictive checks */

/* Summary statistics of a dataset (labels: nconcept x n, -1 = missing):
 * out[0] mean number of classes per concept (concepts with >= 2 observed
 * tips), out[1] share of classes held by one tip, out[2] share of concepts
 * incompatible with the tree (T1 extension fails), out[3 + k] the share of
 * shared concepts for tip pair k. */
static int dataset_stats(gd_t1_ws *ws, const gd_tree *t, const int *labels,
                         size_t nconcept, const size_t *tips, size_t ntip,
                         double *out, gd_error *err)
{
    size_t n = t->nnodes;
    double ksum = 0.0, kcnt = 0.0, single = 0.0, classes = 0.0, incompat = 0.0, used = 0.0;
    size_t npair = ntip * (ntip - 1) / 2;
    for (size_t k = 0; k < npair; k++) out[3 + k] = 0.0;
    double *pden = (double *)gd_calloc_checked(npair ? npair : 1, sizeof(double), err);
    if (pden == NULL) return GD_E_NOMEM;
    gd_t1_options o; gd_t1_options_default(&o);
    o.epsilon = 0.0;
    o.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    for (size_t c = 0; c < nconcept; c++) {
        const int *row = labels + c * n;
        size_t nobs = 0;
        for (size_t a = 0; a < ntip; a++) if (row[tips[a]] >= 0) nobs++;
        if (nobs < 2) continue;
        used += 1.0;
        /* classes and singletons, by brute force over the few tips */
        size_t kc = 0;
        for (size_t a = 0; a < ntip; a++) {
            int la = row[tips[a]];
            if (la < 0) continue;
            bool first = true; size_t size = 0;
            for (size_t b = 0; b < ntip; b++) {
                if (row[tips[b]] != la) continue;
                if (b < a) first = false;
                size++;
            }
            if (first) { kc++; classes += 1.0; if (size == 1) single += 1.0; }
        }
        ksum += (double)kc; kcnt += 1.0;
        gd_cognate_data d = { (int *)row, n };
        gd_t1_concept_result cr;
        if (gd_t1_loglik_concept(ws, 1e-4, &d, &o, &cr, err) != 0) { free(pden); return GD_E_INTERNAL; }
        if (cr.excluded) incompat += 1.0;
        size_t k = 0;
        for (size_t a = 0; a < ntip; a++)
            for (size_t b = a + 1; b < ntip; b++, k++) {
                int la = row[tips[a]], lb = row[tips[b]];
                if (la < 0 || lb < 0) continue;
                pden[k] += 1.0;
                if (la == lb) out[3 + k] += 1.0;
            }
    }
    out[0] = (kcnt > 0.0) ? ksum / kcnt : NAN;
    out[1] = (classes > 0.0) ? single / classes : NAN;
    out[2] = (used > 0.0) ? incompat / used : NAN;
    for (size_t k = 0; k < npair; k++) out[3 + k] = (pden[k] > 0.0) ? out[3 + k] / pden[k] : NAN;
    free(pden);
    return 0;
}

/* Simulate one concept: T1 with per-branch rates, or with probability eps a
 * CRP partition of the observed tips. Missing tips stay missing. */
static void simulate_concept(gd_rng *rng, const gd_tree *t, double rate, const double *bm,
                             double eps, double alpha, const int *mask, int *out)
{
    size_t n = t->nnodes;
    if (gd_rng_uniform_open(rng) < eps) {
        int next = 0;
        size_t seen = 0;
        int *sizes = (int *)calloc(n + 1, sizeof(int));
        for (size_t v = 0; v < n; v++) {
            out[v] = -1;
            if (mask[v] < 0) continue;
            double u = gd_rng_uniform_open(rng) * ((double)seen + alpha);
            int cls = next;
            double acc = 0.0;
            for (int k = 0; k < next && sizes != NULL; k++) {
                acc += (double)sizes[k];
                if (u < acc) { cls = k; break; }
            }
            if (cls == next) next++;
            if (sizes != NULL) sizes[cls]++;
            out[v] = cls;
            seen++;
        }
        free(sizes);
        return;
    }
    int next = 1;
    for (size_t i = 0; i < n; i++) {
        size_t v = t->preorder[i];
        if (v == t->root) { out[v] = 0; continue; }
        size_t p = t->nodes[v].parent;
        double dt = t->nodes[p].age - t->nodes[v].age;
        if (dt < 0.0) dt = 0.0;
        double r = rate * ((bm != NULL) ? bm[v] : 1.0);
        out[v] = (gd_rng_uniform_open(rng) < exp(-r * dt)) ? out[p] : next++;
    }
    for (size_t v = 0; v < n; v++) if (mask[v] < 0) out[v] = -1;
}

/* ------------------------------------------------- calibration variation */

/* Loosen a soft window about its centre. A fixed tip or a hard bound has no
 * soft window to loosen and is left exactly as it was — widening a stated fact
 * would be inventing doubt rather than measuring it. */
static void widen_window(gd_calibration *cal, double f)
{
    gd_cal_widen(cal, f);
}

void gd_analysis_calibrations_free(gd_analysis_calibration *cals, size_t n)
{
    if (cals == NULL) return;
    for (size_t i = 0; i < n; i++) {
        free(cals[i].id);
        free(cals[i].target);
        free(cals[i].type_name);
        free(cals[i].source_key);
    }
    free(cals);
}

int gd_analysis_calibrations_read(const char *calibrations_tsv, size_t len,
                                  gd_analysis_calibration **out, size_t *n,
                                  gd_error *err)
{
    if (calibrations_tsv == NULL || out == NULL || n == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "gd_analysis_calibrations_read: null argument");
    }
    *out = NULL; *n = 0;

    tsv_table ct; memset(&ct, 0, sizeof ct);
    gd_analysis_calibration *cals = NULL;
    int rc = GD_E_INTERNAL;

    if (table_load(calibrations_tsv, len, &ct, err) != 0) goto done;
    size_t cc_id = table_col(&ct, "calibration_id"), cc_t = table_col(&ct, "target_id");
    size_t cc_ty = table_col(&ct, "type"), cc_p1 = table_col(&ct, "p1");
    size_t cc_p2 = table_col(&ct, "p2"), cc_p3 = table_col(&ct, "p3");
    size_t cc_h = table_col(&ct, "is_hard"), cc_s = table_col(&ct, "source_key");
    size_t cc_o = table_col(&ct, "outlier_prob");
    if (cc_t == SIZE_MAX || cc_ty == SIZE_MAX) {
        gd_error_set(err, GD_E_UNKNOWN_COLUMN,
                     "calibrations need the columns target_id and type");
        goto done;
    }

    cals = (gd_analysis_calibration *)gd_calloc_checked(ct.nrows ? ct.nrows : 1, sizeof *cals, err);
    if (cals == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }

    for (size_t r = 0; r < ct.nrows; r++) {
        gd_analysis_calibration *c = &cals[r];
        c->id = dupstr(cc_id != SIZE_MAX ? ct.rows[r][cc_id] : "cal");
        c->target = dupstr(ct.rows[r][cc_t]);
        c->type_name = dupstr(ct.rows[r][cc_ty]);
        if (c->id == NULL || c->target == NULL || c->type_name == NULL) {
            gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
        }
        if (!gd_cal_type_from_string(c->type_name, strlen(c->type_name), &c->type)) {
            gd_error_setf(err, GD_E_PARSE,
                          "unknown calibration type '%s'", c->type_name);
            goto done;
        }
        if (cell_double(&ct, cc_p1, r, 0.0, &c->p1) != 0 ||
            cell_double(&ct, cc_p2, r, 0.0, &c->p2) != 0 ||
            cell_double(&ct, cc_p3, r, 0.0, &c->p3) != 0) {
            gd_error_setf(err, GD_E_BAD_NUMERIC,
                          "non-numeric parameter for calibration '%s'", c->id);
            goto done;
        }
        /* is_hard means something only for uniform_node, where it chooses
         * between a hard window and a soft one with 2.5% tails. Every other
         * type is hard or soft by construction, and a value contradicting
         * that is refused rather than silently ignored, as it used to be. */
        {
            int given = 0; long hv = 0;
            if (cc_h != SIZE_MAX && strcmp(ct.rows[r][cc_h], ".") != 0 &&
                ct.rows[r][cc_h][0] != '\0') {
                if (gd_tsv_parse_int(ct.rows[r][cc_h], strlen(ct.rows[r][cc_h]), &hv) != 0 ||
                    (hv != 0 && hv != 1)) {
                    gd_error_setf(err, GD_E_BAD_NUMERIC,
                                  "is_hard must be 0 or 1 for calibration '%s'", c->id);
                    goto done;
                }
                given = 1;
            }
            bool intrinsic_hard = (c->type == GD_CAL_FIXED_TIP || c->type == GD_CAL_ANCESTOR_TIP ||
                                   c->type == GD_CAL_HARD_MINIMUM ||
                                   c->type == GD_CAL_HARD_MAXIMUM);
            if (c->type == GD_CAL_UNIFORM_NODE) {
                c->is_hard = given ? (hv != 0) : true;
            } else {
                if (given && (hv != 0) != intrinsic_hard) {
                    gd_error_setf(err, GD_E_CALIBRATION,
                                  "calibration '%s': is_hard=%ld contradicts type '%s', which is always %s",
                                  c->id, hv, c->type_name, intrinsic_hard ? "hard" : "soft");
                    goto done;
                }
                c->is_hard = intrinsic_hard;
            }
        }
        c->outlier = 0.0;
        if (cell_double(&ct, cc_o, r, 0.0, &c->outlier) != 0 ||
            !(c->outlier >= 0.0 && c->outlier < 1.0)) {
            gd_error_setf(err, GD_E_BAD_NUMERIC,
                          "outlier_prob must be a number in [0, 1) for calibration '%s'", c->id);
            goto done;
        }
        c->source_key = dupstr(cc_s != SIZE_MAX ? ct.rows[r][cc_s] : ".");
        if (c->source_key == NULL) {
            gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
        }
        c->point = (c->type == GD_CAL_UNIFORM_NODE) ? 0.5 * (c->p1 + c->p2) : c->p1;
        c->is_extant_tip = (c->type == GD_CAL_FIXED_TIP && c->p1 <= 0.0);
    }

    *out = cals; *n = ct.nrows;
    cals = NULL;
    rc = GD_E_OK;

done:
    if (rc != GD_E_OK) rc = fail_code(err);
    gd_analysis_calibrations_free(cals, ct.nrows);
    table_free(&ct);
    return rc;
}

double gd_analysis_trace_draw(const gd_analysis_tree_result *res,
                              size_t param, size_t chain, size_t k)
{
    if (res == NULL || res->trace == NULL) return NAN;
    if (param >= res->nparams || chain >= res->n_chains ||
        k >= res->trace_draws) {
        return NAN;
    }
    return res->trace[(param * res->n_chains + chain) * res->trace_draws + k];
}

const gd_analysis_node_age *gd_analysis_result_node(
        const gd_analysis_tree_result *res, size_t node)
{
    if (res == NULL) return NULL;
    for (size_t i = 0; i < res->nages; i++)
        if (res->ages[i].node == node) return &res->ages[i];
    return NULL;
}

/* ------------------------------------------------------------- date-tree */

void gd_analysis_tree_config_default(gd_analysis_tree_config *cfg)
{
    if (cfg == NULL) return;
    memset(cfg, 0, sizeof *cfg);
    /* The relaxed clock is the default: it predicts held-out concepts as
     * well or better on every family tried, and recovers held-out tip dates
     * where the strict clock is confidently wrong (docs/research/
     * model-uncertainty.md, examples/loco_benchmark). It needs 4,000 +
     * 4,000 to converge on all eight variant runs there; 2,000 + 2,000
     * converged on five. */
    cfg->n_warmup = 4000;
    cfg->n_samples = 4000;
    cfg->n_chains = 4;
    cfg->threads = 0;
    cfg->seed = 42;
    cfg->rate_logmean = -8.80;
    cfg->rate_logsd = 0.80;
    cfg->epsilon = 0.05;
    cfg->sample_epsilon = true;
    cfg->relaxed = true;
    cfg->block_age_scaling = true;
    cfg->root_lex_logmean = -1.2;   /* exp(-1.2) = 0.30 replacements per concept */
    cfg->root_lex_logsd = 1.5;
    cfg->prior_only = false;
    cfg->compare_prior = true;
    cfg->concept_rates = true;
    cfg->loan_threshold = 0.5;
    cfg->concept_rate_prior_logmean = -0.1;
    cfg->concept_rate_prior_logsd = 0.35;
    cfg->hpd_level = 0.90;
    cfg->calibrations.drop_index = GD_ANALYSIS_KEEP_ALL;
    cfg->calibrations.widen = 1.0;
}

/* Per-concept cognate data keyed by tree node; cognate-id strings are interned
 * to integers within each concept, because only equality matters to T1. */
typedef struct {
    size_t  nconcept;
    int    *labels;        /* nconcept * nnodes; class int or -1 */
    size_t  nnodes;
    char  **names;         /* nconcept concept ids, owned */
} cognate_build;

static const char *const TYPE_NAMES[] = {
    "log_rate", "age", "log_sigma_b", "logit_eps", "branch_logmult",
    "concept_dev", "log_sigma_c", "logit_tip_error"
};

struct gd_analysis_tree_job {
    gd_analysis_tree_config cfg;      /* copied: the caller's may not outlive us */
    gd_tree                *tree;     /* owned until a successful run */
    cognate_build           cb;
    gd_cognate_data        *concepts;
    gd_t1_options           opt;      /* must outlive the model adapter */
    gd_calibration_set     *cs;
    gd_t1_ws               *ws;
    gd_mcmc_t1             *model;
    size_t                  n_skipped_rows;
    size_t                  n_unreadable_weights;
    size_t                  n_missing_cognates;
    size_t                  n_loans;
    int                     ran;
};

/* Frees what points at the tree, so the tree can be handed to the result. */
static void job_release_tree_users(gd_analysis_tree_job *job)
{
    if (job->model != NULL) { gd_mcmc_t1_destroy(job->model); job->model = NULL; }
    if (job->ws != NULL) { gd_t1_ws_destroy(job->ws); job->ws = NULL; }
    if (job->cs != NULL) { gd_calset_destroy(job->cs); job->cs = NULL; }
}

void gd_analysis_tree_job_destroy(gd_analysis_tree_job *job)
{
    if (job == NULL) return;
    job_release_tree_users(job);
    free(job->cb.labels);
    if (job->cb.names != NULL)
        for (size_t c = 0; c < job->cb.nconcept; c++) free(job->cb.names[c]);
    free(job->cb.names);
    free(job->concepts);
    gd_tree_destroy(job->tree);       /* NULL once a result took it */
    free(job);
}

/* Tips, counted. (nnodes + 1) / 2 holds only for a binary tree; a polytomy
 * or a unary node makes it wrong. */
static size_t count_leaves(const gd_tree *t)
{
    size_t k = 0;
    if (t == NULL) return 0;
    for (size_t v = 0; v < t->nnodes; v++) k += t->nodes[v].is_leaf ? 1u : 0u;
    return k;
}

size_t gd_analysis_tree_job_languages(const gd_analysis_tree_job *job)
{
    return (job != NULL) ? count_leaves(job->tree) : 0;
}
size_t gd_analysis_tree_job_nodes(const gd_analysis_tree_job *job)
{
    return (job != NULL) ? job->cb.nnodes : 0;
}
size_t gd_analysis_tree_job_concepts(const gd_analysis_tree_job *job)
{
    return (job != NULL) ? job->cb.nconcept : 0;
}
size_t gd_analysis_tree_job_skipped_rows(const gd_analysis_tree_job *job)
{
    return (job != NULL) ? job->n_skipped_rows : 0;
}

int gd_analysis_tree_prepare(const gd_analysis_tree_input *in,
                             const gd_analysis_tree_config *cfg,
                             gd_analysis_tree_job **out, gd_error *err)
{
    if (in == NULL || cfg == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "gd_analysis_tree_prepare: null argument");
    }
    *out = NULL;

    gd_analysis_calibration *cals = NULL;
    size_t ncal = 0;
    gd_forms *forms = NULL;
    size_t  *lang_node = NULL;     /* forms language index -> tree node */
    char ***class_names = NULL;
    size_t *class_count = NULL;
    size_t  nclass_slots = 0;      /* class_names/class_count length */
    gd_analysis_tree_job *job = NULL;
    int rc = GD_E_INTERNAL;

    job = (gd_analysis_tree_job *)gd_calloc_checked(1, sizeof *job, err);
    if (job == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    job->cfg = *cfg;

    /* ---- tree ---- */
    if (in->newick == NULL) {
        gd_error_set(err, GD_E_INVALID_ARG, "no Newick tree supplied");
        goto done;
    }
    size_t newick_len = (in->newick_len == GD_LEN_NUL_TERMINATED)
                      ? strlen(in->newick) : in->newick_len;
    if (gd_newick_parse(in->newick, newick_len, &job->tree, err) != 0) {
        label_error(err, in->newick_name);
        goto done;
    }
    /* Unlabelled internal nodes get stable "mrca:X,Y" names, so every node
     * in every output can be told apart and targeted by a calibration. */
    if (gd_tree_name_internal_nodes(job->tree, err) != 0) goto done;
    size_t n = job->tree->nnodes;
    job->cb.nnodes = n;

    /* ---- forms -> cognate data ---- */
    {
        size_t flen = (in->forms_tsv == NULL) ? 0
                    : ((in->forms_len == GD_LEN_NUL_TERMINATED)
                        ? strlen(in->forms_tsv) : in->forms_len);
        gd_forms_options fo; gd_forms_options_default(&fo);
        fo.analysis_id = in->analysis_id;
        fo.loan_threshold = job->cfg.loan_threshold;
        if (gd_forms_read_opts(in->forms_tsv, flen, &fo, &forms, err) != GD_E_OK) {
            label_error(err, in->forms_name);
            goto done;
        }
    }
    if (forms->nconcept == 0) {
        gd_error_set(err, GD_E_MISSING_FIELD, "forms contains no data rows");
        goto done;
    }
    job->cb.nconcept = forms->nconcept;
    job->cb.names = (char **)gd_calloc_checked(forms->nconcept, sizeof(char *), err);
    if (job->cb.names == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    for (size_t c = 0; c < forms->nconcept; c++) {
        job->cb.names[c] = dupstr(forms->concepts[c]);
        if (job->cb.names[c] == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    }
    job->n_unreadable_weights = forms->n_unreadable_weights;
    job->n_missing_cognates = forms->n_missing_cognates;
    job->n_loans = forms->n_loans;

    /* Rows whose language the tree never mentions are not data for this run.
     * Counted per row, not per language, because that is what the dataset
     * report has always meant by "skipped". */
    lang_node = (size_t *)gd_malloc_array_checked(forms->nlang ? forms->nlang : 1,
                                                  sizeof(size_t), err);
    if (lang_node == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    for (size_t l = 0; l < forms->nlang; l++) {
        lang_node[l] = gd_tree_find_label(job->tree, forms->languages[l]);
        /* Forms describe languages, and languages are tips. A label match on
         * an internal node (a clade named like a language) would give that
         * node observations it cannot have under the model. */
        if (lang_node[l] != GD_NODE_NONE && !gd_tree_is_leaf(job->tree, lang_node[l])) {
            gd_error_setf(err, GD_E_INVALID_ARG,
                          "forms language '%s' names an internal node of the tree, not a tip",
                          forms->languages[l]);
            goto done;
        }
        if (lang_node[l] == GD_NODE_NONE) job->n_skipped_rows += forms->rows_per_lang[l];
    }

    /* concepts x nodes: both come from the input, so the product is checked
     * before it becomes an allocation size. */
    size_t ncell;
    if (!gd_size_mul(job->cb.nconcept, n, &ncell)) {
        gd_error_set(err, GD_E_NOMEM, "forms: concept x node count overflows");
        goto done;
    }
    job->cb.labels = (int *)gd_malloc_array_checked(ncell, sizeof(int), err);
    nclass_slots = job->cb.nconcept;
    class_names = (char ***)gd_calloc_checked(job->cb.nconcept, sizeof(char **), err);
    class_count = (size_t *)gd_calloc_checked(job->cb.nconcept, sizeof(size_t), err);
    if (job->cb.labels == NULL || class_names == NULL || class_count == NULL) {
        gd_error_set(err, GD_E_NOMEM, "out of memory");
        goto done;
    }
    for (size_t i = 0; i < ncell; i++) job->cb.labels[i] = -1;

    /* Replay the assignments in file order, interning each cognate id within
     * its concept the first time it takes a cell. The order is the file's,
     * not the winners': T1 sums lgamma over label-sorted blocks, so numbering
     * the classes differently is a different summation order. */
    for (size_t ai = 0; ai < forms->nassignments; ai++) {
        const gd_forms_assignment *a = &forms->assignments[ai];
        size_t node = lang_node[a->lang];
        if (node == GD_NODE_NONE) continue;
        const char *cog = a->cognate;
        int cls = -1;
        for (size_t k = 0; k < class_count[a->concept]; k++)
            if (strcmp(class_names[a->concept][k], cog) == 0) { cls = (int)k; break; }
        if (cls < 0) {
            char **grown = (char **)gd_realloc_array_checked(class_names[a->concept],
                                            class_count[a->concept] + 1, sizeof(char *), err);
            if (grown == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
            class_names[a->concept] = grown;
            class_names[a->concept][class_count[a->concept]] = dupstr(cog);
            if (class_names[a->concept][class_count[a->concept]] == NULL) {
                gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
            }
            cls = (int)class_count[a->concept]++;
        }
        job->cb.labels[a->concept * n + node] = cls;
    }

    /* A tip with no observed form is a name the forms file does not know
     * (usually a spelling mismatch with the tree). Dating it would sample its
     * branch from the prior and report that as a result. */
    for (size_t v = 0; v < n; v++) {
        if (!gd_tree_is_leaf(job->tree, v)) continue;
        size_t c = 0;
        while (c < job->cb.nconcept && job->cb.labels[c * n + v] < 0) c++;
        if (c == job->cb.nconcept) {
            const char *name = job->tree->nodes[v].label;
            gd_error_setf(err, GD_E_MISSING_FIELD,
                          "tree tip '%s' has no observed form in the forms file",
                          (name != NULL) ? name : "(unnamed)");
            goto done;
        }
    }

    job->concepts = (gd_cognate_data *)gd_malloc_array_checked(job->cb.nconcept, sizeof(gd_cognate_data), err);
    if (job->concepts == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    for (size_t c = 0; c < job->cb.nconcept; c++) {
        job->concepts[c].labels = job->cb.labels + c * n;
        job->concepts[c].nnodes = n;
    }

    /* ---- calibrations ---- */
    /* Through the public reader, not a second decode of the same table.
     * gd_analysis.h says why: two readings of one file risk disagreeing, and
     * before D3 this function was the disagreement it warns about. */
    if (gd_analysis_calibrations_read(in->calibrations_tsv, in->cal_len,
                                      &cals, &ncal, err) != GD_E_OK) {
        label_error(err, in->cal_name);
        goto done;
    }

    /* One row may be withheld (leave-one-calibration-out) and every soft window
     * may be loosened; both are how a caller asks what the chronology owes to
     * its calibrations rather than to its data. */
    size_t drop = (job->cfg.calibrations.drop_index < ncal) ? 1 : 0;
    double widen = job->cfg.calibrations.widen;
    /* Withholding a tip's age must leave that tip free to be estimated. An
     * uncalibrated tip is extant (age 0), so dropping the row outright would
     * not withhold the date, it would replace it with "today". The dropped
     * row is therefore replaced by a bare hard minimum of 0 on the same
     * target, which marks the tip as dated-but-unknown. */
    bool drop_tip = false;
    if (drop) {
        size_t tv = GD_NODE_NONE;
        gd_error ignore = GD_ERROR_INIT;   /* an unresolvable target fails later, with its name */
        (void)gd_tree_resolve_target(job->tree, cals[job->cfg.calibrations.drop_index].target,
                                     &tv, &ignore);
        drop_tip = (tv != GD_NODE_NONE && gd_tree_is_leaf(job->tree, tv) &&
                    cals[job->cfg.calibrations.drop_index].type != GD_CAL_TIP_ERROR);
    }
    job->cs = gd_calset_create(ncal - drop + (drop_tip ? 1 : 0), err);
    if (job->cs == NULL) goto done;
    size_t out_row = 0;
    for (size_t rr = 0; rr < ncal; rr++) {
        if (drop && rr == job->cfg.calibrations.drop_index) {
            if (!drop_tip) continue;
            gd_calibration *free_tip = &job->cs->cals[out_row++];
            free_tip->calibration_id = dupstr("withheld_tip_age");
            free_tip->target_id      = dupstr(cals[rr].target);
            free_tip->source_key     = dupstr(".");
            if (free_tip->calibration_id == NULL || free_tip->target_id == NULL ||
                free_tip->source_key == NULL) {
                gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
            }
            free_tip->type    = GD_CAL_HARD_MINIMUM;
            free_tip->p1      = 0.0;
            free_tip->is_hard = true;
            continue;
        }
        const gd_analysis_calibration *src = &cals[rr];
        gd_calibration *cal = &job->cs->cals[out_row++];
        cal->calibration_id = dupstr(src->id);
        cal->target_id      = dupstr(src->target);
        cal->source_key     = dupstr(src->source_key);
        if (cal->calibration_id == NULL || cal->target_id == NULL ||
            cal->source_key == NULL) {
            gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
        }
        cal->type    = src->type;
        cal->p1      = src->p1;
        cal->p2      = src->p2;
        cal->p3      = src->p3;
        cal->is_hard = src->is_hard;
        cal->outlier = src->outlier;
        if (widen > 1.0) widen_window(cal, widen);
    }
    if (gd_calibrations_resolve(job->cs, job->tree, err) != 0 ||
        gd_calibrations_propagate(job->cs, job->tree, err) != 0)
        goto done;

    /* ---- model ---- */
    job->ws = gd_t1_ws_create(job->tree, err);
    if (job->ws == NULL) goto done;
    gd_t1_options_default(&job->opt);
    job->opt.epsilon = job->cfg.epsilon;
    /* A concept whose pattern no T1 history on this topology can produce
     * (usually ancestral polymorphism sorted differently in the daughters,
     * sometimes borrowing) has an age-independent likelihood term, so dropping
     * it changes no date. It is dropped rather than charged to the noise
     * component, because counting it as noise inflates the sampled epsilon,
     * which then pulls fast but compatible concepts into the noise too. The
     * count is reported (n_incompatible). */
    job->opt.incompatible = GD_INCOMPATIBLE_EXCLUDE;
    job->model = gd_mcmc_t1_create(job->ws, job->concepts, job->cb.nconcept,
                                   &job->opt, err);
    if (job->model == NULL) goto done;

    *out = job;  job = NULL;
    rc = GD_E_OK;

done:
    if (rc != GD_E_OK) rc = fail_code(err);
    if (class_names != NULL) {
        for (size_t c = 0; c < nclass_slots; c++) {
            if (class_names[c] == NULL) continue;
            /* class_count is allocated alongside class_names and grows in step
             * with it. It can still be NULL here: if its own allocation was the
             * one that failed, this cleanup runs with class_names live and
             * class_count absent. Nothing has been added to the per-class
             * arrays at that point, so zero is the right count — reading
             * class_count[c] unconditionally would dereference NULL. */
            size_t nk = (class_count != NULL) ? class_count[c] : 0;
            for (size_t k = 0; k < nk; k++) free(class_names[c][k]);
            free(class_names[c]);
        }
        free(class_names);
    }
    free(class_count);
    free(lang_node);
    gd_forms_destroy(forms);
    gd_analysis_calibrations_free(cals, ncal);
    gd_analysis_tree_job_destroy(job);   /* no-op when the job was handed over */
    return rc;
}

static void free_concepts(gd_analysis_concept *c, size_t n)
{
    if (c == NULL) return;
    for (size_t i = 0; i < n; i++) free(c[i].name);
    free(c);
}

void gd_analysis_tree_result_destroy(gd_analysis_tree_result *res)
{
    if (res == NULL) return;
    free_concepts(res->concepts, res->n_concepts);
    free(res->checks);
    free(res->pointwise);
    if (res->cal_checks != NULL)
        for (size_t i = 0; i < res->n_cal_checks; i++) free(res->cal_checks[i].id);
    free(res->cal_checks);
    free(res->trace);
    free(res->ages);
    free(res->params);
    gd_tree_destroy(res->tree);
    free(res);
}

int gd_analysis_tree_job_run(gd_analysis_tree_job *job,
                             gd_analysis_tree_result **out, gd_error *err)
{
    if (job == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "gd_analysis_tree_job_run: null argument");
    }
    *out = NULL;
    if (job->ran || job->model == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "this analysis job has already run");
    }

    gd_mcmc_result *res = NULL, *prior_res = NULL;
    double *hz = NULL;
    double *col = NULL, *pcol = NULL, *rate_col = NULL, *root_col = NULL;
    gd_analysis_tree_result *r = NULL;
    int rc = GD_E_INTERNAL;

    gd_mcmc_model model; gd_mcmc_t1_model(job->model, &model);
    gd_mcmc_config mc; gd_mcmc_config_default(&mc);
    mc.seed = job->cfg.seed;
    mc.n_warmup = job->cfg.n_warmup;
    mc.n_samples = job->cfg.n_samples;
    mc.thin = 1;
    mc.rate_logmean = job->cfg.rate_logmean;
    mc.rate_logsd = job->cfg.rate_logsd;
    mc.sample_branch = job->cfg.relaxed ? 1 : 0;
    mc.fixed_epsilon = job->cfg.epsilon;
    mc.sample_contamination = job->cfg.sample_epsilon;
    mc.block_age_scaling = job->cfg.block_age_scaling;
    mc.root_lex_logmean = job->cfg.root_lex_logmean;
    mc.root_lex_logsd = job->cfg.root_lex_logsd;
    mc.prior_only = job->cfg.prior_only;
    mc.store_pointwise = job->cfg.keep_pointwise;
    mc.tree_prior = job->cfg.yule ? GD_TREE_PRIOR_YULE : GD_TREE_PRIOR_UNIFORM;
    mc.sample_rate_mix = job->cfg.concept_rates;
    mc.rate_mix_family = job->cfg.concept_rate_gamma ? 1 : 0;
    /* Tip observation-error priors (E2.9), from tip_error rows. */
    double *te_mean = NULL, *te_conc = NULL;
    for (size_t c = 0; c < job->cs->ncals; c++) {
        const gd_calibration *cl = &job->cs->cals[c];
        if (!cl->resolved || cl->type != GD_CAL_TIP_ERROR) continue;
        if (te_mean == NULL) {
            te_mean = (double *)gd_calloc_checked(job->tree->nnodes, sizeof(double), err);
            te_conc = (double *)gd_calloc_checked(job->tree->nnodes, sizeof(double), err);
            if (te_mean == NULL || te_conc == NULL) {
                free(te_mean); free(te_conc);
                gd_error_set(err, GD_E_NOMEM, "out of memory"); return fail_code(err);
            }
        }
        te_mean[cl->node_index] = cl->p1;
        te_conc[cl->node_index] = cl->p2;
    }
    mc.tip_err_mean = te_mean;
    mc.tip_err_conc = te_conc;
    mc.rate_mix_prior_logmean = job->cfg.concept_rate_prior_logmean;
    mc.rate_mix_prior_logsd = job->cfg.concept_rate_prior_logsd;
    mc.yule_kappa = job->cfg.yule_kappa;
    mc.on_progress = job->cfg.on_progress;
    mc.threads = job->cfg.threads;
    mc.progress_user = job->cfg.progress_user;
    mc.progress_every = job->cfg.progress_every;

    job->ran = 1;
    res = gd_mcmc_run(&model, job->tree, job->cs, &mc, job->cfg.n_chains, err);
    if (res == NULL) { free(te_mean); free(te_conc); goto done; }

    /* The prior alone, for the effective-prior report (E1.6). The likelihood
     * is switched off, so this costs a small fraction of the posterior run.
     * Not reported to the progress callback: it is part of this run. */
    if (job->cfg.compare_prior && !job->cfg.prior_only) {
        gd_mcmc_config pc = mc;
        pc.prior_only = true;
        pc.seed = mc.seed ^ 0x9e3779b97f4a7c15ULL;
        pc.on_progress = NULL;
        prior_res = gd_mcmc_run(&model, job->tree, job->cs, &pc, job->cfg.n_chains, err);
        if (prior_res == NULL) { free(te_mean); free(te_conc); goto done; }
    }
    free(te_mean); free(te_conc);
    /* The chains are gone, and the workspace still points at the last
     * chain's tip-error hazards; the summaries below set their own. */
    gd_t1_ws_set_tip_error(job->ws, NULL);
    hz = (double *)gd_calloc_checked(job->tree->nnodes, sizeof(double), err);
    if (hz == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }

    size_t n = job->tree->nnodes;
    /* res->n_chains * res->n_draws was already checked when the sampler
     * allocated its sample matrix from the same two numbers. */
    size_t total = res->n_chains * res->n_draws, np = res->n_param;
    col = (double *)gd_malloc_array_checked(total ? total : 1, sizeof(double), err);
    pcol = (double *)gd_malloc_array_checked(total ? total : 1, sizeof(double), err);
    rate_col = (double *)gd_malloc_array_checked(total ? total : 1, sizeof(double), err);
    root_col = (double *)gd_malloc_array_checked(total ? total : 1, sizeof(double), err);
    r = (gd_analysis_tree_result *)gd_calloc_checked(1, sizeof *r, err);
    if (col == NULL || pcol == NULL || rate_col == NULL || root_col == NULL || r == NULL) {
        gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
    }

    /* Per-draw rate and root age, for the lexical and relative ages. */
    {
        size_t p_rate = SIZE_MAX, p_root = SIZE_MAX;
        for (size_t p = 0; p < np; p++) {
            if (res->param_type[p] == GD_MCMC_P_LOGRATE) p_rate = p;
            if (res->param_type[p] == GD_MCMC_P_AGE && res->param_node[p] == job->tree->root) p_root = p;
        }
        double fixed_root = job->cs->node_lo[job->tree->root];
        for (size_t i = 0; i < total; i++) {
            rate_col[i] = (p_rate != SIZE_MAX) ? exp(gd_mcmc_result_at(res, p_rate, i)) : NAN;
            root_col[i] = (p_root != SIZE_MAX) ? gd_mcmc_result_at(res, p_root, i) : fixed_root;
        }
    }

    /* Thinned trace, when asked for: evenly spaced draws, chains kept apart. */
    size_t td = job->cfg.trace_draws;
    if (td > res->n_draws) td = res->n_draws;
    if (td > 0 && np > 0) {
        size_t ntrace;
        if (!gd_size_mul(np, res->n_chains, &ntrace) ||
            !gd_size_mul(ntrace, td, &ntrace)) {
            gd_error_set(err, GD_E_NOMEM, "trace size overflows");
            goto done;
        }
        r->trace = (double *)gd_malloc_array_checked(ntrace, sizeof(double), err);
        if (r->trace == NULL) { goto done; }
        r->trace_draws = td;
        for (size_t p = 0; p < np; p++)
            for (size_t c = 0; c < res->n_chains; c++)
                for (size_t k = 0; k < td; k++) {
                    size_t draw = (td == 1) ? 0 : k * (res->n_draws - 1) / (td - 1);
                    r->trace[(p * res->n_chains + c) * td + k] =
                        gd_mcmc_result_at(res, p, c * res->n_draws + draw);
                }
    }

    size_t nages = 0;
    for (size_t p = 0; p < np; p++)
        if (res->param_type[p] == GD_MCMC_P_AGE) nages++;
    r->ages = (gd_analysis_node_age *)gd_calloc_checked(nages ? nages : 1, sizeof(gd_analysis_node_age), err);
    r->params = (gd_analysis_param *)gd_calloc_checked(np ? np : 1, sizeof(gd_analysis_param), err);
    if (r->ages == NULL || r->params == NULL) {
        gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
    }

    double level = (job->cfg.hpd_level > 0.0 && job->cfg.hpd_level < 1.0)
                 ? job->cfg.hpd_level : 0.90;
    size_t ai = 0;
    for (size_t p = 0; p < np; p++) {
        (void)gd_mcmc_result_column(res, p, col, total);
        gd_analysis_summary s = summarize(col, total, level);
        size_t node = res->param_node[p];
        const char *label = NULL;
        if (node < n && job->tree->nodes[node].label != NULL)
            label = job->tree->nodes[node].label;

        gd_analysis_param *pp = &r->params[p];
        pp->type = res->param_type[p];
        pp->type_name = TYPE_NAMES[res->param_type[p]];
        pp->node = node;
        pp->label = label;
        pp->stat = s;
        pp->rhat = res->rhat[p];
        pp->ess = res->ess[p];
        pp->ess_tail = res->ess_tail[p];
        pp->mcse = res->mcse[p];
        pp->acceptance = res->acc_rate[p];

        if (res->param_type[p] == GD_MCMC_P_AGE) {
            gd_analysis_node_age *na = &r->ages[ai++];
            na->node = node;
            na->label = label;
            na->age = s;
            na->rhat = res->rhat[p];
            na->ess = res->ess[p];
            na->ess_tail = res->ess_tail[p];
            na->mcse = res->mcse[p];
            for (size_t i = 0; i < total; i++) pcol[i] = rate_col[i] * col[i];
            na->lexical = summarize(pcol, total, level);
            for (size_t i = 0; i < total; i++) pcol[i] = col[i] / root_col[i];
            na->relative = summarize(pcol, total, level);
            na->prior_overlap = NAN;
            na->contraction = NAN;
            memset(&na->prior, 0, sizeof na->prior);
            if (prior_res != NULL) {
                size_t ptotal = prior_res->n_chains * prior_res->n_draws;
                if (ptotal <= total) {
                    (void)gd_mcmc_result_column(prior_res, p, pcol, ptotal);
                    na->prior = summarize(pcol, ptotal, level);
                    na->prior_overlap = overlap_log(col, total, pcol, ptotal);
                    na->contraction = 1.0 - var_log(col, total) / var_log(pcol, ptotal);
                }
            }
        }
    }

    r->nages = nages;
    r->nparams = np;
    r->n_languages = count_leaves(job->tree);
    r->n_nodes = n;
    r->n_concepts = job->cb.nconcept;
    r->n_chains = res->n_chains;
    r->n_warmup = job->cfg.n_warmup;
    r->n_samples = job->cfg.n_samples;
    r->n_draws = res->n_draws;
    r->n_skipped_rows = job->n_skipped_rows;
    r->n_unreadable_weights = job->n_unreadable_weights;
    r->n_missing_cognates = job->n_missing_cognates;
    r->n_loans = job->n_loans;
    r->relaxed = job->cfg.relaxed;
    r->branch_rate_confounded = res->branch_rate_confounded;

    /* Convergence verdict. A NAN diagnostic (a parameter that never moved, or
     * too few draws) counts as a failure: "could not be assessed" is not
     * "converged". */
    {
        double ess_min = GD_CONV_ESS_PER_CHAIN * (double)res->n_chains;
        double mx = 0.0, mb = HUGE_VAL, mt = HUGE_VAL;
        bool nan_r = false, nan_b = false, nan_t = false;
        r->n_unconverged = 0;
        for (size_t p = 0; p < np; p++) {
            double rh = res->rhat[p], eb = res->ess[p], et = res->ess_tail[p];
            bool ok = (rh < GD_CONV_RHAT_MAX) && (eb >= ess_min) && (et >= ess_min);
            if (!ok) r->n_unconverged++;
            if (isnan(rh)) nan_r = true; else if (rh > mx) mx = rh;
            if (isnan(eb)) nan_b = true; else if (eb < mb) mb = eb;
            if (isnan(et)) nan_t = true; else if (et < mt) mt = et;
        }
        /* The extremes report NAN when any parameter could not be assessed,
         * so a summary line cannot look healthier than the worst case. */
        r->max_rhat = nan_r ? NAN : mx;
        r->min_ess_bulk = nan_b ? NAN : mb;
        r->min_ess_tail = nan_t ? NAN : mt;
        r->converged = (np > 0 && r->n_unconverged == 0);
    }
    r->min_branch_duration_corr = res->min_branch_duration_corr;
    if (job->cfg.keep_pointwise && res->pointwise != NULL) {
        r->pointwise = res->pointwise;   /* moves into the result */
        res->pointwise = NULL;
    }

    /* Where the calendar scale comes from (E1.6). */
    {
        bool node_cal = false, tip_cal = false;
        for (size_t c = 0; c < job->cs->ncals; c++) {
            const gd_calibration *cl = &job->cs->cals[c];
            if (!cl->resolved) continue;
            if (cl->type == GD_CAL_TIP_ERROR) continue;   /* not an age */
            if (job->tree->nodes[cl->node_index].is_leaf) {
                if (!(cl->type == GD_CAL_FIXED_TIP && cl->p1 <= 0.0)) tip_cal = true;
            } else {
                node_cal = true;
            }
        }
        r->scale_source = (node_cal && tip_cal) ? "calibration+tip"
                        : node_cal ? "calibration" : tip_cal ? "tip" : "clock_prior";
        r->rate_contraction = NAN;
        for (size_t p = 0; p < np; p++) {
            if (res->param_type[p] != GD_MCMC_P_LOGRATE) continue;
            double sd = r->params[p].stat.sd;
            r->rate_contraction = 1.0 - (sd * sd) / (job->cfg.rate_logsd * job->cfg.rate_logsd);
        }
        r->scale_from_prior = !(r->rate_contraction >= 0.5);
        for (size_t i = 0; i < total; i++) pcol[i] = rate_col[i] * root_col[i];
        r->lexical_root = summarize(pcol, total, level);
    }

    /* Variance decomposition and power-scaling sensitivity (E1.6). */
    for (int k = 0; k < 4; k++) { r->sens_root[k] = NAN; r->sens_lexroot[k] = NAN; }
    r->var_share_data = NAN; r->var_share_clock = NAN;
    if (!isnan(rate_col[0]) && total > 2) {
        double *g_root = (double *)gd_malloc_array_checked(total, sizeof(double), err);
        double *g_lex = (double *)gd_malloc_array_checked(total, sizeof(double), err);
        double *lp = (double *)gd_malloc_array_checked(total, sizeof(double), err);
        double *lr = (double *)gd_malloc_array_checked(total, sizeof(double), err);
        if (g_root == NULL || g_lex == NULL || lp == NULL || lr == NULL) {
            free(g_root); free(g_lex); free(lp); free(lr);
            gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
        }
        for (size_t i = 0; i < total; i++) {
            lr[i] = log(rate_col[i]);
            g_root[i] = log(root_col[i]);
            g_lex[i] = lr[i] + g_root[i];
        }
        /* Covariance split: Var(log a) = Cov(log a, log L) + Cov(log a, -log rho). */
        double vroot = 0, clex = 0, crate = 0, mr = 0, ml = 0, mt = 0;
        for (size_t i = 0; i < total; i++) { mr += g_root[i]; ml += g_lex[i]; mt += lr[i]; }
        mr /= (double)total; ml /= (double)total; mt /= (double)total;
        for (size_t i = 0; i < total; i++) {
            vroot += (g_root[i] - mr) * (g_root[i] - mr);
            clex += (g_root[i] - mr) * (g_lex[i] - ml);
            crate -= (g_root[i] - mr) * (lr[i] - mt);
        }
        if (vroot > 0.0) { r->var_share_data = clex / vroot; r->var_share_clock = crate / vroot; }

        /* Component log densities per draw. */
        double sdr = job->cfg.rate_logsd, sdl = job->cfg.root_lex_logsd;
        for (size_t i = 0; i < total; i++) {
            double z = (lr[i] - job->cfg.rate_logmean) / sdr;
            lp[i] = -0.5 * z * z;
        }
        r->sens_root[GD_SENS_CLOCK] = power_sens(g_root, lp, total);
        r->sens_lexroot[GD_SENS_CLOCK] = power_sens(g_lex, lp, total);
        for (size_t i = 0; i < total; i++) {
            double z = (g_lex[i] - job->cfg.root_lex_logmean) / sdl;
            lp[i] = -0.5 * z * z;
        }
        r->sens_root[GD_SENS_LEXROOT] = power_sens(g_root, lp, total);
        r->sens_lexroot[GD_SENS_LEXROOT] = power_sens(g_lex, lp, total);
        /* Calibrations on sampled nodes; fixed tips contribute a constant. */
        bool any_cal = false;
        for (size_t i = 0; i < total; i++) lp[i] = 0.0;
        for (size_t c = 0; c < job->cs->ncals; c++) {
            const gd_calibration *cl = &job->cs->cals[c];
            if (!cl->resolved) continue;
            size_t pv = SIZE_MAX;
            for (size_t p = 0; p < np; p++)
                if (res->param_type[p] == GD_MCMC_P_AGE && res->param_node[p] == cl->node_index) pv = p;
            if (pv == SIZE_MAX) continue;
            any_cal = true;
            for (size_t i = 0; i < total; i++)
                lp[i] += gd_cal_logprior(cl, gd_mcmc_result_at(res, pv, i));
        }
        if (any_cal) {
            r->sens_root[GD_SENS_CALIBRATIONS] = power_sens(g_root, lp, total);
            r->sens_lexroot[GD_SENS_CALIBRATIONS] = power_sens(g_lex, lp, total);
        }
        if (!job->cfg.prior_only) {
            r->sens_root[GD_SENS_LIKELIHOOD] = power_sens(g_root, res->loglik, total);
            r->sens_lexroot[GD_SENS_LIKELIHOOD] = power_sens(g_lex, res->loglik, total);
        }
        free(g_root); free(g_lex); free(lp); free(lr);
    }

    /* Per-concept summary (E2.2), over at most 200 evenly spaced draws. Each
     * draw's ages and parameters are written back into the tree and the
     * concept likelihoods recomputed: the sampler keeps only totals. */
    if (!job->cfg.prior_only && total > 0) {
        size_t nc = job->cb.nconcept;
        r->concepts = (gd_analysis_concept *)gd_calloc_checked(nc ? nc : 1, sizeof *r->concepts, err);
        double *bm = (double *)gd_malloc_array_checked(n, sizeof(double), err);
        if (r->concepts == NULL || bm == NULL) { free(bm); gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
        size_t p_rate = SIZE_MAX, p_eps = SIZE_MAX, p_sc = SIZE_MAX;
        for (size_t p = 0; p < np; p++) {
            if (res->param_type[p] == GD_MCMC_P_LOGRATE) p_rate = p;
            if (res->param_type[p] == GD_MCMC_P_LOGIT_EPS) p_eps = p;
            if (res->param_type[p] == GD_MCMC_P_LOG_SIGMA_C) p_sc = p;
        }
        size_t ndraw = (total < 200) ? total : 200;
        gd_t1_options opt = job->opt;
        bool ok = true;
        for (size_t k = 0; k < ndraw && ok; k++) {
            size_t i = (ndraw == 1) ? 0 : k * (total - 1) / (ndraw - 1);
            for (size_t v = 0; v < n; v++) {
                job->tree->nodes[v].age = job->cs->node_lo[v];
                bm[v] = 1.0;
            }
            for (size_t p = 0; p < np; p++) {
                if (res->param_type[p] == GD_MCMC_P_AGE)
                    job->tree->nodes[res->param_node[p]].age = gd_mcmc_result_at(res, p, i);
                else if (res->param_type[p] == GD_MCMC_P_BRANCH_LOGMULT)
                    bm[res->param_node[p]] = exp(gd_mcmc_result_at(res, p, i));
            }
            double rate = (p_rate != SIZE_MAX) ? exp(gd_mcmc_result_at(res, p_rate, i)) : exp(job->cfg.rate_logmean);
            opt.epsilon = (p_eps != SIZE_MAX) ? 1.0 / (1.0 + exp(-gd_mcmc_result_at(res, p_eps, i)))
                                              : job->cfg.epsilon;
            double sc = (p_sc != SIZE_MAX) ? exp(gd_mcmc_result_at(res, p_sc, i)) : 0.0;
            for (size_t v = 0; v < n; v++) hz[v] = 0.0;
            for (size_t p = 0; p < np; p++)
                if (res->param_type[p] == GD_MCMC_P_LOGIT_TIPERR)
                    hz[res->param_node[p]] = log1p(exp(gd_mcmc_result_at(res, p, i)));
            gd_t1_ws_set_tip_error(job->ws, hz);
            for (size_t c = 0; c < nc; c++) {
                gd_t1_concept_result cr;
                if (gd_t1_loglik_concept_mixed_family(job->ws, rate, job->cfg.concept_rate_gamma ? 1 : 0,
                                               sc, job->cfg.relaxed ? bm : NULL,
                                               &job->concepts[c], &opt, &cr, err) != 0) { ok = false; break; }
                if (cr.excluded) { r->concepts[c].incompatible = true; continue; }
                r->concepts[c].p_contaminated += cr.p_contaminated / (double)ndraw;
                r->concepts[c].rate_mult += cr.rate_mult_mean / (double)ndraw;
                r->concepts[c].mean_loglik += cr.log_lik / (double)ndraw;
            }
        }
        free(bm);
        if (!ok) goto done;
        for (size_t c = 0; c < nc; c++) {
            if (r->concepts[c].incompatible) {
                r->n_incompatible++;
                r->concepts[c].p_contaminated = NAN;
                r->concepts[c].rate_mult = NAN;
                r->concepts[c].mean_loglik = NAN;
            }
            r->concepts[c].name = dupstr(job->cb.names[c]);
            if (r->concepts[c].name == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
            for (size_t v = 0; v < n; v++)
                if (job->tree->nodes[v].is_leaf && job->concepts[c].labels[v] >= 0) r->concepts[c].n_observed++;
        }
    }

    /* Posterior predictive checks (E2.10). */
    if (!job->cfg.prior_only && total > 0) {
        size_t nc = job->cb.nconcept, ntip = 0;
        size_t *tips = (size_t *)gd_malloc_array_checked(n, sizeof(size_t), err);
        if (tips == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
        for (size_t v = 0; v < n; v++) if (job->tree->nodes[v].is_leaf) tips[ntip++] = v;
        size_t npair = ntip * (ntip - 1) / 2, nstat = 3 + npair;
        size_t R = (total < 200) ? total : 200;
        double *obs = (double *)gd_malloc_array_checked(nstat, sizeof(double), err);
        double *rep = (double *)gd_malloc_array_checked(R * nstat, sizeof(double), err);
        int *sim = (int *)gd_malloc_array_checked(nc * n, sizeof(int), err);
        double *bm = (double *)gd_malloc_array_checked(n, sizeof(double), err);
        gd_rng *rng = gd_rng_create(job->cfg.seed ^ 0x51ed270b27f3a7c5ULL, 7u, err);
        r->checks = (gd_analysis_check *)gd_calloc_checked(nstat, sizeof *r->checks, err);
        bool ok = (obs && rep && sim && bm && rng && r->checks);
        if (ok && dataset_stats(job->ws, job->tree, job->cb.labels, nc, tips, ntip, obs, err) != 0) ok = false;
        size_t p_rate = SIZE_MAX, p_eps = SIZE_MAX, p_sc = SIZE_MAX;
        for (size_t p = 0; p < np; p++) {
            if (res->param_type[p] == GD_MCMC_P_LOGRATE) p_rate = p;
            if (res->param_type[p] == GD_MCMC_P_LOGIT_EPS) p_eps = p;
            if (res->param_type[p] == GD_MCMC_P_LOG_SIGMA_C) p_sc = p;
        }
        for (size_t k = 0; k < R && ok; k++) {
            size_t i = (R == 1) ? 0 : k * (total - 1) / (R - 1);
            for (size_t v = 0; v < n; v++) { job->tree->nodes[v].age = job->cs->node_lo[v]; bm[v] = 1.0; }
            for (size_t p = 0; p < np; p++) {
                if (res->param_type[p] == GD_MCMC_P_AGE)
                    job->tree->nodes[res->param_node[p]].age = gd_mcmc_result_at(res, p, i);
                else if (res->param_type[p] == GD_MCMC_P_BRANCH_LOGMULT)
                    bm[res->param_node[p]] = exp(gd_mcmc_result_at(res, p, i));
            }
            double rate = (p_rate != SIZE_MAX) ? exp(gd_mcmc_result_at(res, p_rate, i)) : exp(job->cfg.rate_logmean);
            double eps = (p_eps != SIZE_MAX) ? 1.0 / (1.0 + exp(-gd_mcmc_result_at(res, p_eps, i)))
                                             : job->cfg.epsilon;
            double sc = (p_sc != SIZE_MAX) ? exp(gd_mcmc_result_at(res, p_sc, i)) : 0.0;
            for (size_t v = 0; v < n; v++) hz[v] = 0.0;
            for (size_t p = 0; p < np; p++)
                if (res->param_type[p] == GD_MCMC_P_LOGIT_TIPERR)
                    hz[res->param_node[p]] = log1p(exp(gd_mcmc_result_at(res, p, i)));
            for (size_t c = 0; c < nc; c++) {
                double lam = 1.0;
                if (sc > 0.0)
                    lam = job->cfg.concept_rate_gamma
                        ? gd_rng_gamma(rng, 1.0 / (sc * sc), sc * sc)
                        : exp(sc * gd_rng_normal(rng, 0.0, 1.0) - 0.5 * sc * sc);
                simulate_concept(rng, job->tree, rate * lam, job->cfg.relaxed ? bm : NULL, eps,
                                 job->opt.noise_alpha, job->concepts[c].labels, sim + c * n);
                /* Observation error: a tip's form is replaced by a new class
                 * with probability e = 1 - exp(-hazard). */
                for (size_t v = 0; v < n; v++)
                    if (hz[v] > 0.0 && sim[c * n + v] >= 0 &&
                        gd_rng_uniform_open(rng) < -expm1(-hz[v]))
                        sim[c * n + v] = 1000000 + (int)v;
            }
            if (dataset_stats(job->ws, job->tree, sim, nc, tips, ntip, rep + k * nstat, err) != 0) ok = false;
        }
        if (ok) {
            double *col_r = (double *)gd_malloc_array_checked(R, sizeof(double), err);
            if (col_r == NULL) ok = false;
            for (size_t s = 0; s < nstat && ok; s++) {
                gd_analysis_check *ck = &r->checks[s];
                if (s == 0) snprintf(ck->name, sizeof ck->name, "classes_per_concept");
                else if (s == 1) snprintf(ck->name, sizeof ck->name, "singleton_share");
                else if (s == 2) snprintf(ck->name, sizeof ck->name, "incompatible_share");
                else {
                    size_t k = s - 3, a = 0, b = 1;
                    for (size_t q = 0; q < k; q++) { if (++b == ntip) { a++; b = a + 1; } }
                    const char *la = job->tree->nodes[tips[a]].label, *lb = job->tree->nodes[tips[b]].label;
                    snprintf(ck->name, sizeof ck->name, "sharing:%s|%s", la ? la : "?", lb ? lb : "?");
                }
                ck->observed = obs[s];
                size_t m = 0; double ge = 0.0, mean = 0.0;
                for (size_t k = 0; k < R; k++) {
                    double x = rep[k * nstat + s];
                    if (isnan(x)) continue;
                    col_r[m++] = x; mean += x;
                    if (x > obs[s]) ge += 1.0; else if (x == obs[s]) ge += 0.5;
                }
                if (m == 0 || isnan(obs[s])) {
                    ck->rep_mean = ck->rep_lo90 = ck->rep_hi90 = ck->ppp = NAN;
                    continue;
                }
                qsort(col_r, m, sizeof(double), cmp_double);
                ck->rep_mean = mean / (double)m;
                ck->rep_lo90 = col_r[(size_t)(0.05 * (double)(m - 1))];
                ck->rep_hi90 = col_r[(size_t)(0.95 * (double)(m - 1))];
                ck->ppp = ge / (double)m;
            }
            free(col_r);
            r->n_checks = nstat;
        }
        free(tips); free(obs); free(rep); free(sim); free(bm);
        gd_rng_destroy(rng);
        if (!ok) { gd_error_set(err, GD_E_INTERNAL, "posterior predictive checks failed"); goto done; }
    }

    /* Calibration check (E2.6): each calibration on a sampled node. */
    {
        size_t nck = 0;
        for (size_t c = 0; c < job->cs->ncals; c++) {
            const gd_calibration *cl = &job->cs->cals[c];
            for (size_t p = 0; p < np; p++)
                if (cl->resolved && res->param_type[p] == GD_MCMC_P_AGE &&
                    res->param_node[p] == cl->node_index) { nck++; break; }
        }
        if (nck > 0) {
            r->cal_checks = (struct gd_analysis_cal_check *)gd_calloc_checked(nck, sizeof *r->cal_checks, err);
            if (r->cal_checks == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
            for (size_t c = 0; c < job->cs->ncals; c++) {
                const gd_calibration *cl = &job->cs->cals[c];
                size_t pv = SIZE_MAX;
                for (size_t p = 0; p < np; p++)
                    if (cl->resolved && res->param_type[p] == GD_MCMC_P_AGE &&
                        res->param_node[p] == cl->node_index) pv = p;
                if (pv == SIZE_MAX) continue;
                struct gd_analysis_cal_check *ck = &r->cal_checks[r->n_cal_checks++];
                ck->id = dupstr(cl->calibration_id);
                if (ck->id == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
                ck->target = job->tree->nodes[cl->node_index].label;
                ck->outlier_prior = cl->outlier;
                (void)gd_mcmc_result_column(res, pv, col, total);
                ck->target_age = summarize(col, total, level);
                double post = 0.0;
                for (size_t i = 0; i < total; i++) post += gd_cal_outlier_posterior(cl, col[i]);
                ck->outlier_post = (cl->outlier > 0.0) ? post / (double)total : NAN;
            }
        }
    }

    /* The tree moves into the result; tear down what still points at it first,
     * so the job and the result can be destroyed in either order. */
    job_release_tree_users(job);
    r->tree = job->tree;
    job->tree = NULL;

    *out = r;  r = NULL;
    rc = GD_E_OK;

done:
    if (rc != GD_E_OK) rc = fail_code(err);
    free(hz);
    if (job != NULL && job->ws != NULL) gd_t1_ws_set_tip_error(job->ws, NULL);
    free(col); free(pcol); free(rate_col); free(root_col);
    gd_mcmc_result_destroy(prior_res);
    if (r != NULL) {
        free_concepts(r->concepts, job->cb.nconcept);
        free(r->checks);
        free(r->pointwise);
        if (r->cal_checks != NULL)
            for (size_t i = 0; i < r->n_cal_checks; i++) free(r->cal_checks[i].id);
        free(r->cal_checks);
        free(r->trace); free(r->ages); free(r->params); free(r);
    }
    if (res != NULL) gd_mcmc_result_destroy(res);
    return rc;
}

int gd_analysis_tree_run(const gd_analysis_tree_input *in,
                         const gd_analysis_tree_config *cfg,
                         gd_analysis_tree_result **out, gd_error *err)
{
    gd_analysis_tree_job *job = NULL;
    int rc = gd_analysis_tree_prepare(in, cfg, &job, err);
    if (rc != GD_E_OK) return rc;
    rc = gd_analysis_tree_job_run(job, out, err);
    gd_analysis_tree_job_destroy(job);
    return rc;
}

/* ------------------------------------------------------------- date-pair */

/* ------------------------------------------------------ pairwise models */

/*
 * The one place a model name is spelled. The CLI and the browser both parse
 * through here, so they cannot come to accept different spellings, and a
 * manifest written by either names the model the same way.
 */
static const struct {
    gd_pair_model  model;
    const char    *name;      /* what a user types */
    const char    *report;    /* what a manifest records */
} PAIR_MODELS[] = {
    { GD_PAIR_P1,                  "p1",            "P1_hierarchical_constant_rate" },
    { GD_PAIR_DECAY_SWADESH,       "swadesh",       "swadesh"       },
    { GD_PAIR_DECAY_STAROSTIN,     "starostin",     "starostin"     },
    { GD_PAIR_DECAY_VASILYEV,      "vasilyev",      "vasilyev"      },
    { GD_PAIR_DECAY_TWO_COMPONENT, "two-component", "two-component" },
};

int gd_pair_model_from_string(const char *name, gd_pair_model *out, gd_error *err)
{
    if (out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG,
                            "gd_pair_model_from_string: null output");
    }
    if (name == NULL || name[0] == '\0') {
        *out = GD_PAIR_P1;
        return GD_E_OK;
    }
    for (size_t i = 0; i < sizeof PAIR_MODELS / sizeof PAIR_MODELS[0]; i++) {
        if (strcmp(name, PAIR_MODELS[i].name) == 0) {
            *out = PAIR_MODELS[i].model;
            return GD_E_OK;
        }
    }
    /* Never a silent fall back to the default: a typo in --decay would
     * otherwise report a P1 date under the name of a law nobody ran. */
    return gd_error_setf(err, GD_E_INVALID_ARG, "unknown pairwise model '%s'", name);
}

const char *gd_pair_model_name(gd_pair_model m)
{
    for (size_t i = 0; i < sizeof PAIR_MODELS / sizeof PAIR_MODELS[0]; i++) {
        if (PAIR_MODELS[i].model == m) return PAIR_MODELS[i].report;
    }
    return "unknown";
}

/* The decay laws map one-to-one onto the internal enum; GD_PAIR_P1 has no
 * decay law and never reaches here. */
static gd_decay_model pair_decay_model(gd_pair_model m)
{
    switch (m) {
        case GD_PAIR_DECAY_STAROSTIN:     return GD_DECAY_STAROSTIN;
        case GD_PAIR_DECAY_VASILYEV:      return GD_DECAY_VASILYEV;
        case GD_PAIR_DECAY_TWO_COMPONENT: return GD_DECAY_TWO_COMPONENT;
        case GD_PAIR_DECAY_SWADESH:
        case GD_PAIR_P1:
        default:                          return GD_DECAY_SWADESH;
    }
}

void gd_analysis_pair_config_default(gd_analysis_pair_config *cfg)
{
    if (cfg == NULL) return;
    memset(cfg, 0, sizeof *cfg);
    cfg->t_lo = 100.0;
    cfg->t_hi = 8000.0;
    cfg->rate_logmean = -8.9;   /* per year (about -2.0 per kiloyear) */
    cfg->rate_logsd = 0.6;
    cfg->age_a = 0.0;
    cfg->age_b = 0.0;
    cfg->loan_threshold = 0.5;
    cfg->model = GD_PAIR_P1;
}

void gd_analysis_pair_result_destroy(gd_analysis_pair_result *res)
{
    if (res == NULL) return;
    free(res->grid_t);
    free(res->grid_density);
    free(res);
}

int gd_analysis_pair_run(const gd_analysis_pair_input *in,
                         const gd_analysis_pair_config *cfg,
                         gd_analysis_pair_result **out, gd_error *err)
{
    if (in == NULL || cfg == NULL || out == NULL ||
        in->lang_a == NULL || in->lang_b == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "gd_analysis_pair_run: null argument");
    }
    *out = NULL;

    gd_forms *forms = NULL;
    int *matches = NULL;
    double *yc = NULL;
    gd_concept_ref *refs = NULL;
    gd_analysis_pair_result *r = NULL;
    int rc = GD_E_INTERNAL;
    int have_result = 0;
    gd_pair_result result;
    memset(&result, 0, sizeof result);

    {
        size_t flen = (in->forms_tsv == NULL) ? 0
                    : ((in->forms_len == GD_LEN_NUL_TERMINATED)
                        ? strlen(in->forms_tsv) : in->forms_len);
        gd_forms_options fo; gd_forms_options_default(&fo);
        fo.analysis_id = in->analysis_id;
        fo.loan_threshold = cfg->loan_threshold;
        if (gd_forms_read_opts(in->forms_tsv, flen, &fo, &forms, err) != GD_E_OK) {
            label_error(err, in->forms_name);
            goto done;
        }
    }

    size_t li_a = gd_forms_language_index(forms, in->lang_a);
    size_t li_b = gd_forms_language_index(forms, in->lang_b);
    /* A language the file never mentions is a typo, not "no evidence": the
     * run would otherwise return the prior and call it a date. */
    if (li_a == SIZE_MAX || li_b == SIZE_MAX) {
        gd_error_setf(err, GD_E_INVALID_ARG, "language '%s' not found in forms",
                     (li_a == SIZE_MAX) ? in->lang_a : in->lang_b);
        goto done;
    }

    /* Concepts present in both languages; a match is an identical cognate id. */
    size_t nc = forms->nconcept;
    matches = (int *)gd_malloc_array_checked(nc ? nc : 1, sizeof(int), err);
    refs = (gd_concept_ref *)gd_malloc_array_checked(nc ? nc : 1,
                                                     sizeof(gd_concept_ref), err);
    if (matches == NULL || refs == NULL) {
        gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
    }
    /* Synonyms are marginalised exactly (E2.3): with candidate weights w,
     * the probability that the two compared forms are cognate is
     * y = sum_i sum_j w_i w_j [class_i == class_j], entered as a soft match
     * (P1's y_c). One form per side gives y = 0 or 1, the plain case. */
    yc = (double *)gd_malloc_array_checked(nc ? nc : 1, sizeof(double), err);
    if (yc == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    size_t m = 0, nmatch = 0, nsyn = 0;
    for (size_t i = 0; i < nc; i++) {
        char *const *ca; char *const *cb; const double *wa; const double *wb;
        size_t na = gd_forms_candidates(forms, i, li_a, &ca, &wa);
        size_t nb = gd_forms_candidates(forms, i, li_b, &cb, &wb);
        if (na == 0 || nb == 0) continue;
        const double *ka, *kb;
        (void)gd_forms_candidates_conf(forms, i, li_a, &ka);
        (void)gd_forms_candidates_conf(forms, i, li_b, &kb);
        /* A shared class counts as a match only as far as both cognate
         * assignments can be trusted (cognacy_confidence, E2.5). */
        double y = 0.0;
        for (size_t x = 0; x < na; x++)
            for (size_t z = 0; z < nb; z++)
                if (strcmp(ca[x], cb[z]) == 0) y += wa[x] * wb[z] * ka[x] * kb[z];
        if (y > 1.0) y = 1.0;
        if (na > 1 || nb > 1) nsyn++;
        matches[m] = 1;
        yc[m] = y;
        nmatch += (y >= 1.0 - 1e-12) ? 1u : 0u;
        /* P1 works per kiloyear internally; the config is per year. */
        refs[m].log_rate_mean = cfg->rate_logmean + log(1000.0);
        refs[m].log_rate_sd = cfg->rate_logsd;
        m++;
    }
    if (m == 0) {
        gd_error_setf(err, GD_E_MISSING_FIELD,
                     "'%s' and '%s' share no observed concept", in->lang_a, in->lang_b);
        goto done;
    }

    /* GD_PAIR_P1 keeps the hierarchical constant-rate posterior; a decay law
     * uses the deterministic decay-law posterior (M12) instead. */
    double unit;                /* multiplier from the result's unit to years */
    int prc;
    if (!(cfg->age_a >= 0.0) || !(cfg->age_b >= 0.0) ||
        !isfinite(cfg->age_a) || !isfinite(cfg->age_b)) {
        gd_error_set(err, GD_E_INVALID_ARG, "tip ages must be finite and >= 0");
        goto done;
    }
    double age_max = (cfg->age_a > cfg->age_b) ? cfg->age_a : cfg->age_b;
    if (!(cfg->t_hi > age_max)) {
        gd_error_setf(err, GD_E_INVALID_ARG,
                      "--t-hi (%.0f) must exceed the older tip age (%.0f)",
                      cfg->t_hi, age_max);
        goto done;
    }
    if (cfg->model != GD_PAIR_P1 && (cfg->age_a > 0.0 || cfg->age_b > 0.0)) {
        gd_error_setf(err, GD_E_INVALID_ARG,
                      "decay law '%s' assumes contemporaneous lists; use P1 with tip ages",
                      gd_pair_model_name(cfg->model));
        goto done;
    }
    if (cfg->model != GD_PAIR_P1) {
        gd_decay_params dp; gd_decay_params_default(&dp);
        prc = gd_decay_posterior(pair_decay_model(cfg->model), &dp, matches, yc,
                                 m, cfg->t_lo, cfg->t_hi, 1024, &result, err);
        unit = 1.0;             /* already in years */
    } else {
        /* d = 2T with T in kyr: P1 works in the rate's time unit, so the bounds
         * go in as kyr and the summaries are rescaled to years below. */
        prc = gd_p1_posterior(refs, m, matches, yc,
                              cfg->age_a / 1000.0, cfg->age_b / 1000.0,
                              cfg->t_lo / 1000.0, cfg->t_hi / 1000.0,
                              0.0, 0.3, 0.0, 1.0, 19.0,
                              GD_SENSITIVITY_STANDARD, 512, 24, &result, err);
        unit = 1000.0;
    }
    if (prc != 0) goto done;
    have_result = 1;

    r = (gd_analysis_pair_result *)gd_calloc_checked(1, sizeof *r, err);
    if (r == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }
    r->model = gd_pair_model_name(cfg->model);
    r->model_id = cfg->model;
    r->n_shared = m;
    r->n_cognate = nmatch;
    r->retention = m ? (double)nmatch / (double)m : 0.0;
    r->n_missing_cognates = forms->n_missing_cognates;
    r->n_synonym_concepts = nsyn;
    r->mean = result.mean * unit;
    r->median = result.median * unit;
    r->sd = result.sd * unit;
    r->et_lo_90 = result.et_lo_90 * unit;
    r->et_hi_90 = result.et_hi_90 * unit;
    r->absolute_anchored = result.absolute_anchored;
    /* The separation is an increasing affine map of T, so its quantiles are
     * the mapped quantiles of T. */
    r->age_a = cfg->age_a;
    r->age_b = cfg->age_b;
    r->sep_median   = 2.0 * r->median   - cfg->age_a - cfg->age_b;
    r->sep_et_lo_90 = 2.0 * r->et_lo_90 - cfg->age_a - cfg->age_b;
    r->sep_et_hi_90 = 2.0 * r->et_hi_90 - cfg->age_a - cfg->age_b;
    r->prior_dominance = result.prior_dominance;

    /* Carry the posterior grid out in years, rather than discarding it. */
    if (result.grid_T != NULL && result.grid_density != NULL && result.grid_n > 0) {
        r->grid_t = (double *)gd_malloc_array_checked(result.grid_n, sizeof(double), err);
        r->grid_density = (double *)gd_malloc_array_checked(result.grid_n, sizeof(double), err);
        if (r->grid_t == NULL || r->grid_density == NULL) {
            gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done;
        }
        for (size_t i = 0; i < result.grid_n; i++) {
            r->grid_t[i] = result.grid_T[i] * unit;
            /* A density in the result's time unit; rescaling the axis by
             * `unit` divides the density by it, so the curve still integrates
             * to one over years. */
            r->grid_density[i] = result.grid_density[i] / unit;
        }
        r->ngrid = result.grid_n;
    }

    *out = r;  r = NULL;
    rc = GD_E_OK;

done:
    if (rc != GD_E_OK) rc = fail_code(err);
    if (r != NULL) { free(r->grid_t); free(r->grid_density); free(r); }
    if (have_result) {
        free((void *)result.grid_T);
        free((void *)result.grid_density);
    }
    free(matches); free(refs); free(yc);
    gd_forms_destroy(forms);
    return rc;
}

int gd_analysis_clade_draws(const gd_analysis_tree_result *res, const char *members,
                            double *draws, size_t cap, size_t *ndraws,
                            bool *monophyletic, gd_error *err)
{
    if (res == NULL || members == NULL || draws == NULL || ndraws == NULL || monophyletic == NULL)
        return gd_error_set(err, GD_E_INVALID_ARG, "clade_draws: NULL");
    size_t mlen = strlen(members);
    char *target = (char *)gd_malloc_checked(mlen + 6, err);
    if (target == NULL) return GD_E_NOMEM;
    memcpy(target, "mrca:", 5);
    memcpy(target + 5, members, mlen + 1);
    size_t v = GD_NODE_NONE;
    int rc = gd_tree_resolve_target(res->tree, target, &v, err);
    free(target);
    if (rc != 0) return rc;
    /* Members, counted; the MRCA's tips, counted. Equal counts mean the
     * clade is exactly the MRCA's tip set (all members are under it). */
    size_t nmem = 1;
    for (const char *p = members; *p; p++) if (*p == ',') nmem++;
    size_t under = 0;
    for (size_t u = 0; u < res->tree->nnodes; u++) {
        if (!res->tree->nodes[u].is_leaf) continue;
        size_t w = u;
        while (w != GD_NODE_NONE && w != v) w = res->tree->nodes[w].parent;
        if (w == v) under++;
    }
    *monophyletic = (under == nmem);
    size_t pidx = SIZE_MAX;
    for (size_t p = 0; p < res->nparams; p++)
        if (res->params[p].type == GD_MCMC_P_AGE && res->params[p].node == v) pidx = p;
    size_t n = 0;
    if (pidx == SIZE_MAX) {
        for (; n < cap; n++) draws[n] = res->tree->nodes[v].age;
    } else {
        if (res->trace == NULL || res->trace_draws == 0)
            return gd_error_set(err, GD_E_INVALID_ARG, "clade_draws: the run kept no trace (trace_draws = 0)");
        for (size_t c = 0; c < res->n_chains && n < cap; c++)
            for (size_t k = 0; k < res->trace_draws && n < cap; k++)
                draws[n++] = gd_analysis_trace_draw(res, pidx, c, k);
    }
    *ndraws = n;
    return 0;
}

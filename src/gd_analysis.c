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

/* ------------------------------------------------- calibration variation */

/* Loosen a soft window about its centre. A fixed tip or a hard bound has no
 * soft window to loosen and is left exactly as it was — widening a stated fact
 * would be inventing doubt rather than measuring it. */
static void widen_window(gd_calibration *cal, double f)
{
    switch (cal->type) {
    case GD_CAL_UNIFORM_NODE: {
        double c = 0.5 * (cal->p1 + cal->p2);
        double h = 0.5 * (cal->p2 - cal->p1) * f;
        double lo = c - h;
        cal->p1 = (lo < 0.0) ? 0.0 : lo;
        cal->p2 = c + h;
        break;
    }
    case GD_CAL_NORMAL_NODE:
    case GD_CAL_LOGNORMAL_NODE:
    case GD_CAL_SOFT_MINIMUM:
    case GD_CAL_SOFT_MAXIMUM:
        cal->p2 *= f;   /* the sd (log-sd for the lognormal) */
        break;
    case GD_CAL_OFFSET_LOGNORMAL_NODE:
        /* p1 is the hard offset and p2 the log-mean: the spread is p3.
         * Multiplying p2 used to move the distribution instead of widening
         * it (by a factor exp(p2 (f - 1)) in its median). */
        cal->p3 *= f;
        break;
    default:
        break;
    }
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
            bool intrinsic_hard = (c->type == GD_CAL_FIXED_TIP ||
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
    cfg->n_warmup = 2000;
    cfg->n_samples = 2000;
    cfg->n_chains = 4;
    cfg->seed = 42;
    cfg->rate_logmean = -8.80;
    cfg->rate_logsd = 0.80;
    cfg->epsilon = 0.05;
    cfg->relaxed = false;
    cfg->block_age_scaling = true;
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
} cognate_build;

static const char *const TYPE_NAMES[] = {
    "log_rate", "age", "log_sigma_b", "logit_eps", "branch_logmult",
    "concept_dev"
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
        if (gd_forms_read(in->forms_tsv, flen, &forms, err) != GD_E_OK) {
            label_error(err, in->forms_name);
            goto done;
        }
    }
    if (forms->nconcept == 0) {
        gd_error_set(err, GD_E_MISSING_FIELD, "forms contains no data rows");
        goto done;
    }
    job->cb.nconcept = forms->nconcept;
    job->n_unreadable_weights = forms->n_unreadable_weights;
    job->n_missing_cognates = forms->n_missing_cognates;

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
        drop_tip = (tv != GD_NODE_NONE && gd_tree_is_leaf(job->tree, tv));
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
    job->opt.incompatible = GD_INCOMPATIBLE_CONTAMINATION;   /* absorb homoplasy */
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

void gd_analysis_tree_result_destroy(gd_analysis_tree_result *res)
{
    if (res == NULL) return;
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

    gd_mcmc_result *res = NULL;
    double *col = NULL;
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
    mc.block_age_scaling = job->cfg.block_age_scaling;
    mc.on_progress = job->cfg.on_progress;
    mc.progress_user = job->cfg.progress_user;
    mc.progress_every = job->cfg.progress_every;

    job->ran = 1;
    res = gd_mcmc_run(&model, job->tree, job->cs, &mc, job->cfg.n_chains, err);
    if (res == NULL) goto done;

    size_t n = job->tree->nnodes;
    /* res->n_chains * res->n_draws was already checked when the sampler
     * allocated its sample matrix from the same two numbers. */
    size_t total = res->n_chains * res->n_draws, np = res->n_param;
    col = (double *)gd_malloc_array_checked(total ? total : 1, sizeof(double), err);
    r = (gd_analysis_tree_result *)gd_calloc_checked(1, sizeof *r, err);
    if (col == NULL || r == NULL) { gd_error_set(err, GD_E_NOMEM, "out of memory"); goto done; }

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

    /* The tree moves into the result; tear down what still points at it first,
     * so the job and the result can be destroyed in either order. */
    job_release_tree_users(job);
    r->tree = job->tree;
    job->tree = NULL;

    *out = r;  r = NULL;
    rc = GD_E_OK;

done:
    if (rc != GD_E_OK) rc = fail_code(err);
    free(col);
    if (r != NULL) { free(r->trace); free(r->ages); free(r->params); free(r); }
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
        if (gd_forms_read(in->forms_tsv, flen, &forms, err) != GD_E_OK) {
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
    size_t m = 0, nmatch = 0;
    for (size_t i = 0; i < nc; i++) {
        const char *ca = (li_a == SIZE_MAX) ? NULL : gd_forms_primary(forms, i, li_a);
        const char *cb = (li_b == SIZE_MAX) ? NULL : gd_forms_primary(forms, i, li_b);
        if (ca == NULL || cb == NULL) continue;
        matches[m] = (strcmp(ca, cb) == 0) ? 1 : 0;
        nmatch += (size_t)matches[m];
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
        prc = gd_decay_posterior(pair_decay_model(cfg->model), &dp, matches, NULL,
                                 m, cfg->t_lo, cfg->t_hi, 1024, &result, err);
        unit = 1.0;             /* already in years */
    } else {
        /* d = 2T with T in kyr: P1 works in the rate's time unit, so the bounds
         * go in as kyr and the summaries are rescaled to years below. */
        prc = gd_p1_posterior(refs, m, matches, NULL,
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
    free(matches); free(refs);
    gd_forms_destroy(forms);
    return rc;
}

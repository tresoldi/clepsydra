/*
 * gd_report.c — tables of reported values, and the two renderers over them.
 *
 * The formats are reproduced exactly as the hand-rolled writers had them,
 * because the committed example outputs are byte-for-byte reproducible and a
 * moved decimal would be a change in the evidence, not in the code.
 */
#include "gd_report.h"
#include "gd_memory.h"
#include "gd_math.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef enum { CELL_EMPTY = 0, CELL_TEXT, CELL_NUM, CELL_UINT } cell_kind;

typedef struct {
    cell_kind kind;
    double    num;
    uint64_t  uint;
    char     *text;      /* owned */
} cell;

struct gd_table {
    const gd_column *cols;
    size_t ncols;
    cell  *cells;        /* nrows * ncols */
    size_t nrows, cap;
};

int gd_table_create(const gd_column *cols, size_t ncols,
                    gd_table **out, gd_error *err)
{
    if (cols == NULL || out == NULL || ncols == 0) {
        return gd_error_set(err, GD_E_INVALID_ARG, "report: bad argument");
    }
    gd_table *t = (gd_table *)gd_calloc_checked(1, sizeof *t, err);
    if (t == NULL) return GD_E_NOMEM;
    t->cols = cols;
    t->ncols = ncols;
    *out = t;
    return GD_E_OK;
}

void gd_table_destroy(gd_table *t)
{
    if (t == NULL) return;
    if (t->cells != NULL) {
        size_t n = t->nrows * t->ncols;
        for (size_t i = 0; i < n; i++) gd_free(t->cells[i].text);
        gd_free(t->cells);
    }
    gd_free(t);
}

int gd_table_add_row(gd_table *t, gd_error *err)
{
    if (t == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "report: null table");
    if (t->nrows == t->cap) {
        size_t ncap = t->cap ? t->cap * 2 : 16;
        size_t ncell;
        if (!gd_size_mul(ncap, t->ncols, &ncell)) {
            return gd_error_set(err, GD_E_NOMEM, "report: table too large");
        }
        cell *grown = (cell *)gd_realloc_array_checked(t->cells, ncell,
                                                       sizeof(cell), err);
        if (grown == NULL) return GD_E_NOMEM;
        t->cells = grown;
        t->cap = ncap;
    }
    memset(&t->cells[t->nrows * t->ncols], 0, t->ncols * sizeof(cell));
    t->nrows++;
    return GD_E_OK;
}

static cell *cell_at(gd_table *t, size_t col)
{
    if (t == NULL || t->nrows == 0 || col >= t->ncols) return NULL;
    return &t->cells[(t->nrows - 1) * t->ncols + col];
}

int gd_table_set_text(gd_table *t, size_t col, const char *s, gd_error *err)
{
    cell *c = cell_at(t, col);
    if (c == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "report: bad cell");
    char *copy = gd_strdup_n(s ? s : "", s ? strlen(s) : 0, err);
    if (copy == NULL) return GD_E_NOMEM;
    gd_free(c->text);
    c->text = copy;
    c->kind = CELL_TEXT;
    return GD_E_OK;
}

int gd_table_set_num(gd_table *t, size_t col, double v, gd_error *err)
{
    cell *c = cell_at(t, col);
    if (c == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "report: bad cell");
    c->kind = CELL_NUM;
    c->num = v;
    return GD_E_OK;
}

int gd_table_set_uint(gd_table *t, size_t col, uint64_t v, gd_error *err)
{
    cell *c = cell_at(t, col);
    if (c == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "report: bad cell");
    c->kind = CELL_UINT;
    c->uint = v;
    return GD_E_OK;
}

size_t gd_table_nrows(const gd_table *t) { return (t != NULL) ? t->nrows : 0; }
size_t gd_table_ncols(const gd_table *t) { return (t != NULL) ? t->ncols : 0; }

static const cell *peek(const gd_table *t, size_t row, size_t col)
{
    if (t == NULL || row >= t->nrows || col >= t->ncols) return NULL;
    return &t->cells[row * t->ncols + col];
}

const char *gd_table_text_at(const gd_table *t, size_t row, size_t col)
{
    const cell *c = peek(t, row, col);
    return (c != NULL && c->kind == CELL_TEXT) ? c->text : NULL;
}

double gd_table_num_at(const gd_table *t, size_t row, size_t col)
{
    const cell *c = peek(t, row, col);
    if (c == NULL) return NAN;
    if (c->kind == CELL_NUM)  return c->num;
    if (c->kind == CELL_UINT) return (double)c->uint;
    return NAN;
}

/* ------------------------------------------------------------------ sinks */

static int emit(const gd_sink *sink, const char *s, size_t len)
{
    if (sink == NULL || sink->write == NULL) return GD_E_INVALID_ARG;
    if (len == 0) return GD_E_OK;
    return sink->write(sink->ctx, s, len);
}

static int emit_str(const gd_sink *sink, const char *s)
{
    return emit(sink, s, strlen(s));
}

/* ---------------------------------------------------------- TSV rendering */

int gd_table_write_tsv(const gd_table *t, const gd_sink *sink, gd_error *err)
{
    if (t == NULL || sink == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "report: bad argument");
    }
    int rc;
    int first = 1;
    for (size_t c = 0; c < t->ncols; c++) {
        if (t->cols[c].tsv_name == NULL) continue;
        if (!first && (rc = emit_str(sink, "\t")) != GD_E_OK) return rc;
        if ((rc = emit_str(sink, t->cols[c].tsv_name)) != GD_E_OK) return rc;
        first = 0;
    }
    if ((rc = emit_str(sink, "\n")) != GD_E_OK) return rc;

    for (size_t r = 0; r < t->nrows; r++) {
        first = 1;
        for (size_t c = 0; c < t->ncols; c++) {
            if (t->cols[c].tsv_name == NULL) continue;
            if (!first && (rc = emit_str(sink, "\t")) != GD_E_OK) return rc;
            first = 0;
            const cell *cl = &t->cells[r * t->ncols + c];
            char buf[64];
            switch (cl->kind) {
            case CELL_TEXT:
                if ((rc = emit_str(sink, cl->text ? cl->text : "")) != GD_E_OK) return rc;
                break;
            case CELL_NUM:
                if (t->cols[c].tsv_style == GD_NUM_GENERAL) {
                    snprintf(buf, sizeof buf, "%.*g", t->cols[c].tsv_precision, cl->num);
                } else {
                    snprintf(buf, sizeof buf, "%.*f", t->cols[c].tsv_precision, cl->num);
                }
                if ((rc = emit_str(sink, buf)) != GD_E_OK) return rc;
                break;
            case CELL_UINT:
                snprintf(buf, sizeof buf, "%llu", (unsigned long long)cl->uint);
                if ((rc = emit_str(sink, buf)) != GD_E_OK) return rc;
                break;
            case CELL_EMPTY:
                break;
            }
        }
        if ((rc = emit_str(sink, "\n")) != GD_E_OK) return rc;
    }
    return GD_E_OK;
}

/* --------------------------------------------------------- JSON rendering */

/* Same escaping and same non-finite handling the shim's emitter used: a NaN
 * or an infinity becomes null, so JSON.parse cannot throw on a degenerate
 * run. */
static int emit_json_string(const gd_sink *sink, const char *s)
{
    int rc;
    if ((rc = emit_str(sink, "\"")) != GD_E_OK) return rc;
    for (const char *p = s ? s : ""; *p != '\0'; p++) {
        unsigned char ch = (unsigned char)*p;
        char esc[8];
        const char *out = NULL;
        size_t len = 1;
        switch (ch) {
        case '"':  out = "\\\""; len = 2; break;
        case '\\': out = "\\\\"; len = 2; break;
        case '\n': out = "\\n";  len = 2; break;
        case '\r': out = "\\r";  len = 2; break;
        case '\t': out = "\\t";  len = 2; break;
        default:
            if (ch < 0x20) {
                snprintf(esc, sizeof esc, "\\u%04x", ch);
                out = esc; len = 6;
            } else {
                out = p; len = 1;
            }
        }
        if ((rc = emit(sink, out, len)) != GD_E_OK) return rc;
    }
    return emit_str(sink, "\"");
}

int gd_table_write_json(const gd_table *t, const gd_sink *sink, gd_error *err)
{
    if (t == NULL || sink == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "report: bad argument");
    }
    int rc;
    if ((rc = emit_str(sink, "[")) != GD_E_OK) return rc;
    for (size_t r = 0; r < t->nrows; r++) {
        if (r > 0 && (rc = emit_str(sink, ",")) != GD_E_OK) return rc;
        if ((rc = emit_str(sink, "{")) != GD_E_OK) return rc;
        int first = 1;
        for (size_t c = 0; c < t->ncols; c++) {
            if (t->cols[c].json_name == NULL) continue;
            if (!first && (rc = emit_str(sink, ",")) != GD_E_OK) return rc;
            first = 0;
            if ((rc = emit_json_string(sink, t->cols[c].json_name)) != GD_E_OK) return rc;
            if ((rc = emit_str(sink, ":")) != GD_E_OK) return rc;
            const cell *cl = &t->cells[r * t->ncols + c];
            char buf[64];
            switch (cl->kind) {
            case CELL_TEXT:
                if ((rc = emit_json_string(sink, cl->text)) != GD_E_OK) return rc;
                break;
            case CELL_NUM:
                if (!gd_is_finite(cl->num)) {
                    if ((rc = emit_str(sink, "null")) != GD_E_OK) return rc;
                } else {
                    snprintf(buf, sizeof buf, "%.10g", cl->num);
                    if ((rc = emit_str(sink, buf)) != GD_E_OK) return rc;
                }
                break;
            case CELL_UINT:
                snprintf(buf, sizeof buf, "%llu", (unsigned long long)cl->uint);
                if ((rc = emit_str(sink, buf)) != GD_E_OK) return rc;
                break;
            case CELL_EMPTY:
                if ((rc = emit_str(sink, "null")) != GD_E_OK) return rc;
                break;
            }
        }
        if ((rc = emit_str(sink, "}")) != GD_E_OK) return rc;
    }
    return emit_str(sink, "]");
}

/* ---------------------------------------------------------- the reports */

const char *gd_report_tree_model_name(const gd_analysis_tree_result *res)
{
    return (res != NULL && res->relaxed) ? "T1_novel_state_relaxed"
                                         : "T1_novel_state";
}

/* posterior_nodes.tsv, and the browser's `nodes` array. The TSV's "node"
 * column carries the label; the JSON reports the label and the node index
 * separately, because the page draws a topology and needs the index. */
static const gd_column NODE_COLS[] = {
    { "node",         "label",   0, GD_NUM_FIXED },
    { NULL,           "node",    0, GD_NUM_FIXED },
    { "mean_bp",      "mean",    1, GD_NUM_FIXED },
    { "median_bp",    "median",  1, GD_NUM_FIXED },
    { "sd",           "sd",      1, GD_NUM_FIXED },
    { "hpd90_lo_bp",  "hpd_lo",  1, GD_NUM_FIXED },
    { "hpd90_hi_bp",  "hpd_hi",  1, GD_NUM_FIXED },
    { "rhat",         "rhat",    3, GD_NUM_FIXED },
    { "ess",          "ess",     0, GD_NUM_FIXED },
};

int gd_report_nodes(const gd_analysis_tree_result *res,
                    gd_table **out, gd_error *err)
{
    if (res == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "report: null argument");
    }
    gd_table *t = NULL;
    int rc = gd_table_create(NODE_COLS, sizeof NODE_COLS / sizeof NODE_COLS[0],
                             &t, err);
    if (rc != GD_E_OK) return rc;

    for (size_t i = 0; i < res->nages; i++) {
        const gd_analysis_node_age *a = &res->ages[i];
        if ((rc = gd_table_add_row(t, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_text(t, 0, a->label ? a->label : "(node)", err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_uint(t, 1, (uint64_t)a->node, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 2, a->age.mean, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 3, a->age.median, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 4, a->age.sd, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 5, a->age.hpd_lo, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 6, a->age.hpd_hi, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 7, a->rhat, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 8, a->ess, err)) != GD_E_OK) goto fail;
    }
    *out = t;
    return GD_E_OK;
fail:
    gd_table_destroy(t);
    return rc;
}

/* diagnostics.tsv, and the browser's `params` array. */
static const gd_column PARAM_COLS[] = {
    { "param",      NULL,         0, GD_NUM_FIXED },
    { "type",       "type",       0, GD_NUM_FIXED },
    { NULL,         "label",      0, GD_NUM_FIXED },
    /* Significant digits, not decimals: these summarise parameters whose
     * magnitudes differ by orders (a log-rate near -9, an age near 2670). */
    { "mean",       "mean",       4, GD_NUM_GENERAL },
    { "sd",         "sd",         4, GD_NUM_GENERAL },
    { "rhat",       "rhat",       3, GD_NUM_FIXED },
    { "ess",        "ess",        0, GD_NUM_FIXED },
    { "acceptance", "acceptance", 3, GD_NUM_FIXED },
};

int gd_report_diagnostics(const gd_analysis_tree_result *res,
                          gd_table **out, gd_error *err)
{
    if (res == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "report: null argument");
    }
    gd_table *t = NULL;
    int rc = gd_table_create(PARAM_COLS, sizeof PARAM_COLS / sizeof PARAM_COLS[0],
                             &t, err);
    if (rc != GD_E_OK) return rc;

    for (size_t p = 0; p < res->nparams; p++) {
        const gd_analysis_param *pp = &res->params[p];
        /* Sized for a node label plus the "age:" prefix. A truncated name is
         * not cosmetic: two long labels can truncate to the same string, and
         * the table would then show one parameter measured twice instead of
         * two parameters. Refuse rather than write that. */
        char pname[512];
        int pn;
        if (pp->type == GD_MCMC_P_AGE && pp->label != NULL)
            pn = snprintf(pname, sizeof pname, "age:%s", pp->label);
        else
            pn = snprintf(pname, sizeof pname, "%s", pp->type_name);
        if (pn < 0 || (size_t)pn >= sizeof pname) {
            rc = gd_error_setf(err, GD_E_RANGE,
                    "report: parameter name is too long to report (%d bytes)", pn);
            goto fail;
        }
        if ((rc = gd_table_add_row(t, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_text(t, 0, pname, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_text(t, 1, pp->type_name, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_text(t, 2, pp->label, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 3, pp->stat.mean, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 4, pp->stat.sd, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 5, pp->rhat, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 6, pp->ess, err)) != GD_E_OK) goto fail;
        if ((rc = gd_table_set_num(t, 7, pp->acceptance, err)) != GD_E_OK) goto fail;
    }
    *out = t;
    return GD_E_OK;
fail:
    gd_table_destroy(t);
    return rc;
}

/* posterior_pair.tsv: a key/value table, so both columns are text and the
 * numbers are pre-formatted to the precision each line has always used. */
static const gd_column KV_COLS[] = {
    { "key",   "key",   0, GD_NUM_FIXED },
    { "value", "value", 0, GD_NUM_FIXED },
};

static int kv(gd_table *t, const char *key, const char *val, gd_error *err)
{
    int rc = gd_table_add_row(t, err);
    if (rc != GD_E_OK) return rc;
    if ((rc = gd_table_set_text(t, 0, key, err)) != GD_E_OK) return rc;
    return gd_table_set_text(t, 1, val, err);
}

static int kv_num(gd_table *t, const char *key, double v, int prec, gd_error *err)
{
    char buf[64];
    snprintf(buf, sizeof buf, "%.*f", prec, v);
    return kv(t, key, buf, err);
}

static int kv_uint(gd_table *t, const char *key, uint64_t v, gd_error *err)
{
    char buf[32];
    snprintf(buf, sizeof buf, "%llu", (unsigned long long)v);
    return kv(t, key, buf, err);
}

int gd_report_pair(const gd_analysis_pair_result *res,
                   const char *lang_a, const char *lang_b,
                   gd_table **out, gd_error *err)
{
    if (res == NULL || out == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "report: null argument");
    }
    gd_table *t = NULL;
    int rc = gd_table_create(KV_COLS, sizeof KV_COLS / sizeof KV_COLS[0], &t, err);
    if (rc != GD_E_OK) return rc;

    if ((rc = kv(t, "model", res->model, err)) != GD_E_OK) goto fail;
    if ((rc = kv(t, "language_a", lang_a, err)) != GD_E_OK) goto fail;
    if ((rc = kv(t, "language_b", lang_b, err)) != GD_E_OK) goto fail;
    if ((rc = kv_uint(t, "n_shared_concepts", (uint64_t)res->n_shared, err)) != GD_E_OK) goto fail;
    if ((rc = kv_uint(t, "n_cognate_matches", (uint64_t)res->n_cognate, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "retention_fraction", res->retention, 4, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "divergence_mean_bp", res->mean, 1, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "divergence_median_bp", res->median, 1, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "divergence_sd_bp", res->sd, 1, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "et90_lo_bp", res->et_lo_90, 1, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "et90_hi_bp", res->et_hi_90, 1, err)) != GD_E_OK) goto fail;
    if ((rc = kv_uint(t, "absolute_anchored", res->absolute_anchored ? 1u : 0u, err)) != GD_E_OK) goto fail;
    if ((rc = kv_num(t, "prior_dominance", res->prior_dominance, 3, err)) != GD_E_OK) goto fail;

    *out = t;
    return GD_E_OK;
fail:
    gd_table_destroy(t);
    return rc;
}

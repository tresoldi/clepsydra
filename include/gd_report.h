/*
 * gd_report.h — the reported tables, described once.
 *
 * A chronology comes out as posterior_nodes.tsv and diagnostics.tsv on the
 * command line, and as `nodes` and `params` arrays in the browser's JSON.
 * Before D5 those were six independent serializers: five hand-rolled
 * out_path/fopen/fprintf blocks in cli/commands.c and a JSON emitter in the
 * shim, none of them callable from a test. The column names lived inside
 * command bodies, and tools/validation/calibration_benchmark.py read them back
 * by name with nothing connecting the two.
 *
 * A table here holds the values; a renderer turns it into bytes. The library
 * still does no I/O: bytes go to a caller-supplied sink, which is a FILE* on
 * the command line, a JSON buffer in the shim, and memory in the tests.
 *
 * TWO NAMES PER COLUMN, ON PURPOSE. The TSV and the JSON are different
 * schemas and always were: the TSV's `node` column holds a label while the
 * JSON has both a numeric `node` and a `label`, and the TSV rounds to a fixed
 * number of decimals where the JSON keeps %.10g. Rather than unify them --
 * which would change a published payload for no gain -- each column states
 * both spellings, so the correspondence is written down in one place for the
 * first time. A NULL name means the column does not appear in that format.
 */
#ifndef GD_REPORT_H
#define GD_REPORT_H

#include "gd_types.h"
#include "gd_error.h"
#include "gd_analysis.h"

#include <stddef.h>

/* Where rendered bytes go. `write` returns GD_E_OK or a GD_E_* code; a
 * failing sink stops the render and the code comes back to the caller. */
typedef struct gd_sink {
    void *ctx;
    int (*write)(void *ctx, const char *buf, size_t len);
} gd_sink;

typedef struct gd_table gd_table;

/* How the TSV renders a number. Both styles are in use and they are not
 * interchangeable: an age at %.1f keeps a tenth of a year, while a posterior
 * mean at %.4g keeps four significant digits whatever its magnitude. */
typedef enum gd_num_style {
    GD_NUM_FIXED = 0,     /* printf %.*f */
    GD_NUM_GENERAL        /* printf %.*g */
} gd_num_style;

/* One column: its spelling in each format, and how the TSV formats it.
 * `tsv_precision` and `tsv_style` are ignored for text cells. */
typedef struct gd_column {
    const char  *tsv_name;    /* NULL: omitted from the TSV */
    const char  *json_name;   /* NULL: omitted from the JSON */
    int          tsv_precision;
    gd_num_style tsv_style;
} gd_column;

/* Build a table. `cols` is borrowed and must outlive the table -- every
 * caller passes a static array. */
int  gd_table_create(const gd_column *cols, size_t ncols,
                     gd_table **out, gd_error *err);
void gd_table_destroy(gd_table *t);

/* Append a row, then fill its cells. A cell never set renders as empty. */
int gd_table_add_row(gd_table *t, gd_error *err);
int gd_table_set_text(gd_table *t, size_t col, const char *s, gd_error *err);
int gd_table_set_num(gd_table *t, size_t col, double v, gd_error *err);
int gd_table_set_uint(gd_table *t, size_t col, uint64_t v, gd_error *err);

/* Read back, which is what a test asserts on: no parsing of rendered text. */
size_t      gd_table_nrows(const gd_table *t);
size_t      gd_table_ncols(const gd_table *t);
const char *gd_table_text_at(const gd_table *t, size_t row, size_t col);
double      gd_table_num_at(const gd_table *t, size_t row, size_t col);

/* Renderers. */
int gd_table_write_tsv(const gd_table *t, const gd_sink *sink, gd_error *err);
int gd_table_write_json(const gd_table *t, const gd_sink *sink, gd_error *err);

/* ------------------------------------------------------------- the reports */

int gd_report_nodes(const gd_analysis_tree_result *res,
                    gd_table **out, gd_error *err);
int gd_report_diagnostics(const gd_analysis_tree_result *res,
                          gd_table **out, gd_error *err);
/* posterior_calibrations.tsv: each calibration on a sampled node, with the
 * target's posterior and, for an outlier-prone one, its outlier probability. */
int gd_report_calibrations(const gd_analysis_tree_result *res,
                           gd_table **out, gd_error *err);
/* adequacy.tsv: the posterior predictive checks. */
int gd_report_adequacy(const gd_analysis_tree_result *res,
                       gd_table **out, gd_error *err);
/* posterior_concepts.tsv: one row per concept (empty for a prior-only run). */
int gd_report_concepts(const gd_analysis_tree_result *res,
                       gd_table **out, gd_error *err);
int gd_report_pair(const gd_analysis_pair_result *res,
                   const char *lang_a, const char *lang_b,
                   gd_table **out, gd_error *err);

/*
 * The model a chronology result was produced by, as both surfaces name it.
 * Borrowed static storage. This was two independent string literals, in
 * cli/commands.c and web/clepsydra_wasm.c, with nothing comparing them.
 */
const char *gd_report_tree_model_name(const gd_analysis_tree_result *res);

#endif /* GD_REPORT_H */

/*
 * gd_forms.c — one reading of the forms table.
 *
 * Streams the TSV: the parser hands back one row at a time with field
 * pointers into its own buffer, and this module interns the language and the
 * concept, then records the row. Winners are resolved in a second pass over
 * the recorded rows, once the table's shape is known.
 *
 * Recording rows rather than growing a concept x language matrix in place is
 * what keeps the two-dimensional growth out of the hot path: neither
 * dimension is known until the last row has been seen.
 */
#include "gd_forms.h"
#include "gd_tsv.h"
#include "gd_memory.h"

#include <stdlib.h>
#include <string.h>

/* One data row, reduced to what choosing a winner needs. */
typedef struct {
    size_t concept;
    size_t lang;
    double weight;
    char  *cognate;      /* owned */
} forms_row;

typedef struct {
    size_t c_lang, c_con, c_cog, c_w;   /* column indices; SIZE_MAX if absent */
    int    header_seen;

    char **concepts;  size_t nconcept, concept_cap;
    char **languages; size_t nlang, lang_cap;
    size_t *rows_per_lang;

    forms_row *rows;  size_t nrows, row_cap;
    size_t n_unreadable_weights;

    int failed;                          /* a GD_E_* code, or 0 */
} forms_builder;

/*
 * A field's effective length. Parser fields are raw bytes and may contain an
 * embedded NUL; the string-table loader this replaced duplicated each cell
 * with gd_strdup_n and compared it with strcmp, which truncates there. Doing
 * the same keeps the reading identical -- and the first version of this file,
 * which used strncmp plus s[len], read past the end of a shorter stored
 * string whenever a field carried an embedded NUL. The fuzzer found it.
 */
static size_t field_len(const char *field, size_t len)
{
    const char *nul = (const char *)memchr(field, '\0', len);
    return (nul != NULL) ? (size_t)(nul - field) : len;
}

/* Compare a field against a C string, both truncated at any embedded NUL. */
static int field_is(const char *field, size_t len, const char *s)
{
    size_t fl = field_len(field, len);
    return strlen(s) == fl && memcmp(s, field, fl) == 0;
}

/* Intern into a growable string table, preserving first-seen order.
 * Returns the index, or SIZE_MAX on allocation failure. */
static size_t intern(char ***tab, size_t *n, size_t *cap,
                     const char *field, size_t len)
{
    for (size_t i = 0; i < *n; i++) {
        if (field_is(field, len, (*tab)[i])) return i;
    }
    if (*n == *cap) {
        size_t nc = *cap ? *cap * 2 : 64;
        char **grown = (char **)gd_realloc_array_checked(*tab, nc,
                                                         sizeof(char *), NULL);
        if (grown == NULL) return SIZE_MAX;
        *tab = grown;
        *cap = nc;
    }
    char *s = gd_strdup_n(field, field_len(field, len), NULL);
    if (s == NULL) return SIZE_MAX;
    (*tab)[*n] = s;
    return (*n)++;
}

static int header_cb(const gd_tsv_row *row, void *user, gd_error *err)
{
    forms_builder *b = (forms_builder *)user;
    b->c_lang = b->c_con = b->c_cog = b->c_w = SIZE_MAX;
    for (size_t i = 0; i < row->nfields; i++) {
        const char *f = row->fields[i];
        size_t      l = row->lengths[i];
        if      (field_is(f, l, "language_id"))    b->c_lang = i;
        else if (field_is(f, l, "concept_id"))     b->c_con  = i;
        else if (field_is(f, l, "cognate_id"))     b->c_cog  = i;
        else if (field_is(f, l, "primary_weight")) b->c_w    = i;
    }
    b->header_seen = 1;
    if (b->c_lang == SIZE_MAX || b->c_con == SIZE_MAX || b->c_cog == SIZE_MAX) {
        b->failed = GD_E_UNKNOWN_COLUMN;
        return gd_error_set(err, GD_E_UNKNOWN_COLUMN,
                "forms needs the columns language_id, concept_id, cognate_id");
    }
    return GD_E_OK;
}

static int row_cb(const gd_tsv_row *row, void *user, gd_error *err)
{
    forms_builder *b = (forms_builder *)user;

    /* A short row is not an error: the missing fields read as empty, which is
     * what the string-table loader this replaced did too. */
    const char *lf = (b->c_lang < row->nfields) ? row->fields[b->c_lang] : "";
    size_t      ll = (b->c_lang < row->nfields) ? row->lengths[b->c_lang] : 0;
    const char *cf = (b->c_con  < row->nfields) ? row->fields[b->c_con]  : "";
    size_t      cl = (b->c_con  < row->nfields) ? row->lengths[b->c_con]  : 0;
    const char *gf = (b->c_cog  < row->nfields) ? row->fields[b->c_cog]  : "";
    size_t      gl = (b->c_cog  < row->nfields) ? row->lengths[b->c_cog]  : 0;

    size_t li = intern(&b->languages, &b->nlang, &b->lang_cap, lf, ll);
    if (li == SIZE_MAX) { b->failed = GD_E_NOMEM; return gd_error_set(err, GD_E_NOMEM, "out of memory"); }
    size_t ci = intern(&b->concepts, &b->nconcept, &b->concept_cap, cf, cl);
    if (ci == SIZE_MAX) { b->failed = GD_E_NOMEM; return gd_error_set(err, GD_E_NOMEM, "out of memory"); }

    /* Weight. Absent column, absent cell and '.' all mean 1.0, and so does a
     * value that will not parse -- counted, not refused. Never atof: it
     * answers 0.0 for anything it cannot read, and a weight that silently
     * becomes zero drops a form from the comparison without saying so. */
    double w = 1.0;
    if (b->c_w != SIZE_MAX && b->c_w < row->nfields) {
        const char *wf = row->fields[b->c_w];
        size_t      wl = row->lengths[b->c_w];
        if (wl != 0 && !gd_tsv_is_missing(wf, wl)) {
            double parsed;
            if (gd_tsv_parse_double(wf, wl, &parsed) == GD_E_OK) {
                w = parsed;
            } else {
                b->n_unreadable_weights++;
            }
        }
    }

    if (b->nrows == b->row_cap) {
        size_t nc = b->row_cap ? b->row_cap * 2 : 256;
        forms_row *grown = (forms_row *)gd_realloc_array_checked(b->rows, nc,
                                                    sizeof(forms_row), NULL);
        if (grown == NULL) { b->failed = GD_E_NOMEM; return gd_error_set(err, GD_E_NOMEM, "out of memory"); }
        b->rows = grown;
        b->row_cap = nc;
    }
    forms_row *fr = &b->rows[b->nrows];
    fr->concept = ci;
    fr->lang    = li;
    fr->weight  = w;
    fr->cognate = gd_strdup_n(gf, field_len(gf, gl), NULL);
    if (fr->cognate == NULL) { b->failed = GD_E_NOMEM; return gd_error_set(err, GD_E_NOMEM, "out of memory"); }
    b->nrows++;
    return GD_E_OK;
}

static void builder_release(forms_builder *b)
{
    for (size_t i = 0; i < b->nrows; i++) gd_free(b->rows[i].cognate);
    gd_free(b->rows);
    for (size_t i = 0; i < b->nconcept; i++) gd_free(b->concepts[i]);
    gd_free(b->concepts);
    for (size_t i = 0; i < b->nlang; i++) gd_free(b->languages[i]);
    gd_free(b->languages);
    gd_free(b->rows_per_lang);
    memset(b, 0, sizeof *b);
}

void gd_forms_destroy(gd_forms *f)
{
    if (f == NULL) return;
    /* Each array is checked against its own pointer, not only its count: this
     * runs on the partially-built object of a failed read, where a count can
     * already be set while the array it counts is still NULL.
     * The assignments own the cognate strings; primary[] only points at them. */
    if (f->assignments != NULL) {
        for (size_t i = 0; i < f->nassignments; i++) {
            gd_free((void *)f->assignments[i].cognate);
        }
    }
    gd_free(f->assignments);
    gd_free((void *)f->primary);
    if (f->concepts != NULL) {
        for (size_t i = 0; i < f->nconcept; i++) gd_free(f->concepts[i]);
    }
    gd_free(f->concepts);
    if (f->languages != NULL) {
        for (size_t i = 0; i < f->nlang; i++) gd_free(f->languages[i]);
    }
    gd_free(f->languages);
    gd_free(f->rows_per_lang);
    gd_free(f);
}

int gd_forms_read(const char *tsv, size_t len,
                  gd_forms **out, gd_error *err)
{
    int rc = GD_E_INTERNAL;
    forms_builder b;
    gd_forms *f = NULL;
    double *best = NULL;

    if (out == NULL || tsv == NULL) {
        return gd_error_set(err, GD_E_INVALID_ARG, "forms: missing input buffer");
    }

    memset(&b, 0, sizeof b);
    b.c_lang = b.c_con = b.c_cog = b.c_w = SIZE_MAX;

    rc = gd_tsv_parse(tsv, len, header_cb, row_cb, &b, err);
    if (rc != GD_E_OK) {
        if (b.failed != 0) rc = b.failed;
        goto done;
    }

    b.rows_per_lang = (size_t *)gd_calloc_checked(b.nlang ? b.nlang : 1,
                                                  sizeof(size_t), err);
    if (b.rows_per_lang == NULL) { rc = GD_E_NOMEM; goto done; }
    for (size_t i = 0; i < b.nrows; i++) b.rows_per_lang[b.rows[i].lang]++;

    f = (gd_forms *)gd_calloc_checked(1, sizeof *f, err);
    if (f == NULL) { rc = GD_E_NOMEM; goto done; }

    /* concepts x languages: both counts come from the input, so the product
     * is checked before it becomes an allocation size. */
    size_t ncell;
    if (!gd_size_mul(b.nconcept, b.nlang, &ncell)) {
        rc = gd_error_set(err, GD_E_NOMEM, "forms: concept x language count overflows");
        goto done;
    }
    f->primary = (const char **)gd_calloc_checked(ncell ? ncell : 1, sizeof(char *), err);
    best = (double *)gd_malloc_array_checked(ncell ? ncell : 1, sizeof(double), err);
    if (f->primary == NULL || best == NULL) { rc = GD_E_NOMEM; goto done; }
    for (size_t i = 0; i < ncell; i++) best[i] = -1.0;

    f->assignments = (gd_forms_assignment *)gd_malloc_array_checked(
                            b.nrows ? b.nrows : 1, sizeof(gd_forms_assignment), err);
    if (f->assignments == NULL) { rc = GD_E_NOMEM; goto done; }

    /* Highest weight wins; the first row wins a tie, which is what both of
     * the intakes this replaces did. Every row that takes over its cell is
     * recorded in order, including ones a later row overrides. */
    for (size_t i = 0; i < b.nrows; i++) {
        forms_row *fr = &b.rows[i];
        size_t slot = fr->concept * b.nlang + fr->lang;
        if (fr->weight <= best[slot]) continue;
        best[slot] = fr->weight;

        gd_forms_assignment *a = &f->assignments[f->nassignments++];
        a->concept = fr->concept;
        a->lang    = fr->lang;
        a->cognate = fr->cognate;   /* ownership moves to the gd_forms */
        fr->cognate = NULL;
        f->primary[slot] = a->cognate;
    }

    /* Hand the interned tables over; the builder must not free them now. */
    f->nconcept      = b.nconcept;  b.nconcept = 0;
    f->nlang         = b.nlang;     b.nlang = 0;
    f->concepts      = b.concepts;  b.concepts = NULL;
    f->languages     = b.languages; b.languages = NULL;
    f->rows_per_lang = b.rows_per_lang; b.rows_per_lang = NULL;
    f->n_unreadable_weights = b.n_unreadable_weights;

    *out = f;
    f = NULL;
    rc = GD_E_OK;

done:
    gd_free(best);
    gd_forms_destroy(f);
    builder_release(&b);
    return rc;
}

size_t gd_forms_language_index(const gd_forms *f, const char *name)
{
    if (f == NULL || name == NULL) return SIZE_MAX;
    for (size_t i = 0; i < f->nlang; i++)
        if (strcmp(f->languages[i], name) == 0) return i;
    return SIZE_MAX;
}

const char *gd_forms_primary(const gd_forms *f, size_t concept, size_t lang)
{
    if (f == NULL || concept >= f->nconcept || lang >= f->nlang) return NULL;
    return f->primary[concept * f->nlang + lang];
}

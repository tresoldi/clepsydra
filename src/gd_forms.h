/*
 * gd_forms.h — one reading of the forms table.
 *
 * The forms TSV answers one question per (language, concept): which cognate
 * id represents that meaning in that language? When a language has synonyms,
 * `primary_weight` picks between them — highest weight wins, and the first
 * row wins a tie.
 *
 * That rule used to be written three times: once for the chronology intake,
 * once for the pairwise intake, and once in tools/analysis/concept_stability.py.
 * The two C copies had drifted apart on what to do with a weight that will not
 * parse -- the chronology ignored it, the pairwise pipeline failed the run --
 * so the same file could be read two ways in the same process. This module is
 * the one reading, and the disagreement is settled: an unreadable weight falls
 * back to 1.0 everywhere.
 *
 * That direction was chosen because the weight only chooses between synonyms:
 * the concept is still usable, and refusing the whole run over one cell of a
 * large file is a harsh answer to a typo. The information is not thrown away
 * -- n_unreadable_weights counts the cells it happened in, so a caller that
 * wants to be loud about it still can. Note what this reverses: the pairwise
 * pipeline used to refuse, deliberately (see 2284c61), because the atof() that
 * preceded it answered 0.0 for anything unparseable and dropped the form
 * silently. Falling back to 1.0 is not that bug -- the form is kept at full
 * weight, not dropped -- but a typo is quiet again unless someone reads the
 * count.
 *
 * What this module deliberately does NOT do is intern cognate ids to integers.
 * The chronology needs per-concept integer classes and the pairwise model
 * needs to compare strings; interning here would have to pick one, and the
 * numbering a caller sees would then depend on rows it does not care about.
 */
#ifndef GD_FORMS_H
#define GD_FORMS_H

#include "gd_types.h"
#include "gd_error.h"

#include <stddef.h>

/*
 * A forms table, read once.
 *
 * `concepts` and `languages` are in first-seen row order, which is
 * load-bearing: the committed example outputs depend on concept order.
 * Both include entries seen anywhere in the file, including languages that
 * no tree mentions — that is what the chronology intake counted rows against.
 *
 * `primary[c * nlang + l]` is the winning cognate id for that cell, or NULL
 * when the language has no form for that concept. Borrowed from this struct
 * and freed with it.
 */
/*
 * One row that took over its cell: at the moment it was read, its weight beat
 * everything seen for that (concept, language) so far.
 *
 * The sequence matters, not just the final winners. A caller that numbers
 * cognate classes itself -- the chronology intake does, because T1 wants
 * per-concept integers -- must assign those numbers in the order the file
 * produced assignments, including assignments later overridden. T1's noise
 * term sums lgamma over label-sorted blocks, so a different numbering is a
 * different summation order, and in a sampler that is a different chain.
 */
typedef struct gd_forms_assignment {
    size_t      concept;
    size_t      lang;
    const char *cognate;      /* borrowed from the gd_forms */
} gd_forms_assignment;

typedef struct gd_forms {
    size_t   nconcept;
    size_t   nlang;
    char   **concepts;        /* nconcept, owned */
    char   **languages;       /* nlang, owned */
    size_t  *rows_per_lang;   /* nlang: data rows the file gave each language */
    const char **primary;     /* nconcept * nlang; borrowed, NULL where absent */

    gd_forms_assignment *assignments;   /* in file order */
    size_t               nassignments;

    /* Cells whose primary_weight was present but unreadable, and so was taken
     * as 1.0. Zero for a clean file. */
    size_t n_unreadable_weights;

    /* Rows that carry no observation (empty, '.' or '?' cognate id, or
     * status missing/excluded/absent). They are not cognate classes and are
     * not counted in rows_per_lang. */
    size_t n_missing_cognates;
} gd_forms;

/*
 * Read a forms TSV. `len` is a byte count; 0 means an empty buffer.
 * Requires the columns language_id, concept_id and cognate_id;
 * primary_weight is optional and defaults to 1.0. An optional status column
 * marks rows as missing, excluded or absent; those rows, and rows whose
 * cognate id is empty, '.' or '?', are skipped and counted.
 *
 * On success *out is an owned gd_forms; on failure *out is untouched and this
 * function owns nothing.
 */
int gd_forms_read(const char *tsv, size_t len,
                  gd_forms **out, gd_error *err);

void gd_forms_destroy(gd_forms *f);

/* Index of a language by name, or SIZE_MAX if the file never mentioned it. */
size_t gd_forms_language_index(const gd_forms *f, const char *name);

/* The winning cognate id for a cell, or NULL. Borrowed from `f`. */
const char *gd_forms_primary(const gd_forms *f, size_t concept, size_t lang);

#endif /* GD_FORMS_H */

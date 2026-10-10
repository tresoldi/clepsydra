/*
 * gd_treeset.h — a sample of trees, as written by MCMC programs (roadmap
 * E3.4).
 *
 * Accepts plain Newick (one tree per ';'), or NEXUS with a trees block:
 * `tree NAME = [&R] (...);` lines, an optional `Translate` table mapping
 * tokens (usually numbers) to taxon names, and `[...]` comments anywhere.
 * BEAST, MrBayes and RevBayes all write the latter. Every tree comes back as
 * a plain Newick string with taxon names restored and comments removed.
 */
#ifndef GD_TREESET_H
#define GD_TREESET_H

#include "gd_error.h"

#include <stddef.h>

typedef struct gd_treeset {
    char  **newick;     /* ntrees owned, NUL-terminated strings */
    size_t  ntrees;
} gd_treeset;

int  gd_treeset_read(const char *buf, size_t len, gd_treeset **out, gd_error *err);

/* The tip names of a Newick tree, comma-separated, in tree order. *out is
 * owned by the caller (free()). */
int  gd_newick_tip_list(const char *newick, char **out, gd_error *err);
void gd_treeset_destroy(gd_treeset *ts);

#endif /* GD_TREESET_H */

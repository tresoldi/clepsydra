/*
 * gd_treeset.c — reading tree samples. See gd_treeset.h.
 */
#include "gd_treeset.h"
#include "gd_memory.h"
#include "gd_newick.h"
#include "gd_tree.h"
#include "gd_tree_build.h"

#include <ctype.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct { char *from, *to; } tr_entry;

static bool ci_prefix(const char *s, const char *end, const char *word)
{
    size_t n = strlen(word);
    if ((size_t)(end - s) < n) return false;
    for (size_t i = 0; i < n; i++)
        if (tolower((unsigned char)s[i]) != tolower((unsigned char)word[i])) return false;
    return true;
}

static char *dup_range(const char *a, const char *b)
{
    size_t n = (size_t)(b - a);
    char *s = (char *)gd_malloc_checked(n + 1, NULL);
    if (s == NULL) return NULL;
    memcpy(s, a, n);
    s[n] = '\0';
    return s;
}

/* Newick text with [...] comments removed and leaf tokens translated. */
static char *clean_tree(const char *a, const char *b, const tr_entry *tr, size_t ntr)
{
    size_t cap = (size_t)(b - a) * 2 + 64, n = 0;
    char *out = (char *)gd_malloc_checked(cap, NULL);
    if (out == NULL) return NULL;
    char prev = '(';
    const char *p = a;
    while (p < b) {
        if (*p == '[') {                          /* comment */
            while (p < b && *p != ']') p++;
            if (p < b) p++;
            continue;
        }
        if (isspace((unsigned char)*p)) { p++; continue; }
        const char *tok = p;
        if (strchr("(),:;", *p) == NULL) {
            while (p < b && strchr("(),:;[", *p) == NULL && !isspace((unsigned char)*p)) p++;
            const char *rep = NULL;
            size_t tl = (size_t)(p - tok);
            if (prev == '(' || prev == ',') {      /* a leaf name position */
                for (size_t i = 0; i < ntr; i++)
                    if (strlen(tr[i].from) == tl && memcmp(tr[i].from, tok, tl) == 0) { rep = tr[i].to; break; }
            }
            const char *src = rep ? rep : tok;
            size_t sl = rep ? strlen(rep) : tl;
            if (n + sl + 2 >= cap) {
                cap = (n + sl + 2) * 2;
                char *g = (char *)gd_realloc_array_checked(out, cap, 1, NULL);
                if (g == NULL) { free(out); return NULL; }
                out = g;
            }
            memcpy(out + n, src, sl);
            n += sl;
            prev = 'x';
            continue;
        }
        if (n + 2 >= cap) {
            cap *= 2;
            char *g = (char *)gd_realloc_array_checked(out, cap, 1, NULL);
            if (g == NULL) { free(out); return NULL; }
            out = g;
        }
        out[n++] = *p;
        prev = *p;
        p++;
    }
    if (n == 0 || out[n - 1] != ';') out[n++] = ';';
    out[n] = '\0';
    return out;
}

static int push_tree(gd_treeset *ts, size_t *cap, char *s, gd_error *err)
{
    if (s == NULL) return gd_error_set(err, GD_E_NOMEM, "treeset: out of memory");
    if (ts->ntrees == *cap) {
        size_t nc = *cap ? *cap * 2 : 64;
        char **g = (char **)gd_realloc_array_checked(ts->newick, nc, sizeof(char *), NULL);
        if (g == NULL) { free(s); return gd_error_set(err, GD_E_NOMEM, "treeset: out of memory"); }
        ts->newick = g;
        *cap = nc;
    }
    ts->newick[ts->ntrees++] = s;
    return 0;
}

int gd_treeset_read(const char *buf, size_t len, gd_treeset **out, gd_error *err)
{
    if (buf == NULL || out == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "treeset: NULL");
    gd_treeset *ts = (gd_treeset *)gd_calloc_checked(1, sizeof *ts, err);
    if (ts == NULL) return GD_E_NOMEM;
    tr_entry *tr = NULL;
    size_t ntr = 0, trcap = 0, cap = 0;
    int rc = 0;
    const char *end = buf + len, *p = buf;
    while (p < end && isspace((unsigned char)*p)) p++;

    if (ci_prefix(p, end, "#nexus")) {
        /* Translate table, if any: "translate a b, c d, ... ;" */
        const char *q = p;
        while (q < end) {
            if (ci_prefix(q, end, "translate") && (q == buf || isspace((unsigned char)q[-1]))) {
                q += 9;
                while (q < end && *q != ';') {
                    while (q < end && (isspace((unsigned char)*q) || *q == ',')) q++;
                    if (q >= end || *q == ';') break;
                    const char *a0 = q;
                    while (q < end && !isspace((unsigned char)*q)) q++;
                    const char *a1 = q;
                    while (q < end && isspace((unsigned char)*q)) q++;
                    const char *b0 = q;
                    while (q < end && *q != ',' && *q != ';' && !isspace((unsigned char)*q)) q++;
                    const char *b1 = q;
                    if (b1 > b0 && (b1[-1] == '\'' || b1[-1] == '"')) b1--;
                    if (b0 < b1 && (*b0 == '\'' || *b0 == '"')) b0++;
                    if (ntr == trcap) {
                        size_t nc = trcap ? trcap * 2 : 32;
                        tr_entry *g = (tr_entry *)gd_realloc_array_checked(tr, nc, sizeof *tr, NULL);
                        if (g == NULL) { rc = gd_error_set(err, GD_E_NOMEM, "treeset: out of memory"); goto done; }
                        tr = g; trcap = nc;
                    }
                    tr[ntr].from = dup_range(a0, a1);
                    tr[ntr].to = dup_range(b0, b1);
                    if (tr[ntr].from == NULL || tr[ntr].to == NULL) {
                        free(tr[ntr].from); free(tr[ntr].to);
                        rc = gd_error_set(err, GD_E_NOMEM, "treeset: out of memory"); goto done;
                    }
                    ntr++;
                }
                break;
            }
            q++;
        }
        /* Tree lines: "tree NAME = ... ;" */
        q = p;
        while (q < end) {
            const char *line = q;
            while (q < end && *q != '\n') q++;
            const char *s = line;
            while (s < q && isspace((unsigned char)*s)) s++;
            if (ci_prefix(s, q, "tree ") || ci_prefix(s, q, "tree\t")) {
                const char *eq = memchr(s, '=', (size_t)(q - s));
                if (eq != NULL) {
                    const char *t0 = eq + 1, *t1 = t0;
                    while (t1 < end && *t1 != ';') t1++;
                    if (t1 < end) t1++;
                    if ((rc = push_tree(ts, &cap, clean_tree(t0, t1, tr, ntr), err)) != 0) goto done;
                    q = t1;
                }
            }
            if (q < end) q++;
        }
    } else {
        const char *a = p;
        for (const char *q = p; q < end; q++) {
            if (*q == '[') { while (q < end && *q != ']') q++; continue; }
            if (*q != ';') continue;
            const char *s = a;
            while (s < q && isspace((unsigned char)*s)) s++;
            if (s < q && (rc = push_tree(ts, &cap, clean_tree(s, q + 1, NULL, 0), err)) != 0) goto done;
            a = q + 1;
        }
    }
    if (ts->ntrees == 0) { rc = gd_error_set(err, GD_E_PARSE, "treeset: no trees found"); goto done; }

done:
    for (size_t i = 0; i < ntr; i++) { free(tr[i].from); free(tr[i].to); }
    free(tr);
    if (rc != 0) { gd_treeset_destroy(ts); return rc; }
    *out = ts;
    return 0;
}

void gd_treeset_destroy(gd_treeset *ts)
{
    if (ts == NULL) return;
    for (size_t i = 0; i < ts->ntrees; i++) free(ts->newick[i]);
    free(ts->newick);
    free(ts);
}

int gd_newick_tip_list(const char *newick, char **out, gd_error *err)
{
    if (newick == NULL || out == NULL) return gd_error_set(err, GD_E_INVALID_ARG, "tip_list: NULL");
    gd_tree *t = NULL;
    int rc = gd_newick_parse(newick, strlen(newick), &t, err);
    if (rc != 0) return rc;
    size_t need = 1;
    for (size_t v = 0; v < t->nnodes; v++)
        if (t->nodes[v].is_leaf && t->nodes[v].label) need += strlen(t->nodes[v].label) + 1;
    char *s = (char *)gd_malloc_checked(need, err);
    if (s == NULL) { gd_tree_destroy(t); return GD_E_NOMEM; }
    size_t k = 0;
    for (size_t v = 0; v < t->nnodes; v++) {
        if (!t->nodes[v].is_leaf || !t->nodes[v].label) continue;
        if (k) s[k++] = ',';
        size_t l = strlen(t->nodes[v].label);
        memcpy(s + k, t->nodes[v].label, l);
        k += l;
    }
    s[k] = '\0';
    gd_tree_destroy(t);
    *out = s;
    return 0;
}

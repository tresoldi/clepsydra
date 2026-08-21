/*
 * clepsydra_wasm.c — the browser's view of the engine: strings in, JSON out.
 *
 * Every entry point returns a caller-owned string that JavaScript must hand
 * back to clepsydra_free. The analysis itself is src/gd_analysis.c, the same
 * code the native CLI runs, so a date computed here and a date computed on the
 * command line come from one implementation (PLAN.md M10).
 *
 * Links with -sFILESYSTEM=0: the engine reads only the buffers it is given.
 *
 * Options arrive as `key=value` lines rather than JSON. Every option this
 * engine takes is a number, so a key/value list needs no quoting, no escaping
 * and no parser worth testing; the strings that could need escaping (language
 * names, the decay law) are separate arguments. An unknown key is an error
 * rather than a silent default — a misspelt `sample=4000` that quietly ran 2000
 * would be a bad way to learn about typos.
 */
#include "gd.h"
#include "gd_options.h"
#include "gd_analysis.h"

#include <math.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

EMSCRIPTEN_KEEPALIVE const char *clepsydra_version(void);
EMSCRIPTEN_KEEPALIVE const char *clepsydra_git_commit(void);
EMSCRIPTEN_KEEPALIVE char *clepsydra_date_tree_json(const char *newick,
                                                    const char *forms_tsv,
                                                    const char *calibrations_tsv,
                                                    const char *opts);
EMSCRIPTEN_KEEPALIVE char *clepsydra_date_pair_json(const char *forms_tsv,
                                                    const char *lang_a,
                                                    const char *lang_b,
                                                    const char *decay,
                                                    const char *opts);
EMSCRIPTEN_KEEPALIVE char *clepsydra_defaults_json(void);
EMSCRIPTEN_KEEPALIVE void clepsydra_free(char *text);

/* ---------------------------------------------------- progress into JS ---- */

#ifdef __EMSCRIPTEN__
/* EM_JS's documented spelling ends in a semicolon, which expands to a stray one
 * at file scope; ISO C forbids that and -Wpedantic says so. The macro is
 * emscripten's, so the suppression is scoped to this one construct. */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wextra-semi"
EM_JS(int, clepsydra_js_progress, (int chain, int iter, int total), {
    if (typeof Module !== "undefined" && typeof Module.onProgress === "function") {
        return Module.onProgress(chain, iter, total) ? 1 : 0;
    }
    return 0;
});
#pragma clang diagnostic pop
#else
static int clepsydra_js_progress(int chain, int iter, int total)
{
    (void)chain; (void)iter; (void)total;
    return 0;
}
#endif

/* Returning nonzero cancels the run; the caller sees ok:false with
 * status "cancelled" rather than a partial chronology. */
static int progress_bridge(void *user, size_t chain, size_t iter, size_t total)
{
    (void)user;
    return clepsydra_js_progress((int)chain, (int)iter, (int)total);
}

/* ------------------------------------------------------ JSON out buffer ---- */

typedef struct { char *data; size_t len, cap; int failed; } jbuf;

static void jb_init(jbuf *b)
{
    b->cap = 1024; b->len = 0; b->failed = 0;
    b->data = (char *)malloc(b->cap);
    if (b->data == NULL) b->failed = 1; else b->data[0] = '\0';
}

static void jb_reserve(jbuf *b, size_t need)
{
    if (b->failed) return;
    if (b->len + need + 1 <= b->cap) return;
    size_t nc = b->cap * 2;
    while (nc < b->len + need + 1) nc *= 2;
    char *nd = (char *)realloc(b->data, nc);
    if (nd == NULL) { free(b->data); b->data = NULL; b->failed = 1; return; }
    b->data = nd; b->cap = nc;
}

static void jb_puts(jbuf *b, const char *s)
{
    size_t n = strlen(s);
    jb_reserve(b, n);
    if (b->failed) return;
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
}

/* Append one byte. A two-element string just to add a character is both
 * roundabout and, as gcc's analyzer read it, indistinguishable from passing an
 * uninitialised buffer. */
static void jb_putc(jbuf *b, char ch)
{
    jb_reserve(b, 1);
    if (b->failed) return;
    b->data[b->len++] = ch;
    b->data[b->len] = '\0';
}

static void jb_printf(jbuf *b, const char *fmt, ...)
{
    /* Initialised before the call, not for correctness — vsnprintf terminates
     * within the buffer whenever it returns >= 0 — but because gcc's analyzer
     * cannot model that and would otherwise read tmp as uninitialised here.
     * Costing one store on an error-free path is cheaper than a suppression
     * that would also hide a real uninitialised read later. */
    char tmp[512] = { '\0' };
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(tmp, sizeof tmp, fmt, ap);
    va_end(ap);
    /* Truncation is a failure, not something to emit: a clipped number would
     * be valid JSON carrying the wrong value, which is worse than no JSON. */
    if (n < 0 || (size_t)n >= sizeof tmp) { b->failed = 1; return; }
    jb_puts(b, tmp);
}

/* JSON has no NaN or Infinity: those become null, which JavaScript can test
 * for, rather than a token that would make JSON.parse throw. */
static void jb_num(jbuf *b, double v)
{
    if (isnan(v) || isinf(v)) jb_puts(b, "null");
    else jb_printf(b, "%.10g", v);
}

static void jb_str(jbuf *b, const char *s)
{
    if (s == NULL) { jb_puts(b, "null"); return; }
    jb_puts(b, "\"");
    for (const unsigned char *p = (const unsigned char *)s; *p != '\0'; p++) {
        switch (*p) {
        case '"':  jb_puts(b, "\\\""); break;
        case '\\': jb_puts(b, "\\\\"); break;
        case '\n': jb_puts(b, "\\n"); break;
        case '\r': jb_puts(b, "\\r"); break;
        case '\t': jb_puts(b, "\\t"); break;
        default:
            if (*p < 0x20) jb_printf(b, "\\u%04x", (unsigned)*p);
            else jb_putc(b, (char)*p);
        }
    }
    jb_puts(b, "\"");
}

/* Hands the buffer to the caller, or a minimal error payload if it never
 * assembled. Always returns something JSON.parse can read. */
static char *jb_finish(jbuf *b)
{
    if (b->failed || b->data == NULL) {
        free(b->data);
        char *fallback = (char *)malloc(64);
        if (fallback != NULL)
            snprintf(fallback, 64, "{\"ok\":false,\"error\":\"out of memory\"}");
        return fallback;
    }
    return b->data;
}

static char *json_error(const char *status, const char *detail, int code)
{
    jbuf b; jb_init(&b);
    jb_puts(&b, "{\"ok\":false,\"error\":");
    jb_str(&b, status);
    jb_puts(&b, ",\"detail\":");
    jb_str(&b, detail);
    jb_printf(&b, ",\"code\":%d}", code);
    return jb_finish(&b);
}

/* --------------------------------------------------------- option lines ---- */

/*
 * Parses `key=value` lines against the option registry in gd_options.h.
 *
 * The table this replaced held its targets as void *, so nothing checked that
 * a slot's declared kind matched the field it pointed at. It was wrong once:
 * two config fields became bool while their slots still said int, and the
 * parser wrote four bytes into a one-byte field, corrupting whatever followed.
 * No compiler saw it; scripts/page-smoke.js did, by noticing the draws had
 * changed. The registry states each field's type next to its offset and
 * asserts the two agree at compile time, so that class of bug is now a build
 * error rather than a browser-only symptom.
 *
 * Blank lines are allowed; anything else -- an unknown key, a missing '=', a
 * value that is not a number or is out of range -- is an error naming the
 * offending text. Returns 0 on success.
 */
static int parse_opts(const char *opts, const gd_option *reg, size_t nreg,
                      void *cfg, char *errbuf, size_t errbuf_size)
{
    if (opts == NULL) return 0;
    const char *p = opts;
    while (*p != '\0') {
        while (*p == '\n' || *p == '\r' || *p == ' ' || *p == '\t') p++;
        if (*p == '\0') break;
        const char *eol = p;
        while (*eol != '\0' && *eol != '\n' && *eol != '\r') eol++;
        const char *eq = p;
        while (eq < eol && *eq != '=') eq++;
        if (eq == eol) {
            snprintf(errbuf, errbuf_size, "option line without '=': %.*s",
                     (int)(eol - p), p);
            return 1;
        }
        size_t klen = (size_t)(eq - p);
        const gd_option *opt = gd_option_find(reg, nreg, p, klen);
        if (opt == NULL) {
            snprintf(errbuf, errbuf_size, "unknown option '%.*s'", (int)klen, p);
            return 1;
        }
        size_t vlen = (size_t)(eol - eq - 1);
        if (vlen == 0) {
            snprintf(errbuf, errbuf_size, "bad value for '%s'", opt->name);
            return 1;
        }
        gd_error err = GD_ERROR_INIT;
        if (gd_option_set(opt, cfg, eq + 1, vlen, &err) != GD_E_OK) {
            snprintf(errbuf, errbuf_size, "%s", err.message);
            return 1;
        }
        p = eol;
    }
    return 0;
}

/*
 * The defaults, as JSON, so the page can populate its controls from the engine
 * instead of restating numbers that live in gd_analysis_*_config_default.
 * Values are formatted by the same code that parses them back.
 */
EMSCRIPTEN_KEEPALIVE
char *clepsydra_defaults_json(void)
{
    gd_analysis_tree_config tcfg;
    gd_analysis_pair_config pcfg;
    gd_analysis_tree_config_default(&tcfg);
    gd_analysis_pair_config_default(&pcfg);

    jbuf jb;
    jb_init(&jb);
    jb_putc(&jb, '{');
    jb_str(&jb, "ok"); jb_putc(&jb, ':'); jb_puts(&jb, "true");

    static const struct { const char *key; int tree; } groups[2] = {
        { "tree", 1 }, { "pair", 0 }
    };
    for (size_t g = 0; g < 2; g++) {
        size_t n = 0;
        const gd_option *reg = groups[g].tree ? gd_options_tree(&n)
                                              : gd_options_pair(&n);
        const void *cfg = groups[g].tree ? (const void *)&tcfg
                                         : (const void *)&pcfg;
        jb_putc(&jb, ',');
        jb_str(&jb, groups[g].key);
        jb_putc(&jb, ':'); jb_putc(&jb, '{');
        for (size_t i = 0; i < n; i++) {
            char val[64];
            gd_error err = GD_ERROR_INIT;
            if (gd_option_format(&reg[i], cfg, val, sizeof val, &err) != GD_E_OK)
                continue;
            if (i > 0) jb_putc(&jb, ',');
            jb_str(&jb, reg[i].name);
            jb_putc(&jb, ':');
            jb_puts(&jb, val);
        }
        jb_putc(&jb, '}');
    }
    jb_putc(&jb, '}');
    return jb_finish(&jb);
}

/* ------------------------------------------------------------- date-tree ---- */

static void emit_topology(jbuf *b, const gd_tree *t)
{
    jb_puts(b, "\"tree\":{\"root\":");
    jb_printf(b, "%zu", t->root);
    jb_puts(b, ",\"nodes\":[");
    for (size_t i = 0; i < t->nnodes; i++) {
        const gd_node *nd = &t->nodes[i];
        if (i > 0) jb_puts(b, ",");
        jb_printf(b, "{\"id\":%zu,\"parent\":", i);
        if (nd->parent == GD_NODE_NONE) jb_puts(b, "null");
        else jb_printf(b, "%zu", nd->parent);
        jb_puts(b, ",\"label\":");
        jb_str(b, nd->label);
        jb_printf(b, ",\"is_leaf\":%s,\"age\":", nd->is_leaf ? "true" : "false");
        /* The sampler's last state. For a fixed tip that is its calibrated age,
         * which is what the time axis needs; for a sampled node the page should
         * prefer the posterior summary in `nodes` below. */
        jb_num(b, nd->age);
        jb_puts(b, "}");
    }
    jb_puts(b, "]}");
}

char *clepsydra_date_tree_json(const char *newick, const char *forms_tsv,
                               const char *calibrations_tsv, const char *opts)
{
    if (newick == NULL || forms_tsv == NULL || calibrations_tsv == NULL)
        return json_error("bad request", "tree, forms and calibrations are all required",
                          GD_E_INVALID_ARG);

    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);

    size_t nreg = 0;
    const gd_option *reg = gd_options_tree(&nreg);
    /* Big enough for a whole gd_error message: gd_option_set's are the
     * longest, and truncating an error is how a report stops naming the
     * option it is about. */
    char optserr[GD_ERROR_MSG_MAX];
    if (parse_opts(opts, reg, nreg, &cfg, optserr, sizeof optserr) != 0)
        return json_error("bad option", optserr, GD_E_INVALID_ARG);

    cfg.on_progress = progress_bridge;

    /* The three buffers arrive from JavaScript through stringToNewUTF8, which
     * terminates them, so their lengths are measured rather than passed. */
    gd_analysis_tree_input in;
    memset(&in, 0, sizeof in);
    in.newick = newick;                     in.newick_len = GD_LEN_NUL_TERMINATED;
    in.forms_tsv = forms_tsv;               in.forms_len  = GD_LEN_NUL_TERMINATED;
    in.calibrations_tsv = calibrations_tsv; in.cal_len    = GD_LEN_NUL_TERMINATED;
    in.newick_name = "tree";
    in.forms_name = "forms";
    in.cal_name = "calibrations";

    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_job *job = NULL;
    if (gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0)
        return json_error("the data could not be read",
                          err.message[0] ? err.message : "?", err.code);

    gd_analysis_tree_result *res = NULL;
    int rc = gd_analysis_tree_job_run(job, &res, &err);
    gd_analysis_tree_job_destroy(job);
    if (rc != 0) {
        if (err.code == GD_E_CANCELLED)
            return json_error("cancelled", "the run was stopped before it finished",
                              err.code);
        return json_error("the sampler failed", err.message[0] ? err.message : "?",
                          err.code);
    }

    jbuf b; jb_init(&b);
    jb_puts(&b, "{\"ok\":true,\"model\":");
    jb_str(&b, res->relaxed ? "T1_novel_state_relaxed" : "T1_novel_state");
    jb_printf(&b, ",\"n_languages\":%zu,\"n_nodes\":%zu,\"n_concepts\":%zu",
              res->n_languages, res->n_nodes, res->n_concepts);
    jb_printf(&b, ",\"n_chains\":%zu,\"n_warmup\":%zu,\"n_samples\":%zu,\"n_draws\":%zu",
              res->n_chains, res->n_warmup, res->n_samples, res->n_draws);
    jb_printf(&b, ",\"n_skipped_rows\":%zu,\"n_unreadable_weights\":%zu,\"relaxed\":%s",
              res->n_skipped_rows, res->n_unreadable_weights,
              res->relaxed ? "true" : "false");
    jb_printf(&b, ",\"branch_rate_confounded\":%s",
              res->branch_rate_confounded ? "true" : "false");
    jb_puts(&b, ",\"min_branch_duration_corr\":");
    jb_num(&b, res->min_branch_duration_corr);
    jb_printf(&b, ",\"seed\":%.0f,\"rate_logmean\":", (double)cfg.seed);
    jb_num(&b, cfg.rate_logmean);
    jb_puts(&b, ",\"rate_logsd\":");
    jb_num(&b, cfg.rate_logsd);
    jb_puts(&b, ",\"epsilon\":");
    jb_num(&b, cfg.epsilon);
    jb_puts(&b, ",\"hpd_level\":");
    jb_num(&b, cfg.hpd_level);
    jb_printf(&b, ",\"cal_widen\":%.4g", cfg.calibrations.widen);

    jb_puts(&b, ",\"nodes\":[");
    for (size_t i = 0; i < res->nages; i++) {
        const gd_analysis_node_age *a = &res->ages[i];
        if (i > 0) jb_puts(&b, ",");
        jb_printf(&b, "{\"node\":%zu,\"label\":", a->node);
        jb_str(&b, a->label);
        jb_puts(&b, ",\"mean\":");   jb_num(&b, a->age.mean);
        jb_puts(&b, ",\"median\":"); jb_num(&b, a->age.median);
        jb_puts(&b, ",\"sd\":");     jb_num(&b, a->age.sd);
        jb_puts(&b, ",\"hpd_lo\":"); jb_num(&b, a->age.hpd_lo);
        jb_puts(&b, ",\"hpd_hi\":"); jb_num(&b, a->age.hpd_hi);
        jb_puts(&b, ",\"rhat\":");   jb_num(&b, a->rhat);
        jb_puts(&b, ",\"ess\":");    jb_num(&b, a->ess);
        jb_puts(&b, "}");
    }
    jb_puts(&b, "],\"params\":[");
    for (size_t p = 0; p < res->nparams; p++) {
        const gd_analysis_param *pp = &res->params[p];
        if (p > 0) jb_puts(&b, ",");
        jb_puts(&b, "{\"type\":");
        jb_str(&b, pp->type_name);
        jb_puts(&b, ",\"label\":");
        jb_str(&b, pp->label);
        jb_puts(&b, ",\"mean\":");       jb_num(&b, pp->stat.mean);
        jb_puts(&b, ",\"sd\":");         jb_num(&b, pp->stat.sd);
        jb_puts(&b, ",\"rhat\":");       jb_num(&b, pp->rhat);
        jb_puts(&b, ",\"ess\":");        jb_num(&b, pp->ess);
        jb_puts(&b, ",\"acceptance\":"); jb_num(&b, pp->acceptance);
        /* The trace belongs to the parameter, not to a parallel array beside
         * it. It used to be emitted as res->trace[p] and paired back up by
         * index in viz.js, which is a silent mislabel the moment anything
         * reorders params. */
        if (res->trace != NULL && res->trace_draws > 0) {
            jb_puts(&b, ",\"trace\":[");
            for (size_t c = 0; c < res->n_chains; c++) {
                if (c > 0) jb_puts(&b, ",");
                jb_puts(&b, "[");
                for (size_t k = 0; k < res->trace_draws; k++) {
                    if (k > 0) jb_puts(&b, ",");
                    jb_num(&b, gd_analysis_trace_draw(res, p, c, k));
                }
                jb_puts(&b, "]");
            }
            jb_puts(&b, "]");
        }
        jb_puts(&b, "}");
    }
    jb_puts(&b, "]");

    if (res->trace != NULL && res->trace_draws > 0)
        jb_printf(&b, ",\"trace_draws\":%zu", res->trace_draws);

    jb_puts(&b, ",");
    emit_topology(&b, res->tree);
    jb_puts(&b, "}");

    gd_analysis_tree_result_destroy(res);
    return jb_finish(&b);
}

/* ------------------------------------------------------------- date-pair ---- */

char *clepsydra_date_pair_json(const char *forms_tsv, const char *lang_a,
                               const char *lang_b, const char *decay,
                               const char *opts)
{
    if (forms_tsv == NULL || lang_a == NULL || lang_b == NULL)
        return json_error("bad request", "forms and two languages are required",
                          GD_E_INVALID_ARG);

    gd_analysis_pair_config cfg;
    gd_analysis_pair_config_default(&cfg);
    size_t nreg = 0;
    const gd_option *reg = gd_options_pair(&nreg);
    /* Big enough for a whole gd_error message: gd_option_set's are the
     * longest, and truncating an error is how a report stops naming the
     * option it is about. */
    char optserr[GD_ERROR_MSG_MAX];
    if (parse_opts(opts, reg, nreg, &cfg, optserr, sizeof optserr) != 0)
        return json_error("bad option", optserr, GD_E_INVALID_ARG);
    /* The registry bounds each option on its own; that the grid's two ends are
     * ordered is a relation between two of them. */
    if (!(cfg.t_hi > cfg.t_lo))
        return json_error("bad option", "t_lo must be less than t_hi", GD_E_INVALID_ARG);

    /* An empty decay string means "no law given", which is how JavaScript will
     * spell an absent option; NULL means the same, and both select P1.
     * Parsed through the library so the page and the CLI accept exactly the
     * same spellings. */
    {
        gd_error derr = GD_ERROR_INIT;
        if (gd_pair_model_from_string(decay, &cfg.model, &derr) != GD_E_OK)
            return json_error("bad option", derr.message, derr.code);
    }

    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = forms_tsv;
    in.forms_len = GD_LEN_NUL_TERMINATED;   /* from stringToNewUTF8 */
    in.forms_name = "forms";
    in.lang_a = lang_a;
    in.lang_b = lang_b;

    gd_error err = GD_ERROR_INIT;
    gd_analysis_pair_result *res = NULL;
    if (gd_analysis_pair_run(&in, &cfg, &res, &err) != 0)
        return json_error("the pair could not be dated",
                          err.message[0] ? err.message : "?", err.code);

    jbuf b; jb_init(&b);
    jb_puts(&b, "{\"ok\":true,\"model\":");
    jb_str(&b, res->model);
    jb_puts(&b, ",\"language_a\":");
    jb_str(&b, lang_a);
    jb_puts(&b, ",\"language_b\":");
    jb_str(&b, lang_b);
    jb_printf(&b, ",\"n_shared\":%zu,\"n_cognate\":%zu", res->n_shared, res->n_cognate);
    jb_puts(&b, ",\"retention\":");   jb_num(&b, res->retention);
    jb_puts(&b, ",\"mean\":");        jb_num(&b, res->mean);
    jb_puts(&b, ",\"median\":");      jb_num(&b, res->median);
    jb_puts(&b, ",\"sd\":");          jb_num(&b, res->sd);
    jb_puts(&b, ",\"et90_lo\":");     jb_num(&b, res->et_lo_90);
    jb_puts(&b, ",\"et90_hi\":");     jb_num(&b, res->et_hi_90);
    jb_printf(&b, ",\"absolute_anchored\":%s",
              res->absolute_anchored ? "true" : "false");
    jb_puts(&b, ",\"prior_dominance\":"); jb_num(&b, res->prior_dominance);

    /* The whole posterior, not just its summary: the page draws this. */
    jb_puts(&b, ",\"grid\":{\"t\":[");
    for (size_t i = 0; i < res->ngrid; i++) {
        if (i > 0) jb_puts(&b, ",");
        jb_num(&b, res->grid_t[i]);
    }
    jb_puts(&b, "],\"density\":[");
    for (size_t i = 0; i < res->ngrid; i++) {
        if (i > 0) jb_puts(&b, ",");
        jb_num(&b, res->grid_density[i]);
    }
    jb_puts(&b, "]}}");

    gd_analysis_pair_result_destroy(res);
    return jb_finish(&b);
}

/* ----------------------------------------------------------------- misc ---- */

const char *clepsydra_version(void)   { return gd_version(); }
const char *clepsydra_git_commit(void) { return gd_git_commit(); }

void clepsydra_free(char *text) { free(text); }

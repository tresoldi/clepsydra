/*
 * commands.c — CLI command handlers for Clepsydra.
 *
 * Option parsing, file I/O, input validation and TSV formatting. The analysis
 * itself lives in src/gd_analysis.c, which the WebAssembly shim calls too, so
 * the page and the command line cannot compute different numbers (PLAN.md M10).
 *
 * Everything that reaches the engine is checked here first: a CLI argument that
 * is not a number, an output path that cannot be created, a file that cannot be
 * written. The engine reports its own errors through gd_error; this file's job
 * is to make sure it is never handed nonsense in the first place.
 */
#include "commands.h"

#include "gd.h"
#include "gd_options.h"
#include "gd_report.h"
#include "gd_analysis.h"
#include "gd_treeset.h"
#include "gd_tree.h"

#include <errno.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define GD_MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define GD_MKDIR(p) mkdir((p), 0777)
#endif

/* ---------------------------------------------------- safe filesystem/parse */

/*
 * Recursively create a directory path (like `mkdir -p`) using POSIX/Win32
 * APIs. Returns 0 on success (including when the path already exists), -1
 * otherwise. This replaces an earlier system("mkdir -p '%s'"), which let an
 * output path containing a single quote inject arbitrary shell commands.
 */
static int mkdir_p(const char *path)
{
    if (path == NULL || *path == '\0') return -1;
    size_t len = strlen(path);
    char *tmp = (char *)malloc(len + 1);
    if (tmp == NULL) return -1;
    memcpy(tmp, path, len + 1);
    /* Drop a single trailing separator so we don't mkdir an empty component. */
    if (len > 1 && (tmp[len - 1] == '/' || tmp[len - 1] == '\\'))
        tmp[len - 1] = '\0';
    for (char *q = tmp + 1; *q != '\0'; q++) {
        if (*q == '/' || *q == '\\') {
            char sep = *q;
            *q = '\0';
            if (GD_MKDIR(tmp) != 0 && errno != EEXIST) { free(tmp); return -1; }
            *q = sep;
        }
    }
    if (GD_MKDIR(tmp) != 0 && errno != EEXIST) { free(tmp); return -1; }
    free(tmp);
    return 0;
}

/* Parse a non-negative count from a CLI argument. Rejects malformed text,
 * negative values, and overflow (atol/atoi silently produced 0 or huge
 * wrapped size_t values). Returns 0 on success, -1 (with a diagnostic) else. */
/*
 * Build "<dir>/<name>" into `out`, or fail loudly.
 *
 * snprintf cannot overflow, so the old unchecked calls were memory-safe — but
 * a truncated path names a *different file*, and the run would then write its
 * results somewhere the user did not ask for and report success. Refusing is
 * the only honest answer: nothing here can shorten a path the user chose.
 */
static bool out_path(char *out, size_t out_size, const char *dir, const char *name)
{
    int n = snprintf(out, out_size, "%s/%s", dir, name);
    if (n < 0 || (size_t)n >= out_size) {
        fprintf(stderr, "error: output path '%s/%s' is too long (limit %zu bytes)\n",
                dir, name, out_size - 1);
        return false;
    }
    return true;
}



static int parse_int_arg(const char *s, const char *name, int *out)
{
    errno = 0;
    char *endp = NULL;
    long v = strtol(s, &endp, 10);
    if (endp == s || *endp != '\0' || errno == ERANGE || v < INT_MIN || v > INT_MAX) {
        fprintf(stderr, "error: %s must be an integer, got '%s'\n", name, s);
        return -1;
    }
    *out = (int)v;
    return 0;
}

/* Parse a finite double from a CLI argument, rejecting junk and nan/inf. */
static int parse_finite(const char *s, const char *name, double *out)
{
    errno = 0;
    char *endp = NULL;
    double v = strtod(s, &endp);
    if (endp == s || *endp != '\0' || !isfinite(v)) {
        fprintf(stderr, "error: %s must be a finite number, got '%s'\n", name, s);
        return -1;
    }
    *out = v;
    return 0;
}

/* FNV-1a 64-bit checksum of an input's bytes, for provenance in the run
 * manifest. Writes a 16-hex-digit string. */
/* 16 hex digits and a terminator, or the word "unavailable". */
#define CLI_CHECKSUM_SIZE 24

static void buffer_checksum(const char *buf, size_t len,
                            char out[CLI_CHECKSUM_SIZE])
{
    if (buf == NULL) { snprintf(out, CLI_CHECKSUM_SIZE, "unavailable"); return; }
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < len; i++) {
        h ^= (uint64_t)(unsigned char)buf[i];
        h *= 1099511628211ULL;
    }
    snprintf(out, CLI_CHECKSUM_SIZE, "%016llx", (unsigned long long)h);
}

static const char *build_compiler(void)
{
#ifdef __VERSION__
    return __VERSION__;
#else
    return "unknown";
#endif
}

static const char *build_platform(void)
{
#if defined(_WIN32)
    return "windows";
#elif defined(__APPLE__)
    return "darwin";
#elif defined(__linux__)
    return "linux";
#else
    return "unknown";
#endif
}

/* ------------------------------------------------------------------ files */

/* Reads a whole file into a NUL-terminated buffer. The analysis runner takes
 * buffers, never paths — the browser has no filesystem — so reading is the
 * caller's job. */
static char *read_file(const char *path, size_t *len_out)
{
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return NULL;
    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return NULL; }
    long end = ftell(fp);
    if (end < 0 || fseek(fp, 0, SEEK_SET) != 0) { fclose(fp); return NULL; }
    size_t len = (size_t)end;
    char *buf = (char *)malloc(len + 1);
    if (buf == NULL) { fclose(fp); return NULL; }
    size_t got = fread(buf, 1, len, fp);
    int bad = ferror(fp);
    fclose(fp);
    /* A short read is an I/O error, not a shorter file: dating a truncated
     * forms table would silently drop its last rows. */
    if (bad || got != len) {
        fprintf(stderr, "error: short read on %s (%zu of %zu bytes)\n", path, got, len);
        free(buf);
        return NULL;
    }
    buf[got] = '\0';
    if (len_out != NULL) *len_out = got;
    return buf;
}

/* The three inputs every tree command takes, read once and passed as buffers. */
typedef struct {
    char  *tree, *forms, *cal;
    size_t tree_len, forms_len, cal_len;
    const char *tree_path, *forms_path, *cal_path;
} tree_inputs;

static void free_tree_inputs(tree_inputs *in)
{
    free(in->tree); free(in->forms); free(in->cal);
    memset(in, 0, sizeof *in);
}

static int load_tree_inputs(tree_inputs *in, const char *tree_path,
                            const char *forms_path, const char *cal_path)
{
    memset(in, 0, sizeof *in);
    in->tree_path = tree_path;
    in->forms_path = forms_path;
    in->cal_path = cal_path;
    in->tree = read_file(tree_path, &in->tree_len);
    in->forms = read_file(forms_path, &in->forms_len);
    in->cal = read_file(cal_path, &in->cal_len);
    if (in->tree == NULL || in->forms == NULL || in->cal == NULL) {
        fprintf(stderr, "error: cannot read %s\n",
                in->tree == NULL ? tree_path :
                in->forms == NULL ? forms_path : cal_path);
        /* Release whatever did load. Callers return straight out on failure
         * without calling free_tree_inputs, so a partial load must not leave
         * them owning anything: on failure this function owns nothing. */
        free_tree_inputs(in);
        return -1;
    }
    return 0;
}

static void fill_analysis_input(const tree_inputs *in, gd_analysis_tree_input *ai)
{
    memset(ai, 0, sizeof *ai);
    ai->newick = in->tree;           ai->newick_len = in->tree_len;
    ai->forms_tsv = in->forms;       ai->forms_len = in->forms_len;
    ai->calibrations_tsv = in->cal;  ai->cal_len = in->cal_len;
    ai->newick_name = in->tree_path; ai->forms_name = in->forms_path;
    ai->cal_name = in->cal_path;
}

static const char *opt_arg(int argc, char **argv, const char *key, const char *def)
{
    for (int i = 2; i + 1 < argc; i++) if (strcmp(argv[i], key) == 0) return argv[i + 1];
    return def;
}

/* The library renders into a sink; on the command line the sink is a file.
 * This is the whole of the CLI's side of the reporting seam. */
static int file_sink_write(void *ctx, const char *buf, size_t len)
{
    FILE *fp = (FILE *)ctx;
    return (fwrite(buf, 1, len, fp) == len) ? GD_E_OK : GD_E_IO;
}

/*
 * Render a table to out_dir/name. Owns the whole open-write-flush-close
 * sequence, including the truncation check on the path: a truncated path
 * names a different file, and date-pair once reported success while writing
 * to "posteri" instead of posterior_pair.tsv.
 */
static int write_table(const gd_table *t, const char *out_dir, const char *name)
{
    char path[1200];
    if (!out_path(path, sizeof path, out_dir, name)) return -1;
    FILE *fp = fopen(path, "w");
    if (fp == NULL) {
        fprintf(stderr, "error: cannot open %s for writing\n", path);
        return -1;
    }
    gd_sink sink = { fp, file_sink_write };
    gd_error err = GD_ERROR_INIT;
    int rc = gd_table_write_tsv(t, &sink, &err);
    if (rc != GD_E_OK) {
        fprintf(stderr, "error: writing %s: %s\n", path,
                err.message[0] ? err.message : "?");
        fclose(fp);
        return -1;
    }
    if (fclose(fp) != 0) {
        fprintf(stderr, "error: failed to flush %s\n", path);
        return -1;
    }
    return 0;
}

static const char *errmsg(const gd_error *err)
{
    return (err != NULL && err->message[0] != '\0') ? err->message : "?";
}

/* --progress: one rewritten stderr line per chain. Never returns nonzero — the
 * CLI has no way to be interrupted mid-run, so it never cancels; the browser
 * uses the same callback to do exactly that. */
struct cli_progress { size_t n_chains; };

static int cli_progress_cb(void *user, size_t chain, size_t iter, size_t total)
{
    const struct cli_progress *p = (const struct cli_progress *)user;
    double frac = total ? (double)iter / (double)total : 0.0;
    fprintf(stderr, "\r  chain %zu/%zu  %3.0f%%", chain + 1, p->n_chains,
            100.0 * frac);
    if (iter >= total && chain + 1 == p->n_chains) fprintf(stderr, "\n");
    fflush(stderr);
    return 0;
}

/*
 * Options that fill a config come from the registry in gd_options.h, so the
 * CLI, the WebAssembly shim and the page all accept the same set with the
 * same bounds. Command-shape flags that are not config fields -- --tree,
 * --out, --mode, --decay -- stay hand-written below, because there is nothing
 * for them to be a row of.
 */
static int apply_registry_opts(int argc, char **argv,
                               const gd_option *opts, size_t nopt, void *cfg)
{
    gd_error err = GD_ERROR_INIT;
    for (size_t i = 0; i < nopt; i++) {
        const char *v = opt_arg(argc, argv, opts[i].flag, NULL);
        if (v == NULL) continue;    /* not given: keep the library default */
        if (gd_option_set(&opts[i], cfg, v, GD_OPT_MEASURE, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            return -1;
        }
    }
    return 0;
}

/*
 * Check the shape of the argument list before any of it is used. Every one of
 * these used to be silent:
 *   - an unknown flag was ignored and the default used;
 *   - a known flag in last position, or followed by another flag, had no value,
 *     so it was skipped and the default used ("--samples --chains 4" ran with
 *     the default sample count);
 *   - a flag given twice took its first value and dropped the second;
 *   - a stray word that is not a flag's value was ignored.
 * Values may start with a single '-' (a negative log-rate), never with "--".
 */
static int reject_unknown_opts(int argc, char **argv,
                               const gd_option *opts, size_t nopt,
                               const char *const *extra, size_t nextra)
{
    for (int i = 2; i < argc; i++) {
        if (argv[i][0] != '-' || argv[i][1] != '-') {
            fprintf(stderr, "error: unexpected argument '%s'\n", argv[i]);
            return -1;
        }
        bool known = (gd_option_find_flag(opts, nopt, argv[i]) != NULL);
        for (size_t k = 0; !known && k < nextra; k++)
            if (strcmp(argv[i], extra[k]) == 0) known = true;
        if (!known) {
            fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
            return -1;
        }
        if (i + 1 >= argc || (argv[i + 1][0] == '-' && argv[i + 1][1] == '-')) {
            fprintf(stderr, "error: option '%s' needs a value\n", argv[i]);
            return -1;
        }
        for (int j = 2; j < i; j++) {
            if (strcmp(argv[j], argv[i]) == 0) {
                fprintf(stderr, "error: option '%s' given more than once\n", argv[i]);
                return -1;
            }
        }
        i++;
    }
    return 0;
}

/* One "  --flag VALUE (default)  help" line per row, with the default read
 * back out of a defaulted config rather than retyped. */
static void print_options(const gd_option *opts, size_t nopt, const void *defaults)
{
    for (size_t i = 0; i < nopt; i++) {
        char val[64];
        gd_error err = GD_ERROR_INIT;
        const char *shown = opts[i].default_note;
        if (shown == NULL) {
            if (gd_option_format(&opts[i], defaults, val, sizeof val, &err) != GD_E_OK) {
                val[0] = '?'; val[1] = '\0';
            }
            shown = val;
        }
        fprintf(stderr, "  %-20s (%s)  %s\n", opts[i].flag, shown, opts[i].help);
    }
}

/* ------------------------------------------------------------- date-tree */

static void usage_tree(void)
{
    gd_analysis_tree_config def;
    gd_analysis_tree_config_default(&def);
    size_t nopt;
    const gd_option *opts = gd_options_tree(&nopt);
    fprintf(stderr,
      "Usage: clepsydra date-tree --tree T.nwk --forms forms.tsv \\\n"
      "         --calibrations cal.tsv --out DIR [options]\n");
    print_options(opts, nopt, &def);
    fprintf(stderr, "  %-20s (0)  report sampling progress on stderr\n",
            "--progress 0|1");
    fprintf(stderr, "  %-20s (0)  exit 0 even if the convergence check fails\n",
            "--allow-unconverged 0|1");
    fprintf(stderr,
      "  Convergence check: every sampled parameter needs rank R-hat < %.2f and\n"
      "  bulk and tail ESS >= %.0f per chain. A failing run still writes its\n"
      "  outputs, then exits with status 3.\n",
      GD_CONV_RHAT_MAX, GD_CONV_ESS_PER_CHAIN);
}

int cmd_date_tree(int argc, char **argv)
{
    gd_error err = GD_ERROR_INIT;
    const char *tree_path = opt_arg(argc, argv, "--tree", NULL);
    const char *forms_path = opt_arg(argc, argv, "--forms", NULL);
    const char *cal_path = opt_arg(argc, argv, "--calibrations", NULL);
    const char *out_dir = opt_arg(argc, argv, "--out", NULL);
    if (!tree_path || !forms_path || !cal_path || !out_dir) { usage_tree(); return 1; }

    int progress = 0;
    if (parse_int_arg(opt_arg(argc, argv, "--progress", "0"), "--progress", &progress) != 0)
        return 1;
    int allow_unconverged = 0;
    if (parse_int_arg(opt_arg(argc, argv, "--allow-unconverged", "0"),
                      "--allow-unconverged", &allow_unconverged) != 0)
        return 1;

    size_t nopt;
    const gd_option *opts = gd_options_tree(&nopt);
    static const char *const tree_extra[] = {
        "--tree", "--forms", "--calibrations", "--out", "--progress",
        "--allow-unconverged", "--analysis"
    };
    if (reject_unknown_opts(argc, argv, opts, nopt,
                            tree_extra, sizeof tree_extra / sizeof tree_extra[0]) != 0) {
        return 1;
    }

    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    if (apply_registry_opts(argc, argv, opts, nopt, &cfg) != 0) return 1;

    struct cli_progress prog = { cfg.n_chains };
    if (progress) {
        cfg.on_progress = cli_progress_cb;
        cfg.progress_user = &prog;
    }

    tree_inputs in;
    if (load_tree_inputs(&in, tree_path, forms_path, cal_path) != 0) return 1;
    gd_analysis_tree_input ai;
    fill_analysis_input(&in, &ai);
    ai.analysis_id = opt_arg(argc, argv, "--analysis", NULL);

    gd_analysis_tree_job *job = NULL;
    if (gd_analysis_tree_prepare(&ai, &cfg, &job, &err) != 0) {
        fprintf(stderr, "error: %s\n", errmsg(&err));
        free_tree_inputs(&in);
        return 1;
    }

    fprintf(stderr, "clepsydra date-tree: %zu tips/%zu nodes, %zu concepts, "
            "%zu chains x (%zu warmup + %zu samples)%s\n",
            gd_analysis_tree_job_languages(job), gd_analysis_tree_job_nodes(job),
            gd_analysis_tree_job_concepts(job), cfg.n_chains,
            cfg.n_warmup, cfg.n_samples,
            cfg.relaxed ? ", relaxed clock" : "");
    if (gd_analysis_tree_job_skipped_rows(job))
        fprintf(stderr, "  note: %zu form rows for languages not in the tree were skipped\n",
                gd_analysis_tree_job_skipped_rows(job));

    gd_analysis_tree_result *res = NULL;
    int rc = gd_analysis_tree_job_run(job, &res, &err);
    gd_analysis_tree_job_destroy(job);
    if (rc != 0) {
        fprintf(stderr, "error: mcmc: %s\n", errmsg(&err));
        free_tree_inputs(&in);
        return 1;
    }

    /* ---- outputs ---- */
    if (mkdir_p(out_dir) != 0) {
        fprintf(stderr, "error: cannot create output directory '%s'\n", out_dir);
        gd_analysis_tree_result_destroy(res);
        free_tree_inputs(&in);
        return 1;
    }
    char path[1200];
    FILE *fp = NULL;

    {
        gd_table *tbl = NULL;
        if (gd_report_nodes(res, &tbl, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            goto done;
        }
        int wrc = write_table(tbl, out_dir, "posterior_nodes.tsv");
        gd_table_destroy(tbl);
        if (wrc != 0) goto done;
    }
    {
        gd_table *tbl = NULL;
        if (gd_report_diagnostics(res, &tbl, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            goto done;
        }
        int wrc = write_table(tbl, out_dir, "diagnostics.tsv");
        gd_table_destroy(tbl);
        if (wrc != 0) goto done;
    }
    if (res->n_checks > 0) {
        gd_table *tbl = NULL;
        if (gd_report_adequacy(res, &tbl, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            goto done;
        }
        int wrc = write_table(tbl, out_dir, "adequacy.tsv");
        gd_table_destroy(tbl);
        if (wrc != 0) goto done;
    }
    if (res->n_cal_checks > 0) {
        gd_table *tbl = NULL;
        if (gd_report_calibrations(res, &tbl, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            goto done;
        }
        int wrc = write_table(tbl, out_dir, "posterior_calibrations.tsv");
        gd_table_destroy(tbl);
        if (wrc != 0) goto done;
    }
    if (res->concepts != NULL) {
        gd_table *tbl = NULL;
        if (gd_report_concepts(res, &tbl, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            goto done;
        }
        int wrc = write_table(tbl, out_dir, "posterior_concepts.tsv");
        gd_table_destroy(tbl);
        if (wrc != 0) goto done;
    }

    if (!out_path(path, sizeof path, out_dir, "run_manifest.tsv")) goto done;
    fp = fopen(path, "w");
    if (fp == NULL) { fprintf(stderr, "error: cannot open %s for writing\n", path); goto done; }
    fprintf(fp, "key\tvalue\n");
    fprintf(fp, "program\tclepsydra %s\n", gd_version());
    fprintf(fp, "git_commit\t%s\n", gd_git_commit());
    fprintf(fp, "compiler\t%s\n", build_compiler());
    fprintf(fp, "platform\t%s\n", build_platform());
    fprintf(fp, "build_mode\tnative\n");
    fprintf(fp, "model\t%s\n", gd_report_tree_model_name(res));
    { char ck[CLI_CHECKSUM_SIZE];
      buffer_checksum(in.tree, in.tree_len, ck);   fprintf(fp, "input_checksum_tree\t%s\n", ck);
      buffer_checksum(in.forms, in.forms_len, ck); fprintf(fp, "input_checksum_forms\t%s\n", ck);
      buffer_checksum(in.cal, in.cal_len, ck);     fprintf(fp, "input_checksum_calibrations\t%s\n", ck);
    }
    fprintf(fp, "n_languages\t%zu\n", res->n_languages);
    fprintf(fp, "n_nodes\t%zu\n", res->n_nodes);
    fprintf(fp, "n_concepts\t%zu\n", res->n_concepts);
    fprintf(fp, "chains\t%zu\n", res->n_chains);
    fprintf(fp, "warmup\t%zu\n", res->n_warmup);
    fprintf(fp, "samples\t%zu\n", res->n_samples);
    fprintf(fp, "seed\t%lu\n", (unsigned long)cfg.seed);
    fprintf(fp, "rate_logmean_per_year\t%.4f\n", cfg.rate_logmean);
    fprintf(fp, "rate_logsd\t%.4f\n", cfg.rate_logsd);
    if (cfg.sample_epsilon) {
        for (size_t p = 0; p < res->nparams; p++)
            if (res->params[p].type == GD_MCMC_P_LOGIT_EPS)
                fprintf(fp, "contamination_epsilon\tsampled, Beta(1, 19); posterior centre %.4f\n",
                        1.0 / (1.0 + exp(-res->params[p].stat.mean)));
    } else {
        fprintf(fp, "contamination_epsilon\t%.4f (fixed)\n", cfg.epsilon);
    }
    fprintf(fp, "branch_rate_confounded\t%d\n", res->branch_rate_confounded);
    if (!isnan(res->min_branch_duration_corr))
        fprintf(fp, "min_branch_duration_corr\t%.3f\n", res->min_branch_duration_corr);
    /* Only when it happened: a clean file should not carry a line saying so. */
    if (res->n_unreadable_weights > 0)
        fprintf(fp, "unreadable_primary_weights\t%zu\n", res->n_unreadable_weights);
    fprintf(fp, "missing_cognate_rows\t%zu\n", res->n_missing_cognates);
    fprintf(fp, "loan_rows_excluded\t%zu (threshold %.2f)\n", res->n_loans, cfg.loan_threshold);
    if (res->concepts != NULL)
        fprintf(fp, "incompatible_concepts\t%zu of %zu (excluded: no T1 history on this topology fits)\n",
                res->n_incompatible, res->n_concepts);
    fprintf(fp, "scale_source\t%s\n", res->scale_source);
    fprintf(fp, "rate_contraction\t%.3f\n", res->rate_contraction);
    fprintf(fp, "scale_from_prior\t%d\n", res->scale_from_prior ? 1 : 0);
    fprintf(fp, "lexical_root_median\t%.4f\n", res->lexical_root.median);
    fprintf(fp, "lexical_root_hpd90\t%.4f-%.4f\n", res->lexical_root.hpd_lo, res->lexical_root.hpd_hi);
    fprintf(fp, "root_lex_prior\tlognormal(%.3f, %.3f)\n", cfg.root_lex_logmean, cfg.root_lex_logsd);
    fprintf(fp, "root_var_share_data\t%.3f\n", res->var_share_data);
    fprintf(fp, "root_var_share_clock\t%.3f\n", res->var_share_clock);
    {
        static const char *const names[4] = { "clock", "lexroot", "calibrations", "likelihood" };
        const double *rows[2] = { res->sens_root, res->sens_lexroot };
        const char *labels[2] = { "sensitivity_root_age", "sensitivity_lexical_root" };
        for (int q = 0; q < 2; q++) {
            fprintf(fp, "%s\t", labels[q]);
            for (int k = 0; k < 4; k++)
                fprintf(fp, "%s%s=%.3f", k ? ";" : "", names[k], rows[q][k]);
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "converged\t%d\n", res->converged ? 1 : 0);
    fprintf(fp, "n_unconverged_params\t%zu\n", res->n_unconverged);
    fprintf(fp, "max_rhat\t%.4f\n", res->max_rhat);
    fprintf(fp, "min_ess_bulk\t%.0f\n", res->min_ess_bulk);
    fprintf(fp, "min_ess_tail\t%.0f\n", res->min_ess_tail);
    fprintf(fp, "note\tabsolute ages depend on the rate prior and calibrations; see docs/scientific-methods.md\n");
    if (fclose(fp) != 0) { fp = NULL; fprintf(stderr, "error: failed to flush %s\n", path); goto done; }
    fp = NULL;

    fprintf(stderr, "wrote posterior_nodes.tsv, diagnostics.tsv, %srun_manifest.tsv to %s/\n",
            res->concepts != NULL ? "posterior_concepts.tsv, " : "", out_dir);

    int status = 0;
    if (!res->converged) {
        fprintf(stderr,
                "%s: NOT CONVERGED: %zu parameter(s) fail (max R-hat %.3f, "
                "min bulk ESS %.0f, min tail ESS %.0f; need R-hat < %.2f and "
                "ESS >= %.0f). See diagnostics.tsv; run longer or more chains.\n",
                allow_unconverged ? "warning" : "error",
                res->n_unconverged, res->max_rhat, res->min_ess_bulk,
                res->min_ess_tail, GD_CONV_RHAT_MAX,
                GD_CONV_ESS_PER_CHAIN * (double)res->n_chains);
        if (!allow_unconverged) status = 3;
    }

    gd_analysis_tree_result_destroy(res);
    free_tree_inputs(&in);
    return status;

done:
    if (fp != NULL) fclose(fp);
    gd_analysis_tree_result_destroy(res);
    free_tree_inputs(&in);
    return 1;
}

/* ------------------------------------------------------- calib-sensitivity */

static void usage_calib(void)
{
    fprintf(stderr,
      "Usage: clepsydra calib-sensitivity --tree T.nwk --forms forms.tsv \\\n"
      "         --calibrations cal.tsv --out DIR [options]\n"
      "  --mode loco|widen (loco)  leave-one-calibration-out, or loosen every\n"
      "                            calibration window\n"
      "  --factor X (2.0)          widen factor for --mode widen\n"
      "  --include-extant 0|1 (0)  in loco, also drop age-0 (modern) tip calibrations\n");
    {
        gd_analysis_tree_config def;
        gd_analysis_tree_config_default(&def);
        size_t nopt;
        const gd_option *opts = gd_options_tree(&nopt);
        print_options(opts, nopt, &def);
    }
    fprintf(stderr, "Writes calibration_sensitivity.tsv to DIR.\n");
}

/* One run of the chronology under a stated variation of the calibrations. */
static gd_analysis_tree_result *run_variant(const gd_analysis_tree_input *ai,
                                            const gd_analysis_tree_config *base,
                                            size_t drop_index, double widen,
                                            const char *what)
{
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_config cfg = *base;
    cfg.calibrations.drop_index = drop_index;
    cfg.calibrations.widen = widen;
    gd_analysis_tree_result *res = NULL;
    if (gd_analysis_tree_run(ai, &cfg, &res, &err) != 0) {
        fprintf(stderr, "error: %s: %s\n", what, errmsg(&err));
        return NULL;
    }
    return res;
}

int cmd_calib_sensitivity(int argc, char **argv)
{
    gd_error err = GD_ERROR_INIT;
    const char *tree_path = opt_arg(argc, argv, "--tree", NULL);
    const char *forms_path = opt_arg(argc, argv, "--forms", NULL);
    const char *cal_path = opt_arg(argc, argv, "--calibrations", NULL);
    const char *out_dir = opt_arg(argc, argv, "--out", NULL);
    const char *mode = opt_arg(argc, argv, "--mode", "loco");
    if (!tree_path || !forms_path || !cal_path || !out_dir) { usage_calib(); return 1; }
    int do_loco = (strcmp(mode, "loco") == 0);
    if (!do_loco && strcmp(mode, "widen") != 0) {
        fprintf(stderr, "error: --mode must be loco or widen\n");
        return 1;
    }

    double factor = 2.0;
    int include_extant = 0;
    if (parse_finite(opt_arg(argc, argv, "--factor", "2.0"), "--factor", &factor) != 0 ||
        parse_int_arg(opt_arg(argc, argv, "--include-extant", "0"), "--include-extant",
                      &include_extant) != 0)
        return 1;
    if (factor <= 1.0) {
        fprintf(stderr, "error: --factor must be greater than 1\n");
        return 1;
    }

    size_t nopt;
    const gd_option *opts = gd_options_tree(&nopt);
    static const char *const calib_extra[] = {
        "--tree", "--forms", "--calibrations", "--out",
        "--mode", "--factor", "--include-extant", "--analysis"
    };
    if (reject_unknown_opts(argc, argv, opts, nopt,
                            calib_extra, sizeof calib_extra / sizeof calib_extra[0]) != 0) {
        return 1;
    }

    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    if (apply_registry_opts(argc, argv, opts, nopt, &cfg) != 0) return 1;

    tree_inputs in;
    if (load_tree_inputs(&in, tree_path, forms_path, cal_path) != 0) return 1;
    gd_analysis_tree_input ai;
    fill_analysis_input(&in, &ai);
    ai.analysis_id = opt_arg(argc, argv, "--analysis", NULL);

    gd_analysis_calibration *cals = NULL;
    size_t ncal = 0;
    gd_analysis_tree_result *base = NULL, *var = NULL;
    FILE *fp = NULL;
    char path[1200];
    int rc = 1;

    if (gd_analysis_calibrations_read(in.cal, in.cal_len,
                                      &cals, &ncal, &err) != 0) {
        fprintf(stderr, "error: calibrations: %s\n", errmsg(&err));
        goto done;
    }

    fprintf(stderr, "calib-sensitivity (%s): baseline over %zu calibrations, "
            "%zu chains x (%zu+%zu)\n", mode, ncal, cfg.n_chains,
            cfg.n_warmup, cfg.n_samples);

    base = run_variant(&ai, &cfg, GD_ANALYSIS_KEEP_ALL, 1.0, "baseline");
    if (base == NULL) goto done;

    /* Divergence ceiling: dropping the only calibration that bounds a node
     * leaves its age unidentifiable, and the sampler wanders to absurd values.
     * Anything beyond 100x the calibrated/baseline scale is reported as
     * "unbounded" rather than as a spurious number. */
    {
        double ceiling = gd_analysis_divergence_ceiling(base, cals, ncal);

        if (mkdir_p(out_dir) != 0) {
            fprintf(stderr, "error: cannot create output directory '%s'\n", out_dir);
            goto done;
        }
        if (!out_path(path, sizeof path, out_dir, "calibration_sensitivity.tsv"))
            goto done;
        fp = fopen(path, "w");
        if (fp == NULL) { fprintf(stderr, "error: cannot open %s\n", path); goto done; }

        if (do_loco) {
            fprintf(fp, "dropped_id\ttarget\ttype\theld_out_bp\tloo_median_bp\t"
                        "loo_hpd_lo_bp\tloo_hpd_hi_bp\trecovered\tabs_err_bp\t"
                        "max_other_shift_bp\n");
            size_t n_cov = 0, n_scored = 0, n_unbounded = 0;
            for (size_t r = 0; r < ncal; r++) {
                const gd_analysis_calibration *c = &cals[r];
                /* A living language's age is not a prediction worth scoring. */
                if (!include_extant && c->is_extant_tip) continue;
                if (c->type == GD_CAL_TIP_ERROR) continue;   /* not an age to recover */

                var = run_variant(&ai, &cfg, r, 1.0, c->id);
                if (var == NULL) goto done;

                gd_sensitivity_row row;
                if (gd_analysis_sensitivity_score(base, var, c, ceiling,
                                                  &row, &err) != GD_E_OK) {
                    fprintf(stderr, "error: %s\n", errmsg(&err));
                    goto done;
                }

                switch (row.outcome) {
                case GD_SENS_SCORED:
                    n_cov += row.covered ? 1u : 0u; n_scored++;
                    fprintf(fp, "%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%d\t%.1f\t%.1f\n",
                            c->id, c->target, c->type_name, c->point,
                            row.median, row.hpd_lo, row.hpd_hi,
                            row.covered ? 1 : 0, row.abs_err, row.max_other_shift);
                    break;
                case GD_SENS_UNBOUNDED:
                    n_unbounded++;
                    fprintf(fp, "%s\t%s\t%s\t%.1f\tunbounded\tunbounded\tunbounded\t"
                                "unbd\t.\t%.1f\n",
                            c->id, c->target, c->type_name, c->point,
                            row.max_other_shift);
                    break;
                case GD_SENS_ABSENT:
                    fprintf(fp, "%s\t%s\t%s\t%.1f\t.\t.\t.\t.\t.\t%.1f\n",
                            c->id, c->target, c->type_name, c->point,
                            row.max_other_shift);
                    break;
                }
                gd_analysis_tree_result_destroy(var);
                var = NULL;
            }
            if (fclose(fp) != 0) { fp = NULL; fprintf(stderr, "error: failed to flush %s\n", path); goto done; }
            fp = NULL;
            fprintf(stderr, "leave-one-calibration-out: %zu/%zu held-out calibrations "
                    "recovered in their 90%% HPD", n_cov, n_scored);
            if (n_unbounded)
                fprintf(stderr, "; %zu load-bearing (target unbounded when dropped)", n_unbounded);
            fprintf(stderr, " -> %s\n", path);
        } else {
            var = run_variant(&ai, &cfg, GD_ANALYSIS_KEEP_ALL, factor, "widen");
            if (var == NULL) goto done;
            fprintf(fp, "node\tbaseline_median_bp\twidened_median_bp\tshift_bp\t"
                        "baseline_width_bp\twidened_width_bp\n");
            double worst = 0.0;
            for (size_t v = 0; v < base->n_nodes; v++) {
                const gd_analysis_node_age *b = gd_analysis_result_node(base, v);
                const gd_analysis_node_age *w = gd_analysis_result_node(var, v);
                if (b == NULL && w == NULL) continue;
                const char *lbl = base->tree->nodes[v].label;
                double sh = (b != NULL && w != NULL)
                          ? fabs(w->age.median - b->age.median) : NAN;
                if (!isnan(sh) && sh > worst) worst = sh;
                fprintf(fp, "%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
                        lbl ? lbl : "(node)",
                        b ? b->age.median : NAN, w ? w->age.median : NAN, sh,
                        b ? (b->age.hpd_hi - b->age.hpd_lo) : NAN,
                        w ? (w->age.hpd_hi - w->age.hpd_lo) : NAN);
            }
            if (fclose(fp) != 0) { fp = NULL; fprintf(stderr, "error: failed to flush %s\n", path); goto done; }
            fp = NULL;
            fprintf(stderr, "widen x%.2f: max node-age shift %.0f BP -> %s\n",
                    factor, worst, path);
        }
    }
    rc = 0;

done:
    if (fp != NULL) fclose(fp);
    gd_analysis_tree_result_destroy(var);
    gd_analysis_tree_result_destroy(base);
    gd_analysis_calibrations_free(cals, ncal);
    free_tree_inputs(&in);
    return rc;
}

/* ------------------------------------------------------------- date-pair */

static void usage_pair(void)
{
    gd_analysis_pair_config def;
    gd_analysis_pair_config_default(&def);
    size_t nopt;
    const gd_option *opts = gd_options_pair(&nopt);
    fprintf(stderr,
      "Usage: clepsydra date-pair --forms forms.tsv --lang-a A --lang-b B \\\n"
      "         --out DIR [options]\n");
    print_options(opts, nopt, &def);
    fprintf(stderr,
      "  %-20s      swadesh|starostin|vasilyev|two-component\n"
      "  Default is the hierarchical constant-rate model (P1). --decay selects an\n"
      "  alternative retention law (Starostin accelerating, Vasilyev two-component)\n"
      "  for model comparison. The rate prior is per YEAR, as in date-tree.\n"
      "  T is reported in years BP; with --age-a/--age-b the lineage separation\n"
      "  2T - age_a - age_b is reported too (for an ancestor-descendant pair it\n"
      "  is the time elapsed between the two lists).\n",
      "--decay");
}

int cmd_date_pair(int argc, char **argv)
{
    gd_error err = GD_ERROR_INIT;
    const char *forms_path = opt_arg(argc, argv, "--forms", NULL);
    const char *la = opt_arg(argc, argv, "--lang-a", NULL);
    const char *lb = opt_arg(argc, argv, "--lang-b", NULL);
    const char *out_dir = opt_arg(argc, argv, "--out", NULL);
    if (!forms_path || !la || !lb || !out_dir) { usage_pair(); return 1; }

    size_t nopt;
    const gd_option *opts = gd_options_pair(&nopt);
    static const char *const pair_extra[] = {
        "--forms", "--lang-a", "--lang-b", "--out", "--decay", "--analysis"
    };
    if (reject_unknown_opts(argc, argv, opts, nopt,
                            pair_extra, sizeof pair_extra / sizeof pair_extra[0]) != 0) {
        return 1;
    }

    gd_analysis_pair_config cfg;
    gd_analysis_pair_config_default(&cfg);
    if (apply_registry_opts(argc, argv, opts, nopt, &cfg) != 0) return 1;

    /* The registry bounds each option on its own; that the grid's two ends are
     * ordered is a relation between two of them, which no single row can say. */
    if (!(cfg.t_hi > cfg.t_lo)) {
        fprintf(stderr, "error: need --t-lo < --t-hi\n");
        return 1;
    }
    /* Parsed through the library so the CLI and the browser accept exactly the
     * same spellings; an unknown one is refused here rather than silently
     * dating the pair under a different model. */
    {
        gd_error perr = GD_ERROR_INIT;
        if (gd_pair_model_from_string(opt_arg(argc, argv, "--decay", NULL),
                                      &cfg.model, &perr) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&perr));
            return 2;
        }
    }

    size_t forms_len = 0;
    char *forms_buf = read_file(forms_path, &forms_len);
    if (forms_buf == NULL) {
        fprintf(stderr, "error: cannot read %s\n", forms_path);
        return 1;
    }

    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = forms_buf;
    in.forms_len = forms_len;
    in.forms_name = forms_path;
    in.lang_a = la;
    in.lang_b = lb;
    in.analysis_id = opt_arg(argc, argv, "--analysis", NULL);

    gd_analysis_pair_result *res = NULL;
    int rc = gd_analysis_pair_run(&in, &cfg, &res, &err);
    free(forms_buf);
    if (rc != 0) {
        fprintf(stderr, "error: %s\n", errmsg(&err));
        return 1;
    }

    if (mkdir_p(out_dir) != 0) {
        fprintf(stderr, "error: cannot create output directory '%s'\n", out_dir);
        gd_analysis_pair_result_destroy(res);
        return 1;
    }
    {
        gd_table *tbl = NULL;
        if (gd_report_pair(res, la, lb, &tbl, &err) != GD_E_OK) {
            fprintf(stderr, "error: %s\n", errmsg(&err));
            gd_analysis_pair_result_destroy(res);
            return 1;
        }
        int wrc = write_table(tbl, out_dir, "posterior_pair.tsv");
        gd_table_destroy(tbl);
        if (wrc != 0) {
            gd_analysis_pair_result_destroy(res);
            return 1;
        }
    }

    fprintf(stderr, "date-pair %s vs %s: %zu shared, %zu cognate; "
            "divergence median %.0f BP; wrote %s/posterior_pair.tsv\n",
            la, lb, res->n_shared, res->n_cognate, res->median, out_dir);

    gd_analysis_pair_result_destroy(res);
    return 0;
}

/* ------------------------------------------------------------- date-trees */

/*
 * Dating over a posterior sample of trees (roadmap E3.4). Each selected tree
 * is dated with the full model; each named clade's MRCA age is pooled over
 * trees with equal weight per tree. When the trees were estimated from the
 * same cognates, this reuses the data (a hybrid, not one joint posterior):
 * the README of the output says so.
 */
static void usage_trees(void)
{
    fprintf(stderr,
      "Usage: clepsydra date-trees --trees SAMPLE.trees --forms forms.tsv \\\n"
      "         --calibrations cal.tsv --out DIR [--clades clades.tsv] \\\n"
      "         [--burnin 0.25] [--ntrees 50] [--draws-per-chain 100] [options]\n"
      "  --trees: Newick (one per ';') or NEXUS with a Translate table (BEAST,\n"
      "  MrBayes, RevBayes). Calibrations must name taxon sets (mrca:A,B) or tips:\n"
      "  sampled trees have no internal labels. --clades: TSV with columns\n"
      "  clade_id and members (comma-separated tips); the root is always reported.\n"
      "  Options of date-tree apply to every tree.\n");
}

static int cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Median and narrowest 90% interval of n values (sorted in place). */
static void summ90(double *x, size_t n, double *med, double *lo, double *hi)
{
    if (n == 0) { *med = *lo = *hi = NAN; return; }
    qsort(x, n, sizeof(double), cmp_dbl);
    *med = (n % 2) ? x[n / 2] : 0.5 * (x[n / 2 - 1] + x[n / 2]);
    size_t w = (size_t)ceil(0.9 * (double)n);
    if (w < 1) w = 1;
    size_t bi = 0; double best = HUGE_VAL;
    for (size_t i = 0; i + w - 1 < n; i++)
        if (x[i + w - 1] - x[i] < best) { best = x[i + w - 1] - x[i]; bi = i; }
    *lo = x[bi]; *hi = x[bi + w - 1];
}

static char *dup_cstr(const char *s)
{
    size_t n = strlen(s) + 1;
    char *d = (char *)malloc(n);
    if (d) memcpy(d, s, n);
    return d;
}

typedef struct {
    char   *id, *members;
    double *all;  size_t nall, capall;     /* MRCA age, every tree */
    double *mono; size_t nmono, capmono;   /* only trees where monophyletic */
    size_t  ntrees_mono;
} clade_acc;

static int acc_push(double **a, size_t *n, size_t *cap, const double *x, size_t k)
{
    if (*n + k > *cap) {
        size_t nc = (*n + k) * 2;
        double *g = (double *)realloc(*a, nc * sizeof(double));
        if (g == NULL) return -1;
        *a = g; *cap = nc;
    }
    memcpy(*a + *n, x, k * sizeof(double));
    *n += k;
    return 0;
}

int cmd_date_trees(int argc, char **argv)
{
    const char *trees_path = opt_arg(argc, argv, "--trees", NULL);
    const char *forms_path = opt_arg(argc, argv, "--forms", NULL);
    const char *cal_path = opt_arg(argc, argv, "--calibrations", NULL);
    const char *out_dir = opt_arg(argc, argv, "--out", NULL);
    const char *clades_path = opt_arg(argc, argv, "--clades", NULL);
    if (!trees_path || !forms_path || !cal_path || !out_dir) { usage_trees(); return 1; }
    size_t nopt;
    const gd_option *opts = gd_options_tree(&nopt);
    static const char *const extra[] = {
        "--trees", "--forms", "--calibrations", "--out", "--clades", "--burnin",
        "--ntrees", "--draws-per-chain", "--allow-unconverged", "--analysis"
    };
    if (reject_unknown_opts(argc, argv, opts, nopt, extra, sizeof extra / sizeof extra[0]) != 0) return 1;
    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    cfg.compare_prior = false;     /* per-tree prior runs would double the cost */
    if (apply_registry_opts(argc, argv, opts, nopt, &cfg) != 0) return 1;
    double burnin = atof(opt_arg(argc, argv, "--burnin", "0.25"));
    int ntrees_req = 50, dpc = 100, allow = 0;
    if (parse_int_arg(opt_arg(argc, argv, "--ntrees", "50"), "--ntrees", &ntrees_req) != 0 ||
        parse_int_arg(opt_arg(argc, argv, "--draws-per-chain", "100"), "--draws-per-chain", &dpc) != 0 ||
        parse_int_arg(opt_arg(argc, argv, "--allow-unconverged", "0"), "--allow-unconverged", &allow) != 0)
        return 1;
    if (!(burnin >= 0.0 && burnin < 1.0) || ntrees_req < 1 || dpc < 1) {
        fprintf(stderr, "error: need 0 <= --burnin < 1, --ntrees >= 1, --draws-per-chain >= 1\n");
        return 1;
    }
    cfg.trace_draws = (size_t)dpc;

    size_t tlen = 0, flen = 0, clen = 0;
    char *tbuf = read_file(trees_path, &tlen);
    char *fbuf = read_file(forms_path, &flen);
    char *cbuf = read_file(cal_path, &clen);
    char *kbuf = clades_path ? read_file(clades_path, NULL) : NULL;
    gd_treeset *ts = NULL;
    gd_error err = GD_ERROR_INIT;
    clade_acc *cl = NULL; size_t ncl = 0;
    FILE *fp = NULL, *ft = NULL;
    int status = 1;
    if (!tbuf || !fbuf || !cbuf || (clades_path && !kbuf)) { fprintf(stderr, "error: cannot read an input file\n"); goto out; }
    if (gd_treeset_read(tbuf, tlen, &ts, &err) != 0) { fprintf(stderr, "error: %s\n", errmsg(&err)); goto out; }

    /* Clades: the root (all tips of the first tree), then the file's rows. */
    {
        char *all = NULL;
        if (gd_newick_tip_list(ts->newick[0], &all, &err) != 0) {
            fprintf(stderr, "error: first tree: %s\n", errmsg(&err)); goto out;
        }
        size_t ncap = 8;
        cl = (clade_acc *)calloc(ncap, sizeof *cl);
        if (cl == NULL) { free(all); fprintf(stderr, "error: out of memory\n"); goto out; }
        cl[0].id = dup_cstr("root"); cl[0].members = all; ncl = 1;
        if (kbuf) {
            int hdr = 1;
            char *line = kbuf;
            while (line != NULL && *line != '\0') {
                char *nl = strchr(line, '\n');
                if (nl) *nl = '\0';
                char *tab = strchr(line, '\t');
                if (!hdr && tab) {
                    *tab = '\0';
                    char *mem = tab + 1;
                    size_t ml = strlen(mem);
                    while (ml && (mem[ml - 1] == '\r' || mem[ml - 1] == ' ')) mem[--ml] = '\0';
                    if (ncl == ncap) {
                        ncap *= 2;
                        clade_acc *g = (clade_acc *)realloc(cl, ncap * sizeof *cl);
                        if (!g) { fprintf(stderr, "error: out of memory\n"); goto out; }
                        cl = g;
                    }
                    memset(&cl[ncl], 0, sizeof cl[ncl]);
                    cl[ncl].id = dup_cstr(line); cl[ncl].members = dup_cstr(mem);
                    ncl++;
                }
                hdr = 0;
                line = nl ? nl + 1 : NULL;
            }
        }
    }

    size_t first = (size_t)floor(burnin * (double)ts->ntrees);
    size_t avail = ts->ntrees - first;
    size_t nuse = ((size_t)ntrees_req < avail) ? (size_t)ntrees_req : avail;
    if (mkdir_p(out_dir) != 0) { fprintf(stderr, "error: cannot create '%s'\n", out_dir); goto out; }
    char path[1200];
    if (!out_path(path, sizeof path, out_dir, "trees.tsv") || (ft = fopen(path, "w")) == NULL) {
        fprintf(stderr, "error: cannot write trees.tsv\n"); goto out;
    }
    fprintf(ft, "tree_index\tconverged\tmax_rhat\tmin_ess_bulk\troot_median_bp\n");
    size_t nconv = 0;
    double *buf = (double *)malloc((size_t)dpc * cfg.n_chains * sizeof(double) + sizeof(double));
    if (!buf) { fprintf(stderr, "error: out of memory\n"); goto out; }
    fprintf(stderr, "clepsydra date-trees: %zu trees in the sample, %zu after burn-in, dating %zu\n",
            ts->ntrees, avail, nuse);
    for (size_t j = 0; j < nuse; j++) {
        size_t ti = first + ((nuse == 1) ? 0 : j * (avail - 1) / (nuse - 1));
        gd_analysis_tree_input ai;
        memset(&ai, 0, sizeof ai);
        ai.newick = ts->newick[ti]; ai.newick_len = GD_LEN_NUL_TERMINATED;
        ai.forms_tsv = fbuf; ai.forms_len = flen;
        ai.calibrations_tsv = cbuf; ai.cal_len = clen;
        ai.newick_name = trees_path; ai.forms_name = forms_path; ai.cal_name = cal_path;
        ai.analysis_id = opt_arg(argc, argv, "--analysis", NULL);
        cfg.seed += 1;
        gd_analysis_tree_result *res = NULL;
        gd_error_clear(&err);
        if (gd_analysis_tree_run(&ai, &cfg, &res, &err) != 0) {
            fprintf(stderr, "error: tree %zu: %s\n", ti, errmsg(&err)); free(buf); goto out;
        }
        if (res->converged) nconv++;
        double rmed = NAN;
        for (size_t c = 0; c < ncl; c++) {
            size_t nd = 0; bool mono = false;
            gd_error_clear(&err);
            if (gd_analysis_clade_draws(res, cl[c].members, buf, (size_t)dpc * cfg.n_chains, &nd, &mono, &err) != 0) {
                fprintf(stderr, "error: clade %s: %s\n", cl[c].id, errmsg(&err));
                gd_analysis_tree_result_destroy(res); free(buf); goto out;
            }
            if (c == 0 && nd > 0) { double *tmp = (double *)malloc(nd * sizeof(double)); if (tmp) { double lo, hi; memcpy(tmp, buf, nd * sizeof(double)); summ90(tmp, nd, &rmed, &lo, &hi); free(tmp); } }
            if (acc_push(&cl[c].all, &cl[c].nall, &cl[c].capall, buf, nd) != 0 ||
                (mono && acc_push(&cl[c].mono, &cl[c].nmono, &cl[c].capmono, buf, nd) != 0)) {
                fprintf(stderr, "error: out of memory\n"); gd_analysis_tree_result_destroy(res); free(buf); goto out;
            }
            if (mono) cl[c].ntrees_mono++;
        }
        fprintf(ft, "%zu\t%d\t%.4f\t%.0f\t%.1f\n", ti, res->converged ? 1 : 0, res->max_rhat,
                res->min_ess_bulk, rmed);
        fprintf(stderr, "  tree %zu/%zu (sample %zu): %s, root median %.0f BP\n", j + 1, nuse, ti,
                res->converged ? "converged" : "NOT converged", rmed);
        gd_analysis_tree_result_destroy(res);
    }
    free(buf);
    fclose(ft); ft = NULL;

    if (!out_path(path, sizeof path, out_dir, "posterior_clades.tsv") || (fp = fopen(path, "w")) == NULL) {
        fprintf(stderr, "error: cannot write posterior_clades.tsv\n"); goto out;
    }
    fprintf(fp, "clade\tmembers\tp_monophyletic\tmrca_median_bp\tmrca_hpd90_lo_bp\tmrca_hpd90_hi_bp\t"
                "if_mono_median_bp\tif_mono_hpd90_lo_bp\tif_mono_hpd90_hi_bp\n");
    for (size_t c = 0; c < ncl; c++) {
        double m, lo, hi, m2, lo2, hi2;
        summ90(cl[c].all, cl[c].nall, &m, &lo, &hi);
        summ90(cl[c].mono, cl[c].nmono, &m2, &lo2, &hi2);
        fprintf(fp, "%s\t%s\t%.3f\t%.1f\t%.1f\t%.1f\t", cl[c].id, cl[c].members,
                (double)cl[c].ntrees_mono / (double)nuse, m, lo, hi);
        if (cl[c].nmono > 0) fprintf(fp, "%.1f\t%.1f\t%.1f\n", m2, lo2, hi2);
        else fprintf(fp, ".\t.\t.\n");
    }
    fclose(fp); fp = NULL;
    fprintf(stderr, "wrote posterior_clades.tsv, trees.tsv to %s/ (%zu of %zu trees converged)\n",
            out_dir, nconv, nuse);
    status = (nconv == nuse || allow) ? 0 : 3;
    if (nconv < nuse)
        fprintf(stderr, "%s: %zu tree(s) did not pass the convergence check\n",
                allow ? "warning" : "error", nuse - nconv);
out:
    if (fp) fclose(fp);
    if (ft) fclose(ft);
    for (size_t c = 0; c < ncl; c++) { free(cl[c].id); free(cl[c].members); free(cl[c].all); free(cl[c].mono); }
    free(cl);
    gd_treeset_destroy(ts);
    free(tbuf); free(fbuf); free(cbuf); free(kbuf);
    return status;
}

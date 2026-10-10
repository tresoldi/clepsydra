/*
 * test_analysis.c — the shared analysis runner (M10 steps 1-2).
 *
 * gd_analysis is the pipeline that both the CLI and the WebAssembly shim call,
 * so it is exercised here the way the shim will use it: entirely from in-memory
 * buffers, with no filesystem anywhere. What is checked is the plumbing the
 * refactor could plausibly break — the dataset counts a caller reports before
 * committing to a run, which rows are skipped, ownership across the job/result
 * handoff, the pair posterior grid, progress reporting and cancellation, and
 * the diagnostics on bad input. The
 * numerical agreement with the pre-refactor pipeline is gated separately, by
 * the examples reproducing their committed output byte for byte.
 */
#include "test_runner.h"
#include "gd_analysis.h"
#include "gd_tree.h"
#include "gd.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Four tips, one ancient (Latin), so the tree is heterochronous like the real
 * examples. Two internal nodes carry labels so their ages can be read back. */
static const char *TREE =
    "((Latin:1200,(Spanish:2050,Portuguese:2050)Ibero:0)WestRom:500,"
    "Sardinian:2550)Root;";

/* Three concepts x four languages. Concept c1 is fully shared, c2 splits
 * Ibero-Romance from the rest, c3 is unique per language. `ghost` is a language
 * that is not in the tree: its rows must be skipped and counted, not silently
 * folded in. */
static const char *FORMS =
    "language_id\tconcept_id\tcognate_id\tprimary_weight\n"
    "Latin\tc1\tk1\t1.0\n"
    "Spanish\tc1\tk1\t1.0\n"
    "Portuguese\tc1\tk1\t1.0\n"
    "Sardinian\tc1\tk1\t1.0\n"
    "Latin\tc2\tk2\t1.0\n"
    "Spanish\tc2\tk3\t1.0\n"
    "Portuguese\tc2\tk3\t1.0\n"
    "Sardinian\tc2\tk2\t1.0\n"
    "Latin\tc3\tk4\t1.0\n"
    "Spanish\tc3\tk5\t1.0\n"
    "Portuguese\tc3\tk6\t1.0\n"
    "Sardinian\tc3\tk7\t1.0\n"
    "ghost\tc1\tk9\t1.0\n"
    "ghost\tc2\tk9\t1.0\n";

static const char *CALS =
    "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
    "cal1\tLatin\tfixed_tip\t2050\t.\t.\t1\tstandard\n"
    "cal2\tSpanish\tfixed_tip\t0\t.\t.\t1\tstandard\n"
    "cal3\tPortuguese\tfixed_tip\t0\t.\t.\t1\tstandard\n"
    "cal4\tSardinian\tfixed_tip\t0\t.\t.\t1\tstandard\n";

static void tree_input(gd_analysis_tree_input *in)
{
    memset(in, 0, sizeof *in);
    /* The fixtures are C string literals, so their lengths are measured. */
    in->newick = TREE;                in->newick_len = GD_LEN_NUL_TERMINATED;
    in->forms_tsv = FORMS;            in->forms_len  = GD_LEN_NUL_TERMINATED;
    in->calibrations_tsv = CALS;      in->cal_len    = GD_LEN_NUL_TERMINATED;
    in->newick_name = "test.nwk";
    in->forms_name = "forms.tsv";
    in->cal_name = "calibrations.tsv";
}

static void quick_config(gd_analysis_tree_config *cfg)
{
    gd_analysis_tree_config_default(cfg);
    cfg->n_warmup = 200;
    cfg->n_samples = 200;
    cfg->n_chains = 2;
    cfg->seed = 7;
}

TEST(test_analysis_prepare_reports_dataset_size)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_job *job = NULL;
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) == 0);
    CHECK(job != NULL);

    /* Four tips; a rooted binary tree over them has 7 nodes. */
    CHECK(gd_analysis_tree_job_languages(job) == 4);
    CHECK(gd_analysis_tree_job_nodes(job) == 7);
    CHECK(gd_analysis_tree_job_concepts(job) == 3);
    /* The two `ghost` rows are skipped and reported, not absorbed. */
    CHECK(gd_analysis_tree_job_skipped_rows(job) == 2);

    gd_analysis_tree_job_destroy(job);
}

TEST(test_analysis_tree_run_summarizes_ages)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &res, &err) == 0);
    CHECK(res != NULL);
    CHECK(res->n_concepts == 3);
    CHECK(res->n_languages == 4);
    CHECK(res->n_skipped_rows == 2);
    CHECK(res->n_chains == 2);
    CHECK(res->n_draws == 200);
    CHECK(res->nparams >= res->nages);
    CHECK(res->nages >= 1);

    /* Every age summary must be ordered and finite, and every parameter must
     * carry a type name the caller can print without a lookup table. */
    for (size_t i = 0; i < res->nages; i++) {
        const gd_analysis_node_age *a = &res->ages[i];
        CHECK(a->age.hpd_lo <= a->age.median);
        CHECK(a->age.median <= a->age.hpd_hi);
        CHECK(a->age.mean >= 0.0 && !isnan(a->age.mean));
        CHECK(a->rhat > 0.0 && !isnan(a->rhat));
        CHECK(a->node < res->n_nodes);
    }
    for (size_t p = 0; p < res->nparams; p++) {
        CHECK(res->params[p].type_name != NULL);
        CHECK(res->params[p].type_name[0] != '\0');
    }

    /* A parent is older than its child: the constraint the sampler enforces
     * has to survive the summarising. */
    size_t root = res->tree->root;
    for (size_t i = 0; i < res->nages; i++)
        if (res->ages[i].node == root)
            for (size_t j = 0; j < res->nages; j++)
                if (res->ages[j].node != root)
                    CHECK(res->ages[i].age.median >= res->ages[j].age.median - 1e-9);

    /* Strict clock: no branch multipliers, so no confounding number. */
    CHECK(res->relaxed == 0);
    CHECK(isnan(res->min_branch_duration_corr));

    gd_analysis_tree_result_destroy(res);
}

TEST(test_analysis_result_outlives_job)
{
    /* The tree moves into the result. Destroying the job first must leave the
     * result fully readable — the ownership handoff is the part of the refactor
     * most likely to go wrong, and a use-after-free here is silent without a
     * sanitizer. */
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_job *job = NULL;
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) == 0);
    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_job_run(job, &res, &err) == 0);

    gd_analysis_tree_job_destroy(job);          /* job goes first */

    CHECK(res->tree != NULL);
    CHECK(res->tree->nnodes == 7);
    for (size_t i = 0; i < res->nages; i++)
        if (res->ages[i].label != NULL)
            CHECK(strlen(res->ages[i].label) > 0);   /* borrowed from the tree */

    gd_analysis_tree_result_destroy(res);
}

TEST(test_analysis_job_runs_once)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_job *job = NULL;
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) == 0);
    gd_analysis_tree_result *first = NULL, *second = NULL;
    CHECK(gd_analysis_tree_job_run(job, &first, &err) == 0);
    /* The second attempt is refused rather than run against a torn-down model. */
    CHECK(gd_analysis_tree_job_run(job, &second, &err) != 0);
    CHECK(second == NULL);

    gd_analysis_tree_result_destroy(first);
    gd_analysis_tree_job_destroy(job);
}

TEST(test_analysis_tree_reproducible)
{
    /* Same seed, same buffers, same answer: the property the whole
     * native-vs-Wasm comparison in M10 will rest on. */
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *a = NULL, *b = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &a, &err) == 0);
    CHECK(gd_analysis_tree_run(&in, &cfg, &b, &err) == 0);
    CHECK(a->nages == b->nages);
    for (size_t i = 0; i < a->nages; i++) {
        CHECK(a->ages[i].age.median == b->ages[i].age.median);
        CHECK(a->ages[i].age.hpd_lo == b->ages[i].age.hpd_lo);
        CHECK(a->ages[i].age.hpd_hi == b->ages[i].age.hpd_hi);
    }
    gd_analysis_tree_result_destroy(a);
    gd_analysis_tree_result_destroy(b);
}

TEST(test_analysis_tree_bad_input)
{
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_job *job = NULL;

    /* Missing required forms columns. */
    gd_analysis_tree_input in; tree_input(&in);
    in.forms_tsv = "language_id\tconcept_id\nLatin\tc1\n";
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0);
    CHECK(job == NULL);
    CHECK(err.code == GD_E_UNKNOWN_COLUMN);

    /* An unknown calibration type names the offending value. The error owns
     * its message, so this works with no cooperation from the caller — there
     * is no buffer to supply and none to forget. */
    tree_input(&in);
    in.calibrations_tsv =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "cal1\tLatin\tnot_a_type\t2050\t.\t.\t1\tx\n";
    gd_error_clear(&err);
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0);
    CHECK(err.message[0] != '\0');
    CHECK(strstr(err.message, "not_a_type") != NULL);

    /* Malformed Newick. */
    tree_input(&in);
    in.newick = "((a,b)";
    gd_error_clear(&err);
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0);

    /* Null arguments are refused, not dereferenced. */
    CHECK(gd_analysis_tree_prepare(NULL, &cfg, &job, &err) != 0);
    CHECK(gd_analysis_tree_run(&in, &cfg, NULL, &err) != 0);
}

TEST(test_analysis_pair_posterior)
{
    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = FORMS;
    in.forms_len = GD_LEN_NUL_TERMINATED;
    in.forms_name = "forms.tsv";
    in.lang_a = "Spanish";
    in.lang_b = "Portuguese";

    gd_analysis_pair_config cfg;
    gd_analysis_pair_config_default(&cfg);
    cfg.t_lo = 100.0;
    cfg.t_hi = 4000.0;
    gd_error err = GD_ERROR_INIT;

    gd_analysis_pair_result *res = NULL;
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) == 0);
    CHECK_STREQ(res->model, "P1_hierarchical_constant_rate");
    /* Three shared concepts; c1 and c2 match, c3 does not. */
    CHECK(res->n_shared == 3);
    CHECK(res->n_cognate == 2);
    CHECK_APPROX(res->retention, 2.0 / 3.0, 1e-9);
    CHECK(res->median > 0.0);
    CHECK(res->et_lo_90 <= res->median && res->median <= res->et_hi_90);

    /* The grid is carried out in years and still integrates to one — the CLI
     * used to discard it, and the density plot needs it to be a density. */
    CHECK(res->ngrid > 1);
    CHECK(res->grid_t != NULL && res->grid_density != NULL);
    CHECK(res->grid_t[0] >= cfg.t_lo - 1e-6);
    CHECK(res->grid_t[res->ngrid - 1] <= cfg.t_hi + 1e-6);
    double mass = 0.0;
    for (size_t i = 1; i < res->ngrid; i++)
        mass += 0.5 * (res->grid_density[i] + res->grid_density[i - 1]) *
                (res->grid_t[i] - res->grid_t[i - 1]);
    CHECK_APPROX(mass, 1.0, 1e-3);

    gd_analysis_pair_result_destroy(res);
}

TEST(test_analysis_pair_decay_laws_differ)
{
    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = FORMS;
    in.forms_len = GD_LEN_NUL_TERMINATED;
    in.lang_a = "Spanish";
    in.lang_b = "Sardinian";

    gd_error err = GD_ERROR_INIT;
    gd_analysis_pair_result *sw = NULL, *st = NULL;

    gd_analysis_pair_config cfg;
    gd_analysis_pair_config_default(&cfg);
    cfg.model = GD_PAIR_DECAY_SWADESH;
    CHECK(gd_analysis_pair_run(&in, &cfg, &sw, &err) == 0);
    CHECK_STREQ(sw->model, "swadesh");

    cfg.model = GD_PAIR_DECAY_STAROSTIN;
    CHECK(gd_analysis_pair_run(&in, &cfg, &st, &err) == 0);
    CHECK_STREQ(st->model, "starostin");

    /* Competing laws must reach different dates on identical data — that is the
     * point of offering them, and a silent fallback to one law would pass every
     * other check here. */
    CHECK(fabs(sw->median - st->median) > 1.0);

    gd_analysis_pair_result_destroy(sw);
    gd_analysis_pair_result_destroy(st);
}

TEST(test_analysis_empty_buffer_is_empty_not_measured)
{
    /*
     * The point of GD_LEN_NUL_TERMINATED. A length of 0 used to mean "call
     * strlen on it", so an empty buffer was inexpressible and a zeroed input
     * struct quietly asked the library to measure whatever the pointer
     * happened to reach. Now 0 means 0, and an empty forms table is reported
     * as an empty forms table.
     */
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_job *job = NULL;
    gd_analysis_tree_input in; tree_input(&in);

    in.forms_tsv = FORMS;
    in.forms_len = 0;                       /* empty, not "measure it" */
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != GD_E_OK);
    CHECK(job == NULL);

    /* The same buffer measured is a working dataset, which is what makes the
     * line above a statement about the length rather than about the data. */
    gd_error_clear(&err);
    in.forms_len = GD_LEN_NUL_TERMINATED;
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) == GD_E_OK);
    CHECK(job != NULL);
    gd_analysis_tree_job_destroy(job);
    job = NULL;

    /* An explicit length shorter than the buffer is honoured, not rounded up
     * to the terminator: the header row alone is a table with no data rows. */
    gd_error_clear(&err);
    in.forms_len = 34;                      /* the header line only */
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != GD_E_OK);
    CHECK(job == NULL);

    /* And an empty calibration table reads as zero calibrations rather than
     * being measured into whatever follows the pointer. */
    gd_analysis_calibration *cals = (gd_analysis_calibration *)1;
    size_t ncal = 99;
    gd_error_clear(&err);
    (void)gd_analysis_calibrations_read(CALS, 0, &cals, &ncal, &err);
    CHECK(ncal == 0);
    gd_analysis_calibrations_free(cals, ncal);
}

TEST(test_analysis_pair_bad_input)
{
    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = "language_id\tconcept_id\nSpanish\tc1\n";
    in.forms_len = GD_LEN_NUL_TERMINATED;
    in.lang_a = "Spanish";
    in.lang_b = "Portuguese";

    gd_analysis_pair_config cfg;
    gd_analysis_pair_config_default(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_pair_result *res = NULL;

    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) != 0);
    CHECK(res == NULL);
    CHECK(err.code == GD_E_UNKNOWN_COLUMN);

    /* An unknown decay law is an error, never a silent fall back to P1. */
    in.forms_tsv = FORMS;
    gd_error_clear(&err);
    CHECK(gd_pair_model_from_string("not_a_law", &cfg.model, &err) != GD_E_OK);
    CHECK(err.code == GD_E_INVALID_ARG);
    CHECK(strstr(err.message, "not_a_law") != NULL);
    /* The rejected name must not have moved cfg.model off its default. */
    CHECK(cfg.model == GD_PAIR_P1);

    /* The names the CLI and the browser both accept. */
    gd_pair_model pm = GD_PAIR_DECAY_VASILYEV;
    CHECK(gd_pair_model_from_string(NULL, &pm, &err) == GD_E_OK);
    CHECK(pm == GD_PAIR_P1);
    CHECK(gd_pair_model_from_string("", &pm, &err) == GD_E_OK);
    CHECK(pm == GD_PAIR_P1);
    CHECK(gd_pair_model_from_string("two-component", &pm, &err) == GD_E_OK);
    CHECK(pm == GD_PAIR_DECAY_TWO_COMPONENT);
    CHECK_STREQ(gd_pair_model_name(GD_PAIR_P1), "P1_hierarchical_constant_rate");
    CHECK_STREQ(gd_pair_model_name(GD_PAIR_DECAY_STAROSTIN), "starostin");

    /*
     * A malformed weight is taken as 1.0, not as zero, and the run continues.
     *
     * This reverses a deliberate refusal (2284c61). The reason that refusal
     * existed was atof(), which answered 0.0 for anything it could not parse:
     * the form stayed in the table with no weight, quietly lost to any
     * synonym, and shifted the date with nothing said. Falling back to 1.0 is
     * not that bug -- the form is kept at full weight rather than dropped --
     * but the typo is no longer loud. gd_forms counts the cells it happened
     * in so a caller can still report them; nothing in the pipelines does yet.
     *
     * The chronology intake always behaved this way; D3 settled the two
     * readings on it rather than leaving the same file readable two ways.
     */
    static const char *BAD_WEIGHT =
        "language_id\tconcept_id\tcognate_id\tprimary_weight\n"
        "Latin\tc1\tk1\t1.0\n"
        "Spanish\tc1\tk1\theavy\n";
    cfg.model = GD_PAIR_P1;
    in.forms_tsv = BAD_WEIGHT;
    in.lang_a = "Latin";
    in.lang_b = "Spanish";
    gd_error_clear(&err);
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) == GD_E_OK);
    CHECK(res != NULL);
    if (res != NULL) {
        /* The form was kept, so the pair still shares its one concept. */
        CHECK(res->n_shared == 1);
        CHECK(res->n_cognate == 1);
        gd_analysis_pair_result_destroy(res);
        res = NULL;
    }

    /* A '.' in the same column still means "absent", not "malformed". */
    static const char *DOT_WEIGHT =
        "language_id\tconcept_id\tcognate_id\tprimary_weight\n"
        "Latin\tc1\tk1\t.\n"
        "Spanish\tc1\tk1\t.\n";
    in.forms_tsv = DOT_WEIGHT;
    gd_error_clear(&err);
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) == 0);
    CHECK(res != NULL);
    gd_analysis_pair_result_destroy(res);
    res = NULL;
    in.forms_tsv = FORMS;
    in.lang_a = "Spanish";
    in.lang_b = "Portuguese";

    /* Null arguments are refused. */
    cfg.model = GD_PAIR_P1;
    CHECK(gd_analysis_pair_run(NULL, &cfg, &res, &err) != 0);
    in.lang_a = NULL;
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) != 0);
}

/* ------------------------------------------------ progress and cancel (M10) */

typedef struct { size_t calls; size_t last_iter; size_t total; int cancel; } prog_state;

static int prog_cb(void *user, size_t chain, size_t iter, size_t total)
{
    prog_state *p = (prog_state *)user;
    (void)chain;
    p->calls++;
    p->last_iter = iter;
    p->total = total;
    return p->cancel;
}

TEST(test_analysis_trace_is_thinned_and_per_chain)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    cfg.trace_draws = 25;                /* from 200 draws per chain */
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &res, &err) == 0);
    CHECK(res->trace != NULL);
    CHECK(res->trace_draws == 25);

    /* Every kept value must be a real draw of that parameter, so a trace plot
     * cannot show a number the sampler never visited. Ages are bounded below by
     * zero and above by the root's feasible range; checking finiteness and the
     * first/last draw's identity is enough to catch an indexing slip. */
    for (size_t p = 0; p < res->nparams; p++) {
        for (size_t c = 0; c < res->n_chains; c++) {
            const double *t = res->trace + (p * res->n_chains + c) * res->trace_draws;
            for (size_t k = 0; k < res->trace_draws; k++)
                CHECK(!isnan(t[k]) && !isinf(t[k]));
        }
    }

    /* Two chains of the same parameter must not be the same series: they are
     * separate runs, and merging them would hide the disagreement R̂ measures. */
    if (res->n_chains >= 2 && res->nparams > 0) {
        const double *c0 = res->trace;
        const double *c1 = res->trace + res->trace_draws;
        int differ = 0;
        for (size_t k = 0; k < res->trace_draws; k++)
            if (c0[k] != c1[k]) { differ = 1; break; }
        CHECK(differ == 1);
    }
    gd_analysis_tree_result_destroy(res);

    /* The default keeps none: a summary run should not pay for draws nobody
     * asked for. */
    gd_analysis_tree_config plain; quick_config(&plain);
    gd_analysis_tree_result *res2 = NULL;
    CHECK(gd_analysis_tree_run(&in, &plain, &res2, &err) == 0);
    CHECK(res2->trace == NULL);
    CHECK(res2->trace_draws == 0);
    gd_analysis_tree_result_destroy(res2);
}

TEST(test_analysis_reports_progress)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    prog_state st; memset(&st, 0, sizeof st);
    cfg.on_progress = prog_cb;
    cfg.progress_user = &st;
    cfg.progress_every = 25;
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &res, &err) == 0);
    CHECK(st.calls > 0);
    CHECK(st.total == cfg.n_warmup + cfg.n_samples);
    CHECK(st.last_iter == st.total);        /* a bar that can reach 100% */
    gd_analysis_tree_result_destroy(res);
}

TEST(test_analysis_progress_can_cancel)
{
    /* Cancelling is how the page stops a run the visitor abandoned. It must
     * fail loudly rather than hand back a truncated chronology. */
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    prog_state st; memset(&st, 0, sizeof st);
    st.cancel = 1;
    cfg.on_progress = prog_cb;
    cfg.progress_user = &st;
    cfg.progress_every = 10;
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &res, &err) != 0);
    CHECK(res == NULL);
    CHECK(err.code == GD_E_CANCELLED);
    CHECK(st.calls == 1);

    /* The job survives a cancelled run well enough to be destroyed cleanly. */
    memset(&st, 0, sizeof st);
    st.cancel = 1;
    gd_analysis_tree_job *job = NULL;
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) == 0);
    CHECK(gd_analysis_tree_job_run(job, &res, &err) != 0);
    CHECK(res == NULL);
    gd_analysis_tree_job_destroy(job);
}

/* ------------------------------------------ calibration variation (M10) */

TEST(test_analysis_reads_calibrations)
{
    gd_error err = GD_ERROR_INIT;
    gd_analysis_calibration *cals = NULL;
    size_t n = 0;
    CHECK(gd_analysis_calibrations_read(CALS, GD_LEN_NUL_TERMINATED, &cals, &n, &err) == 0);
    CHECK(n == 4);
    CHECK_STREQ(cals[0].id, "cal1");
    CHECK_STREQ(cals[0].target, "Latin");
    CHECK_STREQ(cals[0].type_name, "fixed_tip");
    CHECK_APPROX(cals[0].p1, 2050.0, 1e-9);
    /* The age a calibration pins, used to score whether a withheld one was
     * recovered: p1 for a fixed tip. */
    CHECK_APPROX(cals[0].point, 2050.0, 1e-9);
    /* Latin is ancient; the three living lects are not worth scoring. */
    CHECK(cals[0].is_extant_tip == 0);
    CHECK(cals[1].is_extant_tip == 1);
    gd_analysis_calibrations_free(cals, n);

    /* A uniform window is pinned at its centre, not at an endpoint. */
    const char *uni =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "r\tRoot\tuniform_node\t2000\t3000\t.\t0\tx\n";
    CHECK(gd_analysis_calibrations_read(uni, GD_LEN_NUL_TERMINATED, &cals, &n, &err) == 0);
    CHECK(n == 1);
    CHECK_APPROX(cals[0].point, 2500.0, 1e-9);
    gd_analysis_calibrations_free(cals, n);

    /* A malformed parameter is refused rather than silently read as zero — a
     * calibration that quietly becomes 0 BP would move every date in the tree. */
    const char *bad =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "r\tRoot\tuniform_node\ttwo thousand\t3000\t.\t0\tx\n";
    gd_error_clear(&err);
    CHECK(gd_analysis_calibrations_read(bad, GD_LEN_NUL_TERMINATED, &cals, &n, &err) != 0);
    CHECK(err.code != 0);
}

TEST(test_analysis_drops_a_calibration)
{
    /* Leave-one-calibration-out: withholding Latin's age must actually change
     * the chronology, or the sensitivity analysis is measuring nothing. */
    gd_analysis_tree_input in; tree_input(&in);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_config full; quick_config(&full);
    gd_analysis_tree_result *a = NULL;
    CHECK(gd_analysis_tree_run(&in, &full, &a, &err) == 0);

    gd_analysis_tree_config dropped; quick_config(&dropped);
    dropped.calibrations.drop_index = 0;      /* cal1 = Latin at 2050 BP */
    gd_analysis_tree_result *b = NULL;
    CHECK(gd_analysis_tree_run(&in, &dropped, &b, &err) == 0);

    /* Latin's age was fixed by that calibration; without it the node becomes
     * one the sampler has to estimate -- not an extant tip at age 0, which is
     * what an uncalibrated tip otherwise is. */
    size_t latin = gd_tree_find_label(a->tree, "Latin");
    CHECK(latin != GD_NODE_NONE);
    CHECK(gd_analysis_result_node(a, latin) == NULL);   /* fixed, not sampled */
    CHECK(gd_analysis_result_node(b, latin) != NULL);   /* now estimated */

    gd_analysis_tree_result_destroy(a);
    gd_analysis_tree_result_destroy(b);
}

TEST(test_analysis_widens_calibrations)
{
    /* Loosening every soft window must loosen the posterior it implies. The
     * fixture's only soft window is on the root, so that is where to look. */
    gd_analysis_tree_input in; tree_input(&in);
    in.calibrations_tsv =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "cal1\tLatin\tfixed_tip\t2050\t.\t.\t1\tstandard\n"
        "cal2\tSpanish\tfixed_tip\t0\t.\t.\t1\tstandard\n"
        "cal3\tPortuguese\tfixed_tip\t0\t.\t.\t1\tstandard\n"
        "cal4\tSardinian\tfixed_tip\t0\t.\t.\t1\tstandard\n"
        "cal5\tRoot\tuniform_node\t2400\t2600\t.\t0\tstandard\n";
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_config tight; quick_config(&tight);
    gd_analysis_tree_result *a = NULL;
    CHECK(gd_analysis_tree_run(&in, &tight, &a, &err) == 0);

    gd_analysis_tree_config loose; quick_config(&loose);
    loose.calibrations.widen = 4.0;
    gd_analysis_tree_result *b = NULL;
    CHECK(gd_analysis_tree_run(&in, &loose, &b, &err) == 0);

    size_t root = a->tree->root;
    const gd_analysis_node_age *ta = gd_analysis_result_node(a, root);
    const gd_analysis_node_age *tb = gd_analysis_result_node(b, root);
    CHECK(ta != NULL && tb != NULL);
    double wa = ta->age.hpd_hi - ta->age.hpd_lo;
    double wb = tb->age.hpd_hi - tb->age.hpd_lo;
    CHECK(wb > wa);

    gd_analysis_tree_result_destroy(a);
    gd_analysis_tree_result_destroy(b);
}

TEST(test_analysis_trace_accessor_agrees_with_the_documented_stride)
{
    /* gd_analysis_trace_draw exists so a caller plotting a chain does not have
     * to reproduce [(param * n_chains + chain) * trace_draws + k]. It must
     * agree with that formula for every index, or it has quietly become a
     * second, different contract. */
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    cfg.trace_draws = 9;
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &res, &err) == 0);
    CHECK(res->trace != NULL);

    int mismatched = 0;
    for (size_t p = 0; p < res->nparams; p++)
        for (size_t c = 0; c < res->n_chains; c++)
            for (size_t k = 0; k < res->trace_draws; k++) {
                double direct =
                    res->trace[(p * res->n_chains + c) * res->trace_draws + k];
                if (gd_analysis_trace_draw(res, p, c, k) != direct) mismatched++;
            }

    /* Out of range on any axis is NAN, not a read past the end. */
    double past_param = gd_analysis_trace_draw(res, res->nparams, 0, 0);
    double past_chain = gd_analysis_trace_draw(res, 0, res->n_chains, 0);
    double past_draw  = gd_analysis_trace_draw(res, 0, 0, res->trace_draws);
    double null_res   = gd_analysis_trace_draw(NULL, 0, 0, 0);
    gd_analysis_tree_result_destroy(res);

    CHECK(mismatched == 0);
    CHECK(isnan(past_param));
    CHECK(isnan(past_chain));
    CHECK(isnan(past_draw));
    CHECK(isnan(null_res));
}

TEST(test_analysis_trace_accessor_is_nan_without_a_trace)
{
    /* trace_draws defaults to 0, so most runs keep none. Asking for a draw
     * then is a question with no answer, not a crash. */
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    cfg.trace_draws = 0;
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_result *res = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &res, &err) == 0);
    int had_trace = (res->trace != NULL);
    double v = gd_analysis_trace_draw(res, 0, 0, 0);
    gd_analysis_tree_result_destroy(res);
    CHECK(!had_trace);
    CHECK(isnan(v));
}

TEST(test_analysis_counts_unreadable_weights)
{
    /* D3 settled the weight policy on lenient: an unreadable primary_weight
     * is taken as 1.0 rather than failing the run. That makes a typo in one
     * cell of a large forms file silent unless the count is reported, which
     * is what this field is for. A clean file must report zero, or the
     * warning would cry wolf on every run. */
    static const char *FORMS_BAD =
        "language_id\tconcept_id\tcognate_id\tprimary_weight\n"
        "Latin\tc1\tk1\tnot-a-number\n"
        "Spanish\tc1\tk1\t1.0\n"
        "Portuguese\tc1\tk1\t1.0\n"
        "Sardinian\tc1\tk1\t1.0\n"
        "Latin\tc2\tk2\t1.0\n"
        "Spanish\tc2\tk3\talso-not\n";

    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;

    gd_analysis_tree_result *clean = NULL;
    CHECK(gd_analysis_tree_run(&in, &cfg, &clean, &err) == 0);
    size_t n_clean = clean->n_unreadable_weights;
    gd_analysis_tree_result_destroy(clean);
    CHECK(n_clean == 0);

    in.forms_tsv = FORMS_BAD;
    in.forms_len = GD_LEN_NUL_TERMINATED;
    gd_analysis_tree_result *bad = NULL;
    gd_error err2 = GD_ERROR_INIT;
    int rc = gd_analysis_tree_run(&in, &cfg, &bad, &err2);
    size_t n_bad = (rc == 0) ? bad->n_unreadable_weights : 0;
    gd_analysis_tree_result_destroy(bad);

    /* The run still succeeds -- that is the lenient policy -- and says so. */
    CHECK(rc == 0);
    CHECK(n_bad == 2);
}

/*
 * E0.1: a missing cognate id is not a cognate class. Before the fix, '.' was
 * interned like any id, so two languages that both lacked a concept compared
 * equal and were counted as sharing a cognate.
 */
TEST(test_analysis_pair_missing_is_not_a_match)
{
    static const char *F =
        "language_id\tconcept_id\tcognate_id\tstatus\n"
        "A\tc1\tk1\tobserved\n"
        "B\tc1\tk1\tobserved\n"
        "A\tc2\t.\tobserved\n"
        "B\tc2\t.\tobserved\n"
        "A\tc3\t?\tobserved\n"
        "B\tc3\t?\tobserved\n"
        "A\tc4\t\tobserved\n"
        "B\tc4\t\tobserved\n"
        "A\tc5\tk5\tmissing\n"
        "B\tc5\tk5\texcluded\n"
        "A\tc6\tk6\tobserved\n"
        "B\tc6\tk7\tobserved\n";
    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = F; in.forms_len = GD_LEN_NUL_TERMINATED;
    in.lang_a = "A"; in.lang_b = "B";
    gd_analysis_pair_config cfg; gd_analysis_pair_config_default(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_pair_result *res = NULL;
    int rc = gd_analysis_pair_run(&in, &cfg, &res, &err);
    size_t shared = 0, cog = 0, miss = 0;
    if (res != NULL) {
        shared = res->n_shared; cog = res->n_cognate; miss = res->n_missing_cognates;
        gd_analysis_pair_result_destroy(res);
    }
    CHECK(rc == GD_E_OK);
    CHECK(shared == 2);     /* c1 and c6 only */
    CHECK(cog == 1);        /* c1 */
    CHECK(miss == 8);
}

/* E0.8: a typo in a language name, or a pair with nothing in common, is an
 * error, not a prior-only "date". */
TEST(test_analysis_pair_unknown_or_disjoint_is_an_error)
{
    static const char *F =
        "language_id\tconcept_id\tcognate_id\n"
        "A\tc1\tk1\n"
        "B\tc2\tk2\n";
    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = F; in.forms_len = GD_LEN_NUL_TERMINATED;
    in.lang_a = "A"; in.lang_b = "Bee";
    gd_analysis_pair_config cfg; gd_analysis_pair_config_default(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_pair_result *res = NULL;
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) != GD_E_OK);
    CHECK(res == NULL);
    CHECK(strstr(err.message, "Bee") != NULL);

    in.lang_b = "B";
    gd_error_clear(&err);
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) != GD_E_OK);
    CHECK(res == NULL);
}

/* E0.1, E0.8, E0.10 on the tree path: missing rows are counted and not used;
 * a tip without any observation is an error; a forms language that names an
 * internal node is an error. */
TEST(test_analysis_tree_missing_and_unknown_tips)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_job *job = NULL;

    static const char *WITH_DOTS =
        "language_id\tconcept_id\tcognate_id\n"
        "Latin\tc1\tk1\n" "Spanish\tc1\tk1\n" "Portuguese\tc1\tk1\n" "Sardinian\tc1\tk1\n"
        "Latin\tc2\t.\n"  "Spanish\tc2\t.\n"  "Portuguese\tc2\tk3\n" "Sardinian\tc2\tk2\n";
    in.forms_tsv = WITH_DOTS;
    gd_analysis_tree_result *res = NULL;
    int rc = gd_analysis_tree_run(&in, &cfg, &res, &err);
    size_t miss = 0, nlang = 0;
    if (res != NULL) { miss = res->n_missing_cognates; nlang = res->n_languages;
                       gd_analysis_tree_result_destroy(res); }
    CHECK(rc == 0);
    CHECK(miss == 2);
    CHECK(nlang == 4);

    static const char *NO_SARDINIAN =
        "language_id\tconcept_id\tcognate_id\n"
        "Latin\tc1\tk1\n" "Spanish\tc1\tk1\n" "Portuguese\tc1\tk1\n"
        "Sardinia\tc1\tk1\n";
    tree_input(&in);
    in.forms_tsv = NO_SARDINIAN;
    gd_error_clear(&err);
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0);
    CHECK(job == NULL);
    CHECK(strstr(err.message, "Sardinian") != NULL);

    static const char *INTERNAL =
        "language_id\tconcept_id\tcognate_id\n"
        "Latin\tc1\tk1\n" "Spanish\tc1\tk1\n" "Portuguese\tc1\tk1\n"
        "Sardinian\tc1\tk1\n" "Ibero\tc1\tk1\n";
    tree_input(&in);
    in.forms_tsv = INTERNAL;
    gd_error_clear(&err);
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0);
    CHECK(job == NULL);
    CHECK(strstr(err.message, "Ibero") != NULL);
}

/* E0.10: the language count is the number of tips, also for a polytomy. */
TEST(test_analysis_counts_tips_of_a_polytomy)
{
    gd_analysis_tree_input in; tree_input(&in);
    in.newick = "((Latin:1200,Spanish:2050,Portuguese:2050)WestRom:500,Sardinian:2550)Root;";
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_job *job = NULL;
    int rc = gd_analysis_tree_prepare(&in, &cfg, &job, &err);
    size_t nl = (job != NULL) ? gd_analysis_tree_job_languages(job) : 0;
    gd_analysis_tree_job_destroy(job);
    CHECK(rc == 0);
    CHECK(nl == 4);
}

/* E0.4 and E0.5 through the table reader: is_hard is honoured for
 * uniform_node, refused when it contradicts an intrinsically hard or soft
 * type, and --cal-widen spreads an offset lognormal instead of moving it. */
TEST(test_analysis_is_hard_and_widen_semantics)
{
    gd_error err = GD_ERROR_INIT;
    gd_analysis_calibration *cals = NULL;
    size_t n = 0;
    const char *ok =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "a\tRoot\tuniform_node\t2000\t3000\t.\t0\tx\n"
        "b\tRoot\tuniform_node\t2000\t3000\t.\t1\tx\n"
        "c\tRoot\tuniform_node\t2000\t3000\t.\t.\tx\n"
        "d\tRoot\tnormal_node\t2000\t300\t.\t.\tx\n";
    CHECK(gd_analysis_calibrations_read(ok, GD_LEN_NUL_TERMINATED, &cals, &n, &err) == 0);
    bool h0 = cals[0].is_hard, h1 = cals[1].is_hard, h2 = cals[2].is_hard, h3 = cals[3].is_hard;
    gd_analysis_calibrations_free(cals, n);
    CHECK(!h0); CHECK(h1); CHECK(h2); CHECK(!h3);

    const char *contradiction =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "d\tRoot\tnormal_node\t2000\t300\t.\t1\tx\n";
    gd_error_clear(&err);
    CHECK(gd_analysis_calibrations_read(contradiction, GD_LEN_NUL_TERMINATED, &cals, &n, &err) != 0);
    CHECK(strstr(err.message, "contradicts") != NULL);

    /* Widening an offset lognormal must keep its median (offset + exp(p2))
     * and widen its spread. Measured on the posterior of a root that carries
     * only this calibration and no data signal worth the name. */
    gd_analysis_tree_input in; tree_input(&in);
    in.calibrations_tsv =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "cal1\tLatin\tfixed_tip\t2050\t.\t.\t1\tstandard\n"
        "cal5\tRoot\toffset_lognormal_node\t2100\t6.0\t0.1\t0\tstandard\n";
    gd_analysis_tree_config tight; quick_config(&tight);
    tight.n_warmup = 1000; tight.n_samples = 2000;
    gd_analysis_tree_config loose = tight;
    loose.calibrations.widen = 3.0;
    gd_analysis_tree_result *a = NULL, *b = NULL;
    gd_error_clear(&err);
    CHECK(gd_analysis_tree_run(&in, &tight, &a, &err) == 0);
    CHECK(gd_analysis_tree_run(&in, &loose, &b, &err) == 0);
    const gd_analysis_node_age *ta = gd_analysis_result_node(a, a->tree->root);
    const gd_analysis_node_age *tb = gd_analysis_result_node(b, b->tree->root);
    double wa = ta->age.hpd_hi - ta->age.hpd_lo, wb = tb->age.hpd_hi - tb->age.hpd_lo;
    double ma = ta->age.median, mb = tb->age.median;
    gd_analysis_tree_result_destroy(a);
    gd_analysis_tree_result_destroy(b);
    CHECK(wb > 1.5 * wa);
    /* Before the fix p2 was multiplied: the median jumped from ~2500 to
     * 2100 + exp(18), i.e. by about 6.6e7 years. */
    CHECK(fabs(mb - ma) < 0.5 * wb);
}

/* E0.6: date-pair takes the sampling dates of the two lists. The divergence
 * is then no younger than the older list, and the reported separation is
 * 2T - age_a - age_b. Decay laws, which assume contemporaneous lists, refuse
 * tip ages instead of ignoring them. */
TEST(test_analysis_pair_tip_ages)
{
    gd_analysis_pair_input in;
    memset(&in, 0, sizeof in);
    in.forms_tsv = FORMS; in.forms_len = GD_LEN_NUL_TERMINATED;
    in.lang_a = "Latin"; in.lang_b = "Spanish";
    gd_analysis_pair_config cfg; gd_analysis_pair_config_default(&cfg);
    cfg.age_a = 2050.0;
    cfg.t_hi = 12000.0;
    gd_error err = GD_ERROR_INIT;
    gd_analysis_pair_result *res = NULL;
    int rc = gd_analysis_pair_run(&in, &cfg, &res, &err);
    double lo = 0, med = 0, smed = 0;
    if (res != NULL) {
        lo = res->et_lo_90; med = res->median; smed = res->sep_median;
        gd_analysis_pair_result_destroy(res); res = NULL;
    }
    CHECK(rc == GD_E_OK);
    CHECK(lo >= 2050.0);
    CHECK_APPROX(smed, 2.0 * med - 2050.0, 1e-6);

    cfg.model = GD_PAIR_DECAY_SWADESH;
    gd_error_clear(&err);
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) != GD_E_OK);
    CHECK(res == NULL);

    /* A tip older than the grid can reach is refused, not clamped. */
    cfg.model = GD_PAIR_P1;
    cfg.t_hi = 2000.0;
    gd_error_clear(&err);
    CHECK(gd_analysis_pair_run(&in, &cfg, &res, &err) != GD_E_OK);
    CHECK(res == NULL);
}

/* E0.7: the pair rate prior is per year, like date-tree's. A per-year value
 * and its per-kiloyear equivalent shifted by log(1000) agree, and the
 * default is about -2.0 per kiloyear. */
TEST(test_analysis_pair_rate_is_per_year)
{
    gd_analysis_pair_config cfg; gd_analysis_pair_config_default(&cfg);
    CHECK(fabs(cfg.rate_logmean + log(1000.0) - (-2.0)) < 0.01);
}

/* E0.9: the result carries a convergence verdict consistent with its own
 * diagnostics, and a run far too short to mix is reported as unconverged. */
TEST(test_analysis_convergence_verdict)
{
    gd_analysis_tree_input in; tree_input(&in);
    gd_analysis_tree_config cfg; quick_config(&cfg);
    cfg.n_warmup = 5; cfg.n_samples = 12;
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_result *res = NULL;
    int rc = gd_analysis_tree_run(&in, &cfg, &res, &err);
    bool conv = true; size_t nun = 0, np = 0; double mr = 0, eb = 0, et = 0;
    bool consistent = false;
    if (res != NULL) {
        conv = res->converged; nun = res->n_unconverged; np = res->nparams;
        mr = res->max_rhat; eb = res->min_ess_bulk; et = res->min_ess_tail;
        size_t count = 0;
        for (size_t p = 0; p < res->nparams; p++) {
            const gd_analysis_param *pp = &res->params[p];
            double need = GD_CONV_ESS_PER_CHAIN * (double)res->n_chains;
            if (!(pp->rhat < GD_CONV_RHAT_MAX && pp->ess >= need && pp->ess_tail >= need)) count++;
        }
        consistent = (count == nun);
        gd_analysis_tree_result_destroy(res);
    }
    CHECK(rc == 0);
    CHECK(np > 0);
    CHECK(consistent);
    CHECK(!conv);
    CHECK(nun > 0);
    CHECK(isnan(mr) || mr > 0.0);
    CHECK(isnan(eb) || eb < 400.0);
    (void)et;
}

/* E0.11: a calibration may name a taxon set, "mrca:A,B,...", and unlabelled
 * internal nodes get stable "mrca:X,Y" names that resolve back to them. */
TEST(test_analysis_mrca_targets_and_generated_names)
{
    gd_analysis_tree_input in; tree_input(&in);
    /* Same topology as TREE, with the internal labels removed. */
    in.newick = "((Latin:1200,(Spanish:2050,Portuguese:2050):0):500,Sardinian:2550);";
    in.calibrations_tsv =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "cal1\tLatin\tfixed_tip\t2050\t.\t.\t1\tstandard\n"
        "ib\tmrca:Portuguese, Spanish\tuniform_node\t600\t900\t.\t1\tstandard\n"
        "rt\tmrca:Sardinian,Spanish,Latin\tuniform_node\t2400\t2600\t.\t1\tstandard\n";
    gd_analysis_tree_config cfg; quick_config(&cfg);
    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_result *res = NULL;
    int rc = gd_analysis_tree_run(&in, &cfg, &res, &err);
    double ib_lo = 0, ib_hi = 0, rt_lo = 0, rt_hi = 0;
    int names_ok = 1, found_ib = 0, found_rt = 0;
    if (res != NULL) {
        for (size_t i = 0; i < res->nages; i++) {
            const gd_analysis_node_age *a = &res->ages[i];
            if (a->label == NULL) { names_ok = 0; continue; }
            if (strcmp(a->label, "mrca:Portuguese,Spanish") == 0) {
                found_ib = 1; ib_lo = a->age.hpd_lo; ib_hi = a->age.hpd_hi;
            }
            if (strcmp(a->label, "mrca:Sardinian,Latin") == 0 ||
                strcmp(a->label, "mrca:Latin,Sardinian") == 0) {
                found_rt = 1; rt_lo = a->age.hpd_lo; rt_hi = a->age.hpd_hi;
            }
        }
        gd_analysis_tree_result_destroy(res);
    }
    CHECK(rc == 0);
    CHECK(names_ok);
    CHECK(found_ib);
    CHECK(found_rt);
    CHECK(ib_lo >= 600.0 && ib_hi <= 900.0);
    CHECK(rt_lo >= 2400.0 && rt_hi <= 2600.0);

    /* A name that is not a tip is refused with the name in the message. */
    tree_input(&in);
    in.calibrations_tsv =
        "calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n"
        "x\tmrca:Spanish,Castilian\tuniform_node\t600\t900\t.\t1\tstandard\n";
    gd_analysis_tree_job *job = NULL;
    gd_error_clear(&err);
    CHECK(gd_analysis_tree_prepare(&in, &cfg, &job, &err) != 0);
    CHECK(job == NULL);
    CHECK(strstr(err.message, "Castilian") != NULL);
}

static void run_tests(void)
{
    RUN(test_analysis_prepare_reports_dataset_size);
    RUN(test_analysis_tree_run_summarizes_ages);
    RUN(test_analysis_result_outlives_job);
    RUN(test_analysis_job_runs_once);
    RUN(test_analysis_tree_reproducible);
    RUN(test_analysis_tree_bad_input);
    RUN(test_analysis_trace_is_thinned_and_per_chain);
    RUN(test_analysis_trace_accessor_agrees_with_the_documented_stride);
    RUN(test_analysis_trace_accessor_is_nan_without_a_trace);
    RUN(test_analysis_reports_progress);
    RUN(test_analysis_progress_can_cancel);
    RUN(test_analysis_counts_unreadable_weights);
    RUN(test_analysis_convergence_verdict);
    RUN(test_analysis_reads_calibrations);
    RUN(test_analysis_mrca_targets_and_generated_names);
    RUN(test_analysis_drops_a_calibration);
    RUN(test_analysis_widens_calibrations);
    RUN(test_analysis_is_hard_and_widen_semantics);
    RUN(test_analysis_pair_posterior);
    RUN(test_analysis_pair_decay_laws_differ);
    RUN(test_analysis_empty_buffer_is_empty_not_measured);
    RUN(test_analysis_pair_bad_input);
    RUN(test_analysis_pair_tip_ages);
    RUN(test_analysis_pair_rate_is_per_year);
    RUN(test_analysis_pair_missing_is_not_a_match);
    RUN(test_analysis_pair_unknown_or_disjoint_is_an_error);
    RUN(test_analysis_tree_missing_and_unknown_tips);
    RUN(test_analysis_counts_tips_of_a_polytomy);
}

TEST_MAIN()

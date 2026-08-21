/*
 * fuzz_analysis.c — the whole dataset-ingestion path.
 *
 * gd_analysis_tree_prepare is where the three inputs stop being independent:
 * the forms table is interned against the tree's tip labels, and the
 * calibrations are resolved onto node indices. Those cross-references are
 * where a fuzzer earns its keep — the individual parsers can each be correct
 * while an index derived from one and used against the other is not.
 *
 * It stops before the sampler. gd_analysis_tree_prepare loads, builds and
 * resolves; gd_analysis_tree_job_run is what would spend a minute doing MCMC,
 * and fuzzing it would measure the RNG rather than the parsing.
 *
 * Splitting one byte string into three inputs: the first four bytes are two
 * little-endian 16-bit lengths for the tree and the forms table, and whatever
 * is left is the calibration table. The fuzzer controls the split, so it can
 * mutate the boundaries as well as the contents.
 */
#include "gd.h"

#include <stddef.h>
#include <stdint.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    if (size < 4) return 0;

    size_t n_tree  = (size_t)data[0] | ((size_t)data[1] << 8);
    size_t n_forms = (size_t)data[2] | ((size_t)data[3] << 8);
    const uint8_t *p = data + 4;
    size_t remaining = size - 4;

    if (n_tree > remaining) n_tree = remaining;
    remaining -= n_tree;
    if (n_forms > remaining) n_forms = remaining;
    size_t n_cal = remaining - n_forms;

    gd_analysis_tree_input in;
    in.newick           = (const char *)p;              in.newick_len = n_tree;
    in.forms_tsv        = (const char *)p + n_tree;     in.forms_len  = n_forms;
    in.calibrations_tsv = (const char *)p + n_tree + n_forms; in.cal_len = n_cal;
    in.newick_name = "tree";
    in.forms_name  = "forms";
    in.cal_name    = "calibrations";

    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    /* Small enough that a pathological input cannot spend the fuzzer's budget
     * on arithmetic. prepare() does not sample, but the config still sizes
     * some of what it builds. */
    cfg.n_warmup = 1;
    cfg.n_samples = 1;
    cfg.n_chains = 1;

    gd_error err = GD_ERROR_INIT;
    gd_analysis_tree_job *job = NULL;

    if (gd_analysis_tree_prepare(&in, &cfg, &job, &err) == GD_E_OK) {
        /* Read back what prepare promises a caller before it commits to a run;
         * these are the numbers the CLI prints and the page shows. */
        (void)gd_analysis_tree_job_languages(job);
        (void)gd_analysis_tree_job_nodes(job);
        (void)gd_analysis_tree_job_concepts(job);
        (void)gd_analysis_tree_job_skipped_rows(job);
    }
    gd_analysis_tree_job_destroy(job);

    /* The calibration reader is a separate public entry point over the same
     * bytes, and the sensitivity workflows call it directly. */
    gd_analysis_calibration *cals = NULL;
    size_t ncal = 0;
    gd_error cerr = GD_ERROR_INIT;
    if (gd_analysis_calibrations_read((const char *)p + n_tree + n_forms, n_cal,
                                      &cals, &ncal, &cerr) == GD_E_OK) {
        for (size_t i = 0; i < ncal; i++) {
            (void)cals[i].type;
            (void)cals[i].point;
            (void)cals[i].is_hard;
        }
    }
    gd_analysis_calibrations_free(cals, ncal);
    return 0;
}

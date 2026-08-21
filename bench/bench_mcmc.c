/*
 * bench_mcmc.c — a reproducible number for the sampler's cost.
 *
 * `make benchmark` was a stub from M2 until E1, which meant every performance
 * claim about this engine was an argument rather than a measurement. It exists
 * now because E3 puts an indirection in the sampler's inner loop and that
 * change is owed a before and an after.
 *
 * WHAT IS MEASURED, AND WHY IT IS CPU TIME. clock() is C99, needs no POSIX,
 * and reports processor time rather than wall clock -- which is what the
 * question actually is. A machine with a busy scheduler makes wall clock
 * noisy without making the inner loop slower.
 *
 * WHAT IS NOT MEASURED. An earlier draft counted likelihood calls by wrapping
 * the model callback. That was dropped: wrapping the callback inserts exactly
 * the kind of indirection E3 is trying to measure, so the instrument would
 * have changed the reading. Iterations per second is the number to watch --
 * the calls per iteration are fixed for a given configuration, so the two move
 * together, and this one is measured on the path that actually ships.
 *
 * The run is seeded, so repeats differ only in timing, never in result.
 *
 * THE FIRST RUNS ARE DISCARDED. Without that, the first measured repeat came
 * in 40% slower than the third and the spread across three repeats was 47% --
 * enough to swallow any effect worth measuring. The cause is not the code: a
 * processor idling at a low frequency takes a moment to ramp, and whichever
 * dataset happens to run first pays for it. One discarded run brought the
 * spread to 14%, which was still the ramp and not the code; two brings it
 * under a couple of percent. They are timed and thrown away.
 *
 * If the spread reported at the end is large, do not reach for the median and
 * carry on -- the machine is telling you it cannot resolve the question being
 * asked of it.
 */
#define BENCH_DISCARD 2
#include "gd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static char *read_file(const char *path, size_t *out_len)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL) { fprintf(stderr, "bench: cannot open %s\n", path); return NULL; }
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long sz = ftell(f);
    if (sz < 0 || fseek(f, 0, SEEK_SET) != 0) { fclose(f); return NULL; }
    char *buf = (char *)malloc((size_t)sz + 1);
    if (buf == NULL) { fclose(f); return NULL; }
    size_t rd = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    buf[rd] = '\0';
    *out_len = rd;
    return buf;
}

static double cpu_seconds(clock_t a, clock_t b)
{
    return (double)(b - a) / (double)CLOCKS_PER_SEC;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

int main(int argc, char **argv)
{
    if (argc < 8) {
        fprintf(stderr,
            "usage: bench_mcmc LABEL DIR WARMUP SAMPLES CHAINS SEED REPEATS\n");
        return 2;
    }
    const char *label = argv[1];
    const char *dir = argv[2];
    size_t warmup = (size_t)strtoul(argv[3], NULL, 10);
    size_t samples = (size_t)strtoul(argv[4], NULL, 10);
    size_t chains = (size_t)strtoul(argv[5], NULL, 10);
    unsigned long seed = strtoul(argv[6], NULL, 10);
    size_t repeats = (size_t)strtoul(argv[7], NULL, 10);
    if (repeats == 0) repeats = 1;

    char p_tree[1024], p_forms[1024], p_cal[1024];
    snprintf(p_tree, sizeof p_tree, "%s/tree.nwk", dir);
    snprintf(p_forms, sizeof p_forms, "%s/data/forms.tsv", dir);
    snprintf(p_cal, sizeof p_cal, "%s/calibrations.tsv", dir);

    size_t n_tree = 0, n_forms = 0, n_cal = 0;
    char *s_tree = read_file(p_tree, &n_tree);
    char *s_forms = read_file(p_forms, &n_forms);
    char *s_cal = read_file(p_cal, &n_cal);
    if (s_tree == NULL || s_forms == NULL || s_cal == NULL) {
        free(s_tree); free(s_forms); free(s_cal);
        return 1;
    }

    gd_analysis_tree_input in;
    memset(&in, 0, sizeof in);
    in.newick = s_tree;           in.newick_len = n_tree;
    in.forms_tsv = s_forms;       in.forms_len = n_forms;
    in.calibrations_tsv = s_cal;  in.cal_len = n_cal;
    in.newick_name = p_tree; in.forms_name = p_forms; in.cal_name = p_cal;

    gd_analysis_tree_config cfg;
    gd_analysis_tree_config_default(&cfg);
    cfg.n_warmup = warmup;
    cfg.n_samples = samples;
    cfg.n_chains = chains;
    cfg.seed = (uint64_t)seed;

    double *rate = (double *)calloc(repeats ? repeats : 1, sizeof(double));
    if (rate == NULL) { free(s_tree); free(s_forms); free(s_cal); return 1; }

    int rc = 0;
    /* Warm-up, discarded: see the note at the top of this file. */
    for (size_t r = 0; r < repeats + BENCH_DISCARD; r++) {
        int warming = (r < BENCH_DISCARD);
        gd_error err = GD_ERROR_INIT;
        gd_analysis_tree_job *job = NULL;

        clock_t t0 = clock();
        if (gd_analysis_tree_prepare(&in, &cfg, &job, &err) != GD_E_OK) {
            fprintf(stderr, "bench: prepare failed: %s\n", err.message);
            rc = 1; break;
        }
        clock_t t1 = clock();

        gd_analysis_tree_result *res = NULL;
        if (gd_analysis_tree_job_run(job, &res, &err) != GD_E_OK) {
            fprintf(stderr, "bench: run failed: %s\n", err.message);
            gd_analysis_tree_job_destroy(job);
            rc = 1; break;
        }
        clock_t t2 = clock();

        double prep_s = cpu_seconds(t0, t1);
        double run_s = cpu_seconds(t1, t2);
        double iters = (double)chains * (double)(warmup + samples);
        if (!warming) {
            size_t i = r - BENCH_DISCARD;
            rate[i] = (run_s > 0.0) ? iters / run_s : 0.0;
            printf("%s\t%zu\t%zu\t%zu\t%zu\t%.4f\t%.3f\t%.0f\t%.1f\n",
                   label, i + 1, warmup, samples, chains,
                   prep_s, run_s, iters, rate[i]);
            fflush(stdout);
        }

        gd_analysis_tree_result_destroy(res);
        gd_analysis_tree_job_destroy(job);
    }

    if (rc == 0 && repeats > 0) {
        double *sorted = (double *)malloc(repeats * sizeof(double));
        if (sorted != NULL) {
            memcpy(sorted, rate, repeats * sizeof(double));
            qsort(sorted, repeats, sizeof(double), cmp_double);
            double lo = sorted[0], hi = sorted[repeats - 1];
            double med = sorted[repeats / 2];
            /* Spread as a fraction of the median: an effect smaller than this
             * is not distinguishable from scheduling noise, and a report that
             * quotes one anyway is quoting the machine, not the code. */
            double spread = (med > 0.0) ? (hi - lo) / med : 0.0;
            fprintf(stderr,
                "  %-10s median %.1f iters/s  (min %.1f, max %.1f, spread %.1f%%)\n",
                label, med, lo, hi, 100.0 * spread);
        }
        free(sorted);
    }

    free(rate);
    free(s_tree); free(s_forms); free(s_cal);
    return rc;
}

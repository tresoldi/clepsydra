# Clepsydra

Hierarchical Bayesian lexical dating in portable C99, with an optional
Python reference-prior toolchain.

> **Status: v0.2.0.** Milestones 0–12 complete. The C99 engine implements the
> classical and hierarchical pairwise models (P0/P1), the strict novel-state
> tree likelihood (T1) with a relaxed branch clock, the binary cognate-set
> likelihood (T2), a fixed-tree Metropolis-within-Gibbs sampler with block
> age-scaling, convergence diagnostics and checkpoint/resume, and an adequacy +
> cross-family validation suite with a model card per likelihood
> (`docs/model-cards.md`). The `clepsydra` CLI (`date-tree`, `date-pair`) runs
> end-to-end, with reproducible worked examples on real **IE-CoR** cognate data
> (CC-BY-4.0) under `examples/`. `make test` is 190/190, sanitizer-clean. The
> **browser page** (`web/`) runs the same engine via WebAssembly, with a live
> pairwise tab, a worker-driven chronology, a decay-law comparison and the
> validation evidence; it agrees with the CLI to the last digit. See `PLAN.md`
> for the roadmap and `docs/` for the specification.

## What this is

A reusable scientific library and CLI for estimating language divergence
times from basic-vocabulary cognacy data under several explicit
probabilistic models (classical pairwise, hierarchical pairwise, strict
novel-state multistate tree, binary cognate-set tree), with Bayesian
posteriors, external versioned-TSV reference priors, calibrations for
absolute time, robust contamination for incompatible concepts, and a full
validation suite.

The scientific core is strict C99 (C standard library + libm only) and
compiles natively, and also to **WebAssembly**: `web/` is a browser page that
runs the same engine locally, with no server and no upload. Python (under
`tools/`) is optional and never required by the inference engine.

## Layout

```
include/    public C API (gd.h — the umbrella; gd_analysis.h, gd_options.h,
            gd_report.h, gd_tree.h, gd_types.h, gd_error.h)
src/        inference engine (C99)
cli/        command-line program
web/        browser page + WebAssembly shim (runs the engine in a browser)
tools/      Python reference-prior builder, validation, plotting
reference/  raw sources, generated prior packages, manifests
tests/      unit, likelihood, simulation, regression, fixtures
fuzz/       libFuzzer harnesses over the parsers + committed corpora
bench/      `make benchmark`; docs/benchmark-baseline.md states its resolution
STYLE.md    the conventions this code follows, and why
CONTRIBUTING.md  project state, the CI/analysis/fuzzing gates, reading order
examples/   pair_simple, pair_context, romance_tree, five_ie_reproduction
docs/       model-spec, paper-evaluation, input-formats, calibrations,
            inference, validation
literature/ citations for the required reading (the papers themselves are not
            redistributed; see literature/README.md)
```

`CHANGELOG.md` records what changed per release.

## Reading order for the specification

1. `docs/model-spec.md` — mathematics, parameters, generative models,
   identifiability, the exact pairwise posterior, the fixed-tree MCMC
   target, and pseudocode for the strict novel-state likelihood.
2. `docs/paper-evaluation.md` — how the Rönchen–Wiklund–Hammarström 2024
   paper is used (adopted / adapted / rejected / postponed) and notes on
   the other required references.
3. `docs/input-formats.md` — final TSV schemas + Newick.
4. `docs/calibrations.md` — calibration distributions and feasibility.
5. `docs/inference.md` — pairwise integration and tree MCMC.
6. `docs/validation.md` — test plan, brute-force pseudocode, SBC, adequacy.
7. `docs/model-cards.md` — per-likelihood validity, adequacy, failure regimes.
8. `docs/web.md` — the browser page: architecture, build, verification, deploy.
9. `include/gd.h` — the public C API (one include pulls in the whole library).

## Install

```sh
make                                    # build
sudo make install                       # /usr/local by default
make install PREFIX=$HOME/.local        # or anywhere
make install DESTDIR=/tmp/stage         # staged, for packagers
```

That installs `libclepsydra.a`, a versioned `libclepsydra.so`, the public
headers under `include/clepsydra/`, the `clepsydra` command, and a pkg-config
file. `make uninstall` removes exactly those.

The soname is `libclepsydra.so.0.2` — major *and* minor. This is a 0.x series
and the structs in `include/` are young, so a minor release may break ABI, and
the soname says so rather than promising otherwise.

## Use it as a library

```c
#include <clepsydra/gd.h>     /* one header is the whole public API */

gd_analysis_pair_config cfg;
gd_analysis_pair_config_default(&cfg);

gd_analysis_pair_input in = {0};
in.forms_tsv = forms;                     /* the TSV, in memory */
in.forms_len = GD_LEN_NUL_TERMINATED;
in.lang_a = "latin";
in.lang_b = "spanish";

gd_analysis_pair_result *r = NULL;
gd_error err = GD_ERROR_INIT;
if (gd_analysis_pair_run(&in, &cfg, &r, &err) != GD_E_OK) {
    fprintf(stderr, "clepsydra: %s\n", err.message);
    return 1;
}
printf("%zu shared, %zu cognate, median %.0f BP\n",
       r->n_shared, r->n_cognate, r->median);
gd_analysis_pair_result_destroy(r);
```

```sh
cc myprog.c $(pkg-config --cflags --libs clepsydra) -lm -o myprog
```

Three things worth knowing before depending on it:

- **The library does no I/O.** Every pipeline works from buffers, never paths,
  so the caller reads its own files. That is why the same engine runs natively,
  in a fuzzer, and in a browser with no filesystem.
- **Errors are values, never exits.** Every `int`-returning function returns a
  `GD_E_*` code, and it equals `err.code`. The library never calls `exit()` or
  `abort()`. See the top of `gd_error.h`.
- **A seeded run is bit-reproducible**, and the browser engine agrees with the
  CLI to the last digit — `scripts/wasm-verify.sh` checks it on every build.

## Build and run

```
make native     # strict C99, -Wall -Wextra -Wpedantic -Wconversion ...; builds build/clepsydra
make debug      # -O0 -g, assertions on
make release    # -O2 -DNDEBUG, assertions off
make test       # sanitizer-clean unit + regression + simulation tests (261/261)
make analyze    # gcc -fanalyzer + clang --analyze; gate is zero findings
make fuzz-run   # libFuzzer over the parsers, ASan + UBSan, from committed corpora
make benchmark  # iterations/sec on two examples; see docs/benchmark-baseline.md

# Worked examples on real Lexibank cognate data (CC-BY-4.0; see examples/DATA_SOURCE.md):
examples/pair_simple/run.sh            # P1 pairwise (Spanish–Portuguese, IE-CoR)
examples/romance_tree/run.sh           # T1 Romance chronology (IE-CoR, ~30 s)
examples/five_ie_reproduction/run.sh   # IE tip-dating (deep-time, rate-prior-sensitive)
examples/sinotibetan_tree/run.sh       # T1 Sino-Tibetan chronology (Sagart 2019, a 2nd family)
examples/polynesian_tree/run.sh        # T1 Polynesian chronology (Walworth 2018; shallow, clean)
examples/families_pairwise/run.sh      # P1 pairwise across 9 more families (Dravidian, Turkic, …)
```

See `docs/cli-manual.md` for the command reference and
`docs/scientific-methods.md` for honest interpretation.

```
# The same engine in a browser (no build needed — the artifact is committed):
cd web && python3 -m http.server 8080

source ~/emsdk/emsdk_env.sh && ./web/build-wasm.sh   # rebuild it
./scripts/wasm-verify.sh                             # browser engine == CLI
node scripts/page-smoke.js                           # 23 checks on the live page
```

See `docs/web.md` for how the page is built, verified and deployed.

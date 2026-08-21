# Contributing to Clepsydra

What a new contributor needs that is not in `README.md` or `STYLE.md`: where
the project stands, how the gates work, and the handful of decisions that look
arbitrary until you know why they were made.

`STYLE.md` is the code's contract and takes precedence over anything here.

## Current state

**Milestones 0–9, M11 and M12 are complete** (M0 spec, M1 foundation,
M2 classical + pairwise engine, M3 reference-package builder, M4 hierarchical
pairwise, M5 chronology + calibrations, M6 strict novel-state likelihood T1,
M7 binary cognate-set likelihood T2, M8 fixed-tree posterior sampler, M9 relaxed
clocks + robust observations, M11 adequacy + cross-family validation,
M12 reproducible real-data examples + release). **M10 (WebAssembly + browser
page) was reopened on 2026-08-20 and is the current milestone** — the plan is
PLAN.md §M10: a shared `gd_analysis` runner behind both the CLI and a JSON Wasm
shim, a progress/cancel callback in `gd_mcmc`, and a zero-dependency `web/` page
modelled on the sibling projects (merkmal, regulae, yvyra). See `PLAN.md` for
the roadmap and dependency graph.

**The public API is `include/`, and it is the API the shipped consumers use.**
`gd.h` is an umbrella over `gd_analysis.h` (the pipelines), `gd_tree.h` (the
topology a result is reported against, read-only), `gd_types.h` and
`gd_error.h`. Every function it declares is implemented; nothing returns
`GD_E_NOT_IMPLEMENTED`. The CLI and the Wasm shim build against `include/`
alone — **no `-Isrc`** — and that is enforced by the build rules. If either
ever needs a header from `src/`, the fix is to make that header public, not to
add the flag back.

Internal modules live in `src/gd_*.{c,h}`. Tests are the one consumer allowed
into `src/` (they build with `-Isrc`), because a unit test exists to exercise a
module directly rather than only through the pipelines that call it. `make
test` builds every test against the ASan library and must stay at 100%.

The M0/M1 `gd_context` skeleton — a second, half-implemented API that nothing
shipped ever called — was retired in the C1 alignment milestone; see
`docs/c99-alignment-plan.md`.

M6 notes: `gd_model_novel_state.{c,h}` holds the T1 likelihood — data
extension, the `p_v`/`q_v` recursion, the contamination mixture, simulation,
and `gd_t1_extension` for the debug mode. A `gd_t1_ws` workspace is built once
per topology (it re-reads node ages on every call, so it survives an MCMC
chain) and the likelihood itself allocates nothing. Polytomies are resolved
into zero-length ladders inside the workspace; see `docs/model-spec.md` §3.3
for why that is exact for this model and what remains open.

M7 notes: `gd_model_binary.{c,h}` holds the T2 likelihood — the closed-form
asymmetric binary CTMC (`gd_t2_transition_logs`), Felsenstein pruning in log
space (`gd_t2_loglik_class`), the all-absent ascertainment correction, the
King-2026 meaning-weighted partition rates (`gd_t2_partition_weights`, default
`GD_T2_WEIGHT_MEANING`), the within-meaning dependence warning on the dataset
result, and forward simulation (`gd_t2_simulate_class`). A `gd_t2_ws` workspace
re-reads ages every call like T1's, and the pruning allocates nothing. **T2
needs no polytomy resolution** — Felsenstein pruning combines a node's children
directly, so there is no ladder here (contrast T1). `tests/likelihood/
test_t2_likelihood.c` validates it against an independent full-enumeration
brute force (nothing shared with the pruning), plus normalization, simulation
match, and rate recovery. See `docs/model-spec.md` §3.4.

M8 notes: `gd_mcmc.{c,h}` is the fixed-tree Metropolis-within-Gibbs sampler;
`gd_diagnostics.{c,h}` holds rank-normalized split-R̂ and bulk-ESS. The v1
blocks are node ages (constrained single-site RW; free iff `node_hi>node_lo`)
and a global strict-clock log-rate; the likelihood is reached via a model
callback (`gd_mcmc_model`), and `gd_mcmc_t1_create` is the strict-clock T1
adapter. Adaptation is Robbins-Monro during warm-up only, frozen for retained
draws; per-chain PCG32 streams make runs bit-reproducible; the live state
serializes (versioned + checksum) for chunked resume; draws optionally carry
pointwise per-concept loglik. `tests/simulation/test_mcmc_tree.c` runs strict-
clock-tree SBC (rank uniformity), recovery + R̂, reproducibility, warm-up-only
adaptation, and checkpoint round-trip. See `docs/inference.md` §2.5.

M9 notes: the relaxed clock and robust-observation blocks extend the same
sampler. T1 gained `gd_t1_loglik_concept_relaxed` / `_dataset_relaxed`, which
thread a per-branch multiplier array (`m_b`, indexed by child node) through the
recursion — brute-force validated (`test_t1_relaxed_matches_bruteforce`). The
sampler became a typed-parameter engine (`gd_mcmc_ptype`): the M9 blocks are
branch log-multipliers with a `sigma_b` shrinkage hyperparameter
(`log m_b ~ N(-sigma_b^2/2, sigma_b^2)`, `sigma_b ~ HalfNormal`), contamination
`epsilon` (logit RW, Beta prior), and per-concept rate deviations `mu_c`. All
M9 blocks default **off**, so the M8 strict-clock path (and its tests) is byte
-for-byte unchanged. `gd_mcmc_run` reports a **confounding diagnostic** — the
most negative posterior correlation between a branch's log-multiplier and its
log-duration (docs/model-spec.md §7.2); `branch_rate_confounded` fires past a
0.6 magnitude. `sigma_b` is a pure-prior hyperparameter update (no likelihood
recompute). `tests/simulation/test_mcmc_relaxed.c` covers branch-rate recovery
(ages fixed), the confounding warning (ages free), `sigma_b` shrinkage on
strict data, and robustness to injected loans/miscoding via sampled
contamination. Family multipliers `beta_f` reduce to the global scale on a
single fixed tree and need cross-family data (M11); split-pulse stays excluded.
See `docs/model-spec.md` §4 and `docs/inference.md` §2.5.

M11 notes: adequacy + cross-family validation. `docs/model-cards.md` is the
centerpiece — a card per likelihood (P0, P1, T1, T1-relaxed, T2) separating
implementation validity (SBC/brute-force) from empirical adequacy
(posterior-predictive), listing failure regimes with the diagnostic that
catches each, and stating the flexibility-is-not-correctness principle up
front. `tests/simulation/test_posterior_predictive.c` backs it: a well-specified
fit is PP-calibrated (p≈0.56); a homogeneous-rate fit to concept-rate-
heterogeneous data is flagged inadequate (p≈0, the canonical validity≠adequacy
failure regime); and held-out prediction in a shared observable space shows the
relaxed clock NOT beating the strict clock on strict-truth data ("most flexible"
≠ "best"). `tools/validation/leave_family_out.py` (+ `tools/tests/
test_leave_family_out.py`) is the *synthetic* half of cross-family
generalization: a genuine Monte-Carlo hierarchical leave-one-out coverage
experiment (pooled predictive is consistent, flat fallback over-covers, a
τ̂²-dropping predictive undercovers — it can fail, so a pass means something).
`tools/validation/calibration_benchmark.py` (+ tests) is the *real* half:
leave-one-calibration-out through the actual engine — free an attested tip,
re-date, score whether the withheld age lands in the HPD. Manifest-driven across
Indo-European + Sino-Tibetan + Romance (`examples/loco_benchmark/`, 9 held-out
tips): strict-clock coverage 0.33 vs relaxed 0.89, a real held-out model
comparison (relaxed wins on rate-heterogeneous real data, the opposite of the
strict-truth synthetic result). Note: cross-encoding model
comparison is only ever done in a shared observable space, never via raw
T1-vs-T2 marginal likelihoods. See `docs/validation.md` §4.2/§6/§9.

M12 notes: reproducible real-data examples + release. The CLI is now
end-to-end (`cli/main.c` + `cli/commands.c`, built against the public API in
`include/`): `date-tree` runs the T1 chronology sampler over a Newick tree
+ forms TSV + calibrations TSV and writes posterior_nodes/diagnostics/
run_manifest TSV; `date-pair` runs P1; `calib-sensitivity` runs the
calibration-sensitivity workflows (`--mode loco` = leave-one-calibration-out:
drop each calibration, re-date, score recovery of the withheld age + how far the
rest of the tree moves, flagging load-bearing calibrations whose target goes
`unbounded`; `--mode widen` = loosen every window by a factor and re-date). All three commands run on the shared
`src/gd_analysis.c` pipeline; calibration variation is a config field there
(`calibrations.drop_index` / `.widen`, plus `gd_analysis_calibrations_read` to
enumerate the table), so the browser can use it too. Output is byte-identical to
the pre-refactor `date-tree` and to the committed calib-sensitivity artifacts.
Example + committed artifacts: `examples/romance_tree/calib_sensitivity/`.
`tools/datasets/iecor.py` converts the
IE-CoR CLDF release (Heggarty et al., lexibank `iecor` v1.0, **CC-BY-4.0** —
attribution required, see `examples/DATA_SOURCE.md`) into Clepsydra TSV; only
small converted subsets are committed, not the full dataset. Three examples run
from a clean checkout on real cognate data: `examples/pair_simple` (Spanish–
Portuguese P1), `examples/romance_tree` (Romance T1, proto-Romance ~1600 BP,
converges cleanly), `examples/five_ie_reproduction` (IE tip-dating, a deliberately
hard rate-prior-sensitive deep-time case). A **block age-scaling move** was added
to the sampler (`gd_mcmc.c`, off by default, CLI-enabled) to decorrelate nested
ages — validated by `test_mcmc_block_scaling_same_target` (same target as
single-site) and it lifts the 5-IE ESS ~13→170. `docs/cli-manual.md` and
`docs/scientific-methods.md` document usage and honest interpretation. The C
engine still parses no CLDF/JSON (non-goal v1) — conversion is Python-side.
`tools/datasets/lexibank.py` is the family-agnostic converter (any
Lexibank/LexiCore CLDF Wordlist with cognates); `examples/sinotibetan_tree` uses
it on Sagart et al. 2019 Sino-Tibetan (a second family, CC-BY-4.0) — well-
estimated subgroup ages beside an honestly weakly-identified deep root under a
strict clock. Note: the private `Wu-Urbanek/HM-phylogeny` (Hmong-Mien) repo is
**not reachable** from a tresoldi-scoped session (cross-owner add_repo blocked),
and the public `lexibank/chenhmongmien` wordlist has zero cognate coding, so HM
dating needs cognate-coded data supplied via a session scoped to that repo.
Breadth: `examples/polynesian_tree` is a clean shallow full-tree example
(Walworth 2018, node-calibrated Proto-Polynesian ~2900 BP), and
`examples/families_pairwise` (built by `tools/datasets/build_gallery.py`) runs
P1 across nine more families (Dravidian, Turkic, Uralic, Semitic, Polynesian,
Tupían, Arawakan, Alor-Pantar, Koreanic) — eleven families in all across the
examples, every dataset CC-BY-4.0 and cited in `examples/DATA_SOURCE.md`. For
non-CLDF data (e.g. a hand-built cognate matrix), `tools/datasets/matrix.py`
reads a wide (languages×meanings) or long cognate matrix into the same forms TSV
(tested in `tools/tests/test_matrix.py`), so an external cognate matrix — like a
Hmong-Mien one — plugs straight into `date-tree`/`date-pair`.

Ideas ported from the earlier `distglot` toolkit: `src/gd_decay.{c,h}` adds
alternative pairwise **decay laws** (Swadesh constant-rate, Starostin
accelerating `exp(-2λt²)`, Vasilyev + configurable two-component with a stable
floor), each wrapped in a proper Bayesian grid posterior over the divergence
time. `date-pair --decay swadesh|starostin|vasilyev|two-component` selects one;
default stays the hierarchical constant-rate P1. They are added as *competing
models to compare* (held-out prediction), not defaults — same flexibility≠
correctness stance as the model cards. Tested in `tests/regression/test_decay.c`;
documented in `docs/decay-laws.md` and a card in `docs/model-cards.md`. A second
distglot idea is `tools/analysis/concept_stability.py` (opt-in): data-driven
concept stability = mean pairwise cognate agreement, with an optional emitted
empirical `concept_priors.tsv`. It is explicitly NOT a default and NOT the
reference package — the single-`--forms` priors are in-sample (stamped
`evidence_grade=D`, `EMPIRICAL_INSAMPLE`) with a loud circularity caveat
(`docs/priors.md`). The `--leave-out <family>` mode is the disciplined path: it
pools stability over other families (concepts aligned across families by
Concepticon id from each `concepts.tsv`, equal weight per family) and emits an
out-of-sample prior keyed to the held-out family's own `concept_id`s
(`EMPIRICAL_LFO`) — same discipline as `leave_family_out.py`. Tested in
`tools/tests/test_concept_stability.py`.

## CI

`.github/workflows/ci.yml` runs on every push: `make test` under gcc/ASan and
clang/UBSan, a `-Werror` build with both compilers, every `examples/*/run.sh`
from a clean checkout (plus a check that the chronologies produced node ages),
the Python toolchain (`python3 -m unittest discover -s tools/tests`, 46 tests),
and the web job — provenance, `scripts/wasm-verify.sh` (browser engine vs CLI)
and `node scripts/page-smoke.js` in headless Chrome. `.github/workflows/
pages.yml` deploys `web/` and is gated separately; see `docs/web.md`.

## Static analysis

`make analyze` runs **gcc `-fanalyzer`** and **clang `--analyze`** over
`src/`, `cli/` *and* `web/clepsydra_wasm.c`. Both ship with the compilers CI
already pins, so nothing is installed. **The gate is zero findings with no
suppressions anywhere in the tree**, enforced by a CI job. It is not part of
`make test` — analysis is slower than compilation.

Neither tool exits non-zero on a finding, so the targets grep their output. If
you add an analyzer here, make it fail the same way: a target that reports
"clean" while the tool printed a warning is worse than no target.

`docs/static-analysis.md` records the baseline and the four real defects the
first run found — a null dereference on an out-of-memory cleanup path, a leak
on a partial input load, a dead store hiding a checkpoint-layout invariant, and
a JSON writer that was hard to read. Note the shape: **none of them was on a
path the sanitizer-clean suite reaches.** That is what the
analyzers are for.

## Fuzzing

`make fuzz-run` runs three libFuzzer harnesses (`fuzz/fuzz_{newick,tsv,analysis}.c`)
over committed corpora, under ASan + UBSan; a CI job does the same for 45s per
harness. They call the production parsers directly — the page lets a visitor
paste a Newick tree, so a stranger's bytes reach this C.

Two flags are load-bearing and documented in `docs/fuzzing.md`:
`-fno-sanitize-recover` (without it UBSan prints and the fuzzer exits 0, so a
finding scrolls past a green run) and the unsigned-wrap exclusions the RNG
needs. `make fuzz-run` is verified to fail on a finding rather than assumed to.

A crash becomes an **ordinary test** in `tests/regression/`, not just a corpus
entry: a named test says what the bug was.

## Reading order

1. `STYLE.md` — the conventions this code follows; `docs/c99-alignment-plan.md`
   for why each one is what it is (milestones C1–C8).
2. `PLAN.md` — full plan, milestone dependency graph, open questions.
3. `docs/model-spec.md` — the mathematical specification (parameters,
   generative models, identifiability, exact posteriors, T1 pseudocode,
   brute-force pseudocode, provisional assumptions).
4. `docs/paper-evaluation.md` — adopted/adapted/rejected/postponed split for
   Rönchen et al. 2024, plus notes on the other required references and the
   list of missing literature supplements.
5. `docs/input-formats.md`, `docs/calibrations.md`, `docs/inference.md`,
   `docs/validation.md`, `docs/priors.md`.
6. `docs/static-analysis.md`, `docs/fuzzing.md` — the two analysis gates and
   what they have found.
7. `include/gd.h` — the public C API (one include pulls in the whole library).

## Conventions (this project)

**`STYLE.md` is the full statement**; what follows is the short list worth
having in front of you. Where the two disagree, STYLE.md is right and this
should be fixed.

- **C core:** strict C99, stdlib + libm only. Compile with
  `-std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
  -Wstrict-prototypes -Wmissing-prototypes`. Sanitizer-clean native.
- **Never** `exit()`/`abort()` in `src/` — all errors flow through
  `gd_error*`.
- **One error convention, no exceptions** (documented at the top of
  `include/gd_error.h`): an `int`-returning function returns a `GD_E_*` code,
  and that code equals `err->code`. Never a bare `1` or `-1` — `1` is
  `GD_E_NOMEM`. A function that answers a question returns `bool`, so it is
  visibly not a status. A function returning a pointer fails with `NULL`.
  `gd_error_set`/`_setf` return the code, so a failing branch reads
  `return gd_error_set(err, GD_E_PARSE, "...");`.
- **`gd_error` owns its message.** `gd_error_setf` formats a message naming a
  runtime value straight into the error, so nothing needs a caller-supplied
  scratch buffer. `err.message` is an array: never NULL, `""` when unset —
  test it with `err.message[0]`, not against NULL.
- **`assert()` is for programmer invariants only.** Bad data from a file, a
  user or the browser is always a returned code, never an assertion.
- Initialise with `gd_error err = GD_ERROR_INIT;`, never a positional brace
  list.
- **No raw `malloc`/`calloc`/`realloc` in `src/`** outside `gd_memory.c`. Every
  array allocation goes through `gd_malloc_array_checked` /
  `gd_calloc_checked` / `gd_realloc_array_checked`, which take count and
  element size separately so the multiplication is checked — the counts come
  from parsed input, and `malloc(n * sizeof(T))` is the shape that wraps
  silently. When the count is itself a product (chains x draws x parameters),
  multiply it out with `gd_size_mul` first. `grep -n 'malloc(' src/*.c` should
  find nothing outside `gd_memory.c`.
- **`bool` for boolean semantics, enums for closed sets.** A field that is
  true-or-false is `bool`, not `int`; the only `int`s left that look boolean
  are `gd_error.code` (a status) and the T1 conflict class ids (which carry a
  -1 sentinel). A user-facing name from a closed set is parsed to an enum at
  the boundary — `gd_pair_model_from_string` is the **single** place a
  pairwise-model name is spelled, so the CLI and the browser cannot drift into
  accepting different ones, and an unknown name is `GD_E_INVALID_ARG` rather
  than a silent fall back to the default.
- **Buffer lengths are explicit.** A `*_len` on an analysis input is a byte
  count; `GD_LEN_NUL_TERMINATED` asks the library to measure it. **0 means an
  empty buffer**, not "measure it" — it used to mean the latter, which made an
  empty input inexpressible and a zeroed struct silently ask for a `strlen`.
- **Check every `snprintf` that builds a path or a name.** It cannot overflow,
  but a truncated path names a *different file* — `date-pair` used to report
  success while writing `posteri` instead of `posterior_pair.tsv` for a long
  `--out`. Output paths go through `out_path()` in `cli/commands.c`.
- **Iterative** tree traversal only (no deep recursion — browser stack).
- **Log-space** everywhere in the likelihood; helpers `logsumexp`,
  `log1mexp`, `logdiffexp`.
- **Naming:** library `libclepsydra`, public symbols prefixed `gd_` (as the
  spec's API dictates), repo/project `clepsydra`.

## The browser page

Reopened 2026-08-20; steps 1-7 of the PLAN.md §M10 sequence are done. Full
detail in **`docs/web.md`** — read that before touching `web/`. Working rules:

- **`gd_analysis` is the only place the pipeline lives** (`src/gd_analysis.c`,
  declared in the public `include/gd_analysis.h`). The CLI
  (`cli/commands.c`) adds file I/O and TSV formatting; the shim
  (`web/clepsydra_wasm.c`) adds JSON. Neither may grow inference logic, or the
  page and the CLI start to disagree. `gd_analysis_tree_prepare` /
  `_job_run` split so both callers can report the dataset's size before
  committing to a run.
- **Rebuild and re-verify after touching `src/`, `include/` or the shim:**
  `./web/build-wasm.sh` (writes `web/BUILD_INFO`), then
  `./scripts/wasm-verify.sh` (browser engine vs CLI on a seeded run; observed
  difference is 0) and `node scripts/page-smoke.js` (20 checks in headless
  Chrome). `make test` runs `shim-check`, which compiles the shim's
  non-Emscripten branch natively so it cannot rot without emsdk.
- **Options come from the registry**, `include/gd_options.h` + `src/gd_options.c`
  (D4). The shim's own `opt_slot` table is gone. That table reached its fields
  through a `void *` no compiler could check, and it was wrong once — two
  fields became `bool` while their slots still said `OPT_INT`, so every browser
  run wrote four bytes into a one-byte field and the draws changed; only
  `page-smoke.js` caught it. A registry row states the field's type beside its
  offset and asserts they agree at compile time, so that class of bug now fails
  the build. The shim exposes `clepsydra_defaults_json()` so the page reads the
  engine's defaults instead of restating them.
- **The sampler/likelihood seam is `src/gd_mcmc_model.h`** (E3). `gd_mcmc_state`
  carries a rate *vector*, not a scalar log-rate: T2 has a gain and a loss rate,
  and a scalar is why it could not reach the sampler despite being validated
  since M7. A likelihood declares its own adapter beside itself —
  `gd_mcmc_t1_*` lives in `gd_model_novel_state` now, not in `gd_mcmc.c`.
  `gd_mcmc_t2_*` (E4) is the second adapter, which is what makes the seam real.
  **T2 can be called but not yet sampled:** `fill_state` sets `n_rates = 1`,
  and T2 needs a gain and a loss rate. Sampling both needs a new parameter
  block with its own prior — a modelling decision, not a refactoring one.
- **The calibration-sensitivity verdict is library science** (D6), in
  `gd_analysis_sensitivity_score` and `gd_analysis_divergence_ceiling`. The
  100× ceiling, the shift measure and the covered/unbounded/absent decision
  were `static` functions in `cli/commands.c`: policy no test could reach, and
  the browser could ask for the variant runs through `cal_drop` without being
  able to score them. The CLI keeps the loop and the TSV.
- **The reported tables live in `gd_report`** (D5), not in `fprintf` calls in
  `cli/commands.c`. `posterior_nodes.tsv`, `diagnostics.tsv` and
  `posterior_pair.tsv` are built as tables and rendered to a sink, so their
  columns and number formats are testable in C for the first time. Two of
  those formats are `%.4g`, not `%.4f`: `diagnostics.tsv`'s mean and sd carry
  a log-rate near −9 and an age near 2670 in the same column, so they are
  significant digits. `run_manifest.tsv` deliberately stays hand-written — it
  reports build provenance (compiler, platform, input checksums) that the
  library has no business knowing.
- **Never pass a dataset through `ccall("string")`** — it copies onto the wasm
  stack and a real forms file overflows it, failing as an opaque out-of-bounds
  trap. Heap pointers via `stringToNewUTF8`, then `_free`.
- **`gd_mcmc_config.on_progress` defaults to NULL** and must stay
  draw-for-draw invisible when set; `test_mcmc_progress_does_not_change_draws`
  enforces it.
- The committed `web/clepsydra.js` / `.wasm` are artifacts: the deploy compiles
  nothing, and `scripts/wasm-provenance.sh --check` is what stops a stale engine
  shipping. It hashes `src include web/clepsydra_wasm.c` — **not** `cli/`.
- Page defaults per dataset (`web/data/datasets.js`) are measured settings that
  pass the page's own R̂/ESS thresholds. If a prior or a move changes, re-measure
  rather than relaxing the thresholds.

Remaining: step 8 — the GitHub Pages workflow and `CNAME`, on the user's signal,
once the repository is public.

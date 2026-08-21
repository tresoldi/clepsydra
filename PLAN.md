# Clepsydra — implementation plan

Hierarchical Bayesian lexical dating in portable C99/WebAssembly, with a
Python reference-prior toolchain. This plan is derived from
*Hierarchical Bayesian Lexical Dating Specification.pdf* (44 pp.) in this
directory. It restates the brief's decisions, sequences the work, and flags
open questions and risks. The brief's §20 is binding: **Milestone 0 first,
no broad implementation until the spec is internally consistent.**

Open questions are listed in §12. Decisions marked *(assumed)* below are my
defaults pending your answer; they are reversible.

---

## 1. What we are building (one paragraph)

A reusable scientific library + CLI that estimates language divergence times
from basic-vocabulary cognacy data under several explicit probabilistic
models (classical pairwise, hierarchical pairwise, strict novel-state
multistate tree, binary cognate-set tree), with Bayesian posteriors (grid
integration for pairwise, Metropolis-within-Gibbs MCMC for trees), reference
priors loaded from versioned TSV, calibrations for absolute time, robust
contamination for incompatible concepts, and full diagnostics/validation.
Core in strict C99 (stdlib + libm only), compiles native and to Wasm. Python
is optional and never required by the inference engine.

---

## 2. Fixed product decisions (from the brief — not re-litigated)

- Three commands: `date-pair`, `date-pair-context`, `date-tree`. Shared
  definitions of time, rates, concept priors, contamination.
- Fixed topology in v1; topology-ready architecture (no single-tree
  assumption baked into data/output structures). No topology search.
- Posteriors are mandatory; MLE only for init/debug.
- Reference priors are external versioned TSV, never hardcoded in C. Family
  data *updates* priors, not just selects constants.
- Canonical concept IDs (immutable); external mappings (Concepticon) in
  separate tables, not an architectural dependency. No silent label matching.
- Strict TSV parser in C (path or in-memory buffer; mandatory headers; `.`
  = missing; `#` comments; duplicate/unknown-column detection; **never**
  `exit()`/`abort()` in the library).
- Browser-scale target: ≤100 languages, ≤1000 concepts, multiple forms per
  cell, rooted fixed tree, 4 chains. Single-threaded core; browser runs
  chains in separate Web Workers. Bounded memory, chunked/resumable MCMC,
  deterministic seeds, no deep recursion, no POSIX/thread/file/locale/OS
  entropy dependence.
- Distinguish lexical time (identifiable) from calendar time (needs
  calibration/anchor). Report "years not identifiable" honestly.
- Robust contamination: `L_c = (1-ε_c)L_tree,c + ε_c L_noise,c`; report
  posterior P(concept contaminated).
- Models: P0 classical, P1 hierarchical pairwise, T1 strict novel-state, T2
  binary cognate-set, (later) T3 finite-state Mk/CTMC.
- Clocks: strict → family multiplier → relaxed branch → (later) split-pulse.
- Calibrations: fixed_tip, uniform/normal/lognormal/offset_lognormal node,
  hard/soft min/max. Feasibility checks before inference. Sensitivity:
  prior-only, leave-one-out, one-at-a-time, hard-vs-soft, conflict reporting.
- RNG: PCG32 (or equivalent) with streams, test vectors, unbiased bounded
  ints, normal/exp/gamma/beta/categorical. No `rand()`.
- Checkpointing: `gd_chain_serialize/restore/step`; versioned + checksum.
- Diagnostics: rank-normalized split-R̂, bulk/tail ESS, acceptance, MCSE,
  identifiability flags, prior-to-posterior gain, concept diagnostics,
  pointwise loglik (for WAIC/PSIS-LOO in Python).
- Validation is a release requirement: unit, analytical, brute-force
  likelihood, SBC, parameter recovery, posterior predictive, robustness
  injection, calibration validation, leave-family-out.
- Build: strict warnings (`-std=c99 -Wall -Wextra -Wpedantic -Wconversion
  -Wshadow -Wstrict-prototypes -Wmissing-prototypes`); GCC, Clang,
  Emscripten; sanitizer-clean native.
- Non-goals v1: cognate detection from strings, phonological/sound-change
  models, topology search, networks, borrowing-edge inference, geography,
  pop-size models, CLDF/JSON parsing in C, HMC/NUTS, GPU, required BLAS/GSL.

---

## 3. Repository layout (per brief §14, adopted as-is)

```
include/   gd.h gd_analysis.h gd_tree.h gd_types.h gd_error.h
src/       gd_version.c gd_memory.c gd_math.c gd_rng.c gd_tsv.c
           gd_newick.c gd_tree.c gd_data.c gd_calibration.c gd_priors.c
           gd_pairwise.c gd_model.c gd_model_novel_state.c
           gd_model_binary.c gd_contamination.c gd_quadrature.c
           gd_optimizer.c gd_mcmc.c gd_diagnostics.c gd_output.c
           gd_checkpoint.c
cli/       main.c commands.c
wasm/      gd_wasm.c demo/{index.html,app.js,worker.js}
tools/     reference_builder/ validation/ plotting/
reference/ raw/ generated/ manifests/
tests/     unit/ likelihood/ simulation/ regression/ fixtures/
examples/  pair_simple/ pair_context/ romance_tree/ five_ie_reproduction/
docs/      model-spec.md priors.md input-formats.md calibrations.md
           inference.md validation.md wasm.md paper-evaluation.md
```

Naming *(assumed)*: library `libclepsydra`, public symbols prefixed `gd_`
(as the spec's API already dictates), project/repo name `clepsydra`.

---

## 4. Milestone 0 — specification (do this first, do it completely)

The brief's §20 "first action" checklist, mapped to deliverables:

| # | Requirement | Deliverable |
|---|---|---|
| 1 | Math model spec | `docs/model-spec.md` |
| 2 | Every parameter: symbol, unit, support, prior | table in `model-spec.md` |
| 3 | Generative description of each likelihood | `model-spec.md` per-model subsections |
| 4 | Synonyms, missing data, incompatible concepts per model | `model-spec.md` §handling |
| 5 | Nonidentifiable parameter combinations | `model-spec.md` §identifiability |
| 6 | Exact pairwise posterior definition | `model-spec.md` §pairwise |
| 7 | Exact fixed-tree state + MCMC target | `model-spec.md` §tree-target |
| 8 | Pseudocode for strict novel-state likelihood | `model-spec.md` §novel-state-pseudocode |
| 9 | Independent brute-force validation pseudocode | `docs/validation.md` §brute-force |
| 10 | Final C public API + TSV schemas | `docs/input-formats.md` + `include/gd*.h` skeleton |
| 11 | Provisional assumptions needing later simulation | `model-spec.md` §open-assumptions |
| 12 | Milestone dependency graph | this file §6 |

Plus `docs/paper-evaluation.md` with the Rönchen–Wiklund–Hammarström 2024
ideas split into **adopted / adapted / rejected / postponed**, and an
assumptions table for the strict novel-state model (brief §3.2).

**Acceptance gate (brief §15 M0):** no ambiguous rate/age/path-length
definition; every model generative; absolute-time nonidentifiability
documented; no code copied from external implementations.

**Risk:** M0 depends on reading the 2024 paper. If the PDF is unavailable,
`paper-evaluation.md` is written from the spec's summary + DOI and stamped
"pending primary-source verification" — flagged in §12 Q1.

---

## 5. Milestones 1–12 (sequenced, with acceptance gates)

### M1 — Portable foundation
Error system (`gd_error`), checked memory helpers + context arena, TSV
parser, Newick parser (labels, branch lengths, internal labels, polytomies,
optional simple comments), tree struct + iterative (non-recursive)
postorder/preorder traversal, canonical-ID table, dataset validation, PCG32
RNG with test vectors, numerical helpers (logsumexp, log1p, expm1,
log(1−e^{−x}), stable primitives).
- **Accept:** native + Wasm build; sanitizer-clean; deterministic RNG test
  vectors; parser property/fuzz tests where practical; **no inference yet.**
- **Note:** the TSV parser and Newick parser must both accept in-memory
  buffers (path = `(const char*, size_t)`), since Wasm has no filesystem.

### M2 — Classical + hierarchical pairwise engine (deterministic)
Classical baseline (P0) with explicit `--model classical
--legacy-clock swadesh-100|lees-200` (never silent), P1 with
concept-specific rates, literal two-language mode, missing concepts, uncertain
cognacy (soft-observation `L_c = y_c p_c + (1−y_c)(1−p_c)`), known tip ages,
posterior grid integration (adaptive refinement, log-space normalization, CDF
construction, adaptive Simpson/Gauss–Kronrod, Gauss–Hermite for lognormal
rate if validated), threshold probabilities, pairwise simulation.
- **Accept:** agreement with classical formula; posterior normalization;
  analytical/numerical regression tests; **no MCMC** for normal low-dim
  pairwise.

### M3 — Reference-package builder (Python, `tools/reference_builder/`)
Ingest Pagel 2007 rate data, Vejdemo–Hörberg 2016 supporting data, ASJP
pinned release; concept mapping to canonical IDs; center relative rates in
log space; construct broad global priors; emit strict TSV
(`concept_priors.tsv`, `family_priors.tsv`) + provenance manifests +
checksums; leave-family-out packages.
- **Accept:** reproducible build from pinned inputs; no unexplained
  constants; unknown concepts get valid fallback priors; package passes the
  C validator.
- **Hard dependency:** access to the source datasets (§12 Q1/Q6).

### M4 — Hierarchical pairwise inference
Global concept priors (`μ_c ~ N(m_c^ref, s_c^ref²)`), family clock `β_f`,
family–concept deviation `u_{f,c}`, `σ_f, τ ~ HalfNormal`, global
contamination `ε ~ Beta`, exact synonym marginalization, sensitivity presets
(tight/standard/wide), prior-dominance reporting.
- **Accept:** SBC for pairwise dates and rates; recovery across missingness
  and rate heterogeneity; output distinguishes lexical-time vs calendar-time.

### M5 — Fixed-tree chronology + calibrations
Node-age representation, fixed and ancient tips, all calibration
distributions, feasible-age propagation (parent > child, hard-bound
consistency, feasible lower/upper per node, fail informatively on no feasible
chronology), MRCA queries, chronology init, prior-only sampling.
- **Accept:** calibration contradictions caught pre-MCMC; node order
  maintained; correct target-pair MRCA posterior extraction (for
  `date-pair-context`).

### M6 — Strict novel-state likelihood (T1)
Independent implementation: minimal data extension / equivalent structure,
path-property/extendability checks, concept-wise likelihood, missing
marginalization, internal-node observations, polytomies, robust contamination
wrapper, simulation, debug mode emitting extended-state constraints. User
chooses `incompatible = error|exclude|contamination` (default contamination).
- **Accept:** brute-force agreement on tiny trees; ~linear scaling in
  concepts+nodes for compatible data; incompatible concepts get clear
  diagnostics; no automatic wholesale concept removal.
- **Critical:** validate independently (brief §3.4) — enumerate ancestral
  states, compare to optimized C. Agreement with the authors' implementation
  alone is insufficient.

### M7 — Binary cognate-set likelihood (T2)
Binary encoding (each cognate class = binary trait, naturally handles
synonymy/polymorphism), asymmetric binary CTMC `0 --λ_g--> 1 --λ_l--> 0`,
pruning likelihood, ascertainment correction (all-zero unobserved classes),
shared-observable predictions. Stochastic-Dollo birth-and-loss deferred.
- **Accept:** brute-force agreement; simulation recovery; documented
  within-meaning dependence warning; **no direct cross-encoding model
  comparison** without shared observable space.

### M8 — Fixed-tree posterior sampler
Metropolis-within-Gibbs over blocks: internal node ages, global clock scale,
family multiplier, concept rates/deviations, family–concept deviations,
relaxed branch multipliers, contamination, latent primary synonym states,
uncertain calibration hyperparams. Log-scale RW for positive rates,
logit-scale for bounded proportions, age proposals constrained by
parent/child, clade block scaling, slice sampling for 1-D scales where
stable, categorical Gibbs/Metropolis for latent primary forms. Adapt only in
warm-up; freeze before retained samples; deterministic per-chain seeds from
master seed. Multiple chains, checkpoint/resume, samples + pointwise loglik
output.
- **Accept:** SBC for small strict-clock tree; reproducible chains; no
  post-warm-up adaptation; diagnostics detect deliberately stuck chains.

### M9 — Relaxed clocks + robust observations
Family multipliers, branch multipliers (`log m_b ~ N(−σ_b²/2, σ_b²)` so
E[m_b]=1, shrinkage prior on σ_b), hierarchical concept rates, global +
concept contamination, loan modes (exclude/mixture/sensitivity), calibration
leave-one-out, model sensitivity runner.
- **Accept:** simulated branch-rate recovery where identifiable; warn when
  branch rates × dates confounded; robustness tests with injected loans +
  miscoding.

### M10 — WebAssembly interface + browser page
*Reopened 2026-08-20, after M12. Two halves: the Wasm target proper (the
original scope, whose acceptance gate stands) and a public browser page in the
shape of the sibling projects merkmal / regulae / yvyra — `web/`, GitHub Pages,
`clepsydra.tresoldi.org`.*

**Architecture.** The page talks to a Web Worker; the worker `ccall`s a
JSON-in/JSON-out shim (`web/clepsydra_wasm.c`); the shim calls a new shared
runner (`src/gd_analysis.{c,h}`) that the native CLI calls too. The runner is
the ~240 lines now sitting inside `cmd_date_tree` between option parsing and
output writing: TSV → interned concepts → `gd_cognate_data[]`, calibration
resolution, MCMC wiring, quantile summarising. Extracting it leaves the CLI as
I/O + TSV formatting and the shim as I/O + JSON formatting, and makes it
structurally impossible for the page and the CLI to compute different numbers.
`gd_tsv_parse` and `gd_newick_parse` already take `(path, buf, len)` — the path
is only for diagnostics — so running from memory needs no parser work.

**C-side changes** (the only non-trivial part):
- `src/gd_analysis.{c,h}` — `gd_analysis_tree_config` / `_pair_config` in,
  `gd_analysis_tree_result` / `_pair_result` out; no `stdio`, no paths, errors
  through `gd_error` as everywhere else. The tree result carries per-node median
  + 90% HPD + R̂ + ESS + calibration status + topology + the confounding
  diagnostic + a thinned trace; the pair result carries the whole posterior
  grid, not just the median — the density plot draws it.
- `gd_mcmc_config` gains `on_progress` / `progress_user`; a nonzero return
  cancels the run. Default `NULL`, so the M8/M9 path and its SBC tests stay
  byte-for-byte unchanged. The CLI gets `--progress` from the same callback.
- **Gate on the refactor:** `make test` stays at 150/150 *and* every
  `examples/*/run.sh` reproduces its committed output exactly. One moved digit
  means the refactor is wrong.

**`web/`**, mirroring merkmal / regulae file-for-file: `index.html`,
`style.css` (their shared tokens verbatim), `app.js`, `viz.js`, `worker.js`,
`guide-content.js`, `clepsydra_wasm.c`, `build-wasm.sh`, `BUILD_INFO`, and the
committed `clepsydra.js` / `clepsydra.wasm`. **Zero runtime dependencies** — no
CDN, hand-rolled SVG (yvyra's d3-selection is not adopted): the page works
offline and from any static host.

**Four tabs:**

| tab | engine | latency | notes |
|---|---|---|---|
| Pair | `date-pair` (P1) | 0.3 s | debounced, interactive, no worker |
| Tree | `date-tree` (T1) | ~40–60 s in Wasm | worker + progress + cancel; strict/relaxed toggle |
| Decay | `--decay` laws | instant | four retention laws compared on one pair |
| Validation | none | — | model cards, LOCO, posterior-predictive; read-only |

Web defaults for Tree are 2 chains × (2000 warmup + 2000 samples) ≈ 10–15 s,
with a full run (4 × 8000) on request. Two chains is the floor: a single-chain
default cannot detect non-convergence, and shipping one would contradict §9.

**Four visualizations** (`viz.js`):
1. **Dated tree with HPD bars** — the centrepiece, and unlike yvyra's tree it is
   calendar-time and heterochronous: a years-BP axis with the present at right,
   tips at their own ages (Latin at 2050 BP beside Spanish at 0), a translucent
   90% HPD bar on every internal node, calibrated nodes marked distinctly from
   inferred ones. Download SVG with styles inlined.
2. **Pairwise posterior density** — the full grid as a filled curve with HPD
   shading and the median marked; never a bare number (§9, and
   `docs/scientific-methods.md`).
3. **Diagnostics panel** — R̂ and bulk-ESS bars against their thresholds, trace
   sparklines, and the branch-rate × duration confounding warning.
4. **Held-out calibration plot** — LOCO attested vs inferred with HPD whiskers,
   coverage counted per clock.

**Data.** Example datasets served as static TSV under `web/data/`, fetched on
demand rather than inlined as a JS bundle (~700 KB would block first paint).
CC-BY-4.0 attribution carried onto the page from `examples/DATA_SOURCE.md`.
File upload and a paste-in editor as well.

**Build, provenance, deploy.** `web/build-wasm.sh` drives `emcc` over `src/*.c`
plus the shim (plain Makefile project — no CMake step, unlike the siblings;
emsdk is at `~/emsdk`). `scripts/wasm-provenance.sh --write|--check` after
regulae: it hashes `src/ include/ web/clepsydra_wasm.c` into `BUILD_INFO`, so an
engine older than its sources is a failing CI check rather than a silently stale
deploy. (`cli/` is *not* hashed, contrary to the first draft of this plan: it is
not part of the browser build, so including it would fail the check for changes
that cannot reach the page.) `.github/workflows/pages.yml` gates the deploy on that check plus
asset presence. The repository goes public and `CNAME` lands at that point, not
before.

**Honesty guardrails.** The page carries the boundary section the siblings do,
plus caveats specific to dating: the reduced default iterations are for
responsiveness and a real analysis belongs on the CLI; a wide HPD is the answer,
not a failure; the deep-time cases (5-IE) are rate-prior-sensitive by
construction. The Validation tab exists so the page cannot be read as "type
words in, get a date out".

**Order of work:** (1) `gd_analysis` refactor, gated on identical example
output; (2) progress/cancel callback; (3) shim + build + provenance, verifying
Romance in-browser against the CLI; (4) page shell, Pair tab, worker, Tree tab;
(5) `viz.js`; (6) Decay + Validation tabs, guide content; (7) docs
(`docs/web.md`, README, AGENTS); (8) Pages workflow + CNAME on request.

- **Accept:** no filesystem requirement; no blocking for a whole chain (chunked
  stepping + progress + cancel); identical seeded native-vs-Wasm results within
  a documented FP tolerance, checked on the Romance example; medium-scale memory
  benchmark; `make test` unchanged at 150/150 and every example byte-identical
  after the refactor; the deploy fails on a stale engine; every number shown on
  the page reachable from the CLI.

**Status (steps 1-7 done).** Met so far: `-sFILESYSTEM=0`; progress + cancel
with the draws provably unchanged; native-vs-Wasm difference **0**, not merely
within tolerance (`scripts/wasm-verify.sh`); memory 14-19 MB and wasm ~2x native
(`docs/web.md`); the refactor gate held (all seven examples byte-identical, and
the suite has since grown 150 -> 165); the provenance check fails a stale engine.
Steps 1-8 are done, including the Pages workflow and `CNAME`; publishing waits
only on the repository going public. Full write-up: `docs/web.md`.

The `calib-sensitivity` command (leave-one-calibration-out and window widening)
was written against its own factoring of the pipeline while this milestone was
in flight; it now runs on `gd_analysis` like the other two commands, verified by
reproducing its committed artifacts byte for byte. Calibration variation is
therefore available to the browser as well, and the Tree tab uses the widening
half of it.

### M11 — Adequacy + cross-family validation (full suite)
SBC, posterior predictive checks, parameter recovery, leave-family-out
priors, calibration prediction, model comparison, contamination stress
tests. Use King 2026 principles. Model card per likelihood; distinguish
implementation validity from empirical adequacy; document failure regimes.
- **Accept:** model cards; no "most flexible = best" claim; documented
  failure regimes.

### M12 — Reproducible examples + release
Minimal two-language, pairwise-with-context, Romance fixed-tree,
5-language-IE reproduction-style (not copied), custom-concept-list tutorial,
Concepticon-replacement tutorial, native CLI manual, Wasm integration
manual, scientific-methods document.
- **Accept:** every example runs from clean checkout; inputs/outputs
  versioned; results include uncertainty + diagnostics; no posterior median
  presented as a known true date.

---

## 6. Milestone dependency graph

```
M0 (spec) ──► M1 (foundation) ──► M2 (pairwise deterministic) ──► M4 (hier pairwise)
                  │                       │                              │
                  │                       └──────────────────────────────┤
                  │                                                      ▼
                  ├──► M5 (chronology+cal) ──► M6 (novel-state lik) ──► M8 (tree MCMC) ──► M9 (relaxed+robust)
                  │                              │                          │                  │
                  │                              ├──► M7 (binary lik) ──────┤                  ▼
                  │                              │                          │             M10 (Wasm)
                  │                              ▼                          ▼
                  │                           (brute-force validation feeds back to M6/M7/M8)
                  │
                  └──► M3 (ref-prior builder, Python) ──► feeds M4, M9, M11
                                                       │
                                                       ▼
                           M8 + M9 ──► M11 (full validation) ──► M12 (examples/release)
```

Key ordering constraints:
- **M0 gates everything.** No `src/` implementation until `model-spec.md`
  is internally consistent.
- **M3 can run in parallel with M2** once M1's TSV schemas are frozen (they
  are frozen in M0). M3 is Python-only and independent of the C inference.
- **M6/M7 brute-force validation (M11-precursors) must accompany M6/M7**,
  not be deferred — the brief forbids "it runs = it's implemented."
- **M10 depends on M8** (chunked stepping + checkpoint must exist first).

---

## 7. Cross-cutting engineering rules (apply from M1 onward)

- **No `exit()`/`abort()` in `src/`.** All errors flow through `gd_error *`.
- **Checked multiplications before every allocation;** context arena for
  analysis-lifetime memory; reusable scratch buffers; **no allocation in the
  innermost likelihood loop**; no VLAs depending on user input; no compiler
  extensions.
- **Iterative tree traversal only** (no recursion deep enough to overflow a
  browser stack) — parse to a traversal-order array once.
- **Log-space everywhere:** logsumexp, log1p, expm1, stable `log(1−e^{−x})`,
  explicit NaN/Inf checks, no silent underflow, reproducible summation order.
- **Caches:** concept-local structures, branch durations, traversal orders;
  invalidate explicitly; recompute only affected concepts for concept-rate +
  latent-primary proposals.
- **Determinism:** PCG32 with serialized state; per-chain seeds derived from
  a master seed; identical native↔Wasm within FP tolerance.
- **Sanitizers:** ASan/UBSan/MSan on native CI from M1.

---

## 8. Testing strategy (mapped to validation requirements §13)

| Layer | Where | Examples |
|---|---|---|
| Unit (deterministic) | `tests/unit/` | TSV, Newick, canonical map, calibration feasibility, RNG vectors, distribution moments, log-helpers, tree traversal, MRCA, missing data, synonym-weight normalization, pairwise symmetry, zero path length, extreme dates, very slow/fast concepts |
| Analytical pairwise | `tests/regression/` | classical formula, dense Python grid, posterior normalization, known quantiles |
| Brute-force likelihood | `tests/likelihood/` | enumerate ancestral states on tiny trees, direct summation vs pruning, missing leaf, observed internal node, polytomy, compatible + incompatible patterns, 1 and several cognate classes |
| SBC | `tests/simulation/` | pairwise dates, strict-clock tree, concept rates, family mult, relaxed clocks, contamination, latent synonym-primary |
| Parameter recovery | `tests/simulation/` | grids over n_languages, n_concepts, missingness, n_calibrations, calibration precision, rate heterogeneity, synonym freq, loan rate, family depth |
| Posterior predictive | `tools/validation/` (Python) | cognate-class counts, singletons/doubletons, pairwise retention, compatibility violations, synonym counts, clade/concept retention, change concentration |
| Robustness injection | `tools/validation/` | miscoding, false/missed cognacy, within/distant loans, semantic-shift polymorphism, nonrandom missing, bad calibrations |
| Calibration validation | `tests/regression/` | leave-one-out, predict omitted tips/nodes, conflict diagnostics, shape sensitivity |
| Leave-family-out | `tools/validation/` | global vs family-specific vs legacy vs weak prior; coverage + error |

Tolerances documented per test. SBC rank uniformity is the pass criterion
for "model implemented," not "MCMC runs."

---

## 9. Output contracts (frozen in M0, implemented progressively)

Posterior files (strict TSV): `posterior_nodes.tsv`,
`posterior_rates.tsv`, `posterior_hyperparameters.tsv`,
`posterior_concepts.tsv`, `diagnostics.tsv`, `pointwise_loglik.tsv`,
`run_manifest.tsv`, `samples.bin|tsv`. The run manifest carries version, git
commit, compiler, platform, wasm/native, input checksums, prior-package
version, seed, config, model, reference mode, warnings, start/completion.

Posterior summaries: mean, median, sd, mode (where stable), 50/80/90/95%
equal-tail, 80/90/95% HPD, requested threshold probs (e.g. P(age>3000)),
prior-to-posterior comparison, prior-dominance diagnostic.

Scientific-output tone (brief §17): model used, calibrations setting the
scale, prior package, whether family data updated rates, rate-vs-calibration
uncertainty split, strict-inheritance fit, cross-model sensitivity,
influential concepts/calibrations. No false precision; no "diverged in 1470 BP."

---

## 10. Performance plan (brief §16)

Benchmark configurations: 2×100, 30×200, 100×1000. Record prep time, one
likelihood eval, one MCMC iteration, memory, checkpoint size, native vs Wasm.
Principles already in §7 (caches, iterative traversal, affected-concept-only
recompute). **Numerical correctness is never sacrificed for speed before
validation.**

---

## 11. Provisional assumptions flagged for later simulation (M0 deliverable, §11 of brief §20)

These get an explicit list in `model-spec.md §open-assumptions`; each must be
resolved by an SBC/recovery experiment in M4/M8/M11:

1. Contamination-noise model choice (symmetric Dirichlet-multinomial vs
   empirical categorical) — adequacy unverified.
2. Concept-specific contamination with hierarchical shrinkage (deferred from
   global) — identifiability unverified.
3. Relaxed-branch-clock identifiability vs node ages — known to be
   confounded; needs σ_b posterior + calibration-sensitivity diagnostics.
4. Split-pulse model — **not** in v1 unless simulation shows identifiability.
5. Latent primary-synonym Gibbs updates mixing on real polymorphic data.
6. Gauss–Hermite quadrature for lognormal rate integration — only if
   validated against dense grids.
7. Legacy clock presets as fallback anchors only — never universal.
8. Binary-CTMC ascertainment correction adequacy under heavy polymorphism.
9. Cross-encoding model comparability — only via shared observable space
   (held-out pairwise cognacy / held-out cell states).

---

## 12. Open questions for you (blockers / decisions)

1. **2024 paper + supporting data.** Do you have the Rönchen–Wiklund–
   Hammarström PDF and its data? Without it, `paper-evaluation.md` is
   written from the spec summary + DOI and marked pending verification.
   *(Blocks full M0; I can still produce the rest of M0.)*
2. **Other reference sources** (King 2026, Pagel 2007 rate table, Vejdemo–
   Hörberg 2016 data, pinned ASJP release). Available locally? *(Blocks M3.)*
3. **This session's scope:** write PLAN.md only (done), or proceed to
   execute M0 now?
4. **Repo:** *resolved* — `github.com/tresoldi/clepsydra`, currently private;
   it goes public when the M10 page is ready to deploy. Naming confirmed as
   `clepsydra` / `libclepsydra` / `gd_*`.
5. **Emscripten** — *resolved*: emsdk is installed at `~/emsdk`
   (`source ~/emsdk/emsdk_env.sh`). M10 builds from there; CI does not need it,
   because the artifact is committed and the deploy checks its provenance.
6. **Python toolchain** deps okay (numpy/scipy/pandas, reference-builder
   only)?
7. **Real datasets** for examples (Romance, 5-IE) supplied by you, or
   synthetic fixtures until M12?

---

## 13. Immediate next step (pending your answers)

Execute **Milestone 0**: produce `docs/model-spec.md`,
`docs/paper-evaluation.md`, `docs/input-formats.md` (final TSV schemas),
`docs/validation.md` (incl. brute-force pseudocode), `include/gd*.h` API
skeleton, `docs/inference.md`, and the parameter/assumptions/identifiability
tables — then re-confirm internal consistency before any `src/` code in M1.

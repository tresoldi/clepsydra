# Changelog

Notable changes per release. Dates are release dates; the reasoning behind each
structural change is in the commit that made it, and in `docs/deepening-plan.md`
and `docs/c99-alignment-plan.md`.

## Unreleased

Work on `docs/roadmap.md`. Breaking changes to defaults, the CLI, output
columns and example numbers are expected throughout 0.x and are listed here.

### Changed (breaking)

- **Relaxed clock: a rate-shift move** (the log rate up by δ, every branch
  multiplier down by δ). Relaxed runs that stalled at 4,000 + 4,000 now
  converge: Polynesian R̂ 1.018 → 1.007, five-IE 1.025 → 1.009. Relaxed
  draws differ from earlier versions for the same seed. Checkpoint format 8.

- **`root_var_share_data` and `root_var_share_clock` are a covariance
  split.** Each is now Cov(log root, term) / Var(log root), with log L and
  −log ρ as the terms, so the two sum to 1. They were variance ratios that
  left the covariance out. When a calibration pins the root, log L and log ρ
  co-vary strongly, and the old values grew to 12 and 13 (Romance). Either
  share can still fall outside [0, 1] in that case (five-IE: −0.50 / 1.50).

- **Empirical prior on the concept-rate spread.** The prior is now
  `log sigma_c ~ N(-0.1, 0.35)` (`--concept-rate-prior-logmean`, `-logsd`),
  fitted on 17 calibrated subtrees in six families. It replaces
  HalfNormal(1) and `--concept-rate-prior-sd`. Under the vaguer prior, a
  small family's weakly identified spread came out too low and its lexical
  ages too young: in the simulator, 78% coverage at the 90% level
  (`docs/research/concept-rate-prior.md`).

- **Among-concept rate heterogeneity is on by default** (E2.1;
  `--concept-rates`; the spread prior is described below).
  - **The model.** Each concept's rate is `rho · lambda_c`, with
    `log lambda_c ~ N(-sigma_c²/2, sigma_c²)`, integrated out per concept.
    The multipliers are rescaled to mean exactly 1, so `rho` stays the mean
    rate and is comparable across concept lists.
  - **The quadrature.** A 24-point Gauss–Hermite rule, computed by
    Golub–Welsch. Against a dense integral the error is below 0.4% up to
    `sigma_c = 1.5` on long branches; 8 points erred by 6%, in the slow tail
    that dates deep nodes.
  - **The spread.** `sigma_c` is sampled (its prior is now empirical; see below) (`log_sigma_c` in
    `diagnostics.tsv`). On Romance it is about 1.05, and node ages move by up
    to 2%.
  - **Cost.** Concepts with the same cognate partition are now evaluated
    once and weighted by their count, which is exact (Romance: 170 concepts,
    69 patterns). Even so, the mixture makes runs about 12× slower: the
    Romance example takes about 2 minutes.

- **A proper prior on node ages** (E1.1, E1.3, E1.4, E1.7; `docs/inference.md`
  §2.3). Uncalibrated ages were flat inside a hidden `[0, 1e15]` box: the
  root was tilted towards it, and the posterior was improper. On Arawan the
  root ran to about 9e14 years. The prior is now a log-normal lexical root
  `L = rho · a_root` (`--root-lex-logmean -1.2`, `--root-lex-logsd 1.5`)
  with ages uniform below the root, normalised exactly with dated tips. The
  `1e15` sentinel is gone. Moves:
  - age steps are window-relative, and the root's is on the log scale;
  - the scale moves leave free tips alone;
  - a new clock-preserving move shifts the rate against the ages.

  `--prior-only 1` samples the prior. **Every example's dates change.** The
  Romance root is 2,668 BP; the five-IE root falls from about 7,975 to 6,279
  BP, because the old tilt pushed deep roots towards the box; the
  Sino-Tibetan root presses on the top of its soft 4,500–9,000 window. All
  examples now pass the convergence check, five-IE included, and on Arawan
  the bulk ESS of the root went from 9 to over 3,000.
- **`date-pair` rate prior is per year** (E0.7), the same unit as
  `date-tree`. It was per kiloyear. Subtract log(1000) = 6.9078 from an old
  `--rate-logmean`. The default is -8.9 (about the old -2.0 per kiloyear).

### Added

- **`average`: model averaging by stacking** (E2.8). It combines
  `date-tree` runs that fit different models to the same data. It computes
  PSIS-LOO for each run (`gd_loo.h`, now public), with Pareto-k
  diagnostics and `se_diff` from pointwise differences. It then chooses
  stacking weights and writes model-averaged node ages with a
  between-model share of the log-age variance.
  - Tests: the GPD fit recovers known shapes. PSIS-LOO matches exact LOO
    in a conjugate model (−70.28 against −70.30). Stacking hits a
    closed-form optimum.
  - On Romance, the relaxed clock beats the strict one by 10.9 ± 4.6
    elpd and takes all the weight.
- **Checkpoints bound to their inputs, resumable mid-sampling** (E3.3).
  A chain checkpoint (version 7) records a 64-bit binding hash and refuses
  to restore when anything differs. The hash covers the parameter layout,
  tree, calibration windows, the model's data (via a new optional
  `fingerprint` in `gd_mcmc_model`), every trajectory-shaping setting, and
  the chain stream. The checkpoint also carries the draws already stored,
  with their log-likelihoods and pointwise rows. A resume during sampling
  reproduces every draw of an uninterrupted run (tested). So do the
  refusals of other data, another seed and another chain.
- **Parallel chains.** In threaded builds (the default `THREADS=1`; WASM
  is never threaded) each chain runs in its own POSIX thread, on its own
  tree copy and likelihood clone. `--threads N` caps how many run at once.
  Draws are bit-identical to a sequential run (tested), and
  ThreadSanitizer is clean. That needed `gd_lgamma`, because C99's `lgamma`
  writes the global `signgam`. Four Arawan chains run 2.9× faster.
- **Non-centred move for the relaxed clock.** Branch log-multipliers are
  N(−σ_b²/2, σ_b), a funnel in which single-site moves stall as σ_b shrinks:
  on Arawan the relaxed fit gave R-hat 1.16 and bulk ESS 22 for σ_b. A
  joint move now changes log σ_b and maps every multiplier through its
  standardised deviation (Jacobian n_b·δ). A prior-only test recovers the
  hierarchical prior exactly (σ_b bulk ESS 6,589). Without the move the
  same test fails with R-hat 1.04 and ESS 175. Checkpoint version 6.
- **Incremental likelihood evaluation** (E3.1). The T1 adapter caches, per
  cognate pattern, the data extension (computed once) and every node's p/q
  for each mixture node. An evaluation recomputes only the nodes whose
  edges changed, and their ancestors; a change of ε alone recomputes
  nothing. It agrees exactly with full recomputation over a random walk of
  400 state changes. Edge terms (log retention and log change) depend
  on the rate, mixture node, branch multiplier and tip hazard, not on
  the pattern. They are now computed once per mixture node per
  evaluation, instead of once per pattern: 1.65× on Arawan. The
  dirty-node tracking adds little (1.2× on Romance), because the rate,
  spread, block and clock moves change every edge each sweep.
- **`--write-draws 1`** and **`--write-pointwise 1`** (E3.2) write every
  retained draw (`draws.tsv`) and the per-draw, per-concept
  log-likelihoods for PSIS-LOO and WAIC (`pointwise.tsv`).
- **`date-trees`: dating over a posterior sample of trees** (E3.4).
  - **Input.** Newick or NEXUS (BEAST, MrBayes and RevBayes `Translate`
    tables and `[&…]` comments are handled; `gd_treeset.h`). After burn-in,
    `--ntrees` evenly spaced trees are each dated with the full model.
  - **Output.** Named clades (`--clades`) and the root get their MRCA age
    pooled over trees, with `P(monophyletic)` and the age conditional on
    monophyly, in `posterior_clades.tsv`; `trees.tsv` gives per-tree
    convergence.
  - **Library.** `gd_analysis_clade_draws` and `gd_newick_tip_list` are new.
- **Adequacy panel** (E2.10; `adequacy.tsv`). Posterior predictive checks
  on data replicated from 200 posterior draws. The replicates use each
  concept's rate multiplier, CRP contamination with probability ε, and the
  observed missing cells. The checks are classes per concept, the singleton
  share, the incompatible share, and cognate sharing for every pair of tips.
  Each gives the observed value, the replicated mean and 90% range, and the
  posterior predictive p-value.
- **Subtree-scaling move** (E3.0). For each free internal non-root node, the
  free internal ages beneath it are scaled together, so nested nodes with
  close ages mix. On Romance the bulk ESS of those nodes rose from about 530
  to about 3,750 at 1.4× the cost. The Sino-Tibetan and Polynesian examples
  are back to 5,000 samples per chain. A prior-only test on a caterpillar
  checks the move's Jacobian.
- **Codings, loans and synonyms** (E2.3, E2.4).
  - **Codings.** A forms file with several `analysis_id` codings must be
    given `--analysis` (they used to be mixed).
  - **Loans.** Rows at or above `--loan-threshold` (default 0.5) in
    `loan_probability` are read as missing and counted.
  - **Synonyms in `date-pair`** are marginalised exactly, as the soft match
    `y = Σ w_i w_j [same class]`; `n_synonym_concepts` is reported.
    `gd_forms_candidates` and `gd_forms_read_opts` are new.
- **Gamma concept-rate mixing** (E2.8; `--concept-rate-gamma 1`) as the
  alternative to log-normal. It uses a generalised Gauss–Laguerre rule from
  the same Golub–Welsch solver, with 24 to 64 nodes by shape. Against a
  dense integral the worst error is 0.09%, where 24 nodes at cv 1.5 erred
  by 22%.
- **`cognacy_confidence`** (E2.5) weights shared classes in `date-pair`.
- **Tip observation error** (E2.9; `tip_error` calibration rows). A tip's
  list may replace forms "in error" with a sampled probability e, which has
  a Beta prior. Under T1, e is identified only by that prior, so the docs
  present it as a sensitivity analysis. `logit_tip_error` appears in
  `diagnostics.tsv`, and the adequacy replicates include the error process.
- **`ancestor_tip` calibration** (E2.7): a dated tip on a zero-length branch,
  for an old list that samples an ancestor (e.g. Paumari 1833).
- **Calibrations that can be wrong** (E2.6). An optional `outlier_prob`
  column makes a calibration a mixture with a copy ten times as wide. The
  copy sets no sampler window, so a wrong calibration cannot force an
  impossible chronology. `posterior_calibrations.tsv` reports every
  calibrated sampled node's posterior age and, for outlier-prone rows, the
  posterior probability that the calibration is wrong.
  `gd_cal_widen` and `gd_cal_outlier_posterior` are new.
- **Sampled contamination and a per-concept table** (E2.2).
  - **ε is sampled** with a Beta(1, 19) prior (`--sample-epsilon`, default
    on; `--epsilon` is used with `--sample-epsilon 0`).
  - **`posterior_concepts.tsv`** gives, for every concept, the posterior
    probability that it came from the noise component and its posterior
    rate relative to the mean (from the concept-rate mixture).
  - **Incompatible concepts are excluded.** These are concepts no T1 history
    on the topology can produce. They were charged to the noise component
    before, which changed no date (their term does not depend on ages), but
    inflated the sampled ε. They are now counted (`incompatible_concepts` in
    the manifest) and flagged.
  - **On Romance:** 29 of 170 concepts are incompatible, mostly Latin
    synonyms sorted differently in the daughters. ε falls from 0.29 to about
    0.007 once they are excluded.

- **SBC of the proper prior** (Gate E1). It covers five scenarios: free
  root, normal root calibration, offset-lognormal root calibration, dated
  tip, and Yule. Truths are drawn by exact rejection without the sampler,
  ranks are kept per parameter, and the joint log-likelihood is a test
  quantity. The earlier SBC binned 41 ranks into 8 bins, so bin 0 held six
  ranks and expected 17% more counts; it now uses 40.
- **Width decomposition and power-scaling sensitivity** (E1.6).
  `run_manifest.tsv` reports the share of the calendar root's log-variance
  due to the lexical root (the cognates) and to the clock. On uncalibrated
  Arawan they are 1% and ~100%. It also reports first-order power-scaling
  sensitivities (after Kallioinen et al. 2024) of the root and the lexical
  root to the clock prior, the lexical-root prior, the calibrations and the
  likelihood. Per-draw log-likelihoods are now kept in `gd_mcmc_result`.
- **Yule age prior as a sensitivity alternative** (E1.3; `--yule 1`,
  `--yule-kappa`). This is a pure-birth, completely sampled prior conditioned
  on the root. The rate is written as `kappa / root`, so the prior has no
  calendar unit, and the default kappa matches the expected root for the tip
  count. It is defined only when every tip is modern. On the Polynesian
  example it moves nodes by 1–4% against the uniform prior.
- **Three kinds of age, the scale source and the prior beside the posterior**
  (E1.6).
  - `posterior_nodes.tsv` reports every node's lexical age (`rho · age`,
    replacements per concept) and relative age (age / root age) next to the
    calendar age.
  - It also reports the node under a likelihood-free prior run, which is on
    by default (`--compare-prior`), with the prior–posterior overlap.
  - `run_manifest.tsv` gains `scale_source`, `rate_contraction`,
    `scale_from_prior` and the lexical root.
  - On Arawan the lexical root is the same, 0.54 ± 0.05, with or without a
    root window. Only the calendar ages move, and the manifest says why.
  - Missing numbers in output tables are `.`, not `nan`.
- **Taxon-set targets** (E0.11). A calibration may target `mrca:A,B,...`, the
  most recent common ancestor of the listed tips. Unlabelled internal nodes
  are reported as `mrca:X,Y` instead of `(node)`, and those names resolve
  back to the same node.
- **`date-pair --age-a/--age-b`** (E0.6): the sampling dates of the two
  lists, in years BP. They were hard-coded to 0, so an old list (Paumari
  1833, Arua 1869) could only be dated by correcting the output by hand.
  `posterior_pair.tsv` gains `age_a_bp`, `age_b_bp` and the lineage
  separation `2T - age_a - age_b` (median and 90% interval): for an
  ancestor-descendant pair this is the time elapsed between the lists. The
  decay laws assume contemporaneous lists and refuse tip ages.

- **Convergence check** (E0.9). R̂ is now the larger of the bulk and folded
  (tail) rank R̂. Tail ESS and the MCSE of the mean are new. All of them are
  ported from ArviZ 0.23 and golden-tested against it. The old rank
  normalisation used `N - 1/4` where Vehtari et al. use `N + 1/4`, and its
  Geyer truncation differed. `date-tree` now judges convergence: every
  parameter needs R̂ < 1.01 and bulk and tail ESS ≥ 100 per chain. A failing
  run still writes its outputs, then exits with status 3 unless
  `--allow-unconverged 1`. The verdict is in `run_manifest.tsv` and in the
  result struct (`converged`, `n_unconverged`, `max_rhat`, `min_ess_*`).
  `posterior_nodes.tsv` and `diagnostics.tsv` rename `ess` to `ess_bulk` and
  add `ess_tail` (and `mcse_mean`). The Sino-Tibetan example now runs 8,000
  samples per chain to pass. The five-IE example does not pass at any
  practical length under the current sampler; it runs with
  `--allow-unconverged 1` and says so until E1.7.

### Known gaps

- The committed browser engine (`web/`) predates these changes; it is
  rebuilt with `web/build-wasm.sh`, which needs emsdk.

### Fixed

- **Missing cognate ids were counted as matches** (roadmap E0.1). A `.`, `?`
  or empty `cognate_id`, or a row with `status` `missing`, `excluded` or
  `absent`, was interned as an ordinary cognate class. Two languages that both
  lacked a concept therefore "shared a cognate" in the pairwise models, and
  formed a spurious class in T1. Such rows are now skipped and counted:
  `missing_cognate_rows` in `run_manifest.tsv`, `n_missing_cognates` in both
  result structs.
- **Unknown or disjoint languages no longer "date"** (E0.8). `date-pair` with
  a language the forms file never mentions, or with two languages sharing no
  observed concept, is an error instead of returning the prior. `date-tree`
  fails when a tree tip has no observed form, naming the tip.
- **Uncalibrated tips are extant** (E0.2). A tip that no calibration names is
  fixed at age 0. It used to be free on `[0, 1e15]`, so a "modern" language
  could drift into the past and absorb change that belongs to internal
  branches. Leave-one-out on a tip date now replaces the withheld row with a
  bare minimum of 0, so the tip is still estimated rather than made extant.
- **`normal_node` support** (E0.3) is `mean ± 6 sd` in years. It was
  `exp(mean + 5 sd)`, as if the mean were a log-mean. Lognormal windows are
  cut at 6 log-sd (about 1e-9 of the mass) instead of 5.
- **`--cal-widen` on `offset_lognormal_node`** (E0.4) multiplies the log-sd
  `p3`. It multiplied the log-mean `p2`, which moved the calibration (by
  orders of magnitude) instead of widening it. Soft minimum and maximum
  bounds are now widened too.
- **`is_hard` is honoured** (E0.5). It was read and ignored, so every
  `uniform_node` was hard. `uniform_node` with `is_hard=0` is now a soft
  window with 2.5% tails, and `1` or `.` stays hard. For every other type
  `is_hard` must agree with the type or the table is refused. *Changes the
  romance example's calibration-sensitivity outputs, whose two windows are
  soft.*
- **Sampler start-up** (E0.12). A chain no longer marks the caller's
  calibration set feasible on its own (it did, so an unpropagated or
  contradictory set reached the sampler). A chain whose initial state has a
  non-finite likelihood or prior is refused instead of sampled.
- **Strict command line and input** (E0.14). The CLI refuses a flag with no
  value (it used to skip it and use the default), a flag given twice, and
  a stray argument. A short file read is an error, not a shorter file.
  Forms rows with an empty language or concept id, and `primary_weight`
  values that are NaN, infinite or negative, are refused with their line
  number. `tests/cli/test_cli.sh` runs under `make test`.
- **Public API edge cases** (E0.14). `gd_rng_bounded(rng, 0)` returns 0
  instead of dividing by zero, and bound `0x80000000` no longer hits signed
  overflow. `gd_rng_beta` with tiny shapes no longer returns a point mass at
  0.5 when both gamma draws underflow. `gd_newick_parse` refuses a NULL
  output pointer.
- **Language count and tip resolution** (E0.10). `n_languages` counts tips
  instead of assuming a binary tree (`(nnodes + 1) / 2`). A forms language
  that matches an internal node label is an error.

## 0.2.0 — 2026-08-21

First release intended for use outside the repository. The engine's numbers are
unchanged: every one of the 22 worked-example outputs is byte-identical to
0.1.0-m1 apart from the version line, and the four alternative decay laws whose
moments moved are noted below.

### Added

- **Installation.** `make install` places the static and shared libraries, the
  public headers under `include/clepsydra/`, the `clepsydra` command and a
  pkg-config file. `DESTDIR` and `PREFIX` are honoured; `make uninstall`
  removes exactly what was placed. The soname is `libclepsydra.so.0.2` —
  major *and* minor, because a 0.x minor release may still break ABI.
- **`gd_options.h`** — every run option described once: name, CLI flag, type,
  location and bounds. The CLI, the WebAssembly shim and the browser page read
  the same table, so all three now accept the same options with the same
  validation.
- **`gd_report.h`** — the reported tables and the renderers over them. Output
  goes to a caller-supplied sink, so the library still does no I/O.
- **`gd_analysis_sensitivity_score`** and `gd_analysis_divergence_ceiling` —
  the calibration-sensitivity verdict, previously unreachable inside the CLI.
- **`gd_analysis_trace_draw`** and result accessors, so callers no longer
  reproduce the sampler's index arithmetic.
- **`make benchmark`**, with `docs/benchmark-baseline.md` stating how small an
  effect it can resolve (about 3%).
- T2, the binary cognate-set likelihood, can now be attached to the sampler
  through `gd_mcmc_t2_create`. It cannot yet be *sampled*: that needs a second
  rate parameter, which is a modelling decision.

### Changed

- **`gd_mcmc_state` carries a rate vector** rather than one scalar log-rate.
- **The parsers are buffer-only.** `gd_tsv_parse` and `gd_newick_parse` no
  longer take a path; nothing in `src/` opens a file.
- **One reading of the forms table.** The "highest `primary_weight` wins" rule
  had three implementations that disagreed; an unreadable weight is now taken
  as 1.0 everywhere, and counted in `n_unreadable_weights`. The pairwise
  pipeline used to refuse such a file.
- **Alternative decay laws report slightly different moments.** Sharing one
  posterior-grid summariser moved `mean` by at most 1.0e-2 years and `sd` by
  2.2e-2, on values of 3000–7000 years; medians and every equal-tailed interval
  are unchanged. P1, the default, did not move at all.
- `run_manifest.tsv` reports `unreadable_primary_weights` when there were any.
- The browser payload carries each parameter's trace inside the parameter.

### Fixed

- A heap-buffer-overflow reading a TSV field containing an embedded NUL,
  reachable from a pasted file in the browser. Found by libFuzzer.
- `run_manifest.tsv` could record a commit the binary was not built from.
- The Makefile tracked no header dependencies, so editing a header produced a
  mixed build rather than a rebuilt one.
- The WebAssembly shim truncated error messages into a 192-byte buffer.
- `gd.h` did not include the whole public API it claims to.

## 0.1.0-m1

Milestones 0–12: the pairwise and tree likelihoods, the fixed-tree sampler with
relaxed clocks and robust observations, calibrations, convergence diagnostics,
the validation suite and model cards, the CLI, the reproducible worked examples
on IE-CoR data, and the WebAssembly browser page.

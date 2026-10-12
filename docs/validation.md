# Validation

> Milestone 0 deliverable. Validation is a release requirement, not a later
> enhancement (brief §13). This document defines the test plan and gives
> pseudocode for the independent brute-force likelihood check and for
> simulation-based calibration (SBC). Parameter/model references are to
> `docs/model-spec.md`.

**Status.** Written as the Milestone 0 plan, this document mixes the plan with
what was later built. Where a section says "As implemented", that part is
delivered; the rest is the plan. Delivered against planned, at 0.3:

| § | what | status |
|---|---|---|
| 1, 2 | unit and analytical pairwise tests | delivered |
| 3 | brute-force likelihood, T1 and T2 | delivered |
| 4 | simulation-based calibration | delivered for pairwise and the strict-clock tree sampler (§4.1); relaxed clock and robust observations by exact prior moments and recovery, not full rank-based SBC (§4.2); family multipliers not done |
| 5 | parameter-recovery grids | planned as written; `tools/simworld` covers six single-assumption scenarios on 6–9 language trees, not the full grid |
| 6 | posterior-predictive checks | delivered (`adequacy.tsv`) |
| 7 | robustness injection | partly: borrowing, coding error, rate decay, non-random missingness and gamma rates (`docs/research/simworld.md`); not done: false or missed cognacy separately, loans within clades, erroneous calibrations, the T1 / T1+ε / T2 comparison |
| 8 | calibration validation | partly: leave-one-calibration-out on nine tips and `calib-sensitivity`; no conflict diagnostics beyond the prior-versus-posterior overlap |
| 9 | reference-prior validation | pilot scale (`docs/research/clock-pilot.md`), a no-go for sharp dates |
| 10 | model comparison | delivered (PSIS-LOO, `average`) |
| 10b | browser engine equals native | delivered, for the committed build |
| 11 | model cards | delivered |

Validation philosophy (brief §3.4, §11; King 2026): a model is not
"implemented" merely because its MCMC runs. Implementation **validity**
(SBC, brute-force, analytical) is separate from empirical **adequacy**
(posterior-predictive, cross-family). A model card per likelihood records
both (M11).

## 1. Deterministic unit tests (`tests/unit/`)

- TSV parsing (all schemas in `docs/input-formats.md`): headers, `.`
  missing, `#` comments, blank lines, malformed numeric, duplicate IDs,
  unknown required columns, unknown optional columns, in-memory buffer input.
- Newick parsing: labels, branch lengths, internal labels, polytomies,
  comment handling, malformed input.
- Canonical-ID mapping; alias tables; no silent label matching.
- Calibration feasibility: parent>child consistency, hard-bound conflicts,
  dated-tip fit, feasible bound propagation, no-feasible-chronology failure.
- RNG: PCG32 published test vectors (exact sequences per stream); unbiased
  bounded integers; distribution moments (normal, exp, gamma, beta,
  categorical).
- Log-space helpers: `logsumexp`, `log1p`, `expm1`, `log1mexp`,
  `logdiffexp`; NaN/Inf checks; no underflow.
- Tree traversals: iterative postorder/preorder; MRCA/LCA; no deep recursion.
- Missing-data handling; synonym-weight normalization (sum≈1; default
  `1/n_candidates` with warning); pairwise symmetry; zero path length;
  extreme old dates; very slow and very fast concepts.

## 2. Analytical pairwise tests (`tests/regression/`)

For homogeneous pairwise (P0):
- compare against the classical Swadesh formula (`r = −ln 0.86`/kyr etc.);
- compare numerical integration against dense high-precision Python grids;
- test posterior normalization (∫ posterior = 1);
- test quantiles against analytically/numerically known results.

## 3. Brute-force tree likelihood tests (`tests/likelihood/`)

For tiny trees and small state sets, **independently** enumerate ancestral
histories and compare with the optimized C likelihood. Independent means the
brute-force code does **not** reuse the `p_v`/`q_v` recursion
(`docs/model-spec.md` §11). Pseudocode:

```
function brute_force_loglik_concept(m, tree, rates):
    f = data[m]
    total = -inf  (log)
    for g in all_total_assignments(V, K):      # K = #observed classes + few new
        if not has_path_property(g, tree): continue
        if not compatible(g, f): continue
        pr = 0.0  (log)
        for (u,v) in edges(tree):
            t  = age(u) - age(v)
            rr = exp(-rate(v)*t)
            pr += (g[u]==g[v]) ? log(rr) : log(1-rr)    # generative defn, NOT p/q
        total = logsumexp(total, pr)
    return total
```

Cases required (tolerance documented per test):
- tiny binary tree, compatible data;
- missing leaf;
- observed internal node;
- a polytomy (compare against resolved-ladder approximation until the exact
  nonbinary recursion lands — `docs/model-spec.md` §3.3);
- compatible and incompatible patterns (incompatible ⇒ `L_tree=0`, then
  contamination mixture);
- one and several cognate classes;
- unobserved ancestral state.

Agreement with the authors' implementation alone is **not** sufficient.

### 3.1 As implemented for T1 (M6)

`tests/likelihood/test_t1_likelihood.c` carries **two** independent brute
forces, cross-checked against each other as well as against the library, so
that a single misread of the generative process cannot pass:

- `bf_changes` enumerates the `2^|E|` replacement patterns of §3.3 directly;
- `bf_partitions` enumerates every labeling of the nodes (as restricted growth
  strings, i.e. all set partitions of `V`), keeps those with the path property
  (Prop 3) that match the data, and multiplies the per-edge probabilities.

Neither uses `p_v`/`q_v`, the data extension, or LCA preprocessing. Tolerance
is `1e-9` on the log-likelihood throughout. Covered: every partition of the
tips on a 4-tip binary tree at two rates; one and several classes; missing
tips; no observation at all; observed internal nodes including the root;
an unobserved ancestral class; 3-way and 4-way polytomies and a polytomy under
the root with an ancient tip; 60 randomized cases over five topologies with
random ages, random partial data and random rates (many of them not
extendable); and the `r → 0` / `r → ∞` limits.

Two further checks need no reference implementation at all:

- **normalization** — the tree likelihood sums to 1 over all partitions of the
  observed nodes (binary and polytomous trees, with and without observed
  internal nodes), as does the contamination mixture (§6);
- **simulation** — Monte-Carlo frequencies of each tip partition match the
  computed likelihood within 4 standard errors over 40 000 replicates, and a
  profile likelihood over 300 simulated concepts recovers the generating rate.

Extendability is cross-checked against the brute force too: a pattern is
reported non-extendable exactly when no enumerated history produces it.

### 3.2 As implemented for T2 (M7)

`tests/likelihood/test_t2_likelihood.c` validates the binary cognate-set
likelihood against a single independent brute force, `bf_raw`, that enumerates
all `2^n` present/absent assignments to the nodes, keeps those consistent with
the data, and sums `π(root) · ∏_edges P_{parent,child}(t)` — recomputing the
2×2 CTMC transition matrix from the closed form directly, so it shares nothing
with the Felsenstein pruning or `gd_t2_transition_logs`. The ascertainment
correction is applied identically on both sides (divide by
`1 − P(all-absent)`). Tolerance `1e-9`. Covered: every tip present/absent
pattern on a 4-tip binary tree, with and without ascertainment, at two rate
settings; missing tips; observed internal nodes; a fixed (non-stationary) root
prior; a polytomy handled natively (no ladder); and 60 randomized cases over
five topologies with random ages, random partial data and random gain/loss
rates.

Reference-free checks: the ascertainment-corrected likelihood **sums to 1**
over the observable tip patterns (all but all-absent) and the raw likelihood
sums to 1 over all patterns; the transition matrix rows are distributions and
approach the stationary distribution as `t → ∞`; **simulation** frequencies of
each observable tip pattern match the corrected likelihood within 4 standard
errors (all-absent draws discarded, which is exactly the ascertainment
condition), and a profile-likelihood refit over 400 simulated observable
classes recovers the generating rate scale. Meaning-vs-cognate partition
weighting (King 2026) and the within-meaning dependence warning are tested
directly.

## 4. Simulation-based calibration (`tests/simulation/`)

For each supported Bayesian model, following King 2026:

1. draw parameters from the configured prior;
2. simulate a tree dataset (`simulate` vtable method);
3. fit the same model;
4. calculate posterior ranks of the true values;
5. test **rank uniformity** (RUV) — ECDF difference plots within 95% bands;
6. measure credible-interval coverage (e.g. true value in 95% HPD ~95% of
   the time).

Run SBC separately for:
- pairwise dates (P1);
- strict-clock tree dates (T1, strict clock);
- concept rates (`μ_c`, `r_c`);
- family multipliers (`β_f`);
- relaxed clocks (`m_b`, `σ_b`);
- contamination (`ε`);
- latent synonym-primary states (T1).

A model is not considered implemented merely because MCMC runs.

### 4.1 As implemented for the tree sampler (M8)

`tests/simulation/test_mcmc_tree.c` runs SBC for the **strict-clock tree**: the
quartet `((a,b)r1,(c,d)r2)root` with tips at age 0 and the root anchored at `R`,
whose free parameters are the two internal ages (each with a `Uniform(0,R)`
prior) and the global log-rate (`Normal` prior). Over 80 replicates it draws
the parameters from those priors, simulates T1 data (tips observed only),
refits with `gd_mcmc_run`, and records the rank of every true value among the
posterior draws. The pooled ranks are uniform — a χ² over 6 bins sits well
below the 0.999 critical value (a biased sampler, e.g. a missing Jacobian or
mis-scaled acceptance, would spike it). A separate 4-chain run recovers the
generating rate with rank-normalized split-R̂ < 1.1 and healthy ESS.

The diagnostics themselves are unit-tested in `tests/unit/test_diagnostics.c`:
R̂ ≈ 1 for converged ensembles and R̂ ≫ 1 for a deliberately offset chain
(the "detect stuck chains" gate), bulk-ESS near the total draw count for iid
draws and far below it for an AR(1) sequence, and the inverse-normal quantile
against known values. Reproducibility, warm-up-only adaptation, and
checkpoint/resume identity are asserted alongside the SBC.

### 4.2 Relaxed clock and robust observations (M9)

`tests/simulation/test_mcmc_relaxed.c` covers the M9 acceptance gate on the
quartet under T1:

- **Branch-rate recovery where identifiable** — with the node ages fixed
  (an identifiable regime), simulate concepts under known per-branch multipliers
  and refit; the posterior-mean `log m_b` correlates ~0.98 with the truth.
- **Confounding warning** — with the internal ages free instead, the branch
  rate and the node age trade off; `min_branch_duration_corr ≈ −0.77` and
  `branch_rate_confounded` fires. With ages fixed there is no age variance and
  nothing is flagged.
- **σ_b shrinkage** — on strictly-clocked data (`m_b ≡ 1`) the HalfNormal
  shrinkage prior pulls the posterior-mean `σ_b` toward zero.
- **Robustness to loans / miscoding** — injecting incompatible (loan-like)
  concepts raises the sampled contamination `ε` while the rate estimate stays
  near the truth, because the contamination component of §6 absorbs the loans
  rather than biasing the dates.

The relaxed T1 likelihood itself (per-branch multipliers) is checked against the
independent brute force with per-branch rates in
`tests/likelihood/test_t1_likelihood.c` (§3.1), and the all-ones multiplier is
verified to reproduce the strict clock exactly.

## 5. Parameter-recovery experiments (`tests/simulation/`)

Simulation grids varying: number of languages; number of concepts;
missing-data fraction; number of calibrations; calibration precision;
concept-rate heterogeneity; branch-rate heterogeneity; synonym frequency;
loan/error rate; family depth. Report bias, RMSE, interval width, coverage.

## 6. Posterior-predictive checks (`tools/validation/`, Python)

Simulate replicated datasets from posterior draws and compare (King 2026
metrics, adapted):
- number of cognate classes per meaning;
- singleton and doubleton cognate classes;
- pairwise cognate-retention distribution;
- compatibility violations;
- synonym counts;
- clade-specific retention;
- concept-wise retention;
- rate variation;
- missingness-conditioned statistics;
- concentration of changes on branches.

Essential because a model may recover its own parameters while remaining
inadequate for real lexical data.

**As implemented (M11).** `tests/simulation/test_posterior_predictive.c`
computes the posterior-predictive distribution of a discrepancy statistic — the
variance across concepts of the per-concept distinct-class count — by simulating
a replicate dataset from every posterior draw. A well-specified fit is
calibrated (the observed statistic is typical; PP p-value ≈ 0.56), while a
single-rate model fit to strongly concept-rate-heterogeneous data is flagged
inadequate (PP p-value ≈ 0) even though its sampler is validly implemented (SBC
passes). This is the canonical demonstration that **implementation validity is
not empirical adequacy**. The same file compares the strict and relaxed clocks
by held-out predictive log-score in a shared observable space (held-out cell
states): on strictly-clocked data the relaxed clock does not beat the strict
clock, so flexibility is not rewarded out-of-sample.

## 7. Robustness injection tests (`tools/validation/`)

Starting from model-generated clean data, inject: random miscoding; false
cognacy; missed cognacy; loans within clades; loans across distant branches;
semantic-shift-like polymorphism; nonrandom missing data; erroneous
calibrations. Compare strict (T1), contamination (T1+ε), and binary (T2)
models.

## 8. Calibration validation (`tests/regression/`)

- calibration leave-one-out;
- prediction of omitted dated tips;
- prediction of omitted node constraints;
- calibration conflict diagnostics;
- sensitivity to calibration distribution shape (hard vs soft).

## 9. Reference-prior validation (`tools/validation/`)

Leave-family-out experiments. The scientific question is whether an
out-of-sample clock prior for a held-out family is **calibrated** — does its
leave-one-out predictive interval cover the family's true clock offset at the
nominal rate?

Never evaluate a family-specific prior on the same observations used to
construct it without labelling the result in-sample.

**As implemented (M11).** `tools/validation/leave_family_out.py` (tested in
`tools/tests/test_leave_family_out.py`) runs a Monte-Carlo coverage experiment
on the hierarchical normal model behind a family-clock package
(`β_g ~ N(μ, τ²)`, stored estimate `y_g ~ N(β_g, ω_g²)`). Holding each family
out in turn, it builds three priors for the held-out family's true offset and
measures how often the true value lands in each interval, across thousands of
independent draws:

| prior | what it is | measured coverage (τ_true 0.4, 10 families) |
|---|---|---|
| `family_specific` | the family's own estimate (in-sample, **cheating**) | ≈ nominal — a positive control confirming the harness |
| `pooled_loo` | honest hierarchical leave-one-out predictive (grand mean, between-family variance τ̂² by method of moments, Student-t tails) | ~0.82 at nominal 0.90 — **consistent**: climbs to ≈0.89 by 150 families |
| `unknown_fallback` | the flat UNKNOWN row the builder currently emits (fixed SD 0.80) | ~1.0 — over-covers when τ_true < 0.80 |

This is a genuine, falsifiable experiment, **not** a construction that cannot
fail: coverage is a measured Monte-Carlo number that can (and does) miss
nominal. Its honest findings: the flat fallback is over-conservative; the
naive-but-principled pooled predictive *modestly undercovers with few families*
(a real small-sample effect, bounded and documented, that vanishes as families
accumulate — the estimator is consistent), and sits closer to nominal than the
flat fallback. The experiment has **teeth**: a predictive that drops the
between-family variance `τ̂²` collapses to ≈0.50 coverage at nominal 0.90
(`run_coverage(broken=True)`), which the tests assert it must.

**Scope.** This validates the leave-one-out *procedure* on synthetic
hierarchical data — it certifies the uncertainty-propagation arithmetic is
honest. The complementary real-world check is the leave-one-calibration-out
benchmark below, which *can* catch model misspecification.

### Leave-one-calibration-out (LOCO) — real, end-to-end

`tools/validation/calibration_benchmark.py` (tested in
`tools/tests/test_calibration_benchmark.py`) runs the second experiment on real
cognate data through the actual C engine. On a heterochronous tree with attested
ancient tips, each tip calibration is held out in turn: the calibration file is
rewritten with that tip **freed** (fixed age → broad uniform prior, so the
sampler tip-dates it), `date-tree` reruns with the rest fixed, and the withheld
age is scored against the freed node's 90% HPD. The target is a value already in
the calibration file — nothing is fabricated; the test is recovery of a datum
the model was not given.

Run across **three families** (`examples/loco_benchmark/`, N = 9 attested tips:
Indo-European ×5, Sino-Tibetan ×3, Romance ×1), strict vs relaxed clock. Each
withheld date is scored by coverage and by the interval score of the 90%
interval (width plus 20× any miss distance). The baseline is the same run's
prior-only interval for the freed tip; skill = 1 − IS / IS(prior).

| clock | covered | Indo-European | Sino-Tibetan | Romance | RMSE | interval score (prior) | skill |
|---|---|---|---|---|---|---|---|
| strict  | 1/9 | 1/5 | 0/3 | 0/1 | 1,254 BP | 7,864 (3,930) | −1.00 |
| relaxed | 7/9 | 4/5 | 3/3 | 0/1 | 1,087 BP | 4,407 (3,907) | −0.13 |

Findings the synthetic study cannot produce. **(1)** On real data the strict
clock is overconfident: narrow intervals that miss 8 of 9, scoring twice as
badly as the prior-only interval. The relaxed clock covers 7 of 9 with smaller
errors, the opposite of the strictly-clocked synthetic check in
`docs/model-cards.md`, where relaxing does not help. **(2)** The relaxed clock
is less overconfident but only modestly informative: it beats the prior on
Sino-Tibetan (skill 0.47, from a run above the 1.01 R̂ threshold), ties on
Indo-European, and over the eight recoverable targets has skill +0.13. With nine
targets from three families, this shows less overconfidence than the strict
clock, not calibration. **(3)** Romance-Latin is structural: Latin is the family's only
ancient tip, so freeing it leaves no tip-age signal. **(4)** Deep tips (Vedic,
Old Chinese) carry the largest errors under both clocks: freeing them lets
their age trade against the root's.

### Blind dating on GLED against published dates and ASJP

`tools/gled_benchmark/`, report `docs/research/gled-benchmark.md`. Clepsydra
dates 29 known-age nodes (18 with two-sided windows) and 12 family roots from
GLED's 40-concept automatic cognates, with no calibration, and the ASJP method
(Holman et al. 2011 constants, ASJP v21) dates the same nodes. Summary:

| | nodes: MAE | roots: MAE | roots: median within published window |
|---|---|---|---|
| clepsydra, default prior | 823 y | 3,191 y | 3 / 12 |
| ASJP (in-sample constants) | 541 y | 1,535 y | 5 / 12 |

Clepsydra is biased young at depth (log error −0.73 on nodes older than 2,500
BP), its calendar intervals are prior-dominated (rate contraction ≈ 0) and
too wide to inform, a leave-one-family-out rate prior did not help, and
shallow calibrations do not extrapolate to deep roots (Indo-European
1,529 BP against ~7,000 published). ASJP's graded distance orders the nodes
better (ρ = 0.82) than binary automatic cognates (0.51) or clepsydra's
lexical age (0.30). Eighteen nodes and twelve roots, automatic cognates and
UPGMA topologies: a stress test, not a calibration study.

These are **small-N spot checks** (nine held-out points), not coverage
guarantees; the attested dates carry their own windows a point comparison
ignores. Adding a family is a few manifest lines; the constraint is curating
trees with trustworthy attested/published calibrations, tracked in
`docs/paper-evaluation.md`.

## 10. Model comparison

Output pointwise log likelihood by concept so Python can compute WAIC /
PSIS-LOO where diagnostics permit; held-out-concept / held-out-language /
calibration prediction. When data encodings differ, compare predictions in a
**common observable representation** (held-out pairwise cognacy or held-out
language–concept states). Do **not** compare raw AIC/BIC/WAIC/marginal
likelihood across different data codings without explicit justification
(`docs/model-spec.md` §3.4).

## 10b. The browser engine and the page (M10)

Two things the C test suite cannot reach, checked separately.

**The browser engine must equal the native one.** `scripts/wasm-verify.sh` runs
the Romance chronology (2 chains x 500+500, seed 42) and the
Spanish-Portuguese pairwise date through both builds and compares every
reported quantity — every node's mean, median, sd, 90% bounds, R-hat and ESS,
and the whole pairwise summary. The tolerance is relative and documented
(`TOL`, default 1e-9); the observed difference is **0**. This is what licenses
`docs/web.md`'s claim that a date on the page is the date the CLI computes: the
two share `src/gd_analysis.c`, and the check confirms the sharing is real
rather than nominal.

**The page must do what it says.** `node scripts/page-smoke.js` drives headless
Chrome over the DevTools Protocol and asserts 20 properties of the live page.
The ones that carry scientific weight:

- the pairwise date shown equals one obtained by *running the CLI* during the
  test, not a constant pasted into it;
- a deliberately short run is **reported as unconverged**, and the shipped
  per-dataset defaults **do** converge under the page's own R-hat <= 1.01 and
  ESS >= 200 thresholds (so a default that failed its own check would fail the
  test rather than teach visitors to ignore the warning);
- ancient tips are drawn to the left of modern ones on the dated tree — a
  cladogram that ignored the ages would pass every other check;
- the five decay laws produce at least four distinct medians, so a silent
  fallback to one law cannot masquerade as agreement;
- both clocks render in the held-out calibration figure, with coverage reported
  per clock;
- Cancel leaves no partial posterior on screen.

Neither script runs in `make test` (one needs emsdk-built artifacts, the other a
browser). `make test` does run `shim-check`, which compiles the WebAssembly
adapter's non-Emscripten branch with the native compiler, so the browser surface
cannot break unnoticed on a machine without Emscripten.

## 11. Model cards (M11)

One card per likelihood (P0, P1, T1, T1-relaxed, T2, reserved T3): assumptions;
generative description; parameters; identifiability; SBC status;
posterior-predictive adequacy; documented failure regimes; clear distinction
between implementation validity and empirical adequacy; explicit statement that
greater flexibility does **not** imply greater realism.

**Delivered** in `docs/model-cards.md`. Every "validity" claim names a test in
`make test`; every "adequacy" claim is a posterior-predictive or recovery check
with its limits stated; every "failure regime" pairs a data situation with the
diagnostic that catches it (or says none exists yet, and withholds the result —
e.g. the exact non-binary T1 recursion). The flexibility fallacy is stated at
the top and demonstrated by the held-out model comparison (§6).

# Model cards

> Milestone 11 deliverable. One card per likelihood. Following King 2026
> (doi:10.1371/journal.pcbi.1014312), each card separates **implementation
> validity** (does the code target the intended posterior? — established by
> analytic checks, brute force, and SBC) from **empirical adequacy** (is the
> model right for real data? — probed by posterior-predictive checks and
> documented failure regimes). The two are not the same, and neither is
> "flexibility".
>
> **The flexibility fallacy.** A more flexible model is *not* automatically a
> better one. Extra parameters can fit noise, widen posteriors, and hurt
> out-of-sample prediction without describing the data-generating process any
> better. Clepsydra never selects a model on flexibility, marginal likelihood
> across incompatible encodings, or in-sample fit. Model comparison uses
> **held-out prediction in a shared observable space** (held-out pairwise
> cognacy / held-out cell states); see the T1/T2 cards and
> `tests/simulation/test_posterior_predictive.c`
> (`test_model_comparison_shared_observable`), which shows the relaxed clock
> failing to beat the strict clock out-of-sample on strictly-clocked data.

## Legend

- **Validity** — evidence the implementation is correct (test file · check).
- **Adequacy** — what posterior-predictive / recovery checks say, and what they
  *cannot* say.
- **Failure regimes** — data-generating situations where the model is known to
  be wrong or weakly identified; the diagnostic that catches each.

---

## P0 — classical homogeneous baseline

**Purpose.** A regression baseline only (one rate `r` for all concepts,
`P(match) = exp(-r·d)`), for analytic cross-checks. Not a scientific model of
vocabulary evolution. Legacy clocks are selected only with an explicit
`--legacy-clock` flag, never silently.

- **Assumptions.** One global rate; no hierarchy, contamination, synonyms, or
  rate heterogeneity; two contemporaneous daughters ⇒ `d = 2T`.
- **Validity.** `tests/regression/test_pairwise.c` — agreement with the
  closed-form glottochronological posterior; posterior normalization; quantile
  ordering; simulation round-trip.
- **Adequacy.** Deliberately inadequate for real basic-vocabulary data (rates
  vary ~100-fold across meanings, Pagel 2007). Present for calibration and
  teaching, not inference.
- **Failure regimes.** Any concept-rate heterogeneity (the norm) — use P1/T1.
  Absolute dates are non-identifiable without an anchor (reported, not faked).

---

## P1 — hierarchical pairwise cognacy model

**Purpose.** Two-language divergence dating with concept-specific rates,
integrated deterministically (no MCMC in normal low dimensions).

- **Assumptions.** `log r_{f,c} = μ_c + β_f + u_{f,c}`; concepts conditionally
  independent given parameters; uncertain cognacy via `y_c`; global
  contamination `ε`; exact synonym marginalization over candidate pairs.
- **Priors.** `μ_c ~ N(m_c^ref, s_c^ref²)` from the reference package;
  `β_f ~ N(β̄, σ_f²)`; `ε ~ Beta(1,19)`. Sensitivity presets tight/standard/wide
  must all be run for any published result. The shared family clock offset
  `β_f` is **marginalized** (Gauss-Hermite quadrature over its prior), not fixed
  to a point — so the reported intervals include family-clock uncertainty. With
  `σ_f = 0` this reduces exactly to the point-β posterior.
- **Validity.** `tests/simulation/test_sbc_pairwise.c` — SBC rank uniformity
  for the divergence time; parameter recovery; posterior normalization vs a
  dense grid; the sensitivity presets produce distinct posteriors; the
  lexical-vs-calendar-time flag is set honestly. **The SBC marginalizes `β` over
  its prior — it is not told the true offset.** Because `β` and `log T` are
  perfectly confounded in the pairwise likelihood, this makes the test able to
  fail: the point-β posterior it replaced (which ignored `σ_f`) is badly
  miscalibrated here (χ² > 130 vs a ≈ 40 gate), so the passing rank-uniformity
  is genuine evidence the `β`-marginal is calibrated, not an artifact of
  conditioning on the truth.
- **Adequacy.** Reference priors *update* under the family data rather than
  fixing constants; prior-dominance is reported when few concepts inform the
  fit. Leave-family-out coverage is checked in `tools/validation/`.
- **Failure regimes.** Only `r·d` (lexical time) is identifiable without a
  calibration/anchor — reported as "absolute years not identifiable". Two
  languages carry little information about the rate hierarchy, so the posterior
  is often prior-dominated (flagged). Undetected borrowing inflates apparent
  retention; use the contamination component and loan modes.

---

## Pairwise decay laws (alternative retention curves)

**Purpose.** Pairwise divergence dating under **non-constant** replacement, for
comparison with the constant-rate baseline. Ported from the `distglot` toolkit
and wrapped in a Bayesian posterior (`docs/decay-laws.md`).

- **Assumptions.** Constant rate (`swadesh`), accelerating rate
  (`starostin`, `exp(-2λt²)`), or a two-component stable-core + regular mixture
  with a retention floor (`vasilyev`, `two-component`). Each defines the
  expected shared-cognate similarity `s(T)`.
- **Validity.** `tests/regression/test_decay.c` — curve monotonicity/floors,
  round-trip inversion, posterior recovery and normalization, and that the laws
  give materially different dates on the same similarity.
- **Adequacy / the flexibility fallacy.** These are competing hypotheses, not
  improvements. Which fits basic-vocabulary decay is empirical; choose by
  held-out prediction in a shared observable space, never by preference. The
  tool ports them to make the comparison possible.
- **Failure regimes.** Absolute dates rest on the law *and* its parameters
  (`r`, `λ`, the two-component constants); a two-component floor makes very deep
  pairwise dates weakly identified (similarity near the floor → time barely
  constrained). Reported, not hidden.

## T1 — strict novel-state multistate model

**Purpose.** The primary tree model: per meaning, a cognate class is replaced
by a never-before-seen class along the tree (Rönchen et al. 2024, independently
re-derived). One primary cognate per language–meaning.

- **Assumptions (of T1 only, not the project).** One primary state per cell;
  every replacement opens a globally-new class; no back mutation, no parallel
  innovation, no reappearance of extinct classes; meanings evolve
  independently; strictly tree-like inheritance (no undetected borrowing);
  unobserved ancestral states allowed.
- **Parameters/priors.** Concept rates via the §4 hierarchy; node ages from
  calibrations; relaxed branch multipliers `log m_b ~ N(-σ_b²/2, σ_b²)` with a
  HalfNormal shrinkage prior on `σ_b`; contamination `ε ~ Beta(1,19)`.
- **Validity.**
  - *Likelihood:* `tests/likelihood/test_t1_likelihood.c` — agreement with
    **two independent brute forces** (edge-pattern enumeration and
    set-partition enumeration, neither reusing the `p_v/q_v` recursion) to
    `1e-9` on binary trees, polytomies, missing/observed-internal patterns, and
    60 randomized cases; the tree likelihood and the contamination mixture each
    normalize to 1 over the observable partitions; simulation frequencies match
    the likelihood; the relaxed (per-branch) likelihood is brute-force checked
    too.
  - *Inference:* `tests/simulation/test_mcmc_tree.c` — strict-clock-tree SBC
    rank uniformity, seeded reproducibility, warm-up-only adaptation,
    checkpoint identity; `test_mcmc_relaxed.c` — branch-rate recovery where
    identifiable, `σ_b` shrinkage, and contamination robustness.
- **Adequacy.** Posterior-predictive checks
  (`tests/simulation/test_posterior_predictive.c`): a well-specified fit is
  calibrated (moderate PP p-value); a homogeneous-rate fit to strongly
  concept-rate-heterogeneous data is **flagged inadequate** (PP p-value ≈ 0)
  even though it is validly implemented. This is the canonical
  validity-≠-adequacy example.
- **Failure regimes.**
  - *Concept-rate heterogeneity ignored* — detected by the class-count-variance
    PP check; mitigate with per-concept rates (`μ_c` block).
  - *Branch-rate / node-age confounding* — only `r·m_b·t_b` is identifiable per
    branch; `gd_mcmc_run` reports the branch-multiplier / branch-duration
    posterior correlation and flags a strong trade-off. Requires calibrations
    and `σ_b` shrinkage; warn when weakly identified.
  - *Absolute-scale non-identifiability* — without an anchor only lexical time
    is identifiable.
  - *Borrowing / miscoding* — the contamination mixture absorbs incompatible
    concepts and reports `P(concept contaminated)`; incompatible concepts are
    never silently removed.
  - *Polytomies* — resolved to zero-length ladders (exact for this model,
    brute-force verified); the exact non-binary recursion awaits the appendix
    supplement, so **no published polytomy result before it is checked**.
  - *Synonym coexistence* — not modelled (single primary state); use T2 for
    polymorphism.

---

## T1 relaxed clock (extension of T1)

**Purpose.** T1 with per-branch rate multipliers for lineage-specific rate
variation.

- **Assumptions.** As T1, plus `log m_b ~ N(-σ_b²/2, σ_b²)` (so `E[m_b]=1`)
  with shrinkage on `σ_b`; the strict clock is `σ_b = 0`.
- **Validity.** Relaxed likelihood brute-force checked
  (`test_t1_relaxed_matches_bruteforce`); branch-rate recovery where
  identifiable (`test_mcmc_relaxed.c`).
- **Adequacy / the flexibility fallacy.** On strictly-clocked data the relaxed
  clock does **not** improve held-out prediction over the strict clock
  (`test_model_comparison_shared_observable`); the shrinkage prior keeps `σ_b`
  near zero (`test_m9_sigma_b_shrinks_to_strict`). More parameters are not
  automatically a better model. **But the verdict is data-dependent, and the
  discipline is to check.** On *real* data across three families
  (leave-one-calibration-out, `examples/loco_benchmark/`: Indo-European,
  Sino-Tibetan, Romance) the relaxed clock recovers held-out attested tip ages
  far better than the strict clock (90% HPD coverage **0.89 vs 0.33** over nine
  held-out points, with lower RMSE) — because real lineages are rate-
  heterogeneous, so the flexibility earns its keep exactly where the
  strictly-clocked simulation says it should not. Same test, opposite verdict,
  decided by held-out prediction on the data at hand — never by preference for or
  against flexibility.
- **Failure regimes.** Branch rates and node ages are strongly confounded
  (above); without calibrations the multipliers absorb age uncertainty. The
  confounding flag and calibration leave-one-out are required before trusting
  branch-rate estimates.

---

## T2 — binary cognate-set model

**Purpose.** Each cognate class is a binary present/absent trait under an
asymmetric two-state CTMC; a polymorphism-tolerant alternative to T1's
single-primary-state assumption.

- **Assumptions.** Per-class gain/loss CTMC `0 --λ_g--> 1 --λ_l--> 0`; classes
  within a meaning treated as independent binary characters (they are **not**
  truly independent — flagged); all-absent classes unobservable.
- **Key design choices (King 2026).** Ascertainment correction conditions each
  class on being observed somewhere; **partition rates are weighted by the
  number of meanings, not cognate classes** (the cognate-count weighting that
  miscalibrates under ascertainment is opt-in only).
- **Validity.** `tests/likelihood/test_t2_likelihood.c` — Felsenstein pruning
  vs an independent full-enumeration brute force (`2^n` node states) on binary
  trees, polytomies (handled natively — no ladder), missing/observed-internal
  patterns, and randomized cases; the ascertainment-corrected likelihood
  normalizes over the observable patterns; simulation matches; rate recovery.
- **Adequacy.** Within-meaning independence is false; the dataset result raises
  `within_meaning_dependence_warning` whenever a meaning has more than one
  class, so the assumption is never relied on silently.
- **Failure regimes.**
  - *Within-meaning dependence* — meanings with many classes contribute
    disproportionate likelihood; warned.
  - *Cross-encoding comparison* — raw T1-vs-T2 marginal likelihoods are **not**
    comparable; use held-out prediction in a shared observable space only.
  - *Heavy polymorphism* — the ascertainment × partition-rate interaction
    (King 2026) is a known adequacy risk; posterior-predictive checks on
    present-cognate counts and singletons are required.

---

## T3 — finite-state Mk/CTMC (reserved)

Not implemented. Reserved for a finite-state multistate CTMC permitting
homoplasy (relaxing T1's no-parallel-development assumption). Its role would be
robustness comparison, **not** "more flexible = more realistic". It will not be
added before T1 and T2 are validated on real data, and only if held-out
prediction justifies it.

---

## How the cards are kept honest

- Every "Validity" claim names a test that runs in `make test` (currently
  190/190, sanitizer-clean).
- Every "Adequacy" claim is a posterior-predictive or recovery check, with its
  limits stated; passing SBC is never presented as adequacy.
- Every "Failure regime" pairs a data situation with the diagnostic that
  detects it. Where no diagnostic yet exists (e.g. the exact non-binary T1
  recursion), the card says so and the result is withheld.

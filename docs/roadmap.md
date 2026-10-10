# Roadmap: from 0.2 to a reference tool

> Drafted 2026-10-10, revised the same day after two adversarial reviews of
> the first draft (`docs/research/reviews/2026-10-10-roadmap-codex.md`,
> `…-roadmap-claude.md`). Supersedes the open parts of `PLAN.md` (M0–M12 are
> done; this is what comes after). Evidence base, all in `docs/research/`:
> the Arawan field report (`arawan-lessons.md`), two adversarial code reviews
> (`reviews/2026-10-10-codex.md`, `reviews/2026-10-10-claude.md`), the
> lexical-dating literature (`literature-lexical-dating.md`), dating methods
> from other fields (`literature-dating-methods.md`), a mathematical analysis
> with numerical experiments (`mathematical-analysis.md`, `math/`), and a
> 367-row calibration database (`reference/calibrations/`).

## The goal, stated so it can be tested

Given a cognate-coded wordlist for a family, with or without a topology and
with or without calibrations, return for every node (or every named set of
languages) a posterior age distribution. Its intervals are checked out of
sample, on real families whose dates are known, by a pre-registered test
that is allowed to fail.

Four things follow from the research and shape everything below.

1. **Three kinds of age, always reported separately.**
   - *Lexical* age (`ρ·a`, in expected replacements).
   - *Relative* age (`a_v / a_root`).
   - *Calendar* age (needs an external scale).

   Under a strict clock the first two are identified by the cognates
   alone. In the ideal case (200 meanings, known concept rates, correct
   model) lexical age is good to about ±8%. Under a relaxed clock only
   branch lexical lengths are directly identified, so node lexical and
   relative ages depend partly on the prior. The Arawan failure came from
   presenting calendar age as if it were lexical.
2. **For an uncalibrated family, the calendar interval is set by the
   cross-family clock, not by the cognates.**
   - The between-family rate SD σ_f is a floor that no lexical data removes
     (`mathematical-analysis.md` §§3–4).
   - The value 0.3–0.4 (log units) comes from in-sample or non-T1
     estimates. It has not been measured under this model.
   - With σ_f estimated from about 8 families, the predictive is a t₇, not
     a normal. A realistic expectation for a 2,000 BP median is ×/÷ 2.3–3
     (about 700–5,500 BP), not ×/÷ 1.9.
   - The tool promises a method (an interval whose sources are shown and
     whose transport has been tested), not a width.
3. **The calendar promise may turn out to be untenable.** That has to be
   learned early and cheaply (D2, the σ_f pilot), not after three phases of
   engineering. If uncalibrated calendar dating does not beat trivial
   baselines out of sample, the product is lexical and relative chronology
   plus explicit calendar *scenarios*. That outcome can be published.
4. **The current engine has no proper prior on node ages.**
   - Every uncalibrated date it has produced is shaped by a hidden `1e15`
     bound and an `age^(m−1)` volume tilt.
   - The contamination floor makes the posterior improper.
   - Nothing on the engine track is worth doing before this is fixed.

## Status (2026-10-10)

- **E0: done.** All items are fixed and tested. E0.13: no external T1 code
  exists, so the in-repo brute forces stand.
- **E1: done, except the pipeline-level SBC.**
  - E1.1 and E1.2: the lexical-root prior, with calendar calibrations as
    factors.
  - E1.3: uniform prior plus Yule.
  - E1.4: the `1e15` cap removed.
  - E1.5: `--prior-only`.
  - E1.6: three ages, scale source, prior beside posterior, width
    decomposition, power-scaling.
  - E1.7: window-relative moves and the clock-preserving move.

  The gate's SBC runs at the sampler level with the real T1 likelihood
  (`test_mcmc_sbc_proper_prior_scenarios`). It covers a free root, normal and
  offset-lognormal root calibrations, a dated tip, and the Yule prior, with
  ranks per parameter plus the log-likelihood, and all pass. SBC through
  `gd_analysis` (input parsing included) moves to D3, which needs the same
  forms-writing simulator.
- **Arawan at gate E1** (`docs/research/arawan-acceptance.md`):
  - the root is finite and converged without a window;
  - with a window it sits inside it rather than on its top;
  - the lexical root is 0.54 ± 0.05 either way;
  - the clock sets ~100% of the calendar width, and the report says so.

- **M1/M2 at pilot scale: done, and a no-go again**
  (`docs/research/clock-engine.md`).
  - The full engine was run on Arawan-sized Glottolog subtrees around each
    of the 15 two-sided nodes, with no calendar information.
  - It gives σ_tot = 0.70 (95% 0.40–1.43) against the pilot's 0.77. The
    engine's node rates agree with the pilot's (r = 0.81). The crude
    pairwise estimator was not what made the clock spread large.
  - The 95% width at 2 ka is ×/÷ 5.8.
  - Leave one family out, the clock beats the trivial baselines on point
    estimate, but every skill interval includes zero.
- **E2: done apart from its gate.**
  - Done:
    - E2.1: concept rates, with an empirical prior on their spread
      (`concept-rate-prior.md`).
    - E2.2: sampled ε and the concept table.
    - E2.3 and E2.5: synonyms and cognacy confidence in `date-pair`.
    - E2.4: codings and loans.
    - E2.6: outlier calibrations.
    - E2.7: `ancestor_tip`.
    - E2.8: gamma mixing. Lomax and strict/relaxed spike-and-slab remain open.
    - E2.9: `tip_error`, with the empirical prior from D1.3.
    - E2.10: the adequacy panel.
  - The gate's simulator (D3, `tools/simworld`) is running.
- **E3:** done:
  - E3.0: the subtree move.
  - E3.1: incremental likelihood. Edge terms are computed once per mixture
    node, and only dirty nodes are recomputed; this is exact, and 1.65×
    on Arawan.
  - E3.2: draws and pointwise output.
  - E3.4: `date-trees`.
  - Added: parallel chains (`--threads`; bit-identical draws, 2.9× on four
    chains).
  - Added: a non-centred σ_b move for the relaxed clock. Without it, the
    Arawan relaxed fit stalled at R-hat 1.16.

  - E3.3: checkpoints bound to an input hash, and carrying their draws.

  E3 is complete.
- **D1:** done: D1.2 (`reference/calibrations/curated_nodes.tsv`) and D1.3
  (`error-prior.md`).
- **D2, the σ_f pilot: done, and a no-go on its pre-registered rule**
  (`docs/research/clock-pilot.md`).
  - Over 15 nodes in 6 families, σ_tot = 0.77 (95% 0.50–1.40). Most of it
    is within-family node scatter (τ = 0.64), not σ_f (0.33).
  - Out of family, the learned clock does not beat the trivial baselines
    reliably, and is slightly worse than the old default.
  - An uncalibrated 2 ka node gets ×/÷ 6 at 95%. The expected product for
    uncalibrated families is lexical and relative chronology plus labelled
    calendar scenarios.

## Structure: two tracks that meet

The first draft gated every phase on the previous one. Both reviews showed
that this delayed the one experiment that decides the project, the
cross-family clock, behind sampler performance work. The revision has two
tracks:

- **Engine track (E0→E3).** Strictly gated: correctness, a proper target,
  and a frozen observation model.
- **Data and validation track (D0→D5).** Starts on day one. It needs no C
  changes until the tracks merge.
- **Merge (M1→M4).** The clock package, its frozen validation, the Arawan
  result, and the release.
- **Optional expansion (X).** Chosen after M2, on demonstrated need.

**Arawan is a continuous acceptance case**, not a final example. From E0 on,
every gate re-runs KAHD and reports two things:

- whether the root is still glued to anything;
- the Paumari 1833 versus modern Paumari separation, a **negative control**
  that is never used to fit anything.

Calendar dates for Arawan are not published until M2 passes.

Evidence tags: `C#`/`K#` are first-review findings (codex/claude), `RC`/`RK`
are roadmap-review findings, `A#` is an Arawan item, and `§` is a section of
the mathematical analysis.

---

## Engine track

### E0 — Correctness, and an external check on the likelihood

| id | fix | evidence |
|---|---|---|
| E0.1 | Missing cognate IDs (`.`, empty, `?`) and `status ∈ {missing, excluded, absent}` are not observations; count and report skipped rows. *(Verified: today two `.` cells compare equal and count as a cognate match.)* | K-C3, C12 |
| E0.2 | Extant tips default to age 0. A tip with a non-zero or uncertain age needs an explicit tip calibration; an undated tip in a run with ancient tips is an error | K-C4, C4 |
| E0.3 | `normal_node` support uses `mean ± k·sd`, not `exp(mean + 5sd)`; lognormal tails are not silently truncated | K-C5, C5 |
| E0.4 | `--cal-widen` widens the dispersion parameter of every calibration family (log-SD for offset-lognormal), with a test per type | K-M13, C6 |
| E0.5 | `is_hard` is honoured (soft bounds get stated tail mass) or removed from the schema | K-C2, C5 |
| E0.6 | `date-pair` takes tip ages (`--age-a/--age-b`) and flags ancestor–descendant pairs | A7, K-M9 |
| E0.7 | One rate unit (per year, log scale) across all commands; `--rate-unit` for legacy input | A8 |
| E0.8 | Unknown language IDs, zero shared concepts, or a tip with no forms are errors, not prior-only "results" | C20 |
| E0.9 | Convergence gate: folded + rank R-hat, bulk and tail ESS, MCSE (Vehtari et al. 2021, golden-tested against ArviZ vectors); non-zero exit below thresholds unless `--allow-unconverged` | A4, K-M10, C17 |
| E0.10 | Language count from `is_leaf`, not `(n+1)/2`; forms languages resolve to leaves only | C21, C22 |
| E0.11 | Calibrations may target a taxon set (`mrca:lang1,lang2,...`); unlabelled internal nodes get stable generated names | A6 |
| E0.12 | Never mutate `feasible`; the initial state must have finite prior and likelihood or the run fails | K-m1/m2/m3, C23 |
| E0.13 | **Likelihood oracle.** *Status 2026-10-10:* Rönchen et al. published no code; their LDC 2024 appendix (figshare 26089507) holds proofs only, with no worked numbers. The T1 recursion is therefore checked against two in-repo brute forces that do not share it (`tests/likelihood/test_t1_likelihood.c`: enumeration of edge-replacement patterns, and of all node labellings, including missing tips, polytomies and dated tips). An external cross-check moves to X5/M4 (BEAST2/RevBayes on shared fixed-parameter cases). | K §4.5, RC, RK |
| E0.14 | CLI rejects options missing their value; short reads checked; public API null checks; RNG bound arithmetic. *(Needed, but does not block experiments.)* | C27, C28 |
| E0.15 | T2 marked **experimental** in CLI help, docs and README; P1's fixed ε documented. Docs claims match code | C10, C11, C32, K-M6 |

**Gate E0:** regression test per fix; `make test` and `make analyze` clean;
examples re-run with changed outputs explained in the changelog; E0.13
passes; Arawan re-run recorded.

### E1 — A proper target, sampled well, reported honestly

The central change is **two explicit modes**. The first draft's "proper
calendar root prior" and its "free re-dating by convolution" contradict each
other (RC-1, RK-2.1). A proper calendar root prior is not scale-invariant.
The scale-invariant prior `1/a_root` is improper under the contamination
floor.

| id | item | evidence |
|---|---|---|
| E1.1 | **Lexical-time mode (default when there is no calendar calibration).** A proper, broad prior on the *lexical* root age `ρ·a_root`; a prior on the dimensionless ratios; ρ independent with a clock prior. Calendar ages are then *defined* as the lexical posterior convolved with the clock prior. This is exact and proper, so one run can be re-dated under any clock package. Valid with a relaxed clock only if branch multipliers have mean 1, and with concept rates only if they are centred (E2.1) | §4.1, RC-1, RK-2.1 |
| E1.2 | **Calendar time-tree mode** (any calendar calibration, dated tips, or a calendar-unit tree prior): a joint proper prior on calendar ages and ρ, refit whenever a calibration or prior changes; no convolution shortcut. A user's calendar "root window" (e.g. Arawan's U(300, 3000)) is **a calibration** with a justification field and an outlier option (E2.6), never "the root prior". The engine reports which mode ran and why | RC-1, RK-2.1, RK-8 |
| E1.3 | **Ratio priors**: `uniform_ratios` (uniform node ages given the root; the MrBayes uniform-clock construction, *not* MCMCTree's, which is a conditional birth–death kernel) and **birth–death with sampling fraction ≈ 1** as a sensitivity pair. Lexical samples are usually *diversified*, one doculect per subgroup, so the sampling scheme is stated, not assumed. With dated tips the ordered-region volume is piecewise polynomial in the root. Its exact normaliser is computed by a post-order recursion, so the stated root prior is what you get | RC §2 R1.1, RK-2.2, lit-dating §B |
| E1.4 | `1e15` removed. **Propriety checks** refuse improper configurations. Regression test: enlarging any numerical limit leaves the posterior unchanged | C1, C2, K-C1 |
| E1.5 | **Exact prior-only mode** (`--prior-only`: same sampler, likelihood off), replacing `gd_cal_sample_prior` | C3 |
| E1.6 | **Reporting.** Every node gets lexical, relative and calendar age. Calendar ages carry `scale_source ∈ {calibration, tip, clock_package, user_rate_prior}`. The report adds: effective prior beside the posterior (from E1.5); **power-scaling prior sensitivity** (Kallioinen et al. 2024) per component; a **variance decomposition** of each calendar interval into data, branch-rate, clock and model parts; a flag when an HPD endpoint is within 1% of a bound; and a "posterior ≈ prior" flag. These replace `prior_dominance`. In calendar mode, *calibrated* marginals generally differ from what the user wrote (Heled & Drummond 2012), and the report shows this | A2, A10, K-C2, RK-6.1–2, RC-7.5 |
| E1.7 | **Sampler geometry** *(moved from R3, needed to pass this gate)*: ratio + log-root parameterisation with Jacobian; a clock-preserving move; step sizes from the current state with dual-averaging adaptation in warm-up | §5h, K-M11, RK-1.2 |

**Gate E1:** SBC through the full `gd_analysis` pipeline, per parameter, with
the joint log-likelihood as an extra test quantity (Modrák et al. 2025), so
that "posterior = prior" is detectable. Cases: (a) lexical-time mode with a
clock prior; (b) normal and lognormal root calibrations; (c) ancient tips;
(d) both ratio priors. On Arawan, the root no longer sits on any bound, and
the report says the scale is clock-transferred.

### E2 — The frozen observation model

This is the **minimal** model that the clock package will be fitted with
(M1). Everything else is a sensitivity variant with its own package, never a
silent change of default. Every default change after M1 invalidates the
package (RK-1.2).

| id | item | evidence |
|---|---|---|
| E2.1 | **Concept-rate heterogeneity, on by default and centred.** Mean 1 (or geometric mean 1) over a **canonical Concepticon reference scale**, not over each family's own list, so that β_f is transportable and the mean-rate uncertainty sits entirely in ρ (dos Reis et al. 2014). Integrated out per concept with **tail-faithful quadrature**: K ≥ 8 mean-matched categories or generalised Gauss–Laguerre/Hermite, checked against high-K on the deepest training node, because deep ages are dated by the slow tail. Gamma and lognormal both implemented; the choice is a model-uncertainty axis (E2.8) | C8, K-M2, §2, RC-R2.1, RK-2.3 |
| E2.2 | **Contamination sampled** (ε ~ Beta, prior from the package, since 7 languages barely identify it), with per-concept responsibilities in `posterior_concepts.tsv`. The fixed CRP concentration is revisited. P1 gets a real ε and chance-match probability | K-M3, C9, C14, RC-R2.6, RK-8 |
| E2.3 | **Synonyms, first as sensitivity.** Exact marginalisation in pairwise mode. In T1, alternative resolutions (one synonym or all; documented selection rules) run as a coding-sensitivity axis; no latent-primary Gibbs sampler until the latent variable has defined generative semantics | A11, K-M5, C12, RC-R2.3 |
| E2.4 | **Loans and coding columns read** (`loan_probability`, `inheritance_probability`, `cognacy_confidence`, `analysis_id`). Modes: `exclude` (external loans become missing for inheritance) and `sensitivity`. A simulation test with *compatible* sister loans shows the young bias. A general `mixture` is deferred | K-M4, C14, A14, RC-R2.4 |
| E2.5 | **Uncertain-cognacy semantics unified** (one meaning of `y_c`, truth-table tests) | C15 |
| E2.6 | **Calibrations that can be wrong** *(moved from R6, needed by the training data)*: heavy-tailed or outlier-mixture calibration densities; **one-sided** "no later than" and attestation-lag calibrations with stated tail mass; **censored** likelihoods for min/max-only rows; calibration conflict and leave-one-out diagnostics | lit-dating rec. 3, RC-2, RK-1.2 |
| E2.7 | **Ancestor-like tips (minimal).** Allow a zero-length ancestor branch for a dated tip (Latin, Old Norse, Paumari 1833), or exclude such tips from training. Full FBD waits (X3) | K-m16, Chang et al. 2015, RK-1.2 |
| E2.8 | **Sensitivity variants, not defaults:** branch × concept Lomax frailty (not consistent under branch subdivision, so it depends on node density, a transport hazard); source-level coding error (E2.9); mixing family (gamma vs lognormal vs a two-component slow tail); strict vs relaxed clock (a spike-and-slab on σ_b, so the clock is not the user's free choice). Deep ages are reported **model-averaged**, with a printed model-uncertainty term | §2.4–2.5, §11, RC-R2.2, RK-2.4, RK-6.3, RK-6.9 |
| E2.9 | **Source/list-level observation error as sensitivity.** The parameter is attached to a *source* (collector, list, coder), not to a language. Under T1 it cannot be distinguished from terminal-branch time, so it is identified only by its prior. The prior comes from same-language, different-collector doculect pairs in Lexibank (D1.3) and **never from Paumari 1833**, which stays a test | A9, §8, RC-R2.5, RK-2.5 |
| E2.10 | **Adequacy panel** in `date-tree` output: PPCs on class counts (exact `E[K]`), singleton and doubleton frequencies, pairwise sharing residuals, triplet identities; **clock adequacy** (unconstrained branch lexical lengths vs clock-predicted); **bootstrap over meanings** vs the model SD; a **chance-resemblance floor** flag for deep nodes; a **temporal-signal test** (root-to-tip and date randomisation) when tips are dated | §1.4, §12, RC-7.4/7.7/7.9, RK-6.12 |

**Gate E2:** PPC panel runs. The misspecification simulator (D3) reports
bias and coverage for borrowing, split/lump error, Γ heterogeneity,
correlated concept domains, non-random missingness and compatible loans. The
default model is **frozen and versioned**.

### E3 — Throughput (as needed, non-gating)

| id | item | evidence |
|---|---|---|
| E3.0 | **Subtree-scaling move**: nested deep nodes with close ages (Sino-Tibetan inner_TB / eastern_TB) mix slowly under single-site moves; the example needs 12,000 samples per chain to pass the convergence check with concept rates on | E2.1 runs |
| E3.1 | Likelihood caching: data extension, sort keys and the noise term computed once per concept; age moves recompute only the path to the root | K-M11 |
| E3.2 | Full posterior draws, pointwise log-likelihood, rates and hyperparameters, every node (fixed ones too) in the output. *(Needed by M2 scoring, so do it before M1.)* | C18 |
| E3.3 | Checkpoints bound to input hashes, resumable mid-sampling | C25 |
| E3.4 | **Tree-set input** (R5.1 moved forward): date over a supplied posterior tree sample (yvyra, MrBayes, BEAST) and pool; per-clade ages conditional on monophyly. **MRCA-of-named-set ages** with `P(monophyletic)` as the primary quantity when topology varies. When the trees were estimated from the *same* cognates, this is a hybrid approximation: say so, and compare with a split-by-concepts run | C19, K-C6, §5j, RC-R5.1 |

There is no ESS/s gate. Profiling decides when E3.1 is worth doing.

---

## Data and validation track (starts now)

### D0 — Pre-registration

Written and committed **before any real LOFO number is computed**:

- the family list and clustering rules (nested datasets and shared sources
  stay together; leave-one-team-out as a robustness check);
- grade rules, and deduplication by **node identity** (MRCA glottocode set),
  not by row id;
- exclusions:
  - creoles;
  - Bayesian/lexical "calibrations", which are circular (e.g. Pama-Nyungan
    root, Proto-Uto-Aztecan root → grade C);
  - published ASJP dates;
- held-out-family rules: all of that family's node calibrations removed;
  whether dated tips are kept is declared up front, and the two strata are
  reported separately;
- metrics: log score of the calibration *density* under the predictive;
  interval score; PIT against calibration draws; censored scores for
  one-sided rows; Σz² against χ²_F; PIT ECDF with simultaneous bands; depth
  slope; multiplicative width;
- **baselines**: a clock-prior-only predictor; language count / family size
  (Rama & Wichmann 2020); ASJP/Holman; an empirical depth-rank
  ("climatological") predictor; the 0.2.0 default rate;
- **falsifiable pass criteria** (§ M2).

### D1 — Curation

| id | item | evidence |
|---|---|---|
| D1.1 | Lexibank ↔ Glottolog ↔ Concepticon ↔ calibration-node mapping for the 59 local expert-cognate datasets; Concepticon overlap graph; coding protocol, team and source metadata per dataset | §4.2, RC-5, RK-2.7 |
| D1.2 | Calibration database clean-up: merge duplicate nodes (Ma'anyan–Malagasy, Cham, Roglai–Tsat, Malayo-Polynesian, Proto-Oceanic, Proto-Central Pacific, Inuit, Cholan, Afrikaans–Dutch); demote circular rows; mark max-leaning "points" as one-sided; add a **justification schema** (Parham et al. 2012) and uncertain-target-node mixtures for archaeological associations; check whether the Arawakan ceramic correlations were themselves built on lexicostatistics | RK-3.2, RC-3, RK-6.8 |
| D1.3 | **Same-language, different-collector pairs** from Lexibank and KAHD (time ≈ 0), giving the empirical prior for source-level error (E2.9) | RK-2.5 |
| D1.4 | **Realistic F.** The A/B internal-node rows cover 19 nominal families, and after nesting, creoles and circular rows are removed roughly 6–9 independent families remain, several with a single shallow node. Choose among: C-grade rows under the outlier model; tip-calibrated families (e.g. Tupí-Guaraní with Old Tupinambá/Old Guaraní); or accepting F ≈ 7 and the wider t-predictive | RC-3, RK-3.2 |

### D2 — σ_f pilot (go/no-go for the calendar promise)

Use fast estimators: pairwise ML distances with Godambe correction (§5a–b,
Python), and current T1 with the E0 fixes on fixed published trees with a
crude root prior. Fit `b_f ~ N(μ, σ_f² + s_f²)` and run a crude LOFO on the
calibrated families. Also fit σ_f under **one uniform automatic coding**
(LexStat on the Concepticon core, KAHD included). If automatic coding gives
a smaller σ_f, the clock should transport through it. This takes days, not
months, and it sets expectations for everything else.

### D3 — Misspecification world simulator

Built early (R7.2 moved forward). It generates borrowing networks, coding
error, time-dependent rates and split pulses, at target-like sizes (6–9
languages). It calibrates the **power and decision rules** of M2 before real
numbers are seen, and it supplies the E2 gate. Power check at α = 0.05 for
detecting 95% intervals that really cover 80%:

| test | F = 8 | F = 10 |
|---|---|---|
| binomial coverage | ≈ 0.20 | ≈ 0.32 |
| Σz² | ≈ 0.58 | ≈ 0.65 |

### D4 — Arawan data work

| id | item | evidence |
|---|---|---|
| D4.1 | Add the **Natterer 1833** Paumari source to the calibration database (currently only Ehrenreich 1889) with collection and publication dates, lect identity, uncertainty; likewise Chandless 1869 for Arua | RC-4 |
| D4.2 | *(Out of scope, user decision 2026-10-10.)* **Arawan soft calibrations**, curated and graded with experts (Dienst, Dixon, Ferraz Gerardi). Candidates: Tupí-Guaraní / Língua Geral and Portuguese loan strata as upper bounds on late splits; Proto-Arawan ↔ Purus-Arawakan contact loans (`carvalhopurus` is local); *Wörter und Sachen* on cultigens; middle-Purus archaeology as a soft maximum. This is **the only lever that beats the σ_f floor**: even a one-sided soft maximum cuts the long upper tail | RK-4.2 |
| D4.3 | Recode the two historical lists independently (ideally two coders); strict/loose cognacy and alternative resolutions of the 208 polymorphic cells as coding variants; external loans marked | RC-4, A11 |
| D4.4 | *(Out of scope, user decision 2026-10-10.)* **Protocol effect for KAHD.** Double-code about 50 concepts under the IE-CoR protocol, which turns an unknown transport shift into a measured one | RK-2.7, RK-4.2 |

### D5 — Target-like evaluation design

For each held-out family, also evaluate on random **6–9-language subtrees**
that contain the calibrated nodes, with topology from E3.4. That is the
Arawan regime: sparse sampling, long branches, few depths. Report an
**Americas / small-oral-family stratum** separately, even if it is
anecdotal. Out-of-region transport is the assumption the Arawan number
rests on.

---

## Merge

### M1 — The clock package

| id | item | evidence |
|---|---|---|
| M1.1 | Fit the frozen E2 model per training family. Either a **joint multi-family model**, or two stages that carry each family's **likelihood**: draws, or a likelihood grid in β for one-sided cases, prior-corrected. Not posterior means treated as data | RC-R4.2, RK-2.7 |
| M1.2 | Population model: `β_f ~ N(μ + depth term, σ_f²)` with family, protocol and team dispersion (team as a random effect where crossed data exist; otherwise said to be folded into σ_f). **Intercept-only first**, at most one pre-registered shrinkage covariate (population size); covariate selection nested inside LOFO. Concept effects on the canonical scale, with Concepticon covariates for concepts unseen in training | RC-R4.3, RK-2.7, lit-dating rec. 4 |
| M1.3 | A **versioned package the C engine loads**: predictive *draws* (t-like, including σ_f uncertainty), concept effects, depth term, domain manifest (concept list, protocol, model version, training families, D0 plan hash). The engine evaluates the mean concept effect **on the target list** and adds the variance from unobserved concepts to the decomposition (E1.6). Synthetic sources removed | C7, K-M1, RK-8 |

### M2 — Frozen validation (the headline result)

Nested LOFO following D0 and D5, with every learned quantity rebuilt inside
each fold: hyperparameters, concept effects, mixing shape, ε prior,
covariates, calibration supports, thresholds. **Pass criteria, all
pre-registered:**

1. **Skill.** Beats the language-count and climatological baselines in mean
   log score, with a family-clustered bootstrap CI excluding 0.
2. **Sharpness.** The 95% interval for a 2 ka node is no wider than ×/÷ 3.5
   (or the pre-registered value).
3. **Dispersion.** Σz² and the PIT ECDF are reported with the posterior of
   the inflation factor.
4. **No depth slope** in the target depth band.

**On failure:** the result is published as is. For uncalibrated families the
tool then reports lexical and relative ages plus labelled calendar
scenarios. Any recalibrated successor is a new package version, validated
on untouched families or by nested recalibration, never on the folds that
estimated it.

### M3 — Arawan

KAHD, with every layer reported:

- branch lexical lengths;
- strict-clock relative chronology, with relaxed-clock sensitivity;
- topology (E3.4), coding (D4.3) and loan sensitivity;
- runs with and without the historical lists;
- calendar ages under the M1 package (plus D4.2 soft calibrations, should
  anyone curate them later);
- the variance decomposition;
- Paumari 1833 as the held-out check.

The grade C 3,200 BP archaeological maximum appears only as a sensitivity
run. The ASJP and Swadesh values appear as benchmarks only.

### M4 — Documentation and release

| id | item |
|---|---|
| M4.1 | Examples regenerated with prior-only overlays and the new report; LOCO scored by interval score |
| M4.2 | `docs/scientific-methods.md` rewritten around the two modes, the three kinds of age, the scale-source flag and the variance decomposition |
| M4.3 | Model-uncertainty table of sensitivity variants (E2.8, ratio priors, mixing families) |
| M4.4 | Citable release (Zenodo DOI, CITATION.cff). Methods paper whose core claim is the **pre-registered M2 table, including a failure outcome if that is what happens** |

---

## Optional expansion (after M2, on demonstrated need)

| id | item | note |
|---|---|---|
| X1 | **Native topology MCMC** (exchange, subtree-slide, Wilson–Balding), explicit topology prior, CCD summaries (Berling et al. 2025) | Was R5.2. Arawan has BEAST/yvyra tree samples, so E3.4 covers it first |
| X3 | Full FBD / sampled ancestors | E2.7 covers the minimal case |
| X4 | Pairwise composite engine in C | The Python version is built for D2 |
| X5 | Completing T2 (two rates, meaning weights, stable ascertainment, BEAST2/RevBayes cross-check) | Experimental since E0.15 |
| X6 | Decay laws with priors on their constants | Until then, labelled as fixed-constant illustrations |
| X7 | Branch × concept heterogeneity in a subdivision-consistent form (Markov-modulated concept rates) | Replaces the Lomax variant if it earns its keep |

## What is deliberately not on this roadmap

- **Neural SBI or variational inference as primary engines.** VI
  underestimates exactly the variance we must report, and neural SBI is
  unfaithful under misspecification. Either may later serve as a validator
  (§5c–5d).
- **Full reticulate (network) inference.** The PPCs come first (E2.10).
- **Automatic cognate detection in the engine.** LexStat/cognator output
  enters as an alternative coding (E2.3–2.4) and as the D2 uniform-protocol
  experiment.

## What will *not* narrow the Arawan root interval

These may be worth doing for other reasons, but they do not narrow the root
interval much:

- the Lomax branch × concept variant (about 5% at 2 ka);
- source-level error (about +80–180 lexical years);
- topology machinery (the root is the MRCA of everything);
- sampler speed;
- more concepts (§3.3);
- more than about 10 training families (σ̂ 0.315 → 0.305).

What does narrow it, in decreasing value per effort:

1. Arawan soft calibrations (D4.2).
2. Depth-matched, target-like training (D5).
3. A measured protocol effect (D4.4).
4. A covariate that survives nested LOFO.

## Decisions taken with the user (2026-10-10)

- Roadmap first, reviewed, then implementation with tests. Every finished
  step is committed and pushed.
- Topology: supplied tree sets (E3.4) and native topology MCMC (X1).
  **Exhaustive enumeration is dropped** (X2 removed), after both roadmap
  reviews.
- **Default ratio prior: `uniform_ratios`.** It is the λ→0 (critical)
  limit of the birth–death prior conditioned on the root, so it is a
  coherent generative prior on ranked histories and remains valid if
  topology is later sampled. It has no rate parameter, so lexical-time mode
  stays exactly scale-invariant. It is also the transparent baseline the
  reviews asked for. Birth–death with sampling fraction ≈ 1 is the shipped
  sensitivity alternative and appears in the model-uncertainty table.
- **Breaking changes are accepted in 0.x:** defaults, CLI, rate units,
  output columns and example numbers may all change. They are recorded in
  the changelog.
- **Expert-dependent work is out of scope:** D4.2 (Arawan soft
  calibrations with specialists) and D4.4 (protocol double coding). D4.1
  and D4.3 remain only as far as they can be done from published sources.
  Without D4.2, the Arawan calendar interval is expected to sit near the
  σ_f floor.
- Downloading Lexibank/CLDF data is approved. Research outputs live in the
  repository (`docs/research/`, `reference/calibrations/`).

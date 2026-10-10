# D2: the cross-family clock pilot

> Run 2026-10-10. The rules were fixed in advance in `clock-pilot/plan.md`;
> one deviation, made after the first run, is recorded there. Code:
> `tools/clock_pilot/` (`cd tools && python3 -m clock_pilot.run
> --no-lexstat`). Outputs: `clock-pilot/calibration_matches.tsv`,
> `results.tsv`, `lofo.tsv`, `lofo_one_sided.tsv`; first-run outputs are in
> `*_run1.tsv`.

## Question

For a family with no calendar calibration, how wide must an honest calendar
interval be? `docs/roadmap.md` says the width is set by how much the clock
varies between families (σ_f) and between nodes (τ). This pilot measures
both on calibrated splits in expert-coded Lexibank data, then checks by leave
one family out (LOFO) whether a transported clock predicts known ages better
than trivial baselines.

## What was used

- **Calibrations.** 78 grade A/B node-type rows, deduplicated into 53
  Glottolog nodes. 32 of the nodes are two-sided. Creoles, circular
  (Bayesian-derived) dates, and nodes that cannot be measured with living
  languages were excluded.
- **Cognates.** The 59 expert-cognate Lexibank datasets, with Glottolog v5.0
  for classification. Six datasets have about one cognate set per concept
  and were excluded (deviation 1): `carvalhopurus`, `chacontukanoan`,
  `ratcliffearabic`, `simsrma`, `yamfinder`, `zhangrgyalrong`.
- **Result of matching.** 15 two-sided nodes in **F = 6 families**:
  Austronesian, Dravidian, Indo-European, Pama-Nyungan, Tungusic and Uralic.
  Nine one-sided nodes entered a separate censored fit.
- **Estimator.** The pooled cross-pair match proportion p gives a distance
  `d = −ln p = 2rT` (T1, a single rate), and `β = ln(d / 2T)`. The variant
  with gamma among-concept heterogeneity (α = 1) is reported beside it.
  Uncertainty comes from a bootstrap over concepts and from the calibration
  interval.

## Results

| quantity (T1) | estimate | 95% interval |
|---|---|---|
| μ, log rate per year | −8.59 | −9.14 to −7.96 |
| σ_f, between families | 0.33 | 0.01 to 1.18 |
| τ, between nodes within a family | 0.64 | 0.42 to 1.04 |
| **σ_tot = √(σ_f² + τ²)** | **0.77** | **0.50 to 1.40** |

These are posterior medians with HalfNormal(1) priors. Neither the REML
fit nor a HalfNormal(0.5) prior changes the picture: REML gives σ_f = 0 and
τ = 0.64, and the tighter prior gives σ_tot = 0.70. The gamma variant gives
σ_tot = 0.91. The censored fit, which adds the one-sided nodes, gives
σ_f = 0 and τ = 0.59.

**The spread is within families, not between them.** Nodes inside one family
scatter about as much as families do. A single clock rate does not hold even
within a well-coded family. Anything that varies from node to node goes into
τ: lineage-specific rates, contact, the age of the calibration, and sampling
of doculects. The cross-family meta-analysis cannot remove it.

### Leave one family out (family-weighted means)

| predictor | log score | 80% cov. | 95% cov. | ×/÷ width of 95% |
|---|---|---|---|---|
| learned clock (T1) | −1.35 | 0.77 | 1.00 | 7.5 |
| old default clock (−8.8, sd 0.8) | **−1.23** | 0.71 | 0.99 | 4.9 |
| family size | −1.53 | 0.70 | 0.76 | 3.8 |
| climatological (training node ages) | −1.77 | 0.71 | 0.83 | 5.4 |

These are differences in log score, with family bootstrap 95% intervals:

| comparison (clock − baseline) | difference | 95% interval |
|---|---|---|
| vs climatological | +0.42 | −0.29 to +1.56 |
| vs family size | +0.18 | −0.55 to +1.15 |
| vs old default | −0.12 | −0.24 to −0.01 |

For the root of a new family whose median is 2,000 years, the learned
clock's 95% predictive interval is ×/÷ 6.1: about **330 to 12,000 years**.

## Reading against the pre-registered rule

The pre-registered **no-go** condition (plan.md, "Reading the result") is
σ_tot > 0.64 or a clock that does not beat the climatological baseline.
**Both hold:**
- the posterior median of σ_tot is 0.77;
- the advantage over the climatological baseline has an interval that
  includes zero, and over family size likewise.

The learned clock is also slightly worse than the old Swadesh-style default.
That default has a similar centre (−8.8 against −8.6) and a similar spread
(0.8 against 0.77), so the data essentially rediscover it.

**For the calendar promise this is a no-go at the pilot level.** With
today's evidence the honest interval for an uncalibrated node is a factor of
5–8 either way at 95%. That is not decision-relevant, and it is wider than
the ×/÷ 3.5 sharpness bar set for M2. The tool should report lexical and
relative ages as its products for uncalibrated families. Calendar output
there should be labelled scenarios, which is the fallback the roadmap
already specifies.

## What could still change this

The pilot is deliberately crude, and three of its simplifications inflate τ
in ways the full engine addresses:

1. **Distances are pairwise and single-rate.** The full T1 fit with
   integrated concept-rate heterogeneity (E2.1, now in the engine) and a tree
   uses all languages jointly. The gamma variant here did not shrink τ, so
   heterogeneity alone is unlikely to be enough.
2. **Calibration quality.** Several "points" are max-leaning political
   events, and node-to-calibration assignments are uncertain. Part of τ is
   calibration error, not clock variation. E2.6 (outlier and one-sided
   calibrations) and D1.2 (curation) would let the fit separate the two.
3. **F = 6.** Every interval above is wide. With six families σ_f is barely
   constrained, so the claim that the spread is within families rests on the
   point estimates. More families need the tip-calibrated and C-grade routes
   in D1.4.

None of these is likely to bring σ_tot from 0.77 down to the 0.3–0.4 the
first roadmap assumed. That assumption was in-sample and non-T1 (Holman
2011); out of sample and under T1, the number is about twice as large.

## Implications for the roadmap

- M1 and M2 stay. They are the definitive version of this experiment, but
  their expected outcome is now the failure branch: lexical and relative
  chronology plus calendar scenarios.
- The levers that remain are family-specific soft calibrations (D4.2), which
  the user has ruled out of scope for Arawan, and modelling τ explicitly:
  clock variation across nodes, through the relaxed clock with honest priors.
- The Arawan calendar interval from the clock alone should be stated as
  roughly ×/÷ 6 at 95%, not ×/÷ 2.

## Not done

- Uniform automatic coding (LexStat) was skipped (`--no-lexstat`): the
  module was not finished before the earlier session ran out of memory. It
  remains the obvious test of whether one protocol shrinks τ.
- The first run's σ_tot of about 1.4 came almost entirely from Ratcliffe's
  Arabic. That dataset codes cognacy relative to Standard Arabic, with one
  set per concept. Its exclusion is the only deviation from the plan.

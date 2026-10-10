# M1/M2 at pilot scale: the clock with the full engine — pre-registration

Written 2026-10-10, before any subtree was dated. At the time of writing the
only results looked at were those of the D2 pilot (`../clock-pilot.md`), whose
numbers this experiment is meant to check. Code: `tools/clock_engine/`
(`cd tools && python3 -m clock_engine.run`). Changes made after any engine
result was seen are listed at the end under "Deviations".

## Question

The D2 pilot estimated each calibrated node's log clock rate from pooled
pairwise match proportions under a single rate, `β = ln(d / 2T)`, and found
σ_tot = 0.77, most of it within-family node scatter (τ = 0.64): a no-go on its
pre-registered rule. Its report names the estimator as a possible source of
inflated τ. This experiment repeats the pilot with the node's lexical age
estimated by `clepsydra date-tree`: all sampled languages jointly on a tree,
among-concept rate heterogeneity integrated out, contamination (ε) sampled.
Is τ, and σ_tot, smaller?

## Nodes (unchanged from the pilot)

Calibration rows, curation, deduplication, matching, dataset exclusions
(including the pilot's deviation 1, degenerate codings, and its deviations
1–3 on modern "Middle" dialects, rake nodes and the family-primary rule) are
those of the pilot, taken by calling its code (`tools/clock_pilot`:
`calibrations`, `glottolog`, `lexibank`, `match`, and `run.select` with the
`family` rule). ABVD is not used. The result is expected to be the pilot's 15
two-sided nodes in 6 families plus its 9 one-sided nodes; the code checks the
selection against the pilot's `calibration_matches.tsv` and stops if it
differs.

## Subtree per node (roadmap D5, target-like)

From the node's primary dataset only, after the pilot's coding rule (first
form per doculect and concept that has a cognate judgement; ancient
doculects dropped):

1. **Ingroup**: the doculects the pilot's match placed on the node's sides
   (`Match.labels`). At most 6:
   - split and stem nodes: up to 3 per side;
   - crown nodes: round-robin over the Glottolog children of the crown that
     are represented (children ordered by their best doculect's number of
     coded concepts, descending, then glottocode), taking the next best
     doculect of each child in turn, until 6 are chosen or none remain.
   Within a side or child, doculects are ranked by number of coded concepts
   (descending), then doculect id.
2. **Outgroup**: up to 2 doculects of the same dataset (modern, with a
   glottocode in Glottolog, with coded concepts) that are not under the
   node's MRCA and whose Glottolog MRCA with the node is the node's parent.
   Ranked as above. If there is none, the node is the root of the subtree.
   No fallback to more distant outgroups.
3. **Topology**: the Glottolog v5.0 classification pruned to the chosen
   doculects. Each doculect is a tip below its own glottocode (so a doculect
   coded at a node that also has a chosen descendant is a sister of that
   descendant). Unary nodes are collapsed; multifurcations are kept. Internal
   nodes are labelled with their glottocode, and the calibrated node with
   `target`. The code checks that the `target` node is the Glottolog MRCA of
   the ingroup and equals the pilot's node MRCA.
4. **Tips** are named `L_<doculect id>`, with every character outside
   `[A-Za-z0-9_]` replaced by `_` (a numeric suffix if two collide).
5. **forms.tsv**: columns `language_id, concept_id, cognate_id`, one row per
   coded (doculect, concept), concepts being the dataset's own parameters.
6. **calibrations.tsv**: one `fixed_tip` row at 0 (hard) per tip, and nothing
   else. No calendar information enters the run.

## Dating

`clepsydra date-tree` with the default model: concept rates on
(log-normal), ε sampled (Beta(1, 19)), uniform age prior with the default
lexical-root prior LogNormal(−1.2, 1.5), clock prior at its default (it sets
only calendar ages, which are not used). `--warmup 3000 --samples 3000
--chains 4 --seed 20261010`. Runs are strictly sequential, each in its own
memory-limited scope (1 GB).

- If the run's convergence check fails (exit status 3), it is rerun once with
  `--samples 6000` (same seed), and the rerun is used. If that fails too, the
  rerun's result is kept and flagged `converged = 0`; flagged nodes stay in
  the main fit and a sensitivity fit without them is reported.
- The binary used is a fixed copy of `build/clepsydra` (commit 8654480 era,
  0.2.0). It has no `--write-draws`, and `--trace-draws` is not written to a
  file by `date-tree`, so per-draw lexical ages are not available.

**Lexical age and its uncertainty.** From `posterior_nodes.tsv`, row
`target`: `L = lexical_median`, and on the log scale

    ln L ≈ N( ln lexical_median,  s_L² ),   s_L = (ln hpd90_hi − ln hpd90_lo) / (2 · 1.645).

This treats the posterior of ln L as normal and the reported 90% HPD as an
equal-tailed 90% interval on the log scale. The approximation is used
because it is the simplest one available from the outputs; it affects only
the measurement variance, which in the pilot was small next to τ.

## Rate per node

`β_i = ln L_i − ln T_i`. T_i follows the pilot exactly: +50 years to BP
ages; for an interval, ln T centred on the log midpoint with variance
(ln hi − ln lo)²/12; for a point, SD 0.15. There is no factor 2: L is a
node-to-tip quantity, not a pair distance. Measurement variance
s_i² = s_L² + Var(ln T).

Note on scale: under the engine's concept-rate mixture ρ is the mean rate
across concepts, so β_engine is expected to sit above the pilot's single-rate
β_T1 (closer to its gamma variant). The comparison of interest is the spread
(τ, σ_f), not μ. The old default clock (−8.8, sd 0.8) is a single-rate
Swadesh value and is applied to the engine's L as is; this puts it at a
disadvantage that is stated, not corrected.

## Meta-analysis and LOFO (the pilot's code)

- The pilot's `model.fit` (REML with profile intervals; Bayesian grid with
  HalfNormal(1), and HalfNormal(0.5) as sensitivity) on the two-sided nodes:
  μ, σ_f, τ, σ_tot with 95% intervals; depth slope.
- The censored secondary fit (two-sided + one-sided nodes) with the pilot's
  `model.censored_fit`; one-sided rows use the engine's L with variance s_L².
- LOFO with the pilot's `lofo` predictors and `lofo.score`: Bayesian clock
  (grid mixture), REML plug-in t, old default clock, family size,
  climatological; one-sided held-out nodes scored by censored log score.
  Family-weighted means; family bootstrap (2000, the pilot's seed) of the
  log-score difference of the Bayesian clock against each baseline.
- Width factor of the 95% predictive for a new family's node with median
  2000 years, at the median s_L of the two-sided nodes.
- Node by node: β_engine beside the pilot's β_T1 and β_g1.

## Reading the result (the pilot's rule, unchanged)

- σ_tot = √(σ_f² + τ²); ×/÷ 3.5 at 95% corresponds to σ_tot ≈ 0.64.
- **Go**: posterior median σ_tot ≤ 0.5, LOFO 95% coverage ≥ 80%
  (family-weighted), and the clock beats both the climatological and the
  family-size baselines in family-weighted mean log score.
- **No-go**: posterior median σ_tot > 0.64, or the clock does not beat the
  climatological baseline.
- **Inconclusive**: anything else, including F < 6.

"Beats" is read as in the pilot: the family-weighted mean difference is
positive; the bootstrap interval is reported beside it.

## Deviations

Made after the subtrees were built and inspected (doculect lists and
topologies only), **before any dating run**:

1. **Jurchen excluded.** The only outgroup available for TUN-01
   (Manchu–Xibe) was Jurchen (`jurc1239`, oskolskayatungusic), a historical
   language (12th–17th c.) that the pilot's ancient-name rule does not catch.
   As a tip fixed at age 0 it would be dated as a modern language. It is
   now treated as ancient (`subtree.EXTRA_ANCIENT`), so TUN-01 is a two-tip
   subtree with the node at the root. The pilot did not use Jurchen (it was
   on neither side of the split), so the node selection is unaffected.

Noted, not changed: AN-020 (Micronesian) uses `mcd`, a comparative
dictionary with 7,312 "concepts", most coded in one or two doculects. The
pilot used it under the same rules, and the plan keeps the dataset's own
parameters, so it is used as is.

Technical, no effect on the selection: `reference/calibrations/calibrations.tsv`
changed after the pilot (D1.2, commit 7f38b1f, regraded the two circular rows
AU-02 and NA-020, which the pilot excluded anyway). The pilot's loader rejects
the new grades, so the code reads the table as of commit f92ccc1, the version
the pilot used. The node selection is checked against the pilot's
`calibration_matches.tsv` and is identical.

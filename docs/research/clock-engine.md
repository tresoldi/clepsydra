# M1/M2 at pilot scale: the cross-family clock with the full engine

> Run 2026-10-10. Pre-registration: `clock-engine/plan.md`, written before
> any subtree was dated; its two deviations were made before any dating run.
> Code: `tools/clock_engine/` (`cd tools && python3 -m clock_engine.run
> --clepsydra BINARY`). Outputs in `clock-engine/`: `nodes.tsv`,
> `results.tsv`, `lofo.tsv`, and `subtrees/<node>/` (tree, forms,
> calibrations, doculect list, clepsydra output, run log).

## Question

The D2 pilot (`clock-pilot.md`) estimated each calibrated node's log clock
rate from pooled pairwise match proportions under a single rate, and found
σ_tot = 0.77, mostly within-family node scatter (τ = 0.64). It named its
estimator as a possible cause of an inflated τ. This run keeps the pilot's
nodes, datasets, meta-analysis and LOFO, and replaces only the per-node
estimate: each node's lexical age L is estimated by `clepsydra date-tree`
on a small subtree, with all sampled languages jointly, integrated
among-concept rate heterogeneity and sampled contamination. Does τ, and
σ_tot, shrink?

## What was done

- **Nodes.** The pilot's selection, by its own code, checked row by row
  against its `calibration_matches.tsv`: 15 two-sided nodes in 6 families,
  plus 9 one-sided nodes for the censored fit. The calibration table is read
  as of commit f92ccc1, the version the pilot used (D1.2 later regraded two
  rows that the pilot excluded anyway).
- **Subtrees** (roadmap D5). From the node's primary dataset: up to 6
  ingroup doculects (3 per side, or a round-robin over the crown's
  children), ranked by coded concepts, plus up to 2 outgroups whose
  Glottolog MRCA with the node is the node's parent. The topology is
  Glottolog pruned to these tips, with polytomies kept. 18 of 24 nodes have
  an outgroup. Six are subtree roots: Romance, Anglo-Frisian, Modern
  Goidelic, Insular Celtic, Polynesian, and Manchu–Xibe.
- **Deviation 1.** Jurchen, a historical language that the pilot's
  ancient-name rule misses, would have been the Tungusic outgroup as a
  modern tip. It was excluded, so TUN-01 is a two-tip tree.
- **Dating.** `date-tree` with the default model and only `fixed_tip` = 0
  calibrations, so no calendar information enters. Settings: `--warmup 3000
  --samples 3000 --chains 4 --seed 20261010`. Runs were sequential, each
  capped at 1 GB.
- **Lexical age and its uncertainty.** The binary has no per-draw output, so
  ln L ~ N(ln lexical_median, s_L²), with s_L = (ln hi − ln lo)/(2·1.645)
  from the 90% HPD.
- **Rate.** β = ln L − ln T, with T, its +50 y offset and its variance
  exactly as in the pilot.
- **Statistics.** The pilot's code throughout. The engine's estimate enters
  through the pilot's interface (log_d = ln L + ln 2, s_lex = s_L), so every
  formula downstream is unchanged.

**Compute.** 24 nodes and 29 clepsydra runs took 87 min of wall time.
Building the subtrees took 10 s and the analysis 6 min, for about 1.6 h in
all. Most nodes took 1–6 min. AU-01 took 12 min, and AN-020 (mcd, 7,312
dictionary "concepts") took 28 min across its two runs.

**Reruns.** Five runs failed the convergence check and were rerun with
6,000 samples: IE-GER-04, IE-CEL-02, TUN-01, DRA-03 and AN-020. Four of the
five are two-tip trees. Three passed on rerun. IE-GER-04 and IE-CEL-02 still
fail, both on `log_sigma_c` only (R̂ 1.011, ESS 363 and similar). Two tips
barely identify the spread of concept rates. The target age itself has
R̂ ≤ 1.005 and bulk ESS > 1,000 in both. Both nodes are flagged and kept,
and a fit without them is reported.

## Results

### Meta-analysis (15 two-sided nodes, F = 6)

Posterior medians (95% interval) with HalfNormal(1) priors:

| quantity | engine | pilot T1 | pilot gamma α=1 |
|---|---|---|---|
| μ, log rate per year | −8.32 (−8.90, −7.62) | −8.59 (−9.14, −7.96) | −8.28 |
| σ_f | 0.47 (0.03, 1.34) | 0.33 (0.01, 1.18) | 0.37 |
| τ | **0.46** (0.24, 0.84) | 0.64 (0.42, 1.04) | 0.77 |
| **σ_tot** | **0.70** (0.40, 1.43) | 0.77 (0.50, 1.40) | 0.91 |
| REML σ_f / τ | 0.44 / 0.41 | 0.00 / 0.64 | 0.00 / 0.77 |
| ×/÷ 95% width, 2,000 y node | **5.8** (330–11,200 y) | 6.1 (320–11,900 y) | 8.2 |

Sensitivity analyses:

| variant | σ_f | τ | σ_tot |
|---|---|---|---|
| HalfNormal(0.5) priors | 0.37 | 0.44 | 0.62 (0.37, 1.07) |
| converged nodes only (13) | 0.47 | 0.50 | 0.74 (0.42, 1.47) |

Other fits:

- **Censored fit** (adds the 9 one-sided nodes, ML): σ_f = 0.46, τ = 0.35,
  μ = −8.33.
- **SD of family means:** 0.79 (0.33–1.73).
- **Depth slope:** −0.15 (−0.63, 0.34).

### Node by node

β_engine is ln L − ln T. β_pilot is the pilot's ln(d / 2T): single rate
(T1), and gamma α = 1 in brackets.

| node | family | dataset | T (y) | L | β_engine | β_pilot T1 [g1] |
|---|---|---|---|---|---|---|
| IE-ROM-01 Romance | IE | iecor | 1710 | 0.345 | −8.51 | −8.73 [−8.44] |
| IE-GER-04 Anglo-Frisian* | IE | iecor | 1520 | 0.351 | −8.37 | −8.56 [−8.26] |
| IE-GER-05 North Germanic | IE | iecor | 1090 | 0.106 | −9.24 | −9.07 [−8.94] |
| IE-CEL-01 Brythonic | IE | iecor | 1470 | 0.276 | −8.58 | −8.65 [−8.38] |
| IE-CEL-02 Modern Goidelic* | IE | iecor | 929 | 0.154 | −8.71 | −8.79 [−8.64] |
| IE-SLA-02 East Slavic | IE | iecor | 663 | 0.065 | −9.23 | −9.04 [−8.96] |
| URA-FIN-01 Finnic | Uralic | uralex | 1170 | 0.256 | −8.43 | −8.54 [−8.30] |
| URA-SAA-01 Saami | Uralic | uralex | 1724 | 0.309 | −8.63 | −8.79 [−8.51] |
| URA-PER-01 Permian | Uralic | uralex | 901 | 0.476 | −7.55 | −7.90 [−7.54] |
| TUN-01 Manchu–Xibe | Tungusic | oskolskaya | 188 | 0.121 | −7.35 | −7.40 [−7.28] |
| AN-01 Malayo-Polynesian | Austronesian | blust | 4243 | 0.950 | −8.40 | −8.69 [−7.90] |
| AN-014 East Polynesian | Austronesian | walworth | 853 | 0.356 | −7.78 | −7.87 [−7.53] |
| AN-020 Micronesian | Austronesian | mcd | 2040 | 0.339 | **−8.70** | **−10.20** [−10.12] |
| DRA-03 Tamil–Malayalam | Dravidian | dravlex | 980 | 0.086 | −9.34 | −9.33 [−9.24] |
| AU-01 Wati | Pama-Nyungan | bowernpny | 2798 | 1.717 | −7.40 | −8.03 [−6.99] |

\* still unconverged in `log_sigma_c` after the rerun (target age converged).

Engine and pilot agree closely (r = 0.81). Fourteen of the 15 nodes move by
less than 0.65, and most of them sit between the pilot's T1 and gamma
values, as expected when concept rates are integrated out. The one large
change is AN-020 (Micronesian), which moves by +1.5. In `mcd`, a comparative
dictionary, the pilot's pooled cross-pair proportion is dominated by
sparsely coded meanings. The joint fit, with ε sampled and incompatible
concepts excluded, puts the node in line with the rest of Austronesian. The
raw SD of β over nodes barely changes: 0.64 against the pilot's 0.67.

**Where the smaller τ comes from.** The following check was not
pre-registered. Dropping AN-020 alone from the *pilot's* fit gives
σ_f = 0.48, τ = 0.39 and σ_tot = 0.66. Dropping it from the engine's fit
gives 0.49, 0.46 and 0.72. The drop in τ is therefore the correction of one
outlier node. On the other 14 nodes the engine's scatter is no smaller than
the pilot's (slightly larger, in fact, through the two-tip nodes' wider
s_L). Some of the within-family spread has moved into σ_f, so the
partition between the two is not stable at F = 6, but the total is.

### Leave one family out (15 held-out two-sided nodes, family-weighted)

| predictor | log score | 80% cov. | 95% cov. | ×/÷ 95% width |
|---|---|---|---|---|
| learned clock, engine (Bayes) | −1.33 | 0.76 | 1.00 | 7.1 |
| learned clock, engine (REML t) | −1.38 | 0.71 | 1.00 | 7.9 |
| old default (−8.8, sd 0.8) on engine L | −1.37 | 0.63 | 0.94 | 5.3 |
| family size | −1.53 | 0.70 | 0.76 | 3.8 |
| climatological | −1.77 | 0.71 | 0.83 | 5.4 |
| *pilot: learned clock T1 (Bayes)* | *−1.35* | *0.77* | *1.00* | *7.5* |
| *pilot: old default* | *−1.23* | *0.71* | *0.99* | *4.9* |

Differences in log score, with family bootstrap 95% intervals (F = 6):

| clock (engine) − baseline | difference | 95% interval | pilot |
|---|---|---|---|
| vs climatological | +0.45 | −0.49 to +1.89 | +0.42 (−0.29, +1.56) |
| vs family size | +0.21 | −0.58 to +1.40 | +0.18 (−0.55, +1.15) |
| vs old default | +0.04 | −0.31 to +0.35 | −0.12 (−0.24, −0.01) |

Two further results:

- **Calibration.** Σz² = 8.4 over 15 nodes for the engine clock, against an
  expectation of 15. The 95% intervals cover every held-out node. The
  predictive is, if anything, too wide, which is consistent with τ and σ_f
  being poorly separated at F = 6.
- **One-sided held-out nodes** (9, mean censored log score): engine clock
  −0.55, old default −0.30, family size −0.66, climatological −0.87.

The old default is no longer better than the learned clock. This is largely
a scale effect: under concept-rate mixing ρ is the mean rate, which sits
about 0.3 above the single-rate value the −8.8 default was set on. That
puts the default off-centre (family PIT 0.31), as the plan anticipated.

## Reading against the pre-registered rule

The rule is the pilot's, unchanged:
- **Go** needs a posterior median σ_tot ≤ 0.5, LOFO 95% coverage ≥ 80%, and
  a clock that beats both the climatological and the family-size
  baselines.
- **No-go** is σ_tot > 0.64, or a clock that does not beat the
  climatological baseline.

| criterion | value | |
|---|---|---|
| σ_tot posterior median | 0.70 | > 0.64: no-go condition met |
| LOFO 95% coverage (family-weighted) | 1.00 | ≥ 0.80 |
| beats climatological | +0.45 (CI includes 0) | positive point estimate |
| beats family size | +0.21 (CI includes 0) | positive point estimate |

**Reading: no-go**, on σ_tot. It does not depend on the unconverged nodes:
without them σ_tot is 0.74. Only the tighter HalfNormal(0.5) prior brings
the median down to 0.62. That still sits above the 0.5 needed for a go,
which would make it inconclusive. The skill differences have the right
sign but intervals that include zero, as in the pilot.

## What this says

- **The pilot's estimator was not the main source of τ.** With the full
  engine, joint estimation on a tree, integrated concept rates and sampled
  contamination, the per-node rates are essentially the pilot's. The
  exception is one node whose pairwise estimate on a comparative dictionary
  was an outlier.
- **σ_tot ≈ 0.7 stands.** A transported clock gives an uncalibrated
  2,000-year node an interval of about ×/÷ 6 at 95% (330 to 11,000 years),
  against the ×/÷ 3.5 that M2 requires.
- **The clock carries some information.** It is calibrated, if wide, and
  beats the trivial baselines on point estimates, but not reliably at F = 6.
  The expected product for uncalibrated families remains lexical and
  relative chronology plus labelled calendar scenarios, which is the
  roadmap's failure branch for M2.

## Limitations

- **Uncertainty of L.** It comes from the 90% HPD under a log-normal
  approximation, not from per-draw ages. The engine's s_L is small next to
  τ (median 0.14) except on the two-tip nodes (0.40–0.45), so this matters
  little for σ_tot.
- **Two-tip subtrees** (Anglo-Frisian, Modern Goidelic, Manchu–Xibe) are
  pairwise estimates in all but name. Tamil–Malayalam has a single outgroup.
- **Topologies are Glottolog's.** Rake nodes (Malayo-Polynesian, Iranian)
  are left as polytomies. The topology is not estimated (E3.4).
- **F = 6.** σ_f and τ are not separately identified, and σ_f's interval
  runs from 0.03 to 1.34. More families need D1.4.
- **Single-rate rule.** μ is not comparable with the pilot's T1 μ, but the
  spread is, and the spread is what the rule judges.

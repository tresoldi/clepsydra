# Arawan: the M3 result

Roadmap item M3: KAHD Arawan dated with every layer reported. The inputs are
those of `arawan-acceptance.md`: 9 doculects, 253 concepts, the MCC topology
with labelled internal nodes, and tip dates. Arua (Chandless 1869) is fixed
at 131 BP. Paumari 1833 has a soft U(165, 170) tip date and the empirical
list-error prior from D1.3 (Beta with mean 0.088, concentration 13). No
expert calibration exists (D4.2 is out of scope).

Reproduce with `sh tools/arawan_m3/run.sh OUTDIR` (about 2 h on 4 cores,
every run under a 1 GB cap). The fixed-tree runs used 4 chains × (4,000
warm-up + 4,000 draws), seed 42. All runs date from 2026-10-11 and use a
binary built from main. The averaging step ran on 19197ca.

## 1. Lexical ages: the result the data support

Lexical age is the expected number of replacements per concept, at the
mean concept rate, from the node to the present. It needs no clock. Under
every variant below it is stable to the second decimal.

| node | clade | lexical median (90% HPD), strict | relaxed clock | modern lists only |
|---|---|---|---|---|
| root | Arawan | 0.646 (0.544–0.763) | 0.639 (0.496–0.784) | 0.628 (0.526–0.735) |
| n6 | non-Paumari | 0.617 (0.516–0.727) | 0.593 (0.458–0.729) | 0.601 (0.509–0.708) |
| n5 | Madi–Arawa + Deni–Kulina | 0.343 (0.281–0.407) | 0.347 (0.269–0.436) | 0.334 (0.278–0.399) |
| n3 | Madi + Arawa | 0.279 (0.190–0.361) | 0.287 (0.186–0.385) | — |
| n4 | Deni–Kulina | 0.174 (0.135–0.216) | 0.175 (0.118–0.234) | 0.177 (0.137–0.221) |
| n7 | Paumari 1833 + modern | 0.107 (0.038–0.186) | 0.116 (0.033–0.207) | — |
| n2 | Madi | 0.037 (0.022–0.052) | 0.038 (0.021–0.057) | 0.038 (0.024–0.054) |
| n1 | Jamamadi–Jarawara | 0.026 (0.013–0.038) | 0.026 (0.013–0.041) | 0.026 (0.013–0.039) |

The relative chronology is the robust part. The family splits into an early
Paumari / non-Paumari divergence (0.65). This is followed almost at once by
Suruwaha (0.62), then a long gap to the Madi–Arawa / Deni–Kulina split
(0.34). Madi is very shallow (0.04): its three lects are close to dialects
lexically. The relaxed clock widens the intervals but moves no median by
more than 0.02. Dropping the two historical lists changes the root by 0.02.

## 2. Calendar ages: what fixes the scale

The calendar scale comes from the two dated historical lists (165–170 and
131 years of known time). That is little time against a 2,000–4,000-year
tree, so the calendar root is set mostly by the rate prior. The engine
reports this: `scale_from_prior 1`, and 98.6% of the variance in the log of
the root age moves with the clock (`root_var_share_clock` 0.986), 1.4% with
the lexical age (`root_var_share_data` 0.014). With the empirical clock or the
root window the split is the same (0.98 / 0.02).

| scale source | root BP, median (90% HPD) | n5 BP | n4 (Deni–Kulina) BP | converged |
|---|---|---|---|---|
| tip dates, default rate prior | 3,157 (749–7,749) | 1,657 (366–4,112) | 847 (187–2,132) | yes (R-hat ≤ 1.003) |
| tip dates, relaxed clock | 2,930 (630–7,563) | 1,602 (352–4,153) | 803 (154–2,111) | **no** (σ_b: R-hat 1.013, ESS 618)¹ |
| tip dates + empirical clock (M1: log-rate −8.32, sd 0.77) | 2,283 (634–5,152) | 1,210 (321–2,731) | 610 (168–1,396) | yes |
| tip dates + U(300, 3000) root window | 1,982 (1,093–2,995) | 1,043 (488–1,554) | 528 (237–808) | yes |
| modern lists only, default rate prior | 4,139 (498–12,118) | 2,187 (242–6,394) | 1,151 (93–3,385) | yes |

¹ At twice the length (8,000 + 8,000, §5) both relaxed variants converge
(R-hat ≤ 1.003, ESS ≥ 1,088).

The headline calendar estimate is the empirical-clock row. Its prior comes
from the cross-family clock pilot (`clock-pilot.md`), with σ_tot 0.77
pre-registered. On that prior, the root falls at about 2,300 BP with a
90% interval of 600–5,200 BP.
The interval is close to the σ_f floor that `mathematical-analysis.md`
predicts for a family with no internal calibrations. Only expert soft
calibrations (D4.2) would narrow it. The window row adopts the
BEAST scenario's U(300, 3000) as a calibration: its answer is the window's
answer, and it is reported as a scenario, not a finding. The grade C
3,200 BP archaeological maximum would be a sensitivity run of the same
kind, and is not adopted.

## 3. Topology

Uncertainty about the topology comes from the BEAST B1 posterior: 20 trees, evenly spaced after 10%
burn-in, each dated with the tip calibrations (`date-trees`, 4 chains × 3,000
+ 3,000). All 20 converged (max R-hat 1.007, min ESS 1,213).

| clade | P(monophyletic) | MRCA BP, median (90% HPD) |
|---|---|---|
| root | 1.00 | 2,815 (619–7,295) |
| Madi (Banawa, Jamamadi, Jarawara) | 1.00 | 173 (34–460) |
| Paumari (1833 + modern) | 1.00 | 419 (166–1,142) |
| Deni–Kulina | 1.00 | 822 (152–2,123) |
| non-Paumari | 0.70 | 2,660 (521–6,839) |
| Madi + Arawa | 0.35 | 1,692 (311–4,500) |

Per-tree root medians range from 2,352 to 3,557 BP. That spread is small
next to the calendar interval of any one tree, so topology adds little
uncertainty to the root (as the roadmap expected: the root is the MRCA of
everything). It matters for the internal structure. The MCC tree's
non-Paumari clade and its Madi + Arawa grouping are supported by 70% and
35% of the sample. Dates for those two nodes should be read conditional
on the clade (the `if_mono_*` columns of `posterior_clades.tsv`).

## 4. Paumari 1833 as the held-out check

Its date was replaced by a flat U(0, 1000) BP and the run repeated, keeping
its error prior. The model places the list at 232 BP (90% HPD 0–701). The
true 165–170 lies inside. The interval is wide: the separation from modern
Paumari is 0.04 lexical units, and an 1833 list with an error probability of
about 0.19 cannot date itself to better than a few centuries. The check
fails in neither direction. Before the error model (E2.9) the same list
pulled the root to 839 BP (`arawan-acceptance.md`).

Without the 1833 date the root moves from 3,157 to 3,811 BP (560–9,311).
The calendar root is therefore sensitive to the one historical anchor it
has, and the lexical root is not (0.646 to 0.649).

## 5. Model uncertainty (E2.8)

Four models were run on the MCC tree with tip dates, at 8,000 + 8,000 draws
each, with every draw and per-concept log-likelihood written. They were
combined with `clepsydra average`, which uses PSIS-LOO over concepts and
stacking.

| model | elpd_loo (SE) | Δ vs best (SE) | p_loo | max Pareto k | stacking weight |
|---|---|---|---|---|---|
| strict, log-normal concept rates | −763.8 (28.7) | −0.3 (1.3) | 7.7 | 0.22 | 0.34 |
| strict, gamma | −764.0 (28.7) | −0.4 (1.6) | 7.6 | 0.14 | 0.00 |
| relaxed, log-normal | −763.5 (28.8) | 0 | 10.4 | 0.27 | 0.29 |
| relaxed, gamma | −763.7 (28.8) | −0.2 (1.0) | 10.3 | 0.19 | 0.36 |

The data cannot tell the four apart: every difference is well inside one
SE, and the relaxed clock buys about 2.7 effective parameters for nothing.
Every Pareto k is below 0.3, so the LOO estimates are reliable. In the stacked
posterior the root is 2,994 BP (5–95%: 1,108–9,746) and 0.631 lexical
(0.514–0.779). The between-model share of the log-age variance is 0.2% for
the calendar root and at most 9% for any lexical age (n6). The choice of
clock and rate-mixing family does not drive the Arawan result.

## 6. Adequacy

The posterior predictive panel (`adequacy.tsv`, strict run) passes on
classes per concept (ppp 0.49), on singleton share (0.61), and on every
pairwise sharing rate (ppp between 0.13 and 0.97). The exception is
**the share of concepts that no T1 history can produce**. It is 7.5%
observed against 0.3% replicated (90% at most 1.2%; ppp 0.000). These are
concepts where a cognate class reappears in non-sister lects after another
class has intervened. Borrowing, ancestral polymorphism and coding errors
all produce that pattern. The model handles such concepts through the ε
contamination term (posterior about 0.004), but it does not explain them.
This is the main known misfit. The simulation study (`simworld.md`) shows
that coding error is the violation that biases lexical ages most.

## 7. Not done

- **Coding and loan sensitivity (D4.3).** It needs an alternative coding of
  KAHD (automatic cognates or a loan-flagged list). There is none in the
  repository, and D4.4 (double coding) is out of scope.
- **Expert soft calibrations (D4.2).** Out of scope. They are the only lever
  that narrows the calendar interval.

## Summary

- **Lexical ages are robust.** The Arawan root is 0.65 lexical units
  (0.54–0.76). The split order, the relative depths, and the shallowness of
  Madi survive every variant: clock, rate model, topology, and the historical
  lists.
- **Calendar ages are bounded by what fixes the time scale.** The best
  estimate uses the empirical clock prior: root 2,300 BP (600–5,200). With
  only the default rate prior and tip dates, the median is 3,200 BP and the
  upper bound is half again as far (7,700 BP).
- **The checks behave as they should.** The held-out Paumari 1833 date is
  recovered. Model averaging finds no preference between clocks. The
  remaining misfit is the 7.5% of concepts that are incompatible with the
  tree.

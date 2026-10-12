# Arawan: the M3 result

Roadmap item M3: KAHD Arawan dated with every layer reported. The inputs are
those of `arawan-acceptance.md`: 9 doculects, 253 concepts, the MCC topology
with labelled internal nodes, and tip dates. Arua (Chandless 1869) is fixed
at 131 BP. Paumari 1833 has a soft U(165, 170) tip date and the empirical
list-error prior from D1.3 (Beta with mean 0.088, concentration 13). No
expert calibration exists (D4.2 is out of scope).

Reproduce with `sh tools/arawan_m3/run.sh OUTDIR`. The current numbers were
run on main at 80da2a4, under the default model since then: the relaxed
clock at 4 chains × (4,000 warm-up + 4,000 draws), seed 42, every run under
a 1 GB cap (on `mydev`, about 4 h on 2 cores). The empirical-clock row did
not converge at that length (σ_b R-hat 1.015) and was rerun at 8,000 + 8,000.
The strict-clock numbers of the first version of this page (2026-10-11,
before the relaxed default) are kept as the `strict` column and row.

## 1. Lexical ages: the result the data support

Lexical age is the expected number of replacements per concept, at the
mean concept rate, from the node to the present. It needs no clock. Under
every variant below it moves by at most 0.03.

| node | clade | lexical median (90% HPD), default (relaxed) | strict clock | modern lists only |
|---|---|---|---|---|
| root | Arawan | 0.641 (0.507–0.793) | 0.646 (0.544–0.763) | 0.637 (0.503–0.792) |
| n6 | non-Paumari | 0.596 (0.464–0.736) | 0.617 (0.516–0.727) | 0.595 (0.470–0.734) |
| n5 | Madi–Arawa + Deni–Kulina | 0.348 (0.269–0.439) | 0.343 (0.281–0.407) | 0.343 (0.265–0.436) |
| n3 | Madi + Arawa | 0.288 (0.185–0.386) | 0.279 (0.190–0.361) | — |
| n4 | Deni–Kulina | 0.173 (0.119–0.230) | 0.174 (0.135–0.216) | 0.176 (0.122–0.234) |
| n7 | Paumari 1833 + modern | 0.112 (0.033–0.202) | 0.107 (0.038–0.186) | — |
| n2 | Madi | 0.038 (0.022–0.058) | 0.037 (0.022–0.052) | 0.039 (0.022–0.057) |
| n1 | Jamamadi–Jarawara | 0.026 (0.012–0.041) | 0.026 (0.013–0.038) | 0.026 (0.013–0.041) |

The relative chronology is the robust part. The family splits into an early
Paumari / non-Paumari divergence (0.64). This is followed almost at once by
Suruwaha (0.60), then a long gap to the Madi–Arawa / Deni–Kulina split
(0.35). Madi is very shallow (0.04): its three lects are close to dialects
lexically. The relaxed clock widens the intervals but moves no median by
more than 0.02. Dropping the two historical lists changes the root by 0.004.

## 2. Calendar ages: what fixes the scale

The calendar scale comes from the two dated historical lists (165–170 and
131 years of known time). That is little time against a 2,000–4,000-year
tree, so the calendar root is set mostly by the rate prior. The engine
reports this: `scale_from_prior 1`, `rate_contraction` 0.32, and 97.4% of
the variance in the log of the root age moves with the clock
(`root_var_share_clock` 0.974), 2.6% with the lexical age. With the
empirical clock or the root window the split is the same (0.97 / 0.03 and
0.96 / 0.04).

| scale source | root BP, median (90% HPD) | n5 BP | n4 (Deni–Kulina) BP | converged |
|---|---|---|---|---|
| tip dates, default (relaxed clock) | 3,092 (671–7,755) | 1,676 (324–4,230) | 834 (153–2,125) | yes (R-hat ≤ 1.005, ESS ≥ 678) |
| tip dates, strict clock | 3,157 (749–7,749) | 1,657 (366–4,112) | 847 (187–2,132) | yes (R-hat ≤ 1.004) |
| tip dates + empirical clock (M1: log-rate −8.32, sd 0.77), 8,000 + 8,000 | **2,235 (563–5,171)** | 1,212 (317–2,826) | 598 (144–1,431) | yes (R-hat ≤ 1.003, ESS ≥ 1,543) |
| tip dates + U(300, 3000) root window | 1,948 (1,051–3,027) | 1,058 (479–1,655) | 523 (203–843) | yes |
| modern lists only, default rate prior | 4,250 (449–12,225) | 2,291 (213–6,566) | 1,180 (113–3,405) | yes |

The primary result is the lexical chronology of §1. The calendar row that
rests on the most external information is the empirical-clock one, but it is
a scenario conditioned on a clock prior, not a finding of the cognates: the
clock pilot (`clock-pilot.md`, σ_tot 0.77 pre-registered) was a no-go for
sharp uncalibrated dates. On that prior, the root falls at about 2,200 BP
with a 90% interval of 560–5,200 BP, as under the strict clock (2,283,
634–5,152).
The interval is close to the σ_f floor that `mathematical-analysis.md`
predicts for a family with no internal calibrations. Only expert soft
calibrations (D4.2) would narrow it. The window row adopts the
BEAST scenario's U(300, 3000) as a calibration: its answer is the window's
answer, and it is reported as a scenario, not a finding. The grade C
3,200 BP archaeological maximum would be a sensitivity run of the same
kind, and is not adopted.

## 3. Topology

Uncertainty about the topology comes from the BEAST B1 posterior: 20 trees,
evenly spaced after 10% burn-in, each dated with the tip calibrations
(`date-trees`, whose default is the strict clock at 4 chains × 2,000 +
2,000 per tree). 18 of 20 converged and are pooled; the other two (R-hat
1.011 and 1.015) are listed in `trees.tsv` but left out of the summaries.

| clade | P(monophyletic) | MRCA BP, median (90% HPD) |
|---|---|---|
| root | 1.00 | 2,858 (696–7,386) |
| Madi (Banawa, Jamamadi, Jarawara) | 1.00 | 174 (34–461) |
| Paumari (1833 + modern) | 1.00 | 432 (166–1,158) |
| Deni–Kulina | 1.00 | 821 (166–2,119) |
| non-Paumari | 0.67 | 2,712 (625–7,073) |
| Madi + Arawa | 0.39 | 1,653 (318–4,379) |

Per-tree root medians range from 2,313 to 3,423 BP. That spread is small
next to the calendar interval of any one tree, so topology adds little
uncertainty to the root (as the roadmap expected: the root is the MRCA of
everything). It matters for the internal structure. The MCC tree's
non-Paumari clade and its Madi + Arawa grouping are supported by 67% and
39% of the sample. Dates for those two nodes should be read conditional
on the clade (the `if_mono_*` columns of `posterior_clades.tsv`).

## 4. Paumari 1833 as the held-out check

Its date was replaced by a flat U(0, 1000) BP and the default run repeated,
keeping its error prior. The model places the list at 260 BP (90% HPD
0–751). The true 165–170 lies inside. The interval is wide: the separation
from modern Paumari is 0.04 lexical units, and an 1833 list with an error
probability of about 0.19 cannot date itself to better than a few
centuries. The check fails in neither direction. Before the error model
(E2.9) the same list pulled the root to 839 BP (`arawan-acceptance.md`).
Under the strict clock the same check gave 232 BP (0–701).

Without the 1833 date the root moves from 3,092 to 3,729 BP (478–9,441).
The calendar root is therefore sensitive to the one historical anchor it
has, and the lexical root is not (0.641 to 0.645).

Classical glottochronology fails this check by a factor of six to fifteen
(`arawan-glottochronology.md`).

## 5. Model uncertainty (E2.8)

Four models were run on the MCC tree with tip dates, at 4,000 + 4,000 draws
each (all four converged: R-hat ≤ 1.007, ESS ≥ 678), with every draw and
per-concept log-likelihood written. They were combined with `clepsydra
average`, which uses PSIS-LOO over concepts and stacking.

| model | elpd_loo (SE) | Δ vs best (SE) | p_loo | max Pareto k | stacking weight |
|---|---|---|---|---|---|
| strict, log-normal concept rates | −763.8 (28.7) | −0.4 (1.3) | 7.7 | 0.39 | 0.28 |
| strict, gamma | −764.1 (28.6) | −0.6 (1.6) | 7.7 | 0.25 | 0.00 |
| relaxed, log-normal (the default) | −763.4 (28.8) | 0 | 10.3 | 0.30 | 0.50 |
| relaxed, gamma | −763.7 (28.8) | −0.3 (1.0) | 10.4 | 0.20 | 0.22 |

The data cannot tell the four apart: every difference is well inside one
SE, and the relaxed clock buys about 2.6 effective parameters for nothing.
Every Pareto k is below 0.4, so the LOO estimates are reliable. In the stacked
posterior the root is 3,075 BP (5–95%: 1,119–10,022) and 0.636 lexical
(0.519–0.784). The between-model share of the log-age variance is 0.2% for
the calendar root and at most 6% for any lexical age (n6). The choice of
clock and rate-mixing family does not drive the Arawan result.

## 6. Adequacy

The posterior predictive panel (`adequacy.tsv`, default run) passes on
classes per concept (ppp 0.46), on singleton share (0.66), and on every
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
- **Relaxed clock over the tree sample.** §3 uses `date-trees`' strict
  default; a relaxed run over 20 trees costs about 7 h on `mydev`. §5 shows
  the clock choice moves the root by about 2%.

## Summary

- **Lexical ages are robust.** The Arawan root is 0.64 lexical units
  (0.51–0.79). The split order, the relative depths, and the shallowness of
  Madi survive every variant: clock, rate model, topology, and the historical
  lists.
- **Calendar ages are bounded by what fixes the time scale.** As a scenario
  under the empirical clock prior, the root is 2,200 BP (560–5,200). With
  only the default rate prior and tip dates, the median is 3,100 BP and the
  upper bound is half again as far (7,800 BP). Neither is a finding of the
  cognates: the two historical lists move the rate by about a third of its
  prior variance.
- **The checks behave as they should.** The held-out Paumari 1833 date is
  recovered. Model averaging finds no preference between clocks. The
  remaining misfit is the 7.5% of concepts that are incompatible with the
  tree.
- **The move to the relaxed default changes almost nothing here.** Every
  median above is within about 5% of its strict-clock value.
- **Classical glottochronology on the same data** gives roots from 1,750 to
  6,600 BP depending on list and formula, dates a child node older than its
  parent on this tree, and misdates the 1833 list by a factor of six or more
  (`arawan-glottochronology.md`). Its estimates fall inside clepsydra's
  intervals.

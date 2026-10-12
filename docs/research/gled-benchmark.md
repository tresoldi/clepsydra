# Clepsydra on GLED, against published dates and ASJP

Question: when clepsydra is run on public cognate data for families whose age
has been estimated by others, how close does it come to the published dates,
and how does that compare with the ASJP method (Holman et al. 2011)?

This is a blind test. No calibration is given unless stated, tips are modern,
and the rate prior is the default (log r = −8.80, sd 0.80) or one of two
alternatives below. Run: `tools/gled_benchmark/` (`python3 -m
gled_benchmark.runner`); result tables in `docs/research/gled-benchmark/`.

## Setup

- **Cognate data:** GLED (release 20221127; Tresoldi), 40 concepts, 6,572
  doculects, automatic cognate sets. Ancient doculects are dropped. One
  doculect per Glottolog code (the best covered, at least 25 concepts). At
  most 24 tips per node (30 per family root), spread over the subgroups.
- **Topology:** UPGMA on the shared-cognate distance, not a published tree.
  The dated node is the root of that subtree.
- **Clepsydra:** `date-tree`, relaxed clock, 4 chains × 4,000 + 4,000, no
  calibrations. Three rate priors: `default`; `empirical` (the clock from
  the Arawan pilot, −8.32 / 0.77; partly fitted on these families, so not
  independent); `LOFO` (leave-one-family-out: mean and sd of the log rates
  implied by the known-age nodes of *other* families).
- **ASJP:** ASJP v21 (Lexibank CLDF, commit cf4bf0d), the same 40 concepts,
  the same Glottolog codes. LDND, similarity s = 1 − mean LDND over the pairs
  across the node, t = (ln s − ln 0.92) / (2 ln 0.72) millennia, the published
  constants. Reproduction check: Holman et al. report 4,348 (IE), 3,178
  (Uralic), 3,404 (Turkic), 1,564 (Japonic), 5,261 (Sino-Tibetan), 3,633
  (Austronesian), 4,295 (Pama-Nyungan); this run gives 4,212, 3,346, 3,033,
  1,783, 4,174, 3,349, 4,600. **The ASJP side is in-sample**: those
  constants were fitted partly on these nodes.
- **Known ages (29 nodes):** the independent attestation and settlement
  intervals of `reference/calibrations/curated_nodes.tsv` that GLED covers:
  18 two-sided windows, 11 one-sided bounds.
- **Published family roots (12):** Bayesian estimates from
  `literature-lexical-dating.md` Part 2 (status "paper" unless noted). The
  window is the union of their reported 95% intervals, the centre the
  geometric mean of their central values. They are model outputs with their
  own calibrations, a weaker yardstick than the dated nodes. Indo-European
  and Semitic are modern-only here (the published roots include Anatolian
  and Akkadian).

## 1. Known-age nodes

Ages in years before 2000, median and 90% HPD.

| node | known age (BP) | tips | clepsydra default | LOFO prior | ASJP |
|---|---|---|---|---|---|
| IE-ROM-01 | 1500–1950 | 15 | 738 (43–2320) | 642 (88–1814) | 1659 |
| IE-GER-04 | 1400–1650 | 4 | 1272 (110–3979) | 1109 (126–3064) | 1619 |
| IE-GER-05 | 950–1250 | 5 | 979 (55–3214) | 884 (81–2563) | 1173 |
| IE-CEL-03 | ≥ 1700 | 6 | 4501 (396–13061) | 3769 (589–9771) | 3597 |
| IE-SLA-01 | ≥ 1200 | 14 | 982 (72–2939) | 877 (122–2376) | 971 |
| IE-II-01 | ≥ 3000 | 24 | 3940 (206–13713) | 3515 (315–10881) | 3506 |
| IE-II-03 | ≥ 2600 | 18 | 3686 (207–11941) | 3254 (290–9710) | 2662 |
| IE-II-06 | 500–750 | 7 | 788 (66–2502) | 691 (73–1872) | 532 |
| URA-FIN-01 | 720–1900 | 8 | 1148 (67–3514) | 1014 (101–2767) | 918 |
| URA-SAA-01 | 1450–2050 | 6 | 1420 (105–4460) | 1257 (156–3615) | 1568 |
| AN-01 | 4000–4500 | 24 | 816 (60–2482) | 644 (129–1544) | 2625 |
| AN-06 | 1050–1450 | 12 | 608 (45–1971) | 487 (75–1170) | 1991 |
| DRA-01 | ≥ 2250 | 13 | 1856 (145–6074) | 1781 (252–4568) | 2078 |
| SEM-07 | 910–1130 | 7 | 1294 (71–4065) | 1365 (150–3710) | 1606 |
| SEM-08 | ≤ 1370 | 19 | 606 (40–1895) | 654 (84–1674) | 1217 |
| CUS-01 | 350–500 | 5 | 821 (60–2590) | 829 (123–2273) | 718 |
| HOL-003 | 2855–3845 | 24 | 1726 (83–5057) | 1339 (291–2920) | 4330 |
| HOL-004 | ≤ 3200 | 7 | 3670 (315–10576) | 2856 (647–6130) | 3979 |
| HOL-007 | 2047.5–2752.5 | 19 | 622 (47–1866) | 495 (102–1140) | 1965 |
| AN-005 | 3250–3650 | 18 | 3366 (264–10072) | 2658 (599–5994) | 3865 |
| AN-007 | 2850–3250 | 10 | 869 (80–2703) | 696 (135–1590) | 1888 |
| AN-010 | ≤ 3200 | 4 | 4178 (224–12699) | 3403 (660–7859) | 1994 |
| AN-011 | ≤ 2896 | 14 | 380 (22–1258) | 300 (44–737) | 656 |
| AN-014 | 747–975 | 10 | 629 (52–1945) | 512 (104–1190) | 949 |
| AN-017 | ≥ 718 | 7 | 844 (56–2725) | 681 (132–1647) | 671 |
| AN-020 | 1850–2250 | 14 | 1542 (116–4907) | 1196 (221–2893) | 3280 |
| AN-028 | ≤ 1350 | 11 | 675 (45–2154) | 541 (108–1239) | 1032 |
| NA-002 | 750–1150 | 5 | 2711 (218–8093) | 3098 (835–6473) | 1022 |
| AU-01 | 1550–5050 | 5 | 2636 (158–8418) | 2534 (278–6892) | 1477 |

| 18 two-sided windows | median \|log error\| | mean signed log error | MAE (years) | RMSE | 90% HPD overlaps window | median inside window |
|---|---|---|---|---|---|---|
| clepsydra, default | 0.30 | −0.30 | 823 | 1,232 | 15 / 18 | 4 / 18 |
| clepsydra, empirical clock | 0.69 | −0.77 | 1,133 | 1,462 | 14 / 18 | 1 / 18 |
| clepsydra, LOFO prior | 0.53 | −0.43 | 1,005 | 1,385 | 15 / 18 | 3 / 18 |
| ASJP (Holman constants) | 0.24 | +0.02 | 541 | 736 | — | 8 / 18 |

The one-sided bounds are met by all clepsydra runs (11 / 11) and by ASJP 7 / 11.
The clepsydra intervals are wide (median 3,450 years against windows a few
hundred years wide), so meeting a bound or overlapping a window says little.
For reference, Rama & Wichmann (2020) report a leave-one-out MAE of 606 years
for ASJP on 30 calibration points.

The error is not uniform. Nodes younger than 2,500 BP (13) are dated 12%
too young on average (mean log error −0.13); the five older ones are dated
about half the age (−0.73; four Austronesian nodes and Wati).

## 2. Family roots against published estimates

| family | published (centre; windows) | clepsydra default | LOFO prior | with attested calibrations | ASJP (this run) |
|---|---|---|---|---|---|
| IE | 6977; 4870–10410 | 2976 (117–11136) | 2639 (122–9179) | 1529 (1112–2102) | 4212 |
| URA | 5113; 3330–7503 | 3233 (147–10847) | 3000 (269–9054) | 4849 (2681–7682) | 3346 |
| TRK | 2066; 1517–2755 | 881 (93–2620) | 843 (132–2063) | — | 3033 |
| SEM | 5750; 4400–7400 | 1377 (76–4391) | 1291 (151–3538) | 1349 (97–4458) | 3335 |
| JPN | 2182; 1239–4190 | 1187 (108–3561) | 1143 (151–3023) | — | 1783 |
| ST | 6957; 4183–11285 | 5525 (345–18607) | 4990 (579–14040) | — | 4174 |
| DRA | 4540; 3000–6500 | 1244 (91–3910) | 1172 (144–3173) | — | 2619 |
| AN | 5173; 4660–5800 | 1218 (104–3648) | 957 (190–2151) | 1752 (1104–2646) | 3349 |
| PN | 5671; 4455–6966 | 17733 (1123–57064) | 16172 (1143–46533) | 7044 (2639–11965) | 4600 |
| UA | 4095; 3258–5025 | 3681 (180–12046) | 3334 (423–9395) | — | 3463 |
| BNT | 4964; 4500–5770 | 722 (35–2467) | 691 (85–1932) | — | 3477 |
| TG | 1762; 847–2740 | 1296 (76–4084) | 1251 (197–3302) | — | 1378 |

| 12 roots | median \|log error\| | mean signed log error | MAE (years) | published window met by median |
|---|---|---|---|---|
| clepsydra, default | 0.85 | −0.70 | 3,191 | 3 / 12 |
| clepsydra, empirical clock | 1.33 | −1.17 | 3,367 | 0 / 12 |
| clepsydra, LOFO prior | 0.97 | −0.78 | 3,231 | 3 / 12 |
| ASJP (Holman constants) | 0.42 | −0.31 | 1,535 | 5 / 12 |

Clepsydra's 90% interval overlaps the published window for 9 of 12 families
and contains all of it for 6, but at a median width of 4,300 years against
published windows of 1,000–7,000.

Read the posterior beside the prior: with no calibration the rate is not
learned (`rate_contraction` between −0.02 and 0.05 in all twelve runs), the
prior-only median of the root is about 2,000 BP, and the posterior is that
prior scaled by how much lexical change the cognates show. Pama-Nyungan (17.7
thousand years) is the lexical age of 30 languages with an automatically
fragmented cognate structure; it is a data failure, not a dating.

**With attested calibrations (last column).** Where the family has attested
internal nodes in the sample (only the ones that come out monophyletic in the
UPGMA tree are kept: Indo-European 2, Uralic 1, Austronesian 1, Pama-Nyungan
1; Semitic none), clepsydra's calendar scale is no longer prior-only.
Uralic (4,849) and Pama-Nyungan (7,044) land inside the published windows.
Indo-European (1,529, interval 1,112–2,103) and Austronesian (1,752) do not
(Semitic had no usable calibration and stays at its blind 1,350): the calibrations are shallow (500–1,250 BP; 750–975 BP),
they fix a fast rate, and the extrapolation to a root five times older fails.
The interval is narrow because the calibrations are, not because the root is
known. This is the time-dependence of the apparent rate (saturation and
replacement of the more variable concepts first); the relaxed clock has no
term for it, and the review's "depth-trend" model would address it.

## 3. Why ASJP is ahead on this test

The data carry less ordering information in clepsydra's input. Across the
18 windowed nodes, the rank correlation between the lexical measure and the
known age is:

| measure | Spearman ρ |
|---|---|
| ASJP LDND (graded phonetic distance) | 0.82 |
| mean GLED cognate distance across the node | 0.51 |
| clepsydra lexical age of the node (rate-free) | 0.30 |

(Roots: clepsydra lexical age against the published centre, ρ = 0.51.) Binary
automatic cognates on 40 concepts saturate: once most concepts have no shared
cognate, deep nodes look alike; LDND keeps graded information on partially
similar words. The model on top (a tree, UPGMA topology, ≤ 24 tips) then
loses some of the raw ordering (0.51 → 0.30), plausibly from the topology and
sampling choices of the benchmark.

## 4. What this does and does not show

- The uncalibrated calendar ages of clepsydra are not competitive with ASJP's
  on GLED, are biased young at depth, and carry intervals too wide to be
  informative. This is consistent with the design (the primary output is the
  lexical age; the calendar scale is a prior scenario, `calendar_basis`
  column) and with the review.
- A one-time calibration of the rate (LOFO) did not help (MAE 1,005 against
  823 for the default). Rates differ by family more than the pooled prior
  admits (sd of the per-family mean log rate 0.5–0.7).
- Calibrations help only when they are as deep as the target.
- Limits of this test: automatic cognates (GLED) rather than expert ones;
  40 concepts; UPGMA topologies; at most 24–30 tips; the dated-node windows
  are partly "max-leaning" historical events; the published family roots
  are Bayesian outputs under their own assumptions; the ASJP side is
  in-sample; 18 and 12 are small samples; the `empirical` prior is not
  independent. The results are not comparable with the expert-cognate LOCO
  of `examples/loco_benchmark/`, which uses attested tips and 8–9 targets.

## 5. Next steps this points to

1. Use graded information from the forms (LDND or sound-class distance within
   a cognate set) as a second likelihood or as a pre-screen of saturated
   nodes: the largest gap found here.
2. A depth-trend (time-dependent rate) term, tested on the same nodes.
3. Expert-cognate data for the same nodes (IE-CoR, ABVD) to separate the data
   effect from the model effect; GLED + a published tree to remove the
   UPGMA effect.

Files: `nodes.tsv` (29 nodes × 3 priors), `roots.tsv` (12 roots × 3 priors,
with ASJP), `roots_calibrated.tsv` (the attested-calibration variant).

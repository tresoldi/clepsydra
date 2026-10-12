# Multi-family held-out-calibration benchmark (LOCO)

The real-world half of leave-family-out validation, across **three families**.
The synthetic study in `tools/validation/leave_family_out.py` certifies the
uncertainty arithmetic; this runs the actual C engine on real cognate data and
asks whether it recovers calibrations it was **not** given.

## Method

Leave-one-calibration-out. Each family's tree carries age calibrations on
attested ancient tips. For each such tip, the harness
(`tools/validation/calibration_benchmark.py`) rewrites the calibration file with
that tip **freed** (fixed age → broad uniform prior, so the sampler tip-dates
it), reruns `date-tree` with the other calibrations still fixed, and scores
whether the withheld age falls in the freed node's inferred 90% HPD. The target
is a value already in the calibration file — nothing is fabricated; the test is
recovery of a datum the model was not shown.

Driven by `manifest.json` (per-family dataset paths + rate priors). Reproduce:

```
python tools/validation/calibration_benchmark.py --binary build/clepsydra \
    --manifest examples/loco_benchmark/manifest.json --out /tmp/loco_multi   # strict
python tools/validation/calibration_benchmark.py --binary build/clepsydra \
    --manifest examples/loco_benchmark/manifest.json --out /tmp/loco_relaxed --relaxed 1
```

About 20 minutes for both clocks on 4 cores. `loco_multi_report.tsv` is the two
reports stacked with a `clock` column, and `web/data/loco.tsv` must be a
byte-identical copy (the web workflow checks it).

## Held-out targets (attested-language dates)

The held-out values are **attestation dates** — when each ancient language is
historically documented — which are established facts, not contested divergence
estimates, and so make sound independent targets. Midpoints (BP) with windows:

| family | tip | date (BP) | attestation window | reference |
|---|---|---|---|---|
| Indo-European | Vedic (Ṛgvedic Sanskrit) | 3200 | 1500–1000 BCE | Witzel 1995; Mayrhofer, *EWAia* |
| Indo-European | Ancient Greek (Classical) | 2400 | 800–400 BCE | standard |
| Indo-European | Latin (Classical) | 2050 | 100 BCE–100 CE | standard |
| Indo-European | Old Irish | 1250 | 600–900 CE | Thurneysen, *GOI* |
| Indo-European | Old Church Slavonic | 1050 | 860–1100 CE | Lunt, *OCS Grammar* |
| Sino-Tibetan | Old Chinese | 2800 | Western Zhou ~1046–771 BCE | Baxter & Sagart 2014 |
| Sino-Tibetan | Old Tibetan | 1300 | 7th–9th c. CE | imperial documents |
| Sino-Tibetan | Old Burmese | 950 | 11th–12th c. CE inscriptions | Myazedi inscription 1113 CE |
| Romance | Latin (Classical) | 2050 | 100 BCE–100 CE | standard |

Cognate data: Indo-European and Romance from IE-CoR (Heggarty et al. 2023);
Sino-Tibetan from Sagart et al. 2019 (via Lexibank). Each family keeps the rate
prior from its standalone example (`../<family>_tree/run.sh`).

## Scoring

Coverage alone rewards wide intervals. Each held-out date is therefore also
scored by the **interval score** of its 90% interval (Gneiting & Raftery 2007):
the width, plus 20 times the distance by which a miss falls outside. Lower is
better, in years. The baseline is the same run's **prior-only** interval for the
freed tip (`prior_*` columns: every other calibration kept, the cognate
likelihood off). **Skill** is 1 − IS(posterior) / IS(prior). It is above 0 when
the cognates improve on what the calibrations and tree prior already say.

## Results (main after e70bf68, 4 chains × 4,000 + 4,000)

Per-target numbers are in `loco_multi_report.tsv`.

| clock | family (N) | covered | mean width | interval score | prior-only IS | skill |
|---|---|---|---|---|---|---|
| strict | Indo-European (5) | 1/5 | 2,185 | 8,483 | 4,358 | −0.95 |
| strict | Sino-Tibetan (3) | 0/3 | 1,082 | 3,762 | 3,824 | 0.02 |
| strict | Romance (1) | 0/1 | 1,258 | 17,078 | 2,102 | −7.1 |
| strict | **all (9)** | **1/9** | 1,714 | **7,864** | 3,930 | **−1.00** |
| relaxed | Indo-European (5) | 4/5 | 3,295 | 4,528 | 4,355 | −0.04 |
| relaxed | Sino-Tibetan (3)* | 3/3 | 2,004 | 2,004 | 3,784 | **0.47** |
| relaxed | Romance (1) | 0/1 | 1,578 | 11,012 | 2,033 | −4.4 |
| relaxed | **all (9)** | **7/9** | 2,674 | **4,407** | 3,907 | **−0.13** |

\* Sino-Tibetan under the relaxed clock mixes slowly. At 4,000 + 4,000 its
worst R̂ was 1.07 (ESS 51). The rows above come from a 10,000 + 10,000 rerun,
where the worst R̂ is 1.024 (Old Chinese, ESS 212). That is above the project's
own 1.01 threshold, so these three rows are indicative only. The shorter run
gave the same picture (3/3, skill 0.43).

Without Romance-Latin, which cannot be recovered (finding 3), the skill over
the eight remaining targets is **+0.13 relaxed** and **−0.61 strict**.

Findings:

1. **The strict clock is overconfident.** Its intervals are the narrowest, and
   they miss 8 of 9 withheld dates, often by more than 1,000 years. A
   forecaster who ignored the cognates and quoted the prior would have scored
   twice as well. Under the 0.2.0 engine the strict clock covered 3 of 9; its
   intervals were a little wider (mean 1,956 BP) and no closer (RMSE 1,354).
2. **The relaxed clock is less overconfident, and only modestly informative.**
   It covers 7 of 9, with smaller point errors (RMSE 1,087 against 1,254 BP).
   On interval score it beats the prior clearly on Sino-Tibetan, ties on
   Indo-European, and is behind overall (skill −0.13 against the prior-only
   interval). Nine targets from three families, with one structurally
   unrecoverable, show it is less overconfident than the strict clock here;
   they do not measure its calibration. Its gain over the strict
   clock comes mostly from widening the intervals, not from moving them.
   Both clocks are wrong in the same direction for the same tips. Vedic comes
   out too old (5,300–5,600 against 3,200 BP), Greek and Old Irish too young.
3. **Romance-Latin cannot be recovered.** Latin is that family's only ancient
   tip. Freeing it leaves no dated tip, so the tree has no tip-age signal at
   all, and the fit falls back to the node windows (500–720 BP inferred against
   2,050). The prior-only interval covers it only because it spans the whole
   allowed range.
4. **Deep tips are the hard ones.** Freeing the first-diverging tip (Vedic, Old
   Chinese) lets its age trade against the root's, and these carry the largest
   errors under both clocks.

What this means for a user: tip dates predicted from cognates alone are
weak evidence, even with four other dated tips on the tree. Under the default
strict clock they are confidently wrong. Report such a date with the relaxed
clock, alongside its prior-only interval, and read a skill near zero as "the
cognates did not add to the calibrations".

## Reading it honestly

This is a **small-N** benchmark — nine held-out points across three families,
one point each for some. It is a calibrated-ness spot check that (a) confirms
the pipeline predicts held-out real dates sensibly and (b) exposes where it does
not, complementing the large-sample synthetic coverage study. It is **not** a
coverage guarantee, and the attestation midpoints carry their own windows (above)
that a single point comparison ignores. Broadening it further is data curation —
supply more trees with trustworthy attested/published calibrations — not a code
change; the manifest format makes adding a family a few lines.

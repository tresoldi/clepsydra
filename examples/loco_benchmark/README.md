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

## Results

Per-target numbers are in `loco_multi_report.tsv` (both clocks). Coverage of the
held-out attested age by the freed node's 90% HPD:

| clock | overall | Indo-European (5) | Sino-Tibetan (3) | Romance (1) | RMSE | MAE |
|---|---|---|---|---|---|---|
| strict  | 0.33 (3/9) | 0.20 | 0.67 | 0.00 | 1354 BP | 1129 BP |
| relaxed | **0.89 (8/9)** | **1.00** | **1.00** | 0.00 | 1164 BP | 905 BP |

Three findings, consistent across three independent families:

1. **The relaxed clock is decisively better calibrated on real data** — overall
   held-out coverage 0.89 vs 0.33, and with lower error (RMSE 1164 vs 1354 BP).
   It is not merely wider intervals: the point errors shrink too (e.g. IE Latin
   +737→+232 BP, Old Irish −884→−451, ST Old Tibetan −207→+81). Real lineages
   are rate-heterogeneous, so the per-branch multipliers earn their keep. This
   is the **mirror image** of the strictly-clocked synthetic check in
   `docs/model-cards.md`, where relaxing does *not* help — same test, opposite
   verdict, decided by held-out prediction on the data at hand, never by
   preference. That the two-clock comparison flips between synthetic and real
   data across three families is the strongest single piece of adequacy evidence
   in the project.

2. **The one persistent miss is structural, not a clock failure.** Romance-Latin
   is missed under *both* clocks (inferred ~750–930 BP vs 2050): Latin is the
   family's *only* attested ancient tip, so freeing it leaves the tree
   heterochronous nowhere — there is no tip-age signal left, and the fit falls
   back to the internal node bounds. You cannot tip-date the sole ancient tip.
   Honest, and exactly what a benchmark should surface.

3. **Deep outgroup tips stay hard.** Vedic (IE, first-diverging) and Old Chinese
   (ST, deep Sinitic) carry the largest errors under both clocks — freeing a deep
   tip lets its age slide toward the wide root prior. The relaxed clock still
   *covers* them (wider, honester intervals) where the strict clock does not.



## Reading it honestly

This is a **small-N** benchmark — nine held-out points across three families,
one point each for some. It is a calibrated-ness spot check that (a) confirms
the pipeline predicts held-out real dates sensibly and (b) exposes where it does
not, complementing the large-sample synthetic coverage study. It is **not** a
coverage guarantee, and the attestation midpoints carry their own windows (above)
that a single point comparison ignores. Broadening it further is data curation —
supply more trees with trustworthy attested/published calibrations — not a code
change; the manifest format makes adding a family a few lines.

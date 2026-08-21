# Leave-one-calibration-out (LOCO) — real held-out-calibration benchmark

This is the **real-world** half of leave-family-out validation (the synthetic
coverage study lives in `tools/validation/leave_family_out.py`). It answers a
question the synthetic study cannot: does the model, dating real cognate data
through the actual C engine, recover a calibration it was **not** given?

## What it does

The five-language IE tree is heterochronous — five attested ancient tips carry
age calibrations (Vedic 3200, Ancient Greek 2400, Latin 2050, Old Irish 1250,
Old Church Slavonic 1050 BP). For each tip in turn, the harness
(`tools/validation/calibration_benchmark.py`) rewrites the calibration file with
that tip **freed** (its fixed age replaced by a broad uniform prior, so the
sampler tip-dates it), reruns `date-tree` with the other four tips + the root
prior still fixed, and checks whether the withheld age falls in the freed node's
inferred 90% HPD. Nothing is fabricated: the target is a value already in
`calibrations.tsv` — the test is recovery of a datum the model wasn't shown.

## Result (`loco_report.tsv`)

| clock | coverage (90% HPD) | RMSE | MAE | mean width |
|---|---|---|---|---|
| strict   | **0.20** (1/5) | 1615 BP | 1421 BP | 2460 BP |
| relaxed  | **0.80** (4/5) | 1478 BP | 1153 BP | 3819 BP |

Two honest findings:

1. **The relaxed clock is materially better calibrated on this real data.** It
   covers 4/5 held-out tips vs the strict clock's 1/5, and by *better point
   estimates*, not merely wider intervals — every tip's error shrank (Latin
   +699→+177, Old Irish −888→−462, Greek −1439→−1088). Real IE lineages are
   **not** strictly clocked, so the per-branch rate multipliers earn their keep.
   This is the mirror image of the synthetic, strictly-clocked check in
   `docs/model-cards.md`, where the relaxed clock does *not* beat the strict one
   out-of-sample — which is exactly the point: flexibility is judged by held-out
   prediction on the data at hand, never by preference.

2. **Deep-tip / root confounding is real and neither clock escapes it.** Vedic,
   the first-diverging (outgroup) tip, is inferred ~6000 BP under both clocks
   (true 3200): freeing the outgroup tip while the root prior spans 4000–10000 BP
   lets its age slide deep, absorbing the unidentified root age. The tips nested
   inside the tree are far better constrained.

Read this as a **small-N (N=5, one family) calibrated-ness spot check**, not a
coverage guarantee — it certifies the pipeline predicts held-out real dates
sensibly and exposes where it does not, complementing the large-sample synthetic
study. Extending it to more families is a data-curation task (supply trees with
trustworthy attested calibrations), not a code change.

## Reproduce

```
python -m validation.calibration_benchmark \
    --binary build/clepsydra --tree examples/five_ie_reproduction/tree.nwk \
    --forms examples/five_ie_reproduction/data/forms.tsv \
    --calibrations examples/five_ie_reproduction/calibrations.tsv \
    --out /tmp/loco --warmup 5000 --samples 5000 --chains 4 \
    --rate-logmean -8.8 --rate-logsd 0.8 [--relaxed 1]
```

(~80 s strict, ~3.5 min relaxed; one MCMC per held-out tip.) The attested tip
ages here are the illustrative attestation dates from the example; swap in
authoritative published ages to turn this into a citable benchmark.

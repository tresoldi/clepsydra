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

These are the Indo-European rows of the three-family benchmark in
`../../loco_benchmark/`, which also explains the scoring. Interval score is in
years, lower is better; skill compares it with the same run's prior-only
interval for the freed tip.

| clock | covered (90% HPD) | RMSE | mean width | interval score | prior-only IS | skill |
|---|---|---|---|---|---|---|
| strict | 1/5 | 1,439 BP | 2,185 BP | 8,483 | 4,358 | −0.95 |
| relaxed | 4/5 | 1,231 BP | 3,295 BP | 4,528 | 4,355 | −0.04 |

The relaxed clock covers four of the five withheld dates where the strict
clock covers one. Neither beats the prior-only interval: the four other dated
tips already say most of what the cognates can say about the fifth. Both clocks
err the same way on the same tips. Vedic, the first-diverging tip, comes out
far too old (about 5,300–5,600 against 3,200 BP), because freeing it lets its
age trade against the root's. Ancient Greek and Old Irish come out too young.

This is five held-out points from one family, a spot check rather than a
coverage guarantee.

## Reproduce

```
python -m validation.calibration_benchmark \
    --binary build/clepsydra --tree examples/five_ie_reproduction/tree.nwk \
    --forms examples/five_ie_reproduction/data/forms.tsv \
    --calibrations examples/five_ie_reproduction/calibrations.tsv \
    --out /tmp/loco --warmup 4000 --samples 4000 --chains 4 \
    --rate-logmean -8.8 --rate-logsd 0.8 [--relaxed 1]
```

(One MCMC per held-out tip.) The attested tip
ages here are the illustrative attestation dates from the example; swap in
authoritative published ages to turn this into a citable benchmark.

# D1.3: an empirical prior for list-level observation error

> 2026-10-11. Code: `tools/error_prior/run.py`. Data:
> `error-prior/pairs.tsv`.

## Method

Two lists of the same language (the same Glottocode) are separated by
essentially no time. The share of their shared concepts on which the cognate
codes differ therefore measures observation error: different collectors,
transcription, glossing and coding. With per-list error probabilities e₁
and e₂, the mismatch is about e₁ + e₂ − e₁e₂.

The scan covers the expert-coded Lexibank datasets, with the same exclusions
as the clock pilot: ancient doculects, degenerate codings, and pairs with
fewer than 40 shared concepts. It finds 24,169 pairs. 22,968 come from one
dialect survey (`bdpa`), so each (dataset, language) group contributes one
value, its median: 103 groups in 26 datasets.

## Result

| statistic | value |
|---|---|
| pair mismatch, median over groups | 0.139 |
| 10% / 25% / 75% / 90% | 0.048 / 0.073 / 0.226 / 0.286 |
| per-list error e (equal errors: 1 − √(1 − m)), mean | 0.088 |
| sd | 0.075 |
| Beta fit | mean 0.088, concentration 13 |

**This is an upper bound on pure observation error.** A shared Glottocode
often joins distinct dialects, whose real divergence is counted here as
error. It is also a prior for *modern* lists made by trained collectors.
Old lists by non-specialists should be expected to be worse.

## Use

`tip_error` rows (`docs/calibrations.md`) take a Beta mean and
concentration. The default recommendation is mean 0.09 and concentration
13 for an ordinary list. Use a wider, higher prior only with a reason
recorded in the row's `source_key`, e.g. a list known to be non-specialist.

On Arawan the 1833 Paumari list reaches a posterior error of 0.23–0.29
under a looser prior (`arawan-acceptance.md`), about three times the
typical modern value.

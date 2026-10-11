# Minimal pairwise divergence (P1)

The smallest example: a two-language divergence date by deterministic pairwise
integration (P1, no MCMC), from **real IE-CoR cognate data** for Spanish and
Portuguese (170 shared meanings).

```
./run.sh          # instant; writes output/posterior_pair.tsv
```

## Representative output

```
n_shared_concepts     170
n_cognate_matches     142
retention_fraction    0.835
divergence_median_bp  ~619
et90_lo_bp / et90_hi_bp   ~342 / ~1112
absolute_anchored     0
prior_dominance       ~0.217
```

Spanish and Portuguese share ~84% of basic vocabulary; the P1 posterior places
their divergence near ~620 BP with a 90% equal-tail interval of roughly
340–1100 BP. `absolute_anchored = 0` reports honestly that, without a tip-age
anchor, the scale rests on the rate prior (`--rate-logmean/-logsd`, per
year). See `../DATA_SOURCE.md` for data licence/attribution.

Swap in any two languages present in `data/forms.tsv` with `--lang-a/--lang-b`,
or point at another example's `data/forms.tsv` for a different pair.

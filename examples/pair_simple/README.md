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
n_cognate_matches     143
retention_fraction    0.841
divergence_median_bp  ~569
et90_lo_bp / et90_hi_bp   ~405 / ~776
absolute_anchored     0
prior_dominance       ~0.10
```

Spanish and Portuguese share ~84% of basic vocabulary; the P1 posterior places
their divergence near ~570 BP with a 90% equal-tail interval of roughly
400–780 BP. `absolute_anchored = 0` reports honestly that, without a tip-age
anchor, the scale rests on the rate prior (`--rate-logmean/-logsd`, per
year). See `../DATA_SOURCE.md` for data licence/attribution.

Swap in any two languages present in `data/forms.tsv` with `--lang-a/--lang-b`,
or point at another example's `data/forms.tsv` for a different pair.

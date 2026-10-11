# Model uncertainty across the sensitivity variants (M4.3)

Every example family with a tree was fitted under each sensitivity variant
that applies to it:

- strict or relaxed clock;
- log-normal or gamma concept-rate mixing;
- for a family whose tips are all modern, the uniform-ratio or the Yule age
  prior.

The variants were compared on held-out concepts (PSIS-LOO) and stacked with
`clepsydra average`. Arawan's four clock × mixing variants are in
`arawan-result.md` §5.

    sh tools/model_uncertainty/run.sh OUTDIR     # ~17 min on 4 cores

Each run uses the example's own tree, forms and calibrations, with 4 chains ×
(4,000 + 4,000) and seed 42, on main at 3216c67.

## Held-out fit

Δ is elpd_loo against the best variant, with the SE of the difference. No
Pareto k exceeds 0.7 in any run, so every LOO estimate is reliable.

| family (concepts) | strict, log-normal | strict, gamma | relaxed, log-normal | relaxed, gamma | weight on relaxed |
|---|---|---|---|---|---|
| Romance (170) | −11.2 ± 4.7 | −10.7 ± 4.4 | −0.4 ± 1.3 | 0 | 1.00 |
| Polynesian, uniform (210) | −29.5 ± 6.7 | −26.5 ± 6.7 | −3.1 ± 1.2 | 0 | 1.00 |
| Polynesian, Yule | −29.7 ± 6.7 | −26.7 ± 6.7 | −3.1 ± 1.2 | −0.4 ± 0.3 | (in the row above) |
| Five-IE (170) | −4.7 ± 3.2 | −4.0 ± 3.1 | −0.8 ± 0.5 | 0 | 0.90 |
| Arawan (253) | −0.3 ± 1.3 | −0.4 ± 1.6 | 0 | −0.2 ± 1.0 | 0.66 |

Three findings:

1. **The relaxed clock predicts better wherever it can be told apart.** The
   margin is 2.4 SE on Romance and 4 SE on Polynesian. Five-IE has five
   tips, and the margin there is only 1.3 SE. On Arawan the variants tie.
   The strict clock never wins.
2. **Gamma mixing is never worse**, and it is better on Polynesian (3.1 ±
   1.2 under the relaxed clock). The heavier slow tail of the concept rates
   is what the data ask for.
3. **The age prior cannot be chosen from the data.** Uniform ratios and Yule
   differ by 0.4 ± 0.3 on Polynesian. As expected, a prior on the ages is
   not testable on held-out concepts.

## What the choice does to the ages

Per-variant posterior medians in BP. The stacked column is the
`average` result.

**Romance** (Latin fixed at 2,050 BP; root and Proto-Romance windows)

| node | strict, log-normal | relaxed, log-normal | relaxed, gamma | stacked (5–95%) |
|---|---|---|---|---|
| root | 2,686 | 2,571 | 2,575 | 2,574 (2,203–2,794) |
| proto_romance | 1,531 | 1,666 | 1,660 | 1,661 (1,412–2,006) |
| eastern_western | 1,024 | 1,155 | 1,153 | 1,153 (829–1,531) |
| ibero | 374 | 366 | 376 | 373 (194–642) |

**Polynesian** (root window 2,600–3,100 BP; all tips modern)

| node | strict, log-normal | relaxed, log-normal | relaxed, gamma | relaxed, gamma, Yule | stacked (5–95%) |
|---|---|---|---|---|---|
| proto_polynesian | 2,837 | 2,837 | 2,835 | 2,833 | 2,835 (2,611–3,083) |
| Tongic | 1,920 | 2,332 | 2,338 | 2,117 | 2,338 (1,325–2,866) |
| eastern_polynesian | 2,233 | 2,044 | 2,052 | 1,955 | 2,052 (1,577–2,476) |
| central_eastern | 1,747 | 1,260 | 1,293 | 1,178 | 1,293 (783–1,753) |

**Five-IE** (five ancient tips; broad root window)

| node | strict, log-normal | relaxed, log-normal | relaxed, gamma | stacked (5–95%) |
|---|---|---|---|---|
| proto_indo_european | 6,060 | 6,538 | 6,378 | 6,340 (4,413–9,533) |
| italo_celtic | 5,316 | 5,241 | 5,126 | 5,150 (3,547–7,704) |

The calibrated nodes hardly move: each is pinned by its own window or tip.
Internal nodes move by as much as the intervals are wide. Central-Eastern
Polynesian goes from 1,750 BP under the strict clock to about 1,290 BP under
the relaxed one, a quarter younger. Proto-Romance moves 130 years older. On
these families, the choice of clock is a larger source of error in the
internal ages than the posterior width under any single clock suggests.

`between_model_share` in `average`'s `nodes.tsv` is 0.000 in all three
families. This is not a sign of agreement. Stacking put all or nearly all
of the weight on one variant, and the share is computed under those
weights. The per-variant medians above are the sensitivity table.

## Convergence

At 4,000 + 4,000 every strict run converged (R-hat ≤ 1.010). Of the 8
relaxed runs, 4 did not:

| run | max R-hat | min bulk ESS |
|---|---|---|
| Polynesian, relaxed, log-normal, uniform | 1.018 | 509 |
| Polynesian, relaxed, log-normal, Yule | 1.014 | 466 |
| Polynesian, relaxed, gamma, uniform | 1.011 | 459 |
| Five-IE, relaxed, log-normal | 1.025 | 253 |

These misses are small: their medians agree with the converged neighbouring
variants. A result for publication would need longer runs. The relaxed clock
still mixes slowly on 8–9 tips: its branch multipliers and the ages trade
off, which is the confounding that `branch_rate_confounded` flags.

## Consequences

- **The default stays the strict clock for now, with the relaxed clock as
  the standard sensitivity run.** The relaxed clock predicts better on two
  of four families and never worse. As the default, it would have to
  converge at default run lengths, which it does not yet do reliably on 8–9
  tips, and its branch rates are confounded with ages wherever calibrations
  are sparse. This is recorded as an open item: improve relaxed mixing, then
  revisit the default.
- **Gamma concept rates** are a candidate default on the same terms: never
  worse, sometimes better. Their only cost is the heavier tail of the rate
  prior.
- **Report the clock comparison with every tree result.** Use `date-tree
  --relaxed 0/1 --write-pointwise 1` and then `average`, and give the
  per-model medians of the internal nodes.

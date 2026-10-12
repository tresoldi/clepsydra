# Simulated worlds: coverage when the model is wrong (D3)

`tools/simworld/run.py` simulates families of 6–9 modern languages and
150 concepts. Each replicate draws a random topology with uniform node
ages, a root between 1,000 and 4,000 years, a clock rate around 2e-4 per
year, and log-normal concept multipliers (σ 0.8, mean 1). It then breaks
one assumption at a time:

| scenario | what is wrong |
|---|---|
| clean | nothing |
| borrowing | 10% of concepts copied from one tip to another after the split |
| coding | every list replaces 8% of its forms in error |
| time_dep | the rate falls with depth, r(t) = r0 / (1 + t/1000) |
| missing_fast | fast concepts are missing more often (up to 60%) |
| gamma_rates | gamma concept multipliers, CV 1.2 |

Each replicate is fitted with the default `date-tree` (4 chains × 2,000 +
2,000, no calendar information). The score is whether the 90% HPD of each
internal node's lexical age contains the true one (rate × age). There are
20 replicates per scenario and 125–137 nodes. Nodes from one tree are
correlated, so the standard error of a coverage is closer to 0.05 than to
the binomial 0.026.

    cd tools && python3 -m simworld.run --reps 20 --out ../docs/research/simworld

`nodes.tsv` holds the results of the default model, and `nodes_tip_error.tsv`
those of the variant in §2. Per-replicate inputs are written under `runs/`,
which is not committed.

## 1. Default model (2026-10-11, main after 19197ca)

| scenario | 90% coverage | root coverage | mean log bias | mean log width | root runs converged |
|---|---|---|---|---|---|
| clean | 0.81 | 0.85 | −0.06 | 0.61 | 1.00 |
| borrowing | 0.88 | 0.85 | +0.05 | 0.81 | 1.00 |
| coding | **0.32** | 0.60 | **+0.45** | 0.40 | 1.00 |
| time_dep | **0.31** | **0.05** | **−0.38** | 0.65 | 1.00 |
| missing_fast | 0.81 | 0.85 | −0.09 | 0.76 | 1.00 |
| gamma_rates | 0.90 | 0.85 | +0.13 | 0.75 | 0.90 |

With the empirical concept-rate prior (`concept-rate-prior.md`) replacing
the HalfNormal, clean coverage moved from 0.78 to 0.81 and its bias stayed
at −0.06. The remaining shortfall is about 1.5–2 standard errors. The
clean world is not the fitted prior exactly: its topologies come from
random pairing, not from the uniform-ratio prior. The shortfall is
therefore not, on its own, evidence of a bug. Calibration in the strict
sense is tested by the SBC suite, where the generator is the prior.

Borrowing at this rate, non-random missingness and a wrong rate-mixing
family are absorbed: coverage holds within noise. Two violations are not:

- **Coding error** adds terminal change that the model reads as time.
  Every node is too old, by +0.45 in the log on average. The bias is worst
  at shallow nodes (true lexical age < 0.15), where coverage is 0 and
  the bias is +1.0: their age is mostly error.
- **A rate that falls with depth** makes deep nodes look young. The bias is
  −0.38 and the root is covered once in 20. No model with a constant clock
  can detect this from modern lists alone. Only calibrations at more than
  one depth can.

## 2. A list-error term on every list

The obvious fix for coding error is the D1.3 list-error prior (Beta, mean
0.088, concentration 13) on every list, not only on historical ones. Same
seeds, same data:

| scenario | 90% coverage | mean log bias | log width | converged |
|---|---|---|---|---|
| clean | 0.81 → **0.56** | −0.06 → **−0.33** | 0.61 → 1.38 | 1.00 → 0.90 |
| coding | 0.32 → 0.85 | +0.45 → +0.14 | 0.40 → 1.06 | 1.00 → 0.55 |
| borrowing | 0.88 → 0.78 | +0.05 → −0.18 | 0.81 → 1.44 | 1.00 → 0.75 |
| gamma_rates | 0.90 → 0.87 | +0.13 → −0.02 | 0.75 → 1.34 | 0.90 → 0.75 |

It repairs coding error and breaks the clean world. On clean data the
shallow nodes shrink by half (log bias −0.72, coverage 0.44). The deep nodes
are also too young (−0.21). The reason is identifiability. When every
list is modern, a list's error probability and the length of its terminal
branch have the same effect on the data. The prior's mean of about 9% then
takes real terminal change away from the branches. The confound also
slows mixing: convergence falls to 55–90%.

**Decision: per-list error stays opt-in** (a `tip_error` calibration row),
recommended where a list is known to be special, such as a historical,
non-specialist or secondary-source record. That is where it has
independent support. For Arawan's 1833 list, the date separates it from
modern Paumari. The error term cannot be the default until something other
than the prior separates the two. Candidates: same-language list pairs
(D1.3's data used inside the fit), or an error prior concentrated near 0
with a long tail.

## 3. What this means for users

- In the scenarios tried, lexical ages held up when the lists were clean:
  even then coverage was 0.81 against 0.90, and whether coding is "careful"
  cannot be measured for most historical data. Coding noise of a few
  percent per list ages shallow nodes markedly. This is the effect that the
  adequacy panel's singleton share and incompatible share are there to
  flag.
- Calendar ages from modern-only data rest on a constant clock. If rates
  decay with depth, deep nodes will be too young and their intervals will
  not show it.

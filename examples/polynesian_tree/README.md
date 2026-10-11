# Polynesian chronology (a shallow, cleanly-dated family, T1)

Dates a fixed Polynesian topology from **real cognate data** — Walworth (2018) —
for eight doculects across the main Polynesian branches (Tongic, Samoic, Eastern
Polynesian). Polynesian is a *shallow, well-dated* family, so unlike the deep
Indo-European and Sino-Tibetan roots this one converges crisply and its ages are
archaeologically defensible — the clean-case counterpoint in the example set.

```
./run.sh          # ~2 min; writes output/posterior_nodes.tsv, diagnostics.tsv, run_manifest.tsv
```

## Inputs

- `data/` — a LexiCore subset from `lexibank/walworthpolynesian` (CC-BY-4.0; see
  `../DATA_SOURCE.md`), via `tools/datasets/lexibank.py`.
- `tree.nwk` — an **illustrative** fixed topology following the accepted
  Polynesian subgrouping (Tongic vs Nuclear; Samoic and Eastern Polynesian
  within Nuclear; Rapanui basal in Eastern).
- `calibrations.tsv` — the Proto-Polynesian root is anchored by a **node
  calibration** (uniform 2600–3100 BP, the archaeological settlement window),
  rather than by tip-dating a reconstruction; extant tips are at 0. This is the
  clean way to date a shallow family whose members are all modern.

## Representative output (seed 42, 4 chains × 10000)

| node | median BP | 90% HPD | R̂ | ESS |
|---|---|---|---|---|
| proto_polynesian (root) | ~2833 | 2598–3070 | 1.00 | 14164 |
| nuclear_polynesian | ~2778 | 2552–3044 | 1.00 | 12468 |
| Samoic | ~2526 | 2116–2931 | 1.00 | 7511 |
| eastern_polynesian | ~2231 | 1945–2527 | 1.00 | 12303 |
| Tongic | ~1919 | 1380–2435 | 1.00 | 9593 |
| central_eastern (Māori–Tahitian–Hawaiian core) | ~1742 | 1448–2042 | 1.00 | 9443 |

Every node converges (R̂ ≤ 1.001, ESS above 7,500). At 5,000 + 5,000 a few
parameters stopped just short of R̂ 1.01, hence the longer run and the ages line up with the
linguistic/archaeological picture: a Proto-Polynesian break-up ~2900 BP, an
Eastern-Polynesian node in the low thousands, and the Central-Eastern core more
recent. Absolute ages still rest on the root calibration and rate prior, as
`run_manifest.tsv` records — see `docs/scientific-methods.md`.

Contrast with `examples/sinotibetan_tree` and `examples/five_ie_reproduction`,
where the deep root is only weakly identified: the difference is time depth and
calibration, not the tool.

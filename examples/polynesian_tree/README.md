# Polynesian chronology (a shallow, cleanly-dated family, T1)

Dates a fixed Polynesian topology from **real cognate data** — Walworth (2018) —
for eight doculects across the main Polynesian branches (Tongic, Samoic, Eastern
Polynesian). Polynesian is a *shallow, well-dated* family, so unlike the deep
Indo-European and Sino-Tibetan roots this one converges crisply and its ages are
archaeologically defensible — the clean-case counterpoint in the example set.

```
./run.sh          # ~45 s; writes output/posterior_nodes.tsv, diagnostics.tsv, run_manifest.tsv
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

## Representative output (seed 42, 4 chains × 5000)

| node | median BP | 90% HPD | R̂ | ESS |
|---|---|---|---|---|
| proto_polynesian (root) | ~2908 | 2675–3100 | 1.00 | 1601 |
| nuclear_polynesian | ~2849 | 2604–3067 | 1.00 | 1468 |
| eastern_polynesian | ~2250 | 1987–2558 | 1.00 | 1393 |
| central_eastern (Māori–Tahitian–Hawaiian core) | ~1752 | 1450–2050 | 1.01 | 963 |

Every node converges (R̂ ≈ 1.00, ESS ~1000–3900) and the ages line up with the
linguistic/archaeological picture: a Proto-Polynesian break-up ~2900 BP, an
Eastern-Polynesian node in the low thousands, and the Central-Eastern core more
recent. Absolute ages still rest on the root calibration and rate prior, as
`run_manifest.tsv` records — see `docs/scientific-methods.md`.

Contrast with `examples/sinotibetan_tree` and `examples/five_ie_reproduction`,
where the deep root is only weakly identified: the difference is time depth and
calibration, not the tool.

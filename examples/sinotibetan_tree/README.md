# Sino-Tibetan chronology (a different family, T1)

Dates a fixed Sino-Tibetan topology from **real cognate data** — Sagart et al.
(2019), the Sino-Tibetan Database of Lexical Cognates behind the PNAS dated
phylogeny — using eight doculects across five subgroups (Sinitic, Tibetic,
Burmish, Kiranti, rGyalrong), with **Old Chinese, Old Tibetan and Old Burmese**
as attested tip anchors. This complements the Indo-European examples with a
second, unrelated family, driven by the same converter and CLI.

```
./run.sh          # ~60-80 s; writes output/posterior_nodes.tsv, diagnostics.tsv, run_manifest.tsv
```

## Inputs

- `data/` — a LexiCore subset from `lexibank/sagartst` (CC-BY-4.0; see
  `../DATA_SOURCE.md`). Generated with the family-agnostic converter
  `tools/datasets/lexibank.py`.
- `tree.nwk` — an **illustrative** fixed topology grouping the doculects by
  subgroup. Not a phylogenetic claim; Sagart et al. infer the topology.
- `calibrations.tsv` — attested tip ages (Old Chinese 2800, Old Tibetan 1300,
  Old Burmese 950 BP; illustrative), extant tips at 0, and a broad uniform prior
  on the Proto-Sino-Tibetan root.

## Representative output (seed 42, 4 chains × 5000)

| node | median BP | 90% HPD | R̂ | ESS |
|---|---|---|---|---|
| Sinitic (Old Chinese + Mandarin) | ~2979 | 2816–3178 | 1.00 | 3914 |
| Tibetic (Old Tibetan + Lhasa) | ~1407 | 1306–1521 | 1.00 | 3848 |
| Burmish (Old Burmese + Rangoon) | ~1052 | 955–1152 | 1.00 | 3270 |
| Tibeto-Burman core | ~5771 | 5158–6357 | 1.00 | 616 |
| proto_sino_tibetan (root) | ~8793 | 8391–9000* | 1.00 | 2422 |

The **subgroup ages are well estimated** and archaeologically sensible (Sinitic
diversification ~3000 BP; Tibetic ~1400 BP). The **root is not**: under a strict
clock the rate pinned by the recent tip anchors implies a very deep root that
pushes against whatever upper bound the prior sets (`*` the 90% HPD hugs the
9000 BP ceiling). That is an honest identifiability result, not a fitted date —
strict-clock, tip-age-only dating poorly constrains a deep family's root, and
Sagart et al.'s ~7200 BP comes from a **relaxed clock with node calibrations**.
Treat the root here as "deep and weakly identified", and see
`docs/scientific-methods.md`. Adding `--relaxed 1` and internal node
calibrations is the route to a comparable root estimate.

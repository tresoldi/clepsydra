# Sino-Tibetan chronology (a different family, T1)

Dates a fixed Sino-Tibetan topology from **real cognate data** — Sagart et al.
(2019), the Sino-Tibetan Database of Lexical Cognates behind the PNAS dated
phylogeny — using eight doculects across five subgroups (Sinitic, Tibetic,
Burmish, Kiranti, rGyalrong), with **Old Chinese, Old Tibetan and Old Burmese**
as attested tip anchors. This complements the Indo-European examples with a
second, unrelated family, driven by the same converter and CLI.

```
./run.sh          # ~2 min; writes output/posterior_nodes.tsv, diagnostics.tsv, run_manifest.tsv
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
| Sinitic (Old Chinese + Mandarin) | ~2941 | 2810–3096 | 1.00 | 4474 |
| Tibetic (Old Tibetan + Lhasa) | ~1385 | 1306–1486 | 1.00 | 5371 |
| Burmish (Old Burmese + Rangoon) | ~1040 | 964–1129 | 1.00 | 5599 |
| Tibeto-Burman core | ~6483 | 5773–7238 | 1.00 | 1889 |
| proto_sino_tibetan (root) | ~8937 | 8524–9209* | 1.00 | 3750 |

The **subgroup ages are well estimated** and archaeologically sensible (Sinitic
diversification ~3000 BP; Tibetic ~1400 BP). The **root is not**: under a strict
clock the rate pinned by the recent tip anchors implies a very deep root that
pushes against whatever upper bound the calibration sets (`*` the 90% HPD runs
past the soft 4500–9000 BP root window, into its outlier tail). That is an honest identifiability result, not a fitted date —
strict-clock, tip-age-only dating poorly constrains a deep family's root, and
Sagart et al.'s ~7200 BP comes from a **relaxed clock with node calibrations**.
Treat the root here as "deep and weakly identified", and see
`docs/scientific-methods.md`. Adding `--relaxed 1` and internal node
calibrations is the route to a comparable root estimate.

## Prior only against posterior

![Sino-Tibetan: 90% HPD of each node under the prior alone and under the posterior](prior_overlay.svg)

The grey bar is the node's age with the cognate likelihood switched off (the same run's `prior_*` columns: tree prior, calibrations and clock prior only). The blue bar is the posterior. `prior_overlap` in `posterior_nodes.tsv` measures how much the two share.

The data dominate every node (overlaps 0.09–0.35). The subgroups with an
attested tip collapse onto a narrow interval just above that tip. The root is
pushed to 8,937 BP, past the prior's 90% HPD (4,368–8,437). This is a conflict
between the cognates, under a strict clock, and the root window, not
information about the root: the posterior sits in the window's outlier tail.

Regenerate after `./run.sh` with
`python3 ../../tools/plotting/overlay.py output/posterior_nodes.tsv --svg prior_overlay.svg --title "Sino-Tibetan"`, which also prints the numbers as a table.

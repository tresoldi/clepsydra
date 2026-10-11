# Romance chronology (fixed tree, T1)

Dates the internal nodes of a fixed Romance topology from **real IE-CoR cognate
data** (8 doculects, 170 meanings) under the strict novel-state model (T1), with
Classical Latin as an attested tip anchoring the absolute scale.

```
./run.sh          # ~30 s; writes output/posterior_nodes.tsv, diagnostics.tsv, run_manifest.tsv
```

## Inputs

- `data/` — IE-CoR subset (Latin, Italian, Spanish, Portuguese, French, Romanian,
  Catalan, Sardinian). See `../DATA_SOURCE.md` for licence/attribution and how to
  regenerate. CC-BY-4.0 (Heggarty, Anderson & Scarborough, IE-CoR v1.0).
- `tree.nwk` — an **illustrative** fixed Romance topology (Sardinian basal,
  Eastern/Western split, Gallo-Ibero grouping). Not a phylogenetic claim.
- `calibrations.tsv` — Latin fixed at 2050 BP (illustrative Classical-Latin
  attestation), extant tips at 0, broad uniform priors on the root and
  proto-Romance node. Illustrative anchors, not authoritative dates.

## Representative output (seed 42, 4 chains × 4000)

| node | median BP | 90% HPD | R̂ | ESS |
|---|---|---|---|---|
| root (Latin era) | ~2686 | 2480–2825 | 1.00 | 6040 |
| proto_romance | ~1531 | 1377–1761 | 1.00 | 3585 |
| continental | ~1365 | 1156–1631 | 1.00 | 3452 |
| gallo_iberian | ~647 | 472–842 | 1.00 | 4231 |
| ibero (Spanish–Portuguese) | ~374 | 232–523 | 1.00 | 4057 |

Proto-Romance lands near ~1500 BP (Late Latin dispersal). Every node reports a
posterior **distribution** with a credible interval and convergence diagnostics
(R̂, ESS) — never a point date. Absolute ages depend on the rate prior and the
calibrations, as `run_manifest.tsv` records.

Try `--relaxed 1` to add per-branch rate multipliers, and compare the branch-
rate/date confounding diagnostic in `run_manifest.tsv`.

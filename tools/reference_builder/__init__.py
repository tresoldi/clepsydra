"""
Clepsydra reference-prior builder.

A dependency-free (Python stdlib only) toolchain for assembling versioned
reference prior packages consumed by the C inference engine as strict TSV.

Sources (initial evidence):
  - Pagel, Atkinson & Meade 2007 — ~200 IE per-meaning relative rates.
    SUPPLEMENT MISSING: the per-meaning rate table (Supp Table 2) is not in
    our literature/. A synthetic generator stands in until supplied.
  - Vejdemo & Hörberg 2016 — 167 content-concept predictive predictors.
    SUPPORTING INFO MISSING. A synthetic generator stands in.
  - ASJP (pinned release) — worldwide relative concept stability.
    DATABASE MISSING. A synthetic generator stands in.

The builder emits:
  - concept_priors.tsv  (per docs/input-formats.md §5)
  - family_priors.tsv   (per docs/input-formats.md §6)
  - manifest.json       (version, sources, checksums, provenance)
  - leave-family-out variants

Usage:
  python -m reference_builder build --out reference/generated/synthetic_v1
  python -m reference_builder build --out ... --leave-family-out indo_european
"""
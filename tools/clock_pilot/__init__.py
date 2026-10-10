"""D2 σ_f pilot: how much does the lexical clock rate vary between families?

Roadmap item D2 (`docs/roadmap.md`). The pre-registration is
`docs/research/clock-pilot/plan.md`; the report is
`docs/research/clock-pilot.md`.

Usage (deterministic; ~minutes, most of it the optional LexStat run)::

    cd tools
    python3 -m clock_pilot.run              # everything
    python3 -m clock_pilot.run --no-lexstat # skip the uniform automatic coding
    python3 -m clock_pilot.run --no-abvd    # skip the ABVD sensitivity

Inputs (all local, read-only):

- `reference/calibrations/calibrations.tsv` plus `clock_pilot/curation.tsv`
  (per-row decisions: age treatment, Glottolog node, exclusions);
- the 59 Lexibank datasets flagged ExpertCognates in
  `~/work/chl/arcaverborum/datasets.csv`, under
  `~/work/chl/arcaverborum/raw/lexibank/<name>/cldf/` (and `abvd`, for a
  sensitivity analysis only);
- Glottolog CLDF v5.0 as a zip (`--glottolog`; SHA-256 checked).

Outputs, in `docs/research/clock-pilot/`: `calibration_matches.tsv`,
`results.tsv`, `lofo.tsv`, `lofo_one_sided.tsv` (and `lexstat_matches.tsv`
when lingpy is installed).

Modules: `glottolog` (lineages), `calibrations` (curation, node identity,
deduplication), `lexibank` (coding rule), `match` (cross pairs and counts),
`model` (distances, bootstrap, REML/Bayesian random effects, censored ML),
`lofo` (predictors, baselines, scores), `uniform` (LexStat recoding),
`run` (orchestration). Dependencies: Python stdlib, numpy, scipy; lingpy
optional.
"""

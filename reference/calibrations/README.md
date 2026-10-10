# Calibration reference

| file | what it is |
|---|---|
| `oldworld.tsv`, `newworld_pacific.tsv` | the two regional compilations, with their notes (`*.md`) |
| `calibrations.tsv` | both merged into one schema (`tools/reference_builder/merge_calibrations.py`): 367 rows, graded A/B/C/U |
| `curated_nodes.tsv` | the grade A/B node-type rows prepared for fitting (roadmap D1.2) |

## `curated_nodes.tsv`

Each row of `calibrations.tsv` that dates a node gets a decision. The
decisions were made by the D2 clock pilot before any rate was computed
(`docs/research/clock-pilot/plan.md`). The columns:

- `use`: `interval` (two-sided, `lo_bp`–`hi_bp`), `point` (`center_bp`),
  `max` or `min` (one-sided), or `exclude`. Many published "points" are
  political or historical events that bound a split from one side only.
  They are recorded as one-sided, or as an interval when the source gives an
  explicit range.
- `node_spec`: the node in Glottolog terms. `crown:X` is the crown of
  clade X; `split:A|B` is the split between A and B; `stem:X` is X's stem.
- `glottolog_mrca` / `node_name`: the node's Glottolog identity (v5.0).
- `status`: `primary`, `duplicate of …` (the same Glottolog node as an
  earlier row), or `excluded: …` with the reason.

Duplicates are identified by node, not by row id. Two compilations often
date the same split (Ma'anyan–Malagasy, Cham, Inuit, Proto-Oceanic, …), and
counting both would double their weight.

Changes to `calibrations.tsv` itself are recorded in its `notes` column.
On 2026-10-11 two rows were demoted from B to C because their dates are
lexical or Bayesian dating results, which makes them circular as
calibrations: NA-020 (Proto-Uto-Aztecan) and AU-02 (Pama-Nyungan).

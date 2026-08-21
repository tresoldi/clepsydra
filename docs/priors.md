# Reference priors

> Stub (Milestone 0). The reference-prior **builder** is a Python package
> (`tools/reference_builder/`, M3) that emits the strict-TSV prior files the
> C engine loads (`concept_priors.tsv`, `family_priors.tsv`; schemas in
> `docs/input-formats.md`). The C engine **never** hardcodes concept or
> family rates.

## Reference modes

`reference_mode ∈ {uninformative, global, family_specific,
leave_current_family_out}`. Default `family_specific` for prediction;
`leave_current_family_out` for scientific validation (so a family is never
tested with priors trained partly from itself). Every output records the
selected mode.

## Evidence sources (initial)

| Source | Use | Status |
|---|---|---|
| Pagel, Atkinson & Meade 2007 | ~200 IE per-meaning relative rates; centered in log space, stored separately from the clock scale. Not universal absolute rates. | **table missing** (Supp. Table 2 not in PDF) — stub |
| Vejdemo & Hörberg 2016 | 167 content-concept predictive priors (frequency, synonyms, senses, imageability); CV-fit; content vs function separate. | **supporting data missing** — stub |
| ASJP (pinned release) | worldwide relative concept stability; family-balanced; no direct absolute per-year rates from lexical distance. | **database missing** — stub |
| Rama & Wichmann 2020 | 30 calibration groups as benchmark / leave-group-out resource. **Not** our clock. | groups captured in `docs/paper-evaluation.md` |

## Legacy overall-scale priors (fallback anchors)

| prior_id | unit | log_mean | log_sd | description |
|---|---|---|---|---|
| `legacy_swadesh100` | per_kyr | −1.8916490462 | 0.80 | 0.86 retention/kyr |
| `legacy_lees200` | per_kyr | −1.5271140564 | 0.80 | 0.8048 retention/kyr |

Wide SD intentional; fallback anchors, not empirical conclusions.

## Provisional default hyperpriors

See `docs/model-spec.md` §2. Sensitivity presets `tight`/`standard`/`wide`
must all be run for any published result; no published result relies on one
preset.

## Evidence grades

`A` open matrix+topology+calibrations+posterior samples; `B` open
matrix+calibration spec, incomplete posterior; `C` published posterior
summaries only; `D` traditional/approximate. Only A/B train detailed
concept-rate priors by default; C is external validation; D is
qualitative/sensitivity.

## Published family benchmarks (registry, separate from inference priors)

A benchmark row contains: `family_id, study_id, root_definition,
posterior_median, posterior_lower, posterior_upper, interval_probability,
calibration_summary, data_availability, evidence_grade, citation`. Do not
create a numeric row unless the exact root definition and uncertainty
interval are confirmed. Published root ages are benchmark observations,
**not** automatically valid priors; the same family may have several
conflicting rows. (Examples: Indo-European, 5-language Core-IE case study,
Dravidian, Austronesian, Bantu, Sino-Tibetan, Turkic, Uralic, Semitic,
Japonic, selected South American families.)

## Data-driven concept stability (opt-in, in-sample)

`tools/analysis/concept_stability.py` estimates each concept's stability
directly from a dataset as the mean pairwise cognate agreement across the
languages that share it (the `distglot` method), classifies concepts
stable/regular by a threshold, and can emit an empirical `concept_priors.tsv`
(centered relative log-rates: more stable ⇒ lower rate). It is deliberately
**not a default and not the reference package**: the reference builder
(`tools/reference_builder`) draws on independent published sources, whereas this
estimator reads the very data one would then date.

> ⚠️ **Circularity.** The single-`--forms` estimate is *in-sample*. Using it to
> date the same dataset double-counts the evidence and understates uncertainty.
> Those rows are stamped `evidence_grade=D`, `source_key=EMPIRICAL_INSAMPLE`.
> Legitimate uses: exploratory concept ranking; prior-sensitivity checks; and the
> leave-family-out mode below. It is an option because those uses are real, not
> because in-sample priors are sound.

### Leave-family-out (out-of-sample)

`--leave-out <family>` is the disciplined mode: it pools stability across *other*
families and emits a prior for a held-out target family it never saw — the same
discipline as `tools/validation/leave_family_out.py`. Concepts are aligned across
families by their **Concepticon id** (read from each family's `concepts.tsv`), so
a concept lines up even though every family uses a different local `concept_id`
(`C100_yellow` in one, `C57_yellow` in another). Pooling sums cognate agreements
and language pairs across families (**pair-count-weighted**, so a family that
shares a concept over a single language pair — a 0/1 draw — does not count the
same as one with hundreds of pairs). The emitted rows are keyed back to the
*target* family's own local `concept_id`s (each target concept gets its own row,
even two that share a Concepticon id), so the file drops straight into dating that
family, and are stamped `evidence_grade=D`, `source_key=EMPIRICAL_LFO` (still
low-grade empirical, but now legitimately out-of-sample for the target).

> ⚠️ **Depth confounding.** The relative-log-rate transform assumes a single
> divergence depth per concept (`s ≈ exp(-r_c·d_avg)`, centering removes
> `d_avg`). Pooling *raw similarities* across families of different time depths
> violates that, so a concept seen mostly in shallow vs deep families gets a
> depth-biased rate, and agreements among closely-related languages are treated
> as independent when they are not. Leave-family-out priors are a grade-D
> relative **ranking** of concept stability, not a depth-calibrated rate — good
> for exploratory ranking and sensitivity checks, not a load-bearing prior for a
> published date.

```
python -m analysis.concept_stability \
    --families-dir examples/families_pairwise --leave-out dravidian \
    --out pooled_stability.tsv --emit-priors dravidian_priors.tsv
```

Each subdirectory of `--families-dir` with a `forms.tsv` is one family (a
`concepts.tsv` alongside it supplies the Concepticon alignment; without it the
local `concept_id` is used as the key). Families can also be listed explicitly
with repeatable `--family label=forms.tsv[,concepts.tsv]`. The `--out` table
reports pooled per-concept stability with its family count; `--emit-priors`
writes the target-excluded `concept_priors.tsv`. Coverage (how many of the
target's concepts the pool could inform) is printed so partial overlap is never
silent.

## Open items (tracked in `docs/paper-evaluation.md`)

M3 cannot complete the real-data prior build until: Pagel 2007 Supp. Table 2,
Vejdemo–Hörberg 2016 Supporting Information, and a pinned ASJP release are
supplied. Synthetic fixtures unblock M1/M2/M4–M8 in the meantime.

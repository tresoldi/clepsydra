# Input formats

> Milestone 0 deliverable. Final TSV schemas and Newick tree format for
> Clepsydra. The C engine contains a small dependency-free parser for all of
> these (M1). Parser rules (brief §1.6):
> - input from a filesystem path **or** an in-memory buffer `(const char*,
>   size_t)` (Wasm has no filesystem);
> - headers mandatory;
> - fields separated by literal `\t`;
> - UTF-8 stored/passed through, but the engine does not normalize Unicode;
> - `.` = missing;
> - blank lines ignored;
> - lines beginning with `#` are comments;
> - malformed numeric fields produce useful errors;
> - duplicate IDs detected;
> - unknown **required** columns are errors; unknown **optional** columns may
>   be retained or ignored with a warning;
> - no quoted/multiline-field syntax;
> - **never** `exit()`/`abort()` inside the reusable library.

All files are strict TSV (tab-separated). Numeric years are "BP" relative to
`time_origin_year` (default 1950) unless stated.

## 1. `languages.tsv`

```
language_id	name	family_id	is_extinct	metadata
latin	Latin	indo_european	1	.
italian	Italian	indo_european	0	.
```

- `language_id` — unique, stable ASCII identifier.
- `family_id` — used to key `family_priors.tsv` and the `β_f` hierarchy.
- `is_extinct` — `0`/`1`.
- `metadata` — optional opaque JSON-ish or `.`; retained, not parsed by C.

## 2. `concepts.tsv`

```
concept_id	label	external_namespace	external_id	word_class	status
C000948	WATER	Concepticon	948	noun	active
C_CUSTOM_001	CUSTOM_CONCEPT	.	.	unknown	active
```

- `concept_id` — canonical, immutable (integer or ASCII). The inference code
  uses **only** this ID.
- `external_namespace` / `external_id` — e.g. Concepticon; separate table,
  not an architectural dependency. `.` if none.
- `word_class` — noun/verb/adj/…/unknown (used to separate content vs
  function concepts in priors).
- `status` — `active`/`unknown`/`deprecated`. Unknown concepts receive broad
  fallback priors. **Never** match concepts silently by approximate label.

## 3. `forms.tsv`

```
analysis_id	language_id	concept_id	form_id	form	cognate_id	primary_weight	cognacy_confidence	inheritance_probability	loan_probability	status
expert_a	italian	C000948	f001	acqua	COG_WATER_1	1.0	1.0	0.99	0.01	observed
expert_a	english	C_SMALL	f002	small	COG_SMALL_1	0.75	0.95	0.98	0.02	observed
expert_a	english	C_SMALL	f003	little	COG_SMALL_2	0.25	0.90	0.95	0.05	observed
```

- `analysis_id` — names an alternative cognacy coding; the CLI runs them
  separately or samples among complete codings (§5.4 of `model-spec.md`).
- `cognate_id` — cognate class label (arbitrary unique ID; only equality
  structure matters for T1).
- `primary_weight` — within-cell weight; non-missing weights must sum to ~1
  (validated; if all missing, set `1/n_candidates` with a warning).
- `cognacy_confidence` — `y_c` for pairwise soft-observation likelihood.
- `inheritance_probability`, `loan_probability` — feed loan modes
  (`exclude`/`mixture`/`sensitivity`).
- `status` — `observed`/`missing`/`excluded`/`absent`. Missing means unknown.

Multiple forms per (analysis, language, concept) cell are retained (synonyms
/ polymorphism). A missing cell has no rows (or rows with `status=missing`).

## 4. `calibrations.tsv`

```
calibration_id	target_id	type	p1	p2	p3	is_hard	source_key
latin_age	latin	fixed_tip	1900	.	.	1	HIST_001
root_prior	romance_root	offset_lognormal	3000	7.5	0.4	0	HIST_002
```

See `docs/calibrations.md` for the exact meaning of `p1/p2/p3` per `type`,
feasibility checking, and sensitivity procedures. Each calibration names a
labelled node or tip. `is_hard` `1` ⇒ hard bound; `0` ⇒ soft.

## 5. `concept_priors.tsv`

```
concept_id	scope_id	unit	log_rate_mean	log_rate_sd	n_families	n_observations	evidence_grade	source_key	version
C000948	GLOBAL	per_kyr	-2.15	0.55	20	400	B	META_001	1.0.0
C000948	indo_european	per_kyr	-2.32	0.32	1	87	A	PAGEL_2007	1.0.0
```

- `scope_id` — `GLOBAL` or a `family_id` (for `family_specific` /
  `leave_current_family_out` modes).
- `unit` — `per_kyr` (per millennium) is canonical; the loader converts.
- `evidence_grade` — `A`/`B`/`C`/`D` (§9.6 of the brief). Only A/B train
  detailed concept-rate priors by default; C is external validation; D is
  qualitative/sensitivity.
- `version` — prior-package version (recorded in the run manifest).

These are the reference priors the C engine loads (never hardcoded). Family
data *updates* them via the hierarchy (`μ_c ~ N(m_c^ref, s_c^ref²)` etc.).

## 6. `family_priors.tsv`

```
family_id	log_clock_mean	log_clock_sd	evidence_grade	source_key	version
indo_european	0.00	0.40	B	META_IE_001	1.0.0
UNKNOWN	0.00	0.80	D	DEFAULT	1.0.0
```

`UNKNOWN` is the fallback family prior (broad). Used for `β_f`.

## 7. `analysis.tsv` (run configuration)

One key–value pair per row:

```
key	value
model	novel_state
reference_mode	family_specific
clock	hierarchical
contamination	global
chains	4
warmup	5000
samples	10000
seed	123456789
time_origin_year	1950
```

Unknown keys are errors unless namespaced `x_custom_*`. Every output records
the selected mode.

## 8. Tree input (Newick)

Rooted Newick with labelled leaves and labelled internal nodes:

```
((italian,spanish)romance_west,romanian)romance_root;
```

The minimal independent parser (M1) supports:
- leaf and internal **labels**;
- **branch lengths** (after `:`);
- **polytomies** (nonbinary);
- **comments** in `[...]` only if simple to support safely (otherwise
  rejected with an error; never silently misparsed).

In dating mode, supplied branch lengths may be **ignored**, used as
**starting values**, or interpreted per an explicit configuration key
(`branch_lengths = ignore|start|duration`). **Never** guessed. Ages come
from calibrations and the MCMC, not from branch lengths, unless configured.

## 9. Output files (strict TSV)

`posterior_nodes.tsv`, `posterior_rates.tsv`,
`posterior_hyperparameters.tsv`, `posterior_concepts.tsv`, `diagnostics.tsv`,
`pointwise_loglik.tsv`, `run_manifest.tsv`, `samples.tsv` (or `samples.bin`
with a documented binary layout). The run manifest carries: program version,
git commit, compiler, platform, native/wasm mode, input checksums,
prior-package version, seed, configuration, model, reference mode, warnings,
start/completion metadata. See `docs/inference.md` §outputs.

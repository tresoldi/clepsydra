# CLI manual

The `clepsydra` binary (built by `make native` into `build/clepsydra`) runs the
validated inference engine over the strict-TSV inputs of `docs/input-formats.md`.

```
clepsydra <command> [options]
```

Running a command with no options prints its usage. `clepsydra version` prints
the library version; `clepsydra --help` lists the commands.

## `date-tree` — fixed-topology chronology (T1)

Estimates internal node ages on a fixed rooted tree from cognate data, by
Metropolis-within-Gibbs over the strict novel-state likelihood.

```
clepsydra date-tree --tree T.nwk --forms forms.tsv --calibrations cal.tsv --out DIR \
    [--warmup 2000] [--samples 2000] [--chains 4] [--seed 42] \
    [--rate-logmean -8.80] [--rate-logsd 0.80] [--epsilon 0.05] [--relaxed 0|1] \
    [--progress 0|1]
```

| option | meaning | default |
|---|---|---|
| `--tree` | rooted Newick; tip labels are `language_id`s, internal nodes labelled | — |
| `--forms` | forms TSV (`language_id, concept_id, cognate_id, primary_weight`) | — |
| `--calibrations` | calibrations TSV (`docs/input-formats.md` §4) | — |
| `--out` | output directory (created) | — |
| `--warmup` / `--samples` | adapting / retained iterations per chain | 2000 / 2000 |
| `--chains` | independent chains (per-chain seeds off `--seed`) | 4 |
| `--rate-logmean` / `--rate-logsd` | prior on `log r` (per **year**) | −8.80 / 0.80 |
| `--epsilon` | contamination weight for incompatible (homoplastic/loan) concepts | 0.05 |
| `--relaxed` | `1` adds per-branch rate multipliers + `σ_b` | 0 |
| `--block-age-scaling` | decorrelate nested ages; helps on deep trees | 1 |
| `--hpd-level` | credible level for the reported interval | 0.90 |
| `--trace-draws` | evenly spaced draws kept per chain per parameter | 0 |
| `--progress-every` | progress interval in iterations; `0` is about every 1% | 0 |
| `--cal-drop` | withhold this calibration row, by index | keep all |
| `--cal-widen` | multiply every soft calibration window by this | 1.0 |
| `--progress` | `1` reports sampling progress on stderr, one rewritten line per chain | 0 |

The five options above `--progress` were reachable only from the browser until
D4; they and the CLI's flags are now one table (`include/gd_options.h`), so the
two surfaces accept the same set with the same bounds. `clepsydra date-tree`
with no arguments prints the table, with each default read out of the library
rather than restated.

`--progress` uses the same callback the browser page uses to drive its progress
bar and its Cancel button (`docs/web.md`); it does not alter the draws, which
`test_mcmc_progress_does_not_change_draws` checks bit-for-bit.

The rate prior default (`−8.80`, sd `0.80`) is the legacy Swadesh-100 anchor per
year — a **weak fallback**, not a universal constant. Incompatible concepts are
absorbed by the contamination component (never silently dropped). The block
age-scaling move is enabled to decorrelate nested ages on deep trees. Absolute
ages depend on the rate prior and calibrations — the tool says so.

**Outputs** (strict TSV in `DIR/`):
- `posterior_nodes.tsv` — per free node:
  - **calendar age**: mean, median, sd and 90% HPD in years BP;
  - R̂, bulk ESS and tail ESS;
  - **lexical age** (`lexical_*`, `rho · age`): expected replacements per
    concept since the node, which is what the cognates measure;
  - **relative age** (`relative_*`): the age as a fraction of the root's,
    which needs no calendar scale;
  - the node under the **prior alone** (`prior_*`, from a likelihood-free run;
    `--compare-prior 0` skips it), and `prior_overlap`, the overlap of prior
    and posterior on log age (1 means the data said nothing about the node).

  Missing values are written as `.`.
- `diagnostics.tsv` — every sampled parameter: mean, sd, R̂, bulk and tail ESS,
  MCSE of the mean, acceptance.
- `draws.tsv` (with `--write-draws 1`) — every retained draw of every
  sampled parameter, one row per (chain, draw).
- `run_manifest.tsv` — program/commit, model, sizes, seed, priors, the
  branch-rate confounding flag, where the calendar scale comes from
  (`scale_source`: `calibration`, `tip`, both, or `clock_prior`), how far the
  clock moved off its prior (`rate_contraction`; `scale_from_prior` is 1
  below 0.5), the lexical root (median and 90% HPD) and its prior, where the
  calendar root's width comes from (`root_var_share_data` and
  `root_var_share_clock`, the shares of Var(log root) due to the lexical root
  and to the clock; Arawan uncalibrated is 0.01 and 1.0), the power-scaling
  sensitivities, the convergence verdict (`converged`,
  `n_unconverged_params`, `max_rhat`, `min_ess_bulk`, `min_ess_tail`), and the
  scale-dependence note.

**Concepts** (`posterior_concepts.tsv`). One row per concept. `incompatible`
marks a concept whose pattern no single-state T1 history on the given
topology can produce: typically an ancestral synonym pair kept by different
daughters, sometimes a loan. Such concepts are excluded. Their share,
reported in the manifest, measures how much of the data the model cannot
use; on Romance it is 17%. For the rest, `p_contaminated` is the posterior
probability that the concept came from the noise component, ε, which is
sampled with a Beta(1, 19) prior. `relative_rate` is the concept's posterior
rate relative to the mean, under the concept-rate mixture.

**Among-concept rates** (`--concept-rates 1`, the default). Basic meanings
are replaced at very different rates. A single rate makes slow meanings look
young and fast ones look old, and it biases deep nodes. Each concept's rate is
integrated over a log-normal with mean 1 and a sampled spread `sigma_c`.
`--concept-rates 0` gives the homogeneous model. `--concept-rate-gamma 1`
swaps the log-normal for a gamma with mean 1 and coefficient of variation
`sigma_c`. Its slow tail is heavier: gamma mixing gives power-law retention.
Deep dates depend on that tail, so the two families are a model-uncertainty
axis to report, not a choice to make silently.

**Adequacy** (`adequacy.tsv`). Posterior predictive checks compare
statistics of the observed data with data replicated from 200 posterior
draws. The statistics are the classes per concept, the singleton share, the
share of incompatible concepts, and the cognate sharing for every pair of
tips. A `ppp` near 0 or 1 marks a feature the model does not reproduce. A
pair of tips sharing more than predicted suggests borrowing between them, or
a misplaced tip.

**Power-scaling sensitivity.** `sensitivity_root_age` and
`sensitivity_lexical_root` measure, for each component (the clock prior, the
lexical-root prior, the calibrations, the likelihood), how far the posterior
of log root age or log lexical root moves when that component is raised to a
power α, at α = 1. This is the first-order version of Kallioinen et al.
(2024): a covariance with the component's log density, as the larger shift
of mean or sd, in posterior sd. Above 0.2 the component is influential.
Prior-sensitive and not likelihood-sensitive means prior-dominated; sensitive
to both means prior-data conflict. A flat-topped calibration (a uniform
window) is invisible to power-scaling, because its log density is constant
inside the window. Judge those from the `prior_*` columns instead.

**Convergence check.** Every sampled parameter must have rank R̂ < 1.01 (the
larger of the bulk and folded R̂) and bulk and tail ESS ≥ 100 per chain
(Vehtari et al. 2021). A failing run still writes its outputs, prints which
criteria failed, and exits with status **3**. `--allow-unconverged 1` makes it
a warning; use it to inspect a failure, not to publish one.

## `date-pair` — pairwise divergence (P1)

Deterministic pairwise divergence date for two languages (no MCMC).

```
clepsydra date-pair --forms forms.tsv --lang-a A --lang-b B --out DIR \
    [--t-lo 100] [--t-hi 8000] [--rate-logmean -8.9] [--rate-logsd 0.6] \
    [--age-a 0] [--age-b 0]
```

`--t-lo/--t-hi` bound the divergence time in years BP; the rate prior
(`--rate-logmean/-logsd`) is the log replacement rate **per year**, as in
`date-tree` (before 0.3 it was per kiloyear: subtract 6.9078 from an old
value). `--age-a/--age-b` are the sampling dates of the two lists in years BP
(0 = modern); the decay laws refuse them. Concepts present in both languages
count as a cognate match iff they share a `cognate_id`; a missing cognate id
(`.`, `?`, empty, or `status` missing/excluded/absent) is not a match.

**Output** `DIR/posterior_pair.tsv`: shared-concept and cognate-match counts,
retention fraction, divergence mean/median/sd and 90% equal-tail interval (years
BP), the two tip ages and the lineage separation `2T - age_a - age_b` (median
and 90% interval, in years; for an ancestor-descendant pair it is the time
elapsed between the lists), the `absolute_anchored` flag, and the
prior-dominance diagnostic.

## `date-trees` — dating over a sample of trees

```
clepsydra date-trees --trees SAMPLE.trees --forms forms.tsv \
    --calibrations cal.tsv --out DIR [--clades clades.tsv] \
    [--burnin 0.25] [--ntrees 50] [--draws-per-chain 100] [date-tree options]
```

This dates a posterior sample of topologies, which `date-tree` cannot do.
`--trees` takes plain Newick (one tree per `;`) or NEXUS with a
`Translate` table and `[&…]` annotations, as written by BEAST, MrBayes and
RevBayes. After burn-in, `--ntrees` trees evenly spaced through the sample
are dated with the full model, one after another. Sampled trees have no
internal labels, so calibrations must target taxon sets (`mrca:A,B`) or
tips. Tip names must match the forms.

**Output.**
- `posterior_clades.tsv`: for the root and each row of `--clades`
  (`clade_id`, comma-separated `members`), the fraction of trees in which
  the clade is monophyletic, and its MRCA age pooled over trees with equal
  weight. The MRCA of a named set is defined in every tree. The age is also
  reported conditional on monophyly.
- `trees.tsv`: one row per dated tree, with its convergence and root
  median.

The command exits with status 3 if any tree fails the convergence check,
unless `--allow-unconverged 1`. It skips the per-tree prior run
(`--compare-prior`) to halve the cost.

**Caveat.** When the trees were estimated from the same cognates, the result
is a hybrid, not one joint posterior: the data choose the topologies and
then date them again. Report it as such, and compare with a split-by-concept
run.

## `calib-sensitivity` — calibration-sensitivity workflows

How much does the T1 chronology depend on the calibrations you supplied? This
command re-dates the tree under systematically varied calibrations and reports
how the node ages move. Two modes:

```
clepsydra calib-sensitivity --tree T.nwk --forms forms.tsv --calibrations cal.tsv \
    --out DIR --mode loco|widen [--factor 2.0] [--include-extant 0] \
    [--warmup N] [--samples N] [--chains N] [--seed S] \
    [--rate-logmean X] [--rate-logsd X] [--epsilon X] [--relaxed 0|1]
```

- **`--mode loco`** (leave-one-calibration-out, default). Drops each calibration
  in turn, re-dates with the rest, and for the dropped calibration's target node
  reports the withheld value, the inferred age + 90% HPD, whether the value is
  **recovered** (in the HPD), the absolute error, and `max_other_shift_bp` (how
  far the rest of the chronology moves). If dropping a calibration leaves its
  target age **unbounded** (the calibration is the only thing constraining it),
  the row is marked `unbounded`/`unbd` — a load-bearing calibration, not a
  failure. Age-0 (modern) tip calibrations are skipped unless `--include-extant 1`.
- **`--mode widen`**. Loosens every soft calibration window by `--factor` (uniform
  bounds expand about their centre; normal/lognormal SDs scale) and re-dates once,
  reporting each node's baseline vs widened median, the shift, and both interval
  widths — i.e. how much the dates depend on calibration *tightness*.

**Output** `DIR/calibration_sensitivity.tsv` (schema depends on the mode) plus a
one-line summary on stderr. Each mode runs the full sampler once per variation,
so it costs roughly (number of calibrations + 1) dating runs in `loco` mode.
Example on the committed Romance data:

```
clepsydra calib-sensitivity --tree examples/romance_tree/tree.nwk \
    --forms examples/romance_tree/data/forms.tsv \
    --calibrations examples/romance_tree/calibrations.tsv \
    --out /tmp/rom_cs --mode loco
```

## Notes

- The engine never calls `exit()`/`abort()`; malformed inputs print an error and
  return non-zero.
- Node ages are years before `time_origin_year` (default 1950, "BP"); the C
  engine never parses calendar strings.
- For reproducible worked examples on real cognate data (eleven families), see
  `examples/` and `docs/scientific-methods.md`.
- To date **your own** data: convert it to the forms TSV first — from a Lexibank
  CLDF release with `tools/datasets/lexibank.py`, or from a plain cognate matrix
  (wide languages×meanings, or a long table) with `tools/datasets/matrix.py` —
  then point `date-tree` / `date-pair` at the result.

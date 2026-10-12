# Scientific methods and honest interpretation

How to read a clepsydra result without over-claiming. The worked examples in
`examples/` and the Arawan analysis (`docs/research/arawan-result.md`)
follow it.

## What the number is (and is not)

A clepsydra date is a **posterior distribution** under a stated model,
tree, rate prior and set of calibrations. It is reported with a credible
interval and convergence diagnostics, never as a single "true" date.

"Proto-Romance diverged in 1662 BP" is wrong. This is right: "the posterior
median for Proto-Romance is about 1,660 BP (90% HPD 1,383–1,959; R̂ 1.00,
ESS 2,928), with Latin fixed at 2,050 BP, the default T1 model with a
relaxed clock, and the default rate prior".

## Three kinds of age

Cognate data measure change, not time. Every node therefore gets three ages
in `posterior_nodes.tsv`, and they are not equally well known.

| age | columns | what it is | what sets it |
|---|---|---|---|
| **lexical** | `lexical_*` | ρ · a: expected replacements per concept, at the mean concept rate, from the node to the present | the cognate data |
| **relative** | `relative_*` | a / a_root: the node's depth as a fraction of the root | the cognate data and the tree prior |
| **calendar** | `median_bp`, `hpd90_*` | a, in years before present | the lexical age divided by the clock rate ρ |

Lexical and relative ages are what the data identify. A calendar age is
only as good as whatever fixes ρ: dated tips, node calibrations, or the
clock prior. For Arawan the lexical root is 0.64 (0.51–0.79) under every
variant tried. The calendar root ranges from 1,950 to 4,250 BP, depending
on which source fixes the time scale.

## Two modes, and the scale-source flag

The engine runs in one of two modes and records which in `run_manifest.tsv`.

- **Lexical-time mode** (`scale_source clock_prior`). No calibration, no
  dated tip. The ages are given a proper prior on the lexical root and on the
  ratios below it (the default `uniform_ratios`). The clock rate is
  independent of the data, so the calendar posterior is the lexical posterior
  rescaled by the clock prior: the data cannot narrow it beyond the clock's
  own spread. The cross-family pilot puts that spread at σ ≈ 0.7–0.8 in
  log rate (`docs/research/clock-pilot.md`), ×/÷ 5–6 at 95%. For an
  uncalibrated family, report lexical and relative ages, plus calendar
  scenarios labelled as such.
- **Calendar time-tree mode** (`calibration`, `tip`, or `calibration+tip`).
  At least one calendar fact enters the fit. A dated tip or a calibrated
  node informs ρ jointly with the cognates. A user's "root window" is a
  calibration like any other, with a source key and an optional outlier
  probability, never "the root prior".

Being in calendar mode does not mean the calendar fact did much. Two fields
say how much:

- `rate_contraction`: how far the posterior of log ρ narrowed from its
  prior (1 − posterior variance / prior variance). `scale_from_prior 1` when it is below 0.5,
  meaning the calendar ages still rest mostly on the clock prior. Arawan's
  two historical lists give 0.32: they move the scale, but do not fix it.
- **The variance split** `root_var_share_data` / `root_var_share_clock`.
  Since log a_root = log L − log ρ, the variance of the log root age splits
  exactly into the part that moves with the lexical root
  (Cov(log a, log L) / Var(log a)) and the part that moves with the clock
  (Cov(log a, −log ρ) / Var(log a)). The two sum to 1. Arawan with tip dates
  gives 0.03 / 0.97: almost all of the calendar width is the clock's.
  When a calibration pins the root, L and ρ co-vary, and one share can leave
  [0, 1] (five-IE: −0.45 / 1.45). That means the calibration, through ρ,
  is what holds the root, and extra lexical spread is cancelled by the rate.

## The choices that set the answer

Every date rests on choices `run_manifest.tsv` records:

1. **The model.** T1 (novel-state, with concept-rate mixing and ε
   contamination) is the engine. T2 is experimental. P0 and P1 are pairwise.
   See `docs/model-cards.md`.
2. **The calibrations.** The examples' calibrations are illustrative and
   flagged as such. A calibration with no stated source is a guess, and the
   output will carry it.
3. **The rate prior.** In lexical-time mode it *is* the calendar answer. In
   calendar mode it still dominates whenever `scale_from_prior` is 1.
4. **The topology.** `date-tree` uses one fixed tree. `date-trees` dates a
   sample of trees (BEAST, MrBayes, yvyra) and reports clade probabilities,
   with ages conditional on each clade. A clade present in 35% of trees has
   an age only in those trees. This is a modular sensitivity analysis, not a
   joint posterior over topology and ages: each tree is dated as given, and
   when the sample was estimated from the same cognates the data are used
   twice. Only trees that converge are pooled; the output says how many.
5. **The observation model.** Concepts that no single-state history can
   produce are absorbed by ε and listed in `posterior_concepts.tsv`. They are
   never silently removed. Lists known to be unusual (historical,
   non-specialist, secondary) can carry a `tip_error` term.

## Reading the diagnostics first

Before interpreting any age:

- **Convergence.** Every sampled parameter needs rank R̂ < 1.01 and bulk and
  tail ESS ≥ 400. A failing run exits with status 3 and says which
  parameters failed.
- **Prior against posterior.** The `prior_*` columns give each node's age
  under the prior alone, from the same run. `prior_overlap` near 1 means
  the data did little for that node. `tools/plotting/overlay.py` draws the
  two intervals side by side; every tree example has one. In five-IE the
  root's overlap is 0.85: its width is the priors', not the cognates'.
- **Power-scaling** (`sensitivity_root_age`, `sensitivity_lexical_root`).
  This measures how far the root moves when one component (clock prior,
  lexical-root prior, calibrations, likelihood) is up- or down-weighted.
  Above 0.2 the component is influential. If the lexical root is not mostly
  likelihood-sensitive, the data are not driving it.
- **`branch_rate_confounded`** (relaxed clock). Branch rates and node ages
  trade off, so a flagged run needs calibrations before its branch rates
  mean anything.
- **`absolute_anchored`** (pairwise). A value of 0 means the scale is the
  rate prior's.

## Validity, adequacy, model choice

Three distinct questions, kept separate:

- **Validity:** does the code target the intended posterior? It is
  established by analytic checks, independent brute force, and
  simulation-based calibration (SBC) with the real likelihood
  (`docs/validation.md`).
- **Adequacy:** is the model right for these data? It is checked with the
  posterior-predictive panel in `adequacy.tsv`: classes per concept,
  singleton share, the share of incompatible concepts, and sharing for
  every pair of tips. A ppp near 0 or 1 is a feature the model does not
  reproduce. Arawan fails one: 7.5% incompatible concepts against at most
  1.2% replicated.
- **Model choice:** a more flexible model is not automatically better.
  Compare models on held-out prediction over concepts. `date-tree
  --write-pointwise 1` writes the per-concept log-likelihoods, and
  `clepsydra average` computes PSIS-LOO and stacking weights. It also
  reports the between-model share of each node's uncertainty. On Arawan
  the four strict/relaxed × log-normal/gamma models differ by less than
  their SE, and model choice accounts for 0.2% of the root's calendar
  uncertainty. On Romance the relaxed clock wins by 10.9 ± 4.6.

## Against classical glottochronology

`tools/glottochronology/dating.py` dates a tree by Swadesh (1955), Lees
(1953) and Starostin (2000), the methods a reader may compare clepsydra
with. `docs/research/arawan-glottochronology.md` runs them on Arawan and on
the Romance control. Three things to expect when making that comparison:

- **They disagree with each other.** On the same Arawan cognates the root
  ranges from 1,750 to 6,600 BP by word list and formula. Their bootstrap
  intervals look narrow because they treat the rate as known.
- **They read list error as time.** The Paumari 1833 list comes out
  1,000–2,500 years old instead of 167. Clepsydra's list-error term covers
  the truth.
- **They do not give a chronology of a tree.** Pooled over pairs, they can
  date a child older than its parent.

Where they are right (shallow nodes; the Romance root under Lees),
clepsydra agrees. Their point estimates sit inside clepsydra's intervals.

## Known weak points, measured

The simulation study (`docs/research/simworld.md`) breaks one assumption at a
time and scores the 90% interval of each lexical age:

- **Absorbed:** borrowing in 10% of concepts, fast concepts missing more
  often, gamma instead of log-normal rate variation. Coverage stays at
  0.81–0.90.
- **Not absorbed:** coding error. With 8% wrong forms per list, ages are
  ~55% too old (coverage 0.32), worst at shallow nodes. Careful coding
  matters more than any model choice.
- **Not absorbed:** a rate that falls with depth. Deep ages come out ~30%
  too young, and the root is covered once in 20. Only calibrations at more
  than one depth can reveal this.
- **Measured on public data (GLED, 40 concepts, no calibration):** calendar
  ages are prior-scaled and biased young at depth, and ASJP's published
  formula is closer to the published family dates (roots MAE 1,535 against
  3,191 years). Shallow calibrations do not extrapolate to deep roots
  (`docs/research/gled-benchmark.md`). Read uncalibrated calendar ages as
  scenarios (`calendar_basis`).
- **Clean data:** coverage is 0.81 against 0.90. That is within about
  two standard errors, given that nodes from one tree are correlated.

## How to report a result

State the following:
- the model;
- the mode and `scale_source`;
- the calibrations, with their sources;
- the rate prior;
- `rate_contraction` and the variance split;
- the adequacy failures, if any;
- the model comparison or sensitivity runs;
- the topology source.

Lead with lexical and relative ages when the calendar scale is
prior-driven. Give intervals, not point dates. When the data cannot answer
the question, say so and show the width. The deep IE and Sino-Tibetan roots
are examples.

**A run from the browser page is not a result to report.** The page
(`docs/web.md`) computes exactly what the CLI computes, and is checked to the
last digit. But its default chain lengths are chosen so a run finishes while
someone waits, and it keeps no manifest. Use it to explore and to teach. Do
the run you intend to cite with `clepsydra date-tree`, which writes
`run_manifest.tsv` beside the posterior.

## The worked examples

- `examples/pair_simple` — P1 pairwise (Spanish–Portuguese), the smallest case.
- `examples/romance_tree` — T1 chronology of Romance, anchored by attested
  Latin. Proto-Romance at about 1,660 BP. `calib_sensitivity/` shows which
  calibrations the chronology leans on.
- `examples/five_ie_reproduction` — IE tip-dating from five ancient
  languages. The root is about 6,580 BP (4,190–9,430): wide,
  rate-prior-sensitive, and reported with its full uncertainty.
- `examples/sinotibetan_tree` — Sino-Tibetan (Sagart et al. 2019).
  Subgroup ages are well estimated. The root is deep and wide (about 8,250
  BP, 6,810–9,170), a weakly identified deep root reported as such; under a
  strict clock it presses past its soft window.
- `examples/polynesian_tree` — a shallow family (Walworth 2018) with a
  node-calibrated root (about 2,840 BP) and well-converged internal ages.
  The clock choice moves the Central-Eastern core by a quarter (relaxed
  ~1,270 BP, strict ~1,740 BP).
- `examples/families_pairwise` — a nine-family pairwise gallery, reported as
  unanchored lexical time.

All run from a clean checkout (`make`, then each `run.sh`) on real,
CC-BY-4.0 cognate data from Lexibank / LexiCore (`examples/DATA_SOURCE.md`).

# Scientific methods and honest interpretation

> Milestone 12 deliverable (product brief §17). How to read a Clepsydra result
> without over-claiming. The worked examples in `examples/` follow this to the
> letter.

## What the number is (and is not)

A Clepsydra date is a **posterior distribution**, reported with a credible
interval and convergence diagnostics — never a single "true" date. A posterior
median is a summary of that distribution, not a fact. "Proto-Romance diverged in
1590 BP" is wrong; "the posterior median for the proto-Romance node is ~1590 BP,
90% HPD 1400–1810 BP, R̂ 1.00, ESS 740, under this tree, rate prior, and
calibrations" is right.

## The five things that set the scale

Every absolute date rests on choices the output records in `run_manifest.tsv`:

1. **The model.** T1 (strict novel-state), T2 (binary), or pairwise P0/P1 encode
   different assumptions (see `docs/model-cards.md`). A date is a date *under a
   model*.
2. **The calibrations.** Absolute time is not identifiable from cognacy alone —
   only lexical time (`r·d`) is (identifiability §7). Tip ages and node priors
   supply the scale; the examples' calibrations are **illustrative**, flagged as
   such, not authoritative.
3. **The rate prior.** With few or no calibrations, the clock-rate prior largely
   determines the absolute ages. The five-language IE example shows this
   starkly: the PIE-root posterior shifts materially with `--rate-logmean/-logsd`,
   and its 90% HPD spans the Steppe and Anatolian hypotheses. That width is the
   honest answer, not a defect to be tuned away.
4. **The topology.** v1 uses a fixed rooted tree; it does not infer or search
   topologies. A different tree is a different analysis.
5. **Contamination / loans.** Incompatible concepts (homoplasy, borrowing) are
   absorbed by the contamination component and reported, never silently removed.

## Reading the diagnostics first

Before interpreting any age:

- **R̂** near 1 (≤ ~1.01) and **ESS** large (hundreds+) — otherwise the chains
  have not converged and the summaries are unreliable. A short run says so.
- **`branch_rate_confounded`** (relaxed clock) — branch rates and node ages
  trade off; a flagged run needs calibrations before branch-rate estimates mean
  anything (identifiability §7.2, `test_mcmc_relaxed.c`).
- **`absolute_anchored`** (pairwise) — `0` means the scale is prior-driven.
- **prior-dominance** — a posterior close to the prior carries little data
  signal (common with few concepts or languages).

## Validity is not adequacy is not flexibility

Three distinct questions, kept separate (King 2026; `docs/model-cards.md`,
`docs/validation.md`):

- **Validity** — does the code target the intended posterior? Established by
  analytic checks, independent brute force, and SBC.
- **Adequacy** — is the model right for *these* data? Probed by
  posterior-predictive checks; a validly-implemented homogeneous-rate model can
  still be inadequate for rate-heterogeneous data, and the PP check catches it.
- **Flexibility** — a more flexible model is **not** automatically better. On
  strictly-clocked data the relaxed clock does not beat the strict clock on
  held-out prediction. Model choice uses held-out prediction in a shared
  observable space, never in-sample fit or cross-encoding marginal likelihoods.

## How to report a result

State the model, the calibrations that set the scale, the prior package /
rate prior, whether family data updated the rates, the split between rate and
calibration uncertainty, the strict-inheritance fit (contamination /
compatibility), cross-model sensitivity, and the influential concepts and
calibrations. Give intervals, not point dates. When the data cannot answer the
question (deep-time IE root), say so and show the width — as the five-language
example does.

**A run from the browser page is not a result to report.** The page
(`docs/web.md`) computes exactly what the CLI computes — that is checked to the
last digit — but its default chain lengths are chosen so a run finishes while
someone waits, and it keeps no manifest. Use it to explore, to teach, and to see
whether a dataset behaves; do the run you intend to cite with `clepsydra
date-tree`, which writes `run_manifest.tsv` beside the posterior so the result
can be reproduced.

## The worked examples

- `examples/pair_simple` — P1 pairwise (Spanish–Portuguese), the smallest case.
- `examples/romance_tree` — T1 chronology of Romance, anchored by attested Latin;
  converges cleanly and dates proto-Romance to ~1600 BP with honest intervals.
- `examples/five_ie_reproduction` — reproduction-style IE tip-dating; a hard,
  rate-prior-sensitive deep-time case reported with its full uncertainty.
- `examples/sinotibetan_tree` — a **different family** (Sino-Tibetan, Sagart et
  al. 2019): well-estimated subgroup ages (Sinitic, Tibetic, Burmish) alongside a
  root that is honestly reported as weakly identified under a strict clock — the
  same validity/identifiability lessons on unrelated data.
- `examples/polynesian_tree` — a **shallow, cleanly-dated** family (Polynesian,
  Walworth 2018): a node-calibrated Proto-Polynesian root and crisply converged
  internal ages (~2900 BP root, Eastern Polynesian ~2250 BP), the clean-case
  counterpoint to the deep IE/ST roots.
- `examples/families_pairwise` — a **nine-family** pairwise gallery (Dravidian,
  Turkic, Uralic, Semitic, Polynesian, Tupían, Arawakan, Alor-Pantar, Koreanic):
  the same engine and pipeline run unchanged across unrelated families, with
  divergences honestly reported as unanchored lexical time.

All run from a clean checkout (`make native`, then each `run.sh`) on real,
CC-BY-4.0 cognate data from Lexibank / LexiCore (`examples/DATA_SOURCE.md`).

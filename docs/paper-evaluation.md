# Paper evaluation

> Milestone 0 deliverable. Evaluates the required references for Clepsydra.
> The primary one — Rönchen, Wiklund & Hammarström 2024 — is split into
> **adopted / adapted / rejected / postponed** as the brief (§3) demands. The
> others are summarized for how Clepsydra uses them. No code is copied from
> any implementation; the model is independently re-derived in
> `docs/model-spec.md`.

## Primary paper

**Philipp Rönchen, Tilo Wiklund & Harald Hammarström.** *Likelihood
calculation in a multistate model of vocabulary evolution for linguistic
dating.* Language Dynamics and Change 14(1), 2024, 1–41.
doi:10.1163/22105832-bja10032.

The "Multistate model" (Warnow et al. 2004 / Mossel & Steel 2004): per
meaning, vocabulary evolves along a tree by replacement of the current
cognate class with a **new, never-before-seen** class. No back mutation, no
parallel innovation. One primary cognate per language–meaning. The paper's
main contribution is an exact, ~O(n log n) likelihood algorithm: (1) extend
the partial data `f_m` to its unique extension `f̃_m` (path property /
extendability); (2) a bottom-up double recursion on `p_v` (subtree
likelihood) and `q_v` (likelihood × preserves-a-data-value-to-v), with four
cases (i, iia, iib, iic) for binary trees (Theorem 10), generalized to
nonbinary trees in appendix Theorems A.8/A.9.

### Adopted (independently implemented)

- **Meaning-wise multistate cognate representation** (Table 3 coding): one
  phylogenetic character per meaning with multiple states.
- **Lexical replacement as transition to a novel cognate state** (no back
  mutation, no parallel innovation) — the generative process of
  `docs/model-spec.md` §3.3.
- **Support for cognate states not observed at the leaves** (unseen ancestral
  states), handled natively by the `p_v`/`q_v` recursion — no ascertainment
  correction needed for T1, unlike binary/Stochastic-Dollo.
- **The path-property interpretation of strict vertical inheritance** (Def 2,
  Prop 3): an outcome has positive probability iff it has the path property.
- **Extendability / compatibility checking** (Defs 4–6, Prop 7): incompatible
  data ⇒ `P(f_m)=0`. Used as the basis for T1's `incompatible` handling.
- **The `p_v`/`q_v` double recursion (Theorem 10)** as the T1 likelihood,
  re-implemented in C with log-space numerics (`docs/model-spec.md` §10).
- **Efficient tree likelihood** — O(n log n) extension + O(n) recursion,
  linear in concepts; the complexity target for compatible data.
- **Missing observations** marginalized via the extension/recursion (a node
  not in `dom(f_m)` is unconstrained).
- **Lexical observations on internal nodes** (ancient languages as observed
  internal/leaf nodes) — Case ii of the recursion.
- **Nonbinary trees** — supported in principle (the paper covers them in the
  appendix); see *Postponed* for the current status.
- **Meaning-specific change rates** `c_m` — directly maps to our concept
  rates `r_c` and the hierarchy of `docs/model-spec.md` §4.
- **Simulation from the same model** for validation (their generative process
  = our `simulate`).

### Adapted (used, but changed for this project)

- **Rates.** The paper uses a single per-meaning rate `c` (or a fixed lexical
  clock). Clepsydra replaces this with the **hierarchical log-rate**
  `log r_{f,c} = μ_c + β_f + u_{f,c}` and a relaxed branch clock
  (`docs/model-spec.md` §4), so family data *updates* reference priors and
  branch-rate variation is explicit. The paper itself notes relaxed clocks
  are straightforward to add Bayesianly (citing Drummond et al. 2006).
- **Inference.** The paper does ML grid search + bootstrap on a fixed small
  tree. Clepsydra is **Bayesian** (posteriors mandatory; MCMC for trees,
  deterministic integration for pairwise), with calibrations, contamination,
  and full diagnostics.
- **Compatibility failures.** The paper excludes incompatible meanings
  (likelihood 0). Clepsydra wraps T1 in a **robust contamination mixture**
  (`L_c = (1−ε_c)L_tree + ε_c L_noise`) so one incompatible concept does not
  zero the dataset, and reports `P(c contaminated)` — directly addressing the
  paper's central stated weakness (§7: "central weakness … lack of robustness
  against non-genealogical developments").
- **Synonyms.** The paper assumes exactly one primary cognate per
  language–meaning. Clepsydra keeps this for T1 but adds (a) uncertain primary
  state (latent in MCMC / marginalized in pairwise) and (b) the **T2 binary
  model** as a polymorphism-tolerant alternative, preserving the paper's
  input-data semantics where possible.
- **Calibration.** The paper fixes leaf ages and grid-searches internal ages.
  Clepsydra supports the full calibration suite (fixed/ancient tips,
  uniform/normal/lognormal/offset-lognormal nodes, hard/soft bounds) per
  `docs/calibrations.md`.

### Rejected (not carried into Clepsydra)

- **Do not make their model the sole basis.** The brief is explicit: T1 is
  "only one model module." P0, P1, T2 (and reserved T3) coexist.
- **Do not reproduce their software architecture / repository.** Their
  Python implementation (gitlab `py-cognateevolutionmodels`) is not copied or
  mirrored; the C engine is independent.
- **Do not copy code** from their implementation.
- **Do not inherit the limitations silently** (brief §3.2). The strict
  novel-state assumptions (one primary state; every replacement creates a new
  class; no parallel innovation / back mutation / undetected borrowing /
  synonym coexistence; independent meanings; tree-like inheritance) are
  **explicit model assumptions of T1 only**, documented as such in
  `docs/model-spec.md` §3.3 and §12, not project-wide truths.
- **Do not validate only against their implementation** (brief §3.4).
  Clepsydra validates T1 by **independent brute-force enumeration**
  (`docs/model-spec.md` §11); agreement between two implementations of the
  same misunderstood recurrence is not sufficient.

### Postponed (reserved for later, not blocking v1)

- **Nonbinary exact recursion (Theorems A.8/A.9).** The appendix supplement
  (figshare, doi:10.6084/m9.figshare.26089507) is **not in `literature/`**.
  Until supplied, T1 supports binary trees exactly and resolves polytomies to
  zero-length ladders as a **documented approximation**, with brute-force
  checks on small polytomies. Logged in `docs/model-spec.md` §3.3/§12.
  > **Action:** obtain the appendix to implement the exact nonbinary
  > recursion before any published polytomy result.
- **`k`-forms-per-meaning extension** (paper §7): store probabilities of
  preserving 0..k cognate classes per meaning, to handle synonymy/semantic
  shift within the multistate framework. Promising but not in v1; T2 covers
  polymorphism meanwhile.
- **Explicit borrowing/contact edges** (Nakhleh et al. 2005; Kelly & Nicholls
  2017, cited by the paper). NP-hard in general; not in v1. Clepsydra uses
  contamination + loan modes (`docs/model-spec.md` §5.3, §6) instead.
- **Homoplastic-states extension** (Warnow et al. 2006). Deferred to the
  reserved T3 finite-state CTMC (`docs/model-spec.md` §3.5).
- **ML topology search / EM branch-length optimization.** Out of scope
  (fixed topology in v1; topology-ready architecture only).

---

## Other required references

### Benedict King 2026 — *Testing the validity and adequacy of linguistic
phylogenetic analyses.* PLOS Comput Biol 22(5): e1014312.
doi:10.1371/journal.pcbi.1014312.

**Used for:** the validation framework (M11) and one T2 design decision.

- **Simulation-Based Calibration (SBC) + Rank Uniformity Validation (RUV)**
  are the central validity tests: draw from the prior, simulate, refit, check
  rank uniformity (ECDF difference plots) and HPD coverage. Clepsydra's SBC
  suite (`docs/validation.md`) follows this directly.
- **Validation ≠ adequacy.** A model can be validly implemented yet
  *inadequate*. Clepsydra separates the two: SBC tests implementation
  validity; posterior-predictive checks test adequacy. A model card per
  likelihood records both (M11).
- **Ascertainment × partition-rate reweighting miscalibration.** Weighting
  partition rates by the number of *cognate classes* (the BEAST2 default)
  yields poorly-calibrated partition rates, birth rate, and tree height.
  Weighting by the number of *meanings* is well-calibrated and fits better.
  ⇒ Clepsydra's T2 defaults to **meaning-weighted** partition rates; the
  cognate-count weighting is opt-in only for reproducing old analyses
  (`docs/model-spec.md` §3.4).
- **Posterior-predictive inadequacy of the covarion model**: lower proportion
  of present cognates, higher mean pairwise distance, fewer singletons,
  excess variance in #present cognates — attributed to semantic shift and
  non-independent evolution of cognate sets within a meaning. ⇒ Clepsydra's
  posterior-predictive metrics (`docs/validation.md` §6) include these, and
  the multistate (T1) approach is adopted partly because it avoids the
  within-meaning independence assumption King flags.
- King explicitly flags the **multistate model with infinite state space**
  (i.e. T1 / Rönchen et al.) as "most promising" but not yet implemented at
  scale — motivating this project.

### Pagel, Atkinson & Meade 2007 — *Frequency of word-use predicts rates of
lexical evolution throughout Indo-European history.* Nature 449, 717–720.
doi:10.1038/nature06176.

**Used for:** relative per-concept rate evidence in the reference-prior
builder (M3).

- ~200 IE meanings across 87 languages; rates estimated as expected cognate
  replacements per 10,000 years, **scaled to an 8,700-year IE age** (so they
  are *relative*, not universal absolute). ~100-fold variation across
  meanings (mean 3.05 ± 1.82, median 2.74, min 0.09, max 9.27 per 10⁴ yr).
- Word-use frequency predicts ~50% of the variance; part of speech adds
  structure (numbers/pronouns slowest; conjunctions/prepositions fastest).
- **Clepsydra use:** map the ~200 meaning-specific rates to canonical IDs;
  center them in log space
  (`m̃_c = log r_c − Σ w_c log r_c / Σ w_c`) and store the centered
  multiplier **separately** from the overall clock scale. Treated as
  relative-rate evidence, **not** universal absolute rates.
- > **Action / open item:** the per-meaning rate **table** is in the paper's
  > Supplementary Table 2, which is **not** in the PDF text in `literature/`.
  > Need the supplement (or the Dyen–Kruskal–Black source database + a
  > re-fit) to build the Pagel-based prior. Until then the Pagel prior is a
  > stub. Logged in `docs/priors.md` (to be written at M3).

### Vejdemo & Hörberg 2016 — *Semantic Factors Predict the Rate of Lexical
Replacement of Content Words.* PLOS ONE 11(1): e0147924.
doi:10.1371/journal.pone.0147924.

**Used for:** predictive priors for concepts lacking direct historical rate
estimates (M3).

- 167 **content** (open-class) Swadesh items (treat content and function
  words separately — function-word rates are dominated by word class).
- Predictors: frequency, word class, age of acquisition, mutual information,
  arousal, imageability, senses (polysemy), synonyms. Significant in the
  regression: **frequency, synonyms, senses, imageability** (synonyms
  especially).
- **Clepsydra use:** fit predictive priors (with cross-validation — *not*
  in-sample fit) for concepts without direct historical rates; map concepts
  to canonical IDs. Content vs function concepts handled separately.
- > **Action / open item:** the predictor table is in the paper's Supporting
  > Information file, **not** in the PDF text. Need that file (or the
  > underlying corpus/lexicon sources) to build the Vejdemo–Hörberg prior.
  > Stub until supplied.

### Rama & Wichmann 2020 — *A test of Generalized Bayesian dating.* PLOS ONE
15: e0236522. doi:10.1371/journal.pone.0236522.

**Used for:** the 30 calibration groups as a benchmark / leave-group-out
validation resource (M3, M11). **Not** as Clepsydra's clock.

- 30 phylogenetically independent calibration groups (ASJP db v18, 40-item
  Swadesh lists), each with a known age. A Gamma regression of median tree
  length → age gives a *generalized* (no family-specific priors) dating
  formula, validated leave-one-out.
- The 30 groups (Table 1, captured): Benue-Congo 6500, Brythonic 1450,
  Chinese 2000, Cholan 1600, Dardic 3550, East Polynesian 1050, East Slavic
  760, English-Frisian 1550, Ethiopian Semitic 2450, Hmong-Mien 2500, Inuit
  800, Iranian 3900, Ma'anyan-Malagasy 1350, Malayo-Chamic 2400,
  Maltese-Maghreb Arabic 910, Mississippi Valley Siouan 2475, Mongolic 750,
  Ongamo-Maa 1150, Oromo 460, Pama-Nyungan 4500, Romance 1729, Romani 650,
  Saami 1750, Scandinavian 1100, Southern Nilotic 2500, Southern Songhai 550,
  Temotu 3200, Tupi-Guarani 1750, Turkic 2500, Wakashan 2500 (years BP).
- **Clepsydra use:** as benchmark/validation calibrations and for
  leave-group-out experiments — **never** copying their fitted Gamma formula
  as the project's clock. Stored with provenance (ASJP v18, retrieval date,
  checksum, concept mapping, families, exclusion rules).

### Maurits, de Heer, Honkola, Dunn & Vesakoski 2020 — *Best practices in
justifying calibrations for dating language families.* J. Language Evolution
5(1):17–38. doi:10.1093/jole/lzz009. (The `lzz009.pdf` in `literature/`.)

**Used for:** calibration methodology (M5, `docs/calibrations.md`).

- Calibrations are priors on **tree nodes** (MRCA points), not on
  protolanguages-as-branches; the protolanguage-as-gradual-entity vs
  node-as-instantiation tension is documented and reflected in our
  `calibrations.tsv` `target_id` semantics.
- The **clock/age invariance** (double branch lengths, halve the rate ⇒
  identical likelihood) is the formal statement of our identifiability §7.1;
  it is why an anchor is required for absolute time.
- Hard vs soft bounds, calibration justification from loanword analysis and
  historical scholarship, and the Uralic case study (revised calibrations
  shifted the root age ~850 yr and widened the HPD ~3000 yr) motivate our
  required sensitivity procedures (prior-only, leave-one-out,
  one-at-a-time, hard-vs-soft, conflict reporting — `docs/calibrations.md`).
- "Do not silently treat a historical date as exact merely because the input
  gave one number" — enforced: every calibration carries a `type` and
  `is_hard` flag and a `source_key`.

### ASJP (database)

**Used for:** worldwide relative concept stability and broad cross-family
behaviour (M3). A pinned version is stored with version, retrieval date,
checksum, concept mapping, families included, exclusion rules. Precise
absolute per-year rates are **not** derived directly from ASJP lexical
distance; the database is family-balanced so large families do not dominate.

> **Action / open item:** the ASJP database itself is **not** in
> `literature/`. Need a pinned release (the Rama–Wichmann paper used v18) to
> build the ASJP-derived prior. Stub until supplied.

### Stochastic Dollo / missing data / lateral transfer literature
(Nicholls & Gray 2006/2008; Ryder & Nicholls 2011; Kelly & Nicholls 2017)

**Used for:** context, not implementation. The paper compares the Multistate
model to Stochastic Dollo (Table 5): both disallow homoplasy; Multistate is
meaning-based with one rate per meaning, Stochastic Dollo is trait-based with
birth+death. Clepsydra does **not** implement Stochastic Dollo in v1 but
preserves the extension point; the reserved T2 stochastic-Dollo-like
birth-and-loss module is noted in `docs/model-spec.md` §3.4.

---

## Summary of open literature items

| Item | Needed for | Status |
|---|---|---|
| Rönchen 2024 appendix supplement (figshare 26089507) | exact nonbinary T1 recursion (Thm A.8/A.9) | **missing** — polytomies resolved to zero-length ladders (M6). That resolution reproduces our §3.3 generative process exactly (brute-force verified on unresolved polytomies), so what is open is whether A.8/A.9 read a polytomy the same way, plus the direct `O(k)` recursion. No published polytomy result before the appendix is checked. |
| Pagel 2007 Supplementary Table 2 (per-meaning rates) | M3 Pagel-based concept prior | **missing** — stub |
| Vejdemo–Hörberg 2016 Supporting Information (predictor table) | M3 predictive priors | **missing** — stub |
| ASJP database (pinned release) | M3 ASJP-derived stability prior | **missing** — stub |
| IE-CoR / Romance / 5-language datasets | M12 real-data examples | deferred (synthetic fixtures until then, per user) |
| Held-out-family calibration benchmark | real-world leave-family-out validity (infer a held-out family's dates end-to-end, compare to its historical calibrations) | **partial (N=9, three families)** — `tools/validation/calibration_benchmark.py` runs manifest-driven leave-one-calibration-out across Indo-European, Sino-Tibetan, and Romance through the real engine (`examples/loco_benchmark/`): strict-clock 90% coverage 1/9 vs relaxed 7/9 over nine attested-tip held-out points, scored also by interval score against a prior-only baseline (relaxed skill −0.13 overall, +0.13 without the unrecoverable Romance-Latin), a real model-comparison finding the synthetic study cannot produce. Broadening further is data curation (more trees with authoritative attested/published calibrations), a few manifest lines, not a code change. |

None of these block Milestone 0 or M1/M2 (pairwise + foundation, synthetic
fixtures). They block the *real-data* portions of M3 (reference-prior
builder) and M12 (reproducible examples), and the exact-nonbinary portion of
M6. They are tracked here and in `PLAN.md` §12.

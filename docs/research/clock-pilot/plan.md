# D2 σ_f pilot: pre-registration

Written 2026-10-10, before any match count, distance or rate was computed.
At the time of writing the only things looked at were the calibration rows,
the Glottolog structure around each calibrated node, and per-dataset
metadata (number of forms, cognate coverage, language tables). Changes made
after results were seen are listed at the end under "Deviations", with the
reason. The code that implements this plan is `tools/clock_pilot/`.

## Question

For a family with no calendar calibration, the calendar interval is set by
how much the log clock rate varies between families (σ_f), plus how much it
varies between nodes within a family (τ). This pilot estimates both from
calibrated splits in expert-coded Lexibank data, and checks out of sample
(leave one family out) whether a transported clock predicts calibration ages
better than trivial baselines.

## Data, pinned

- Calibrations: `reference/calibrations/calibrations.tsv` at the repository
  HEAD of 2026-10-10. Rows with grade A or B and `target_kind` in
  {node, root, stem, settlement_bound}: 78 rows.
- Per-row decisions: `tools/clock_pilot/curation.tsv` (one line per row,
  with the reason). That file is part of this plan.
- Cognate data: the 59 Lexibank datasets flagged `ExpertCognates` in
  `~/work/chl/arcaverborum/datasets.csv`, read from
  `~/work/chl/arcaverborum/raw/lexibank/<name>/cldf/`. The git commit of each
  checkout is written to `results.tsv`. No other cognate data are used
  (in particular not `abvd`, which is local but not flagged).
- Classification: Glottolog CLDF v5.0 (`glottolog-cldf-4dbf078`), local zip,
  SHA-256 `e954bd77…a433`.

## Inclusion and exclusion

1. **Excluded outright** (curation `use = exclude`):
   - creoles (CR-02, CR-03, CR-05, CR-11);
   - circular rows, whose date is itself a lexical or Bayesian dating result
     (NA-020 Uto-Aztecan root, AU-02 Pama-Nyungan root). No ASJP-derived date
     is used: Holman (HOL-*) rows are used only for their non-lexical
     "calibration point", never for the ASJP estimate;
   - nodes that cannot be measured with contemporary languages: Latin +
     Romance (Latin is an ancestor), Tocharian A/B (both ancient), Semitic
     root (East Semitic is extinct), the Hungarian stem (Ugric is not a
     Glottolog node; the stem would be the Uralic root, and the min is
     uninformative).
2. **Age treatment**, from `bound_type` and the notes, decided per row:
   - `interval`: two-sided; log T uniform on [lo, hi]. Max-leaning "points"
     are used as two-sided intervals only where the row's
     `suggested_distribution` gives an explicit two-sided range.
   - `point`: log T normal, centred on the point, SD 0.15 (Holman's 15%).
   - `max` / `min`: one-sided. Excluded from the main σ_f fit; used in a
     separate censored fit and scored separately.
3. **Tip age.** Modern word lists were collected around 1950–2010. All BP
   ages get +50 years, so T is time from the split to about 2000 CE.
4. **Ancient doculects** are dropped before matching: a dataset `historical`
   flag, an `Is_Proto` flag, or a name or Glottolog name matching Old,
   Middle, Ancient, Classical, Archaic, Proto, Early, Late, Epigraphic,
   Medieval, Vedic, or a known ancient language (Latin, Gothic, Hittite,
   Tocharian, Sanskrit, Avestan, Akkadian, Luvian, Oscan, Umbrian, Gaulish,
   Pali, Prakrit, Sogdian, Khotanese, Parthian, Bactrian, Ge'ez).
5. **Deduplication by node identity.** Each usable row gets a node: its
   Glottolog MRCA (crown: the code itself; split A|B: MRCA of A and B; stem
   X: the parent of X). Rows with the same MRCA are one node. The primary row
   is the first by (two-sided interval > point > one-sided; grade A > B; file
   order). If a node has no two-sided row but has both a min and a max row,
   they are combined into an interval [min, max].
6. **Family** = Glottolog top-level family of the node. Nested or
   overlapping datasets therefore always stay with their family.

## Matching

- Node spec forms: `crown:X` (pairs of doculects under X whose lineages part
  at X), `split:A|B` (one doculect under A, one under B; B minus A), and
  `stem:X` (one under X, one under X's parent but not under X). Side lists
  may name several codes and `-code` exclusions.
- A (node, dataset) match needs at least one doculect on each side and a
  median of at least 40 shared concepts per cross pair.
- **Primary dataset.** For each family, the dataset that matches the most of
  the family's two-sided nodes (ties: larger median shared concepts, then
  name) is used for every node it matches. Other nodes use their own best
  dataset (largest median shared concepts, then more cross pairs, then
  name). Rationale: keep the coding protocol fixed within a family where
  possible, so τ is not inflated by protocol mixing. Sensitivity: the
  per-node "longest list" rule for every node.

## Coding

- Per (doculect, concept), the first form in `forms.csv` order that has a
  cognate judgement is used; per form, the first row in `cognates.csv`
  order. Synonyms beyond the first are ignored. The same rule everywhere.
- A form with no row in `cognates.csv` is treated as missing, not as a
  singleton. Datasets below 95% cognate coverage are flagged.
- Concepts are the dataset's own parameters (comparisons are within a
  dataset).

## Estimator

For a node and a dataset, over all cross pairs (one doculect per side) and
the concepts both members of a pair have:

- K_c, N_c = matches and comparisons for concept c, summed over cross pairs;
  p̂ = ΣK / ΣN (pooled). Under a strict clock every cross pair has the same
  expected match probability, so the pooled ratio is the composite
  (pairwise) maximum-likelihood estimate of that common p. The median of
  per-pair estimates is reported as a check.
- **T1, single rate**: d = −ln p̂ = 2rT, β = ln(d / 2T) (β = log rate per
  year, the scale of clepsydra's −8.8 default).
- **Gamma among-concept heterogeneity**: p = (1 + 2rT/α)^(−α), so
  2rT = α(p̂^(−1/α) − 1). α = 1 is primary (§2.1, Lomax); α = 0.5 and 2 as
  sensitivity. α is not tuned on the results.
- If p̂ = 1, one half non-match is added (flagged).
- **Lexical uncertainty**: bootstrap over concepts (1000 replicates, fixed
  seed), which keeps the correlation between cross pairs that share a
  concept. s_lex = SD of ln d over replicates.
- **Age uncertainty**: interval: Var(ln T) = (ln hi − ln lo)²/12; point:
  0.15². ln T is centred on the log midpoint (interval) or the point.
- Measurement variance s_i² = s_lex² + Var(ln T).

## Meta-analysis

b_i = μ + u_f + e_i, u_f ~ N(0, σ_f²), e_i ~ N(0, τ² + s_i²), family f(i).

- **REML** with profile likelihood intervals for σ_f (χ²₁ / 2 = 1.92 drop).
- **Bayesian** grid posterior: flat μ (integrated analytically),
  σ_f ~ HalfNormal(1), τ ~ HalfNormal(1); posterior median and 95% interval.
  Sensitivity: HalfNormal(0.5) on both (clepsydra's current scale).
- **F** = number of families with at least one two-sided node in the fit.
- The predictive SD for a node of a new family is
  √(σ_f² + τ² + Var μ̂ + s_lex²). τ counts: a new family contributes a new
  node, not a family mean.
- Depth: the slope of b_i on ln T is reported (not used for prediction).

## Censored secondary fit

Max rows give b ≥ ln(d / 2·max); min rows give b ≤ ln(d / 2·min). Fitted by
maximum likelihood with the same random-effects structure, integrating u_f
by Gauss-Hermite quadrature, together with the two-sided rows.

## Leave one family out

For each family with a two-sided node: drop all of its nodes (both
two-sided and one-sided), refit on the rest, and predict each dropped
two-sided node from its own observed d:

ln T_pred = ln d − ln 2 − β, β drawn from the predictive.

Predictors, all scored the same way:

1. **Clock, Bayesian**: mixture over the grid posterior of
   N(μ̂_θ, σ_f² + τ² + 1/A_θ), plus s_lex².
2. **Clock, REML plug-in**: Student t with F_train − 1 df, scale
   √(σ̂_f² + τ̂² + Var μ̂ + s_lex²).
3. **Old default clock** (baseline a): ln r ~ N(−8.8, 0.8²), P(match) =
   exp(−2rT), plus s_lex².
4. **Family size** (baseline b): ln T = a + b·ln(n), n = Glottolog
   languages under the node's MRCA (for a stem, under the parent); OLS on
   the training nodes, t predictive.
5. **Climatological** (baseline c): ln T ~ t(mean, SD·√(1+1/n)) of the
   training nodes' log midpoints.

Both clock predictors are run for T1 and for gamma α = 1.

Metrics, per held-out node, with the calibration as a density q on ln T:

- log score: ln ∫ p_pred(x) q(x) dx;
- expected coverage of the central 50/80/95% predictive intervals:
  ∫ 1[x in interval] q(x) dx;
- expected PIT: ∫ F_pred(x) q(x) dx;
- interval score (95%, Gneiting & Raftery) at the calibration's log centre,
  in log units;
- z = (x_centre − mean)/SD; Σz² reported against its expectation;
- sharpness: ×/÷ half-width factor of the 95% interval, and that factor
  for a 2000-year node with median s_lex.

Means are reported per node and family-weighted (each family's mean first).
One-sided held-out nodes get a censored log score, ln P_pred(T ≤ max) or
ln P_pred(T ≥ min), reported separately.

## Uniform automatic coding (if lingpy is present)

lingpy 2.6.13 is installed. For each primary (node, dataset) match:
the doculects in the cross pairs (at most 10 per side, first in
`languages.csv` order), the `Segments` column, and the concepts whose
Concepticon ID occurs in at least half of the primary datasets. LexStat
scorer with runs = 1000, `cluster(method='lexstat', threshold=0.55,
cluster_method='infomap')`, fixed random seed. Expert and automatic codings
are then compared on the same concept subset: σ_f and τ from both.

## Reading the result (decided now)

- σ_tot = √(σ_f² + τ²) is the number that sets the calendar interval of an
  uncalibrated node. ×/÷ 3.5 at 95% (the M2 sharpness criterion)
  corresponds to σ_tot ≈ 0.64.
- **Go** (for continuing the calendar track to M1/M2): the posterior median
  of σ_tot ≤ 0.5 under α = 1, LOFO 95% coverage ≥ 80% (family-weighted),
  and the clock beats both the climatological and the family-size
  baselines in family-weighted mean log score.
- **No-go**: posterior median of σ_tot > 0.64, or the clock does not beat
  the climatological baseline.
- **Inconclusive**: anything else, including any case with F < 6.

## Deviations

Made after the coverage check (which doculects fall on each side, and the
median shared concepts per cross pair), **before any match count, distance
or rate was looked at**:

1. **Two modern dialects** were caught by the ancient-name rule: "Middle
   Lozva Mansi" (bdpa, zhivlovobugrian) and "Middle Chulym"
   (savelyevturkic). They are kept by name.
2. **Rake nodes.** The Malayo-Polynesian crown has 35 children in Glottolog.
   For a crown node, a regional sample (e.g. Borneo only) gives cross pairs
   whose real MRCA is younger than the node. So the node-level choice of
   dataset now ranks first by the **number of Glottolog children of the
   crown that are represented** (split and stem nodes count 2), then by
   median shared concepts, then by number of cross pairs, then by name.
3. **Family-primary dataset** is used only when it matches at least two of
   the family's two-sided nodes; otherwise each node takes its own best
   dataset by the rule above. (With one node per dataset, "most nodes" is
   a tie and carries no information.)
4. **ABVD sensitivity.** `abvd` is local, expert-coded (291k cognate rows,
   84% of forms) but not flagged `ExpertCognates`. It is not in the primary
   analysis. A pre-declared sensitivity adds it as a candidate for
   Austronesian nodes under the same rules.

1. **Degenerate cognate codings excluded** (added after the first run). The
   first run gave σ_tot ≈ 1.4, driven almost entirely by one node: Maltese
   vs Maghrebi Arabic in `ratcliffearabic`, with 99.2% matches and an
   implied log rate of −12.5. That dataset's README says its cognate sets
   "have only been assigned with respect to Standard Arabic": it has 100
   sets for 100 concepts, so every pair of dialects "shares" every concept by
   construction. It carries no replacement signal. The rule, applied to
   every dataset alike, excludes a dataset whose mean number of distinct
   cognate sets per concept (over concepts coded in two or more doculects)
   is below 1.05. The first run's results are kept in `results_run1.tsv` for
   comparison.

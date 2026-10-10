# Adversarial review of `docs/roadmap.md` (clepsydra, drafted 2026-10-10)

Reviewer stance: hostile referee in Bayesian phylogenetics, computational statistics and historical
linguistics. The repository was read only and not modified. I checked the roadmap against
`docs/research/` (arawan-lessons, both reviews, both literature surveys, mathematical-analysis
with its `math/` scripts) and `reference/calibrations/calibrations.tsv`. The numbers that are not
quoted from those files come from small calculations I ran (binomial and χ² power, t-predictive
intervals; `scipy`, outside the repo).

---

## 0. Verdict in one paragraph

The diagnosis is right: no proper age prior, calendar scale silently borrowed, no validated clock.
The research is good, and R0–R1 are the right first moves. The problem is the **ordering after R1**.
The roadmap puts roughly three phases of engineering (R2 observation model, R3 sampler
performance) in front of the one experiment that decides whether the project's promise can be kept
at all (R4: is σ_f ≈ 0.3, or 0.5+, out of sample, on the families we actually have?). The gates
are strictly sequential, so a 5× ESS/s performance target (Gate R3) blocks the scientific
centrepiece. Meanwhile the slowest task, R4.1 data curation, cannot start until R3 is done. Three
of the scientific claims also need tightening before anyone writes code:

1. The convolution identity (R1.4) is inconsistent with the default prior chosen in R1.1.
2. The R4 gate as written cannot fail: it either passes or recalibrates, and with ~8–10
   independent families it has roughly 20–35% power to detect 95% intervals that really cover 80%.
3. The calibration database supports fewer independent training families than the "≥ 10" that
   R4.1 needs, and it contains rows that would leak or be circular.

The Arawan promise (~1000–4000 BP for a 2000 BP median) is the **best case** from
`mathematical-analysis.md` §3.3. It assumes σ_f = 0.3 is *known*. A realistic estimate from ~8
families gives something like **×/÷ 2.3–2.9, i.e. ~700–5700 BP**, unless Arawan-specific soft
calibrations are added. The roadmap does not contain a single item that produces those.

---

## 1. Ordering and the critical path

### 1.1 What is actually on the critical path to the owner's goal

The goal is calibrated calendar posteriors for an uncalibrated ~7-language family. By
`mathematical-analysis.md` §3.2–3.3, the calendar width there is set ~80–90% by σ_f (root log-SD
0.33 total against 0.07 from data under a strict clock, or 0.13 under a relaxed clock). Hence:

- **Critical:** a proper prior (R1.1–R1.3), honest scale reporting (R1.4–R1.6), a curated and
  leak-free training set (R4.1), a frozen dating model (R4.2), the meta-analysis (R4.3), the clock
  package (R4.4), and LOFO (R4.5). Plus the minimum observation model without which β_f is
  biased: missing data (R0.1), concept heterogeneity (R2.1), sampled ε (R2.6), synonyms (R2.3),
  and calibration semantics in the *training* families (R0.3–R0.5, R6.1, R6.2; see 1.3).
- **Second order for the Arawan root:** R2.2 (branch×concept Γ, worth ~5% at 2 ka per §2.4(d)),
  R2.5 (tip offset, worth +80–180 y at most, §8), R5 (the root is the MRCA of every taxon set, so
  its age hardly depends on topology), R3.4–R3.5 (speed and checkpoints).

### 1.2 Specific misplacements

| Problem | Items | Fix |
|---|---|---|
| R4.1 (match 59 Lexibank datasets to calibrations, Glottolog/Concepticon mapping, dedup, protocol metadata) is the longest-lead task, and it is pure data work. It waits for Gates R1–R3. | R4.1 | Start it **now**, in parallel with R0. It needs none of the C changes. |
| The σ_f go/no-go information comes last. If out-of-sample σ_f is 0.5–0.6 (plausible: Rama & Wichmann 2020's GBD LOO MAE of 770 y barely beats the "number of languages" baseline at 822 y, `literature-lexical-dating.md` Part 3), the Arawan deliverable and the paper framing change completely. | R4.3/R4.5 | Add **R4.0, a σ_f pilot**, right after R0. Use fast estimators (pairwise ML distances §5b, or the current T1 on fixed published trees with R0 fixes and a crude root prior) on the calibrated families. Fit `b_f ~ N(μ, σ_f² + s_f²)` and run a crude LOFO. This takes days in Python and sets expectations for everything else. |
| Gate R1(a) is "SBC, uncalibrated root with a clock prior". The ρ·a ridge will mix badly with the current single-site sampler, so the gate will fail for sampler reasons. | R3.1–R3.3 vs Gate R1 | Move R3.1 (ratio + log-root), R3.2 (clock-preserving move) and R3.3 (state-relative steps) **into R1**. Or better (1.4 below): for uncalibrated runs, sample in lexical time, where the ridge does not exist. |
| Gate R3 (5× ESS/s) blocks R4. Performance is not a scientific prerequisite. | Gate R3 | Make R3.4/R3.5 "as needed", not gating. |
| R6.1 (outlier model) and R6.2 (attestation-lag / one-sided) come *after* R4, yet the R4 training calibrations are mostly "point (max-leaning)" political or historical events and min-only bounds (see 3.2). Fitting β_f to these as Normal points biases β_f and inflates σ_f. | R6.1, R6.2 | Do them **before R4.2**. They are small (literature-dating-methods Part C rec. 3 rates them S–M). |
| R6.3 (sampled ancestors): training families have 118 A/B `tip` rows (Latin, Old Norse, Old Irish …). Dating ancestor-like tips as sisters makes trees too old (Chang et al. 2015; lit-lexical bias table), which biases β_f. Arawan's Paumari-1833 is itself essentially an ancestor of modern Paumari. | R6.3 | Before R4: exclude ancestor-like tips from training, or treat them as ancestor calibrations (zero-length branch). Full FBD can wait. |
| R4.2 says "fit the same model used to date", but R2 (and later R6, R7) keep changing defaults. Every default change invalidates the package. | R2 vs R4.2 | **Freeze a minimal default model** (T1 + concept Γ + sampled ε + synonym marginalisation + missing handling + strict/relaxed choice) before R4.2. Put R2.2, R2.5 and other mixing families in as sensitivity variants with their own refitted package, not as new defaults. |
| R7.2 (misspecified "world" simulator) is placed after R4, but it is what tells you, *before* reading real LOFO numbers, what coverage and power to expect. | R7.2 | Run it alongside R4 and use it to fix the R4 gate's decision rule (section 3). |
| R7.5 (identical likelihoods to Rönchen et al.'s T1 code) is cheap and checks the foundation of every number. | R7.5 | Pull it into R0/R1. |

### 1.3 Gating philosophy

"No phase starts before the previous gate passes" is right for correctness (R0→R1) and wrong after
that. Use two tracks: **engine** (R0→R1→minimal R2→R3 as needed) and **data/validation**
(R4.1 curation, R4.0 pilot, R7.2 simulator, Arawan calibration curation). The tracks meet at R4.2.

---

## 2. Items that are scientifically wrong or likely to fail

### 2.1 R1.1 default prior vs the R1.4 convolution identity (they contradict each other)

- `mathematical-analysis.md` §4.1: the identity *log a = log(ρa) − log ρ, with the lexical
  posterior independent of the rate prior* is exact only if the prior on ages is
  **scale-invariant** and ρ is a priori independent of it.
- R1.1's default is "an explicit root-age prior × uniform-conditional". A proper calendar root
  prior (lognormal, or a uniform window like Arawan's U(300, 3000)) is *not* scale-invariant. With
  it, the calendar posterior is the product of two informative scale priors (root window × clock
  package), not a convolution. That double counting is exactly how Arawan's root got glued to
  3000 (A2).
- The obvious scale-invariant choice, p(a_root) ∝ 1/a_root, gives an **improper** posterior when
  the contamination floor is on. As lexical length goes to ∞ the likelihood tends to a positive
  constant (codex finding 2), and ∫ dL/L diverges. R1.2 would then have to refuse the very
  configuration R1.4 relies on.

**Correct construction.** Put a proper, broad prior on the **lexical** root age (e.g. a lognormal
on ρ·a_root in expected-replacements units), uniform-conditional (or BD) on the ratios, and an
independent clock prior on ρ. Calendar ages are then *defined* as the convolution, which is exact
and proper. Calendar-unit root windows supplied by a user are treated as **calibrations**, with an
outlier model (R6.1) and an effective-prior display (R1.5), never as "the root prior".

The identity **breaks** under:

- any calendar calibration, including internal nodes, root windows and **dated tips**. Arawan has
  two dated tips, so the identity is only approximate for the motivating case;
- BD/Yule with an informative λ prior in calendar units. It holds if λ has a 1/λ prior and the
  turnover and sampling fraction are unitless;
- concept-rate or branch-rate priors whose *means* are in absolute units. Relative effects are
  fine, and so is a relaxed clock with mean-1 multipliers: the identity survives a relaxed clock;
- covariate models that put calendar time into the rate, e.g. epoch or time-dependent clocks.

R1.4 should state this list. The engine should keep a joint calendar-time path for those cases and
use convolution only where it is exact. A self-check is cheap: importance-reweight the lexical
draws by the calendar-prior ratio and report a Pareto-k̂.

### 2.2 Uniform-conditional prior: three technical traps

1. **Dated tips.** With non-zero tip ages, the ordered region given the root is not
   `c(τ)·root^(m−1)`. Nodes ancestral to dated tips have lower limits, so the volume is a
   piecewise polynomial in the root. Dividing by `root^(m−1)` does *not* make the root marginal
   equal the stated root prior, and Gate R1(c) (ancient tips) will fail SBC by a small but
   detectable amount. Compute the exact normaliser by a post-order polynomial recursion (cheap on
   a fixed tree), or use MCMCTree-style conditional construction.
2. **Internal calibrations.** Multiplying calibration densities into the uniform-conditional
   density gives effective priors that differ from what the user wrote (Heled & Drummond 2012;
   lit-dating Part B §2). R1.5 *displays* this, which is good. Say explicitly whether the target is
   "display" or "condition so that marginals match" (calibrated-Yule style). The current text is
   ambiguous.
3. **Variable topology (R5).** Uniform-conditional on each topology times a uniform topology prior
   is not a coherent generative prior. The normaliser c(τ) (the share of ranked histories) must
   enter every per-topology evidence, or topologies get silently reweighted (§5j already warns
   about this; R5.3 does not inherit the warning).

For small, almost completely sampled families (Arawan: ~7 of ~8 languages), BD with sampling
fraction ≈ 1 is cleaner than a uniform-conditional prior and is a coherent topology prior.
Relative ages under a relaxed clock are prior-sensitive (§3.2: ratio SD 0.16 → 0.27), so a
uniform-conditional vs BD sensitivity panel belongs in R1, not "later".

### 2.3 Discrete-gamma concept rates (R2.1)

- Identification: a 7-language family cannot estimate the shape (§2.5: 6 node depths do not
  identify curvature). In Arawan the hyperprior on the SD will be prior-dominated, so its prior
  **must come from the package** (R4.3), fitted under the same parameterisation.
- **Tails:** deep nodes are dated by the slow tail (§3.1, §2.3). Yang's equal-probability K = 4
  with category means or medians truncates exactly that tail. Use K ≥ 8 with mean-matched
  categories, or generalised Gauss–Laguerre for gamma and Gauss–Hermite for lognormal. Check
  against high-K quadrature on the deepest node of each training family.
- **Normalisation:** the concept-rate prior must be centred (mean-1 or sum-to-zero on the log
  scale), with all mean-rate uncertainty in ρ. Otherwise the prior variance of the effective mean
  rate shrinks like 1/C (dos Reis et al. 2014; lit-dating Part C rec. 4). The roadmap never says
  this, and it is precisely what makes β_f a transportable quantity. Note also that with mean-1
  normalisation, β_f is the mean rate *of that concept list*. Transport to KAHD's 253 concepts
  needs the Concepticon concept effects of R4.3, evaluated on the target list.

### 2.4 Lomax branch×concept heterogeneity (R2.2)

- Exactness: within one branch, Lomax `(1+rt/k)^(−k)` is the exact keep-probability for a gamma
  multiplier drawn independently per branch and per concept, and the recursion is unchanged.
  That part is right.
- **It is not consistent under branch subdivision.** Two segments t₁ + t₂ with independent draws
  give a different retention from one branch of length t₁ + t₂. So the implied process depends on
  how many nodes the tree happens to have: a **node-density artefact**. Paths through densely
  sampled trees (IE-CoR, ABVD) average more multipliers and look more exponential. Arawan's long,
  sparse branches look more heavy-tailed. This is a **transport hazard for β_f**, because the
  training trees are dense and the target is sparse. It also breaks under sampled ancestors (R6.3),
  where a degree-2 node splits a branch.
- If you keep it, implement a subdivision-consistent version: a Markov-modulated (covarion-like)
  concept rate that switches categories at rate s, so the state is (class, category); or a
  per-lineage-per-concept multiplier inherited down the tree (autocorrelated). Failing that,
  **subsample training families to target-like sizes** (see 3.4).
- Payoff for Arawan: about 5% at a 2 ka root (§2.4(d) table, CV 0.5). This makes it a sensitivity
  variant, not a default. Its real value is in deep training families, not the target.

### 2.5 Per-language coding error (R2.5): identifiability

- Under T1 a split error is *exactly* a terminal-branch replacement (§8). With a relaxed clock, a
  per-tip offset δ_i is **completely confounded** with that terminal branch's multiplier m_b, so
  it is identified only by priors. Under a strict clock it is identified only through departures
  from ultrametricity, which are also what a relaxed clock explains.
- So everything depends on where the prior on δ comes from, and the roadmap does not say. Two
  empirical sources exist in the local data:
  1. **Same language, different collector** doculect pairs in Lexibank/KAHD. These have time ≈ 0,
     so their lexical distance is pure elicitation, transcription and coding error.
  2. **Known ancestor–descendant pairs** (Paumari 1833 → modern; old vs modern lists elsewhere).
- Hold Paumari-1833 out of the prior fit so it remains a *test*. (A9 is the only external check
  Arawan has. Don't train on it.)
- Payoff for the Arawan **root**: small (+80–180 y on lexical length, §8). The two old lists have
  41 and 70 concepts and bring essentially no scale information: about 2 expected replacements in
  170 years, swamped by noise. The simplest honest option for root dating is to run with and
  without the old lists. R2.5 matters for shallow nodes and pair dating, not the headline.

### 2.6 Exhaustive topology enumeration with Laplace evidences (R5.3)

This is expensive and fragile, with low payoff for this goal:

- KAHD has **9** doculects (7 modern + Arua + Paumari-1833), which is 2,027,025 rooted binary
  topologies, not 10,395. At §5j's 10–50 ms each that is 6–28 h per model variant, before
  sampled-ancestor placements.
- Each topology needs a mode search and Hessian over ages, 2L−2 log-multipliers, ε, concept-Γ
  shape and σ_b, with INLA-style nesting over the hyperparameters. Near-polytomies, which §5j says
  are common, put modes on boundaries where Laplace is biased *in a topology-dependent way*. For
  BMA, a 1-nat relative error multiplies a weight by e.
- The payoff for the root age is near zero (it is the MRCA of all taxa). For internal named sets,
  R5.1 (date over the existing BEAST/yvyra tree sample) and R5.2 (topology MCMC) deliver the same
  thing with standard, testable machinery.
- **Cut R5.3 from the C roadmap.** Keep an *oracle* in Python for L ≤ 6 (945 topologies): MCMC per
  topology plus stepping-stone, used only to validate R5.2. Reorder R5 as R5.1 → R5.4 → R5.2.

### 2.7 The cross-family meta-analysis (R4.3): is β_f transportable?

β_f is not a property of a family. It is a property of (family × concept list × coding protocol ×
model × calibration depth × sampling density). The roadmap handles the first and the model
(R4.2, R4.4). The rest are not handled:

| Source of non-transport | Status in roadmap | What to do |
|---|---|---|
| Concept list | Concepticon pooling (R4.3); refuse mismatched lists (R4.4) | Good. But "refuse a mismatched list" plus heterogeneous training lists (IE-CoR 170, ABVD 210, DravLex 100, KAHD 253) contradicts itself. Evaluate μ + mean concept effect *on the target list*, and report the overlap-induced variance. |
| Coding protocol (root vs word cognacy, lumper/splitter, loans, synonyms) | Only "coding protocol" in the manifest | Protocol is **perfectly confounded with dataset and team**, so it lands in σ_f. That is honest only if the target's protocol is exchangeable with the training protocols. Record protocol and team per dataset, model team as a random effect (or leave-one-team-out), and estimate the protocol effect directly: double-code a KAHD subset under the IE-CoR protocol, or the reverse. |
| Calibration depth (time dependence from a mis-specified mixing tail, §2.4(b)) | "bias vs depth" reported in R4.5 | Make depth a *covariate* in R4.3 (lit-dating Part C rec. 4 says so). Most A/B node calibrations are under 2.5 ka. Arawan's root (~1.5–3 ka) lies inside that band, which is good news: depth-matched training is achievable. Say so and exploit it. |
| Sampling density / node count (bursts §10, Lomax subdivision 2.4) | Absent | Subsample training families to 6–9-language subtrees (3.4). |
| Population / sociolinguistic shift (Bromham 2015 vs Greenhill 2018) | Covariate "tested by LOFO" | Fine, but covariate selection must be **nested** inside LOFO (3.1). |
| A uniform protocol as a transport device | "Automatic cognate detection out of scope" | Reconsider *as an experiment*. ASJP's edge in Rama & Wichmann (MAE 606 vs 695 GBD, 728 glottochronology) may come partly from applying one protocol everywhere. Fit σ_f under expert codings and under one automatic coding (LexStat/cognator on the same Concepticon core) across all families including KAHD. If the automatic σ_f is smaller, transport through it. |

**Two-stage normal approximation:** `b_f ~ N(μ, σ_f² + s_f²)` is poor for families with one
min-only or max-leaning calibration, where p(β_f | D_f) is one-sided. Carry the family posteriors
as samples (or as a likelihood-ratio grid in β), not as (b_f, s_f).

**Small-F hyperparameter uncertainty:** with F = 8, the 95% CI for σ_f is (0.66, 2.04)·σ̂_f.
The predictive is a t₇ distribution, not a normal. This is not in the "≈ 0.3–0.4" story. The
package must ship the *predictive distribution* (draws), not a point σ_f.

### 2.8 LOFO with ~10 families (R4.5): power

- Within-family nodes share the β error, so the effective sample size for coverage is the number
  of **families**, not nodes.
- Binomial test of nominal 95% coverage, α = 0.05:
  - F = 8: reject if ≤ 5 covered. Power against true coverage 80% is **0.20**.
  - F = 10: reject if ≤ 7. Power **0.32** against 80%, **0.18** against 85%.
  - F = 15: 0.60. F = 30: 0.74.
- Nominal 80%: F = 10 has power 0.37 against a true 60%.
- A sharper statistic, the variance of the standardised log errors (Σz² ~ χ²_F), detects the
  case "95% intervals actually cover 80%" (variance ratio ≈ 2.3) with power **0.58 at F = 8 and
  0.65 at F = 10**. Better, but not decisive.
- Worse, σ_f is *fitted on the same F families*. In-distribution LOFO coverage is therefore
  roughly right **by construction**. LOFO mainly tests exchangeability among the training
  families. It **cannot** test the step that matters, transport to a small, oral, Amazonian,
  KAHD-coded family.

---

## 3. Is the R4 headline validation able to fail? Leakage

### 3.1 As written, it cannot fail

Gate R4 reads: "coverage within sampling error of nominal; if not, recalibrate and report the
inflation". With the power above, "within sampling error" almost always holds, and when it does
not, the remedy is automatic. Neither branch is a failure. Pre-register falsifiable criteria
instead:

1. **Skill against trivial baselines** (proper scores: log score and interval score of the
   *calibration density* under the predictive, see 3.3). Baselines:
   - (a) the number of languages / Glottolog family size (Rama & Wichmann: 822 y MAE);
   - (b) the Holman/ASJP formula;
   - (c) a "climatological" predictive, i.e. the empirical distribution of all A/B node ages given
     node depth rank;
   - (d) the 0.2.0 default Swadesh rate prior.

   **Fail** if clepsydra does not beat (a) and (c) in mean log score with a family-clustered
   bootstrap CI excluding 0. In that case the tool reports only lexical and relative ages for
   uncalibrated families. That outcome is a legitimate, publishable result.
2. **Sharpness ceiling:** fail if the LOFO-calibrated 95% interval for a 2 ka node is wider than
   (say) ×/÷ 3.5. At that width the calendar output is not decision-relevant, and the docs must
   say so.
3. **Dispersion test:** Σz² against χ²_F, plus a PIT ECDF with simultaneous bands (Säilynoja et
   al. 2022). Report the posterior of the inflation factor, not a pass/fail.
4. **Depth slope:** regress log(â/a) on depth. §2.4 predicts a negative slope if the mixing tail is
   wrong. Fail if the slope's CI excludes 0 *in the depth band of the target*.

Pre-register these, together with family list, grade rules, model and baselines, **before**
computing any real LOFO numbers. Any model choice made after looking at LOFO (mixing family,
covariates, grade cut-offs, which families to drop) needs **nested** LOFO, or it optimistically
biases the headline.

### 3.2 Leakage and circularity found in `calibrations.tsv`

My tabulation: 367 rows. Grades A = 114, B = 84, C = 151, U = 18. 45 families have some A/B row.
But A/B **internal-node/root** calibrations are concentrated in IE (~16), Austronesian (~16 before
dedup), Arawakan (8), Uralic (3), Afroasiatic (4), Eskimo-Aleut (3–4 near-duplicates of one Inuit
node), and a scatter of single nodes. Concrete problems:

- **Duplicates inflate weight.** The same node appears more than once across `oldworld.tsv` and
  `newworld_pacific.tsv`, or as node + settlement_bound: Ma'anyan–Malagasy (×2), Eastern/Western
  Cham (×2), Northern Roglai–Tsat (×2), Malayo-Polynesian root (×2), Proto-Oceanic (node +
  settlement), Proto-Central Pacific (node + settlement), Inuit (×3 variants), Cholan (×2),
  Afrikaans–Dutch (×2, once under "Germanic (Low Franconian)"). R4.1's "one row per node" has to be
  enforced by node identity (MRCA glottocode set), not by row id.
- **Circular "calibrations".** The Pama-Nyungan root (B, "Bayesian estimate") and the
  Proto-Uto-Aztecan root (B, "Bayesian inference … benchmark") are outputs of lexical dating.
  Using them in training or testing measures agreement with BEAST, not accuracy. Demote them to C.
- **Nested families as separate "families".** "Romance (French)" and "Germanic (Low Franconian)"
  sit inside IE. Lexibank has `meloniromance`, `saenkoromance`, `walworthpolynesian` and
  `servamalagasy` nested inside IE-CoR/ABVD families. Creoles ("English-lexified creole",
  "French-lexified creole", "Suriname creoles") are derived from IE lexifiers, and their lexical
  history is creolisation/relexification, not tree-like replacement. **Exclude creoles from
  training entirely**, and group nested datasets with their parent for LOFO.
- **Max-leaning points.** Many A/B rows are political events that bound a split from one side
  (e.g. IE-ROM-01 "point (as published); better treated as max-leaning"; East Slavic, Goidelic,
  Romani, Cham, Roglai, Oromo, Mongolic, Tungusic). Scoring them as two-sided points miscalibrates
  the *test*, not the method. Use censored scoring for one-sided bounds.
- **Archaeological correlations for Arawakan** (ceramic traditions → clades, SA-001…SA-011, mostly
  B). Check whether these associations were themselves built with lexicostatistical dates. If so,
  they are circular. They are also among the few South American calibrations, so they deserve
  extra scrutiny, not less.
- **Same coders across families.** Several South American datasets (e.g. `gerarditupi`, `kahd`,
  Chacón's Arawakan sets) share authors or coding teams. Shared protocol is good for transport but
  correlates LOFO errors. Treat team as a clustering unit (leave-one-team-out as a robustness
  check).
- **Everything learned has to be refit inside the loop:** concept effects, mixing shape, k, σ_b
  hyperprior, ε prior, covariate coefficients, and the calibration support/upper bound. Codex C30
  found that the current LOCO computes support from all rows, including the held-out one. R4.5
  should say "every data-derived quantity is refit without family g", not only "the
  hyperparameters".

**Realistic F.** Cross the A/B internal nodes with the local expert-cognate Lexibank sets and
remove nesting, creoles, circular rows and min-only singletons. What is left is roughly IE
(`iecor`), Austronesian (`abvd`), Uralic (`uralex`), Arawakan (`walkerarawakan`/`chacon*`),
Semitic/Arabic (`kitchensemitic`, `ratcliffearabic`), Dravidian (`dravlex`, min bounds only),
Tungusic (`oskolskayatungusic`, one 186 y node), Ainu (`leeainu`, one node), and perhaps
Pama-Nyungan (`bowernpny`, Wati 1500–5000, almost uninformative). That is **~6–9 families, several
with one shallow node**. The "≥ 10 genuinely independent families" target in R4.1 is probably not
reachable at grade A/B. Decide now whether to (i) admit C-grade rows with an outlier model (R6.1)
and a higher outlier prior, (ii) add tip-calibrated families (Tupí-Guaraní with Old Tupinambá/Old
Guaraní, a regionally apt one), or (iii) accept F ≈ 7 and the wider t-predictive.

### 3.3 How to score against an uncertain truth

A calibration is a distribution, not a number. PIT against a point "truth" mixes up calibration
error with method error. Score by ∫ p_pred(a) p_cal(a) da (log score of the calibration under the
predictive), or compute PIT against draws from p_cal. Use censored likelihoods for min/max-only
rows.

### 3.4 Make LOFO look like the target

Dating a held-out 100-language IE tree "as if uncalibrated" is not the Arawan problem. For each
held-out family, also evaluate on **random 6–9-language subsamples** that contain the calibrated
node(s), with **topology not given** (R5.1/R5.2) if that is how Arawan will be run. This gives:

- the right regime: sparse sampling, long branches, few depths, relaxed-clock floor;
- a check on the node-density and Lomax transport hazards (2.4);
- more (correlated) test cases per family, though not more independent units.

Report a **stratified** result for the Americas / small oral families (Arawakan, Tupí-Guaraní,
Mayan if data can be found), even if it is only anecdotal. Out-of-region transport is the
assumption the Arawan number rests on, and the report should say it was not validated.

---

## 4. Arawan: what would make the outcome concretely better, and is ~1000–4000 BP credible?

### 4.1 Credibility of the promise

`arawan-lessons.md` and §3.3 give ~1050–3800 BP (×/÷ 1.9) for σ_f = 0.33 **known**, a correctly
specified model, and a perfectly transportable clock. Adding what the roadmap itself implies:

| Assumption | 95% for a 2000 BP median |
|---|---|
| σ_f = 0.30 known (the promise) | 1110–3600 (×/÷ 1.8) |
| σ_f = 0.33 estimated from F = 8 (t₇ predictive) | 870–4580 (×/÷ 2.3) |
| … + relaxed-clock data floor 0.13 | 830–4840 (×/÷ 2.4) |
| σ_f = 0.40, F = 8, + 0.13 | 700–5710 (×/÷ 2.9) |
| σ_f = 0.50, F = 8, + 0.13 | 550–7270 (×/÷ 3.6) |

Holman's 29% mean absolute discrepancy is **in-sample**, which corresponds to a log-SD of ~0.36.
Out of sample it should be larger. Rama & Wichmann's LOO RMSE of 936 y on ages of mostly 1–6 ka
is a relative error of ~40%+. Then add model uncertainty (mixing family, §11), protocol shift, and
the unvalidated Amazonian transport. **The honest expectation is ×/÷ 2.3–3, i.e. roughly
700–5500 BP.** The roadmap should promise the method (an LOFO-backed interval whose sources are
shown), not the width.

### 4.2 What would actually narrow it

In decreasing value per effort:

1. **Arawan-specific soft calibrations with an outlier model.** This lever is absent from the
   roadmap and is the only one that beats the σ_f floor. The database has only ARW-003/005
   (grade C) and ARW-004 (U) for the root. Candidates:
   - loanword strata: Tupí-Guaraní/Língua Geral and Portuguese loans give upper bounds on late
     splits (lit-lexical idea 14; Maurits et al. 2020);
   - contact loans between Proto-Arawan and Arawakan, especially Purus Arawakan
     (`carvalhopurus` is local);
   - *Wörter und Sachen* on reconstructible cultigens and technology;
   - middle-Purus archaeology as a soft maximum.

   Even a one-sided soft maximum (e.g. root < 4000 BP, 5% tail mass) cuts the long upper tail of a
   lognormal predictive considerably. Add an item: **R6.5, curate and grade Arawan soft
   calibrations, with an expert (Dienst, Dixon, Ferraz Gerardi), before R8.1.**
2. **Restricting transport to depth-matched, target-like training data** (2.7, 3.4). This reduces
   bias more than variance, but bias is what destroys coverage.
3. **Estimating the protocol effect for KAHD** by double-coding ~50 concepts under the training
   protocol(s). This turns an unknown shift into a measured one with an SE.
4. **Covariates that LOFO actually supports** (population size first). This may shift Arawan rather
   than narrow it (Bromham 2015 predicts faster loss in small populations, so younger dates).
   Report either outcome.
5. **The cheap things that are already planned:** sample in lexical time, show the three kinds of
   age, power-scaling sensitivity on the rate prior (2.9 below). These do not narrow the interval,
   but they make it defensible.

Things that will **not** materially narrow the Arawan root interval: R2.2, R2.5, R5.2/R5.3, R3.x,
R7.3, more concepts (§3.3), more training families beyond ~10 (§4.2: 0.315 → 0.305).

---

## 5. Scope and effort: cut, defer, or demote

| Item | Recommendation | Reason |
|---|---|---|
| R5.3 exhaustive enumeration + Laplace/AGHQ in C | **Cut** (Python oracle for L ≤ 6 only) | 2.6: 2 M topologies for KAHD's 9 doculects, boundary-biased Laplace, near-zero payoff for the root |
| R5.2 native topology MCMC | **Defer** until after R4/R8.1; do R5.1 first | Arawan already has BEAST/yvyra posteriors; R5.1 + R5.4 cover named-set ages |
| R7.4 completing T2 | **Mark experimental now; don't complete** | King 2026a and Rönchen favour T1; T2 needs its own package and validation |
| R7.6 decay-law priors | **Label as illustrations** | Off the critical path |
| R7.3 Godambe composite engine in C | **Prototype in Python for the R4.0 pilot**; C later if at all | Most of its value is as a fast estimator for σ_f screening |
| R3.4 caching, R3.5 hashed checkpoints, R3 gate | **Non-gating; do when profiling demands** | Engineering, not science |
| R1.1 BD with hyperpriors | Keep uniform-conditional + BD(sampling fraction ≈ 1) as a sensitivity pair; defer full hyperprior machinery | Small families: sensitivity matters more than richness |
| R2.2 Lomax | **Sensitivity variant**, not default; fix subdivision consistency first | 2.4 |
| R2.4 loan mixture | Keep `exclude` + read the columns; defer `mixture` | KAHD loan columns are empty (A14) |
| R6.3 FBD | Defer; meanwhile exclude ancestor-like tips from training and allow zero-length ancestor branches | 1.2 |
| R0.14 docs | Do it, but last within R0 | Cheap, not urgent |

---

## 6. What the research recommends that the roadmap leaves out

1. **Power-scaling prior sensitivity** (Kallioinen et al. 2024). This is `mathematical-analysis.md`
   §6.4 and ranked recommendation 8, and it does not appear in the roadmap. It is cheap
   post-processing on stored draws and gives exactly the "prior-dominated, no conflict" diagnosis
   Arawan needs. Put it in R1.5/R1.6.
2. **Variance decomposition** (σ_f vs data vs branch-rate vs model): §11, §3.3 corollary,
   lit-lexical idea 2. R1.4 reports three kinds of age but not *where the width comes from*.
3. **Model-averaged deep ages under ≥ 2 mixing families** (§2.5 practical rule; §11 law of total
   variance), with a printed model-uncertainty term. The roadmap only mentions a gamma *or*
   lognormal option.
4. **Concept-rate centring / normalisation** (dos Reis et al. 2014; lit-dating rec. 4). See 2.3.
5. **Depth (time-dependence) covariate in the cross-family model** (lit-dating Part B §5 and
   rec. 4/11). Shallow-calibrated rates applied to deeper targets are "the most likely source of a
   large systematic error".
6. **Posterior SBC on the user's own dataset** (King 2026b; lit-lexical idea 19). Cheap with a
   fixed tree, and directly addresses "is *my* Arawan posterior calibrated under the model".
7. **SBC with data-dependent test quantities** (Modrák et al. 2025: the joint log-likelihood as a
   test quantity). Without it, Gates R1/R7.1 cannot detect "posterior = prior", which is exactly
   the Arawan failure mode.
8. **Calibration-justification schema** (Parham et al. 2012; lit-dating rec. 24), feeding R6.1
   outlier priors. Also **uncertain target-node mixtures** (Rutschmann 2007, Pyron 2010). Several
   archaeological correlations (Arawakan ceramics, Indo-Iranian) are uncertain about *which node*
   they date.
9. **Strict-vs-relaxed spike-and-slab** (lit-lexical idea 21). The Arawan paper rejected a strict
   clock. With 7 tips, P(strict | data) is interpretable and decides whether the 0.13 relaxed floor
   applies.
10. **Two-component borrowing (internal vs external)** (lit-lexical idea 8). External loans should
    be *missing*, internal loans bias young. This is relevant to Amazonia, even if KAHD's loan
    columns are empty today.
11. **Polytomy / near-zero branch handling** (lit-lexical idea 24; §5j). Small Amazonian families
    often have rapid radiations.
12. **Bootstrap over meanings as a misspecification check** (§12). If the bootstrap SD is much
    larger than the model's data-only SD, meanings are more heterogeneous than modelled.

---

## 7. Revised ordering

**Track A: engine**

- **A0 (= R0, plus R7.5).** Correctness fixes. Likelihood identity against Rönchen's T1 code.
- **A1 (= R1, revised).**
  - Proper prior on **lexical** root age, uniform-conditional/BD on ratios, clock prior on ρ.
    Calendar-unit root windows become calibrations.
  - Exact normaliser with dated tips.
  - Prior-only mode, effective-prior and power-scaling report, variance decomposition, three
    kinds of age.
  - A lexical-time sampler for uncalibrated runs; ratio + log-root parameterisation and the
    clock-preserving move for calibrated runs (R3.1–R3.3 moved here).
  - Gate: full-pipeline SBC incl. ancient tips, with a log-likelihood test quantity.
- **A2 (minimal R2, frozen as the package model).**
  - R2.1 with centred rates and tail-faithful quadrature; R2.6; R2.3; R2.7.
  - R6.1 + R6.2 (needed by training).
  - Ancestor-tip handling (minimal R6.3).
  - PPC gate (as in R2).
- **A3 (as needed).** R3.4–R3.6.

**Track B: data and validation (starts on day 1)**

- **B0.** Pre-registered R4 analysis plan: families, grade rules, dedup by MRCA identity,
  exclusions (creoles, circular Bayesian rows, nested datasets), metrics, baselines, falsifiable
  criteria (3.1).
- **B1 (= R4.1).** Curation: Lexibank ↔ Glottolog ↔ Concepticon ↔ calibration node mapping;
  protocol/team metadata; same-language-different-collector pairs for the δ prior.
- **B2 (new R4.0).** σ_f pilot with fast estimators. Decide whether the calendar promise is
  viable, and whether a uniform automatic coding transports better than expert codings.
- **B3 (= R7.2, early).** World simulator to calibrate the power and decision rules of the R4
  gate.
- **B4 (new R6.5).** Arawan soft-calibration curation with experts.

**Merge**

- **M1 (= R4.2–R4.4).** Fit on the frozen A2 model. Package ships predictive *draws*, the concept
  effects, depth and covariate terms, protocol and team metadata.
- **M2 (= R4.5, revised).**
  - Nested LOFO.
  - Target-like subsampled evaluation, topology unknown.
  - Censored and density-based scoring.
  - Baselines.
  - Americas stratum.
  - Pre-registered pass/fail.
- **M3 (= R8.1).** Arawan:
  - lexical and relative chronology;
  - calendar ages under the package ± Arawan soft calibrations;
  - with and without the old lists;
  - R5.1 tree-set integration for internal nodes;
  - Paumari-1833 kept as the held-out check.
- **M4 (= R8.2–R8.4).** Paper and release. Sensitivity variants (R2.2 consistent version, R2.5,
  BD vs uniform-conditional, mixing families) appear as a model-uncertainty table.

**Deferred:** R5.2, R6.3 full FBD, R2.4 mixture, R7.3 in C, R7.4 completion.
**Cut:** R5.3 in C, R7.6 as inference.

---

## 8. Smaller points

- R1.6 should also flag *calendar root windows supplied by the user*. Arawan's U(300, 3000) was
  "copied from the BEAST prior". The tool should warn that a window without a justification row is
  treated as an outlier-prone calibration.
- R0.9: also add a gate on **prior-sensitivity** (power-scaling) for calendar outputs, not only on
  convergence. A converged posterior that equals the prior should be flagged (the Polynesian root
  in K§4.1 is exactly this).
- R2.6: with 7 languages ε is weakly identified. Give it an empirically based Beta prior from the
  package, or the improperness issue (codex 2) comes back through the ε → 1 corner.
- R4.4 "warns on low concept overlap": quantify it as the predictive variance added by
  unobserved concept effects, and print it in the variance decomposition.
- R8.4: the paper's core claim should be the *pre-registered* R4.5 table **including a failure
  outcome**, e.g. "uncalibrated calendar dating does not beat the node-count baseline for small
  families". Otherwise the headline is unfalsifiable.

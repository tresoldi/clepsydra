# Adversarial review of `docs/roadmap.md`

## Verdict

The roadmap has the right diagnosis of the current software—R0 and a proper time-tree prior are genuine release blockers—but the scientific critical path is wrong after R1. The project’s central claim is R4, yet the calibration model, topology uncertainty, model-adequacy tests, external likelihood checks, and the validation design needed to make R4 credible are mostly postponed to R5–R7. That ordering can produce a polished clock package whose uncertainty is calibrated to defects in its training calibrations and observation model.

I would not approve the roadmap as a sequence of gated phases. I would approve a narrower programme: establish a mathematically coherent lexical-time model; validate its implementation; curate a genuinely independent, protocol-characterised calibration panel; integrate supplied topology uncertainty; then run a frozen, family-level validation. Native topology MCMC, exhaustive topology evidences, full FBD, a second inference engine, and completed T2 should not be on the path to the first defensible Arawan result.

The proposed Arawan interval of roughly 1,000–4,000 BP around a 2,000 BP median is a useful illustration of what a **best-case conditional** answer might look like. It is not yet a credible product promise. The calculation assumes a between-family log-rate SD near 0.3 under the same model, list, coding protocol, and target depth; none of those conditions has yet been established for Arawan.

## Highest-priority changes

| Priority | Roadmap item(s) | Finding | Required change |
|---|---|---|---|
| 1 | R1.1, R1.4 | The proposed proper calendar-root prior and the claimed exact rate-prior convolution are not the same model and generally cannot both hold. | Split calendar-time-tree mode from lexical-time/ratio mode; specify the joint density and when re-dating is exact. |
| 2 | R4 versus R6 | Calibration reliability is downstream of clock fitting, although the clock is learned from those calibrations. | Move one-sided calibration semantics, provenance/dependence, conflict checks, and robust calibration likelihoods before R4. |
| 3 | R4 versus R5.1/R7 | Family-rate posteriors will condition on one selected topology and an only partly validated likelihood. | Move supplied tree-set integration, matched-likelihood checks, and misspecification coverage before the production R4 fit. |
| 4 | R4.5 | With about ten independent families, 95% coverage and PIT shape are nearly untestable; post-hoc inflation on the same folds is leakage. | Pre-register family-level folds and pass criteria; freeze the method; use proper scores and an untouched later test set or nested recalibration. |
| 5 | R4.2–R4.4 | `beta_f` is not transportable without a fixed rate normalisation and control for list composition, coding protocol, coder/source, ascertainment, and depth. | Define a canonical concept-rate scale and package domain; obtain crossed coder/protocol data or admit these effects are folded into predictive variance. |
| 6 | R2.5 | A per-language terminal offset is not identified separately from terminal branch length/rate. | Treat it as externally informed measurement-error sensitivity, preferably at source/list/coder level, not a freely learned default. |
| 7 | R5.2–R5.3 | Two new topology engines are expensive and do not address the dominant Arawan scale uncertainty. Laplace screening is especially fragile near polytomies. | Implement R5.1/R5.4 first; defer R5.2 and R5.3 until the fixed-tree target and evidences have independent oracles. |
| 8 | R8.1 | Arawan is postponed to the end even though it exposes the project’s decisive failure modes. | Make Arawan a continuous acceptance case from R0 onward, with calendar publication blocked until R4 passes. |

## 1. Ordering and the actual critical path

### What is ordered correctly

The opening diagnosis is sound. The current engine initializes unconstrained ages to `[0, 1e15]` (`src/gd_calibration.c:198-202`), gives an uncalibrated age no density beyond calibration factors (`src/gd_mcmc.c:196-259`), and by default mixes the tree likelihood with a positive CRP noise component (`src/gd_analysis.c:360-371`; `src/gd_model_novel_state.c:660-669`). Thus the present uncalibrated result is cap-defined, and without a proper age/lexical-scale prior the positive contamination floor makes the unbounded target improper. R1 must precede substantive dating.

The silent-input fixes in R0.1–R0.5, R0.8–R0.11, and R0.13 are also genuine blockers, not polish. Missing cognate IDs, free “modern” tips, wrong calibration support, zero-data successes, and internal-node/language confusion can change the scientific answer. R0.9 belongs early because no validation score should count an unconverged fold.

### What is misplaced

1. **R6 is a prerequisite of R4, not a later enhancement.** The clock is inferred from historical, archaeological, epigraphic, and tip-date evidence. If first attestations are encoded as symmetric dates, shared sources are multiplied as independent, or contested events are treated as exact, the error is absorbed into `beta_f` and shipped to Arawan. At minimum R6.1, R6.2, and R6.4—and a provenance/dependence model absent from the roadmap—must precede clock fitting. `literature-dating-methods.md:397-404,424-450` makes robust calibration and predictive validation part of the scale-transfer problem itself.

2. **R7.1, R7.2, and R7.5 are prerequisites of the R4 headline, not postscript validation.** SBC checks computation, richer-world coverage checks robustness, and agreement with an independent implementation checks the likelihood. A production clock learned under a model that has not passed these tests is not a reference clock. The roadmap’s R2 gate has some misspecification simulation, but it is not a substitute for matched fixed-parameter likelihood tests or depth- and protocol-transfer experiments.

3. **R5.1/R5.4 should precede the final R4 fit.** Training-family rates are affected by topology and ancient-language placement. Conditioning all training fits on one selected tree understates `beta_f` uncertainty and can shift it. Supplied posterior tree sets are much cheaper than native tree inference and are enough to propagate much of this uncertainty. The caveat in `literature-dating-methods.md:510-514`—reusing the same cognates once to estimate the tree and again to date it is an approximation—must be reported, but it is still preferable to pretending the topology is certain.

4. **R3 is too monolithic.** R3.1–R3.4 and R3.6 are needed before dozens of cross-family fits. R3.5 checkpoint binding is worthwhile engineering, but it should not gate the statistical feasibility study. Conversely, the global joint-prior interface and transformed parameterisation are partly R1 work, not merely later sampler optimisation: the current single-age update only adds a node-local calibration density (`src/gd_mcmc.c:499-520`), whereas R1.1 requires a joint time-tree density.

5. **Arawan must be an acceptance test throughout.** Waiting until R8.1 invites the tool to become excellent on large Eurasian training families while still mishandling sparse historical lists, disputed topology, and a genuinely uncalibrated small family. The Paumarí historical/modern comparison should be a negative control from the observation-model phase onward.

6. **The strict “no phase starts before the previous gate” rule is inefficient.** Calibration curation, Concepticon overlap analysis, coder/protocol metadata, fold design, and power calculations can and should start while R0–R1 are implemented. They are read-only/data-design work and may reveal that R4 is infeasible under the planned package contract before years are spent on inference machinery.

## 2. Scientific problems by item

### R1.1: “uniform-conditional” is proper, but not neutral and not the general MCMCTree construction

An explicit proper root prior times a correctly normalised conditional distribution for the remaining ages can remove the current volume tilt. That makes it a valid baseline. It does not make it scientifically uninformative:

- It induces a particular distribution on node ratios and branch spacings, which can be influential with five to seven languages and a relaxed clock.
- With ancient tips, multiple calibrations, and inherited order constraints, the normalising volume is piecewise and topology-dependent; it is not always just `root^-(m-1)`.
- If topology is sampled, a topology prior plus a ranking/age density is required. “Uniform on ages conditional on the root” does not define a neutral topology prior.
- Calling this “the MCMCTree construction” is inaccurate. The review’s own literature notes describe MCMCTree as a conditional **birth–death-with-sampling kernel** (`literature-dating-methods.md:41-50`), whereas uniform node ages conditional on the root are associated with the MrBayes uniform-clock prior (`literature-dating-methods.md:107-109`).

R1.1 should therefore offer `uniform_ratios` as a transparent baseline, not an unquestioned default. It should require prior-predictive plots and sensitivity against Yule/birth–death and alternative sampling schemes. Language datasets are usually diversified samples—one doculect per subgroup—not random complete samples; the literature explicitly says the sampling scheme changes the node-age prior (`literature-dating-methods.md:49-57,111-121`).

There is a second subtlety. Multiplying a tree prior by calibration densities does not generally preserve the user-specified root or node marginal. R1.1’s promise that the root marginal is “exactly the stated root prior” is true only for the specified conditional construction before other calibration factors, or if all those factors are included in the conditional normalisation. R1.5 is therefore mandatory even after R1.1: effective priors can still be distorted (`literature-dating-methods.md:70-80`; `literature-lexical-dating.md:480-493`).

### R1.4: the convolution identity is overgeneralised

The exact identity in `mathematical-analysis.md:501-516` is conditional on a scale-invariant age/ratio prior and no calendar calibration:

\[
\log a_v=\log(\rho a_v)-\log\rho.
\]

It does **not** follow from R1.1’s proper calendar-root prior. A calendar root prior couples `rho` and lexical age and supplies a second scale source. Nor is `p(a_root) proportional to 1/a_root` a proper rescue by itself; with the contamination floor it remains nonintegrable without additional proper structure.

The roadmap should define two explicit modes:

1. **Lexical-time mode:** put a proper prior directly on lexical root scale and on dimensionless ratios; fit lexical branch lengths; independently draw a clock-package rate to obtain calendar scenarios. Here post-hoc re-dating can be exact if all concept rates are dimensionless deviations around the common scale.
2. **Calendar time-tree mode:** use a proper calendar Yule/BD/FBD/root prior and calibrations. Refit when that prior or the calibrations change; do not advertise free convolution.

The identity breaks, or changes interpretation, in the cases the roadmap most needs:

- any node or root calibration and non-zero calendar tip age couples rate and time;
- ordinary calendar-scale Yule/BD/FBD priors couple shape/scale to calendar time;
- under a relaxed clock, even relative node chronology is separated from branch multipliers by their prior (`mathematical-analysis.md:495-499`), so `rho*a_v` and `a_v/a_root` are not simply “identified by data alone” as the roadmap’s opening claim says;
- concept heterogeneity does not break scale invariance only if `r_fc = rho_f * lambda_c` has an explicit normalisation and `lambda_c` is dimensionless and independent of `rho_f`;
- time-dependent/epoch rates do not generally admit one scalar convolution.

Thus the headline statement at `docs/roadmap.md:23-27` needs a strict-clock qualifier. Under a relaxed clock the directly data-supported objects are branch lexical lengths; node “lexical ages” and relative ages are model/prior-dependent decompositions.

### R2.1: discrete-gamma concept rates are reasonable, but the scale and tail are the model

Integrating a per-concept latent rate over 4–8 categories is a reasonable C99 design and preferable to sampling roughly 200 weakly identified scalar effects. But the roadmap omits three essentials:

- constrain the relative-rate distribution to a fixed mean (or fixed geometric mean), or `rho_f` and the concept-rate centre are nonidentifiable;
- define the constraint over a canonical concept universe, not separately over each family’s observed list, or list composition silently changes `beta_f`;
- compare gamma, lognormal, and at least one heavier/flexible slow-tail model. The slow tail controls deep dates, and the research shows materially different retention at depth (`mathematical-analysis.md:270-349`).

In a seven-language family the mixing distribution is weakly identified from only about six node depths. Turning a gamma/lognormal choice “on by default” does not make that uncertainty disappear. R2 should treat the mixing family as model uncertainty, scored by held-out family/depth prediction, and propagate between-model variance. It can also compete with contamination, borrowing, coding lumps, and time dependence, so recovery under the fitted model is insufficient.

### R2.2: Lomax survival is exact for a narrow frailty model, not a general branch correction

Replacing `exp(-r t)` by `(1+r t/k)^(-k)` is exact if a fresh, independent gamma multiplier is drawn for every branch × concept and integrated out (`mathematical-analysis.md:303-332`). It does not represent a persistent fast lineage, an autocorrelated clock, a persistent fast concept, or a borrowing process.

It also lacks the exponential semigroup property: subdividing a branch and independently redrawing frailty changes the survival distribution. That makes the model depend on where lineage boundaries are placed and deserves special scrutiny when comparing topologies or treating near-zero branches. R2.2 should therefore say “branch-boundary × concept frailty sensitivity,” not promise to “remove” deep-node bias. Compare it against lognormal/two-component mixing and test topology/evidence behaviour under branch subdivision before making it a default.

### R2.3–R2.4: “exact” synonym marginalisation and a generic loan mixture need generative semantics

R2.3 is exact computationally only after specifying what the latent variable means. If multiple forms are genuine coexisting synonyms, choosing one latent “primary” state is not observation error; it is a false single-state model. If the alternatives are uncertainty about which elicited form should be coded, their `primary_weight` values need a calibrated probabilistic interpretation before they can be used as prior probabilities. For the first release, alternative-coding sensitivity (strict/loose cognacy, one/all synonyms, documented selection rules) is safer than a Gibbs sampler that can turn arbitrary weights into apparently precise inference. The research itself calls a genuine polymorphism dynamics model expensive and not core (`mathematical-analysis.md:1071-1075`).

R2.4 similarly conflates distinct processes. An external loan is usually uninformative about change within the family and should be missing; an internal loan can create a perfectly tree-compatible shared class and bias a node young; uncertain cognacy is an observation model; and an alternative `analysis_id` is a model/coding sensitivity, not another independent observation. `loan_probability`, `inheritance_probability`, and `cognacy_confidence` cannot simply be treated as calibrated probabilities unless their provenance supports that interpretation. Specify these likelihoods separately, begin with exclude/sensitivity modes, and validate the compatible-loan case before advertising a general “mixture” (`literature-lexical-dating.md:1003-1008`; `mathematical-analysis.md:940-994`).

### R2.5: the proposed coding-error parameter is not identified by the cognates

Under T1, a homogeneous split error on one list is observationally equivalent to extra lexical evolution on that terminal branch. The mathematical note states this explicitly and says the tip offset needs a prior from double-coding studies (`mathematical-analysis.md:913-919`). A wider prior for a “low-quality” list does not identify error; it merely gives the model permission to explain away a long branch.

The unit is also wrong. The problem is usually **source/list/coder/orthography**, not an inherent per-language error rate. Paumarí 1833 and modern Paumarí should have separate observation-error effects even though they represent the same language lineage. Conversely, several languages coded by one lumper/splitter may share an error effect.

Make this a sensitivity term until there are independent duplicate codings, overlapping sources, or an informative split/lump prior. Prefer an explicit latent observation model that distinguishes splits from lumps where the patterns permit it. The known Paumarí pair is a validation datum, not evidence that a free offset successfully “found” coding error.

### R2.6: contamination and error models will compete

Sampling epsilon and exporting responsibilities is an improvement, but whole-concept CRP contamination is not a borrowing model. A compatible sister-to-sister loan remains tree-compatible and biases ages young; sending an entire concept to the CRP discards valid vertical information. The roadmap should also revisit the fixed CRP concentration/noise partition law noted in `reviews/2026-10-10-claude.md:293`, and it should state how epsilon, concept-rate tails, source error, and borrowing are kept distinguishable. Empirical concept-specific borrowability priors are suggested in `literature-lexical-dating.md:149-167,999-1008` and may be more defensible than a single exchangeable epsilon.

### R4.2–R4.4: `beta_f` is not automatically transportable

R4 correctly says the rate must be learned under the same model and concept list, but that condition is more restrictive than the roadmap acknowledges.

`beta_f` changes with:

- concept-list composition and Concepticon mapping;
- word- versus root-cognacy and lumper/splitter practice;
- treatment of loans, synonyms, partial cognacy, missingness, and singleton ascertainment;
- coder/source and automatic versus expert coding;
- topology/tree prior and relaxed-clock model;
- calibration depth when the slow concept-rate tail is wrong.

If each family has one dataset and one coding team, family rate, coder effect, and protocol effect are perfectly or nearly confounded. A random-effects regression cannot separate them merely by naming covariates. The training design needs crossed evidence—some families coded under multiple protocols and/or common coders spanning families—or the new-family predictive must include the combined family+protocol+coder dispersion. A manifest warning is not a statistical correction.

The factorisation also needs an anchor. If relative concept effects are centred within each family’s observed list, a list rich in slow meanings and a list rich in fast meanings define different `rho_f`. Define `lambda_c` on a canonical reference scale with a sum-to-zero/geometric-mean constraint and infer missing concepts by a covariate or hierarchical model. The literature specifically recommends frequency/semantic covariates for transfer to unseen concepts (`mathematical-analysis.md:867-878`; `literature-lexical-dating.md:991-998`).

Finally, do not “meta-analyse posteriors” naively. R4.2 produces a posterior for each `beta_f`; R4.3 then writes a population model for latent `beta_f` but omits the family estimation likelihood/uncertainty shown in `mathematical-analysis.md:520-531`. A two-stage analysis must carry each family likelihood (or a prior-corrected approximation with `s_f`), not treat posterior means/draws as raw observations and count their first-stage priors twice. A joint multi-family model is cleaner.

### R4.3: covariate selection is underpowered

“Speaker population first” is not justified as a near-term model. The cited evidence is contradictory and family-specific (`literature-lexical-dating.md:247-267`). With roughly ten independent families, even one slope plus residual variance is fragile, and trying population, contact, geography, subsistence, and writing tradition invites selection optimism. Pre-register at most one strongly motivated covariate with a shrinkage prior, evaluate it in fully nested LOFO, and expect the intercept-only predictive to win. Defer a serious meta-regression until the effective family count is much larger.

### R5.3: exhaustive topology enumeration with Laplace evidences is oversold

The topology counts are correct, but the feasibility claim assumes machinery that does not exist and benign posteriors that Arawan is unlikely to have. The current likelihood recomputes data extension/noise terms; there is no optimiser, Hessian, AGHQ, evidence oracle, or topology iterator. Each topology also requires a rebuilt tree/workspace. R3.4 is a hard prerequisite.

More importantly, the mathematical analysis already warns that Laplace evidence is topology-dependently biased when internal branches approach zero (`mathematical-analysis.md:767-773`). Screening on biased Laplace mass and refining only the winners can permanently discard precisely the near-polytomy topologies for which Laplace performed worst. Hierarchical concept mixtures, epsilon near a boundary, and relaxed branch effects further undermine Gaussian evidence approximations.

There is an Arawan mismatch too. The field case has nine doculects—seven modern plus two historical lists (`arawan-lessons.md:9-16`)—whereas R5.3 stops at eight. Nine rooted binary topologies number 2,027,025 (`mathematical-analysis.md:735-742`). The promised R8.1 topology-averaged Arawan result therefore cannot use R5.3 as written unless historical lists are removed or represented as observations on existing lineages. R5.1 is the practical route.

### R5.1: supplied tree sets are useful, but naive pooling is not a joint analysis

If the supplied posterior trees were estimated from independent evidence, pooling conditional age posteriors is straightforward. If they were estimated from the same cognate matrix—as the proposed BEAST/MrBayes route usually would be—combining external posterior topology frequencies with a second clepsydra likelihood fit is a hybrid of two models, not automatically a sample from either model's joint posterior. It can reuse the lexical evidence in incompatible ways and inherit the external program's tree prior, ascertainment, and coding assumptions. The roadmap should require one of: independent tree evidence; a genuinely joint clepsydra target; or a documented importance/reweighting approximation with simulation showing its bias. At minimum report the approximation and compare it with data-splitting by concepts. The source review already flags this same-data reuse (`literature-dating-methods.md:510-514`).

## 3. R4 validation: can it genuinely fail?

Not yet. The stated outputs are good—PIT, multiple coverage levels, log and interval scores, and baselines—but the gate and data design allow several false victories.

### Effective sample size and power

At ten independent families, nominal 95% coverage has binomial standard error about 0.069. Zero misses gives 100%; one miss gives 90%. Nominal 80% coverage has standard error about 0.126. A PIT histogram of ten family-level values has essentially no shape resolution. Scoring several nested nodes per family does not solve this because their clock error, topology, calibrations, cognates, and coders are correlated.

The 367-row calibration table is therefore misleading as a sample-size statement. Direct tabulation of `reference/calibrations/calibrations.tsv` gives 114 grade A, 84 B, 151 C, and 18 U rows. Among A/B rows, 120 are tips and only 67 are node/root targets. Those 67 occur under just 19 nominal family labels, and 40/67 are Indo-European, Austronesian, or Arawakan. Some remaining labels are nested descendants or creoles rather than independent draws from a worldwide family population. R4.1’s “at least ten” may be enough for a feasibility study; it is not enough to establish nominal 95% coverage or a population covariate.

### Leakage controls that must be explicit

For each outer LOFO fold, remove the held-out family from **all** of the following:

- family-rate hyperparameters and any interval inflation;
- concept-rate means, rate-tail/heterogeneity parameters, and contamination priors;
- covariate selection, preprocessing thresholds, protocol mappings, and package normalisation;
- calibration-support construction and any empirical root/time bounds;
- model selection and tuning if they were informed by that family’s score.

Within the held-out family, remove **all** node calibrations if the claim is uncalibrated-family prediction. Removing only the target node while retaining sister/root calibrations is LOCO, not LOFO. Decide in advance whether dated historical tips are retained. Retaining them tests heterochronous tip-dating, which is relevant to Arawan but is not the same as a fully contemporaneous uncalibrated family.

The repository already contains the leakage pattern to avoid: the existing LOCO helper computes a free upper support from all rows, including the subsequently held-out row (`reviews/2026-10-10-codex.md:274-280`). R4 needs a new fold-isolated pipeline, not a relabelled version of that benchmark.

### Calibration values are measurements, not truth

Most internal linguistic “dates” are uncertain historical or archaeological proxies. Grade B is not ground truth. PIT against the midpoint of a soft calibration asks whether one model reproduces another prior assumption. Conversely, using the same historical event both to train the family clock and to score a nested node is direct leakage.

Predeclare a small number of evaluation targets per family whose evidence is independent of the cognate data and of the calibrations used elsewhere in that fold. Carry their measurement uncertainty in the score. Secure tip dates are abundant but mainly test tip handling, not root-age transport; secure internal split dates are rare. Report those two validation strata separately.

Published ASJP or Bayesian lexical dates are not eligible truth. The Arawan ASJP value is explicitly labelled circular for lexical benchmarking in `calibrations.tsv:157`, and the literature warns that secondary phylogenetic calibrations can be too narrow (`literature-dating-methods.md:219-221`).

### The gate rewards width and post-hoc fitting

R4’s gate requires only coverage “within sampling error.” Arbitrarily wide intervals pass coverage, especially at n≈10. Although R4.5 lists proper scores, the gate does not require improvement over ASJP, Rama–Wichmann, or the node-count baseline. Require both calibration and utility: pre-specified interval/log-score improvement, reported sharpness, and maximum tolerable multiplicative width in the relevant depth range.

The fallback—fit inflation when coverage fails and report the inflated coverage—is not validation if the same folds estimate and assess inflation. It is training error. Preserve the failed result, estimate inflation inside a nested loop, and evaluate it on an outer fold, or freeze a new package and test it on later untouched families. With only ten families, honest nested validation will be noisy; that is evidence to report, not a reason to reuse folds.

### A design that can fail

1. Freeze the observation model, clock normalisation, concept mapping, calibration types, and one covariate specification before looking at outer-fold results.
2. Define independent family clusters; all nested subfamilies, shared data sources, and duplicated events stay in one fold.
3. Choose one primary internal target per family where possible, plus a separately reported tip-date panel. Carry target-date uncertainty.
4. Rebuild every learned artefact inside each fold.
5. Report family-clustered coverage with exact/binomial or bootstrap uncertainty, log score, interval score, median absolute log error, multiplicative width, and depth-specific bias.
6. Require improvement over a clock-prior-only predictor, ASJP where available, and a simple language-count/tree-length baseline. A likelihood that adds no predictive value should fail even if its intervals cover.
7. Publish the initial frozen result even if it fails. Any recalibrated successor gets a new version and new or nested validation.

## 4. What would concretely improve Arawan

### First fix the target dataset, not only the sampler

The Arawan field report contains nine doculects, including Aruá/Chandless 1869 and Paumarí/Natterer 1833, with only 41 and 70 concepts in the two old lists (`arawan-lessons.md:9-18`). Yet the calibration database contains the A-grade 1867 Arawá list and an **unverified 1889 Ehrenreich Paumarí/Jamamadi entry**, not the Natterer 1833 list (`reference/calibrations/calibrations.tsv:155-159`). Curating the actual 1833 source, its collection/publication date, lect identity, and uncertainty is a missing roadmap task.

The highest-payoff Arawan work is:

1. Independently re-transcribe/recode the two historical lists, ideally by at least two coders, and record source/coder uncertainty. This supplies the external information R2.5 otherwise lacks.
2. Produce alternative strict/loose cognacy and synonym resolutions rather than pretending one coding is true. The 208 manually collapsed polymorphic cells in the field report make this material (`arawan-lessons.md:38-41`).
3. Mark external loans as missing for internal inheritance and test plausible internal loans separately. Report concept and doculect influence, especially the slow concepts driving the root.
4. Decide whether historical lists are serial observations of existing lineages, possible ancestors, or separate tips. Do not let a full FBD implementation block a fixed-topology “observation on branch/ancestor sensitivity” for this case.
5. Use a supplied posterior tree set and named-set MRCA ages. Aruá’s position was fixed in the original analysis despite being described as unresolved (`arawan-lessons.md:15-16`).
6. Train/package a clock on data genuinely comparable to KAHD: similar Concepticon overlap, expert coding protocol, sparse/archival lists, and preferably South American low-documentation families. A Eurasian written-history clock may be out of domain even if its global average covers.
7. Make the Paumarí 1833–modern separation a release-blocking negative control. The observed 1,011–1,418-year estimate versus 165–170 years (`arawan-lessons.md:44-54`) must be predicted without giving a free source-error parameter the answer.

The Arawan outputs should be layered: branch lexical lengths; strict-clock relative chronology; relaxed-clock relative sensitivity; topology/coding/loan sensitivity; then several calendar scenarios from explicitly named clock packages. The weak C-grade 3,200-BP archaeological maximum may be shown as a sensitivity analysis, not folded silently into the “uncalibrated” headline. The C-grade ASJP root and U-grade Swadesh calculation are benchmarks only, exactly as `calibrations.tsv:157-159` says.

### Is 1,000–4,000 BP credible?

As arithmetic, yes under favourable assumptions. `mathematical-analysis.md:447-453` combines data/branch uncertainty with a cross-family log-rate SD of 0.3 to obtain root log-SD around 0.33, a 95% factor of about 1.9: 2,000 BP becomes about 1,050–3,800 BP.

As an empirical promise, no:

- `sigma_f=0.3` has not been estimated under T1 on matched expert-coded concept lists. Holman’s 29% discrepancy is in-sample and from 40-item ASJP distances (`literature-lexical-dating.md:517-528,921-930`), not an estimate of T1 family-rate variance. Rama–Wichmann’s 770-year LOOCV MAE is age-dependent and not a log-SD.
- With rate SD 0.4 plus branch/data uncertainty, the factor is roughly 2.3, already about 870–4,600 BP before topology, coding, borrowing, calibration, and model uncertainty. At 0.5 it is much wider.
- The ±6–8% “lexical data” figure assumes 200 meanings, known concept rates, correct model, and a seven-tip tree. KAHD is heavily incomplete, its historical lists have 41/70 meanings, and its coding/topology are uncertain. The roadmap states this idealised lower bound as though it described the actual data (`docs/roadmap.md:23-26`; compare `mathematical-analysis.md:434-474,1158-1162`).
- A 2,000-BP median is itself not data-guaranteed. The only published Arawan absolute estimate is the circular ASJP value around 1,764 years; Dixon-based arithmetic points older (`literature-lexical-dating.md:860-890`; `calibrations.tsv:157-158`).

The defensible promise is therefore conditional: “If a frozen, comparable clock package validates out of family with predictive log-SD near 0.3, then a factor-of-two 95% interval is plausible.” Until then, advertise a family of calendar scenarios, not 1,000–4,000 BP as a performance target.

## 5. Revised ordering

I would replace the phase order with the following.

### Phase A — Define the estimand and validation contract

- Rewrite the opening claim to distinguish strict-clock lexical ages from relaxed-clock branch lexical lengths.
- Specify lexical-time versus calendar-time-tree modes and the exact rate normalisation.
- Freeze the LOFO unit, leakage rules, target-date uncertainty, package compatibility rules, and success criteria.
- Audit the calibration panel, Concepticon overlap graph, coder/protocol metadata, and independent-family count now.

### Phase B — Minimal correctness and proper target

- R0.1–R0.5, R0.8–R0.11, R0.13, plus strict finite/range validation needed by the scientific path.
- R1.1–R1.6, but with separate lexical/calendar modes and a global joint tree-prior interface.
- R3.1–R3.3 enough to sample the new target reliably.
- Analytic/tiny-tree prior marginals, cap-invariance tests, exact prior-only runs, and full-pipeline SBC.

R0.12 and documentation cleanup should still be done, but need not block every statistical experiment.

### Phase C — Independent likelihood and data-model validation

- Move R7.5 here: match T1 fixed-parameter likelihoods with the Rönchen reference code and exact small-tree enumeration.
- Implement R2.1, R2.3, R2.4, and R2.6 plus the R2 PPC panel.
- Treat R2.2 and R2.5 as sensitivity models until they win predictive tests and have defensible priors.
- Run richer-world misspecification coverage (R7.2) now, including time depth, correlated concept domains, nonrandom missingness, coder splits/lumps, and compatible loans.

### Phase D — Efficient fixed-tree production engine

- R3.1–R3.4 and R3.6; checkpoint binding R3.5 can follow once the target format is stable.
- Add concept/bootstrap influence and power-scaling prior sensitivity from the literature notes.

### Phase E — Calibrations and supplied topology uncertainty

- Move R6.1, R6.2, R6.4, calibration-source dependence/provenance, and tip-age uncertainty here.
- Implement R5.1 and R5.4. Validate pooling and report the between-topology variance share.
- Use a limited ancient-observation/ancestor sensitivity for Arawan; full R6.3 FBD can remain later.

### Phase F — Cross-family clock and frozen validation

- Curate R4.1 with family/source/coder/protocol clusters, not row counts.
- Fit R4.2–R4.4 jointly or with prior-corrected family likelihoods and a canonical concept-rate scale.
- Run the frozen, nested R4.5 experiment. Gate on proper scores and sharpness as well as coverage.
- If it fails, report failure; design a new package version and validate it independently.

### Phase G — Arawan release result

- R8.1 immediately after a successful R4, with every uncertainty layer and negative control.
- Documentation/release items R8.2–R8.4.

### Phase H — optional expansion

- R5.2 native topology MCMC or R5.3 exhaustive evidences, chosen after comparative feasibility work rather than committing to both.
- Full R6.3 FBD, R7.3 composite likelihood, and a completed T2 only if a demonstrated user/scientific need justifies them.

## 6. Scope and effort: cut or defer

### Defer

- **R5.2 native topology MCMC.** Large dynamic-tree, calibration, caching, Hastings-ratio, and tree-diagnostic subsystem. R5.1 supplies most immediate value.
- **R5.3 exhaustive L=8 Laplace/AGHQ.** High implementation and validation cost; fragile evidence screening; it does not cover the nine-doculect Arawan case as represented in the field report.
- **R6.3 full FBD sampled ancestors.** Scientifically important for ancient-language-rich families, but a narrower ancestor/serial-observation sensitivity can address Arawan sooner.
- **R7.3 pairwise composite/Godambe engine.** Attractive independent check, but not before the primary target is validated. It does not solve the cross-family scale problem.
- **Full R7.4 T2 completion.** Mark T2 experimental now and stop making product claims. Completing it is a separate project with ascertainment and meaning-weighting hazards.
- **R3.5 polished resumability.** Useful for production, low payoff before model freeze.

### Narrow

- **R4.3 meta-regression:** intercept-only first; one pre-registered shrinkage covariate at most.
- **R2.3 synonym latent-state MCMC:** begin with alternative-coding sensitivity and exact pairwise handling. A full within-T1 synonym dynamics model is expensive and weakly identified.
- **R2.5 source error:** sensitivity with external priors, not free per-language estimation.
- **R6.1 calibration outliers:** begin with heavy-tailed/mixture densities and shared-source clusters; do not build a large OxCal clone.
- **R7.1 “every executable configuration”:** define supported publication profiles and parameter-block coverage. The Cartesian product of all options is not a finite practical gate.

### Keep high priority

- R1 proper target/prior-only/effective-prior work;
- R2.1 concept heterogeneity with explicit scale anchoring and mixing-family sensitivity;
- R2 PPCs and realistic misspecification tests;
- R3 ratio geometry, clock-preserving moves, caching, complete draws/pointwise output;
- calibration provenance, one-sided likelihoods, and conflict diagnostics;
- R5.1/R5.4;
- the R4 fold-isolated package and a continuously exercised Arawan benchmark.

## 7. Literature lessons the roadmap omits or weakens

1. **Tree sampling design.** Birth–death priors depend on complete/random versus diversified sampling. Lexical datasets deliberately sample one doculect per subgroup (`literature-dating-methods.md:49-57,111-121`). A sampling fraction alone is not enough.

2. **Rate-scale prior collapse over many concepts.** Independent per-concept rate priors can make the average rate artificially precise as concept count grows. The global scale must retain cross-family uncertainty and concept effects need a centring/gamma-Dirichlet construction (`literature-dating-methods.md:200-217,435-444`). R2.1/R4 should state this explicitly.

3. **Calibration dependence and uncertain target assignment.** Several dates for one event should not simply be multiplied; calibration sources can share evidence, and archaeology may attach to one of several nodes. Hierarchical multi-source likelihoods and target-node mixtures are suggested in `literature-dating-methods.md:239-251,340-351`. R6 models only independent outliers.

4. **Clock adequacy, not only cognate-pattern PPC.** The literature recommends comparing unconstrained lexical branch lengths with clock-predictive branch lengths (`literature-dating-methods.md:381-386,516-520`). R2’s PPC panel does not directly test whether the clock can reproduce branchwise rate structure.

5. **Prior sensitivity and information limits.** Power-scaling by component and amplified-data/infinite-sites reports are cheap and directly explain what more cognates cannot improve (`mathematical-analysis.md:821-841`; `literature-dating-methods.md:470-483`). R1.5 overlap alone is not a substitute.

6. **Concept influence and dependence.** Meanings are not necessarily independent because of semantic domains, derivation, concerted sound change, and coding practice. A Bayesian bootstrap/concept influence report is explicitly recommended (`literature-dating-methods.md:298-300,562-565`), while the code review warns that model uncertainty must include concept-domain dependence (`reviews/2026-10-10-codex.md:282-288`).

7. **Chance-resemblance/saturation floor.** Deep nodes should be flagged when shared-state signal is indistinguishable from chance (`literature-lexical-dating.md:190-209,1069-1074`). CRP contamination responsibilities are not the same diagnostic.

8. **Strict-versus-relaxed model uncertainty.** The literature proposes a strict-clock spike/model averaging and notes clock choice is particularly important for tip dating (`literature-lexical-dating.md:473-478,1059-1063`). Merely offering `--relaxed` leaves users to select the answer-producing model.

9. **Ancient-list time signal should be tested before dating.** Root-to-tip temporal-signal checks and date randomisation are cheap analogues of TempEst (`literature-dating-methods.md:139-145`). For Arawan’s 165–170-year shallow span and noisy old list, failure to detect temporal signal would be important evidence.

10. **Polytomies and ILS.** Near-zero branches, rapid radiations, and persistence of polymorphism are not just topology nuisances (`literature-lexical-dating.md:587-599,1075-1081`). A binary-topology enumerator with a continuous prior that assigns zero mass to true polytomies can be confidently wrong.

## Bottom line

The roadmap should not promise that a sufficiently elaborate C99 model will manufacture calibrated calendar dates for Arawan. The achievable result depends overwhelmingly on whether a model- and protocol-matched cross-family predictive distribution transports to a sparse Amazonian dataset. That is an empirical question, and the current R4 design lacks the independence, power, and leakage controls to answer it.

The shortest path to a defensible result is: proper lexical/calendar target definitions; exact prior diagnostics and independent likelihood checks; realistic observation-model stress tests; robust calibration curation; supplied-tree integration; then a frozen family-level validation scored for both calibration and usefulness. If that experiment fails, the correct Arawan product is still valuable—lexical branch lengths, relative chronology, and explicit calendar scenarios—but it is not a single well-calibrated BP posterior.

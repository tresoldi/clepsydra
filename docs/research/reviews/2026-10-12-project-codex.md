# Independent review of `clepsydra`

Review date: 2026-10-12. Repository state: `0a671f3`.

I reviewed the requested documentation, empirical result pages, earlier 2026-10-10 reviews, recent history, examples, principal C implementation, public headers, and tests. I ran `build/clepsydra --help` and all compiled native unit, likelihood, regression, and simulation executables. They passed. LeakSanitizer had to be disabled because it cannot run under this sandbox’s tracing; Python tests could not start because the globally read-only environment provides no usable temporary directory. I made no changes.

“Verified” below means checked directly in code, data, or an executed binary. “Inferred” means a scientific interpretation of those materials.

## 1. What the project is and what it has done

Clepsydra is a dependency-light C99 engine for estimating divergence times from lexical cognacy data. Its serious model is T1: a fixed or externally sampled tree whose cognate states evolve by a novel-state replacement process. The current default adds:

- an uncorrelated lognormal relaxed branch clock;
- integrated lognormal concept-rate heterogeneity;
- a sampled whole-concept contamination mixture;
- optional list-specific observation error;
- dated tips and several node-calibration distributions;
- a proper prior on lexical root depth and conditional node ages.

It also supplies deterministic pairwise dating, an implemented but non-runnable binary gain/loss likelihood (T2), dating over a posterior tree sample, calibration sensitivity, posterior predictive checks, PSIS-LOO model averaging, and prior-versus-posterior reporting.

The central inferential distinction is correct and valuable: cognates directly identify lexical branch length, approximately `rate × time`; calendar time additionally requires a clock or calendar calibration. The output accordingly separates lexical, relative, and calendar ages (`docs/scientific-methods.md:18-33`). In an uncalibrated analysis, the calendar posterior is substantially a transformation of the chosen clock prior rather than information created by cognates (`docs/scientific-methods.md:35-52`).

Inference is Metropolis-within-Gibbs with adaptive warm-up, four chains, robust rank-normalized/folded R-hat, bulk and tail ESS, checkpointing, parallel chains, and several ridge-following proposals for the difficult relaxed-clock geometry. Concept-rate integration uses Gauss-Hermite or generalized Gauss-Laguerre quadrature (`src/gd_model_novel_state.c:924-995`). The default is four chains with 4,000 warm-up and 4,000 retained iterations, relaxed clock, sampled contamination, concept-rate mixture, and broad log-clock and lexical-root priors (`src/gd_analysis.c:527-559`). **Verified.**

The validation portfolio is now substantial:

- two structurally independent brute-force implementations for the T1 likelihood;
- brute-force and normalization checks for T2;
- several strict-clock SBC scenarios, including a free root, calibrations, ancient tips, Yule ages, and integrated concept rates;
- exact prior-moment checks for relaxed-clock transformations;
- recovery/confounding tests for branch rates;
- posterior predictive checks;
- deliberately misspecified simulated worlds;
- a real-data leave-one-calibration-out benchmark.

The native test binaries all passed in this review. In particular, T1 reported 34/34 tests, T2 17/17, tree MCMC 17/17, relaxed MCMC 6/6, and all regression/unit suites passed. **Verified by rerun.**

The main empirical findings are more cautious than the project’s headline aspiration:

- The cross-family clock pilot is a pre-registered no-go for sharp uncalibrated calendar dates. Estimated total log-rate spread is about 0.77; a nominal 2,000-year node has an approximate 95% range of 330–12,000 years (`docs/research/clock-pilot.md:38-51,76-98`). The learned clock does not demonstrably beat simple baselines.
- Arawan lexical chronology is stable across the investigated model variants: root lexical age is about 0.64 expected replacements per concept (`docs/research/arawan-result.md:18-40`). Its calendar root is not stable: approximately 1,950–4,250 BP in median depending on what supplies the scale, with very wide intervals (`docs/research/arawan-result.md:42-70`).
- Arawan model variants differ little in held-out concept prediction, so model choice contributes little to the already large root uncertainty. Topology affects internal clades more than the family root.
- Classical formulas give Arawan root estimates from roughly 1,750 to 6,600 BP and turn the noisy 1833 Paumari list into 1,000–2,500 years of false divergence (`docs/research/arawan-glottochronology.md:40-59,107-140`).
- The real LOCO benchmark is unfavorable to the strict clock: 1/9 coverage and interval-score skill −1.00. The relaxed clock covers 7/9, but its aggregate skill remains −0.13 against prior-only prediction (`examples/loco_benchmark/README.md:64-85`).
- Misspecification simulations identify serious failure regimes: coding error gives 32% coverage with positive age bias, while a rate decreasing with depth gives 31% node coverage and 5% root coverage (`docs/research/simworld.md:31-60`). A general list-error parameter repairs coding-error worlds but damages clean ones because terminal rate and error are confounded (`docs/research/simworld.md:62-90`).

### Recent trajectory

The last 80 commits show a rapid transition from a promising but scientifically unsafe 0.2 engine toward a considerably more defensible 0.3:

1. The October 10 reviews exposed an improper/hidden age prior, erroneous treatment of missing cognates, weak calibration semantics, incomplete diagnostics, unmodelled heterogeneity and observation error, and no topology propagation.
2. Subsequent commits introduced a proper lexical-root/node-age prior, strict input validation, corrected calibration handling, missing-data semantics, modern diagnostics, concept-rate integration, error/contamination models, tree-sample dating, checkpoints bound to inputs, pointwise likelihoods, and parallel chains.
3. The most recent work has concentrated on empirical calibration, negative results, prior-versus-posterior decomposition, relaxed-clock mixing, and honest presentation rather than adding another likelihood.

That is the right trajectory. The project has responded unusually well to adversarial review. It has not, however, accumulated enough independent empirical validation to turn engineering maturity into reference-method status.

## 2. Does it achieve its goal?

### Short answer

It now achieves a narrower and scientifically useful goal:

> A transparent, reproducible engine for conditional lexical and relative chronologies, and for broad calendar scenarios when calibration information is weak.

It does **not yet** achieve the stronger stated goal of producing generally well-calibrated posterior calendar ages from cognate lists or serving as a field reference. Its own strongest external tests say those calendar distributions are often prior-dominated, uninformative, or miscalibrated under plausible misspecification.

### Where the evidence is strong

1. **The T1 likelihood is unusually well checked.** Two independent enumerations, normalization, random cases, polytomies, ancient/internal observations, boundary rates, and simulation recovery provide strong evidence that the recursion computes the intended likelihood (`docs/validation.md:79-109`). **Verified in source and rerun tests.**

2. **The former improper-age-prior failure has been repaired.** The code now constructs a proper lexical-root prior and conditional chronology volume (`src/gd_tree_prior.c:65-104`), and the simulation suite checks prior moments and SBC with a free root. The October 10 critical finding about the `1e15` ceiling is no longer current. **Verified.**

3. **Sampler transformations have targeted validation.** The scale moves include explicit Jacobians—for example, subtree scaling uses `k × delta` (`src/gd_mcmc.c:1198-1243`)—and the relaxed-clock prior-only tests reproduce exact moments. **Verified.**

4. **Diagnostics are genuine rather than cosmetic.** The implementation computes the maximum of bulk and folded rank R-hat, plus bulk and tail ESS, and rejects runs failing R-hat `<1.01` or ESS `≥100` per chain (`src/gd_diagnostics.c:323-400`; `src/gd_analysis.c:1140-1160`). **Verified.**

5. **The software reports prior dependence unusually clearly.** Prior-only node summaries, rate contraction, scale-source flags, and the root-width covariance split are exactly the information most linguistic dating papers omit. **Verified in implementation and examples.**

6. **Negative empirical results are preserved.** The clock pilot, LOCO scores and misspecified-world simulations are not buried. That substantially increases confidence in the project’s scientific culture.

### Where the evidence is weak

1. **External calibration is tiny.** The decisive end-to-end benchmark has nine tips from three families, with one structurally unrecoverable Romance case and a Sino-Tibetan fit still failing the stated convergence threshold at 10,000 + 10,000 iterations: R-hat 1.024, ESS 212 (`examples/loco_benchmark/README.md:68-85`). Seven-of-nine coverage is therefore not an estimate of calibration in any conventional sense. **Verified.**

2. **Clean-world coverage is already below nominal.** The misspecification study obtains 0.81 rather than 0.90 node coverage under its “clean” scenario (`docs/research/simworld.md:31-48`). The document reasonably notes that the topology/age generator differs from the fitted prior, but that means this is a model-robustness experiment, not confirmation of operational calibration.

3. **The strongest SBC remains small-tree, in-model validation.** Passing SBC says the implementation samples its declared posterior. It does not establish that the posterior describes lexical evolution. Current relaxed-clock tests emphasize moments and recovery in identifiable synthetic configurations rather than full parameter-wise SBC over realistic trees, missingness, ancient tips, coding errors, loans and calibrations. **Inferred from `docs/validation.md:138-210` and test coverage.**

4. **The clock prior is only provisionally empirical.** The pilot uses 15 two-sided nodes from six families and a simplified pairwise estimator, finds wide uncertainty, and explicitly fails its pre-registered usefulness criterion (`docs/research/clock-pilot.md:19-37,79-98`). That is good evidence against sharp universal clocks, not a mature reference-prior package.

5. **The observation model is still coarse.** Thresholding `loan_probability`/`inheritance_probability` excludes rows rather than modelling their uncertainty (`src/gd_forms.c:189-217`). Cognacy confidence is parsed, but the scientific meaning of supplied confidence values remains dependent on how they were produced. Whole-concept contamination handles incompatible patterns; it is not a borrowing process and cannot identify compatible internal loans. **Verified in code; scientific consequence inferred.**

6. **The relaxed clock is only weakly identified without multidepth calibration.** The program diagnoses branch-duration confounding, but diagnosis does not restore information. Branch-specific multipliers can absorb age, lineage, coding, borrowing, and model inadequacy.

### Claims that outrun the evidence

- README’s “full validation suite” (`README.md:58-65`) is fair as a description of breadth but can be read as “fully validated.” The suite validates implementation much better than empirical calibration.
- “The relaxed clock is roughly honest” (`examples/loco_benchmark/README.md:94-98`) is too strong. Seven of nine coverage, aggregate negative interval-score skill, three families, correlated design choices, and one unconverged family establish only that it is less overconfident than the strict clock in this small benchmark.
- The claim that relaxed prediction is “never worse” supports a default only weakly. Four families and concept-wise PSIS-LOO do not establish transportability; more importantly, the real LOCO aggregate score is still worse than prior-only.
- “Lexical ages are trustworthy when coding is careful” (`docs/research/simworld.md:92-100`) should be narrowed to “showed tolerable behavior in the simulated scenarios tested.” Even clean coverage was 0.81, and “careful coding” is not an empirically measurable assumption for most historical datasets.
- Arawan’s “headline calendar estimate” of 2,235 BP (`docs/research/arawan-result.md:61-67`) is questionable wording. The same document establishes that the calendar scale is overwhelmingly clock-driven and that the clock pilot was a no-go. The scientifically primary result is the lexical chronology; the calendar distribution is a named external-clock scenario.

## 3. Correctness and validation risks

### Critical scientific risks

1. **Time-dependent rates can produce severely miscalibrated deep ages with no internal warning.**  
   The simulated declining-rate world covers the root once in 20 and biases ages young (`docs/research/simworld.md:53-60`). The default uncorrelated branch clock is not a time-dependent clock: it permits branch scatter but does not encode systematic rate decay with depth. No modern-only dataset can identify this cleanly. **Verified result; inferred general severity.**

2. **Coding/list error is confounded with shallow time.**  
   Eight percent list error produces zero coverage for shallow nodes in the reported simulation. Adding a generic list-error term shifts real terminal evolution into error and reduces clean coverage to 0.56 (`docs/research/simworld.md:62-81`). This directly threatens the node ages users are most likely to regard as precise.

3. **Calendar calibration is not demonstrated across enough independent families.**  
   The strict clock clearly fails, while the relaxed clock mainly widens intervals and has negative aggregate skill against the prior (`examples/loco_benchmark/README.md:68-105`). The goal specifically promises well-calibrated posteriors; this is the central unmet requirement.

### High implementation/statistical risks

4. **`date-trees` pools draws from unconverged per-tree analyses.**  
   The loop appends clade draws regardless of `res->converged` (`cli/commands.c:1195-1210`) and only returns status 3 after writing the pooled posterior (`cli/commands.c:1220-1240`). The Arawan topology result includes 18 converged and two unconverged trees (`docs/research/arawan-result.md:73-79`). A user can therefore read an apparently ordinary aggregate containing unreliable components. **Verified.**  
   At minimum, pooled summaries should default to converged trees only or be withheld entirely, with an explicit sensitivity comparison if exclusions alter topology weights.

5. **The model combines several weakly identifiable flexibility terms.**  
   Branch variation, concept-rate spread, contamination, tip error, node ages and clock rate can explain overlapping features. The project has good one-feature recovery tests, but not a full factorial SBC/coverage design for their joint default configuration on realistic trees. **Inferred.**

6. **Topology propagation is post hoc, not joint Bayesian inference.**  
   `date-trees` runs the lexical model over externally sampled trees and pools equally selected tree draws. Unless the supplied topology sample is already conditional on exactly the same cognate data/model, this is not a joint posterior; if it is based on the same cognates under another model, it risks an incoherent modular update or double use of evidence. The output is useful sensitivity analysis but should not automatically be called full topology uncertainty.

7. **Calibration dependence is not modelled.**  
   Multiple historical/archaeological calibrations may share sources, assumptions or reconstruction chains. Their likelihood factors are multiplied as independent evidence. Provenance is recorded, but statistical dependence is not represented. This can make apparently rich calibration sets overconfident.

8. **Whole-concept contamination misses the most dangerous loans.**  
   Tree-incompatible patterns can be absorbed. A compatible transfer within a clade resembles inheritance and can make that clade too young. The existing borrowing simulation is only one design and its acceptable coverage should not be generalized.

### Medium engineering/numerical risks

9. **Quadrature accuracy varies by parameter regime.**  
   The rerun likelihood test reported gamma-rule relative errors of 6.2% at eight nodes, 1.6% at 16 and 0.39% at 24, while the adaptive rule uses more nodes for high CV (`src/gd_model_novel_state.c:973-982`). Defaults appear safe, but user-set extreme spread should trigger a numerical accuracy diagnostic or adaptive convergence check.

10. **Malformed primary weights can silently become 1.0.**  
    Nonfinite or negative numeric weights are rejected, but a non-numeric value increments `n_unreadable_weights` and retains the default weight 1 (`src/gd_forms.c:220-244`). That is inconsistent with the otherwise strict scientific-input policy. **Verified.**

11. **Fixed default seed encourages accidental pseudoreplication.**  
    Seed 42 is good for reproducibility (`src/gd_analysis.c:537-542`), but users may mistake repeated default runs for independent robustness checks. The manifest makes this visible, so severity is moderate.

12. **T2 has extensive likelihood validation but no usable inference path.**  
    Keeping it clearly labelled experimental is correct. Its presence nevertheless enlarges maintenance and documentation surface without yet contributing to empirical validation.

### Tests that would not catch important errors

- In-model SBC cannot reveal time-dependent clocks, compatible loans, systematic coding conventions, nonrandom semantic change, or circular calibrations.
- Small quartets do not stress large-tree chronology-volume numerics, deep unbalanced trees, extensive polytomies, or severe calibration conflicts.
- Aggregate chi-square rank tests have limited power and can miss local tail defects; ECDF-based SBC plots and parameter-specific coverage would be stronger.
- The LOCO benchmark’s tip dates test heterochronous prediction, not divergence-node calibration.
- Posterior predictive checks based on class counts/incompatibility can miss a compatible wrong history.
- Passing convergence thresholds does not diagnose nonidentifiability or prior domination.

## 4. Low-hanging fruit

### Engineering

1. **Do not pool unconverged `date-trees` runs** — 2–4 hours, high payoff.  
   Refuse the aggregate by default, or emit separate all-tree and converged-only summaries with retained topology mass.

2. **Reject unreadable weights instead of substituting 1.0** — 1–2 hours, medium payoff.  
   Preserve an explicit opt-in lenient mode only if legacy imports require it.

3. **Add quadrature convergence diagnostics** — 1–2 days, medium/high payoff.  
   Compare `K` and `2K` for extreme concept-rate spreads, or reject parameter ranges where the bound exceeds a declared tolerance.

4. **Expose per-block acceptance and effective jump diagnostics** — 1–2 days, medium payoff.  
   Current R-hat/ESS catch outcomes but not why difficult fits fail.

### Statistics

5. **Make the principal estimand lexical age in unanchored runs** — hours, high payoff.  
   Calendar columns can remain, but should be named and visually grouped as clock-prior scenarios.

6. **Add a time-dependent-clock sensitivity model** — several days to two weeks, very high payoff.  
   Even a simple epoch/depth trend is more important than another contamination variant because the present failure is large and undetectable.

7. **Run joint-default SBC and simulation coverage** — 2–5 days of implementation plus compute, high payoff.  
   Exercise relaxed clock, integrated concept rates, sampled epsilon, dated tips, calibration outliers and observation error together.

8. **Report calibration windows as targets, not only midpoint points** — about 1 day, medium payoff.  
   Score predictive distributions against the historical-date distribution/window.

### Data

9. **Expand LOCO before making calibration claims** — days to weeks of curation, highest empirical payoff.  
   Aim for at least 30–50 withheld observations across 8–10 families and several depths, with family-clustered uncertainty.

10. **Create replicated-coding benchmarks** — several days, high payoff.  
    Use independently coded cognate analyses to measure real rather than simulated coding error and to calibrate `cognacy_confidence`.

11. **Curate compatible borrowing cases** — several days, high payoff.  
    The current adequacy machinery is best at visibly incompatible patterns; it needs examples where contact preserves tree compatibility.

### Documentation

12. **Replace “roughly honest” with the exact LOCO result** — under 1 hour, high payoff.  
    Say “less overconfident in this nine-target benchmark; aggregate interval-score skill remained −0.13.”

13. **Demote the Arawan calendar headline** — under 1 hour, high payoff.  
    Headline the 0.64 lexical root and present 2,235 BP explicitly as the cross-family-clock scenario.

14. **Separate completed validation from planned validation in `docs/validation.md`** — 2–4 hours, medium payoff.  
    Sections containing prospective lists can otherwise look like delivered coverage.

15. **Describe `date-trees` as modular sensitivity analysis unless its input semantics justify a joint posterior** — 2 hours, high payoff.

## 5. Comparison with other methods and tools

| Method/tool | What it offers | Clepsydra advantage | Clepsydra deficit / likely reviewer objection |
|---|---|---|---|
| Swadesh / Lees glottochronology | Simple retention formula with a fixed universal rate | Full posterior, tree-consistent chronology, dated tips, rate uncertainty, concept/branch heterogeneity, observation-error terms | Clepsydra’s uncalibrated calendar scale still ultimately comes from a transferable clock prior; its gain is mainly honest uncertainty, not clock-free dating |
| Starostin | Attempts to correct changing replacement behavior and emphasizes loan removal/stable lists | Joint tree chronology and uncertainty propagation; explicit empirical failure checks | Clepsydra lacks Starostin’s substantive lexical curation protocol and does not yet model systematic time dependence, the failure that its simulation shows is most damaging |
| ASJP dating, Holman et al. 2011 | Uniform 40-item data, automated lexical distance, large cross-family calibration and regression | Cognate-based generative model, full node posteriors, explicit priors/calibrations, better auditability | ASJP’s breadth is far larger; Clepsydra’s empirical clock study covers only six families and does not demonstrate improved out-of-family prediction |
| BEAST 2 / BEASTling | Joint topology, branch lengths, clock and calibrations; mature MCMC ecosystem; language-oriented setup in BEASTling | Small portable binary, fast specialized likelihood, no Java stack, direct prior/posterior decomposition, browser potential | No joint topology inference, fewer clock/tree priors, no mature operator ecosystem, marginal-likelihood machinery, XML-level extensibility, or decades of cross-checking |
| RevBayes | Programmable graphical models, rich tree/clock priors and custom moves | Much simpler reproducible interface for this one task; easier deployment and auditing | Far less model flexibility; novel models cannot be assembled compositionally; fewer inference algorithms |
| MrBayes | Mature Bayesian phylogenetics, topology inference, established clock/tree-prior machinery | Lexical-specific observation and reporting layer; explicit list error and cognacy-oriented outputs | MrBayes is far more established for tree-space inference; Clepsydra’s conditional-topology strategy is not equivalent |
| TraitLab / stochastic Dollo | Models cognate gain/loss, borrowing-related extensions, stochastic Dollo process tailored to lexical traits | T1 works on multistate meaning partitions rather than treating every cognate class as an independent binary character; avoids some within-meaning pseudoreplication | T1 assumes novel replacement states and has no explicit lateral-transfer history; TraitLab-type models address processes Clepsydra currently relegates to contamination/sensitivity |
| Gray & Atkinson; Bouckaert et al.; Chang et al. 2015 | Full Bayesian language phylogenies, relaxed clocks, archaeological/historical calibrations; influential empirical demonstrations | Stronger audit trail around prior domination, calibration prediction and negative results; lexical-age output is more honest than a single historical date | Those analyses jointly estimate topology and often use richer clock/tree models. Clepsydra has not yet shown that its specialized likelihood changes empirical conclusions or improves predictive calibration |
| Rönchen–Wiklund–Hammarström 2024 and recent lexical-dating work | Renewed focus on pairwise/tree-aware lexical dating, heterogeneous stability, empirical evaluation | Clepsydra operationalizes many recommendations in a reproducible C engine and separates validity from adequacy | A reviewer will ask for an external head-to-head benchmark on identical data, not merely architectural comparison or in-repository simulations |
| Recent large lexical phylogenies such as Heggarty et al. | Broad datasets, sophisticated historical calibration, topology and date inference | Clepsydra is easier to inspect and explicitly displays how much age width comes from data versus clock | Data breadth and historical calibration depth are not comparable; Clepsydra’s empirical evidence remains pilot-scale |

The distinctive contribution is not presently “more accurate dates.” It is the combination of:

- a multistate lexical likelihood;
- an explicit lexical/calendar separation;
- prior-only comparison and width attribution;
- robust, inspectable C99 implementation;
- candid predictive failures.

That is already publishable. Claims of superior calendar accuracy are not.

## 6. Paper angles

### A. Software paper

**Claim:** Clepsydra is a portable, reproducible C99 engine for Bayesian lexical chronology with unusually explicit prior-dependence and convergence reporting.

**Existing evidence:** Public API/CLI, native and WebAssembly builds, examples, extensive tests, checkpointing, parallel deterministic chains, likelihood brute forces, diagnostics and model cards.

**Missing:** A stable 0.3 release, archived artifact/DOI, fuller user tutorial, API stability statement, independent installation/reproduction on multiple systems, and performance comparison with BEAST/RevBayes workflows.

**Venues:** *Journal of Open Source Software*, *SoftwareX*, possibly *Methods in Ecology and Evolution* if paired with stronger methodological evidence.

**Strongest objection:** It is feature-rich but not yet used independently, and its flagship T2 likelihood is not end-to-end runnable.

### B. Methods paper: lexical time before calendar time

**Claim:** Separating lexical, relative and calendar ages, with explicit prior-width attribution, prevents false precision in lexical dating.

**Existing evidence:** Proper prior construction, analytic prior checks, Arawan decompositions, cross-family clock no-go, prior-only overlays, width split, sensitivity outputs.

**Missing:** Formal derivation in a paper-quality supplement, broader examples, comparison with corresponding BEAST/BEASTling analyses, and a decision-theoretic evaluation showing that the reporting changes conclusions.

**Venues:** *Systematic Biology*, *Journal of Language Evolution*, *PLOS ONE*.

**Strongest objection:** The separation is conceptually correct but may be seen as reporting discipline rather than a new inferential method.

### C. T1 likelihood/methods paper

**Claim:** A novel-state multistate likelihood models cognate partitions directly and avoids treating within-meaning cognate classes as independent binary characters.

**Existing evidence:** Mathematical specification, two independent brute forces, normalization, simulation recovery, missing/internal observations, polytomy support, concept-rate integration.

**Missing:** Direct comparison with stochastic Dollo/binary models under realistic data-generating processes; power, bias and runtime studies; explanation of the relationship to infinite-alleles/no-homoplasy models; empirical model comparison in a shared observable space.

**Venues:** *Systematic Biology*, *Bioinformatics*, *Journal of Language Evolution*.

**Strongest objection:** The no-recurrence/novel-state assumption and lack of explicit borrowing may be less realistic than the dependence problem it solves.

### D. Empirical Arawan paper

**Claim:** Arawan cognates robustly support a relative lexical chronology, but available historical lists and cross-family clocks do not sharply identify calendar ages.

**Existing evidence:** Nine doculects, stable lexical root around 0.64, historical-list error analysis, tree-sample sensitivity, model variants, adequacy checks, classical comparisons.

**Missing:** Linguist-led cognate and loan audit, multiple defensible topologies/codings, source-critical validation of historical lists, influence analysis by concept/doculect, and external expert treatment of archaeological constraints.

**Venues:** *Diachronica*, *Journal of Historical Linguistics*, *Language Dynamics and Change*.

**Strongest objection:** The calendar conclusion is mostly a negative identifiability result, while the topology and cognate coding are inherited rather than independently established.

### E. Validation/benchmark paper

**Claim:** Standard strict-clock lexical dating is badly overconfident on held-out historical tips, while relaxed clocks improve coverage mainly by widening uncertainty and often add little beyond the prior.

**Existing evidence:** LOCO strict 1/9 versus relaxed 7/9, interval scores against prior-only, synthetic coverage, model-uncertainty runs, convergence diagnostics.

**Missing:** An order-of-magnitude more families/targets, target-window scoring, cluster-aware intervals, preregistered analysis, comparison with BEASTling/RevBayes/TraitLab and ASJP predictions, and fully converged fits.

**Venues:** *Journal of Language Evolution*, *Systematic Biology*, *Royal Society Open Science*.

**Strongest objection:** Nine targets are too few, not independent, and one relaxed analysis fails the project’s own convergence criterion.

### F. Negative-results paper on universal lexical clocks

**Claim:** Neither classical retention constants nor a learned cross-family clock support useful sharp dates for an uncalibrated family; apparent precision comes from conditioning on a rate as known.

**Existing evidence:** Pre-registered clock pilot, factor-of-four disagreement among classical Arawan estimates, failure on Paumari 1833, broad cross-family predictive distributions, LOFO comparisons.

**Missing:** Uniform coding across many families, stronger historical calibration curation, uncertainty in node/calibration assignments, direct ASJP reanalysis, and robust treatment of phylogenetic/nonindependent calibration nodes.

**Venues:** *Journal of Language Evolution*, *Language Dynamics and Change*, *PLOS ONE*.

**Strongest objection:** The clock pilot’s simplified pairwise estimator may inflate heterogeneity that the full tree model could explain.

### G. Coding error and identifiability paper

**Claim:** Small rates of lexical coding error can dominate shallow divergence estimates, while generic list-error parameters are themselves confounded with terminal branch length.

**Existing evidence:** Dramatic simulated-world failures and the Paumari historical-list case.

**Missing:** Real replicated codings, identifiable validation datasets, analytic treatment of the confounding, hierarchical error models informed by duplicated lists, and comparisons with coder-aware Bayesian models.

**Venues:** *Journal of Language Evolution*, *Linguistic Typology* methods venues, *Open Research Europe*.

**Strongest objection:** The current result relies largely on one stylized simulation mechanism.

## 7. Recommendation

1. **Make the empirical calibration benchmark the project’s critical path.** Expand LOCO to more families, more depths and properly scored historical windows; compare against prior-only, classical methods, ASJP, and at least one BEAST/BEASTling or RevBayes workflow. Until then, avoid general “well-calibrated” claims.

2. **Address the two demonstrated model failures before adding breadth:** time-dependent rates and coding/list error. Add at least a depth-trend sensitivity model and replicated-coding experiments. These failures are much more consequential than incremental sampler or UI features.

3. **Harden topology propagation.** Stop pooling unconverged tree runs, state the modular interpretation of external tree samples, and show converged-only sensitivity. Longer term, define precisely when topology samples can be combined without incoherent reuse of lexical evidence.

4. **Publish the software/methodological transparency contribution first.** The implementation and validation of the T1 likelihood, lexical/calendar separation, prior overlays, and negative-results culture are already the strongest story. A software paper or focused methods paper is supportable before a claim of superior calendar dating.

5. **Reframe the flagship scientific promise.** The defensible promise today is: “well-audited posterior lexical chronologies and honest calendar scenarios, with diagnostics that reveal when calendar information is weak.” If broader validation later demonstrates nominal coverage and useful skill, promote calendar dating back into the headline.

Overall: this is now a serious, unusually transparent research prototype and a plausible future reference implementation. It is not yet a reference empirical dating method. Its greatest strength is that its own evidence makes that distinction visible.
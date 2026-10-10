# Methods from outside linguistics for calibrated node-age estimation: annotated bibliography and recommendations for clepsydra

Compiled 2026-10-10 for `~/work/chl/clepsydra` (v0.2.0).

**Scope.** These are methods from molecular-clock dating, phylodynamics, archaeological
chronology, stemmatology and Bayesian computation that bear on one goal: a
well-calibrated posterior for every node age from cognate data, on a given tree or
without one, including families that have no calibrations of their own.

**What clepsydra already does, as read from `docs/` for this review.** It runs on a
fixed topology. It has T1 (strict novel-state) and T2 (binary CTMC) likelihoods, a
concept hierarchy `log r_{f,c} = μ_c + β_f + u_{f,c}`, and an iid lognormal relaxed
clock with a HalfNormal `σ_b`. Calibrations are priors on named nodes and tips. It also
has prior-only sampling, leave-one-calibration-out (LOCO), SBC, posterior predictive
checks (PPC) in Python, and contamination and loan mixtures. Uncalibrated node ages
get a flat prior bounded only by order constraints, and there is no tree prior.

**How the DOIs were checked.** Every DOI below was resolved against the Crossref API on
2026-10-10, and the title, journal, year and authors matched. Items that have no DOI
(arXiv, JMLR, ICLR, J-SFdS) are marked **[no DOI; verified by web search]**. Statements
about what a paper says come from the abstract (Crossref or PubMed) where it was
available. Everything else comes from my working knowledge of these papers and is marked
*(from knowledge; check before quoting)* wherever the detail is specific enough to matter.

---

## Part A — Annotated bibliography

### A1. Dating software and tree priors

**Drummond, A.J., Ho, S.Y.W., Phillips, M.J. & Rambaut, A. (2006). Relaxed phylogenetics and dating with confidence. *PLoS Biology* 4(5): e88. doi:10.1371/journal.pbio.0040088**
- *What it says.* Introduces the uncorrelated lognormal (UCLN) and uncorrelated exponential relaxed clocks. Branch rates are iid draws, and time and topology are co-estimated in BEAST.
- *Transfer.* This is clepsydra's existing relaxed clock. The paper parameterizes the UCLN with a mean-rate parameter and a separate SD, as clepsydra does with `E[m_b]=1`. Its UCLN simulations also show that a relaxed clock, unlike a strict one, gives interval coverage close to nominal when rates really do vary. That makes it the baseline clock for the coverage studies in Recommendation 7.

**Drummond, A.J. & Rambaut, A. (2007). BEAST: Bayesian evolutionary analysis by sampling trees. *BMC Evolutionary Biology* 7: 214. doi:10.1186/1471-2148-7-214**; **Bouckaert, R. et al. (2019). BEAST 2.5. *PLoS Computational Biology* 15: e1006650. doi:10.1371/journal.pcbi.1006650**; **Suchard, M.A. et al. (2018). BEAST 1.10. *Virus Evolution* 4: vey016. doi:10.1093/ve/vey016**
- *Transfer.* These are the reference implementations of most priors and operators discussed below. When clepsydra needs an independent cross-check (for example, binary T2 against BEAST 2 with a binary CTMC and a fixed tree), these are the tools to use.

**Höhna, S. et al. (2016). RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language. *Systematic Biology* 65: 726–736. doi:10.1093/sysbio/syw021**
- *Transfer.* RevBayes treats the time-tree prior, the clock and the calibrations as separate nodes of a graphical model, and its tutorials always sample the prior first. Clepsydra can adopt that discipline: every prior component should be swappable and testable on its own.

**Yang, Z. (2007). PAML 4. *Molecular Biology and Evolution* 24: 1586–1591. doi:10.1093/molbev/msm088** (MCMCTree)
**Yang, Z. & Rannala, B. (2006). Bayesian estimation of species divergence times under a molecular clock using multiple fossil calibrations with soft bounds. *MBE* 23: 212–226. doi:10.1093/molbev/msj024**
- *What it says.* Defines soft bounds (a small tail probability beyond the minimum and maximum), and defines the time prior by **conditional construction**. Calibrated nodes get the user's densities. Uncalibrated node ages are then drawn from a birth–death-with-species-sampling kernel *conditional* on the calibrated ages, and the result is truncated to respect node order.
- *Transfer.* This is the most direct template for replacing clepsydra's flat prior. Under it, uncalibrated nodes have a proper, interpretable prior, while calibrated nodes stay close to (but not exactly at) the user's densities.

**Yang, Z. & Rannala, B. (1997). Bayesian phylogenetic inference using DNA sequences: a Markov chain Monte Carlo method. *MBE* 14: 717–724. doi:10.1093/oxfordjournals.molbev.a025811**
- *Transfer.* This is the source of the birth–death-with-sampling kernel density for node ages. Conditional on the root age, the uncalibrated node ages are order statistics of iid draws from a closed-form density with parameters (λ, μ, ρ). For a fixed topology that density is all clepsydra needs, and it is a few lines of C *(from knowledge; check before quoting)*.

**Gernhard, T. (2008). The conditioned reconstructed process. *Journal of Theoretical Biology* 253: 769–778. doi:10.1016/j.jtbi.2008.04.005**
- *Transfer.* Gives node-age densities under birth–death conditioned on the number of extant taxa and on the root (or origin) age. It is the theoretical basis for "node ages given root age are iid order statistics", which makes a fixed-topology birth–death prior cheap to evaluate.

**Stadler, T. (2009). On incomplete sampling under birth–death models and connections to the sampling-based coalescent. *J. Theor. Biol.* 261: 58–66. doi:10.1016/j.jtbi.2009.07.018**; **Höhna, S., Stadler, T., Ronquist, F. & Britton, T. (2011). Inferring speciation and extinction rates under different sampling schemes. *MBE* 28: 2577–2589. doi:10.1093/molbev/msr095**
- *Transfer.* Language samples are rarely complete, and they are often *diversified*: one language per subgroup is chosen deliberately. Under-sampled trees put more prior mass near the tips. Höhna et al. show that the sampling scheme (uniform or diversified) changes the implied node-age prior. Clepsydra's tree prior therefore needs a sampling-fraction ρ and a choice of scheme, and should not silently assume complete sampling.

**Stadler, T. (2010). Sampling-through-time in birth–death trees. *J. Theor. Biol.* 267: 396–404. doi:10.1016/j.jtbi.2010.09.010**
- *What it says.* Gives the birth–death–sampling (BDS) tree density when tips are sampled through time, which is what is needed for dated non-contemporaneous tips.
- *Transfer.* Ancient languages (Latin, Hittite, Old Irish, Sanskrit, Tocharian) are serially sampled tips. Once dated tips enter the analysis, BDS or FBD is the coherent tree prior.

**Heath, T.A., Huelsenbeck, J.P. & Stadler, T. (2014). The fossilized birth–death process for coherent calibration of divergence-time estimates. *PNAS* 111: E2957–E2966. doi:10.1073/pnas.1319091111**
- *What it says.* The fossilized birth–death (FBD) process treats fossils as samples from the same process that generated the tree, so they are not ad hoc calibration densities. Calibration then follows from the model without user-chosen shapes.
- *Transfer.* Attested ancient languages and dated inscriptions can be modelled as FBD samples. This avoids the "which lognormal do I put on the Italic node" problem whenever the evidence is actually an attested descendant.

**Gavryushkina, A., Welch, D., Stadler, T. & Drummond, A.J. (2014). Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. *PLoS Comput Biol* 10: e1003919. doi:10.1371/journal.pcbi.1003919**
- *What it says.* Defines sampled-ancestor (SA) trees, in which a sampled lineage can be a direct ancestor (a zero-length pendant edge), and gives the MCMC moves that jump between the ancestor and tip states.
- *Transfer.* Latin is plausibly ancestral, or nearly so, to the Romance languages. Forcing it to be a sister tip pushes the Romance root older than Latin's date. On a fixed topology, an "ancestor indicator" per ancient tip is a reversible-jump move on one pendant edge.

**Gavryushkina, A., Heath, T.A., Ksepka, D.T., Stadler, T., Welch, D. & Drummond, A.J. (2017). Bayesian total-evidence dating reveals the recent crown radiation of penguins. *Syst Biol* 66: 57–73. doi:10.1093/sysbio/syw060**
- *Transfer.* A worked FBD and sampled-ancestor analysis. It shows how the crown age moves once fossils can be ancestors.

**Heled, J. & Drummond, A.J. (2012). Calibrated tree priors for relaxed phylogenetics and divergence time estimation. *Syst Biol* 61: 138–149. doi:10.1093/sysbio/syr087**
- *What it says (abstract).* Calibration densities as implemented in BEAST 1.5 "do not represent the marginal prior distribution of the calibration node". The paper shows this analytically on small trees. It constructs a calibrated Yule prior whose marginal on one calibrated node equals the user's density, with or without monophyly, by computing the Yule prior conditional on the calibrated height. The multi-calibration case was left "elusive", and the paper concludes that "direct estimation of the prior induced by specifying multiple calibration densities should be a prerequisite of any divergence time dating analysis".
- *Transfer.* This is the central reference for Recommendation 2. Multiplying a tree prior by calibration densities (the "multiplicative construction") does not give the calibrations you wrote down. Clepsydra must always estimate and display the effective prior.

**Heled, J. & Drummond, A.J. (2015). Calibrated birth–death phylogenetic time-tree priors for Bayesian inference. *Syst Biol* 64: 369–383. doi:10.1093/sysbio/syu089**
- *What it says (abstract).* Defines a class of multi-calibration birth–death priors in which the calibrated nodes have exactly the user-specified marginal, while the birth–death density ratio is kept for uncalibrated nodes. One formulation factorizes times from ranked topologies and stays efficient for up to about six calibrations.
- *Transfer.* For a *fixed* topology the ranked-topology part drops out, which makes the factorized formulation easier to implement. It is the "exact marginal" option to set beside MCMCTree-style conditional construction.

**Rannala, B. (2016). Conceptual issues in Bayesian divergence time estimation. *Phil. Trans. R. Soc. B* 371: 20150134. doi:10.1098/rstb.2015.0134**
- *What it says (abstract).* Times are not identifiable without both calibrations and data. Marginal calibration priors can conflict with node order, which changes the effective prior for a given topology, and this must not be mistaken for signal in the data. The paper defines a topology-consistent prior that preserves the marginals. Its most important warning is that conflicts between calibrations and relative branch lengths can give **grossly incorrect dates with narrow posteriors**, so overly narrow posteriors must be scrutinized.
- *Transfer.* Clepsydra should flag "narrow but conflicted" nodes, where the posterior is narrow and calibration tension is high (Recommendations 9 and 24). Its truncation diagnostics should also separate order-induced prior shifts from shifts caused by the data.

**Inoue, J., Donoghue, P.C.J. & Yang, Z. (2010). The impact of the representation of fossil calibrations on Bayesian estimation of species divergence times. *Syst Biol* 59: 74–89. doi:10.1093/sysbio/syp078**
- *What it says (abstract).* Of all the factors examined, the time prior (including the calibrations) had the largest impact on posterior times. MULTIDIVTIME and MCMCTree implement minimum and maximum bounds differently, which produced very different time priors and explained large differences between published estimates.
- *Transfer.* The same nominal calibration table can produce different results in different implementations. Clepsydra should document its construction exactly and report the effective prior.

**Barba-Montoya, J., dos Reis, M. & Yang, Z. (2017). Comparison of different strategies for using fossil calibrations to generate the time prior in Bayesian molecular clock dating. *Mol. Phylogenet. Evol.* 114: 386–400. doi:10.1016/j.ympev.2017.07.005**
- *What it says (abstract).* Evaluates three strategies for turning minimum and maximum bounds into a time prior. They differ in whether a node borrows the maximum bounds of its ancestors and the minimum bounds of its descendants. The comparison covers MCMCTree, MrBayes and BEAST 2, on an analytic example and two real datasets. *(Detail from knowledge; check before quoting:* the strategies give markedly different effective priors, and therefore different posteriors.)
- *Transfer.* Clepsydra's feasibility propagation (`docs/calibrations.md` §2) already computes inherited bounds. The open question is whether those bounds enter the *prior* or only the support for proposals. That should be decided explicitly and documented.

**Warnock, R.C.M., Yang, Z. & Donoghue, P.C.J. (2012). Exploring uncertainty in the calibration of the molecular clock. *Biology Letters* 8: 156–159. doi:10.1098/rsbl.2011.0710**
- *Transfer (from knowledge; check before quoting).* Different but defensible ways of coding the same fossil evidence as densities change the posterior a lot. This argues for reporting sensitivity to how calibrations are encoded, and clepsydra's `calib-sensitivity` already does part of this.

**Warnock, R.C.M., Parham, J.F., Joyce, W.G., Lyson, T.R. & Donoghue, P.C.J. (2015). Calibration uncertainty in molecular dating analyses: there is no substitute for the prior evaluation of time priors. *Proc. R. Soc. B* 282: 20141013. doi:10.1098/rspb.2014.1013**
- *What it says (abstract).* A posteriori consistency methods (cross-validation) "can lead to the selection of erroneous calibrations". Posteriors are "extremely sensitive to the probabilistic interpretation of temporal constraints". Effective priors for one calibration differ depending on which others are used with it, which breaks the assumption behind every consistency method that a calibration has the same impact alone or in combination.
- *Transfer.* This is an important counterweight to Recommendation 9. Cross-validation of calibrations is a *diagnostic*, never an automatic filter. The default response to an inconsistent calibration should be the outlier model (Recommendation 3), not deletion.

**dos Reis, M., Donoghue, P.C.J. & Yang, Z. (2016). Bayesian molecular clock dating of species divergences in the genomics era. *Nature Reviews Genetics* 17: 71–80. doi:10.1038/nrg.2015.8**
- *What it says (from knowledge; check before quoting).* A review covering rate–time confounding, the infinite-sites limit, the dominance of the time prior, the difference between specified and effective priors, partitioning and the rate prior, and approximate likelihood.
- *Transfer.* The best single overview for designing clepsydra's prior and diagnostics.

**Bromham, L., Duchêne, S., Hua, X., Ritchie, A.M., Duchêne, D.A. & Ho, S.Y.W. (2018). Bayesian molecular dating: opening up the black box. *Biological Reviews* 93: 1165–1191. doi:10.1111/brv.12390**; **Ho, S.Y.W. & Duchêne, S. (2014). Molecular-clock methods for estimating evolutionary rates and timescales. *Molecular Ecology* 23: 5947–5965. doi:10.1111/mec.12953**
- *Transfer.* These are practitioner checklists: prior-only runs, clock-model comparison, calibration sensitivity, and checking the dates against independent information. Clepsydra's "required sensitivity procedures" largely match them already.

**Brown, J.W. & Smith, S.A. (2018). The past sure is tense: on interpreting phylogenetic divergence time estimates. *Syst Biol* 67: 340–353. doi:10.1093/sysbio/syx074**
- *Transfer (from knowledge; check before quoting).* Dates are deeply prior-conditioned, and their precision is often illusory. The paper is a good reference for the wording of clepsydra's model cards and output caveats.

**Ronquist, F., Klopfstein, S., Vilhelmsen, L., Schulmeister, S., Murray, D.L. & Rasnitsyn, A.P. (2012). A total-evidence approach to dating with fossils, applied to the early radiation of the Hymenoptera. *Syst Biol* 61: 973–999. doi:10.1093/sysbio/sys058**
- *What it says (from knowledge; check before quoting).* Tip-dating with morphological characters, and the **uniform node-age (uniform clock) tree prior** in MrBayes. Node ages are uniform conditional on the root age, and the root age has its own explicit prior.
- *Transfer.* The uniform node-age prior is the closest proper relative of clepsydra's current flat prior. The decisive difference is that the root age needs an explicit proper prior, and node ages are uniform *given* the root. See Recommendation 1.

**Zhang, C., Stadler, T., Klopfstein, S., Heath, T.A. & Ronquist, F. (2016). Total-evidence dating under the fossilized birth–death process. *Syst Biol* 65: 228–249. doi:10.1093/sysbio/syv080**
- *Transfer (from knowledge; check before quoting).* Adds diversified-sampling FBD. The extant-sampling scheme strongly affects dates, so the representative-sampling design of a cognate database (one doculect per subgroup) should be modelled as diversified sampling.

**Matzke, N.J. & Wright, A. (2016). Inferring node dates from tip dates in fossil Canidae: the importance of tree priors. *Biology Letters* 12: 20160328. doi:10.1098/rsbl.2016.0328**
- *Transfer (from knowledge; check before quoting).* When dating relies on tip dates, the tree prior has a strong effect on node ages, and naive or uniform priors can be badly biased. This argues for testing tip-dating in clepsydra under BDS or FBD rather than under a flat prior.

**Stadler, T. & Yang, Z. (2013). Dating phylogenies with sequentially sampled tips. *Syst Biol* 62: 674–688. doi:10.1093/sysbio/syt030**
- *Transfer.* Gives the theory of dating from tip dates alone and the conditional construction with BDS. It supports adding dated ancient-language tips as a calibration source.

**Ritchie, A.M., Lo, N. & Ho, S.Y.W. (2017). The impact of the tree prior on molecular dating of data sets containing a mixture of inter- and intraspecies sampling. *Syst Biol* 66: 413–425. doi:10.1093/sysbio/syw095**
- *Transfer.* Mis-chosen tree priors bias node ages when sampling density is uneven. Language datasets mix densely sampled dialect clusters with sparsely sampled branches, so the tree prior needs a sensitivity check (Yule against birth–death against uniform).

### A2. Fast and penalized-likelihood dating (non-Bayesian cross-checks)

**Sanderson, M.J. (2002). Estimating absolute rates of molecular evolution and divergence times: a penalized likelihood approach. *MBE* 19: 101–109. doi:10.1093/oxfordjournals.molbev.a003974**; **Smith, S.A. & O'Meara, B.C. (2012). treePL: divergence time estimation using penalized likelihood for large phylogenies. *Bioinformatics* 28: 2689–2690. doi:10.1093/bioinformatics/bts492**
- *What they say.* Rates are allowed to change between adjacent branches, with a roughness penalty whose smoothing weight is chosen by cross-validation.
- *Transfer.* Penalized likelihood is cheap and dependency-free, and makes a good initializer and sanity check (Recommendation 20). Its cross-validated smoothing is a precursor of fossil-based model cross-validation (Near & Sanderson 2004).

**Paradis, E. (2013). Molecular dating of phylogenies by likelihood methods: a comparison of models and a new information criterion. *Mol. Phylogenet. Evol.* 67: 436–444. doi:10.1016/j.ympev.2013.02.008** (`chronos` in ape)
- *Transfer.* Offers strict, correlated and discrete clock models under penalized likelihood, with an information criterion for choosing between them. The discrete model (a few rate categories) is a cheap middle ground between strict and fully relaxed clocks.

**To, T.-H., Jung, M., Lycett, S. & Gascuel, O. (2016). Fast dating using least-squares criteria and algorithms. *Syst Biol* 65: 82–97. doi:10.1093/sysbio/syv068**; LSD2 ships in **Minh, B.Q. et al. (2020). IQ-TREE 2. *MBE* 37: 1530–1534. doi:10.1093/molbev/msaa015**
- *What it says.* Dating by least squares, with weighting proportional to the expected variance of the branch lengths, in near-linear time. It handles tip dates and node constraints and gives confidence intervals by resampling branch lengths.
- *Transfer.* LSD is the right algorithm for a C99 "quick date" mode and for initial node ages. Its variance weighting (variance ∝ branch length plus a constant) carries over directly to lexical branch lengths in units of expected replacements per meaning.

**Sagulenko, P., Puller, V. & Neher, R.A. (2018). TreeTime: maximum-likelihood phylodynamic analysis. *Virus Evolution* 4: vex042. doi:10.1093/ve/vex042**
- *Transfer (from knowledge; check before quoting).* Joint marginal maximum-likelihood dating using message passing over discretized node-time distributions, which gives approximate marginal posteriors in linear time. The messages are convolutions of branch-length likelihood with time. Clepsydra already has a `gd_grid` module, and a grid-based belief-propagation pass on a fixed tree with a strict clock could give near-exact marginals for small trees. That is useful as an oracle for MCMC validation.

**Tamura, K. et al. (2012). Estimating divergence times in large molecular phylogenies. *PNAS* 109: 19333–19338. doi:10.1073/pnas.1213199109** (RelTime)
**Battistuzzi, F.U., Billing-Ross, P., Murillo, O., Filipski, A. & Kumar, S. (2015). A protocol for diagnosing the effect of calibration priors on posterior time estimates. *MBE* 32: 1907–1912. doi:10.1093/molbev/msv075**
- *What it says (abstract).* Uses calibration-free relative times (RelTime) to check calibration priors. In the metazoan case, two maximum-only calibrations alone were found to force young dates that "severely violate minimum boundaries of almost all of the other 22 calibration constraints".
- *Transfer.* This protocol needs only relative ages. Clepsydra can produce those itself by running with a strict-scale anchor and no calibrations. Each calibration can then be compared against the relative chronology anchored by the others (Recommendation 9).

**Rambaut, A., Lam, T.T., Carvalho, L.M. & Pybus, O.G. (2016). Exploring the temporal structure of heterochronous sequences using TempEst. *Virus Evolution* 2: vew007. doi:10.1093/ve/vew007**
- *Transfer.* Root-to-tip regression against tip date as a check for clock-like signal before Bayesian analysis. With dated ancient languages, the same regression (lexical distance from the root against the tip's attested date) is a cheap pre-check for temporal signal.

### A3. Clock models

**Thorne, J.L., Kishino, H. & Painter, I.S. (1998). Estimating the rate of evolution of the rate of molecular evolution. *MBE* 15: 1647–1657. doi:10.1093/oxfordjournals.molbev.a025892**; **Kishino, H., Thorne, J.L. & Bruno, W.J. (2001). Performance of a divergence time estimation method under a probabilistic model of rate evolution. *MBE* 18: 352–361. doi:10.1093/oxfordjournals.molbev.a003811**
- *Transfer.* The autocorrelated (geometric Brownian) clock. In language, conservatism is plausibly inherited (a conservative ancestor has conservative daughters, as with Icelandic or Lithuanian), which argues for testing an autocorrelated option.

**Rannala, B. & Yang, Z. (2007). Inferring speciation times under an episodic molecular clock. *Syst Biol* 56: 453–466. doi:10.1080/10635150701420643**
- *What it says (abstract).* Adds autocorrelated (geometric Brownian) and independent lognormal clocks to MCMCTree. Uncalibrated nodes get a birth–death-with-species-sampling prior. It develops the **infinite-sites theory**: as the amount of sequence data goes to infinity, the posterior CI width and posterior mean of divergence times form a **perfect linear relationship**, and the slope measures the uncertainty that sequence data alone cannot reduce. "Posterior time estimates typically involve considerable uncertainties even with an infinite amount of sequence data."
- *Transfer.* This sets the theoretical ceiling on what a cognate list can tell you (Recommendation 8).

**Lepage, T., Bryant, D., Philippe, H. & Lartillot, N. (2007). A general comparison of relaxed molecular clock models. *MBE* 24: 2669–2680. doi:10.1093/molbev/msm193**
- *Transfer (from knowledge; check before quoting).* Compares uncorrelated, lognormal-autocorrelated and CIR clocks using Bayes factors. Autocorrelated models were often favoured on real data. Clepsydra should not hard-wire the uncorrelated clock without testing alternatives.

**Lartillot, N., Phillips, M.J. & Ronquist, F. (2016). A mixed relaxed clock model. *Phil. Trans. R. Soc. B* 371: 20150132. doi:10.1098/rstb.2015.0132**
- *What it says (abstract).* Combines Brownian long-term trends with white-noise short-term fluctuations. On mammals the mixed model beats both pure clocks with roughly equal contributions from each. Tip-dating "is particularly sensitive to the choice of the relaxed clock model", and the pure Brownian clock is "overly rigid".
- *Transfer.* The mixed clock is cheap to add on a fixed tree: a log-rate random walk along the tree plus an iid per-branch term. It covers both the inherited-conservatism story and episodic contact-driven bursts. Defer it until tip-dating is in place.

**Drummond, A.J. & Suchard, M.A. (2010). Bayesian random local clocks, or one rate to rule them all. *BMC Biology* 8: 114. doi:10.1186/1741-7007-8-114**
- *Transfer.* Rate shifts on a few branches, with the strict clock as the simplest case, sampled by a spike-and-slab or indicator scheme. This suits a language subgroup with a known episode of heavy replacement, such as a taboo-driven or contact-driven lexical overhaul.

**Fisher, A.A., Ji, X., Nishimura, A., Baele, G., Lemey, P. & Suchard, M.A. (2023). Shrinkage-based random local clocks with scalable inference. *MBE* 40(11): msad242. doi:10.1093/molbev/msad242**
- *What it says (abstract).* An autocorrelated model of heritable rates. Increments between neighbouring branches get heavy-tailed, mean-zero shrinkage priors (a Bayesian bridge). An HMC sampler with closed-form gradients makes it faster than the random local clock (RLC), and it recovers known local clocks.
- *Transfer.* A continuous, easier-to-sample replacement for random local clocks. In a C engine without HMC it can still be sampled by MH on increments (a fixed tree has few branches). Medium priority.

**Douglas, J., Zhang, R. & Bouckaert, R. (2021). Adaptive dating and fast proposals: revisiting the phylogenetic relaxed clock model. *PLoS Comput Biol* 17: e1008322. doi:10.1371/journal.pcbi.1008322**
- *What it says (abstract).* Compares three rate parameterizations (real rates, categories, quantiles) and introduces an adaptive operator that learns operator weights. Adds AVMVN (adaptive multivariate normal) and Bactrian proposal kernels and relaxed-clock-aware narrow exchange. Together these give "up to 65 times more efficient" mixing of relaxed-clock parameters. Implemented as ORC in BEAST 2.
- *Transfer.* Recommendations 6 and 12. These operators are cheap, dependency-free, and attack the same rate–time ridge that clepsydra's documentation (§7.2) already flags.

**Zhang, R. & Drummond, A.J. (2020). Improving the performance of Bayesian phylogenetic inference under relaxed clock models. *BMC Evol Biol* 20: 54. doi:10.1186/s12862-020-01609-4**
- *What it says (abstract).* A proposal kernel that changes a node's age and its three adjacent branch rates *together*, keeping the implied genetic distances constant. Three root strategies (Simple Distance, Small Pulley, Big Pulley) are discussed. It gives large efficiency gains.
- *Transfer.* This is the single most valuable operator for clepsydra (Recommendation 6). The likelihood is unchanged by the move, so only the prior ratio and the Jacobian need computing, which makes it nearly free.

**Bouckaert, R.R. & Drummond, A.J. (2017). bModelTest: Bayesian phylogenetic site model averaging and model comparison. *BMC Evol Biol* 17: 42. doi:10.1186/s12862-017-0890-6**; **Li, W.L.S. & Drummond, A.J. (2012). Model averaging and Bayes factor calculation of relaxed molecular clocks in Bayesian phylogenetics. *MBE* 29: 751–761. doi:10.1093/molbev/msr232**
- *Transfer.* Model averaging across clock models by reversible jump, so dates are not conditioned on one clock choice. In clepsydra, the HalfNormal on `σ_b` already shrinks toward a strict clock. A point mass at `σ_b=0` (spike-and-slab) would give an explicit posterior probability for the strict clock.

**Heath, T.A., Holder, M.T. & Huelsenbeck, J.P. (2012). A Dirichlet process prior for estimating lineage-specific substitution rates. *MBE* 29: 939–955. doi:10.1093/molbev/msr255**
- *Transfer.* Clusters branch rates into an unknown number of rate classes. The same machinery can cluster *families* into rate classes when rate information is shared across families (Recommendation 4).

#### Time-dependent rates (relevant to lexical saturation at deep time)

**Ho, S.Y.W., Phillips, M.J., Cooper, A. & Drummond, A.J. (2005). Time dependency of molecular rate estimates and systematic overestimation of recent divergence times. *MBE* 22: 1561–1568. doi:10.1093/molbev/msi145**; **Ho, S.Y.W. et al. (2011). Time-dependent rates of molecular evolution. *Mol. Ecol.* 20: 3087–3101. doi:10.1111/j.1365-294X.2011.05178.x**; **Ho, S.Y.W., Duchêne, S., Molak, M. & Shapiro, B. (2015). Time-dependent estimates of molecular evolutionary rates: evidence and causes. *Mol. Ecol.* 24: 6007–6012. doi:10.1111/mec.13450**
- *Transfer.* Rates calibrated on shallow nodes are faster than rates calibrated on deep nodes. Applying a shallow-calibrated rate to deep nodes underestimates their ages, and the reverse overestimates shallow ages. Most language calibrations are shallow (historically attested splits under 2.5 kyr), while many questions of interest are deep. This is the main systematic risk in transferring rates across families and depths.

**Aiewsakun, P. & Katzourakis, A. (2016). Time-dependent rate phenomenon in viruses. *Journal of Virology* 90: 7184–7195. doi:10.1128/JVI.00593-16**
- *Transfer (from knowledge; check before quoting).* Across virus datasets, the log rate falls linearly with the log of the timescale of measurement (a power law). The paper proposes correcting dates with this empirical relationship. Clepsydra can test the same relationship across its own calibrations: compute the implied rate for each calibration and regress the log rate on the log calibration age (Recommendation 11).

**Membrebe, J.V., Suchard, M.A., Rambaut, A., Baele, G. & Lemey, P. (2019). Bayesian inference of evolutionary histories under time-dependent substitution rates. *MBE* 36: 1793–1803. doi:10.1093/molbev/msz094**
- *What it says (abstract).* A formal Bayesian clock in which the rate is a function of time (epochs). On foamy viruses and lentiviruses it finds a strong time-dependent effect, with rates varying over four orders of magnitude.
- *Transfer.* An epoch clock (a piecewise-constant rate multiplier by absolute-time interval) is simple to put into a fixed-tree likelihood: each branch's effective length is the integral of the rate over the branch's time span. It is identifiable only with calibrations at several depths.

**Ghafari, M., Simmonds, P., Pybus, O.G. & Katzourakis, A. (2021). A mechanistic evolutionary model explains the time-dependent pattern of substitution rates in viruses. *Current Biology* 31: 4689–4696.e5. doi:10.1016/j.cub.2021.08.020** (preprint title: "Prisoner of War dynamics explains the time-dependent pattern…", bioRxiv doi:10.1101/2021.02.09.430479)
- *What it says (abstract).* A mechanistic model (Prisoner of War, PoW) that reproduces the power-law rate decay with slope about −0.65. Standard substitution models "fail to correctly estimate divergence times once the most rapidly evolving sites saturate". Re-estimated dates are much older, for example an HCV genotype diversification about 423 kya and a sarbecovirus MRCA about 21 kya, "nearly thirty times older".
- *Transfer.* The mechanism is **among-site rate heterogeneity combined with saturation of the fast sites**, which maps directly onto lexical meanings with very different replacement rates (Pagel et al. 2007). A model with a well-specified *distribution* of per-concept rates (clepsydra's `μ_c` hierarchy) already captures the PoW mechanism in principle. The risks are (i) a truncated or misshapen rate distribution, especially a missing very slow tail, and (ii) borrowing, which plays the role of homoplasy. Recommendation 11 turns this into tests.

### A4. Precision limits: infinite sites

**dos Reis, M. & Yang, Z. (2013). The unbearable uncertainty of Bayesian divergence time estimation. *J. Syst. Evol.* 51: 30–43. doi:10.1111/j.1759-6831.2012.00236.x**
- *Transfer (from knowledge; check before quoting).* Gives a finite-sites theory, quantifies how fast uncertainty approaches the infinite-sites limit, and introduces the "infinite-sites plot" (posterior CI width against posterior mean, with the R² of the fit indicating closeness to the limit).

**Zhu, T., dos Reis, M. & Yang, Z. (2015). Characterization of the uncertainty of divergence time estimation under relaxed molecular clock models using multiple loci. *Syst Biol* 64: 267–280. doi:10.1093/sysbio/syu109**
- *What it says (abstract and knowledge).* Times and rates are confounded, so "posterior time estimates will not approach point values even if an infinite amount of sequence data are used". Under relaxed clocks, uncertainty falls with the number of *loci*, each with its own rate draws, not with the number of sites per locus. In the infinite-sites-and-infinite-loci limit only the calibration uncertainty remains *(the multi-locus point is from knowledge; check before quoting)*.
- *Transfer.* This matters for how clepsydra's relaxed clock is structured. If one set of branch multipliers `m_b` is shared by all concepts, adding concepts **cannot** remove the uncertainty from that one realization of branch-rate noise. Each concept adds sites, not loci. So under clepsydra's present structure, the floor on precision is set by calibrations plus `σ_b`, however long the list. Any claim that "more concepts give tighter dates" must be checked against this (Recommendation 8).

**dos Reis, M. et al. (2012). Phylogenomic datasets provide both precision and accuracy in estimating the timescale of placental mammal phylogeny. *Proc. R. Soc. B* 279: 3491–3500. doi:10.1098/rspb.2012.0683**
- *Transfer.* An applied infinite-sites plot showing that genome-scale data reached the calibration-limited regime. This is the template for clepsydra's report.

**dos Reis, M., Zhu, T. & Yang, Z. (2014). The impact of the rate prior on Bayesian estimation of divergence times with multiple loci. *Syst Biol* 63: 555–565. doi:10.1093/sysbio/syu020**
- *What it says (abstract).* With iid priors on locus rates, the prior variance of the *average* rate shrinks as 1/L. The rate prior then dominates, and if it is misspecified the times "converge to wrong values with very narrow credibility intervals". The fix is a gamma-Dirichlet prior: a gamma on the mean rate, split across loci by a Dirichlet.
- *Transfer.* **This applies directly to the concept hierarchy.** If clepsydra's per-concept rates `μ_c` (about 100–200 concepts) were iid around a *fixed* centre, the implied prior on the mean lexical rate would collapse, and the reference-prior centre would silently fix absolute time for uncalibrated families. The global clock scale must carry the uncertainty of the mean rate, the `μ_c` must be deviations with a centring constraint or a Dirichlet-type split, and the global scale's prior must be the honest cross-family predictive (Recommendation 4).

**Rannala, B., Zhu, T. & Yang, Z. (2012). Tail paradox, partial identifiability, and influential priors in Bayesian branch length inference. *MBE* 29: 325–335. doi:10.1093/molbev/msr210**
- *Transfer.* iid priors on many branch parameters can induce a highly informative, unintended prior on their sum (tree length). The same induced-prior logic applies to clepsydra's flat node ages and to the `m_b` (Recommendation 1).

**Graur, D. & Martin, W. (2004). Reading the entrails of chickens: molecular timescales of evolution and the illusion of precision. *Trends in Genetics* 20: 80–86. doi:10.1016/j.tig.2003.12.003**; **Schenk, J.J. (2016). Consequences of secondary calibrations on divergence time estimates. *PLoS ONE* 11: e0148228. doi:10.1371/journal.pone.0148228**
- *What Schenk says (abstract).* When a primary posterior is reused as a secondary calibration, the 95% CIs come out "significantly younger and narrower", primary and secondary estimates differ in 97% of replicates, and normal-shaped secondary priors are worse than uniform ones.
- *Transfer.* Clepsydra must never take a published glottochronological or phylogenetic date as a node calibration without inflating its uncertainty. It is better to model the families jointly (Recommendation 4) than to chain posteriors.

### A5. Tip-dating and total evidence

**Ronquist, F., Lartillot, N. & Phillips, M.J. (2016). Closing the gap between rocks and clocks using total-evidence dating. *Phil. Trans. R. Soc. B* 371: 20150136. doi:10.1098/rstb.2015.0136**
- *Transfer (from knowledge; check before quoting).* Tip-dating often gives older ages than node-dating, and the cause is often the clock model and the tree prior rather than the data. This is the warning to heed when clepsydra adds dated ancient languages.

**Barido-Sottani, J., van Tiel, N.M.A., Hopkins, M.J., Wright, D.F., Stadler, T. & Warnock, R.C.M. (2020). Ignoring fossil age uncertainty leads to inaccurate topology and divergence time estimates in time calibrated tree inference. *Frontiers in Ecology and Evolution* 8: 183. doi:10.3389/fevo.2020.00183**
- *Transfer (from knowledge; check before quoting).* Fixing tip ages at point values when they are really intervals biases node ages and topology. Clepsydra's `fixed_tip` should get a `uniform_tip` or `normal_tip` sibling that samples the tip age (Recommendation 10). Old Irish, Hittite and Tocharian B are attested across centuries, not at a point in time.

Linguistic applications of these ideas, given here for orientation only: **Chang, W., Cathcart, C., Hall, D. & Garrett, A. (2015). Ancestry-constrained phylogenetic analysis supports the Indo-European steppe hypothesis. *Language* 91: 194–244. doi:10.1353/lan.2015.0005** (forcing ancient languages to be ancestors moves the root age a lot), and **Heggarty, P. et al. (2023). Language trees with sampled ancestors support a hybrid model for the origin of Indo-European languages. *Science* 381: eabg0818. doi:10.1126/science.abg0818** (sampled ancestors via FBD in a lexical analysis).

### A6. Archaeological Bayesian chronology (robust calibrations)

**Bronk Ramsey, C. (2009a). Bayesian analysis of radiocarbon dates. *Radiocarbon* 51: 337–360. doi:10.1017/S0033822200033865**
- *What it says (from knowledge; check before quoting).* OxCal's framework. Dates are grouped into **Phases** (unordered events between a start boundary and an end boundary) and **Sequences** (ordered events). The boundaries get uniform-span priors, and the paper stresses that a naive uniform prior on each event's date biases the span of a group outwards. That bias is why boundaries are needed: they give a prior on the *span* that does not grow with the number of dated events.
- *Transfer.* This is the same "induced prior grows with n" problem as flat node ages (Recommendation 1). OxCal's boundary construction is the archaeological counterpart of the uniform node-age prior conditional on the root, and it is a useful explanatory analogy for users. Phase and Sequence logic also maps onto linguistic "event sets": a set of attestations that must predate a split, or sound changes ordered in relative chronology.

**Bronk Ramsey, C. (2009b). Dealing with outliers and offsets in radiocarbon dating. *Radiocarbon* 51: 1023–1045. doi:10.1017/S0033822200034093**
- *What it says (from knowledge; check before quoting).* Each measurement has a prior outlier probability (often 0.05). If it is an outlier, its date is shifted by a scaled offset drawn from a Student-t (general model), a normal, or a one-sided exponential. The exponential is the "Charcoal" model: material older than the event, as with old wood. The posterior outlier probability is reported for each date, and the offset scale can be shared across dates.
- *Transfer.* **Directly transferable to linguistic calibrations (Recommendation 3).** (i) A general t-outlier for each node calibration, for a contested historical identification. (ii) A one-sided exponential "attestation lag" model: the first attestation of a language postdates its split by an unknown, exponentially distributed lag. This is the formal version of "the first attestation is a minimum age", without needing an ad hoc lognormal offset. (iii) Posterior outlier probabilities per calibration, which are more honest than drop-and-refit cross-validation (Warnock et al. 2015).

**Christen, J.A. (1994). Summarizing a set of radiocarbon determinations: a robust approach. *Applied Statistics* 43: 489–503. doi:10.2307/2986273**; **Christen, J.A. & Pérez E., S. (2009). A new robust statistical model for radiocarbon data. *Radiocarbon* 51: 1047–1059. doi:10.1017/S003382220003410X**
- *Transfer.* The original Bayesian outlier-indicator model, followed by a model with t-distributed errors. The t-error version is the simplest robust calibration density: replace `normal_node` with `t_node(ν)`, which keeps a heavy tail. That is one line of C.

**Lanos, P. & Philippe, A. (2017). Hierarchical Bayesian modeling for combining dates in archeological context. *Journal de la Société Française de Statistique* 158(2): 72–88.** **[no DOI; verified by web search, Numdam JSFS_2017__158_2_72_0]**; **Lanos, P. & Philippe, A. (2018). Event date model: a robust Bayesian tool for chronology building. *Communications for Statistical Applications and Methods* 25: 131–157. doi:10.29220/CSAM.2018.25.2.131** (ChronoModel)
- *What they say.* Several dates of one "event" are combined hierarchically. Each date gets its own individual variance, which automatically down-weights outlying measurements (an "event model") without the user choosing outlier probabilities.
- *Transfer.* When one linguistic node has *several* independent calibration sources (a historical record, an archaeological correlate, an epigraphic date), clepsydra should combine them hierarchically with individual error variances. Multiplying their densities together is wrong, because conflicting sources then produce a falsely narrow product.

**Haslett, J. & Parnell, A. (2008). A simple monotone process with application to radiocarbon-dated depth chronologies. *JRSS C* 57: 399–418. doi:10.1111/j.1467-9876.2008.00623.x**; **Parnell, A.C., Haslett, J., Allen, J.R.M., Buck, C.E. & Huntley, B. (2008). A flexible approach to assessing synchroneity of past events using Bayesian reconstructions of sedimentation history. *Quaternary Science Reviews* 27: 1872–1885. doi:10.1016/j.quascirev.2008.07.009** (Bchron)
- *Transfer.* Monotone compound Poisson–gamma processes for depth–age models. The idea that carries over is joint chronology uncertainty: answer questions such as "was split A before split B?" from joint posterior samples, never from marginal intervals. Clepsydra's outputs should support ordering probabilities for pairs of nodes.

**Bayliss, A. (2009). Rolling out revolution: using radiocarbon dating in archaeology. *Radiocarbon* 51: 123–147. doi:10.1017/S0033822200033750**
- *Transfer.* Describes archaeology's practice of model reporting: prior-only runs, sensitivity analyses, agreement indices, and full publication of model code. That is the standard to match for clepsydra's run manifest.

### A7. Stemmatology and cultural-evolution dating

**Barbrook, A.C., Howe, C.J., Blake, N. & Robinson, P. (1998). The phylogeny of The Canterbury Tales. *Nature* 394: 839. doi:10.1038/29667**; **Spencer, M., Davidson, E.A., Barbrook, A.C. & Howe, C.J. (2004). Phylogenetics of artificial manuscripts. *J. Theor. Biol.* 227: 503–511. doi:10.1016/j.jtbi.2003.11.022**; **Roos, T. & Heikkilä, T. (2009). Evaluating methods for computer-assisted stemmatology using artificial benchmark data sets. *Literary and Linguistic Computing* 24: 417–433. doi:10.1093/llc/fqp002**
- *Transfer.* **Artificial traditions with known truth** are built by having volunteers copy a text, which gives real copying processes with a known stemma. This is the cultural-evolution analogue of SBC with misspecification: data produced by a real process the model does not exactly describe. Clepsydra's validation would gain from a "known-history" benchmark built from documented language histories: Romance with Latin, Germanic with Gothic and Old Norse, Sinitic with Middle Chinese. That is partly what its cross-family validation does already.

**McCollum, J. & Turnbull, R. (2024). Using Bayesian phylogenetics to infer manuscript transmission history. *Digital Scholarship in the Humanities* 39: 258–279. doi:10.1093/llc/fqad089** (Crossref issue year 2024, online 2023)
- *Transfer (abstract not read; check before quoting).* Applies Bayesian phylogenetic dating (BEAST-type, with tip dates from dated manuscripts) to textual traditions. It is a direct analogue of tip-dating with attested ancient languages.

**Kestemont, M. et al. (2022). Forgotten books: the application of unseen species models to the survival of culture. *Science* 375: 765–769. doi:10.1126/science.abl7655**
- *Transfer.* Unseen-species estimators (Chao1 and similar) estimate how many witnesses were lost. For linguistics, the analogue is estimating how many lineages were never sampled, which feeds the sampling-fraction ρ of a birth–death tree prior instead of leaving it as a guess.

**Tehrani, J.J. (2013). The phylogeny of Little Red Riding Hood. *PLoS ONE* 8: e78871. doi:10.1371/journal.pone.0078871**; **da Silva, S.G. & Tehrani, J.J. (2016). Comparative phylogenetic analyses uncover the ancient roots of Indo-European folktales. *Royal Society Open Science* 3: 150645. doi:10.1098/rsos.150645**; **Ross, R.M., Greenhill, S.J. & Atkinson, Q.D. (2013). Population structure and cultural geography of a folktale in Europe. *Proc. R. Soc. B* 280: 20123065. doi:10.1098/rspb.2012.3065**
- *What they say.* Ancestral-state reconstruction of tale presence on a language tree. It dates traits *relative to a dated tree*, with one tale traced to the Bronze Age, and controls for horizontal transmission with spatial autologistic terms.
- *Transfer.* (i) Dates are only as good as the tree's dating, so clepsydra's node-age posteriors should be exportable as *samples* that downstream trait analyses can use. (ii) The test that separates spatial from phylogenetic signal is a model of borrowing that clepsydra's loan mixture could borrow as a diagnostic: does contamination correlate with geography?

### A8. Topology uncertainty and summaries

**Wilson, I.J. & Balding, D.J. (1998). Genealogical inference from microsatellite data. *Genetics* 150: 499–510. doi:10.1093/genetics/150.1.499**; **Drummond, A.J., Nicholls, G.K., Rodrigo, A.G. & Solomon, W. (2002). Estimating mutation parameters, population history and genealogy simultaneously from temporally spaced sequence data. *Genetics* 161: 1307–1320. doi:10.1093/genetics/161.3.1307**; **Lakner, C., van der Mark, P., Huelsenbeck, J.P., Larget, B. & Simon, D. (2008). Efficiency of Markov chain Monte Carlo tree proposals in Bayesian phylogenetics. *Syst Biol* 57: 86–103. doi:10.1080/10635150801886156**; **Höhna, S. & Drummond, A.J. (2012). Guided tree topology proposals for Bayesian phylogenetic inference. *Syst Biol* 61: 1–11. doi:10.1093/sysbio/syr074**
- *Transfer.* The standard time-tree moves: Wilson–Balding prune-and-regraft, subtree-slide (Drummond et al. 2002), narrow and wide exchange, and guided proposals from clade frequencies. For a ultrametric or tip-dated tree with fixed tip ages, **narrow exchange plus subtree-slide plus Wilson–Balding** is the minimal set that works. Lakner et al. show that local moves dominate efficiency.

**Altekar, G., Dwarkadas, S., Huelsenbeck, J.P. & Ronquist, F. (2004). Parallel Metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference. *Bioinformatics* 20: 407–415. doi:10.1093/bioinformatics/btg427**
- *Transfer.* Metropolis-coupled MCMC (MC3, heated chains with swaps). It is easy to add to a 4-chain engine and helps with multimodality caused by calibration conflict or outlier indicators (Recommendation 16).

**Whidden, C. & Matsen, F.A. IV (2015). Quantifying MCMC exploration of phylogenetic tree space. *Syst Biol* 64: 472–491. doi:10.1093/sysbio/syv006**
- *Transfer.* Tree-space "peaks" can trap chains even when standard diagnostics look fine, so clade-split frequencies must be compared across chains (ASDSF) once topology is free.

**Larget, B. (2013). The estimation of tree posterior probabilities using conditional clade probability distributions. *Syst Biol* 62: 501–511. doi:10.1093/sysbio/syt014**
- *Transfer.* The conditional clade distribution (CCD) estimates the probabilities of trees never seen in the sample. It is the basis for the summary below.

**Berling, L., Klawitter, J., Bouckaert, R., Xie, D., Gavryushkina, A. & Drummond, A.J. (2025). Accurate Bayesian phylogenetic point estimation using a tree distribution parameterized by clade probabilities. *PLoS Comput Biol* 21(2): e1012789. doi:10.1371/journal.pcbi.1012789**
- *What it says (abstract).* A tractable tree distribution (a CCD) built from a posterior sample, and a point estimator that "performs at least as well and often better than standard methods", where which method is best depends on sample size and dimensionality. *(From knowledge; check before quoting:* the CCD0 and CCD1 variants and their MAP tree are proposed as replacements for MCC trees; implemented in the BEAST 2 `ccd` package.)*
- *Transfer.* If clepsydra ever samples topologies, the summary tree should be the CCD MAP tree, with node ages summarised *conditional on the clade*, over only the samples that contain it. An MCC tree should not be used.

**Heled, J. & Bouckaert, R.R. (2013). Looking for trees in the forest: summary tree from posterior samples. *BMC Evol Biol* 13: 221. doi:10.1186/1471-2148-13-221**
- *Transfer (from knowledge; check before quoting).* Shows that MCC trees with "common ancestor heights" (CAH) are better than mean or median heights over clade-containing samples. CAH assigns each summary-tree node the age of the MRCA of its taxa *in every sample*, whether or not that clade is present. This gives a topology-independent age, which is exactly the "MRCA age of a fixed taxon set" quantity (Recommendation 13).

**Holder, M.T., Sukumaran, J. & Lewis, P.O. (2008). A justification for reporting the majority-rule consensus tree in Bayesian phylogenetics. *Syst Biol* 57: 814–821. doi:10.1080/10635150802422308**; **O'Reilly, J.E. & Donoghue, P.C.J. (2018). The efficacy of consensus tree methods for summarizing phylogenetic relationships from a posterior sample of trees estimated from morphological data. *Syst Biol* 67: 354–362. doi:10.1093/sysbio/syx086**
- *Transfer.* The majority-rule consensus is the Bayes estimator under a symmetric loss, and MCC trees can over-resolve. This is background for choosing a summary.

**Dinh, V., Darling, A.E. & Matsen, F.A. IV (2018). Online Bayesian phylogenetic inference: theoretical foundations via sequential Monte Carlo. *Syst Biol* 67: 503–517. doi:10.1093/sysbio/syx087**
- *Transfer.* Adding taxa one at a time with SMC. Not a priority. It would matter only if clepsydra needed to add a language to an existing analysis without rerunning.

**Rubin, D.B. (1981). The Bayesian bootstrap. *Annals of Statistics* 9: 130–134. doi:10.1214/aos/1176345338**
- *Transfer.* Dirichlet-weighted resampling of concepts gives a cheap robustness check: which node ages depend on a few concepts? This complements `pointwise_loglik.tsv` and PSIS-LOO.

### A9. Modern samplers

**Ji, X., Zhang, Z., Holbrook, A., Nishimura, A., Baele, G., Rambaut, A., Lemey, P. & Suchard, M.A. (2020). Gradients do grow on trees: a linear-time O(N)-dimensional gradient for statistical phylogenetics. *MBE* 37: 3047–3060. doi:10.1093/molbev/msaa130**
- *Transfer.* Branch-length gradients of the pruning likelihood, computed by one postorder and one preorder traversal. They can be written for both T1 and T2 in C, but they take careful bookkeeping (pre-order partials for each concept).

**Ji, X., Fisher, A.A., Su, S., Thorne, J.L., Potter, B., Lemey, P., Baele, G. & Suchard, M.A. (2023). Scalable Bayesian divergence time estimation with ratio transformations. *Syst Biol* 72(5): 1136–1153. doi:10.1093/sysbio/syad039**
- *What it says (abstract).* Maps the N−1 internal heights to one height and N−2 ratios, gives linear-time gradients and Jacobian terms, and uses HMC. This "resolves a mixing issue" in West Nile virus and gives at least 5-fold efficiency gains.
- *Transfer.* The **ratio parameterization is valuable even without HMC**. Each internal node's age is a ratio within (lower bound, parent age), with the lower bound set by the oldest tip below it. This is essentially clepsydra's "transformed node proportions" option. Making it the default and adding a joint random walk on the root height with fixed ratios (which rescales the whole tree proportionally) is cheap. HMC on top is a larger project (Recommendation 25).

**Nishimura, A., Zhang, Z. & Suchard, M.A. (2025). Zigzag path connects two Monte Carlo samplers: Hamiltonian counterpart to a piecewise deterministic Markov process. *JASA* 120: 1077–1089. doi:10.1080/01621459.2024.2395587** (arXiv:2104.07694; DOI from publisher listing and IDEAS/RePEc, not Crossref-checked by me) **[DOI partly verified]**
- *Transfer.* Zig-zag HMC suits truncated and constrained targets such as the order constraints on node ages, but the phylogenetic applications so far are probit trait models, not node ages. It does not pay off for clepsydra.

**Bouchard-Côté, A., Sankararaman, S. & Jordan, M.I. (2012). Phylogenetic inference via sequential Monte Carlo. *Syst Biol* 61: 579–593. doi:10.1093/sysbio/syr131**; **Wang, L., Bouchard-Côté, A. & Doucet, A. (2015). Bayesian phylogenetic inference using a combinatorial sequential Monte Carlo method. *JASA* 110: 1362–1374. doi:10.1080/01621459.2015.1054487**
- *Transfer.* SMC builds ultrametric trees by successive coalescence and also returns an estimate of the marginal likelihood. It is attractive for topology-free *clock* trees, which is clepsydra's eventual case, and needs no dependencies. Its cost is many particles times the likelihood evaluations. Lower priority than MCMC moves.

**Zhang, C. & Matsen, F.A. IV (2019). Variational Bayesian phylogenetic inference. *ICLR 2019*; extended as Zhang, C. & Matsen, F.A. IV (2024). A variational approach to Bayesian phylogenetic inference. *JMLR* 25 (paper 22-0348). [no DOI; verified by web search]**
- *Transfer.* VBPI uses subsplit Bayesian networks and stochastic-gradient training. It needs automatic differentiation and is known to under-cover. It does **not** suit a goal of well-calibrated posteriors in a small dependency-free engine. Not recommended.

### A10. Hierarchical and cross-clade calibration, and calibration cross-validation

**Near, T.J. & Sanderson, M.J. (2004). Assessing the quality of molecular divergence time estimates by fossil calibrations and fossil-based model selection. *Phil. Trans. R. Soc. B* 359: 1477–1483. doi:10.1098/rstb.2004.1523**
- *What it says (abstract).* Introduces "fossil cross-validation", which ranks calibrations by their effect on overall estimation and finds those with exceptionally large error, and "fossil-based model cross-validation", which uses fossils to choose the rate-smoothing model. Both were applied to monocots (8 calibrations) and placentals (9 calibrations).
- *Transfer.* Clepsydra's LOCO already does the first. The second carries over directly: choose the clock model (strict, iid lognormal, autocorrelated) and the tree prior by **predictive accuracy on held-out calibrations**.

**Near, T.J., Meylan, P.A. & Shaffer, H.B. (2005). Assessing concordance of fossil calibration points in molecular clock studies: an example using turtles. *Am. Nat.* 165: 137–146. doi:10.1086/427734**
- *Transfer (from knowledge; check before quoting).* For each calibration in turn, date the tree from it alone and compute the deviation D_i of the other calibrated nodes from their fossil ages. Then iteratively remove the calibration with the largest SS, testing whether the variance in D changes significantly. Warnock et al. (2015) qualify this.

**Marshall, C.R. (1990). Confidence intervals on stratigraphic ranges. *Paleobiology* 16: 1–10. doi:10.1017/S0094837300009672**
- *What it says (from knowledge; check before quoting).* If H fossil horizons are distributed at random over a true range, the confidence interval on the true origin extends the observed range R by α = (1 − C)^{−1/(H−1)} − 1 (in units of R).
- *Transfer.* A *data-driven soft maximum* for a split, computed from the density of attestation. Example: the number of independent datable attestations of a branch over its observed documented span. This gives a principled upper tail where the user currently guesses a lognormal shape (Recommendation 19). Its randomness assumption is strained by the documentary record, which is a caveat to state.

**Marshall, C.R. (2008). A simple method for bracketing absolute divergence times on molecular phylogenies using multiple fossil calibration points. *Am. Nat.* 171: 726–742. doi:10.1086/587523**
- *What it says (abstract).* Uses relative branch lengths from an ultrametric tree built without fossils. Finds the lineage whose true temporal range is best covered by its fossil record, uses its oldest fossil as the minimum, and adds a confidence interval to obtain the *maximum* for all nodes. The method can also identify misdated or misassigned fossils.
- *Transfer.* Run clepsydra without calibrations to get relative ages. For each calibration, compute the implied scale. The calibration implying the oldest ages is the least-gapped minimum. This gives both a bracketing estimate and a misassignment detector that is independent of the Bayesian machinery (Recommendation 9).

**Lukoschek, V., Keogh, J.S. & Avise, J.C. (2012). Evaluating fossil calibrations for dating phylogenies in light of rates of molecular evolution: a comparison of three approaches. *Syst Biol* 61: 22–43. doi:10.1093/sysbio/syr075**
- *What it says (abstract).* Compares Near et al. (2005) cross-validation, Marshall (2008) coverage, and the Sanders & Lee (2007) Bayesian multicalibration approach, explicitly including saturation by data type.
- *Transfer.* Saturation (fast concepts here) interacts with calibration assessment: a deep calibration looks "inconsistent" under a saturating model. Calibration diagnostics should therefore be run under the best-fitting rate-heterogeneity model.

**Pyron, R.A. (2010). A likelihood method for assessing molecular divergence time estimates and the placement of fossil calibrations. *Syst Biol* 59: 185–194. doi:10.1093/sysbio/syp090**; **Rutschmann, F., Eriksson, T., Salim, K.A. & Conti, E. (2007). Assessing calibration uncertainty in molecular dating: the assignment of fossils to alternative calibration points. *Syst Biol* 56: 591–608. doi:10.1080/10635150701491156**
- *Transfer.* Uncertainty about *which node* a calibration belongs to. Linguistic analogue: does an archaeological culture date Proto-Indo-Iranian, or Proto-Indo-Aryan? Clepsydra could represent this as a **mixture over target nodes**, a discrete latent with prior weights, and report the posterior assignment.

**Dornburg, A., Beaulieu, J.M., Oliver, J.C. & Near, T.J. (2011). Integrating fossil preservation biases in the selection of calibrations for molecular divergence time estimation. *Syst Biol* 60: 519–527. doi:10.1093/sysbio/syr019**
- *Transfer.* Calibration quality depends on how good preservation (here, documentation) is around the node. Calibrations from poorly documented regions should get heavier tails.

**Heath, T.A. (2012). A hierarchical Bayesian model for calibrating estimates of species divergence times. *Syst Biol* 61: 793–809. doi:10.1093/sysbio/sys032**
- *What it says (abstract).* Puts a Dirichlet process hyperprior on the rate parameters of exponential calibration densities offset at fossil minima, which clusters them. This avoids specifying each calibration's tail by hand and is "more accurate and precise" than fixed exponentials.
- *Transfer.* For linguistic minima such as first attestations, a *shared, estimated* lag distribution is better than per-calibration guesses. A simpler version is one estimated exponential or gamma lag rate shared by all attestation-type calibrations (Recommendation 3). The DPP is optional.

**Parham, J.F. et al. (2012). Best practices for justifying fossil calibrations. *Syst Biol* 61: 346–359. doi:10.1093/sysbio/syr107**
- *Transfer.* A five-point justification checklist (specimen, phylogenetic placement, age, and so on). Clepsydra's calibration TSV has a `source_key`, and a parallel "justification" schema with fields for identification, chronology and node-assignment confidence could feed the prior outlier probabilities.

**Ho, S.Y.W. & Phillips, M.J. (2009). Accounting for calibration uncertainty in phylogenetic estimation of evolutionary divergence times. *Syst Biol* 58: 367–380. doi:10.1093/sysbio/syp035**
- *Transfer.* A taxonomy of calibration densities (uniform, normal, lognormal, exponential, gamma) and what each encodes. It is already reflected in `docs/calibrations.md`, and the bibliography should cite it.

**Ho, S.Y.W. et al. (2008). The effect of inappropriate calibration: three case studies in molecular ecology. *PLoS ONE* 3: e1615. doi:10.1371/journal.pone.0001615**; **Sauquet, H. et al. (2012). Testing the impact of calibration on molecular divergence times using a fossil-rich group: the case of *Nothofagus*. *Syst Biol* 61: 289–313. doi:10.1093/sysbio/syr116**; **Duchêne, S., Lanfear, R. & Ho, S.Y.W. (2014). The impact of calibration and clock-model choice on molecular estimates of divergence times. *Mol. Phylogenet. Evol.* 78: 277–289. doi:10.1016/j.ympev.2014.05.032**
- *Transfer (from knowledge; check before quoting).* (a) Borrowing rates from deep-calibrated clades for shallow questions, or the reverse, is biased because of time dependency. (b) Several well-spread calibrations beat a single one, and the *position* of calibrations (root against internal, deep against shallow) matters more than their number. Concretely, clepsydra should report which calibrations are deep or shallow relative to the queried node.

### A11. Validation methodology

**Cook, S.R., Gelman, A. & Rubin, D.B. (2006). Validation of software for Bayesian models using posterior quantiles. *JCGS* 15: 675–692. doi:10.1198/106186006X136976**; **Talts, S., Betancourt, M., Simpson, D., Vehtari, A. & Gelman, A. (2018/2020). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788 [no DOI; verified by web search]**
- *Transfer.* Clepsydra already uses SBC. Keep the rank-based version with thinning to roughly independent draws.

**Modrák, M., Moon, A.H., Kim, S., Bürkner, P., Huurre, N., Faltejsková, K., Gelman, A. & Vehtari, A. (2025). Simulation-based calibration checking for Bayesian computation: the choice of test quantities shapes sensitivity. *Bayesian Analysis* 20(2): 461–488. doi:10.1214/23-BA1404**
- *What it says (abstract).* Standard SBC "could never detect large classes of problems including when the posterior is equal to the prior". Adding data-dependent test quantities, especially the **joint log-likelihood of the data**, makes SBC sensitive in principle to any problem. It also clarifies the difference between SBC and checks based on the data-averaged posterior.
- *Transfer.* Add log L(data | θ) as an SBC test quantity, along with node-age *ratios* and tree height (Recommendation 7). Otherwise a sampler that ignores the cognate data entirely (for example, because of a bug that zeroes the likelihood gradient in an update) would still pass SBC on node ages.

**Säilynoja, T., Bürkner, P.-C. & Vehtari, A. (2022). Graphical test for discrete uniformity and its applications in goodness-of-fit evaluation and multiple sample comparison. *Statistics and Computing* 32: 32. doi:10.1007/s11222-022-10090-6**
- *Transfer.* Simultaneous confidence bands for ECDF-difference plots of SBC ranks. They replace histogram eyeballing and multiple per-parameter χ² tests.

**Mendes, F.K., Bouckaert, R., Carvalho, L.M. & Drummond, A.J. (2025). How to validate a Bayesian evolutionary model. *Syst Biol* 74(1): 158–175. doi:10.1093/sysbio/syae064**
- *What it says (abstract and knowledge; check details).* Good practice for checking that a model *implementation* is correct, centred on coverage of 95% HPD intervals across many simulated replicates, plus direct simulation against MCMC prior sampling and similar checks, with automated tooling (LPhy and BEAST 2).
- *Transfer.* Use it as the checklist for clepsydra's validation doc. In particular, report **coverage per parameter class** (root age, internal node ages, `σ_b`, the clock scale), and compare prior samples from the MCMC against direct prior simulation. That second check is the effective-prior test.

**Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: an improved R̂ for assessing convergence of MCMC. *Bayesian Analysis* 16: 667–718. doi:10.1214/20-BA1221**
- *Transfer.* Rank-normalized split-R̂, bulk ESS and tail ESS. Tail ESS matters most here because people quote 95% HPDs of node ages. Implementing these takes a sort and some arithmetic in C.

**Huelsenbeck, J.P. & Rannala, B. (2004). Frequentist properties of Bayesian posterior probabilities of phylogenetic trees under simple and complex substitution models. *Syst Biol* 53: 904–913. doi:10.1080/10635150490522629**
- *Transfer.* The template for "coverage under misspecification": simulate under a richer model than the one fitted, then measure coverage. Under-modelled data gives overconfident posteriors.

**Bollback, J.P. (2002). Bayesian model adequacy and choice in phylogenetics. *MBE* 19: 1171–1180. doi:10.1093/oxfordjournals.molbev.a004175**; **Brown, J.M. (2014). Detection of implausible phylogenetic inferences using posterior predictive assessment of model fit. *Syst Biol* 63: 334–348. doi:10.1093/sysbio/syu002**; **Höhna, S., Coghill, L.M., Mount, G.G., Thomson, R.C. & Brown, J.M. (2018). P3: Phylogenetic posterior prediction in RevBayes. *MBE* 35: 1028–1034. doi:10.1093/molbev/msx286**
- *Transfer.* Data-based statistics (Bollback's multinomial and pattern-frequency statistics) and *inference-based* statistics (Brown: refit trees to posterior-predictive datasets and compare tree length and similar quantities). For cognates, the data-based statistics would be counts of singleton classes per language, the distribution of cognate-class sizes, and the frequency of patterns shared between pairs.

**Duchêne, D.A., Duchêne, S., Holmes, E.C. & Ho, S.Y.W. (2015). Evaluating the adequacy of molecular clock models using posterior predictive simulations. *MBE* 32: 2986–2995. doi:10.1093/molbev/msv154**
- *What it says (abstract).* A posterior-predictive test of clock adequacy that compares per-branch lengths estimated from simulated data with those from the observed data. It is sensitive to the branch-length bias caused by under-parameterized clocks. The multinomial statistic had low power for clock adequacy.
- *Transfer.* Directly implementable (Recommendation 15). For each branch, compare the posterior-predictive distribution of lexical branch length with the length estimated from the real data under an unconstrained (non-clock) fit. The proportion of branches with extreme predictive p-values is the adequacy index.

**Xie, W., Lewis, P.O., Fan, Y., Kuo, L. & Chen, M.-H. (2011). Improving marginal likelihood estimation for Bayesian phylogenetic model selection. *Syst Biol* 60: 150–160. doi:10.1093/sysbio/syq085**; **Baele, G., Lemey, P., Bedford, T., Rambaut, A., Suchard, M.A. & Alekseyenko, A.V. (2012). Improving the accuracy of demographic and molecular clock model comparison while accommodating phylogenetic uncertainty. *MBE* 29: 2157–2167. doi:10.1093/molbev/mss084**
- *Transfer.* Stepping-stone and path sampling (power posteriors). Clepsydra's tempered chains (Recommendation 16) can give these almost for free. Caution: proper priors are required, and the current flat node-age prior makes marginal likelihoods undefined.

### A12. Context from linguistics (cited only for the transfer arguments)

**Pagel, M., Atkinson, Q.D. & Meade, A. (2007). Frequency of word-use predicts rates of lexical evolution throughout Indo-European history. *Nature* 449: 717–720. doi:10.1038/nature06176** (rates differ by orders of magnitude across meanings, which supports the PoW-type heterogeneity argument). **Greenhill, S.J. et al. (2017). Evolutionary dynamics of language systems. *PNAS* 114: E8822–E8829. doi:10.1073/pnas.1700388114** (rates of lexical and structural change vary across families, which bears on cross-family rate sharing). **Gray, R.D. & Atkinson, Q.D. (2003). *Nature* 426: 435–439. doi:10.1038/nature02029**; **Bouckaert, R. et al. (2012). *Science* 337: 957–960. doi:10.1126/science.1219669**.

---

## Part B — Lessons that cut across the sources

1. **Clepsydra's current flat prior is improper, and it is not uninformative.** Rate and time are confounded (Zhu et al. 2015), and node ages have only order constraints with no upper bound. A family with no calibration and no proper clock-scale prior therefore has an improper posterior: ages go to ∞ as the rate goes to 0. Even with a bounded root, a uniform density over the order-constrained region implies a root-age marginal proportional to the *volume* of the region below it. For contemporaneous tips this volume grows polynomially in root age, with degree equal to the number of other uncalibrated internal nodes. The prior therefore leans toward old roots more strongly the larger the tree. MrBayes's uniform prior, OxCal's boundary construction and MCMCTree's birth–death kernel all exist to avoid exactly this.
2. **Calibration densities are not marginal priors.** Every construction distorts them (Heled & Drummond 2012, 2015; Inoue et al. 2010; Barba-Montoya et al. 2017; Rannala 2016). The only honest practice is to sample the joint prior and show its marginals next to the specified ones.
3. **Cognate data cannot beat the calibrations.** In the infinite-data limit, node-age uncertainty is set by the calibrations and, under a relaxed clock, by the single realized pattern of branch-rate noise (Rannala & Yang 2007; Zhu et al. 2015). Concepts are sites, not loci, so more concepts do not average away `m_b`. For uncalibrated families, the limit is the spread of the *between-family* rate predictive. If that spread has a coefficient of variation of 30%, no cognate list can date such a family more precisely than about ±30% (1σ). The output should say so.
4. **Calibrations are fallible measurements and should be modelled that way.** Archaeology learned to fit outliers rather than delete them (Bronk Ramsey 2009b; Christen 1994; Lanos & Philippe 2017). Phylogenetics learned that a posteriori deletion can pick the wrong calibrations (Warnock et al. 2015). Mixture or outlier models with reported posterior outlier probabilities satisfy both lessons.
5. **Rates depend on timescale whenever rate heterogeneity is under-modelled** (Ho et al. 2005, 2011; Aiewsakun & Katzourakis 2016; Ghafari et al. 2021). Language calibrations are mostly shallow. Taking shallow-calibrated rates to deep, uncalibrated families is the most likely source of a large systematic error, and it calls for a test of rate against calibration depth.
6. **Calibration has to be validated predictively.** SBC checks the code. Coverage under *misspecified* simulation checks robustness. Leave-one-calibration-out and leave-one-family-out predictive coverage checks whether the claimed uncertainty is honest on real data. These are three different questions, and all three need answering.

---

## Part C — Ranked recommendations

Ranking is by expected improvement in the *calibration* of node-age posteriors per unit of implementation effort. Difficulty assumes clepsydra's current dependency-free C99 code base: **S** = small (≤ 1–2 days, local change), **M** = medium (≤ 1–2 weeks, new module and tests), **L** = large (multi-week, new subsystem).

### 1. Replace the flat node-age prior with a proper tree prior, and require a proper root or clock-scale prior — **M**
- **What.** Add a `tree_prior` option. (a) `uniform_root`: node ages uniform given the root age, with an explicit root prior. Without a root calibration, use a broad but proper default, for example a lognormal tied to the clock-scale predictive, and say so in the manifest. (b) `birth_death`: the Yang & Rannala (1997, 2006) kernel with (λ, μ, ρ), in which uncalibrated ages are order statistics of iid kernel draws conditioned on the root and calibrated nodes. Give λ−μ a hyperprior and ρ a user value or one estimated from unseen-lineage counts (Kestemont et al. 2022). (c) Keep `flat` only for backwards-compatibility tests, and warn that it is improper.
- **Why.** Fixes the improper posterior for uncalibrated families and the hidden prior pressure toward old roots. It also makes marginal likelihoods (Recommendation 16) defined.
- **Sources.** Yang & Rannala 1997, 2006; Rannala & Yang 2007; Gernhard 2008; Stadler 2009; Höhna et al. 2011; Ronquist et al. 2012; Rannala et al. 2012; Bronk Ramsey 2009a (the boundary analogy).
- **C99 notes.** On a fixed topology there are no ranked-topology terms. The conditional-construction density is a product over uncalibrated nodes of the kernel density, renormalized within each calibrated interval, so it is O(n) per evaluation. Update it incrementally on single-node moves.

### 2. Make the effective prior a first-class output — **S**
- **What.** Run prior-only automatically, or on `--effective-prior`. For every calibrated node, report the specified density against the effective marginal: quantiles, the total-variation or KS distance, and a "distortion" flag. For every uncalibrated node, report the prior 95% interval next to the posterior. Also report the prior/posterior overlap (posterior interval width ÷ prior interval width) as a measure of how much the data informed each node.
- **Why.** Calibration densities are never the marginals (Heled & Drummond 2012). Effective priors shift with other calibrations (Warnock et al. 2015). Truncation by node order can look like data signal (Rannala 2016).
- **Sources.** Heled & Drummond 2012, 2015; Inoue et al. 2010; Barba-Montoya et al. 2017; Warnock et al. 2012, 2015; Rannala 2016; Mendes et al. 2025.
- **C99 notes.** Prior-only sampling exists already. This adds a summary table and validates the MCMC prior sample against *direct* simulation where possible: for the birth–death prior, simulate directly and compare.

### 3. Robust calibrations: outlier mixtures, t-densities and one-sided attestation lags — **S–M**
- **What.**
  - (a) A `p_outlier` column per calibration (default 0, recommended 0.05–0.1 for contested identifications). The density becomes `(1−p)·f(t) + p·g(t)`, where `g` is either a broad t-shifted copy of `f` (the OxCal General model) or a vague density over the feasible interval. Marginalise analytically and compute the posterior outlier probability per draw as `p·g/((1−p)f + p·g)`; no latent indicator is needed.
  - (b) A `t_node(mean, scale, ν)` type.
  - (c) An `attestation_lag` type: age = attestation date + Lag, with Lag ~ Exponential(κ) or Gamma, and κ *shared and estimated* across all attestation-type calibrations, as in Heath's (2012) hierarchical approach without the DPP.
  - (d) A *multi-source* calibration, in which several independent dates for one node are combined hierarchically with individual variances (the ChronoModel event model) instead of being multiplied.
  - (e) Optionally, an uncertain target node: a mixture over two or three candidate nodes (Rutschmann et al. 2007; Pyron 2010).
- **Why.** Linguistic calibrations (archaeological correlates, historical identifications) are often wrong in ways a lognormal cannot express. Deleting them after the fact is unreliable. Multiplying conflicting sources gives false precision.
- **Sources.** Bronk Ramsey 2009b; Christen 1994; Christen & Pérez 2009; Lanos & Philippe 2017, 2018; Heath 2012; Warnock et al. 2015; Rutschmann et al. 2007; Dornburg et al. 2011.
- **C99 notes.** (a)–(c) are new density functions plus a column. (d) adds a latent per source. (e) needs a discrete latent with a Gibbs step. Report the outlier posteriors in `posterior_calibrations.tsv`.

### 4. Hierarchical cross-family clock-scale prior, with an honest predictive for uncalibrated families — **M**
- **What.**
  - Model the family log clock scale as `β_f ~ Normal(μ_β, τ²)` or Student-t, with hyperpriors on `μ_β` and `τ`.
  - Fit it either jointly over several calibrated families in one run (a multi-tree likelihood, a sum over families) or in two stages (empirical Bayes: fit calibrated families, export the posterior predictive for a new family's `β`, as a t-mixture or as samples).
  - For an uncalibrated family, use the *predictive*, not a plug-in.
  - Fix the concept-rate parameterization so that the prior variance of the mean lexical rate does not shrink as 1/C. Either centre the `μ_c` (sum-to-zero on the log scale) with all mean-rate uncertainty in the global scale, or use a gamma-Dirichlet split (dos Reis et al. 2014).
  - Add an optional time-depth covariate (Recommendation 11) so that deep, uncalibrated families do not inherit shallow-calibrated rates unadjusted.
- **Why.** Removes absolute-time non-identifiability for uncalibrated families with honest, wide uncertainty. Prevents the iid-rate-prior collapse that gives "wrong values with very narrow credibility intervals". Avoids the false precision of secondary calibrations.
- **Sources.** dos Reis et al. 2014; Heath et al. 2012 (a DPP over families as an extension); Greenhill et al. 2017 (variation across families); Schenk 2016; Graur & Martin 2004; Ho et al. 2008.
- **C99 notes.** Joint multi-family MCMC means a vector of tree states sharing hyperparameters. The Gibbs update for `μ_β, τ` given the `β_f` is conjugate or slice-sampled. The two-stage version is easier: export the predictive as a TSV reference prior, which fits the existing `reference/` package design.

### 5. Leave-one-family-out (LOFO) and leave-one-calibration-out predictive coverage as the headline validation — **S**
- **What.** For each calibrated family f, remove all its calibrations, date it using only the cross-family predictive from Recommendation 4 (fitted *without* f), and record the PIT value of the removed calibration ages under the resulting posterior. Do the same within families with LOCO. Report empirical coverage of the 50%, 80% and 95% intervals and a PIT histogram or ECDF band across all held-out calibrations.
- **Why.** This is the only test that measures whether *the uncalibrated-family pipeline* is calibrated on real data. SBC cannot do it, because SBC assumes the model is true. It is fossil-based model cross-validation applied to the whole method.
- **Sources.** Near & Sanderson 2004; Near et al. 2005; Säilynoja et al. 2022; Huelsenbeck & Rannala 2004.
- **C99 notes.** Orchestration already exists (LOCO). The new piece is the PIT bookkeeping. Treat PIT values from overlapping families as dependent and say so.

### 6. Rate–time ridge operators (constant-distance moves) — **S**
- **What.**
  - **Node-age with constant distance** (Zhang & Drummond 2020): propose a new age for internal node v and rescale the multipliers of v's parent branch and its two child branches so that each `m_b·t_b` is unchanged. The likelihood is unchanged, so the acceptance ratio is the prior ratio (ages, `m_b`, tree prior) times the Jacobian, which is the product of the `t_old/t_new` ratios.
  - **Global scale move**: multiply all ages by c and the global rate by 1/c.
  - At the root, use the "simple distance" variant.
- **Why.** Branch multipliers and ages are confounded (clepsydra's docs §7.2 already flag this). Single-site MH crawls along the ridge, which deflates tail ESS and makes HPDs untrustworthy. Douglas et al. (2021) report up to 65× efficiency gains from operators of this kind.
- **Sources.** Zhang & Drummond 2020; Douglas et al. 2021.
- **C99 notes.** About 100 lines. The likelihood cache stays valid because path lengths are unchanged, which makes this the cheapest high-payoff item on the list.

### 7. SBC and coverage upgrades — **S**
- **What.**
  - Add the joint log-likelihood, tree height, node-age ratios and `σ_b` as SBC test quantities (Modrák et al. 2025).
  - Use ECDF-difference plots with simultaneous bands (Säilynoja et al. 2022).
  - Add a *misspecification* coverage suite: simulate with borrowing, an autocorrelated clock, time-dependent rates and wrong calibration centres, fit the standard model, and report 95% HPD coverage per parameter class.
  - Report rank-normalized R̂ and tail ESS (Vehtari et al. 2021).
- **Why.** Standard SBC cannot detect a posterior that equals the prior (Modrák et al. 2025). Calibration under misspecification is what users actually experience.
- **Sources.** Cook et al. 2006; Talts et al. 2018; Modrák et al. 2025; Säilynoja et al. 2022; Mendes et al. 2025; Huelsenbeck & Rannala 2004; Vehtari et al. 2021.

### 8. Precision-limit (infinite-sites) report — **S**
- **What.** Output the infinite-sites plot (posterior 95% width against posterior mean for all nodes) with the regression slope and R². Run an "amplified data" diagnostic, refitting with every concept replicated k = 2, 4, 8 times, and report how much each node's interval shrinks. Report the decomposition: if widths barely shrink, the node is calibration-limited or rate-prior-limited. For uncalibrated families, report the theoretical floor, the CV of the clock-scale predictive.
- **Why.** Tells users what a cognate list *can* tell them. It prevents demands for "more concepts" when calibrations are the bottleneck, and makes the "concepts are sites, not loci" point concrete.
- **Sources.** Rannala & Yang 2007; dos Reis & Yang 2013; Zhu et al. 2015; dos Reis et al. 2012.
- **C99 notes.** Replication is a data-weight multiplier on the per-concept log-likelihood (a power posterior with exponent k), which is trivial.

### 9. Calibration-consistency diagnostics that do not delete anything — **S**
- **What.**
  - (a) A calibration-free relative chronology: run with the root fixed to 1 and no calibrations (RelTime-like).
  - (b) For each calibration, the implied absolute scale and its interval, shown in one forest plot. Calibrations whose scale intervals do not overlap are in conflict.
  - (c) The Near et al. (2005) SS statistic and the Marshall (2008) "most complete lineage" bracket as tables.
  - (d) A "narrow-but-conflicted" flag: nodes with small posterior width *and* high calibration tension (Rannala 2016).
  - The output text must state that these are diagnostics, not filters (Warnock et al. 2015), and point to the outlier model in Recommendation 3.
- **Sources.** Near & Sanderson 2004; Near et al. 2005; Marshall 2008; Lukoschek et al. 2012; Battistuzzi et al. 2015; Rannala 2016; Warnock et al. 2015.

### 10. Tip-age uncertainty for ancient languages — **S**
- **What.** Add `uniform_tip`, `normal_tip` and `attested_span_tip` (uniform over the documented span), with tip ages sampled by MH inside their bounds.
- **Why.** Fixing tip ages biases node ages and topology (Barido-Sottani et al. 2020). Ancient languages are attested over centuries.
- **Sources.** Barido-Sottani et al. 2020; Stadler & Yang 2013; Ronquist et al. 2012.

### 11. Time-dependence and saturation diagnostics, then an optional rate-heterogeneity fix — **M**
- **What.**
  - (a) **Diagnostic.** For each calibration (and LOCO fit), compute the implied family or global rate and regress log rate on log calibration age, across families too. A significantly negative slope (power-law decay) is a red flag for deep dating (Aiewsakun & Katzourakis 2016).
  - (b) **Model check.** Posterior-predictive comparison of the proportion of non-cognate pairs against pair divergence time across all calibrated pairs. A saturating mismatch at depth suggests the concept-rate distribution has too thin a slow tail (the PoW mechanism).
  - (c) **Remedy.** Make the concept-rate distribution flexible: a mixture of two lognormals, or a discretized gamma with a "very slow" component. Optionally add an epoch clock (Membrebe et al. 2019), identifiable only with calibrations at several depths.
- **Why.** Most linguistic calibrations are shallow and many targets are deep. Rate decay with depth causes systematic underestimation of deep ages, as in viruses (Ghafari et al. 2021: up to 30× too young).
- **Sources.** Ho et al. 2005, 2011, 2015; Aiewsakun & Katzourakis 2016; Membrebe et al. 2019; Ghafari et al. 2021; Pagel et al. 2007; Lukoschek et al. 2012.
- **C99 notes.** (a) and (b) are reporting. (c) changes the concept-rate prior; the likelihood code already takes per-concept rates.

### 12. Adaptive multivariate and Bactrian proposals during warm-up — **S–M**
- **What.** During warm-up, estimate the covariance of (log node-age ratios, log `σ_b`, log global rate, the `μ_c` block). Then use an AVMVN proposal (Haario-type adaptive Metropolis, frozen after warm-up as the docs require) for the continuous block, plus Bactrian random-walk kernels, which have bimodal step distributions and better acceptance-versus-step trade-offs. Add adaptive operator weights.
- **Why.** Node ages are highly correlated, so single-site updates mix slowly.
- **Sources.** Douglas et al. 2021.
- **C99 notes.** A Cholesky factor of a dense (n+k)² matrix. Fine for trees with fewer than about 200 nodes.

### 13. Topology-independent targets: MRCA of taxon sets, clade-conditional ages and node-ordering probabilities — **S**
- **What.** Let calibrations and queries name *taxon sets* (an MRCA) instead of labelled nodes, with an optional `monophyly` flag. Report node ages as CAH-style MRCA ages of named sets, plus joint probabilities such as P(age A > age B). This is a prerequisite for Recommendations 14 and 22.
- **Why.** Once topology varies, a "named node" may not exist in every sample. The MRCA age of a fixed taxon set always does (Heled & Bouckaert 2013). Joint ordering questions need the joint samples (Parnell et al. 2008).
- **Sources.** Heled & Drummond 2012 (calibrations with and without monophyly); Heled & Bouckaert 2013; Parnell et al. 2008.

### 14. Dating over a supplied posterior tree set (integrating topological uncertainty in two steps) — **M**
- **What.** Accept a Newick or NEXUS *set* of topologies, such as the posterior from a BEAST or MrBayes topology run. Run clepsydra's fixed-tree sampler on each of N sampled topologies (N = 100–1000, parallel and seeded), then pool the samples of MRCA ages of named taxon sets. Report the between-topology share of the variance.
- **Why.** Clepsydra gets most of the benefit of topology uncertainty without writing tree-space MCMC. Report node ages conditional on a clade only over the topologies that contain it.
- **Caveat.** This conditions each topology on the same data twice: once in the topology posterior and once in dating. The bias is usually small but should be documented. It is an approximation, not the joint posterior.
- **Sources.** Larget 2013; Berling et al. 2025; Heled & Bouckaert 2013.

### 15. Posterior-predictive clock adequacy — **M**
- **What.** Following Duchêne et al. (2015): per branch, compare the posterior-predictive lexical branch length (re-estimated on simulated data under a free, non-clock model) with the length from the observed data. Report the share of branches outside the 95% predictive interval, and also cognate-specific statistics (singleton counts per language, the class-size distribution).
- **Why.** Shows whether the clock model can generate the data. Relative model comparison cannot.
- **Sources.** Bollback 2002; Brown 2014; Höhna et al. 2018; Duchêne et al. 2015.
- **C99 notes.** The simulator exists for SBC. The free branch-length fit can use the penalized-likelihood or LSD routine from Recommendation 20.

### 16. Metropolis-coupled MCMC and stepping-stone marginal likelihoods — **M**
- **What.** Heated chains `L^β·prior` with swaps, which reuse the 4-chain infrastructure. Use a power-posterior ladder for stepping-stone marginal likelihoods to compare tree priors and clocks, valid only with proper priors (Recommendation 1).
- **Why.** Multimodality from conflicting calibrations or outlier mixtures, and principled model comparison alongside PSIS-LOO.
- **Sources.** Altekar et al. 2004; Xie et al. 2011; Baele et al. 2012.

### 17. Sampled ancestors and a fossilized-birth–death-style prior for dated ancient languages — **M–L**
- **What.** Treat dated ancient tips as serial samples under a BDS or FBD prior (Stadler 2010; Heath et al. 2014). Give each ancient tip an "is ancestor" indicator: a reversible jump between a zero-length pendant edge (the tip lies on the branch below its attachment point) and a positive-length one. On a fixed topology the ancestor state is defined only where the attachment node is binary-compatible. Report P(ancestor) per ancient language.
- **Why.** Treating Latin as a sister of Romance forces older Romance roots. Sampled ancestors resolve this coherently and make ancient tips calibrations through the tree process rather than through ad hoc densities.
- **Sources.** Gavryushkina et al. 2014, 2017; Heath et al. 2014; Zhang et al. 2016; Stadler 2010; Chang et al. 2015; Heggarty et al. 2023.

### 18. Ratio-transform node parameterization as the default — **S–M**
- **What.** Parameterize the ages as root height × ratios, with each internal node's ratio in (0,1) between its lower bound (the oldest descendant tip or calibration minimum) and its parent's age. Include the Jacobian. Use random walks on logit(ratio), and a root move with the ratios held fixed, which rescales the whole subtree proportionally.
- **Why.** Order constraints become box constraints and correlated heights are decorrelated, which improves mixing (Ji et al. 2023). It is also the prerequisite for any later gradient sampler.
- **Sources.** Ji et al. 2023.
- **Note.** Clepsydra's "transformed node proportions" may already be close to this, so check it before reimplementing.

### 19. Attestation-density-based soft maxima (Marshall CI tool) — **S**
- **What.** A Python helper under `tools/` that takes a list of dated attestations or horizons for a lineage and outputs a soft-maximum calibration using Marshall's (1990) range-extension formula, with a caveat flag when H is small or sampling is clearly non-random.
- **Why.** Replaces guessed lognormal tail shapes with a reproducible, data-derived bound.
- **Sources.** Marshall 1990, 2008; Dornburg et al. 2011.

### 20. Fast non-Bayesian dating as cross-check and initializer (LSD/penalized likelihood) — **S–M**
- **What.** A least-squares dating routine (LSD-style, with variance ∝ branch length + c) on posterior-mean lexical branch lengths from an unconstrained fit, with node constraints from the calibrations. Optionally add a treePL-style roughness penalty with cross-validated smoothing. Use it to initialize the MCMC and print a "fast estimate" beside the Bayesian one.
- **Why.** Gives a quick sanity check: large disagreements point to prior dominance. It speeds burn-in and supports the adequacy test in Recommendation 15.
- **Sources.** To et al. 2016; Sanderson 2002; Smith & O'Meara 2012; Paradis 2013; Sagulenko et al. 2018.

### 21. Clock-model breadth: autocorrelated, mixed and shrinkage-local clocks, plus a strict-clock spike — **M**
- **What.**
  - A geometric-Brownian autocorrelated option (Thorne & Kishino) and the Lartillot et al. (2016) mixed clock (Brownian plus white noise).
  - Later, the Fisher et al. (2023) shrinkage local clock, sampled by MH on increments.
  - A spike at `σ_b = 0` to give P(strict clock).
- **Why.** Lexical conservatism is plausibly inherited, and the iid clock can mis-attribute it to time. Tip-dating results are especially sensitive to the clock (Lartillot et al. 2016).
- **Sources.** Thorne et al. 1998; Kishino et al. 2001; Lepage et al. 2007; Lartillot et al. 2016; Drummond & Suchard 2010; Fisher et al. 2023; Li & Drummond 2012.

### 22. Full topology MCMC with clock-tree moves, and CCD summaries — **L**
- **What.** Narrow exchange, subtree-slide and Wilson–Balding on time trees, compatible with tip dates, monophyly constraints and calibrations on taxon sets (Recommendation 13). Summarise with a CCD MAP tree and clade-conditional ages, and diagnose with ASDSF and tree-space checks.
- **Why.** Gives the joint posterior of tree and dates, the "without a given tree" use case.
- **Sources.** Wilson & Balding 1998; Drummond et al. 2002; Lakner et al. 2008; Höhna & Drummond 2012; Larget 2013; Berling et al. 2025; Whidden & Matsen 2015.
- **C99 notes.** The hard parts are incremental likelihood recomputation after pruning and regrafting (partials invalidated along two paths), the Hastings ratios, and keeping calibration and feasibility propagation dynamic.

### 23. Bayesian-bootstrap concept-influence report — **S**
- **What.** Reweight the per-concept log-likelihood with Dirichlet(1,…,1) weights, using a short refit or importance reweighting of existing draws via `pointwise_loglik`. Report how much each node age moves, and which concepts drive it.
- **Why.** Detects dates driven by a handful of concepts, for example a few loans that were missed.
- **Sources.** Rubin 1981. It complements PSIS-LOO.

### 24. Calibration provenance and justification schema — **S**
- **What.** Extend the calibration TSV with Parham-style justification fields: identification confidence, chronology source, node-assignment confidence. Map them to default `p_outlier` and tail weights.
- **Why.** Makes calibration quality explicit and feeds Recommendation 3.
- **Sources.** Parham et al. 2012; Dornburg et al. 2011; Bayliss 2009.

### 25. HMC with ratio transforms: only for large trees — **L**
- **What.** Linear-time gradients of the T1 and T2 log-likelihoods with respect to the ratio-transformed heights, then HMC or NUTS.
- **Why.** Pays off only above roughly 100–200 tips (Ji et al. 2020, 2023). Most lexical datasets are 20–150 languages, where Recommendations 6, 12 and 18 recover most of the mixing gain.
- **Verdict.** Defer. Variational inference (VBPI) and zig-zag samplers are not recommended for a calibration-first, dependency-free engine. SMC (Bouchard-Côté et al. 2012; Wang et al. 2015) is worth reconsidering only if Recommendation 22 proves to mix badly.

### 26. Grid / message-passing exact marginals as a test oracle for small trees — **M**
- **What.** For strict-clock fixed trees with 4–8 tips, compute node-age marginals by discretized message passing, in the style of TreeTime, using the existing `gd_grid`. Compare them with the MCMC.
- **Why.** An independent, near-exact check of the sampler and the tree prior that does not depend on MCMC.
- **Sources.** Sagulenko et al. 2018.

### Suggested order of work
Phase 1 (priors and validation, about 3–4 weeks): Recommendations 1, 2, 6, 3, 7, 8, 10. Phase 2 (calibrations across families): 4, 5, 9, 11, 19, 24. Phase 3 (mixing and adequacy): 12, 18, 15, 16, 20, 23, 26. Phase 4 (beyond a fixed tree): 13, 14, 17, 21, 22. Defer 25.

---

## Part D — Verification notes and unverified items

- Every DOI in Part A was resolved on Crossref on 2026-10-10 except: Nishimura et al. 2025 (DOI from the publisher listing and RePEc, so treat it as verified but not Crossref-checked); Talts et al. (arXiv:1804.06788, no DOI); Zhang & Matsen (ICLR 2019, JMLR 25 paper 22-0348, no DOI); and Lanos & Philippe 2017 (J-SFdS, Numdam, no DOI).
- Content was checked against abstracts for: Heled & Drummond 2012, 2015; Rannala 2016; Rannala & Yang 2007; Zhu et al. 2015; dos Reis et al. 2014; Inoue et al. 2010; Barba-Montoya et al. 2017; Warnock et al. 2015; Near & Sanderson 2004; Marshall 2008; Heath 2012; Lukoschek et al. 2012; Battistuzzi et al. 2015; Schenk 2016; Ghafari et al. 2021; Membrebe et al. 2019; Douglas et al. 2021; Zhang & Drummond 2020; Lartillot et al. 2016; Fisher et al. 2023; Ji et al. 2023; Berling et al. 2025; Modrák et al. 2025; Mendes et al. 2025; Duchêne et al. 2015; da Silva & Tehrani 2016.
- Statements marked *(from knowledge; check before quoting)* were not checked against the papers in this session. That includes specific formulas: Marshall's 1990 α, the Yang–Rannala kernel, the details of OxCal's outlier models, and the exact LSD weighting.
- The claim in Part B §1 about root-age pressure under a flat prior is a derivation (the volume of the order-constrained region), not a citation. Confirm it on clepsydra itself by running prior-only on an uncalibrated tree with a bounded root.
- Two findings that do not support the user's framing: (i) the Ghafari et al. PoW paper was published in *Current Biology* (2021), not *Nature Communications*; (ii) the Modrák SBC paper's published title is "Simulation-based calibration checking for Bayesian computation: the choice of test quantities shapes sensitivity" (*Bayesian Analysis* 2025), not "SBC can fail to detect…".
- Not located in this session, so not cited: a Bayesian stemmatology study that uses tip-dated birth–death priors on manuscripts (McCollum & Turnbull 2024 is the closest). Also Sanders & Lee 2007 (*Biol. Lett.* 3: 275–279), which is cited by Lukoschek et al., but whose DOI I did not verify.

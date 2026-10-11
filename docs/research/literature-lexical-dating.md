# Lexical and linguistic dating: literature survey for clepsydra

Compiled 2026-10-10. Scope: lexical (basic-vocabulary) dating of language-family nodes, from
classical glottochronology to Bayesian phylolinguistics, with emphasis on what clepsydra
(pairwise exponential retention with concept-specific rates; strict novel-state multistate tree
likelihood after Rönchen, Wiklund & Hammarström 2024; binary gain/loss CTMC; fixed topology;
calibrations; relaxed clock; contamination mixture) can adopt, test against, or expose.

**Verification conventions.** Every DOI in this file was resolved against Crossref (or the
publisher page) during compilation unless marked **[DOI unverified]**. Numbers quoted from
abstracts or publisher pages were read during compilation; numbers marked **[unverified]**
are from memory or secondary sources and must be checked before being used as benchmarks.
"BP" = before present (conventionally 1950, though most linguistic papers use "years ago"
loosely); "BCE" is kept where a paper uses it.

---

## Part 1. Annotated bibliography

### 1.1 Classical lexicostatistics and glottochronology

**Swadesh, M. (1950).** Salish internal relationships. *International Journal of American
Linguistics* 16(4): 157–167. https://doi.org/10.1086/464084
First explicit use of shared-cognate percentages on a fixed meaning list to subgroup and date a
family. *Relevance:* historical origin of the pairwise model (P0) that clepsydra still ships as its
"classical" baseline.

**Swadesh, M. (1952).** Lexico-statistic dating of prehistoric ethnic contacts. *Proceedings of the
American Philosophical Society* 96(4): 452–463. (JSTOR stable 3143802; no Crossref DOI found.)
The 200-item list and the formula t = ln C / (2 ln r). *Relevance:* the constant-rate,
equal-rate-across-meanings, independent-lineages assumptions that every later model relaxes; a
useful didactic "decay law" to keep for comparison in the web page.

**Swadesh, M. (1955).** Towards greater accuracy in lexicostatistic dating. *IJAL* 21(2): 121–137.
https://doi.org/10.1086/464321
Reduced the list to 100 items, chosen for stability and universality, with retention r ≈ 0.86 per
millennium (commonly cited figure; **[value unverified against the paper's text]**). *Relevance:*
the 100-item list is still the core of most datasets clepsydra will ingest; meaning-list choice
interacts with rate priors (a prior calibrated on a 200-item list is not transferable to a
100-item or 40-item list without per-concept rates, which clepsydra already has).

**Lees, R. B. (1953).** The basis of glottochronology. *Language* 29(2): 113–127.
https://doi.org/10.2307/410164
Estimated retention from 13 language pairs with documented histories, giving the canonical
r = 0.805 ± 0.0176 per millennium for the 215-item list (**[figures as commonly cited;
unverified against text]**). Also first to articulate the statistical model as a binomial
process. *Relevance:* Lees' 13 control cases are, in effect, the first "calibration set"; the
problem he faced — rate heterogeneity across pairs larger than binomial noise — is exactly
what clepsydra's hierarchical rate prior must absorb.

**Bergsland, K. & Vogt, H. (1962).** On the validity of glottochronology. *Current Anthropology*
3(2): 115–153. https://doi.org/10.1086/200264
Showed that Icelandic vs. Norwegian (Riksmål), Georgian, and Armenian retain far more (or less)
than predicted: Icelandic's retention since Old Norse implies an age far too young. The
canonical demonstration of *lineage-specific rate variation*. Starostin later argued that much of
the Riksmål discrepancy disappears once Danish/German loans are removed (Vasilyev 2009 handout,
secondary). *Relevance:* (i) motivates a branch-rate (relaxed clock) model and reporting its
variance; (ii) is a clean held-out test case: run clepsydra on Old Norse/Icelandic/Norwegian with
and without the loan-exclusion option and report whether the contamination mixture recovers
Starostin's correction.

**Dyen, I., James, A. T. & Cole, J. W. L. (1967).** Language divergence and estimated word
retention rate. *Language* 43(1): 150–171. https://doi.org/10.2307/411390
Introduced meaning-specific retention rates estimated from Austronesian data. *Relevance:* the
direct ancestor of clepsydra's concept-specific rates.

**Sankoff, D. (1973).** Mathematical developments in lexicostatistic theory. In T. A. Sebeok (ed.),
*Current Trends in Linguistics* 11: 93–114. The Hague: Mouton.
https://doi.org/10.1515/9783111418797-007
Formal treatment of synonymy, borrowing, meaning-specific rates (gamma-distributed), and the
consequent nonlinearity of log-retention with time: with heterogeneous rates, the
expected shared fraction decays more slowly than exponential, so a single-rate formula
*underestimates* deep ages. *Relevance:* the theoretical justification for clepsydra's gamma
(or per-concept) rate model; also provides closed forms for the pairwise expectation under
gamma rate heterogeneity, E[C] = (1 + 2θ t/α)^(-α) for shape α (standard result;
**[cite exact equation number after checking]**).

**Kruskal, J. B., Dyen, I. & Black, P. (1973).** Some results from the vocabulary method of
reconstructing language trees. In I. Dyen (ed.), *Lexicostatistics in Genetic Linguistics*,
30–55. https://doi.org/10.1515/9783110880847-003
Meaning-specific rates, tree-based (not pairwise) estimation. *Relevance:* early tree likelihood
using meaning rates; historical precedent for T1/T2.

**Embleton, S. M. (1985).** Lexicostatistics applied to the Germanic, Romance, and Wakashan
families. *WORD* 36(1): 37–60. https://doi.org/10.1080/00437956.1985.11435862
**Embleton, S. M. (1986).** *Statistics in Historical Linguistics.* Bochum: Brockmeyer
(Quantitative Linguistics 30). No DOI.
Simulation-based treatment of lexicostatistics with borrowing modelled as a function of
geographic adjacency, and meaning-specific rates; shows that borrowing biases ages
*downwards*. *Relevance:* first explicit borrowing-aware lexicostatistic model; the clepsydra
contamination mixture is a statistical cousin. Embleton's simulations are a template for the
"borrowing → too young" validation scenario.

**Dyen, I., Kruskal, J. B. & Black, P. (1992).** An Indoeuropean classification: a
lexicostatistical experiment. *Transactions of the American Philosophical Society* 82(5).
https://doi.org/10.2307/1006517
The 200-meaning, ~95-variety IE dataset with cognate judgments, used (via Dyen et al. 1997
digital version) by Gray & Atkinson 2003, Pagel et al. 2007 and many others. *Relevance:* a
standard benchmark dataset; known issues (loan coding, lect selection, missing ancient
languages) make it a good test for robustness of clepsydra's contamination model.

**Starostin, S. A. (1999).** Comparative-historical linguistics and lexicostatistics. In V.
Shevoroshkin & P. Sidwell (eds), *Historical Linguistics and Lexicostatistics*, 3–50.
Melbourne: Association for the History of Language. No DOI. (Russian original 1989.)
Starostin's "recalibrated glottochronology": (a) borrowings are removed before counting, not
treated as replacements; (b) the replacement rate itself grows with time since a word entered
the list (words "age"), giving C(t) = exp(−λ C t²)-type decay, i.e. t = √(−ln C /(λ C)), with
λ ≈ 0.05 per millennium² for the 100-item list (the formula and constant as reported by
secondary sources — Novotná & Blažek 2007, *Baltistica* 42; Vasilyev 2009 Swadesh-centenary
handout — **[exact functional form and constant unverified against the primary text]**);
(c) dates are computed from *reconstructed* protolanguage lists step by step, not from modern
leaves. *Relevance:* Starostin's formula is a non-exponential "decay law" producing older dates
at low C than Swadesh's; clepsydra's `docs/decay-laws.md` should include it as a named
alternative and compare it against the gamma-mixture closed form (which produces a similar
slowing of decay for a principled reason). Point (c) is a different estimator altogether.

**Kassian, A., Starostin, G., Dybo, A. & Chernov, V. (2010).** The Swadesh wordlist: an attempt at
semantic specification. *Journal of Language Relationship* 4: 46–89. **[DOI unverified;
JLR DOIs are 10.31826/…]**
**Starostin, G. (2010).** Preliminary lexicostatistics as a basis for language classification: a new
approach. *Journal of Language Relationship* 3: 79–116. **[DOI unverified]**
The Moscow (Global Lexicostatistical Database) protocol: strict semantic definitions for each
Swadesh slot, one "most basic" item per slot, synonym rules, and etymological notes. *Relevance:*
clepsydra's contamination and multistate likelihoods both assume a "one-word-per-slot" state;
the Moscow semantic specifications are the best documented way to reduce polymorphism and
semantic drift — worth citing in input-format guidance.

**Kassian, A. S., Zhivlov, M., Starostin, G., Trofimov, A. A., Kocharov, P. A., Kuritsyna, A. &
Saenko, M. N. (2021).** Rapid radiation of the inner Indo-European languages: an advanced
approach to Indo-European lexicostatistics. *Linguistics* 59(4): 949–979.
https://doi.org/10.1515/ling-2020-0060
Uses 13 110-item *reconstructed* protolanguage/ancient lists (Proto-Germanic, Hittite, etc.),
with "onomasiological reconstruction", "derivational drift elimination" and "homoplastic
optimization"; Bayesian, parsimony and NJ give a consensus with Inner IE multifurcating into four
branches c. 3357–2162 BCE. *Relevance:* (i) shows the effect of reducing input noise vs.
modelling it — a data-cleaning alternative to the contamination mixture; (ii) a benchmark for
inner-IE node ages; (iii) a star-like polytomy is an important case for clepsydra's fixed-topology
mode: the tool should permit multifurcations (or near-zero internal branches) without forcing
artefactual resolution.

**Egorov, I. M., Dybo, A. V. & Kassian, A. S. (2022).** Phylogeny of the Turkic languages inferred
from basic vocabulary: limitations of the lexicostatistical methods in an intensive contact
situation. *Journal of Language Evolution* 7(1): 16–39. https://doi.org/10.1093/jole/lzac006
Same Moscow protocol on Turkic; Chuvash first outlier, then an eight-way multifurcation of
Nuclear Turkic; Siberian Turkic is areal, not genealogical. *Relevance:* explicitly a family where
contact defeats tree signal; a test case for clepsydra's contamination mixture and for a
"tree-likeness" diagnostic.

**Tadmor, U., Haspelmath, M. & Taylor, B. (2010).** Borrowability and the notion of basic vocabulary.
*Diachronica* 27(2): 226–246. https://doi.org/10.1075/dia.27.2.04tad
**Haspelmath, M. & Tadmor, U. (eds) (2009).** *Loanwords in the World's Languages: A Comparative
Handbook.* Berlin: De Gruyter Mouton. https://doi.org/10.1515/9783110218442
WOLD (41 languages, 1,460 meanings): nouns more borrowable than verbs/adjectives, content more
than function words; combines borrowability, age, simplicity and representativeness into the
100-item Leipzig–Jakarta list. *Relevance:* per-concept *borrowability* is an empirically
available covariate. Clepsydra's contamination mixture could take a concept-specific prior on the
contamination probability from WOLD borrowing rates rather than a single global weight.

**Bowern, C., Epps, P., Gray, R., Hill, J., Hunley, K., McConvell, P. & Zentz, J. (2011).** Does
lateral transmission obscure inheritance in hunter-gatherer languages? *PLoS ONE* 6(9): e25195.
https://doi.org/10.1371/journal.pone.0025195
122 languages (N. Australia, NW Amazonia, California/Great Basin), 204-item list. Mean loan rate in
basic vocabulary 5.06% (median 2.49%, range 0–48%). Amazonia: no language above 4%, most 1–2%.
Australia highest (mean 9.4%; Gurindji ~48%, Mudburra 38%, Gooniyandi 33%). WOLD comparison mean
10.24% (median 5.3%). *Relevance:* the best empirical prior for the contamination weight: a
default Beta prior centred near 0.03–0.05 with a heavy right tail; region-specific defaults
(Amazonia low, N. Australia high). Directly relevant to Arawan.

### 1.2 Rates of lexical change: meanings, lineages, time

**Pagel, M., Atkinson, Q. D. & Meade, A. (2007).** Frequency of word-use predicts rates of lexical
evolution throughout Indo-European history. *Nature* 449: 717–720.
https://doi.org/10.1038/nature06176
200 meanings × 87 IE languages; per-meaning replacement rates estimated on a posterior sample of
trees; corpus frequency (English, Spanish, Russian, Greek) explains ~50% of variance in rates;
frequent words change slowly. (Commonly cited half-lives range from ~750 to ~10,000 years
**[range unverified against the paper]**.) *Relevance:* (i) the per-concept rate table is the most
widely reused rate prior; (ii) frequency is a covariate that could make concept rates
transferable to families without calibrations: a regression prior log θ_c = β0 + β1 log f_c + ε.

**Calude, A. S. & Pagel, M. (2011).** How do we use language? Shared patterns in the frequency of
word use across 17 world languages. *Phil. Trans. R. Soc. B* 366: 1101–1107.
https://doi.org/10.1098/rstb.2010.0315
Frequencies of the 200 Swadesh meanings in 17 languages from six families correlate at 0.73 on
average; the first principal component explains 70% of variance; worldwide frequency loadings plus
part of speech predict 46% of IE rate variance and 30% of variance in a 110-word stability ranking
across 14 families. *Relevance:* evidence that concept-rate *rankings* transfer across families —
the empirical basis for a cross-family hierarchical concept-rate prior.

**Pagel, M., Atkinson, Q. D., Calude, A. S. & Meade, A. (2013).** Ultraconserved words point to deep
language ancestry across Eurasia. *PNAS* 110(21): 8471–8476. https://doi.org/10.1073/pnas.1218726110
Claims a ~14,450-year "Eurasiatic" tree from 23 putatively ultraconserved meanings in LWED
reconstructions; words used >1/1000 in speech 7–10× more likely to show deep ancestry.
**Heggarty, P. (2013).** Ultraconserved words and Eurasiatic? The "faces in the fire" of language
prehistory. *PNAS* 110(35): E3254. https://doi.org/10.1073/pnas.1309114110
Critique: input "cognates" are proto-form look-alikes of doubtful status; chance resemblance
and coding choice drive the result. *Relevance:* a cautionary example of the deep-time regime
where the signal is at chance level. Clepsydra should report a *chance-resemblance floor* and
refuse (or flag) node ages whose shared-state fraction is indistinguishable from it.

**Pagel, M. & Meade, A. (2018).** The deep history of the number words. *Phil. Trans. R. Soc. B*
373: 20160517. https://doi.org/10.1098/rstb.2016.0517 (Correction: Pagel 2025,
https://doi.org/10.1098/rstb.2025.0364 — **content of the correction not checked**.)
Low numerals (1–5) replaced 3.5–20× more slowly than average Swadesh items, replicated in IE,
Bantu, Austronesian; number words stable for 10,000 to >100,000 years. *Relevance:* extremely slow
concepts contribute almost no dating information but strongly anchor deep relatedness; in a
contamination model they are the concepts where a mismatch is most likely due to borrowing or
coding error (e.g. numeral borrowing in Amazonian/Australian languages — numerals are frequently
borrowed or absent above 3 in Amazonia). Suggests concept-specific contamination priors.

**Vejdemo, S. & Hörberg, T. (2016).** Semantic factors predict the rate of lexical replacement of
content words. *PLoS ONE* 11(1): e0147924. https://doi.org/10.1371/journal.pone.0147924
Regression on Pagel et al. 2007 IE rates: frequency (−), number of synonyms (+, strongest semantic
factor), imageability (−), number of senses (weak −); no word-class effect once controlled;
R² = 0.34 (N = 117). Unique ΔR²: frequency ~16%, synonyms ~12.5%, imageability ~6%, senses ~3%.
*Relevance:* additional covariates for a concept-rate regression prior; synonym count also
predicts *polymorphism*, which clepsydra must handle (multiple states per slot).

**Monaghan, P. (2014).** Age of acquisition predicts rate of lexical evolution. *Cognition* 133(2):
530–534. https://doi.org/10.1016/j.cognition.2014.08.007
Age of acquisition predicts IE replacement rates beyond frequency. *Relevance:* another
covariate (contested by Vejdemo & Hörberg's controlled model).

**Atkinson, Q. D., Meade, A., Venditti, C., Greenhill, S. J. & Pagel, M. (2008).** Languages evolve
in punctuational bursts. *Science* 319: 588. https://doi.org/10.1126/science.1149683
Bantu, IE, Austronesian: 10–33% of vocabulary divergence arises at language-splitting events
(path-length/node-count correlation). *Relevance:* a direct model extension: a per-node "burst"
parameter (instantaneous replacement probability at each split). Without it, a clock calibrated on
shallow, node-dense parts of the tree overestimates rate per unit time and makes node-sparse
deep branches too young (or vice versa, depending on calibration placement). Note the
node-density artefact must be corrected (Venditti et al.).

**Greenhill, S. J., Atkinson, Q. D., Meade, A. & Gray, R. D. (2010).** The shape and tempo of
language evolution. *Proc. R. Soc. B* 277: 2443–2450. https://doi.org/10.1098/rspb.2010.0051
Typological features evolve at rates similar to basic vocabulary but much less tree-like; lexical
rates correlate between Austronesian and IE, structural rates do not. *Relevance:* supports
transferring *lexical* concept-rate priors across families (rates correlate), not structural ones.

**Greenhill, S. J., Wu, C.-H., Hua, X., Dunn, M., Levinson, S. C. & Gray, R. D. (2017).** Evolutionary
dynamics of language systems. *PNAS* 114(42): E8822–E8829. https://doi.org/10.1073/pnas.1700388114
Dirichlet-process mixture over rates of 81 Austronesian languages: most grammatical features
change *faster* than basic vocabulary; grammar shows less schismogenesis, more homoplasy and
contact bursts; a stable core exists in both. *Relevance:* (i) a DP mixture over concept rates is a
flexible alternative to a parametric gamma — clusters of "slow/medium/fast" concepts; (ii)
undercuts the hope of using grammar to date deeper than vocabulary.

**Bromham, L., Hua, X., Fitzpatrick, T. G. & Greenhill, S. J. (2015).** Rate of language evolution is
affected by population size. *PNAS* 112(7): 2097–2102. https://doi.org/10.1073/pnas.1419704112
Polynesian (20 closely related languages): larger populations gain new words faster; smaller
populations lose words faster. *Relevance:* gain and loss respond differently to demography —
supports clepsydra's binary gain/loss CTMC with *separate* branch-specific gain and loss
multipliers, and a branch-rate covariate (log population) as an option.

**Greenhill, S. J., Hua, X., Welsh, C. F., Schneemann, H. & Bromham, L. (2018).** Population size and
the rate of language evolution: a test across Indo-European, Austronesian, and Bantu languages.
*Frontiers in Psychology* 9: 576. https://doi.org/10.3389/fpsyg.2018.00576
153 sister pairs: word loss faster in smaller IE languages; no significant pattern in Austronesian
or Bantu. *Relevance:* population effects are weak and family-specific — do not hard-wire them;
offer as optional covariate with a shrinkage prior centred on zero.

**Wichmann, S. & Holman, E. W. (2009).** Population size and rates of language change. *Human
Biology* 81(2–3): 259–274. https://doi.org/10.3378/027.081.0308
ASJP-based test: no detectable effect of population size on lexical change rates
(**[finding summarized from memory; verify]**). **Nettle, D. (1999).** Is the rate of linguistic
change constant? *Lingua* 108: 119–136. https://doi.org/10.1016/S0024-3841(98)00047-3 — argued
small populations change faster (simulation + typological reasoning). *Relevance:* contradictory
literature → treat as an explorable controversy, not a default.

**Cathcart, C. A. (2024).** Rate variation in language change: toward distributional phylogenetic
modeling. In *Biocultural Evolution: An Agenda for Integrative Approaches* (MPG/MPI-EVA volume),
179–202. https://doi.org/10.51315/9783935751384.008 (**content not read; title/venue from
Crossref**). *Relevance:* argues for modelling rate *distributions* across characters and lineages
jointly — the same design space as clepsydra's hierarchical rate prior.

**Cathcart, C. A. (2018).** Modeling linguistic evolution: a look under the hood. *Linguistics
Vanguard* 4(1). https://doi.org/10.1515/lingvan-2017-0043 (abstract not retrieved). Tutorial-level
critique of model assumptions in linguistic phylogenetics.

### 1.3 Bayesian phylolinguistic dating: models and software

**Gray, R. D. & Atkinson, Q. D. (2003).** Language-tree divergence times support the Anatolian
theory of Indo-European origin. *Nature* 426: 435–439. https://doi.org/10.1038/nature02029
Binary cognate presence/absence (Dyen et al. 1992 + Hittite, Tocharian), restriction-site
(gain/loss) model with gamma rate heterogeneity, MrBayes, then penalized-likelihood rate
smoothing with 14 calibrations. Root 7,800–9,800 BP (see Part 2). *Relevance:* the template for
clepsydra's T2 (binary CTMC); the sensitivity analyses (coding of loans, removing calibrations)
are a model for the "controversy explorer" outputs.

**Atkinson, Q. D. & Gray, R. D. (2005).** Curious parallels and curious connections — phylogenetic
thinking in biology and historical linguistics. *Systematic Biology* 54(4): 513–526.
https://doi.org/10.1080/10635150590950317
**Atkinson, Q., Nicholls, G., Welch, D. & Gray, R. (2005).** From words to dates: water into wine,
mathemagic or phylogenetic inference? *Transactions of the Philological Society* 103(2): 193–219.
https://doi.org/10.1111/j.1467-968X.2005.00151.x
Answers glottochronology objections using the Ringe et al. data (24 languages, 12 node constraints,
20 dated ancient languages) and Dyen et al. (87 languages), under a finite-sites model with r8s rate
smoothing and under stochastic Dollo (TraitLab); both support 7,800–9,800 BP. Validation is **by
simulation**, not calibration leave-one-out: with a true root of 8,680, 20% borrowing gave slight
underestimation and 100% borrowing about 30% underestimation. *Relevance:* the earliest quantitative
"borrowing → too young" calibration curve; a template scenario for clepsydra's simulation suite.

**Atkinson, Q. D. & Gray, R. D. (2006).** How old is the Indo-European language family? Illumination
or more moths to the flame? In P. Forster & C. Renfrew (eds), *Phylogenetic Methods and the
Prehistory of Languages*, 91–109. Cambridge: McDonald Institute. No DOI.
Re-analysis with alternative data and assumptions; the root date is given in secondary sources only as
"ca. 7000 BCE" **[unverified; one listing gives the subtitle as "Progress or…"]**.

**Nicholls, G. K. & Gray, R. D. (2008).** Dated ancestral trees from binary trait data and their
application to the diversification of languages. *J. R. Stat. Soc. B* 70(3): 545–566.
https://doi.org/10.1111/j.1467-9868.2007.00648.x
Stochastic Dollo (SD): each cognate class born once (Poisson birth) and lost independently;
correct conditioning on removal of singletons; uniform marginal root prior; predictive checks
against external data. Finding: node ages robust, topology posteriors not. *Relevance:* SD is the
binary analogue of the novel-state assumption in clepsydra's T1 (each state arises once). The
"ages robust, topology less so" result supports clepsydra's fixed-topology design.

**Ryder, R. J. & Nicholls, G. K. (2011).** Missing data in a stochastic Dollo model for binary trait
data, and its application to the dating of Proto-Indo-European. *J. R. Stat. Soc. C* 60(1): 71–92.
https://doi.org/10.1111/j.1467-9876.2010.00743.x (arXiv:0908.1735)
Missing-at-random data handled properly (Hittite 12% missing was poorly predicted when missing
was coded as absent); "catastrophes" (instantaneous multi-trait change on a branch) for
punctuational heterogeneity; Bayes-factor tests of 30 historical constraints, three rejected;
PIE posterior centred 8,400 BP, 95% HPD 7,100–9,800 BP. Implemented in TraitLab. *Relevance:*
(i) clepsydra must treat missing cells as marginalized, never as absent — a common user error with
binary encodings; (ii) catastrophes = branch-located bursts, an alternative to per-node bursts;
(iii) **calibration consistency testing** (BF per calibration) is directly adoptable.

**Kelly, L. J. & Nicholls, G. K. (2017).** Lateral transfer in stochastic Dollo models. *Annals of
Applied Statistics* 11(2): 1146–1168. https://doi.org/10.1214/17-AOAS1040
Borrowing among contemporaneous lineages in SD; likelihood requires solving large ODE systems;
improved fit on Eastern Polynesian. *Relevance:* the principled but expensive alternative to
clepsydra's contamination mixture; the mixture is a cheap approximation that ignores *who*
borrowed from whom. A middle ground: contamination probability that scales with the duration
that two lineages coexist.

**Bouckaert, R., Lemey, P., Dunn, M., Greenhill, S. J., Alekseyenko, A. V., Drummond, A. J., Gray,
R. D., Suchard, M. A. & Atkinson, Q. D. (2012).** Mapping the origins and expansion of the
Indo-European language family. *Science* 337: 957–960. https://doi.org/10.1126/science.1219669
Joint phylogeographic (relaxed random walk) + lexical model, BEAST, binary covarion / SD variants,
relaxed clock, ancient languages as dated tips, calibrations. Anatolian homeland supported. Root
date in Part 2. *Relevance:* standard set-up (binary + covarion + UCLN + calibrations) that King
2026 later tested; also the target of Rönchen, Billing & Wiklund's simulation critique.

**Chang, W., Cathcart, C., Hall, D. & Garrett, A. (2015).** Ancestry-constrained phylogenetic
analysis supports the Indo-European steppe hypothesis. *Language* 91(1): 194–244.
https://doi.org/10.1353/lan.2015.0005
Re-curated IE dataset; when ancient languages are constrained to be direct ancestors (e.g. Latin
of Romance, Vedic of Indic, Old Irish of Irish), the root becomes much younger (steppe-compatible),
because unconstrained trees place ancient languages on side branches, inflating branch lengths
and rates' apparent slowdown. Argues "advergence" (parallel retention/borrowing among
sister lineages) biases dates. *Relevance:* the single most important design point for clepsydra's
calibrations: **ancestor-as-node calibrations (sampled ancestors)** vs. tip-dated sister calibrations
give systematically different results. With a fixed topology clepsydra can implement
ancestor constraints exactly: an ancient language's data attached at an internal node at a known
time (branch length zero to the ancestor).

**Heggarty, P., Anderson, C., Scarborough, M., King, B., Bouckaert, R., et al. (33 authors, ending
… Greenhill, S. J., Kühnert, D. & Gray, R. D.) (2023).** Language trees with sampled ancestors support
a hybrid model for the origin of Indo-European languages. *Science* 381(6656): eabg0818.
https://doi.org/10.1126/science.abg0818
IE-CoR: 161 languages (109 modern, 52 ancient/historical, tip-dated), 170 meanings, explicit
cognate/loan coding protocol; sampled-ancestor (FBD) prior so attested ancient languages *may* sit on
the line of descent (few turn out to be direct ancestors); binary covarion with rates by meaning
group, relaxed clock. Root ~8,120 BP (6,740–9,610). Chang et al.'s steppe result attributed to
their constraints rather than data. Critiqued by **Kassian, A. & Starostin, G. (2025).**
Do 'language trees with sampled ancestors' really support a 'hybrid model' for the origin of
Indo-European? *Humanities and Social Sciences Communications* 12: 682.
https://doi.org/10.1057/s41599-025-04986-7 (content not read here). *Relevance:* IE-CoR is the dataset clepsydra already uses
(`examples/`); its meaning-level partitioning, sampled ancestors and loan-flagging are directly
reusable. The Chang vs. Heggarty disagreement is the canonical "controversy" to expose.

**Ritchie, A. M. & Ho, S. Y. W. (2019).** Influence of the tree prior and sampling scale on Bayesian
phylogenetic estimates of the origin times of language families. *Journal of Language Evolution*
4(2): 108–123. https://doi.org/10.1093/jole/lzz005
Four lexical datasets: root ages robust to tree prior (less so than for genes); speciation priors
preferred over constant coalescent; narrower sampling fits birth–death worse. *Relevance:* with a
fixed topology clepsydra still needs a prior on node ages; birth–death (with sampling fraction) is
the defensible default, and sensitivity to it should be reported.

**Rama, T. (2018).** Three tree priors and five datasets: a study of the effect of tree priors in
Indo-European phylogenetics. *Language Dynamics and Change* 8(2): 182–218.
https://doi.org/10.1163/22105832-00802005
Uniform, fossilized birth–death (FBD) and coalescent priors on five IE datasets; Steppe-compatible
root ages for two priors; uniform or FBD preferred by Bayes factors. *Relevance:* contrasts with
Ritchie & Ho; IE root depends on tree prior by >1 millennium in some datasets. Clepsydra should
expose the node-age prior explicitly (uniform-on-ages vs. birth–death vs. FBD-like with dated
ancestors).

**Hoffmann, K., Bouckaert, R., Greenhill, S. J. & Kühnert, D. (2021).** Bayesian phylogenetic
analysis of linguistic data using BEAST. *Journal of Language Evolution* 6(2): 119–135.
https://doi.org/10.1093/jole/lzab005
Tutorial: substitution models (binary CTMC, covarion, SD), clocks, tree priors, path/nested
sampling. *Relevance:* defines the "standard" against which clepsydra is compared; documentation
of BEAST defaults that King 2026 shows to be problematic.

**Maurits, L., de Heer, M., Honkola, T., Dunn, M. & Vesakoski, O. (2020).** Best practices in
justifying calibrations for dating language families. *Journal of Language Evolution* 5(1): 17–38
(online 2019). https://doi.org/10.1093/jole/lzz009
Calibration theory, critique of previously used calibrations, recommended justification practices,
worked Uralic calibration set using loanword evidence (e.g. contact with dated Indo-Iranian or
Germanic stages gives bounds). *Relevance:* (i) clepsydra's calibration file should require a
justification field and a calibration *type* (attestation, archaeology, loanword-stratum,
geography); (ii) loanword-stratum bounds are the main calibration resource for families like
Arawan that lack written history (e.g. dated contact layers).

**Neureiter, N., Ranacher, P., Efrat-Kowalsky, N., Kaiping, G. A., Weibel, R., Widmer, P. &
Bouckaert, R. R. (2022).** Detecting contact in language trees: a Bayesian phylogenetic model with
horizontal transfer. *Humanities and Social Sciences Communications* 9: 205.
https://doi.org/10.1057/s41599-022-01211-7
contacTrees (BEAST 2): contact edges between clades along which words are borrowed; simulation
shows ignoring contact biases tree height, rates and topology; IE case study recovers known
loans and gives a tree height "significantly closer to well-attested estimates". *Relevance:*
quantitative evidence that unmodelled borrowing distorts *age*, and a model design (sparse contact
edges with per-edge borrowing probability) for a future clepsydra extension.

**Rönchen, P., Wiklund, T. & Hammarström, H. (2024).** Likelihood calculation in a multistate model
of vocabulary evolution for linguistic dating. *Language Dynamics and Change* 14(1): 1–41.
https://doi.org/10.1163/22105832-bja10032
Infinite-state ("novel-state") replacement model on a tree; general recursive likelihood
algorithm; argues greater realism than SD; case study. *Relevance:* the basis of clepsydra's T1.

**Rönchen, P. (2026).** *Evaluating models of vocabulary evolution for linguistic dating and
subgrouping.* PhD thesis, Uppsala University, defended 7 January 2026 (urn:nbn:se:uu:diva-571691).
Compares finite-state Markov chains (covarion etc.), the multistate model and SD on simulated and
empirical data; abstract: finite-state chains such as the covarion "can overfit linguistic data
sets and infer incorrect language trees." (**Thesis body not read.**)
**Rönchen, P., Billing, O. & Wiklund, T. (2023).** Using simulated data to evaluate models of
Indo-European vocabulary evolution. Abstract, ICHL 26 (Heidelberg), workshop 3. Simulated datasets
tailored to the Chang et al. and Bouckaert et al. pipelines (with loans, near-consensus topology,
varying rates): "both … methodologies fail to correctly infer the age of Indo-European that was
used to create our simulated data sets." *Relevance:* the strongest current argument for
clepsydra's design choice of a multistate likelihood + simulation-based validation; clepsydra
should reproduce these simulation scenarios as a validation suite.

**King, B. (2026).** Testing the validity and adequacy of linguistic phylogenetic analyses. *PLOS
Computational Biology* 22(5): e1014312. https://doi.org/10.1371/journal.pcbi.1014312
SBC (100 replicates, LinguaPhylo → BEAST 2.7.8) of binary cognates + ascertainment correction +
binary covarion + UCLN + Yule. (1) all-absent retained: calibrated; (2) ascertainment-corrected,
partition rates unweighted: calibrated; (3) BEAST default reweighting of partition rates by cognate
count: **fails** — partition rates and birth rate underestimated, *tree height and length
overestimated*. Empirical (Indo-Iranic subset of IE-CoR, 170 meanings, 5 partitions by
cognate-set count): meaning-weighted partition rates best by path sampling, cognate-weighted worst;
on real data tree height unchanged across schemes. Posterior predictive checks (35
Germanic/Italic/Celtic languages, 1,251 cognate sets, 8 statistics): 4 of 8 show p < 0.05 —
misfit attributed to semantic shift and non-independence of cognates within a meaning (one lexeme
per slot). Recommends meaning-weighting, explicit Dirichlet priors, and the infinite-state
multistate model. *Relevance:* very high. (i) Independent validation that the binary encoding's
non-independence is a real misfit — the case for T1; (ii) clepsydra's T2 must avoid the
cognate-count weighting pitfall; (iii) the eight posterior-predictive statistics are an
off-the-shelf adequacy panel (number of prevalent cognates, singletons, variance in prevalence,
pairwise-distance max/min/mean, mean presence, variance in presence).

**King, B. (2026).** Posterior simulation-based calibration tests of phylogenetic dating methods.
bioRxiv preprint. https://doi.org/10.64898/2026.04.14.718426
Posterior SBC (simulate from the posterior, re-infer) of BEAST 2 tip-dating on an IE vocabulary
dataset: good calibration, but posterior-predictive replicates add *no* precision to node ages,
consistent with theoretical limits on node-age identifiability (rate–time confounding).
*Relevance:* clepsydra already does prior SBC; posterior SBC is a cheap addition specific to a
user's dataset. The identifiability result means the posterior width on ages is floored by
calibration/rate-prior information, not by data size — users must be told this.

**Häuser, L., Jäger, G. & Stamatakis, A. (2026).** A systematic exploration of current limitations
of cognate-based phylogenetic inference. *Open Research Europe* 5: 258 (v3).
https://doi.org/10.12688/openreseurope.20351.3
Lexibench cognate data vs 10,406 TreeBASE DNA matrices: similar per-matrix information; binary
encoding violates independence within a concept; multi-valued per-concept models (MK, GTR, new
COG/COGs) overparameterized for current data sizes (κ = 3 subsets: median relative CV error
0.21–0.40 vs 0.03 for binary); recommend binary with simple models until data grow.
*Relevance:* a counterpoint to King/Rönchen: *finite* multistate per-concept models overfit —
supports clepsydra's choice of the infinite-state model (one parameter per concept, no state
frequencies) rather than finite multistate.

**Xu, Y., Douglas, J., Bouckaert, R. R. & Drummond, A. J. (2026).** A hierarchical clock-mixture
model for Bayesian phylogenetic dating. bioRxiv. https://doi.org/10.64898/2026.07.15.738815
Analytical averaging over strict vs relaxed, then UCLN vs autocorrelated on a shared branch-rate
vector (RelaxClockAveraging, BEAST 2). *Relevance:* clepsydra's relaxed clock could adopt a
spike-and-slab on the clock variance (strict ↔ relaxed) to let data decide, reporting the
posterior probability of a strict clock — a natural "glottochronology test".

**Bouckaert, R., Heled, J., Kühnert, D., Vaughan, T., Wu, C.-H., Xie, D., Suchard, M. A., Rambaut,
A. & Drummond, A. J. (2014).** BEAST 2: a software platform for Bayesian evolutionary analysis.
*PLoS Comput. Biol.* 10(4): e1003537. https://doi.org/10.1371/journal.pcbi.1003537
**Drummond, A. J., Ho, S. Y. W., Phillips, M. J. & Rambaut, A. (2006).** Relaxed phylogenetics and
dating with confidence. *PLoS Biology* 4(5): e88. https://doi.org/10.1371/journal.pbio.0040088
UCLN clock. **Heath, T. A., Huelsenbeck, J. P. & Stadler, T. (2014).** The fossilized birth–death
process for coherent calibration of divergence-time estimates. *PNAS* 111(29): E2957–E2966.
https://doi.org/10.1073/pnas.1319091111 . **Gavryushkina, A., Welch, D., Stadler, T. &
Drummond, A. J. (2014).** Bayesian inference of sampled ancestor trees for epidemiology and fossil
calibration. *PLoS Comput. Biol.* 10(12): e1003919. https://doi.org/10.1371/journal.pcbi.1003919
*Relevance:* the biological machinery behind sampled ancestors (Heggarty 2023) and coherent
calibration (avoid the "calibration densities interact with tree prior" problem: the effective
prior on a calibrated node is not the user-specified density; clepsydra should *display the
effective joint prior* by sampling from the prior).

**Ho, S. Y. W., Phillips, M. J., Cooper, A. & Drummond, A. J. (2005).** Time dependency of molecular
rate estimates and systematic overestimation of recent divergence times. *Mol. Biol. Evol.* 22(7):
1561–1568. https://doi.org/10.1093/molbev/msi145 — biological analogue of time-dependent rates
(saturation + purifying selection); in linguistics, the analogue is Starostin's rate acceleration
and the saturation of shared fraction at chance level. No equivalent linguistic study of rate
time-dependency across calibration depths was found during this survey (**gap**).

### 1.4 Distance-based automated dating (ASJP and relatives)

**Serva, M. & Petroni, F. (2008).** Indo-European languages tree by Levenshtein distance. *EPL* 81:
68005. https://doi.org/10.1209/0295-5075/81/68005
**Petroni, F. & Serva, M. (2008).** Language distance and tree reconstruction. *J. Stat. Mech.*
P08012. https://doi.org/10.1088/1742-5468/2008/08/P08012
Normalized Levenshtein distance over 200-item Swadesh lists for 50 IE languages; tree similar to
Gray & Atkinson; glottochronological conversion of distance to time with a fitted constant.
**Serva, M. (2025).** Evolution of the lexicon: a probabilistic point of view. arXiv:2510.22220.
Argues that Swadesh-type estimators have intrinsic probabilistic accuracy limits even under ideal
assumptions, and that modelling gradual phonetic modification of words (not just replacement)
substantially improves precision. *Relevance:* a "sub-cognate" signal (phonetic distance within
cognate sets) as an extra clock is a potential later extension; the accuracy-limit argument should
be reproduced as a theoretical floor in clepsydra's pairwise documentation.

**Holman, E. W., Brown, C. H., Wichmann, S., Müller, A., Velupillai, V., Hammarström, H., Sauppe, S.,
Jung, H., Bakker, D., Brown, P., Belyaev, O., Urban, M., Mailhammer, R., List, J.-M. & Egorov, D.
(2011).** Automated dating of the world's language families based on lexical similarity. *Current
Anthropology* 52(6): 841–875. https://doi.org/10.1086/662127
ASJP (40-item lists, LDND distance) calibrated on 52 language groups with historical,
epigraphic or archaeological dates; discrepancies between estimated and calibration dates are on
average 29% of the estimated dates — **an in-sample figure** (the same 52 points fit the formula; the
authors call it an upper bound on accuracy). Formula t = (log s − log s₀)/(2 log r), s₀ = 0.92,
r = 0.72 per millennium. Dates for most of the world's families. (Further details in
Part 2.) *Relevance:* the main *uncalibrated-family* dating baseline, with an empirically measured
error. Clepsydra's cross-family validation should be benchmarked against this 29% figure using
the same calibration set where possible.

**Holman, E. W., Wichmann, S., Brown, C. H., Velupillai, V., Müller, A. & Bakker, D. (2008).**
Explorations in automated language classification. *Folia Linguistica* 42(3–4): 331–354.
https://doi.org/10.1515/flin.2008.331 — stability ranking of Swadesh items from ASJP; origin of the
40-item list. **Wichmann, S., Holman, E. W., Bakker, D. & Brown, C. H. (2010).** Evaluating
linguistic distance measures. *Physica A* 389: 3632–3639. https://doi.org/10.1016/j.physa.2010.05.011
**Wichmann, S. & Holman, E. W. (2013).** Languages with longer words have more lexical change. In
L. Borin & A. Saxena (eds), *Approaches to Measuring Linguistic Differences*, 249–282.
https://doi.org/10.1515/9783110305258.249 — word length as a rate covariate.

**Rama, T. (2013).** Phonotactic diversity predicts the time depth of the world's language families.
*PLoS ONE* 8(5): e63238. https://doi.org/10.1371/journal.pone.0063238 — alternative, non-lexical
clock trained on ASJP calibrations.

**Rama, T. & Wichmann, S. (2020).** A test of Generalized Bayesian dating: a new linguistic dating
method. *PLoS ONE* 15(8): e0236522. https://doi.org/10.1371/journal.pone.0236522
Automated cognates; Bayesian trees per family *without* family-specific priors or calibrations;
a gamma regression of known time depth on inferred tree length fitted to 30 phylogenetically
independent calibration points; predicts 1,287 root and subgroup dates for 116 families. Results
"similar to published Bayesian studies"; but "in terms of errors and correlations with known dates,
ASJP works better than the new method and both work better than automated glottochronology."
*Relevance:* the closest published analogue to what clepsydra wants for Arawan: transfer of a
time-scale learnt across families. Two lessons: (i) a tree-length → time regression is a calibration
transfer that ignores within-family rate variation; (ii) it did *not* beat ASJP distances, so
gains must come from better rate modelling, not only from Bayesian machinery.

**Rama, T., List, J.-M., Wahle, J. & Jäger, G. (2018).** Are automatic methods for cognate detection
good enough for phylogenetic reconstruction in historical linguistics? *Proc. NAACL-HLT 2018*,
393–400. https://doi.org/10.18653/v1/N18-2063 — automated cognates (LexStat, SVM, etc.) give
trees and dates close to expert cognates on several families. *Relevance:* clepsydra could accept
automatically clustered cognates with a cognate-error rate parameter (another contamination source).

**Jäger, G. (2018).** Global-scale phylogenetic linguistic inference from lexical resources.
*Scientific Data* 5: 180189. https://doi.org/10.1038/sdata.2018.189
PMI sound-class scores and SVM cognate clustering over ~7,000 ASJP lists; binary characters for
character-based inference. *Relevance:* a ready-made automated cognate set for any family
(including Arawan) to run clepsydra when expert cognates are lacking; useful as a cross-check.

**Bouckaert, R., Redding, D., Trebski, A., Sheehan, O., Gray, R., Jones, K. E. & Atkinson, Q.
(2026).** Global language diversification is linked to socio-ecology and threat status. SocArXiv
preprint. https://doi.org/10.31235/osf.io/f8tr6_v2
A dated supertree of 6,635 languages integrating prior knowledge and uncertainty. *Relevance:*
a source of published (prior-heavy) root-age distributions for every family; comparisons with
clepsydra's uncalibrated estimates are useful but not independent (the supertree's dates are
largely priors).

### 1.5 Data infrastructure and conceptual reviews

**Forkel, R., List, J.-M., Greenhill, S. J., Rzymski, C., Bank, S., Cysouw, M., Hammarström, H.,
Haspelmath, M., Kaiping, G. A. & Gray, R. D. (2018).** Cross-Linguistic Data Formats, advancing data
sharing and re-use in comparative linguistics. *Scientific Data* 5: 180205.
https://doi.org/10.1038/sdata.2018.205
**List, J.-M., Forkel, R., Greenhill, S. J., Rzymski, C., Englisch, J. & Gray, R. D. (2022).**
Lexibank, a public repository of standardized wordlists with computed phonological and lexical
features. *Scientific Data* 9: 316. https://doi.org/10.1038/s41597-022-01432-0
*Relevance:* clepsydra should read CLDF Wordlist (FormTable + CognateTable + ParameterTable with
Concepticon IDs) directly, mapping concepts via Concepticon to its per-concept rate priors.

**List, J.-M. (2016).** Beyond cognacy: historical relations between words and their implication for
phylogenetic reconstruction. *Journal of Language Evolution* 1(2): 119–136.
https://doi.org/10.1093/jole/lzw006 — partial cognacy, root vs. word cognacy, and how coding choice
changes trees. *Relevance:* clepsydra's input must declare cognacy type (root/loose vs
word/strict); strict coding ⇒ higher apparent rates ⇒ different priors.

**Jacques, G. & List, J.-M. (2019).** Save the trees: why we need tree models in linguistic
reconstruction (and when we should apply them). *Journal of Historical Linguistics* 9(1): 128–167.
https://doi.org/10.1075/jhl.17008.mat
Incomplete lineage sorting explains much "non-tree" signal; trees remain necessary. *Relevance:*
ILS in a dialect continuum = cognate states persisting polymorphically across splits; the
novel-state model assumes immediate fixation — a known misspecification near rapid radiations.

**Jacques, G. & Pellard, T. (2021).** Phylogenies based on lexical innovations refute the Rung
hypothesis. *Diachronica* 38(1): 1–24. https://doi.org/10.1075/dia.19058.jac — lexical
innovation-based subgrouping of Sino-Tibetan.
**Pellard, T., Ryder, R. J. & Jacques, G. (2026; online 2024?).** Family tree model. In *The Wiley
Blackwell Companion to Diachronic and Historical Linguistics*, 1–33.
https://doi.org/10.1002/9781119898023.wbcdl056 — review of tree vs wave models, sources of
non-tree signal (borrowing, mixing, parallel innovation, ILS), inference of topology and dates, and
validity assessment. (The "Pellard, Ryo & Jacques 2024" item in the task brief could not be found;
this chapter — Pellard, *Ryder* & Jacques — is the likely intended reference. **[Flag]**)

**Bowern, C. (2018).** Computational phylogenetics. *Annual Review of Linguistics* 4: 281–296.
https://doi.org/10.1146/annurev-linguistics-011516-034142 — review of methods, data and critiques.
**Bowern, C. (2017).** *The Indo-European Controversy* and Bayesian phylogenetic methods
(review article). *Diachronica* 34(3): 421–436. https://doi.org/10.1075/dia.34.3.04bow — a
measured rebuttal of Pereltsvaig & Lewis.
**Pereltsvaig, A. & Lewis, M. W. (2015).** *The Indo-European Controversy: Facts and Fallacies in
Historical Linguistics.* Cambridge University Press. https://doi.org/10.1017/CBO9781107294332 —
critique of Bouckaert et al. 2012 (dates, geography, Romance anomalies). *Relevance:* the critique
list doubles as a list of diagnostics clepsydra should be able to produce (e.g. are attested
sub-family dates reproduced when held out?).
**Dunn, M. (2015).** Language phylogenies. In C. Bowern & B. Evans (eds), *The Routledge Handbook
of Historical Linguistics*, 190–211 (Crossref gives 208–229 for the e-book chapter).
https://doi.org/10.4324/9781315794013-19 **[page range inconsistent between sources]**.
**Gray, R. D., Atkinson, Q. D. & Greenhill, S. J. (2011).** Language evolution and human history:
what a difference a date makes. *Phil. Trans. R. Soc. B* 366: 1090–1100.
https://doi.org/10.1098/rstb.2010.0378 — programmatic review.
**Greenhill, S. J., Currie, T. E. & Gray, R. D. (2009).** Does horizontal transmission invalidate
cultural phylogenies? *Proc. R. Soc. B* 276: 2299–2306. https://doi.org/10.1098/rspb.2008.1944
Simulated borrowing: topologies robust to realistic borrowing levels; "inferences about divergence
dates are slightly less robust and show a tendency to underestimate dates." *Relevance:* the
canonical quantitative statement of "borrowing → too young"; clepsydra's simulation suite should
reproduce it and show how much of the bias the contamination mixture removes.
**List, J.-M., Nelson-Sathi, S., Geisler, H. & Martin, W. (2014).** Networks of lexical borrowing and
lateral gene transfer in language and genome evolution. *BioEssays* 36(2): 141–150.
https://doi.org/10.1002/bies.201300096
**Ringe, D., Warnow, T. & Taylor, A. (2002).** Indo-European and computational cladistics.
*Transactions of the Philological Society* 100(1): 59–129. https://doi.org/10.1111/1467-968X.00091
— dataset used by Nicholls & Gray / Ryder & Nicholls; lexical + phonological + morphological
characters. **Nakhleh, L., Ringe, D. & Warnow, T. (2005).** Perfect phylogenetic networks.
*Language* 81(2): 382–420. https://doi.org/10.1353/lan.2005.0078

### 1.6 Family case studies (methods relevant to clepsydra)

(Dates and intervals are collected in Part 2.)

**Gray, R. D. & Jordan, F. M. (2000).** Language trees support the express-train sequence of
Austronesian expansion. *Nature* 405: 1052–1055. https://doi.org/10.1038/35016575 — parsimony
topology test, no dating.

**Gray, R. D., Drummond, A. J. & Greenhill, S. J. (2009).** Language phylogenies reveal expansion
pulses and pauses in Pacific settlement. *Science* 323: 479–483.
https://doi.org/10.1126/science.1166858
**Greenhill, S. J., Drummond, A. J. & Gray, R. D. (2010).** How accurate and robust are the
phylogenetic estimates of Austronesian language relationships? *PLoS ONE* 5(3): e9573.
https://doi.org/10.1371/journal.pone.0009573 — high congruence with traditional subgroups;
Austronesian ~5,200 BP robust when only the most robust calibrations are used.
*Relevance:* calibration-subset robustness is a validation clepsydra should run automatically
(drop-one, keep-most-robust).
**Greenhill, S. J., Blust, R. & Gray, R. D. (2008).** The Austronesian Basic Vocabulary Database.
*Evolutionary Bioinformatics* 4: 271–283. https://doi.org/10.4137/EBO.S893

**Bowern, C. & Atkinson, Q. (2012).** Computational phylogenetics and the internal structure of
Pama-Nyungan. *Language* 88(4): 817–845. https://doi.org/10.1353/lan.2012.0081 — Bayesian subgrouping (four
major divisions); no dating found (**not confirmed from full text**).
**Bouckaert, R. R., Bowern, C. & Atkinson, Q. D. (2018).** The origin and expansion of
Pama–Nyungan languages across Australia. *Nature Ecology & Evolution* 2: 741–749.
https://doi.org/10.1038/s41559-018-0489-3 — see Part 2 for the prior/calibration strategy in a
family with almost no historical calibrations.

**Grollemund, R., Branford, S., Bostoen, K., Meade, A., Venditti, C. & Pagel, M. (2015).** Bantu
expansion shows that habitat alters the route and pace of human dispersals. *PNAS* 112(43):
13296–13301. https://doi.org/10.1073/pnas.1503793112 — 424 languages, 100 meanings, covarion,
variable-rates clock; four node calibrations.
**Koile, E., Greenhill, S. J., Blasi, D. E., Bouckaert, R. & Gray, R. D. (2022).** Phylogeographic
analysis of the Bantu language expansion supports a rainforest route. *PNAS* 119(32):
e2112853119. https://doi.org/10.1073/pnas.2112853119

**Kitchen, A., Ehret, C., Assefa, S. & Mulligan, C. J. (2009).** Bayesian phylogenetic analysis of
Semitic languages identifies an Early Bronze Age origin of Semitic in the Near East. *Proc. R. Soc.
B* 276: 2703–2710. https://doi.org/10.1098/rspb.2009.0408
**Lee, S. & Hasegawa, T. (2011).** Bayesian phylogenetic analysis supports an agricultural origin of
Japonic languages. *Proc. R. Soc. B* 278: 3662–3669. https://doi.org/10.1098/rspb.2011.0518
**Honkola, T., Vesakoski, O., Korhonen, K., Lehtinen, J., Syrjänen, K. & Wahlberg, N. (2013).**
Cultural and climatic changes shape the evolutionary history of the Uralic languages. *J. Evol.
Biol.* 26(6): 1244–1253. https://doi.org/10.1111/jeb.12107
**Savelyev, A. & Robbeets, M. (2020).** Bayesian phylolinguistics infers the internal structure and
the time-depth of the Turkic language family. *J. Lang. Evol.* 5(1): 39–53.
https://doi.org/10.1093/jole/lzz010
**Kolipakam, V., Jordan, F. M., Dunn, M., Greenhill, S. J., Bouckaert, R., Gray, R. D. & Verkerk, A.
(2018).** A Bayesian phylogenetic study of the Dravidian language family. *R. Soc. Open Sci.* 5:
171504. https://doi.org/10.1098/rsos.171504
**Sagart, L., Jacques, G., Lai, Y., Ryder, R. J., Thouzeau, V., Greenhill, S. J. & List, J.-M.
(2019).** Dated language phylogenies shed light on the ancestry of Sino-Tibetan. *PNAS* 116(21):
10317–10322. https://doi.org/10.1073/pnas.1817972116
**Zhang, M., Yan, S., Pan, W. & Jin, L. (2019).** Phylogenetic evidence for Sino-Tibetan origin in
northern China in the Late Neolithic. *Nature* 569: 112–115.
https://doi.org/10.1038/s41586-019-1153-z
**Robbeets, M., Bouckaert, R., Conte, M., et al. (2021).** Triangulation supports agricultural spread
of the Transeurasian languages. *Nature* 599: 616–621. https://doi.org/10.1038/s41586-021-04108-8
**Hruschka, D. J., Branford, S., Smith, E. D., Wilkins, J., Meade, A., Pagel, M. & Bhattacharya, T.
(2015).** Detecting regular sound changes in linguistics as events of concerted evolution. *Current
Biology* 25(1): 1–9. https://doi.org/10.1016/j.cub.2014.10.064 — Turkic; concerted (correlated)
change as a model component.
**Greenhill, S. J., Haynie, H. J., Ross, R. M., Chira, A. M., List, J.-M., Campbell, L., Botero, C. A.
& Gray, R. D. (2023).** A recent northern origin for the Uto-Aztecan family. *Language* 99(1):
81–107. https://doi.org/10.1353/lan.2023.0006 (Crossref also lists a pre-publication DOI
10.1353/lan.0.0276). 34 UA + 2 Kiowa-Tanoan varieties, 121 meanings, 1,125 cognate sets; covarion +
relaxed clock best (BF 36); root 4,095 BP; homeland in southern California. *Relevance:* a clear,
well-documented calibration set (offset lognormals from culture-vocabulary arguments — e.g.
mesquite terms for Northern UA) that is a model for families with no written history.
**Auderset, S., Greenhill, S. J., DiCanio, C. T. & Campbell, E. W. (2023).** Subgrouping in a
'dialect continuum': a Bayesian phylogenetic analysis of the Mixtecan language family. *J. Lang.
Evol.* 8(1): 33–63. https://doi.org/10.1093/jole/lzad004
**Buckley, C. D., Kopp, E., Pellard, T., Ryder, R. J. & Jacques, G. (2025).** Contrasting modes of
cultural evolution: Kra-Dai languages and weaving technologies. *Evolutionary Human Sciences* 7.
https://doi.org/10.1017/ehs.2025.10008
**King, B., Greenhill, S. J., Reid, L. A., Ross, M., Walworth, M. & Gray, R. (2023).** Bayesian
phylogenetic analysis of Philippine languages supports a rapid migration of Malayo-Polynesian
languages. SocArXiv. https://doi.org/10.31235/osf.io/re8m6 — pervasive geographic (diffusion)
signal; language shift in "Negrito" groups leaves no discernible effect on basic vocabulary.

**South American families (Amazonia).**
**Walker, R. S. & Ribeiro, L. A. (2011).** Bayesian phylogeography of the Arawak expansion in lowland
South America. *Proc. R. Soc. B* 278(1718): 2562–2567. https://doi.org/10.1098/rspb.2010.2579
(the task brief's "Biology Letters" is incorrect). 60 Arawak languages, 100-item lists, 694 binary
sites; BEAST 1.6.1, binary covarion, Yule, strict clock, relaxed random-walk phylogeography;
**no calibrations, no absolute ages**; Western Amazonian origin favoured. *Relevance:* the typical
Amazonian situation — relative chronology only — which clepsydra's cross-family prior is meant to
turn into calibrated absolute ages with honest intervals.
**Michael, L., de Carvalho, F., Chacon, T., Rybka, K., Sabogal, A., Chousou-Polydouri, N. &
Kaiping, G. (2023; online 2022).** Deriving calibrations for Arawakan using archaeological evidence.
*Interface Focus* 13(1): 20220049. https://doi.org/10.1098/rsfs.2022.0049 — 12 candidate
calibrations, no dated tree (**content summarized from secondary reading**).
*Relevance:* a ready list of candidate calibrations for an Arawakan test of clepsydra.
**Michael, L., Chousou-Polydouri, N., Bartolomei, K., Donnelly, E., Meira, S., Wauters, V. &
O'Hagan, Z. (2015).** A Bayesian phylogenetic classification of Tupí-Guaraní. *LIAMES* 15(2):
193–221. https://doi.org/10.20396/liames.v15i2.8642301 — 30 TG languages + Awetí and Mawé, 543
meanings, MrBayes restriction-site; not time-calibrated.
**Ferraz Gerardi, F. & Reichert, S. (2021).** The Tupí-Guaraní language family: a phylogenetic
classification. *Diachronica* 38(2): 151–188. https://doi.org/10.1075/dia.18032.fer
**Ferraz Gerardi, F., Tresoldi, T., Coelho Aragon, C., Reichert, S., Gregorio de Souza, J. & Noelli, F. S. (2023).** Lexical phylogenetics of the Tupí-Guaraní family:
language, archaeology, and the problem of chronology. *PLOS ONE* 18(6): e0272226.
https://doi.org/10.1371/journal.pone.0272226 — 40 doculects,
415 concepts; BEAST 2.6.6 binary covarion, UCLN, birth–death; **no archaeological calibrations**
(uniform root prior + tip dates; Tupinambá and Old Guaraní judged unsuitable as node calibrations);
TG breakup 1,762 BP (847–2,740); marginal likelihood favoured a strict clock but a relaxed clock was
chosen on expert grounds. *Relevance:* the closest methodological precedent for clepsydra's
Amazonian use case (and an in-house dataset): rerun with clepsydra's T1, the cross-family clock
prior and strict/relaxed averaging, and compare intervals.
**Ferraz Gerardi, F., Wientzek, T., Roksandic, I., Gregorio de Souza, J. & Orphão de Carvalho, F.
(2025).** A phylogenetic classification of the Je language family. *Open Research Europe* 5: 29 (v3).
https://doi.org/10.12688/openreseurope.19346.3 — 14 languages, 303 concepts, BEAST 2, uniform root
prior 500–6,000 BP plus three node calibrations; root interval reaches the prior bound (per a
reviewer's report; **not checked in the figure**). *Relevance:* an explicit example of a prior-bounded root — the failure mode
idea 16 is designed to make visible.
**Galucio, A. V., Meira, S., Birchall, J., Moore, D., Gabas Júnior, N., Drude, S., Storto, L.,
Picanço, G. & Rodrigues, C. R. (2015).** Genealogical relations and lexical distances within the
Tupian linguistic family. *Bol. Mus. Para. Emílio Goeldi Ciênc. Hum.* 10(2): 229–274.
https://doi.org/10.1590/1981-81222015000200004 — distance-based, undated.
**Chacon, T. C. (2014).** A revised proposal of Proto-Tukanoan consonants and Tukanoan family
classification. *IJAL* 80(3): 275–322. https://doi.org/10.1086/676393
— comparative-method classification; no Bayesian dating of Tukanoan found.
**Ferraz Gerardi, F., Aragon, C. C. & Reichert, S. (2022).** KAHD: Katukinan-Arawan-Harakmbut
Database (pre-release). *Journal of Open Humanities Data* 8. https://doi.org/10.5334/johd.80
— **directly relevant for Arawan**: an open lexical/cognate-coded database of Katukinan, Arawan and
Harakmbut lists that can feed clepsydra. (Coding and coverage not checked.)
**Birchall, J., Dunn, M. & Greenhill, S. J. (2016).** A combined comparative and phylogenetic analysis
of the Chapacuran language family. *IJAL* 82(3): 255–284. https://doi.org/10.1086/687383 —
3 extant + 9 historically attested lects, 285 cognate sets; trees calibrated with ethnohistorical
information and data-collection dates (root age not retrieved). *Relevance:* uses colonial word lists
as dated tips — the same resource Arawan has.
**Zhang, H., Ji, T., Pagel, M. & Mace, R. (2020).** Dated phylogeny suggests early Neolithic origin of
Sino-Tibetan languages. *Scientific Reports* 10: 20792. https://doi.org/10.1038/s41598-020-77404-4
**Sicoli, M. A. & Holton, G. (2014).** Linguistic phylogenies support back-migration from Beringia to
Asia. *PLoS ONE* 9(3): e91722. https://doi.org/10.1371/journal.pone.0091722 — typological, undated.
(Further Amazonian items in Part 2 and the verification notes.)

---

## Part 2. Published node-age estimates (benchmarks)

Conventions: "BP" = years before present as reported (ASJP dates are years before 2000). Interval =
95% HPD unless noted. **Source status:** "paper" = read in the paper/abstract/SI; "sec." = taken from
a secondary source or a figure/tree file and still to be checked; "unv." = unverified.

### 2.1 Indo-European

| Node | Point estimate | Interval | Method | Citation | Status |
|---|---|---|---|---|---|
| PIE root (Hittite split) | ~8,700 BP | 7,800–9,800 (95% CI) | binary restriction-site + Γ, MrBayes; r8s penalized likelihood, 14 constraints; 87 langs, 2,449 cognate sets | Gray & Atkinson 2003 | paper |
| PIE root (exp. branch-length prior) | — | 7,100–9,200 | same | Gray & Atkinson 2003 | paper |
| PIE root | centred 8,400 BP | 7,100–9,800 | stochastic Dollo + missing data + catastrophes (TraitLab), Ringe et al. data, 24 langs | Ryder & Nicholls 2011 | paper (arXiv abstract) |
| PIE root | median not given in text | 7,116–10,410 | BEAST, covarion, phylogeography, 103 langs | Bouckaert et al. 2012 | paper (HPD); the often-quoted ~8,700 median unv. |
| PIE root (corrected 2013 analysis, as replicated by Chang) | 7,580 | 5,970–9,350 | as above | Bouckaert et al. 2013 erratum, via Chang et al. 2015 | sec. |
| PIE root, ancestry-constrained (A1 narrow, 52 langs) | 5,950 | 4,870–7,190 | BEAST, RSC, UCLN, skyline coalescent, 8 ancestry constraints | Chang et al. 2015 | paper |
| PIE root, A1 corrected for "advergence" | 6,480 | 5,370–7,790 | same | Chang et al. 2015 | paper |
| PIE root, A2 medium (82 langs, 143 meanings) | 6,410 | 5,190–7,910 | same | Chang et al. 2015 | paper |
| PIE root, A3 broad (94 langs, 197 meanings) | 5,930 | 4,860–7,250 | same | Chang et al. 2015 | paper |
| PIE root, C1 with constraints vs C3 without | 6,470 vs 7,670 | 4,920–8,700 vs 6,190–9,460 | same; isolates the effect of constraints (~+1,200 y) | Chang et al. 2015 | paper |
| PIE root, replication of Bouckaert 2013 data (B1) | 7,750 | 6,180–9,720 | same | Chang et al. 2015 | paper |
| PIE root | ~8,120 BP | 6,740–9,610 | IE-CoR, 161 langs (52 ancient), 170 meanings; BEAST 2 sampled ancestors (FBD), covarion, relaxed clock | Heggarty et al. 2023 | paper (root); sub-nodes sec. |
| Indo-Iranian split | ~6,980 BP | 5,650–8,400 | same | Heggarty et al. 2023 | sec. |
| Indic/Iranian split | ~5,520 BP | 4,540–6,800 | same | Heggarty et al. 2023 | sec. |
| PIE root (Anatolian split), 13 proto/ancient 110-item lists | mean 3686 BCE (~5,650 BP) | 4139–3450 BCE | MrBayes 3.2.6 + NJ + MP; Moscow protocol | Kassian et al. 2021 | paper |
| Tocharian split | mean 3011 BCE | 3727–2262 BCE | same | Kassian et al. 2021 | paper |
| Inner-IE radiation (4-way) | mean 2717 BCE | 3357–2162 BCE | same | Kassian et al. 2021 | paper |
| Balto-Slavic–Indo-Iranian | mean 2241 BCE | 2723–1790 BCE | same | Kassian et al. 2021 | paper |
| Indo-Iranian | mean 1740 BCE | 2044–1458 BCE | same | Kassian et al. 2021 | paper |
| Balto-Slavic | mean 1250 BCE | 1686–855 BCE | same | Kassian et al. 2021 | paper |
| PIE root, uniform / FBD / coalescent priors, dataset B1 | 6,789 / 7,512 / 9,821 | FBD: 6,244–8,766 | MrBayes, five IE datasets; medians across all 5,935–9,821 | Rama 2018 | sec. (arXiv preprint numbers) |
| IE (excl. Anatolian, Tocharian) | 4,348 BP | — (calibration 5,500; −26%) | ASJP LDND, 40 items | Holman et al. 2011 | paper |
| IE root | 4,826 | — | Generalized Bayesian Dating (automated cognates, GBD regression) | Rama & Wichmann 2020 | sec. |
| IE root (Levenshtein, 50 langs) | ~5,800–5,900 BP (read from figure) | — | T = −ε ln(1 − γD), ε = 1750, γ = 1.09, two calibration pairs, UPGMA | Serva & Petroni 2008 | sec. (figure) |

### 2.2 Asia, Near East, Eurasia

| Family / node | Point | Interval | Method | Citation | Status |
|---|---|---|---|---|---|
| Uralic root | median ≈5,113 BP (mean ≈5,266) | ≈3,330–7,503 | BEAST, CTMC+Γ, relaxed clock; 17 langs, Ura100 | Honkola et al. 2013 | sec. (parsed from published MCC tree) |
| Finno-Ugric | ~3,900 BP | 2,695–5,371 | same | Honkola et al. 2013 | sec. |
| Uralic | 3,178 | — | ASJP | Holman et al. 2011 | paper |
| Turkic root (Bulgharic vs Common Turkic) | median 2,066 BP | 1,517–2,755 | BEAST 2.5.1, covarion, strict clock, FBD; 32 varieties, 254 concepts, tip dates Old Turkic, Cuman | Savelyev & Robbeets 2020 | paper |
| Common Turkic | 474 CE | 7 BCE–809 CE | same | Savelyev & Robbeets 2020 | paper |
| Turkic | 3,404 (calibration 2,500; +27%) | — | ASJP | Holman et al. 2011 | paper |
| Semitic root | ~5,750 BP | 4,400–7,400 | BEAST; 25 langs; Akkadian, Ugaritic, Hebrew, Aramaic, Ge'ez calibrations | Kitchen et al. 2009 | root paper; HPD sec. |
| Japonic root | median 2,182 BP (mean 2,398) | 1,239–4,190 | BEAST, covarion, relaxed; 59 lects, 210 meanings | Lee & Hasegawa 2011 | sec. (thesis restating paper) |
| Japonic | 1,564 | — | ASJP | Holman et al. 2011 | paper |
| Sino-Tibetan root | mean 7,184 BP | 5,093–9,568 | BEAST, covarion relaxed; 50 langs, 180 concepts | Sagart et al. 2019 | paper |
| Sino-Tibetan root | mean 5,871 BP | 4,183–7,817 | BEAST; 109 langs, 949 root-meanings | Zhang et al. 2019 | paper (SI) |
| Tibeto-Burman | 4,665 BP | 3,716–5,688 | same | Zhang et al. 2019 | paper (SI) |
| Sino-Tibetan root | mean 7,983 BP | 4,778–11,285 | BEAST 2; 131 langs, 110 items | Zhang, Ji, Pagel & Mace 2020, *Sci. Rep.* 10: 20792, https://doi.org/10.1038/s41598-020-77404-4 | paper |
| Sino-Tibetan | 5,261 | — | ASJP | Holman et al. 2011 | paper |
| Dravidian root | mean 4,650, median 4,433 | ~3,000–6,500 (across models) | BEAST 2 relaxed covarion; 20 langs, 100 concepts | Kolipakam et al. 2018 | paper |
| Dravidian | 2,055 | — | ASJP | Holman et al. 2011 | paper |
| Transeurasian root | 9,181 BP | 5,595–12,793 | BEAST 2; 98 langs, 254 concepts | Robbeets et al. 2021 | paper |
| Proto-Altaic | 6,811 | 4,404–10,166 | same | Robbeets et al. 2021 | paper |
| Japano-Koreanic | 5,458 | 3,335–8,024 | same | Robbeets et al. 2021 | paper |
| Mongolo-Tungusic | 4,491 | 2,599–6,373 | same | Robbeets et al. 2021 | paper |
| "Eurasiatic" | 14.45 ± 1.75 kya | 11.72–18.38 kya | ultraconserved words, LWED proto-forms | Pagel et al. 2013 | paper; strongly disputed (Heggarty 2013) |
| Afro-Asiatic | 6,016 | — | ASJP | Holman et al. 2011 | paper |
| Austro-Asiatic | 3,635 (Mon-Khmer 3,406; Munda 2,574) | — | ASJP | Holman et al. 2011 | paper; no peer-reviewed Bayesian date found |

### 2.3 Pacific, Australia, Africa, North America

| Family / node | Point | Interval | Method | Citation | Status |
|---|---|---|---|---|---|
| Austronesian root | mean 5,230 BP | 4,750–5,800 | 400 langs, 210 meanings, covarion; r8s PL with 10 node + 4 tip calibrations | Gray, Drummond & Greenhill 2009 | paper |
| Austronesian root, 5 robust calibrations only | mean 5,117 BP | 4,660–5,680 | same trees | Greenhill, Drummond & Gray 2010 | paper |
| Austronesian | 3,633 | — | ASJP | Holman et al. 2011 | paper |
| Austronesian | 5,647 | 3,314–8,831 | GBD | Rama & Wichmann 2020 | sec. |
| Austronesian (Levenshtein) | ~4,000 (from figure) | — | Petroni & Serva | Petroni & Serva 2008 | sec. (figure) |
| Pama-Nyungan root | mean 5,671 | 4,455–6,966 | BEAST, BD-skyline, UCLN; 306 langs, 200 meanings; **single archaeological calibration (Wati)** | Bouckaert, Bowern & Atkinson 2018 | paper |
| Pama-Nyungan | 3,905 | 2,676–5,403 | GBD | Rama & Wichmann 2020 | sec. |
| Pama-Nyungan | 4,295 (calibration 4,500) | — | ASJP | Holman et al. 2011 | paper |
| Bantu in-group | 4,846.5 ± 138.1 | — | BayesPhylogenies covarion; 424 langs, 100 meanings | Grollemund et al. 2015 | sec. (PMC summary) |
| Bantoid+Bantu root | 6,929.7 ± 418.6 | — | same | Grollemund et al. 2015 | sec. |
| Bantu root (lexical only / combined) | 4,940 / 5,110 | 4,500–5,400 / 4,640–5,770 | BEAST 2, CTMC+Γ, relaxed, Grollemund calibrations as lognormals | Koile et al. 2022 | paper |
| Bantu rainforest-crossing node | 4,420 | 4,040–5,000 | same | Koile et al. 2022 | paper |
| Uto-Aztecan root | mean 4,095 BP | 3,258–5,025 (range over models 3,825–5,395) | BEAST 2.4.1, covarion + relaxed; 34 UA + 2 Kiowa-Tanoan, 121 meanings; 5 node calibrations + Classical Nahuatl tip | Greenhill et al. 2023 | paper/preprint |

### 2.4 South America (the regime clepsydra most needs)

| Family / node | Point | Interval | Method | Citation | Status |
|---|---|---|---|---|---|
| Tupí-Guaraní breakup | 1,762 BP | 847–2,740 | BEAST 2.6.6 binary covarion, UCLN, BD; 40 doculects, 415 concepts; **no archaeological calibration** (uniform root prior + tip dates) | Ferraz Gerardi et al. 2023 | paper |
| Awetí + TG | 2,603 | 1,404–4,037 | same | Ferraz Gerardi et al. 2023 | paper |
| Mawé/Awetí–TG | 3,312 | 2,500–4,620 | same | Ferraz Gerardi et al. 2023 | paper |
| Tupi-Guarani coastal | 2,043 (calibration 1,750) | — | ASJP | Holman et al. 2011 | paper |
| Tupi (whole family) | 3,585 | — | ASJP | Holman et al. 2011 | paper; no Bayesian Tupian-wide date found |
| Arawakan (Maipuran) | 4,134 | — | ASJP (no family-wide date given) | Holman et al. 2011 | paper |
| Arawak | no absolute date | — | BEAST, binary covarion, strict clock, uncalibrated | Walker & Ribeiro 2011 | paper |
| **Arauan (Arawan)** | **1,764** | — | ASJP | Holman et al. 2011 | paper; **no other published date found** |
| Tucanoan (Eastern / Western) | 2,699 (1,241 / 2,156) | — | ASJP | Holman et al. 2011 | paper; no Bayesian date found |
| Chapacura-Wanham | 1,931 | — | ASJP | Holman et al. 2011 | paper (Birchall et al. 2016 Bayesian root unv.) |
| Panoan / Tacanan | 1,853 / 1,590 | — | ASJP | Holman et al. 2011 | paper |
| Carib | 2,362 | — | ASJP | Holman et al. 2011 | paper |
| Macro-Ge (Ge-Kaingang) | 7,266 (4,989, said to be exaggerated) | — | ASJP | Holman et al. 2011 | paper |
| Jê root | figure only; interval reaches the 6,000 BP prior bound | uniform 500–6,000 prior | BEAST 2; 14 langs, 303 concepts; node calibrations | Ferraz Gerardi et al. 2025, *Open Res. Europe* 5: 29, https://doi.org/10.12688/openreseurope.19346.3 | sec. |
| Quechuan / Aymaran | 1,717 / 1,057 | — | ASJP | Holman et al. 2011 | paper |
| Mayan | 2,220 (paper says too young) | — | ASJP | Holman et al. 2011 | paper |
| Oto-Manguean | 6,591 | — | ASJP | Holman et al. 2011 | paper |

**Reading the table for Arawan.** The only published absolute figure for Proto-Arawan is ASJP's 1,764
years (in-sample error ~29%, likely more in Amazonia where Holman's calibrations are sparse). The
comparative-method literature (Dixon 2004 "Proto-Arawá phonology", *Anthropological Linguistics*
46(1): 1–83 **[verify]**; Dienst 2008 "The internal classification of the Arawan languages",
*LIAMES* 8 **[verify]**) gives structure (~460 cognate sets, no demonstrated external relatives)
but no date. Tupí-Guaraní, the best-studied comparable Amazonian family, shows that even a careful
uncalibrated Bayesian analysis gives a ~3× interval. A clepsydra Arawan result should therefore be
presented with (i) the cross-family prior-predictive age, (ii) the ASJP figure for comparison, and
(iii) the prior→posterior contraction. That is how it is now reported (`arawan-result.md`): root
2,235 BP (563–5,171) on the empirical clock, `rate_contraction` 0.32–0.38, with ASJP and classical
glottochronology alongside (`arawan-glottochronology.md`).

**Calibrations used for Uto-Aztecan (template for prior-light families).** Offset lognormals: Northern
UA median 3,230 (2,960–3,970; mesquite vocabulary), Takic 2,980 (2,690–3,730), Tübatulabal 2,450
(2,390–2,900), Numic 1,350 (1,290–1,800), Nahuan 1,740 (1,570–2,410); Classical Nahuatl tip 445 BP.
**Pama-Nyungan:** single Wati calibration, Gamma(α = 2, β = 359) offset 3,000, 95% between 3 and 5 ka;
the authors note the whole chronology scales roughly linearly with it.

**Other ASJP figures (Holman et al. 2011, Table 1 calibration vs ASJP estimate):** Benue-Congo
6,500 / 4,940; IE 5,500 / 4,348; Pama-Nyungan 4,500 / 4,295; Tupi-Guarani coastal 1,750 / 2,043;
Cholan 1,600 / 1,148; Mongolic 750 / 2,267. Formula t = (log s − log s₀)/(2 log r), s = 1 − LDND,
s₀ = 0.92, r = 0.72 per millennium (if all error is attributed to similarities instead: s₀ = 0.62,
r = 0.79).

**Rama & Wichmann 2020 accuracy (leave-one-out on 30 calibrations):** MAE 770 y, RMSE 936 y; baseline
"number of languages only" MAE 822 y. In their head-to-head table: MAE GBD 695 / ASJP 606 /
glottochronology 728; RMSE 858 / 823 / 1,039; r 0.76 / 0.78 / 0.61 (differences not significant).
Worst: Benue-Congo (6,500 vs GBD 4,496), Mongolic (750 vs 2,683).

---

## Part 3. Validation evidence: how accurate are lexical dates?

What is actually known, with the caveats that matter for a "well-calibrated posterior" claim.

1. **Classical glottochronology.** Lees (1953) derived r from 13 documented pairs; Bergsland & Vogt
   (1962) showed lineage deviations large enough (Icelandic) to invalidate a universal constant.
   No modern cross-validated error rate for Swadesh's formula exists in the sources read; Rama &
   Wichmann (2020) found automated glottochronology the worst of three methods
   (leave-one-out MAE 728 y vs 606 y for ASJP and 695 y for their Bayesian method;
   r = 0.61 vs 0.78 vs 0.76; differences not significant by Wilcoxon).
2. **ASJP (Holman et al. 2011).** Mean absolute discrepancy 29% of the estimated date over 52
   calibration points; 5/52 off by >50%, 1/52 by >100%; correlation of log similarity with date
   −0.84. **This is in-sample** (same 52 points fit s0 = 0.92 and r = 0.72 per millennium), so it
   is an optimistic figure; Rama & Wichmann's leave-one-out comparison is the better number.
   Known big misses: Benue-Congo (calibration 6,500 vs ASJP 4,940), Mongolic (750 vs 2,267).
3. **Generalized Bayesian Dating (Rama & Wichmann 2020).** Leave-one-out over 30 calibrations:
   MAE 770 y, RMSE 936 y; a trivial baseline using only the number of languages in the family has
   MAE 822 y. *Implication for clepsydra:* a no-calibration method that does not clearly beat the
   "number of languages" baseline is not yet informative; clepsydra's cross-family validation
   should report this baseline.
4. **Bayesian with calibrations.** Austronesian: root essentially unchanged when restricted to the 5
   least controversial calibrations (Greenhill et al. 2010; 5,117 vs 5,230 BP means). IE: Nicholls &
   Gray report leave-one-calibration-out "congruence" without numbers; Ryder & Nicholls (2011) reject
   3/30 historical constraints by Bayes factor; a leave-one-out test cited by Rama & Wichmann (Table 7)
   found 5 of 8 ancient IE languages with inferred medians outside their known ranges
   (**source of that test unconfirmed; probably Chang et al. 2015**).
5. **Simulation studies.** Atkinson et al. (2005): with true root 8,680, borrowing at 20%
   produced slight underestimation; at 100% about 30% underestimation. Greenhill, Currie & Gray
   (2009): topologies robust, dates "tend to underestimate". Neureiter et al. (2022): ignoring
   contact biases tree height, rates and topology. Rönchen, Billing & Wiklund (2023): the
   Chang et al. and Bouckaert et al. pipelines both fail to recover the simulated IE root age.
   King (2026a): BEAST's default cognate-count rate reweighting inflates tree height in SBC.
6. **Uncalibrated Amazonian studies.** Ferraz Gerardi et al. (2023) Tupí-Guaraní used no
   archaeological calibrations (uniform root prior + tip dates); TG breakup 1,762 BP (847–2,740)
   — a factor-of-three interval, i.e. what honest no-calibration uncertainty looks like.
   Walker & Ribeiro (2011) reported no absolute ages for Arawak at all.

**Systematic bias directions reported in the literature** (and what clepsydra can do):

| Cause | Direction | Evidence | Clepsydra lever |
|---|---|---|---|
| Unrecognized borrowing (between sister lineages) | too young | Embleton 1986; Atkinson et al. 2005; Greenhill et al. 2009; Neureiter et al. 2022 | contamination mixture; concept-specific borrowability prior |
| Borrowing from outside the family | too old (extra apparent replacement) | logic of the model; Bowern et al. 2011 rates | separate "external loan" flag → treat as missing |
| Saturation / multiple hits not modelled (single-rate pairwise) | too young at depth | Sankoff 1973; gamma-mixture argument; Starostin's time-squared law is an empirical response | gamma or DP concept rates; report the depth beyond which the posterior is prior-dominated |
| Ancient languages placed on side branches when they are ancestors | too old | Chang et al. 2015 | ancestor (sampled-ancestor) calibrations at internal nodes |
| Partition-rate reweighting by cognate count (binary model) | too old (tree height overestimated) | King 2026a (SBC) | weight by meaning; prefer T1 |
| Punctuational change at splits not modelled | depends on calibration placement (node-dense shallow calibrations ⇒ fast rate ⇒ deep nodes too young) | Atkinson et al. 2008 | per-node burst parameter |
| Cognate coding granularity (root vs word cognacy) | strict coding ⇒ more change ⇒ older if rates transferred | List 2016; Chang et al. 2015 data curation | declare coding type; separate rate priors |
| Polymorphism / synonyms coded as multiple binary presences | inflates shared states ⇒ too young; or misfit | King 2026a PPC; Vejdemo & Hörberg (synonyms ↔ rate) | treat polymorphic slots explicitly in T1 |
| Tree prior misfit (birth–death vs coalescent; sampling fraction) | root shifts up to >1 kyr | Rama 2018; Ritchie & Ho 2019 | expose node-age prior; sensitivity report |
| Rate–time non-identifiability | posterior width floored by priors | King 2026b; molecular dating theory | report prior→posterior contraction per node |

---

## Part 4. Concrete, actionable ideas for clepsydra

Each item names its sources. "T1" = strict novel-state multistate tree likelihood; "T2" = binary
gain/loss CTMC; "P0/P1" = classical/hierarchical pairwise.

1. **Cross-family hierarchical clock prior for uncalibrated families.** Fit a hierarchical model over
   all calibrated families (IE-CoR, ABVD Austronesian, Bantu, Pama-Nyungan, Uto-Aztecan, Tupí-Guaraní
   tip dates, Semitic, Japonic, Turkic, Uralic): family-level mean replacement rate
   log ρ_f ~ N(μ, τ²), per-concept rates shared across families through Concepticon IDs. For an
   uncalibrated family, use the posterior predictive for ρ_f (including τ) as its clock prior. This
   is the principled version of Rama & Wichmann's (2020) tree-length regression and of Holman et
   al.'s (2011) single global constant, and it propagates between-family rate variance (which those
   methods ignore) into the node-age posterior. *Sources:* Holman 2011; Rama & Wichmann 2020;
   Greenhill et al. 2010 (lexical rates correlate across families); Calude & Pagel 2011.
2. **Report the between-family rate variance τ honestly**, and show users how much of their
   uncalibrated node-age interval comes from τ alone (a variance decomposition: τ vs. data vs.
   branch-rate variance). If τ ≈ 0.3 on the log scale, no lexical data can give better than ~±30%
   (1σ) on an uncalibrated root — consistent with Holman's 29% in-sample error. *Sources:* Holman
   2011; Rama & Wichmann 2020; King 2026b (identifiability).
3. **Leave-one-family-out (LOFO) and leave-one-calibration-out validation as a built-in report**,
   with coverage of 50/80/95% intervals and MAE in years and in log-ratio; always include the
   "number of languages only" baseline (MAE 822 y in Rama & Wichmann) and ASJP (MAE 606 y) as the
   bars to beat. *Sources:* Rama & Wichmann 2020; Holman 2011; Greenhill et al. 2010.
4. **Use Holman et al.'s 52 calibration points (Table 1) as a shared external benchmark**, so that
   clepsydra's error can be compared directly to ASJP's on the same targets (where cognate data
   exist). Flag the known pathological points (Mongolic, Benue-Congo, Mayan "too young").
5. **Concept-rate regression prior.** log θ_c = β·x_c + ε_c, with x_c = (log frequency from Calude &
   Pagel's cross-linguistic frequencies, synonym count, imageability, part of speech, ASJP stability
   rank, WOLD borrowability). For concepts absent from the reference families, this gives a
   principled prior instead of the global mean. *Sources:* Pagel et al. 2007; Calude & Pagel 2011;
   Vejdemo & Hörberg 2016; Monaghan 2014; Holman et al. 2008; Tadmor et al. 2010.
6. **Dirichlet-process (or finite) mixture over concept rates** as an alternative to the gamma:
   slow/medium/fast clusters with cluster labels shared across families. *Source:* Greenhill et al.
   2017 (DP mixture on Austronesian).
7. **Concept-specific contamination priors.** Replace the single contamination weight with
   π_c ~ Beta(a_c, b_c) centred on WOLD borrowability for concept c, scaled by a region factor
   (Amazonia ~1–4% basic-vocabulary loans; N. Australia mean 9.4%, outliers to 48%). *Sources:*
   Tadmor et al. 2010; Haspelmath & Tadmor 2009; Bowern et al. 2011.
8. **Distinguish internal from external borrowing.** An external loan (from outside the family)
   should be treated as *missing* (no information about the internal clock), whereas an internal
   loan between sister lineages creates spurious shared states (bias young). The contamination
   mixture should have two components, and inputs should allow a `loan_source` column (IE-CoR
   records loans explicitly). *Sources:* Heggarty et al. 2023 (IE-CoR loan coding); Starostin 1999
   (exclude loans); Embleton 1986.
9. **Coexistence-weighted contamination.** Make the probability that a match between two lineages is a
   loan proportional to the time they coexisted as separate languages (and optionally geographic
   adjacency), a cheap approximation to Kelly & Nicholls (2017) and contacTrees (Neureiter et al.
   2022) that keeps the fixed-tree likelihood tractable.
10. **Per-node punctuational bursts.** Add an instantaneous replacement probability β at each
    split (each daughter independently replaces each concept with probability β at birth). In T1
    this is a simple modification of the transition at the top of each branch. Report the
    posterior fraction of change attributable to bursts and compare with the 10–33% of Atkinson
    et al. (2008). Alternative: Ryder & Nicholls's (2011) "catastrophes" on branches.
11. **Separate gain and loss multipliers per branch in T2**, with an optional covariate (log
    speaker population) under a shrinkage prior centred on zero effect. *Sources:* Bromham et al.
    2015; Greenhill et al. 2018; Wichmann & Holman 2009; Nettle 1999.
12. **Ancestor calibrations (sampled ancestors) at internal nodes.** Allow an attested ancient
    language (Latin, Vedic, Old Tupinambá, Old Guaraní, colonial Arawan word lists) to be attached
    either as a tip with a date or as the ancestor at an internal node (zero-length branch), and
    let the user compare. This is exactly the Chang et al. (2015) vs Heggarty et al. (2023) axis;
    Heath et al. 2014 and Gavryushkina et al. 2014 give the FBD/sampled-ancestor machinery.
13. **Historical word lists as dated tips with heavy missingness.** Colonial-era Amazonian lists
    (e.g. the extinct Arawá itself, last attested c. 1877, Arawan wordlists from 19th-century
    travellers) are the only time anchors for families like Arawan. They must enter with
    missing-at-random marginalization (Ryder & Nicholls 2011) and with an orthography-uncertainty
    contamination term. Ferraz Gerardi et al. (2023) and the Je study (Ferraz Gerardi et al. 2025)
    used tip dates this way.
14. **Loanword-stratum calibrations.** Support calibration type "loan_stratum": e.g. if a node's
    reconstructed vocabulary contains a loan from a dated source (Portuguese, Lingua Geral, Quechua,
    Spanish), the daughters must postdate contact. Maurits et al. (2020) show this is the most
    under-used calibration source; in Amazonia, Lingua Geral Amazônica and Portuguese loans give
    upper bounds on late splits.
15. **Calibration consistency diagnostics.** For each calibration, report (a) the posterior of the
    node with that calibration removed vs. the calibration density, (b) a Bayes-factor/posterior
    predictive p-value for conflict (Ryder & Nicholls 2011 rejected 3/30), and (c) the
    *effective* prior on calibrated nodes after interaction with the tree prior (Heath et al. 2014).
16. **Prior-only runs and prior→posterior contraction per node.** Always produce the joint prior on
    node ages implied by calibrations + tree prior + clock prior (by running the sampler with the
    likelihood off) and report a per-node contraction statistic. This addresses King 2026b's
    identifiability finding and makes "data-informed vs. prior-informed" ages visible.
17. **Posterior predictive adequacy panel from King (2026a).** Implement his eight statistics
    (max/min/mean pairwise distance; mean presence; variance of presence per language; number of
    prevalent sets >80%; singletons; variance of set prevalence), adapted for T1 (per-concept number of
    distinct states, fraction of concepts with all-identical states, state-size distribution).
    4/8 failed for binary covarion on IE — a natural headline comparison for T1.
18. **Avoid the cognate-count weighting pitfall in T2.** If T2 uses partition/meaning rates, weight
    by meaning and put an explicit Dirichlet on relative rates; include an SBC test that would catch
    the King (2026a) failure (tree height overestimated).
19. **Posterior SBC on the user's own dataset.** Simulate replicates from the fitted posterior,
    re-infer, and check rank uniformity for node ages (King 2026b). Cheap with a fixed tree and
    worth exposing as `clepsydra check --posterior-sbc`.
20. **Re-run Rönchen, Billing & Wiklund-style simulations** (loans + rate variation + near-consensus
    IE topology) and report clepsydra's bias next to the Chang/Bouckaert pipelines. This directly
    tests the claim that T1 + contamination is more robust than binary covarion.
21. **Clock-model averaging (strict vs relaxed) with a spike-and-slab on the branch-rate variance.**
    Report P(strict clock | data) — a direct, interpretable "is glottochronology adequate here?"
    test. Ferraz Gerardi et al. (2023) found the strict clock favoured by marginal likelihood for
    Tupí-Guaraní yet chose relaxed — exactly the decision this would formalise. *Source:* Xu et al.
    2026; Drummond et al. 2006.
22. **Starostin's decay law and the gamma-mixture law in `decay-laws.md`.** Add Starostin's
    C = exp(−λ C t²)-type law (constant to verify) and show that the gamma-heterogeneity law
    (1 + 2θt/α)^(−α) reproduces the same "slowing of decay" with fewer ad hoc assumptions;
    let users see how each converts a given C to t. *Sources:* Starostin 1999; Sankoff 1973;
    Novotná & Blažek 2007.
23. **Chance-resemblance / saturation floor.** For each node, compute the probability of the
    observed number of shared states under "unrelated" (state frequencies from the family's own
    lists, or ASJP-based chance similarity) and flag nodes whose posterior is prior-dominated
    because the signal is at the floor. Guards against "ultraconserved words" style over-reach.
    *Sources:* Pagel et al. 2013 vs Heggarty 2013; Holman 2011 (ASJP distances near 100% beyond
    ~6–8 kyr).
24. **Allow multifurcations / near-zero internal branches in fixed-topology mode** and report the
    posterior of internal branch lengths so rapid radiations (Kassian et al. 2021 inner IE;
    Egorov et al. 2022 Nuclear Turkic) are not forced into artefactual dated splits.
25. **Incomplete lineage sorting (polymorphism persistence) option in T1.** Allow a concept's
    ancestral slot to carry two states that sort after a split (a short "polymorphism lifetime"
    parameter). Addresses Jacques & List (2019) and the dialect-continuum results of Auderset et al.
    (2023) and King et al. (2023, Philippines).
26. **Automated-cognate input with a cognate-error rate.** Accept LexStat/SVM cognates (Rama et al.
    2018; Jäger 2018) with an extra contamination component for clustering errors; provide a
    side-by-side expert vs automatic report.
27. **CLDF-native input with Concepticon mapping**, so per-concept priors and borrowability
    covariates attach automatically. Lexibank (List et al. 2022) and KAHD (Ferraz Gerardi et al.
    2022, Katukinan-Arawan-Harakmbut) are the immediate targets for an Arawan demonstration.
28. **Coding-protocol declaration.** Input metadata states root vs word cognacy, loan treatment and
    synonym policy (Moscow protocol vs IE-CoR vs ABVD); rate priors are protocol-specific because
    coding granularity changes apparent rates (List 2016; Kassian et al. 2010; Heggarty et al. 2023).
29. **Benchmark table as a regression test.** Ship Part 2's table (with sources) as data; a
    validation command reruns clepsydra on the public datasets (IE-CoR, ABVD, Bantu, Uto-Aztecan,
    Tupí-Guaraní) and reports where its posterior medians fall relative to published intervals —
    not as truth, but to flag surprising divergences.
30. **Node-age prior choice exposed** (uniform-on-ages / birth–death with sampling fraction /
    FBD-like with dated ancestors), with an automatic sensitivity panel (Rama 2018; Ritchie & Ho
    2019). For small families like Arawan (5–7 languages), the tree prior matters more than for
    large ones.

---

## Part 5. Open controversies the tool should let users explore

1. **Indo-European root age: Anatolian (~8–9.5 kyr) vs Steppe (~5.5–6.5 kyr) vs "hybrid".** Levers:
   ancestor constraints vs tip-dated ancient languages (Chang et al. 2015 vs Heggarty et al. 2023);
   dataset curation (Dyen vs IE-CoR vs Chang); tree prior (Rama 2018); borrowing handling
   (Neureiter et al. 2022); model family (covarion vs SD vs multistate — Rönchen 2026). Clepsydra on
   IE-CoR with a toggle for "Latin/Vedic/Old Irish/… as ancestors" would let a user reproduce the
   main swing.
2. **Is there a universal lexical clock?** Constant glottochronological rate (Swadesh, Lees) vs
   strong lineage variation (Bergsland & Vogt 1962) vs relaxed clocks; how large is the
   between-family variance τ? (Holman 2011; Rama & Wichmann 2020.) The strict-vs-relaxed
   posterior probability (idea 21) and the cross-family τ (idea 2) make this inspectable.
3. **Time-dependence of rates / saturation.** Starostin's accelerating replacement vs gamma-mixture
   heterogeneity vs exponential; is there a depth limit (~6–10 kyr; Greenhill et al. 2010, 2017;
   Nichols 1992)? Users should be able to compare decay laws on their data.
4. **Ultraconserved words and deep families.** Pagel et al. 2013 vs Heggarty 2013; Pagel & Meade
   2018 on numerals. Controversy over whether slow concepts carry signal beyond chance.
5. **How much does borrowing bias dates, and in which direction?** Greenhill et al. 2009 ("slight"
   underestimation) vs Atkinson et al. 2005 simulations (up to ~30% at extreme borrowing) vs
   Neureiter et al. 2022 (systematic height errors) vs Rönchen et al. 2023 (pipelines fail).
   Hunter-gatherer borrowing (Bowern et al. 2011): not higher than agriculturalists within region.
6. **Punctuational bursts: real or artefact?** Atkinson et al. 2008 (10–33%) vs node-density
   artefact concerns; whether bursts change node ages materially.
7. **Population size and rates.** Bromham et al. 2015 (Polynesian, effect) vs Greenhill et al.
   2018 (weak/absent in Austronesian, Bantu) vs Wichmann & Holman 2009 vs Nettle 1999.
8. **Binary vs multistate encodings.** King 2026a (binary covarion misfit 4/8 statistics),
   Rönchen 2026 (covarion overfits), Häuser et al. 2026 (finite multistate overparameterized;
   stick with binary for now). The infinite-state model sits between these; clepsydra can run
   T1 and T2 on the same data and compare node ages and adequacy.
9. **Calibration practice.** Which calibrations are defensible (Maurits et al. 2020)? Narrow
   archaeological priors (Pama-Nyungan's single Wati calibration, which scales the whole tree
   ~linearly — Bouckaert et al. 2018) vs broad historical ones. The leave-one-out and
   effective-prior diagnostics expose this.
10. **Families without calibrations (Arawan, Tukanoan, Panoan, Cariban).** Is it better to (a)
    borrow a global clock with its between-family variance, (b) regress on tree length (Rama &
    Wichmann), (c) use ASJP distances (Holman), or (d) report only relative ages? Clepsydra should
    make (a) the default and always print (d) alongside.
11. **Trees vs waves / dialect continua.** Jacques & List 2019; Pellard, Ryder & Jacques
    (Wiley companion); Auderset et al. 2023 (Mixtecan); Egorov et al. 2022 (Turkic areal clusters):
    when should a fixed tree be refused?
12. **Data processing vs modelling.** Moscow-school curation (Kassian et al. 2021: reconstructed
    proto-lists, derivational drift elimination) vs model-based handling of noise (contamination,
    lateral transfer). Do cleaned inputs and robust models converge on the same dates?

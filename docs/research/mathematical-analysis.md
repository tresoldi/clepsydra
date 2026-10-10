# Lexical dating from cognate partitions: a mathematical and statistical review

*Prepared 2026-10-10 for the clepsydra project. Scope: what the problem is in
mathematical terms, how much it can deliver, and which other formulations are
worth building. It includes numerical experiments (scripts in
`docs/research/math/`, run with `python3 -I`; numbers below are their output).*

**Reference legend.** **[v]** = bibliographic details (and DOI) checked by web
search in this session. **[m]** = standard reference given from memory; the
DOI is believed correct but was not re-checked here. **[u]** = details
uncertain (year, pages or DOI); check before citing. No reference here was
invented. Where I could not recall a reference I give none.

---

## 0. The problem in one paragraph of notation

Languages $i=1..L$, meanings $c=1..C$. For each meaning the data form a
**set partition** $\pi_c$ of the languages that have a form for it (missing
cells drop out; synonyms make it a covering rather than a partition). The
family is a rooted tree $\tau$ with node ages $a_v$ (tips at 0 or at an
attested age). Branch $b$ has duration $t_b$. Under the **strict novel-state
model** (T1) meaning $c$ is replaced on branch $b$ at rate
$r_{bc}=\rho\, m_b\, \lambda_c$, where $\rho$ is the family clock (1/yr),
$m_b$ a relaxed-clock multiplier and $\lambda_c$ a relative concept rate.
Every replacement creates a class never seen before. The target is
$p(a_v \mid \pi_{1..C})$ for every $v$, marginalised over $\tau$ if needed.

Everything below follows from three facts:

1. **T1 is the infinite-alleles model** (Kimura–Crow; Ewens) run on a fixed
   genealogy that has deterministic branch lengths, instead of a coalescent.
2. **The likelihood depends on $(\rho, a)$ only through the products
   $\rho\, m_b\, t_b$.** So calendar time is identifiable only through prior
   information on $\rho$ (calibrations, or a rate prior transferred from other
   families). This is the same structure as the molecular clock, and the
   "infinite-sites" theory of Yang and Rannala applies directly.
3. **Each meaning is a current-status (case-I interval-censored) survival
   observation.** We never see when a word was replaced, only whether the
   class is still shared. This is why the shape of the retention curve is
   hard to learn nonparametrically and has to be supplied by a parametric
   family, or by a mixing distribution over rates.

---

## 1. Population-genetics connection: infinite alleles, Ewens, coalescent

### 1.1 Dictionary

| Population genetics | Lexical dating |
|---|---|
| Gene genealogy (coalescent tree) | Language family tree (fixed, or from a Yule/BD prior) |
| Locus | Meaning (concept) |
| Allele | Cognate class |
| Infinite-alleles mutation, rate $u$ | Novel-state replacement, rate $r_c$ |
| Identity by descent (IBD) = identity by state (IBS) under infinite alleles | "Cognate" = same class |
| Homozygosity $F$ | Pairwise shared-cognate probability |
| Ewens sampling formula (ESF), $\theta=4Nu$ | CRP noise component in clepsydra |
| Hoppe urn | Sequential CRP construction |
| Number of alleles $K$ (sufficient for $\theta$ under ESF) | Number of cognate classes per meaning |
| Watterson's estimator | "Lexical tree length" estimator (below) |
| Griffiths–Tavaré importance sampling over genealogies | Integration over topologies/ages |

Key references: Ewens (1972) *Theor. Popul. Biol.* 3:87–112,
doi:10.1016/0040-5809(72)90035-4 **[m]**; Karlin & McGregor (1972) "Addendum to
a paper of W. Ewens", *TPB* 3:113–116, doi:10.1016/0040-5809(72)90036-6
**[m]**; Kingman (1982) "The coalescent", *Stoch. Proc. Appl.* 13:235–248,
doi:10.1016/0304-4149(82)90011-4 **[m]**; Kingman (1978) "Random partitions in
population genetics", *Proc. R. Soc. A* 361:1–20, doi:10.1098/rspa.1978.0089
**[m]**; Hoppe (1984) "Pólya-like urns and the Ewens sampling formula", *J.
Math. Biol.* 20:91–94, doi:10.1007/BF00275863 **[m]**; Watterson (1975)
*TPB* 7:256–276, doi:10.1016/0040-5809(75)90020-9 **[m]**; Griffiths &
Tavaré (1994) "Ancestral inference in population genetics", *Stat. Sci.*
9:307–319, doi:10.1214/ss/1177010378 **[m]**; Pitman (2006) *Combinatorial
Stochastic Processes*, LNM 1875, Springer **[m]**.

### 1.2 The CRP noise component is literally "infinite alleles on a random tree"

ESF: for $n$ sampled lineages, the probability of a partition with block
sizes $n_1,\dots,n_K$ is

$$
P(\pi)=\frac{\theta^K\,\Gamma(\theta)}{\Gamma(\theta+n)}\prod_{j=1}^{K}(n_j-1)!,
$$

which is the Chinese-restaurant-process (CRP) law. ESF is exactly the
marginal partition law of infinite-alleles mutation on a **Kingman coalescent**
genealogy integrated out. So the clepsydra contamination component has a clean
generative reading: *this meaning evolved on some random genealogy that has
nothing to do with the family tree.* That makes it a principled null rather
than an ad-hoc noise model. Three consequences follow.

- **$K$ is sufficient for $\theta$ under ESF.** This is a classical result
  (Ewens 1972). So the contamination component's own parameter is learnt
  from class counts only. A cheap diagnostic follows: plot each meaning's
  $K_c$ against its expected value under the tree model (§1.4). Meanings far
  above it are contamination candidates.
- **Better noise models are available off the shelf.** The two-parameter
  Pitman–Yor (Ewens–Pitman) partition, with
  $P\propto \prod_{k=1}^{K-1}(\theta+k\alpha)\,\prod_j (1-\alpha)_{n_j-1}$, gives
  heavier tails of singletons. That is what coding noise and look-alike
  replacements produce. It is one extra parameter and still closed-form.
- **The comparison "tree vs CRP" is a test of tree-likeness per meaning.**
  The per-meaning posterior mixture weight is a likelihood ratio between
  the tree-genealogy and random-genealogy hypotheses. The ESF is
  exchangeable over languages, so it cannot know the tree. That is exactly
  what you want from a null.

### 1.3 Exact peeling for infinite alleles on a fixed tree (the p/q recursion, re-derived)

On a fixed tree the infinite-alleles likelihood is **exact and cheap**.
Importance sampling à la Griffiths–Tavaré is needed only when the genealogy
is integrated out. Here is a self-contained recursion. Clepsydra's
`docs/model-spec.md` §10 implements an equivalent p/q form from Rönchen,
Wiklund & Hammarström (2024), *Language Dynamics and Change* 14(1):1–41,
doi:10.1163/22105832-bja10032 **[v]** (the Brill appendix confirms the title
"Likelihood calculation in a multistate model of vocabulary evolution for
linguistic dating").

For node $v$ let $\Lambda(v)$ be the labels at the tips below $v$. A label is
**confined** to $v$ if all its tips lie below $v$. Define $f_v(s)$ as the
probability of the tip data below $v$ given that the state at $v$ is $s$. Here
$s\in\Lambda(v)$, or $s=\varnothing$ meaning "a class not seen below $v$". For a
child $c$ on a branch of length $t$, let $e=e^{-rt}$ and

$$
g_c = f_c(\varnothing)+\sum_{a\in\Lambda(c)\ \text{confined to }c} f_c(a)
\qquad\text{(a brand-new class is born on the branch)}.
$$

Then

$$
f_v(s)=\prod_{c}\begin{cases}
e\,f_c(s) & s\in\Lambda(c)\quad\text{(tips of class $s$ must inherit it: no change allowed)}\\[2pt]
e\,f_c(\varnothing)+(1-e)\,g_c & s\notin\Lambda(c)
\end{cases}
$$

with tips $f(\text{own label})=1$, $f(\cdot)=0$ otherwise ($f\equiv 1$ if the
cell is missing). At the root, $P(\pi)=\sum_{a\in\Lambda}f_\text{root}(a)+f_\text{root}(\varnothing)$.

**The pitfall** (I hit it while writing the experiment). The obvious version
uses $e f_c(s)+(1-e)g_c$ when $s\in\Lambda(c)$. That lets a lost class be
"re-born" further down, which is impossible under novel states. It
overcounts badly: with 7 tips and $r=10$, the total probability over all
877 partitions came to 85.8 instead of 1. After the correction above the
total is 1 to machine precision, and it matches a 200k-replicate brute-force
simulation per pattern. Pairwise marginals reproduce $e^{-r d_{ij}}$ exactly.
Clepsydra's brute-force validator (§11 of the spec) should already guard
against this. It is worth checking that its test trees include a class that is
confined to an internal clade *and* shares a label with the root state.

Cost is $O(\sum_v |\Lambda(v)|)\le O(L^2)$ per meaning. The derivative
with respect to every branch length comes from a standard pre-order (outside)
pass in the same complexity. This is the "linear-time gradient" idea of Ji et
al. (2020) *MBE* 37:3047, doi:10.1093/molbev/msaa130 **[v]**, and it is
what makes HMC/Laplace (§5) cheap.

It also helps to note that T1 is the $K\to\infty$ limit of the Mk/Jukes–Cantor
model: $P(\text{same})=\frac1K+\left(1-\frac1K\right)e^{-\frac{K}{K-1}rt}\to e^{-rt}$.
So clepsydra's reserved T3 (finite-state Mk) is a strict generalisation, with
homoplasy $1/K$.

### 1.4 Moment statistics from population genetics

**Identity by state (pairwise).** For tips $i,j$ with MRCA age $T$:
$F_{ij}^{(c)}=P(\text{same class})=\exp(-r_c d_{ij})$, with $d_{ij}=2T$ under
a clock. This is Malécot-style IBD with a deterministic genealogy. The
sample mean over meanings is a sufficient-ish statistic for $T$ when rates
are homogeneous. The moment estimator is the infinite-state Jukes–Cantor
analogue:

$$
\hat d_{ij}=-\frac{1}{\bar r}\log\hat F_{ij},\qquad
\operatorname{Var}(\hat d_{ij})\approx\frac{1-F_{ij}}{C\,\bar r^{2}\,F_{ij}}
\quad\text{(delta method, homogeneous rates)}.
$$

With heterogeneous known rates, solve $\sum_c (y_c-e^{-r_c d})\,\frac{r_c e^{-r_c d}}{e^{-r_cd}(1-e^{-r_cd})}=0$
(the score equation). Feeding $\hat d_{ij}$ to UPGMA/BIONJ/least-squares
dating gives a tree-and-date estimate in milliseconds (§5b).

**Number of classes (Watterson analogue).** Let $\rho(v)$ be the probability that
a class present at $v$ survives in at least one tip below $v$:
$\rho(\text{tip})=1$, $\rho(v)=1-\prod_{c}(1-e^{-rt_c}\rho(c))$. A class born at
distance $u$ above child $c$ survives with probability $e^{-ru}\rho(c)$.
Integrating births over the branch gives the exact identity

$$
\mathbb E[K_c]=\rho(\text{root})+\sum_{b}\bigl(1-e^{-r_c t_b}\bigr)\rho(\text{child}_b)
\;\approx\;1+r_c\,\mathcal L\quad(r_c\mathcal L\ll1),
$$

where $\mathcal L=\sum_b t_b$ is the tree length. Verified by simulation on the
7-tip test tree (`e7_watterson.py`): $r=0.5$ gives 2.8038 from the formula and
2.8095 from simulation; $r=1$ gives 4.1151 and 4.1189. The linear
approximation $1+r\mathcal L$ already fails at $r\mathcal L\approx 2$ (3.05
vs 2.80). As in population genetics, class counts estimate the **total
lexical tree length** $\rho\mathcal L$, while pairwise identities estimate
path lengths. They are cheap posterior-predictive test statistics: a tree
model that gets $F_{ij}$ right but $K_c$ wrong is signalling borrowing,
synonym handling problems or coding splits.

**Higher-order identity.** For a triple $((A,B),C)$ with ages $t_1<t_2$,
$P(\text{all three same})=e^{-r(t_1+2t_2)}$ and
$P(AB)\,P(AC)/P(ABC)=e^{-rt_1}$. Three-way identity coefficients carry
topology information that pairwise ones do not. This is the
population-genetic route to a "triplet composite likelihood" (§5a).

**What does not transfer.** Griffiths–Tavaré IS, Kingman coalescent
priors and ESF-based estimators assume a *random genealogy with exchangeable
lineages*. The family tree is not exchangeable (it is the parameter). The
transferable pieces are (i) the exact fixed-tree recursion, (ii) the moment
statistics, (iii) the ESF as a null, and (iv) the IS/SMC machinery if
topology is integrated out (§5i).

---

## 2. Rate heterogeneity among meanings: power laws, not exponentials

### 2.1 Gamma frailty gives a Lomax (Pareto II) retention law

If a randomly chosen meaning has rate $r\sim\text{Gamma}(k,\theta)$ (mean
$\bar r=k\theta$, CV $1/\sqrt k$), the expected retention of a lineage is the
Laplace transform:

$$
S(t)=\mathbb E[e^{-rt}]=(1+\theta t)^{-k},\qquad
h(t)=-\frac{d\log S}{dt}=\frac{k\theta}{1+\theta t}=\frac{\bar r}{1+\theta t}.
$$

So **the population hazard declines with time even though every meaning has
a constant hazard**. This is the "frailty" effect of survival analysis
(Vaupel, Manton & Stallard 1979, *Demography* 16:439–454,
doi:10.2307/2061224 **[m]**): fast meanings are used up first, and the
survivors are the slow ones. For a pair, $F(T)=S(2T)$. The asymptotics are
a power law $S\sim(\theta t)^{-k}$, so cognate percentages never "run out"
the way the exponential says they should. This apparent slowdown is pure
selection on rates, with no change in any process.

Glottochronology met this early. van der Merwe (1966) "New mathematics
for glottochronology", *Current Anthropology* 7:485–500 **[v, no DOI
found]**, gave different classes of the list different rates. Dyen, James &
Cole (1967) "Language divergence and estimated word retention rate",
*Language* 43:150–171 **[v, DOI not found]**, gave each meaning its own rate.

### 2.2 Starostin's law is also heavy-tailed, in closed form

Starostin's revised glottochronology uses $c=e^{-\lambda c t^2}$
($\lambda\approx0.05$ ky$^{-2}$). The search confirmed the form and the
constant from secondary sources (e.g. the Vasilyev handout at the MPI-EVA
Swadesh centenary); the primary source, Starostin (1989; English 1999) **[u]**,
was not checked here. Solving the implicit equation with $cy\,e^{cy}=y$ for
$y=\lambda t^2$ gives

$$
c(t)=\frac{W(\lambda t^2)}{\lambda t^2}\qquad (W=\text{Lambert }W).
$$

Hence $c\approx 1-\lambda t^2$ at shallow depth, a *quadratic* start with zero
initial hazard, and $c\sim \ln(\lambda t^2)/(\lambda t^2)$ at depth: a
$t^{-2}$ power law, up to a log. Numerically (`e7_watterson.py`) $c=0.953,
0.845, 0.613, 0.447, 0.339, 0.265, 0.162, 0.110$ at
$t=1,2,4,6,8,10,15,20$ ky. Read as a pair-level retention, the tail behaves
like a gamma mixture with $k\approx2$ on the path length. The difference
is at the shallow end: a gamma mixture has a *linear* start (nonzero initial
hazard). These two are separable only with dated shallow pairs.

### 2.3 Lognormal vs gamma mixing

For log-normal rates (log-SD $s$) there is no closed form. The tail is
$S(t)\approx\exp(-(\ln t)^2/2s^2)$: slower than any exponential, but faster
than any power law eventually. At equal mean and CV, the lognormal keeps
*less* retention at depth than the gamma. From `e3_hetero.py`, with
$\bar r=0.15$/ky/lineage and CV 0.8, pair shared fractions at
$T=4,8,15$ ky are:

| | exp | gamma CV .5 | gamma CV .8 | gamma CV 1 | lognormal CV .8 |
|---|---|---|---|---|---|
| $T=4$ | 0.301 | 0.350 | 0.410 | 0.455 | 0.393 |
| $T=8$ | 0.091 | 0.153 | 0.234 | 0.294 | 0.200 |
| $T=15$ | 0.011 | 0.049 | 0.120 | 0.182 | 0.080 |

**The choice of mixing family matters only in the slow tail.** That is
exactly the region that dates deep nodes.

### 2.4 Bias from misspecification (numbers)

*(a) Heterogeneity ignored entirely* (a single mean rate, Swadesh-style). The
naive date is $-\ln F/(2\bar r)$. At true 4 / 8 / 15 ky with gamma CV 0.8 it
gives 2.97 / 4.85 / 7.06 ky: **too young, by more the deeper the node**.

*(b) Calibrated at 2 ky, then extrapolated with an exponential.* The calibration
absorbs the slowdown at 2 ky ($\hat r=0.127$ vs true mean 0.150). The
date/true ratio is then 1.08 at 1 ky, 1.00 at 2 ky, 0.88 at 4 ky, 0.79 at
6 ky and 0.66 at 10 ky. Shallow calibrations do not fix deep bias.

*(c) Relative chronology is biased too.* With a true ratio $t_1/t_2=0.5$, the
naive estimate is 0.54 (at $t_2=2$), 0.57 (4) and 0.61 (8): **compressed
toward the root**.

*(d) The case clepsydra is actually in.* Concept rates are modelled, but
there is **unmodelled concept × lineage heterogeneity**: the same meaning
replaced at different rates on different branches. Per-concept rates are
known exactly (200 meanings, lognormal across meanings), and each lineage's
rate for each meaning is multiplied by an independent Gamma($k$,$1/k$). The
pair MLE (`e4_bias.py`) is:

| $k$ (CV) | T=1 | 2 | 4 | 6 | 8 ky |
|---|---|---|---|---|---|
| 1 (1.0) | 0.90 | 1.69 | 3.07 | 4.30 | 5.43 |
| 2 (0.71) | 0.95 | 1.83 | 3.47 | 5.01 | 6.48 |
| 4 (0.5) | 0.97 | 1.91 | 3.71 | 5.46 | 7.17 |
| 10 (0.32) | 0.99 | 1.96 | 3.88 | 5.78 | 7.65 |

So **even with perfect concept rates**, moderate lineage-by-concept
variation (CV≈0.5) costs about 9% at 6–8 ky. Clepsydra's $u_{f,c}$ is a
*family-level* concept deviation shared by all branches, so it does not
absorb this. The branch multiplier $m_b$ is shared by all meanings, so it
does not either. A **"+Γ across branch × concept"** term is the lexical
analogue of Yang's (1994) +Γ among-site heterogeneity (*J. Mol. Evol.*
39:306–314, doi:10.1007/BF00160154 **[m]**). It is cheap: replace
$e^{-r t}$ on each branch by $(1+r t/k)^{-k}$. That is exact for a
per-branch, per-meaning gamma multiplier drawn independently on each branch,
and it keeps the peeling recursion unchanged because it is still a
per-branch survival probability. It adds one parameter $k$. With discrete
categories (Yang's approach) you would instead average whole-tree
likelihoods. That is the right model when the multiplier is per meaning
(shared across branches) rather than per branch. Both variants are cheap and
both are worth having. They answer different questions: "meaning-level
residual heterogeneity" vs "branch × meaning heterogeneity".

### 2.5 Can the mixing distribution be estimated, and across families?

- **Within one family.** Meanings appear only through current-status data, so
  the mixing distribution is weakly identified. Its effect only appears
  through the *curvature* of $\log F$ against depth across many node depths.
  In a 7-language family there are 6 node depths, and that is not enough.
- **Across families (the right place).** Pool calibrated families (IE-CoR and
  others) in a hierarchical model with
  $\log\lambda_c\sim\mathcal N(\mu_c,\sigma_c^2)$ for concept effects and a
  candidate family for the residual (gamma $k$, lognormal $s$, or a
  two-component mixture). Then compare families via leave-one-family-out
  predictive density (§4.3). The deep-time behaviour is then *learnt from
  deep calibrated families* rather than assumed.
- **Practical rule.** Report deep ages under at least two mixing families
  (gamma/Lomax and lognormal, matched on mean and variance from the
  reference data). Present the spread as a model-uncertainty component (§11).

---

## 3. Information-theoretic limits: what is achievable

### 3.1 One pair, one meaning, known rate

$y\sim\text{Bernoulli}(p)$ with $p=e^{-x}$ and $x=2rT$. The Fisher information about
$\log T$ is

$$
\mathcal I_{\log T}=\frac{(\partial p/\partial\log T)^2}{p(1-p)}=\frac{x^2 e^{-x}}{1-e^{-x}}=\frac{x^2}{e^{x}-1}.
$$

It is **maximised at $x^\star=1.594$**, i.e. when the two languages share
**20.3%** of the meaning's classes in expectation, with value **0.648
nats**. The Cramér–Rao bound on the CV of $T$ from $C$ meanings at the same
rate is therefore never below $1/\sqrt{0.648\,C}$: 0.124 for $C=100$ and
0.088 for $C=200$. From `e1_pairinfo.py`:

| $x=2rT$ | shared $p$ | info/meaning | CV, C=100 | CV, C=200 |
|---|---|---|---|---|
| 0.2 | 0.82 | 0.181 | 0.235 | 0.166 |
| 1.0 | 0.37 | 0.582 | 0.131 | 0.093 |
| 1.59 | 0.20 | 0.648 | 0.124 | 0.088 |
| 3.0 | 0.05 | 0.472 | 0.146 | 0.103 |
| 6.0 | 0.0025 | 0.090 | 0.334 | 0.236 |
| 8.0 | 0.0003 | 0.022 | 0.682 | 0.483 |

Interpretation:

- **Shallow saturation (too little change) is as real as deep saturation.**
  At $x=0.2$ (two languages ~600 years apart at average rates) each
  meaning carries only 28% of the optimum.
- **Each meaning is informative only in a depth window around
  $T\approx0.8/r_c$.** Deep nodes are dated by the slow meanings, shallow nodes
  by the fast ones. *So the slow-rate tail of the concept-rate prior is what
  determines deep dates.* Its misspecification (§2) is the main source of
  deep bias.
- **Heterogeneous rates widen the window.** With 200 meanings whose rates are
  lognormal around 0.15/ky (log-SD 0.8, known), the pair CV stays at
  0.10–0.15 for every $T$ from 1 to 10 ky (0.148, 0.118, 0.103, 0.101, 0.102,
  0.106 at 1, 2, 4, 6, 8, 10 ky). Rate diversity is a resource: a list made
  only of ultra-stable or only of fast meanings would collapse at one end.

### 3.2 The whole tree: exact Fisher information for a 7-language family

`e2_tree_fisher.py` enumerates all **877 set partitions** of 7 tips and
computes the exact T1 likelihood per partition. Patterns sum to 1 exactly. It
then computes the expected Fisher information for: 6 log internal ages; 12
log branch multipliers with prior $\mathcal N(0,\sigma_b^2)$; and the log
family rate with prior $\mathcal N(0,\sigma_f^2)$. Priors enter as precision
(the Bayesian information/van Trees bound, a lower bound on posterior SD for
large $C$). The tree is
`(((((A,B):0.2,C):0.4,(D,E):0.3):0.7,(F,G):0.5):1.0)` in relative ages, the
root has age $T$, and there are 200 meanings with lognormal rates (mean 0.15/ky,
log-SD 0.8).

**SD of log root age (≈ CV):**

| Root $T$ | data only (rates known, strict clock) | + relaxed clock σ_b=0.3 | + family-rate prior σ_f=0.3 | σ_b=0.3 & σ_f=0.3 | σ_b=0.3 & σ_f=0.5 |
|---|---|---|---|---|---|
| 1 ky | 0.111 | 0.157 | 0.320 | 0.339 | 0.524 |
| 2 ky | 0.084 | 0.139 | 0.312 | 0.331 | 0.519 |
| 4 ky | 0.068 | 0.130 | 0.308 | 0.327 | 0.517 |
| 6 ky | 0.064 | 0.128 | 0.307 | 0.326 | 0.516 |
| 8 ky | 0.062 | 0.127 | 0.306 | 0.326 | 0.516 |

With $C=100$, the data-only column is 0.157 / 0.119 / 0.096 / 0.090 / 0.088.
The σ_b=0.3, σ_f=0.3 column becomes 0.356 / 0.341 / 0.334 / 0.332 / 0.332.
Making each concept's rate uncertain (independent log-SD 0.5 per concept)
barely matters: the root SD goes from 0.068 to 0.080 at 4 ky, because
independent concept errors average out over 200 meanings.

**SD of log(node age / root age), i.e. relative chronology:**

| Root 4 ky, C=200 | AB (shallowest) | ABC | DE | ABCDE | FG |
|---|---|---|---|---|---|
| strict clock | 0.162 | 0.117 | 0.142 | 0.086 | 0.114 |
| relaxed σ_b=0.3 | 0.273 | 0.206 | 0.258 | 0.156 | 0.219 |

These ratios are **independent of σ_f**: relative chronology does not need
calibrations.

### 3.3 Reading the numbers: the irreducible CV for Arawan

1. **Data alone is not the bottleneck.** With 200 meanings, a 7-tip tree and
   known concept rates, the root age is pinned to ±6–8% (1 SD). The 7-tip tree
   is worth about 1.4–1.6× the information of a single pair at the same depth
   (0.068 vs 0.103 at 4 ky).
2. **Branch-rate variation roughly doubles it** (σ_b=0.3 → ±13%). This is the
   finite-lineage floor described by Zhu, dos Reis & Yang (2015) "Characterization of the
   uncertainty of divergence time estimation under relaxed molecular clock
   models using multiple loci", *Syst. Biol.* 64:267–280,
   doi:10.1093/sysbio/syu109 **[v]**. More meanings do *not* remove it,
   because all meanings share the branch multipliers. Only more independent
   lineages do (more languages, or more calibrated nodes).
3. **The family-rate prior dominates everything.** With σ_f=0.3 the root
   log-SD is about 0.33, so a 95% interval is a factor of
   $e^{\pm1.96\times0.33}\approx\times/\div1.9$. A median of 2000 BP becomes
   roughly 1050–3800 BP. With σ_f=0.5 (clepsydra's current HalfNormal(0.5)
   hyperprior scale) it is ×/÷2.8. **That is the honest answer for Arawan
   with zero calibrations, and no amount of modelling sophistication changes
   it.** Only a better-justified σ_f does.
4. Empirical anchors for σ_f: Holman et al. (2011) "Automated dating of the
   world's language families based on lexical similarity", *Current
   Anthropology* 52:841–875, doi:10.1086/662127 **[v]**, report a mean
   absolute discrepancy of **29%** of the estimated date across 52 calibration
   groups. That is consistent with a log-scale SD of about 0.3–0.4.
   Rama & Wichmann (2020) "A test of Generalized Bayesian dating", *PLoS ONE*
   15:e0236522, doi:10.1371/journal.pone.0236522 **[v]**, regress known ages
   on tree length over 30 independent calibration groups. They report
   leave-one-out MAE 770 y and RMSE 936 y, and some known ages fall outside
   their 95% HPDs. Both are direct precedents for "transfer a clock across
   families", and both say the cross-family floor is around 30% (1 SD).
5. **This matches Rannala & Yang's infinite-sites theorem.** With infinite
   data, posterior age uncertainty converges to the uncertainty of the
   calibrations or rate prior, scaled by the tree: Yang & Rannala (2006)
   *MBE* 23:212–226, doi:10.1093/molbev/msj024 **[m]**; Rannala & Yang (2007)
   "Inferring speciation times under an episodic molecular clock", *Syst.
   Biol.* 56:453–466, doi:10.1080/10635150701420643 **[m]**; dos Reis & Yang
   (2013) "The unbearable uncertainty of Bayesian divergence time
   estimation", *J. Syst. Evol.* 51:30–43, doi:10.1111/j.1759-6831.2012.00236.x
   **[m]**. The table shows Arawan is already near that limit at
   $C\approx100$–200.

**Corollary for presentation.** Report three nested numbers per node:
(i) *lexical age* $\rho a_v$ with its data-only interval, (ii) *relative age*
$a_v/a_\text{root}$, and (iii) *calendar age*. The first two can be sharp
even when the third is wide. Users then see where the width comes from.

---

## 4. Identifiability and transferring an absolute scale

### 4.1 What is identified with zero calibrations

The likelihood $\ell(\rho, m, a, \lambda)$ depends on $(\rho,a)$ only through
$\{\rho m_b t_b\}$. Under a **strict clock** ($m_b\equiv1$):

- identified: $\rho a_v$ for all $v$ (lexical ages), and so all ratios
  $a_v/a_w$;
- not identified: $\rho$ and $a_\text{root}$ separately (a one-dimensional
  ridge $\rho a=\text{const}$).

Under a **relaxed clock**, the $2L-2$ products $\rho m_b t_b$ (branch lexical
lengths) are identified. The $L-1$ ages and $2L-2$ multipliers are not: they
are separated only by the prior on $m_b$. So relative chronology becomes
prior-dependent too. The §3.2 numbers show by how much (ratio SD 0.16 → 0.27
for the shallowest node).

**A useful exact factorisation.** Suppose the prior on ages is
**scale-invariant**: e.g. $p(a_\text{root})\propto 1/a_\text{root}$ with an
arbitrary prior on the ratios $a_v/a_\text{root}$, and $\rho$ is independent.
Then in log coordinates

$$
\log a_v = \log(\rho a_v) - \log\rho,
$$

and the posterior of the lexical ages does not depend on the rate prior. The
calendar posterior is the **convolution** of the lexical-age posterior with
the prior on $-\log\rho$. In practice: run the MCMC once in lexical time,
then produce calendar-time posteriors for *any* rate prior by adding
independent draws of $-\log\rho$. Rate-prior sensitivity analysis (§6) then
costs nothing. It is exact for T1/T2 without calibrations. With a Yule/BD
age prior whose rate has its own prior, the factorisation is approximate.

### 4.2 Transferring the scale: hierarchical cross-family rate model

Let $\beta_f=\log\rho_f$. For calibrated families $f=1..F$, the posterior
$p(\beta_f\mid D_f)$ is roughly $\mathcal N(b_f,s_f^2)$. The random-effects
meta-analysis is

$$
b_f\sim\mathcal N(\mu,\ \sigma_f^2+s_f^2),
\qquad \beta_\text{new}\mid\text{data}\sim\mathcal N\!\left(\hat\mu,\ \sigma_f^2+\operatorname{Var}(\hat\mu)\right),
$$

with $\operatorname{Var}(\hat\mu)\approx(\sigma_f^2+\bar s^2)/F$. The
between-family SD $\sigma_f$ is the honest floor. More calibrated families
shrink $\operatorname{Var}(\hat\mu)$ but never $\sigma_f$.

- With $F=10$ families and $\sigma_f=0.3$: $\sqrt{0.09+0.009}=0.315$. With
  $F=30$: 0.305. **Beyond about 10 well-calibrated families, adding
  families no longer narrows Arawan's interval.** What helps is *explaining*
  $\sigma_f$ with covariates.
- Covariates for $\beta_f$: speaker population, writing tradition,
  contact intensity, geography and subsistence. Bromham et al. (2015) "Rate
  of language evolution is affected by population size", *PNAS*
  112:2097–2102, doi:10.1073/pnas.1419704112 **[m]**, found that word loss
  is faster in smaller populations. That is directly relevant to small
  Amazonian families: an Arawan-specific rate prior may need to be shifted,
  not just widened. A meta-regression
  $\beta_f=\mu+\gamma^\top z_f+\eta_f$ reduces $\sigma_f$ only if
  $\gamma$ is real. This must be tested by leave-one-family-out (LOFO).
- The rate scale must be estimated in the **same model** used to date the
  target. A $\rho$ learnt under T2 with ascertainment is not a $\rho$ for T1.
  Concept-rate definitions and synonym handling must match too.
  Otherwise the transfer carries a systematic offset that no interval
  reports.

### 4.3 Leave-one-family-out validation (the decisive check)

For each calibrated family $g$: refit the cross-family hyperparameters
without $g$, date $g$'s calibrated nodes *as if uncalibrated*, and record
the PIT value $u_g=P(a_\text{node}\le a_\text{true}\mid D_{-g})$. Under a
calibrated procedure the $u$'s are Uniform(0,1). Use the coverage of the
50/80/95% intervals, the log score of the true ages, and the bias of
$\log(\hat a/a_\text{true})$ as a function of depth (§2.4 predicts a
negative slope if the mixing tail is wrong). This is the only test that
measures what users of Arawan dates actually rely on. With ~10–30 families
it is coarse but informative (§6).

---

## 5. Alternative inference schemes

Throughout, "C99 fit" means suitability for a small, dependency-free C99
engine.

### 5a. Composite (pairwise) likelihood with an ultrametric constraint

$\ell_\text{CL}(a)=\sum_{i<j}\sum_c \log P(y_{ijc}\mid a_{\text{mrca}(i,j)})$.
This is Lindsay (1988) **[m]**; for the review see Varin, Reid & Firth (2011)
*Statistica Sinica* 21:5–42 **[m]**. It is unbiased (each term is a true
marginal), but its curvature $H$ overstates information because pairs
sharing a branch are dependent. The correct asymptotic variance is the
**Godambe/sandwich** $G^{-1}=H^{-1}JH^{-1}$, with $J=\operatorname{Var}(\text{score})$.

Exact computation on the 7-tip tree (`e5_composite.py`), strict clock, known
rates, $C=200$:

| Root $T$ | exact full-likelihood SD(log root) | pairwise CL, Godambe | pairwise CL, naive $H^{-1}$ |
|---|---|---|---|
| 1 ky | 0.111 | 0.116 | 0.047 |
| 2 ky | 0.084 | 0.090 | 0.037 |
| 4 ky | 0.068 | 0.076 | 0.033 |
| 8 ky | 0.062 | 0.072 | 0.032 |

Averaged over all 6 nodes, the CL/exact SD ratio is 1.06–1.07.

**Pairwise CL is about 88% efficient (in variance) relative to the full T1
likelihood**, so very little information is lost. **But a naive pairwise
posterior** — one that multiplies the 21 pair likelihoods as if independent —
**is overconfident by a factor of about 2.2–2.4 in SD.** The overconfidence
grows roughly like $\sqrt{L}$ (more pairs share each branch). If clepsydra
ever combines pairwise likelihoods over a tree, it must curvature-adjust
(Ribatet, Cooley & Davison 2012 for the Bayesian composite-likelihood
adjustment **[u]**) or temper by $\approx H/J$. Cheap and robust; C99 fit:
excellent.

A natural middle ground is a **triplet composite likelihood**, which recovers
topology information (§1.4) at $O(L^3)$ terms.

### 5b. Distance methods

The infinite-state "Jukes–Cantor" distance per pair is
$\hat d_{ij}=\arg\max_d\sum_c[y_c\log e^{-\lambda_c d}+(1-y_c)\log(1-e^{-\lambda_c d})]$
(rates known up to $\rho$), with variance from §1.4. Then:

- **UPGMA** on $\hat d$ gives an ultrametric tree and lexical node ages. It
  is consistent under a strict clock.
- **NJ/BIONJ** (Saitou & Nei 1987, doi:10.1093/oxfordjournals.molbev.a040454
  **[m]**; Gascuel 1997, doi:10.1093/oxfordjournals.molbev.a025808 **[m]**)
  followed by **least-squares dating** (LSD/LSD2: To, Jung, Lycett & Gascuel
  2016, *Syst. Biol.* 65:82–97, doi:10.1093/sysbio/syv068 **[m]**). LSD
  minimises weighted squared deviations between branch lengths and the
  clock-implied lengths, under calibration constraints, in near-linear time.
  The variance weights from §1.4 slot in directly.

*What they buy:* an instant point estimate, a starting tree for MCMC, a
topology screen for §5j, and an independent check. *Costs:* no principled
uncertainty (a parametric bootstrap over meanings is the fix, with meanings
as the resampling unit), and dependence between distances is ignored.
C99 fit: excellent (UPGMA/BIONJ are ~200 lines).

### 5c. ABC and neural simulation-based inference (SBI)

ABC: Beaumont, Zhang & Balding (2002) *Genetics* 162:2025–2035,
doi:10.1093/genetics/162.4.2025 **[m]**. Neural SBI (NPE/NLE/NRE):
Cranmer, Brehmer & Louppe (2020) *PNAS* 117:30055,
doi:10.1073/pnas.1912789117 **[m]**; `sbi` toolkit (Tejero-Cantero et al.
2020, *JOSS*, doi:10.21105/joss.02505 **[m]**).

*Honest assessment.*

- For the core models the likelihood is **exact and cheap** (§1.3). SBI
  then adds approximation error and buys nothing.
- SBI's real selling point is models with **no tractable likelihood**:
  borrowing with realistic contact geometry, synonym dynamics
  (gain-then-loss of secondary forms), word-level coding error, and
  semantic-shift chains. For a 7-tip family, summary statistics are natural
  ($F_{ij}$ matrix, $K_c$ histogram, triplet identities), and an amortised
  NPE could be trained once over a wide prior and reused across families.
- **Risks are serious.** Hermans et al. (2022) "A trust crisis in
  simulation-based inference? Your posterior approximations can be
  unfaithful", *TMLR*, arXiv:2110.06581 **[v]**, show that standard SBI
  methods can give overconfident posteriors and that ensembling helps.
  Cannon, Ward & Schmon (2022), arXiv:2209.01845 **[v]**, show that model
  misspecification badly degrades neural SBI. A real family is never exactly
  simulable, so this is the core problem, not a corner case.
- Best use: **as a validation tool rather than the engine.** Simulate from a
  richer "world model" (borrowing, coding error, heterogeneity) and run
  clepsydra on the simulations. Measure bias and coverage (§6), and train an
  NPE on the same simulations only as a *second opinion*. C99 fit: none
  (Python side only), and that is fine.

### 5d. Variational inference

Mean-field or full-rank Gaussian VI on (log ages via ratio transform, log
rates). Fast, but VI is known to **underestimate posterior variance**, and
here the main deliverable is the variance. Not recommended as the primary
engine. Acceptable for initialisation or for screening topologies. C99 fit:
moderate (needs gradients and a stochastic optimiser).

### 5e. Laplace / INLA

The posterior in (log root age, logit ratios, log rates, log multipliers)
is smooth and, for $C\ge100$, close to Gaussian (§3: SDs of 0.07–0.3 on the
log scale). INLA (Rue, Martino & Chopin 2009, *JRSS-B* 71:319–392,
doi:10.1111/j.1467-9868.2008.00700.x **[m]**) uses nested Laplace: Laplace
over the many latent parameters ($m_b$, concept deviations) plus a grid or
adaptive quadrature over a few hyperparameters ($\sigma_b$, $\sigma_f$,
$\varepsilon$). Each Laplace step needs a mode search and a Hessian. With the
outside pass of §1.3, gradients cost one extra traversal. Hessians come
from finite-differenced gradients ($O(P)$ gradient calls for
$P\approx20$–40). **Approximate marginal likelihood per topology** follows as
a by-product, and that is what makes §5j feasible. C99 fit: good (BFGS plus a
dense Cholesky of a 40×40 matrix).

### 5f. Deterministic quadrature for small trees

For 7 tips there are 6 internal ages. Use **adaptive Gauss–Hermite
quadrature (AGHQ)** in the Laplace-whitened coordinates. With 5 nodes per
dimension that is $5^6=15{,}625$ likelihood evaluations. At roughly
2×10⁴ flops per evaluation (200 meanings × ~100 peeling ops) that is about
0.2 s in C per topology, giving **essentially exact** marginal posteriors
of every age and the evidence. Smolyak sparse grids cut this to a few
hundred points at level 3–4. Branch multipliers (12 dims) cannot be
tensor-gridded; integrate them by Laplace inside each grid point
(INLA-style). C99 fit: good. This is the most attractive "exact"
alternative to MCMC for families with $L\le8$.

### 5g. Nested sampling

Skilling (2006) *Bayesian Analysis* 1:833–860, doi:10.1214/06-BA127
**[m]**. It gives the evidence $Z$ and posterior samples together, and copes
with multimodality (e.g. ties between near-equal clades). Cost is ~10⁴–10⁵
likelihood calls for 10–20 dimensions, i.e. 0.1–1 s per topology in C. It
suits model comparison (T1 vs T2 vs +Γ vs contamination) and is a
**cross-check of the Laplace evidences** on the top topologies. Too slow to
run on all 10,395 topologies naively (about 1–3 h), which is still
feasible overnight. C99 fit: good (simple constrained-prior sampler via
slice steps).

### 5h. HMC with the ratio transform

Ji et al. (2023) "Scalable Bayesian divergence time estimation with ratio
transformations", *Syst. Biol.* 72:1136–1153, doi:10.1093/sysbio/syad039
**[v]**, map the $N-1$ node heights to one height plus $N-2$ ratios in
$(0,1)$ and run HMC. The ratio parametrisation is older (Kishino, Thorne &
Bruno 2001 **[u]**). This is the right parametrisation even without HMC:
it decouples the scale ($\rho a_\text{root}$, strongly determined) from the
shape (ratios), which is exactly the identifiability structure of §4.1.
Clepsydra's "block age-scaling" move is a step toward this. A
**ratio + log-scale reparametrisation of the MWG sampler** would remove
most of the remaining posterior correlation. Full HMC needs gradients (§1.3)
and a leapfrog integrator: C99 fit moderate, worthwhile if larger trees
matter.

### 5i. SMC over topologies

Bouchard-Côté, Sankararaman & Jordan (2012) "Phylogenetic inference via
sequential Monte Carlo", *Syst. Biol.* 61:579–593,
doi:10.1093/sysbio/syr131 **[m]**. It builds trees by successive mergers
(a natural fit with ultrametric, clock-like trees: merge the closest pair,
proposing the merge time), gives an unbiased evidence estimate, and
parallelises well. For $L\le8$, enumeration (§5j) is simpler and exact. For
$L=10$–30 it is the natural successor. C99 fit: moderate.

### 5j. Exhaustive topology enumeration with Bayesian model averaging

Counts (`e6_topo.py`):

| L | rooted binary $(2L-3)!!$ | ranked histories $L!(L-1)!/2^{L-1}$ | rooted, any degree (OEIS A000311) |
|---|---|---|---|
| 5 | 105 | 180 | 236 |
| 6 | 945 | 2,700 | 2,752 |
| 7 | **10,395** | 56,700 | 39,208 |
| 8 | **135,135** | 1,587,600 | 660,032 |
| 9 | 2,027,025 | 57,153,600 | 12,818,912 |
| 10 | 34,459,425 | 2.57×10⁹ | 282,137,824 |

**For Arawan ($L\approx7$) full Bayesian model averaging over topologies by
brute force is cheap.** Laplace/INLA evidence per topology at ~10–50 ms in C
gives 2–9 minutes for L=7 and 20 minutes to 2 hours for L=8. AGHQ refinement
(§5f) then runs on the topologies holding 99.9% of the Laplace mass. Then:

$$
p(a_{\text{MRCA}(S)}\mid D)=\sum_\tau p(\tau\mid D)\,p(a_{\text{MRCA}_\tau(S)}\mid\tau,D),
\qquad p(\tau\mid D)\propto p(\tau)\,Z_\tau .
$$

Subtleties:

- **The topology prior matters.** A Yule process is uniform on *ranked*
  histories, so it is not uniform on topologies: balanced trees get more
  weight in proportion to their number of rankings. Pick a prior on
  purpose and state it. The evidences $Z_\tau$ must include the
  topology-conditional age prior density, normalised (e.g. the conditioned
  Yule/BD density: Gernhard 2008, *J. Theor. Biol.* 253:769–778,
  doi:10.1016/j.jtbi.2008.04.005 **[m]**).
- **MRCA age of a named set is defined in every topology** (it may be the
  root). This makes BMA for user-named clades well defined, unlike "age of
  clade X", which needs X to be monophyletic. Report $P(S\text{ monophyletic})$
  alongside.
- **Laplace evidences are biased in a topology-dependent way** when a
  topology's posterior presses against a boundary, e.g. a near-zero internal
  branch. Use the ratio/log parametrisation (no hard boundaries) and refine the
  top set with AGHQ or nested sampling.
- Multifurcations: either allow them as topologies (A000311 counts), or
  rely on near-zero branches. With ~200 meanings, near-polytomies are common
  in small families. A soft polytomy prior is the honest treatment.
- C99 fit: **excellent**. A topology iterator is ~100 lines, and the rest
  reuses the fixed-tree likelihood. A fixed-tree engine becomes a
  topology-integrating engine for small families with no tree-MCMC moves
  at all.

---

## 6. Calibration of uncertainty

### 6.1 Simulation-based calibration (computational correctness)

SBC (Talts et al. 2018, arXiv:1804.06788 **[m]**; Modrák et al. 2023
*Bayesian Analysis*, doi:10.1214/23-BA1404 **[m]**): draw parameters from the
prior, simulate data, fit, and check that posterior ranks are uniform. This
checks the *sampler and likelihood code*, not the model. Clepsydra's
simulation validation already does something close. Include $\log a_\text{root}$,
ratios and $K$-based statistics as test quantities (Modrák et al. stress that
the choice of quantity shapes sensitivity).

### 6.2 Coverage under misspecification (the real question)

Draw data from a *richer* simulator than the fitted model: gamma
branch×concept heterogeneity, borrowing, coding error and split pulses
(§§7–10). Measure frequentist coverage of the 50/80/95% intervals as a
function of depth. The §2.4/§8/§9 bias tables predict where coverage will
fail: deep nodes, high borrowing, lumping.

### 6.3 Recalibrating intervals with calibrated families

- **Conformal-style recalibration.** On $F$ held-out calibrated families,
  compute the nonconformity score $s_g=|\log a_\text{true}-\log\hat a_g|/\hat
  \sigma_g$ (posterior median and SD on the log scale). The
  $\lceil(1-\alpha)(F+1)\rceil$-th smallest score gives an inflation factor for
  the target family. Validity needs exchangeability of families, which is
  doubtful when calibrated families are mostly Eurasian, written-tradition
  families and the target is Amazonian. Granularity is $1/(F+1)$: 90% needs
  $F\ge9$. See Angelopoulos & Bates, "A gentle introduction to conformal
  prediction", arXiv:2107.07511 **[m]**.
- **Bayesian alternative (preferred).** Put the inflation inside the model
  as a "dating error" random effect with its own hyperprior, learnt through
  the LOFO likelihood. This is the same as the σ_f meta-analysis of §4.2, so
  in practice **a well-built cross-family hierarchy *is* the recalibration**.
  Conformal is then a check that the hierarchy is not overconfident.
- Rama & Wichmann (2020) is a warning sign. Known ages for some groups fall
  outside the 95% HPDs and the paper does not report a coverage rate. A
  LOFO coverage table should be a standard clepsydra output.

### 6.4 Prior sensitivity / robust Bayes

- **Power-scaling sensitivity** (Kallioinen, Paananen, Bürkner & Vehtari
  2023/24, "Detecting and diagnosing prior and likelihood sensitivity with
  power-scaling", *Stat. Comput.* 34:57, doi:10.1007/s11222-023-10366-5
  **[v]**). Raise a prior component to the power $\alpha$ and reweight
  existing posterior draws by $p(\theta)^{\alpha-1}$ with Pareto-smoothed
  IS. The derivative of a distance between posteriors at $\alpha=1$ flags
  prior-data conflict versus prior dominance. Do this per component: rate
  prior, $\sigma_b$, tree prior, contamination. For Arawan, expect the rate
  prior to show "prior-dominated, no conflict". That is the correct and
  reportable diagnosis. *Implementation:* post-processing on stored draws, in
  Python or C, with no refits. Cheap.
- **ε-contamination classes** (Berger 1994, *Test* 3:5–124,
  doi:10.1007/BF02562676 **[m]**): report the range of a posterior
  quantity over $\{(1-\epsilon)\pi_0+\epsilon q\}$. For a single positive
  scale parameter this becomes a closed-form bound via §4.1 (the
  convolution), so it is computable exactly for the rate prior.
- **Tree-prior sensitivity.** Yule vs BD vs scale-invariant "uniform
  ratios". With few tips the age prior is not swamped by data on relative
  ages under a relaxed clock (§3.2), so this must be reported.

---

## 7. Survival-analysis framing

**Each (meaning, lineage) is a survival unit, and we observe current status
only.** Pair cognacy is $\mathbb 1\{T_i>t_i,\ T_j>t_j\}$ (both kept the
ancestral class), with the twist that a class *shared* via a deeper
ancestor is still a "survival". That is case-I interval censoring. Its
nonparametric MLE converges at $n^{1/3}$, not $n^{1/2}$ (Groeneboom &
Wellner 1992, *Information Bounds and Nonparametric Maximum Likelihood
Estimation*, Birkhäuser **[m]**). So **the shape of the retention curve
cannot be learnt nonparametrically from cognate data at any useful
precision.** It must come from a parametric family or from many calibrated
families.

Useful constructs:

- **Frailty = rate mixing.** Gamma frailty gives Lomax (§2.1). Positive
  stable frailty gives $S=\exp(-(\theta t)^{\alpha})$, a Weibull with
  $\alpha<1$: a *stretched exponential*, another closed form for the
  apparent slowdown, with one parameter. Inverse-Gaussian frailty gives
  $S=\exp\{\frac{1}{\kappa}(1-\sqrt{1+2\kappa\bar r t})\}$ **[derivation:
  Laplace transform of IG; standard]**. All three are drop-in replacements
  for $e^{-rt}$ per branch (§2.4).
- **Proportional hazards with meaning covariates.**
  $\log\lambda_c=\mu+\gamma_1\log\text{freq}_c+\gamma_2\,\text{POS}_c+\dots+u_c$.
  Pagel, Atkinson & Meade (2007) "Frequency of word-use predicts rates of
  lexical evolution throughout Indo-European history", *Nature*
  449:717–720, doi:10.1038/nature06176 **[m]**, found that frequency explains
  much of the among-meaning rate variance, and the effect replicated in other
  families. Calude & Pagel (2011) *Phil. Trans. R. Soc. B* 366:1101–1107,
  doi:10.1098/rstb.2010.0315 **[m]**, found cross-linguistically shared
  frequency patterns. *Why it matters here:* a covariate model transfers to a
  new family better than a per-concept random effect estimated in IE,
  because it only needs the covariate. It also shrinks unstable estimates
  for meanings that are rare in the reference data.
- **Lineage covariates** (population size, contact, writing) enter as
  branch-level PH multipliers: $m_b=\exp(\delta^\top w_b+\eta_b)$. This is a
  regression on the relaxed clock, identifiable only across many families.
- **Time-varying hazards** (an initial "lag", as in Starostin's quadratic
  start, or acceleration) are hazards $h(t-\text{age of word})$. They are
  *not* Markov in the class state, so the peeling recursion needs the age of
  the current word. That expands the state space. Avoid unless calibrated
  evidence demands it.

---

## 8. Measurement error in cognacy

Let the true partition be $\pi_c$ and the coded partition $\tilde\pi_c$. The
simplest word-level error model has two parameters:

- **split** ($\varepsilon$): a word is wrongly placed in its own new class;
- **lump** ($\eta$): a non-cognate pair is wrongly joined (chance resemblance,
  undetected loan, over-eager etymology).

Pairwise consequences: $\tilde F=F(1-\varepsilon)$ with $\varepsilon$ the per-*pair*
split probability (≈ twice the per-word rate), and $\tilde F=F+(1-F)\eta$.

Biases from `e4_bias.py` (pair MLE, 200 meanings, lognormal known rates):

| Error | 0.5 ky | 1 | 2 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|
| split 2% | 0.57 | 1.08 | 2.08 | 4.10 | – | 8.13 | – |
| split 5% | 0.68 | 1.19 | 2.22 | 4.26 | – | 8.33 | – |
| split 10% | 0.88 | 1.40 | 2.45 | 4.53 | – | 8.67 | – |
| lump 1% | – | 0.99 | 1.97 | 3.90 | 5.79 | 7.62 | 9.38 |
| lump 3% | – | 0.96 | 1.90 | 3.71 | 5.41 | 6.96 | 8.38 |
| lump 5% | – | 0.94 | 1.84 | 3.54 | 5.07 | 6.41 | 7.58 |

- **Splitting is an additive offset** (about +80–180 y per 5% at all depths),
  so it is a *relative* disaster for shallow nodes (+36% at 500 y with 5%)
  and negligible at depth. Under T1 a word-level split is *exactly*
  indistinguishable from a replacement on the terminal branch. So it is
  absorbed by a **tip-offset parameter** (extra lexical length on each
  terminal branch), which is weakly identified, confounded with shallow
  ages, and needs a prior from double-coding studies.
- **Lumping is a floor.** The deep curve saturates at $\eta$, so deep nodes
  come out too young (−16% at 10 ky with 3%). Like borrowing (§9),
  it imitates a slowing clock.
- **Latent-partition models.** Put a prior on true partitions given the tree
  (T1 itself) and a likelihood $P(\tilde\pi\mid\pi)$, then sum over $\pi$.
  Exact summation is over all partitions consistent with the coded ones.
  For small $L$ this is feasible (Bell(7)=877 per meaning), so for Arawan
  **coding error can be marginalised exactly**. A "soft cognacy" input
  ($P(\text{pair cognate})$ from an expert or an automatic detector) fits the
  same frame. Rama, List, Wahle & Jäger (2018, NAACL) **[u: venue/DOI not
  checked]** discuss automated cognate detection for phylogenetics. List
  (2016) "Beyond cognacy", *J. Lang. Evol.* 1:119–136,
  doi:10.1093/jole/lzw006 **[m]**, discusses partial cognacy.
- **Coder-level effects.** Split-happy vs lump-happy coders shift *all*
  dates of a family in opposite directions. A cross-family rate prior learnt
  from coders with different habits picks up this bias in $\sigma_f$. This
  is one reason the cross-family floor is ~0.3 and not smaller.

---

## 9. Borrowing

**Two-lineage model.** Each replacement is, with probability $b$, a copy of
the other lineage's current class. The identity process satisfies
$\dot F=2rb(1-F) - 2r(1-b)F$ (borrowing makes a non-identical pair identical; a novel replacement breaks identity), which simplifies to
$\dot F=2rb-2rF$, so

$$
F(T)=b+(1-b)e^{-2rT}.
$$

Same floor form as lumping. Dates from `e4_bias.py`:

| b | 1 | 2 | 4 | 6 | 8 ky |
|---|---|---|---|---|---|
| 5% | 0.94 | 1.84 | 3.54 | 5.07 | 6.41 |
| 10% | 0.87 | 1.69 | 3.14 | 4.36 | 5.35 |
| 20% | 0.76 | 1.42 | 2.51 | 3.33 | 3.94 |

This agrees with simulation studies. Greenhill, Currie & Gray (2009) "Does
horizontal transmission invalidate cultural phylogenies?", *Proc. R. Soc.
B* 276:2299–2306, doi:10.1098/rspb.2008.1944 **[v]**: topologies are robust,
but dates "show a tendency to underestimate". Neureiter et al. (2022)
"Detecting contact in language trees: a Bayesian phylogenetic model with
horizontal transfer" (contacTrees, BEAST 2), *Humanit. Soc. Sci. Commun.*
9:205, doi:10.1057/s41599-022-01211-7 **[v]**: ignoring contact distorts
tree height, rate and topology. Kelly & Nicholls (2017) "Lateral transfer
in Stochastic Dollo models", *Ann. Appl. Stat.* 11:1146–1168,
doi:10.1214/17-AOAS1040 **[v]**: an exact-likelihood lateral-transfer model
for binary cognate traits.

Options for clepsydra, in increasing cost:

1. **Detection by posterior predictive.** Pairwise residuals
   $F^\text{obs}_{ij}-F^\text{pred}_{ij}$. Borrowing creates positive residuals
   for geographically adjacent non-sisters. Also look at the $K_c$
   distribution and triplet-identity residuals. Cheap, and it should be
   default output.
2. **Per-meaning mixture (already present).** Contamination soaks up
   heavily borrowed meanings. Too blunt for diffuse borrowing.
3. **Floor model.** Replace $e^{-r d}$ by $b_{ij}+(1-b_{ij})e^{-rd}$ in a
   pairwise model, with $b_{ij}$ a function of geographic distance. Only a
   pairwise CL (with Godambe correction) gives it a tractable likelihood.
4. **Tree-plus-transfer.** Kelly & Nicholls-type models: exact but complex.
   With $L=7$, a *network* with one or two reticulation edges is tractable as
   a mixture of displayed trees (each meaning follows one displayed tree).
   The T1 recursion runs on each tree and the results are mixed per meaning.
   This is the cleanest borrowing model compatible with the existing engine.
   Cost: ×2 per reticulation.

Dated borrowing *among sisters* makes the sisters look younger. Borrowing
between non-sisters pulls the topology. Both are visible in pairwise
residuals.

---

## 10. Non-constant processes

- **Punctuational (split-pulse) change.** Atkinson et al. (2008) "Languages
  evolve in punctuational bursts", *Science* 319:588,
  doi:10.1126/science.1149683 **[m]**, building on Pagel, Venditti & Meade
  (2006) *Science* 314:119–121, doi:10.1126/science.1129647 **[m]**. Model:
  each branch gets extra lexical length $\kappa$ at its start, so the branch
  length becomes $r(t_b+\kappa)$. Bias if ignored (`e4_bias.py`, κ in
  mean-rate-ky): with κ=0.2, true 0.5 → 0.70, 1 → 1.20, 4 → 4.20. With
  κ=0.5: 0.5 → 1.00, 1 → 1.50. **An additive offset, like split coding
  error.** *Identifiability:* κ is separated from the rate only by variation
  in the **number of nodes along root-to-tip paths** (the node-density
  test). In a 7-tip tree that is 2–4 nodes per path. It is barely
  identified, so it should be a prior-driven sensitivity term, not a free
  estimate.
- **Time-dependent rates** (rates that appear faster over short timescales,
  cf. Ho et al. 2005, *MBE* 22:1561–1568, doi:10.1093/molbev/msi145
  **[m]**): in lexical data these are mostly produced by mixing (§2) and
  coding error (§8). Model the mechanisms rather than a time-dependency
  curve.
- **Autocorrelated rates** (Thorne, Kishino & Painter 1998, *MBE*
  15:1647–1657, doi:10.1093/oxfordjournals.molbev.a025892 **[m]**;
  OU on log-rates): child log-multiplier $\sim\mathcal N(\phi\,\text{parent},
  \sigma^2(1-\phi^2))$. For $L=7$ there are 12 branches. The data cannot
  distinguish autocorrelated from uncorrelated rates (Drummond et al. 2006,
  *PLoS Biol.* 4:e88, doi:10.1371/journal.pbio.0040088 **[m]**, motivate
  model comparison, and it is weak at this size). It matters only through
  the prior it induces on relative ages. Include as a sensitivity case.
- **Rate shifts in a clade** (e.g. a contact-heavy subgroup): a
  local-clock model with a few shift points. Enumeration over shift
  locations is possible for $L\le8$ (12 branches → 12 single-shift models).

General principle: **every non-constant process enters the likelihood
through branch lexical lengths**. Without calibrations, it is separated from
time only by the prior. Its main effect on reported dates is through the
*relative* ages, and it should be shown as model uncertainty (§11).

---

## 11. Model averaging vs selection

- **The data unit is the meaning.** Leave-one-meaning-out predictive
  densities are natural and cheap via PSIS-LOO (Vehtari, Gelman & Gabry
  2017, *Stat. Comput.* 27:1413–1432, doi:10.1007/s11222-016-9696-4 **[m]**)
  from stored draws. Each draw's per-meaning likelihood is already computed
  in the sampler.
- **Stacking** (Yao, Vehtari, Simpson & Gelman 2018, *Bayesian Analysis*
  13:917–1007, doi:10.1214/17-BA1091 **[m]**) of T1, T1+Γ, T1+tip-offset and
  T1-with-contamination by maximising LOO log score. It is preferable to BMA
  when all models are wrong. **Caveat:** stacking optimises prediction of
  *cognate partitions*, not of *ages*. Two models with equal predictive
  fit can give different ages (that is the identifiability problem). So
  stacking weights should be complemented by LOFO age-prediction scores
  (§4.3), which target what you report.
- **Comparable data level.** T2 (binary characters with ascertainment)
  predicts a different data object than T1 (partitions). Map T2 to
  partitions per meaning, or compare only on pairwise identity, before any
  weighting.
- **Honest reporting.** Report the stacked/averaged age posterior, plus the
  per-model medians in a small table. Decompose variance by the law of total
  variance:
  $\operatorname{Var}(\log a)=\mathbb E_M[\operatorname{Var}(\log a\mid M)]+\operatorname{Var}_M(\mathbb E[\log a\mid M])$.
  The second term is "model uncertainty" and should be printed. For deep
  nodes, expect the mixing-family choice (§2.3) to dominate it.

---

## 12. Other directions found along the way

- **Optimal meaning lists.** Per-meaning information peaks where
  $2r_cT\approx1.6$ (§3.1). For a target depth one can *weight* or *select*
  meanings by expected Fisher information. More importantly, report which
  meanings drive each node date (influence via PSIS weights). For Arawan the
  root date will be driven by the ~30 slowest meanings, so their rate
  priors deserve scrutiny.
- **Ancestral-state uncertainty and synonyms.** Synonymy can be modelled as
  a birth–death of *secondary* forms within a slot (polymorphism),
  analogous to allele frequency dynamics. That is the step from
  infinite-alleles to a Moran-type model with drift inside each slot.
  Expensive. Not recommended for the core.
- **Rate-prior convolution (§4.1) as a product feature.** Ship lexical-age
  posteriors and let users (or the web page) convolve them with alternative
  rate priors instantly.
- **Bootstrap over meanings** as a model-free uncertainty check of the
  data-only component. It should reproduce §3.2's "data-only" SD. If it is
  much larger, the model is misspecified (meanings are more heterogeneous
  than the model allows).

---

## 13. Ranked recommendations (for small uncalibrated families like Arawan)

Ranked by value per unit of effort for the stated goal: well-calibrated
calendar ages for every node, integrating over topology, for L≈7, C≈100–200,
no calibrations.

1. **Make the cross-family rate prior the centrepiece and validate it by
   leave-one-family-out (LOFO).** §§3.3 and 4.2–4.3. It is the dominant
   uncertainty (≈0.3 log-SD vs ≈0.07 from data). Build the meta-analysis with
   σ_f as the reported floor, add family covariates (population size first),
   and publish a LOFO coverage table and PIT histogram. *C99:* the engine
   only consumes a prior. The fitting is Python. **Highest value.**
2. **Exhaustive topology enumeration with Laplace/INLA evidences, refined by
   AGHQ.** §§5e, 5f and 5j. 10,395 topologies at L=7 and 135,135 at L=8 take
   minutes in C. It gives true BMA, MRCA-of-named-set ages and
   $P(\text{monophyly})$, with no tree MCMC to tune. *C99:* excellent.
3. **Report lexical age, relative age and calendar age separately, and use
   the rate-prior convolution identity.** §§3.3 and 4.1. It is exact under a
   scale-invariant age prior and makes rate-prior sensitivity free. It also
   shows users that relative chronology is much sharper than absolute.
   *C99:* trivial.
4. **Add branch×concept heterogeneity (+Γ via Lomax branch survival) and
   compare mixing families.** §2.4. It fixes a 5–15% deep-node young bias
   that per-concept rates do not remove. One parameter, recursion unchanged.
   *C99:* trivial.
5. **Posterior predictive checks with population-genetic statistics.**
   §§1.4 and 9: $F_{ij}$ residuals, the $K_c$ distribution via the exact
   $\mathbb E[K]$ formula, and triplet identities. These detect borrowing,
   lumping and splitting, which are the biases that dates cannot reveal
   otherwise. *C99:* easy.
6. **Coverage-under-misspecification simulation suite** (SBC plus a richer
   "world" simulator with borrowing, coding error, Γ heterogeneity and
   pulses). §§6.1–6.2. It turns the bias tables in §§2, 8 and 9 into
   empirical coverage for clepsydra's own intervals. *Python side.*
7. **Exact marginalisation of coding error for small L** (latent partitions,
   split/lump rates with priors from double-coding studies) plus a
   tip-offset parameter. §8. Shallow nodes are where split error dominates.
   *C99:* moderate (Bell(L) enumeration per meaning, ≤877 at L=7).
8. **Power-scaling prior sensitivity** on stored draws (rate prior, σ_b,
   tree prior, contamination). §6.4. Cheap, standard and reportable.
   *C99 or Python.*
9. **Ratio + log-height reparametrisation of the sampler** (and HMC later
   via the linear-time gradient). §§5h and 1.3. It decorrelates scale from
   shape, which is the posterior's natural geometry under zero calibrations.
   *C99:* moderate.
10. **Pairwise composite likelihood with Godambe adjustment, plus BIONJ/LSD
    point estimates** as a fast independent second engine and a starting
    point. §§5a–5b. About 88% efficient. Note the 2.2–2.4× overconfidence of
    any naive pairwise product: if pairwise posteriors are ever combined
    over a tree, they *must* be curvature-adjusted. *C99:* easy.

Not recommended as primary engines: variational inference (underestimates
variance, which is the deliverable), and neural SBI (unfaithful under
misspecification. Use it as a validator only). Worth keeping in reserve:
nested sampling for model evidences on the top topologies, SMC when
families exceed ~9 languages, and one- or two-reticulation network mixtures
if posterior-predictive checks show borrowing.

---

## Appendix A. Scripts (all in `docs/research/math/`, run with `python3 -I`)

| Script | What it computes |
|---|---|
| `e1_pairinfo.py` | Per-meaning Fisher info $x^2/(e^x-1)$, its maximum, and pair CRBs with heterogeneous rates and a family-rate prior |
| `e2_tree_fisher.py` | Exact T1 likelihood on a 7-tip tree by enumerating the 877 partitions (corrected recursion); Bayesian information matrix for ages, branch multipliers and family rate; output in `out_e2_*.txt` |
| `e3_hetero.py` | Gamma/lognormal mixture retention, naive-date bias, shallow-calibration extrapolation, ratio bias |
| `e4_bias.py` | Pair-MLE bias under branch×concept Γ, split/lump coding error, borrowing, split pulses |
| `e5_composite.py` | Godambe efficiency of pairwise composite likelihood vs exact T1 |
| `e6_topo.py` | Topology counts |
| `e7_watterson.py` | Exact $\mathbb E[K]$ formula vs simulation; Starostin–Lambert-W curve |

Assumptions in all experiments: mean concept rate 0.15 /ky per lineage
(retention ≈0.86 per millennium), concept log-SD 0.8, strict clock unless
stated, and Fisher-information (large-$C$) approximations. The numbers are
lower bounds on posterior SD under a correctly specified model. Real
posteriors are wider by the misspecification terms in §§2, 8, 9 and 10.

## Appendix B. Additional sources consulted (verified this session)

- Kelly & Nicholls 2017: https://doi.org/10.1214/17-AOAS1040
- Zhu, dos Reis & Yang 2015: https://pmc.ncbi.nlm.nih.gov/articles/PMC4380039/
- Holman et al. 2011: https://pure.mpg.de/pubman/item/item_2395214_2/component/file_2432001/shh768.pdf
- Rama & Wichmann 2020: https://pmc.ncbi.nlm.nih.gov/articles/PMC7423060/
- Greenhill, Currie & Gray 2009: https://pmc.ncbi.nlm.nih.gov/articles/PMC2677599/
- Neureiter et al. 2022: https://ideas.repec.org/a/pal/palcom/v9y2022i1d10.1057_s41599-022-01211-7.html
- Hermans et al. 2022: https://arxiv.org/abs/2110.06581
- Cannon, Ward & Schmon 2022: https://arxiv.org/abs/2209.01845
- Ji et al. 2020: https://arxiv.org/abs/1905.12146 ; Ji et al. 2023: https://pmc.ncbi.nlm.nih.gov/articles/PMC10636426
- Kallioinen et al. (priorsense): https://cran.r-project.org/package=priorsense
- Rönchen, Wiklund & Hammarström 2024 appendix: https://brill.figshare.com/articles/journal_contribution/Appendix_with_the_article_Likelihood_calculation_in_a_multistate_model_of_vocabulary_evolution_for_linguistic_dating_by_Philipp_R_nchen_Tilo_Wiklund_and_Harald_Hammarstr_m_published_in_LDC_14_2024_/26089507/1
- Starostin formula (secondary): https://www.eva.mpg.de/lingua/conference/09_SwadeshCentenary/pdf/Handouts/Vasiljev_Handout.pdf
- van der Merwe 1966 (citation via JSTOR TOC): https://www.jstor.org/stable/i327882

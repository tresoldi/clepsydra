# Model specification

> Milestone 0 deliverable. Defines every parameter, unit, support and prior;
> gives a generative description of each likelihood; states how synonyms,
> missing data and incompatible concepts enter each model; identifies
> nonidentifiable combinations; defines the exact pairwise posterior and the
> exact fixed-tree MCMC target; gives pseudocode for the strict novel-state
> likelihood; and lists provisional assumptions requiring later simulation.
>
> Notation follows the product brief (`Hierarchical Bayesian Lexical Dating
> Specification.pdf`) and, for the novel-state model, Rönchen, Wiklund &
> Hammarström 2024 (*Language Dynamics and Change* 14(1):1–41,
> doi:10.1163/22105832-bja10032), independently re-derived here. No code is
> copied from any external implementation.

## Contents

1. [Time, ages and path lengths](#1-time-ages-and-path-lengths)
2. [Parameter table](#2-parameter-table)
3. [Models](#3-models)
   - 3.1 [P0 — classical homogeneous baseline](#31-p0--classical-homogeneous-baseline)
   - 3.2 [P1 — hierarchical pairwise cognacy model](#32-p1--hierarchical-pairwise-cognacy-model)
   - 3.3 [T1 — strict novel-state multistate model](#33-t1--strict-novel-state-multistate-model)
   - 3.4 [T2 — binary cognate-set model](#34-t2--binary-cognate-set-model)
   - 3.5 [T3 — finite-state Mk/CTMC (reserved)](#35-t3--finite-state-mkctmc-reserved)
4. [Clocks and rate variation](#4-clocks-and-rate-variation)
5. [Synonyms, missing data, loans, incompatible concepts](#5-synonyms-missing-data-loans-incompatible-concepts)
6. [Contamination](#6-contamination)
7. [Identifiability analysis](#7-identifiability-analysis)
8. [Exact pairwise posterior](#8-exact-pairwise-posterior)
9. [Exact fixed-tree state and MCMC target](#9-exact-fixed-tree-state-and-mcmc-target)
10. [Strict novel-state likelihood — pseudocode](#10-strict-novel-state-likelihood--pseudocode)
11. [Independent brute-force validation — pseudocode](#11-independent-brute-force-validation--pseudocode)
12. [Provisional assumptions requiring later simulation](#12-provisional-assumptions-requiring-later-simulation)

---

## 1. Time, ages and path lengths

All internal times are **numeric years before a configurable origin**
(default `time_origin_year = 1950`, i.e. "BP"). The C engine never parses
calendar strings; calibrations supply numeric years.

- **Age** `a(v)` of a node `v`: years before origin. Larger = older. The root
  is the oldest node. Extant tips have age 0 (or a small positive age for
  ancient/attested tips). Constraint on every edge `(parent p → child c)`:
  `a(p) > a(c)` whenever the edge has positive duration.
- **Branch duration** `t_b = a(parent_b) − a(child_b) > 0` (years).
- **Path length** between two nodes: sum of branch durations along the unique
  path. For a pair of tips with MRCA age `T` and tip ages `a_1, a_2`:
  `d = (T − a_1) + (T − a_2) = 2T − a_1 − a_2`, with `T ≥ max(a_1, a_2)`.
- **Lexical time vs calendar time.** Rates have units *replacements per slot
  per year* (calendar). If no absolute calibration or clock-scale prior is
  supplied, only the product `r·d` (dimensionless "lexical time") is
  identifiable; outputs then report *lexical-time units* and state that
  absolute years are not identifiable. See §7.

Rates are stored internally on the **log scale** as `log r` to keep them
positive and to make the hierarchy additive (§4).

---

## 2. Parameter table

| Symbol | Meaning | Unit | Support | Prior (default; all configurable) |
|---|---|---|---|---|
| `a(v)` | age of node `v` | years BP | `≥ 0`, tree-ordered | from calibrations / root prior |
| `T` | MRCA age of a target pair | years BP | `≥ max(a_1,a_2)` | calibrated or root prior |
| `r_c` | replacement rate of concept `c` | 1/year | `> 0` | via `log r_{f,c}` (below) |
| `μ_c` | global reference log-rate for concept `c` | log(1/year) | ℝ | `N(m_c^ref, (s_c^ref)²)` |
| `β_f` | family-wide log clock multiplier for family `f` | log-scale | ℝ | `N(0, σ_f²)` |
| `u_{f,c}` | family–concept interaction (deviation) | log-scale | ℝ | `N(0, τ²)` |
| `σ_f` | SD of `β_f` across families | log-scale | `> 0` | `HalfNormal(0, 0.50)` |
| `τ` | SD of `u_{f,c}` | log-scale | `> 0` | `HalfNormal(0, 0.50)` |
| `m_b` | relaxed branch multiplier (branch `b`) | multiplier | `> 0` | `log m_b ~ N(−σ_b²/2, σ_b²)` so `E[m_b]=1` |
| `σ_b` | SD of `log m_b` | log-scale | `> 0` | `HalfNormal(0, 0.30)` (shrinkage) |
| `r_{f,c}` | realized concept rate in family `f` | 1/year | `> 0` | `exp(μ_c + β_f + u_{f,c})` |
| `r_{b,f,c}` | realized concept rate on branch `b` | 1/year | `> 0` | `r_{f,c} · m_b` |
| `ε` | global contamination probability | prob | `(0,1)` | `Beta(1, 19)` |
| `ε_c` | concept contamination probability | prob | `(0,1)` | `Beta(a_ε, b_ε)` (hierarchical, later) |
| `q_c` | match prob. under noise model (concept `c`) | prob | `[0,1]` | empirical / Dirichlet-multinomial |
| `y_c` | observed cognacy probability (pairwise) | prob | `[0,1]` | **input**, not learned |
| `λ_g, λ_l` | binary CTMC gain / loss rates (T2) | 1/year | `> 0` | weak lognormal / hierarchical |
| `w_{i}` | primary weight of form `i` in a cell | weight | `[0,1]`, sums to 1 | **input** |
| `π_obs, π_h` | T2 equilibrium freqs (present/absent) | prob | simplex | prior / fixed uniform |

**Legacy scale anchors (fallback presets, not universal constants):**

| Preset | Unit | log_mean | log_sd | Meaning |
|---|---|---|---|---|
| `legacy_swadesh100` | per_kyr | −1.8916490462 | 0.80 | 0.86 retention/kyr (`r = −ln 0.86 = 0.1508228897`/kyr) |
| `legacy_lees200` | per_kyr | −1.5271140564 | 0.80 | 0.8048 retention/kyr (`r = −ln 0.8048 = 0.2171614796`/kyr) |

The wide SD (0.80) is intentional; these are weak fallback anchors. They are
selected only with an explicit `--model classical --legacy-clock ...` flag.

**Provisional default hyperpriors** (clearly marked provisional; sensitivity
presets `tight` / `standard` / `wide` must all be run for any published
result):

| Hyperparameter | Default |
|---|---|
| global clock log-SD | 0.80 |
| minimum concept prior log-SD | 0.35 |
| family clock SD | `HalfNormal(0, 0.50)` |
| family–concept SD | `HalfNormal(0, 0.50)` |
| branch-clock SD | `HalfNormal(0, 0.30)` |
| global contamination | `Beta(1, 19)` |
| unknown concept log-rate SD | 1.00 |

---

## 3. Models

A common vtable interface lets each model expose `prepare`, `log_likelihood`,
`log_likelihood_concept`, `simulate`, `destroy` (see `include/gd.h`). All
models share the time/age definitions of §1 and the clocks of §4.

### 3.1 P0 — classical homogeneous baseline

Regression baseline only. One rate `r` for all concepts. For a concept with
total path length `d`:

```
P(same inherited cognate) = exp(−r·d).
```

For two contemporaneous daughters diverging `t` years ago, `d = 2t`. The
classical Swadesh-100 retention 0.86/kyr gives `r = −ln 0.86`/kyr; Lees-200
gives `r = −ln 0.8048`/kyr (§2). Selection of a legacy clock is **never**
silent: `--model classical --legacy-clock swadesh-100|lees-200`.

P0 supports no contamination, no hierarchy, no synonyms. It exists for
analytical regression tests against the closed-form glottochronological
posterior (§8, `docs/validation.md`).

### 3.2 P1 — hierarchical pairwise cognacy model

Two languages, concepts `c = 1..C`, family `f`. For concept `c` and total
path length `d`:

```
s_{f,c}(d) = exp(−r_{f,c}·d),        log r_{f,c} = μ_c + β_f + u_{f,c}.
```

Hierarchy (defaults in §2):

```
μ_c     ~ N(m_c^ref, (s_c^ref)²)
β_f     ~ N(0, σ_f²)
u_{f,c} ~ N(0, τ²)
σ_f, τ  ~ HalfNormal(0, s)
```

The family data **updates** the reference priors: `m_c^ref, s_c^ref` come
from the loaded prior package (reference mode `global`/`family_specific`/
`leave_current_family_out`), and `β_f, u_{f,c}` let the observed family shift
the log-rate. With no calibration or clock-scale prior, only `r_{f,c}·d`
(lexical time) is identifiable; the absolute scale `r` vs `d` is not (§7).

**Uncertain cognacy.** Each concept supplies an observed cognacy probability
`y_c ∈ [0,1]` (expert/upstream uncertainty, **not** learned). With
`p_c = P(observed match | T, r_c, ε)`,

```
L_c = y_c·p_c + (1 − y_c)·(1 − p_c).
```

**Pairwise contamination** (§6):

```
p_c = (1 − ε_c)·exp(−r_c·d) + ε_c·q_c,
```

with `contamination ∈ {none, global, concept}`. Concepts are conditionally
independent given parameters, so `L = ∏_c L_c`.

**Synonyms (pairwise).** If a cell has several candidate forms with primary
weights `w_i` summing to 1, the match probability is marginalized exactly
over candidate pairs (§5.1). No MCMC is needed for P1 in normal low
dimensions; the posterior is obtained by deterministic integration (§8,
`docs/inference.md`).

### 3.3 T1 — strict novel-state multistate model

The model most closely related to Rönchen et al. 2024, independently
re-implemented. One model module among several.

**Generative process (per concept `c`, rate `r_{b,f,c}` on branch `b`).**
Let the tree be rooted with node set `V`, edges `E`. For an edge `e` of
duration `t_e` and concept rate `r_e`:

```
ch(e) = 1 − exp(−r_e·t_e)     # prob of ≥1 replacement on e
ret(e) = exp(−r_e·t_e)        # prob of no replacement on e
```

1. Assign the root an arbitrary cognate class.
2. Downward, for each child `v` of `u` on edge `(u,v)`: with prob `ch(u,v)`
   assign `v` a **new** cognate class never seen before in the tree; with
   prob `ret(u,v)` assign `v` the same class as `u`.
3. Continue to the leaves.

Consequences (model assumptions, **not** project-wide truths):
- one primary cognate state per language–concept cell;
- every replacement creates a class never used before;
- no parallel innovations, no back mutation, no reappearance of extinct
  classes;
- unobserved ancestral states are permitted (a class may exist only at
  internal nodes);
- meanings evolve independently;
- inheritance is strictly tree-like (no undetected borrowing).

**Data, compatibility, path property, extension.** For concept `m`, data is
a partial function `f_m : V ⇀ ℕ` (cognate-class labels; labels are arbitrary
unique IDs, only equality structure matters). An outcome `f : V → ℕ` is
*compatible* with `f_m` iff one of (equivalent): (i)
`f_m(u)=f_m(w) ⟺ f(u)=f(w)` for all observed `u,w`; (ii) some `f'~f`
extends `f_m`; (iii) some `f_m'~f_m` is extended by `f`.

- **Path property**: `g` has it iff whenever `g(u)=g(u')`, every node on the
  path `u…u'` is in `dom(g)` and maps to the same class. (Prop 3 of the
  paper: a total `f` has `P([f])>0` iff it has the path property.)
- **Determined node**: `v` is determined to class `n` if some pair
  `u_1,u_2 ∈ dom(f_m)` with `f_m(u_1)=f_m(u_2)=n` has `v` on the path
  between them.
- **Extendability**: `f_m` is extendable iff no node is determined to two
  different classes. If not extendable, `P(f_m)=0` (data is incompatible).
- **Extension** `f̃_m`: assign every determined node its determined class.
  `f̃_m ⊇ f_m`, has the path property, and `P(f_m) = P(f̃_m)` when extendable.

**Likelihood recursion (binary tree; nonbinary generalization below).**
For node `v`, subtree `T_v`:

```
p_v = P(F compatible with f̃_m ↾ T_v)              # subtree likelihood
q_v = P(F compatible with f̃_m ↾ T_v  AND  F preserves a data value up to v)
```

`PresData(f, v)` holds iff some determined `w ∈ T_v` has `f(v)=f(w)`. Leaves:
`p_v=1`; `q_v=1` iff `Det(v)` else `0`. For an internal `u` with children
`v_1,v_2`, `ch_i=ch(u,v_i)`, `ret_i=ret(u,v_i)`, and `Det(v)` meaning
`v ∈ dom(f̃_m)`:

- **Case i** (`¬Det(u)`):
  ```
  p_u = p_{v1}·p_{v2} − q_{v1}·q_{v2}·ret_1·ret_2
  q_u = q_{v1}·ret_1·(p_{v2} − q_{v2}·ret_2)
      + q_{v2}·ret_2·(p_{v1} − q_{v1}·ret_1)
  ```
- **Case ii** (`Det(u)`; `q_u = p_u`), let `prob_i = ret_i` if
  `f̃_m(v_i)=f̃_m(u)` else `ch_i`:
  - **iia** (`Det(v_1)` and `Det(v_2)`): `p_u = p_{v1}·prob_1·p_{v2}·prob_2`
  - **iib** (`¬Det(v_1)` and `¬Det(v_2)`):
    `p_u = (p_{v1} − q_{v1}·ret_1)·(p_{v2} − q_{v2}·ret_2)`
  - **iic** (WLOG `Det(v_1)`, `¬Det(v_2)`):
    `p_u = p_{v1}·prob_1·(p_{v2} − q_{v2}·ret_2)`

The likelihood is `p_root`. Complexity: O(n log n) for the extension, O(n)
for the recursion, linear in the number of concepts. All quantities stay in
`[0,1]`; computed in **log space** internally with `logsumexp`/`log1mexp`
helpers (see `docs/inference.md` §numerical stability) to avoid underflow;
the subtraction `p_{v1}·p_{v2} − …` is evaluated stably as
`logdiffexp` of the two log terms.

**Nonbinary trees.** The binary recursion generalizes to polytomies via the
paper's Theorems A.8/A.9 (appendix supplement, doi:10.6084/m9.figshare.26089507).
> ⚠️ **Open item:** that appendix supplement is **not** in our `literature/`
> folder. Until it is obtained, the C implementation will (a) support binary
> trees fully, and (b) handle polytomies by **deterministically resolving**
> them into zero-length ladders (a documented approximation) **OR** by
> brute-force enumeration on small polytomies in tests. The exact nonbinary
> recursion will be implemented once the appendix is supplied. This is logged
> in `docs/paper-evaluation.md` and §12.

**Ladder resolution as implemented (M6).** `gd_t1_ws_create` rewrites each
node with `k > 2` children into a ladder of `k − 2` inserted nodes that carry
the polytomy's *own age*, so every inserted branch has duration zero and every
original child keeps its own branch duration. Under the generative process
above this rewriting is not an approximation at all: each of the `k` children
still, independently, either keeps the parent's class with probability
`exp(−r·t_child)` or opens a class new to the whole tree, which is exactly the
polytomy's own joint distribution. The inserted nodes are unobserved and
always carry the polytomy's class. `tests/likelihood/test_t1_likelihood.c`
confirms this numerically: for 3-way and 4-way polytomies (including one below
the root with an ancient tip), the ladder likelihood agrees with brute-force
enumeration over the **unresolved** polytomy for every partition of the tips,
to `1e-9`.

That is a statement about the model of §3.3, not a claim about the paper. If
Rönchen's Theorems A.8/A.9 give a polytomy a different reading — e.g. an
unresolved binary split rather than `k` independent edges — then that is a
different model and this section must be revisited when the appendix arrives.
The appendix therefore remains an open item (§12 item 10), and the exact
`O(k)` nonbinary recursion is still unimplemented; the ladder costs `k − 2`
extra nodes instead.

**Robust contamination wrapper** (§6): `L_c = (1−ε_c)·L_{tree,c} + ε_c·L_{noise,c}`.
The user chooses `incompatible ∈ {error, exclude, contamination}` (default
`contamination`); incompatible concepts never silently disappear. The
posterior `P(concept c contaminated)` is reported.

**Synonyms in T1** (§5.1): one latent primary state per cell; in pairwise
mode marginalize exactly over candidate pairs; in tree mode sample the latent
primary state within MCMC, recomputing only the affected concept likelihood.

**Simulation** from T1 is the generative process above; used for SBC and
posterior-predictive checks (`docs/validation.md`).

### 3.4 T2 — binary cognate-set model

Each observed cognate class is a binary trait (present/absent) per language.
Naturally permits multiple simultaneous cognate classes for one meaning
(synonymy/polymorphism) — an alternative to T1's single-primary-state
assumption.

**Asymmetric binary CTMC** per cognate class `g`, on branch `b`:

```
0 --λ_g--> 1      (gain)        1 --λ_l--> 0      (loss)
```

Transition matrix exponentiated by `Q·t_b` (or the closed-form 2×2 solution)
and pruned with Felsenstein's algorithm. Rates `λ_g, λ_l` may be concept- and
family-scaled via the same hierarchy as T1 (§4).

**Ascertainment correction.** All-zero (absent-everywhere) cognate classes
are unobservable in the input, so the likelihood is conditioned on
observability:

```
P(x | observable) = P(x) / (1 − P(all-absent))
```

Following King 2026, **partition rates must be weighted by the number of
meanings, not the number of cognate classes**, to avoid the
ascertainment × reweighting miscalibration he documents
(see `docs/paper-evaluation.md`). This is the default; the cognate-count
weighting is available only as an explicit opt-in for reproducing old
analyses.

**Warnings** (enforced in output):
- cognate sets within one meaning are **not** truly independent;
- meanings with many cognate classes contribute disproportionate likelihood;
- raw marginal likelihoods from binary and multistate encodings are **not**
  directly comparable. Cross-encoding model comparison uses predictions in a
  shared observable space (held-out pairwise cognacy / held-out cell states),
  never raw AIC/BIC/WAIC across codings (§9, `docs/validation.md`).

**T2 as implemented (M7).** `gd_model_binary.{c,h}` holds the T2 likelihood.
`gd_t2_transition_logs` is the closed-form 2×2 CTMC solution
`P_ij(t) = π_j + (δ_ij − π_j)e^{−μt}` with `μ = λ_g + λ_l` and stationary
`π = (λ_l, λ_g)/μ`, in log space; `gd_t2_loglik_class` is Felsenstein pruning
over an iterative postorder, storing `log L_v(0)`, `log L_v(1)` per node with a
`logsumexp` over the two child states. Because pruning combines a node's
children directly, **a polytomy needs no resolution — there is no ladder in
T2** (contrast §3.3). A `gd_t2_ws` workspace is built once per topology and
re-reads node ages every call (it survives an MCMC chain); the pruning
allocates nothing. Observed internal nodes (ancient languages) and missing
cells are handled uniformly: a missing node is marginalized (both child
messages kept), an observed one fixes its state.

The ascertainment correction divides each class likelihood by
`1 − P(absent at every observed node)`, the all-absent term computed by the
same pruning with every observed node forced to state 0 and missing nodes still
marginalized. This makes the ascertainment-corrected likelihood a distribution
over exactly the **observable** tip patterns (all patterns except all-absent),
which `tests/likelihood/test_t2_likelihood.c` confirms sums to one; the raw
(uncorrected) likelihood sums to one over all patterns. Both are checked
against an **independent** brute force that enumerates all `2^n` node-state
assignments — sharing nothing with the pruning recursion — on binary trees,
polytomies, missing/observed-internal patterns, a fixed root prior, and
randomized trees, to `1e-9`. Simulation (`gd_t2_simulate_class`) draws the root
from the prior and propagates the CTMC downward; its observable tip-pattern
frequencies match the corrected likelihood, and a profile-likelihood refit
recovers the rate scale.

Per King 2026, `gd_t2_partition_weights` defaults to **meaning weighting**
(`GD_T2_WEIGHT_MEANING`: every meaning counts once, regardless of its cognate
count); cognate-count weighting is opt-in only. `gd_t2_loglik_dataset` scales a
meaning's `λ_g, λ_l` by its concept rate and raises
`within_meaning_dependence_warning` whenever any meaning contributes more than
one cognate class (so the false within-meaning independence assumption is
actually being relied on); `cross_encoding_incomparable` is a standing 1. The
relaxed/hierarchical rate scaling of §4 and the robust contamination wrapper of
§6 are applied by the caller and the sampler (M8/M9), not baked into the class
likelihood.

A stochastic-Dollo-like birth-and-loss module is reserved for later; it does
not block the initial Bayesian system.

### 3.5 T3 — finite-state Mk/CTMC (reserved)

Interface reserved for a finite-state multistate CTMC that permits homoplasy
(relaxing T1's no-parallel-development assumption). **Not implemented** before
T1 and T2 are validated. Its role is robustness comparison, not "more flexible
= more realistic."

---

## 4. Clocks and rate variation

A progression; each is a strict special case of the next.

1. **Strict clock**: `r_{b,f,c} = r_{f,c}`. Baseline tree model.
2. **Family clock multiplier**: `r_{f,c} = exp(β_f)·r_c`. Relevant when priors
   were trained across families.
3. **Relaxed branch clock**: `r_{b,f,c} = r_{f,c}·m_b`,
   `log m_b ~ N(−σ_b²/2, σ_b²)` so `E[m_b]=1`, shrinkage prior on `σ_b`.
   Strict clock = `σ_b=0`. Branch multipliers and branch durations are highly
   confounded (§7); require calibration-sensitivity diagnostics (the
   `clepsydra calib-sensitivity` command: leave-one-calibration-out and
   window-widening re-dating, `docs/cli-manual.md`) and warn when `σ_b` is
   weakly identified.
4. **Split-pulse model** (reserved): `H_{b,c} = r_{b,f,c}·t_b + δ_v` (an
   extra lexical-change pulse at a split). **Excluded from v1** unless
   simulation demonstrates identifiability (§12).

The concept hierarchy `log r_{f,c} = μ_c + β_f + u_{f,c}` is shared by P1, T1
and T2; P0 uses a single global `r`.

**As implemented (M9).** The relaxed branch clock (case 3) and the per-concept
deviation `μ_c` are inference blocks in the tree sampler (`gd_mcmc.c`), reaching
the likelihood through `gd_t1_loglik_concept_relaxed`, whose per-branch rate is
`exp(logr + μ_c) · m_b` (the multiplier `m_b` is indexed by the branch's child
node; inserted polytomy-ladder edges are zero-length, so their multiplier never
matters). `σ_b` carries the HalfNormal shrinkage hyperprior and is updated
against the branch-multiplier prior alone (no likelihood recompute). Because
`r·m_b·t_b` is the only per-branch identifiable combination, branch multipliers
and node ages are confounded (§7.2); the sampler reports the branch-multiplier /
branch-duration posterior correlation and flags a strong trade-off. The **family
multiplier `β_f`** (case 2) is meaningful only across families: on a single
fixed tree it is absorbed into the global clock scale, so it is documented as a
reduction here and activated with cross-family reference priors in M11. The
split-pulse model (case 4) stays excluded from v1 (§12).

---

## 5. Synonyms, missing data, loans, incompatible concepts

### 5.1 Synonyms and polymorphism

Input retains **all** forms per language–concept cell. Each form carries a
`primary_weight`, `cognacy_confidence`, `inheritance_probability`,
`loan_probability`, `status`. Weights apply within a cell; non-missing
primary weights must sum to ~1 (validated; if absent, set
`1/n_candidates` with a warning).

- **T1 (strict multistate)**: one latent primary cognate state per cell.
  - *Known primary*: one form has `primary_weight=1`.
  - *Uncertain primary*: weights sum to 1. In **pairwise** mode, marginalize
    exactly over candidate pairs. In **tree** mode, treat the primary state
    as a latent categorical, sample it in MCMC, recompute only the affected
    concept. **No** combinatorial enumeration over synonym configurations.
- **T2 (binary)**: all cognate classes present in the cell are coded
  present; no form is discarded for not being primary.
- **Genuine coexistence** (birth/competition/narrowing/loss of synonyms) is
  **not** modeled in v1; documented as a limitation. Enough input data is
  preserved to support a future model. The paper's suggested `k`-forms-per-
  meaning extension (store probabilities of preserving 0..k classes) is
  recorded as a postponed idea (`docs/paper-evaluation.md`).

### 5.2 Missing data

"Missing" = unknown, **not** absence. A missing cell contributes no direct
lexical observation. Distinguished statuses (preserved in diagnostics though
v1 likelihood treats all non-observed as missing): `missing`,
`concept absent/not lexicalized`, `form unknown`, `form intentionally
excluded`.

In T1, missing leaves/internal observations are marginalized by the
extension + recursion (a node not in `dom(f_m)` is simply not constrained).
In T2, missing cells contribute `1` (the trivial factor) to the pruning
partial likelihood.

### 5.3 Known and suspected loans

`loan_mode ∈ {exclude, mixture, sensitivity}`.
- **exclude**: forms marked loans above a threshold are treated as missing
  in the vertical-inheritance likelihood.
- **mixture**: `loan_probability` increases the cell's probability of
  belonging to the contamination component (§6).
- **sensitivity**: run several analyses (e.g. exclude `p≥0.50`, exclude
  `p≥0.80`, include all with contamination).

Explicit borrowing networks are **not** implemented in v1.

### 5.4 Alternative cognacy analyses

A dataset may carry several named codings (`analysis_id = expert_a`,
`expert_b`, `automated_01`). The CLI runs them separately or samples among
complete alternative codings. C does not infer arbitrary cognate partitions
in v1.

### 5.5 Incompatible concepts

`incompatible ∈ {error, exclude, contamination}` (default `contamination`).
In T1, an unextendable `f_m` gives `P(f_m)=0`; under `contamination` the
mixture (§6) keeps the dataset likelihood nonzero and reports
`P(c contaminated)`. Under `error` the run aborts with a diagnostic; under
`exclude` the concept is dropped with a recorded warning. **No** automatic
wholesale concept removal.

---

## 6. Contamination

Robust contamination component so incompatible concepts do not zero the whole
likelihood:

```
L_c = (1 − ε_c)·L_{tree,c} + ε_c·L_{noise,c}.
```

- **Global**: `ε ~ Beta(a_ε, b_ε)` (default `Beta(1,19)`).
- **Concept-specific** (later): `ε_c` with hierarchical shrinkage toward the
  global `ε`; identifiability unverified (§12).
- **Noise model** `L_{noise,c}`: an explicit, documented statistical error
  model — a symmetric Dirichlet–multinomial or empirical categorical over
  observed states. It is **not** a claim that the real event was random.

**Noise model as implemented for T1 (M6).** The observable of a concept under
T1 is the *partition* of the observed nodes into cognate classes (the label
integers are arbitrary), and `L_{tree,c}` is a probability distribution over
exactly those partitions. The noise component must live on the same space or
the mixture is not a distribution at all, so it is a Chinese-restaurant
process with concentration `α` (default `1.0`):

```
L_noise = α^B · Π_b Γ(n_b) / [Γ(α+n)/Γ(α)]
```

for a partition of `n` observed nodes into `B` blocks of sizes `n_b`. This is
the `K → ∞` limit of the symmetric Dirichlet–multinomial named above, and the
form of it that is exchangeable over nodes and normalized over partitions
rather than over labelings. Both properties are checked numerically
(`test_t1_noise_model_normalized`, `test_t1_mixture_normalized`: the mixture
sums to 1 over all partitions of the observed nodes, including the ones that
are not extendable). Adequacy of this choice against alternatives remains
§12 item 1 — normalization is not adequacy.

The posterior probability that each concept came from the contamination
component is reported in `posterior_concepts.tsv`. In `mixture` loan mode,
`loan_probability` feeds `ε_c` for that cell.

---

## 7. Identifiability analysis

1. **Absolute scale vs clock rate.** Without a dated tip, node calibration,
   or clock-scale prior, the likelihood depends only on `r·d` (dimensionless
   lexical time). Scaling all rates by `k` and all durations by `1/k` leaves
   the likelihood invariant. ⇒ Output reports lexical-time units and states
   "absolute years not identifiable" unless an anchor is supplied. (Cf.
   Maurits et al. 2020 on the clock/age invariance.)
2. **Branch multiplier vs branch duration.** `r_{b,f,c}·t_b = r_{f,c}·m_b·t_b`:
   `m_b` and `t_b` are confounded per branch. Identifiable only through (a)
   calibrations pinning some node ages and (b) the shrinkage prior on `σ_b`.
   Warn when `σ_b` posterior is weakly identified; require calibration
   leave-one-out sensitivity.
3. **`β_f` vs `μ_c` vs overall scale.** `μ_c + β_f + u_{f,c}`: a constant
   added to all `β_f` can be absorbed into the global scale; anchored by the
   reference priors on `μ_c` and by calibrations.
4. **Concept rate vs root age** (tree): posterior correlation reported;
   concepts with posterior ≈ prior are flagged as uninformative.
5. **Contamination vs rate.** High `ε_c` can mask a mis-specified rate;
   report `P(c contaminated)` and prior-to-posterior gain.
6. **Nodes determined mainly by one calibration** are flagged.

Diagnostics (§12.2 of the brief): absolute-scale-anchored flag;
posterior correlations (clock–root, branch-mult–duration); prior≈posterior
concepts; single-calibration-dominated nodes; bound-mass; prior-to-posterior
information gain.

---

## 8. Exact pairwise posterior

For `date-pair` (P0/P1), deterministic integration is preferred over MCMC.

**P0 closed form.** With a single rate `r` and a prior `π(r)`, the posterior
over divergence time `T` (path length `d=2T−a_1−a_2`) for concepts with
match indicators `x_c ∈ {0,1}` and contamination off:

```
p_c(T) = exp(−r·d)            if x_c = 1 (match)
p_c(T) = 1 − exp(−r·d)        if x_c = 0 (mismatch)
L(T, r) = ∏_c [ y_c p_c + (1−y_c)(1−p_c) ]     # uncertain cognacy
posterior(T, r) ∝ L(T, r) · π(r) · π(T)
```

`π(T)` incorporates `T ≥ max(a_1,a_2)` and any pair calibration. The marginal
posterior of `T` is `∫ L(T,r) π(r) π(T) dr`. We compute it by:

- adaptive grid refinement in `T` (and in `log r` if `r` is not fixed);
- log-space normalization (`logsumexp`);
- numerical CDF construction from the grid;
- adaptive Simpson / Gauss–Kronrod where suitable;
- Gauss–Hermite quadrature for lognormal rate integration **if validated**
  against dense grids (else fall back to grid; §12).

**P1 hierarchy.** With `log r_{f,c} = μ_c + β_f + u_{f,c}`, the latent
vector `(μ_c, β_f, u_{f,c}, ε)` is integrated. For a single family pair
`β_f` and `u_{f,c}` enter linearly; deterministic integration is practical up
to a few latent dimensions (grid / quadrature product). **Switch to MCMC
only** when pairwise mode includes enough latent variables that deterministic
integration is no longer practical (e.g. many concept-specific `u_{f,c}`
simultaneously). The decision is recorded in the run manifest.

**Outputs**: evaluated posterior grid (if requested), and summaries — mean,
median, sd, mode (where stable), 50/80/90/95% equal-tail, 80/90/95% HPD,
threshold probabilities `P(T > τ)`, prior-to-posterior comparison,
prior-dominance diagnostic.

---

## 9. Exact fixed-tree state and MCMC target

For `date-tree` (T1/T2) and `date-pair-context`, the target is the Bayesian
posterior over the fixed rooted topology:

```
posterior(θ, z) ∝ L(data | θ, z, tree) · π(θ)
```

where **parameters** `θ` = {internal node ages `a(v)`; global clock scale;
family multipliers `β_f`; concept log-rates/deviations `μ_c, u_{f,c}`;
relaxed branch multipliers `m_b`; contamination `ε` (and `ε_c` later);
binary CTMC rates `λ_g, λ_l` for T2; uncertain calibration hyperparameters
when explicitly modelled}, and **latent variables** `z` = {primary synonym
state per uncertain cell (T1)}.

Constraints: `a(parent) > a(child)` on every positive-duration edge;
calibration distributions act as priors on the named nodes/tips (see
`docs/calibrations.md`); `T ≥ max(a_1,a_2)` for any target pair.

Sampler: **Metropolis-within-Gibbs** (dependency-free). Blocks and updates
(per `docs/inference.md`): log-scale random walks for positive rates;
logit-scale for bounded proportions; age proposals constrained by parent/child
ages; block scaling of a clade's internal ages; slice sampling for 1-D scale
parameters where stable; categorical Gibbs/Metropolis for latent primary
forms (recompute only the affected concept). Adaptation **only during
warm-up**; frozen before retained samples. Per-chain deterministic seeds
derived from a master seed; 4 chains by default. Checkpoint/resume via
`gd_chain_serialize/restore` (versioned + checksum).

**Node-age parameterization.** Either transformed node proportions within
feasible intervals (Jacobian documented) **or** constrained proposals that
reject invalid ages. No hidden hard bounds beyond configuration.

**`date-pair-context`**: run the tree analysis, then extract the posterior of
the MRCA age of the target pair from the joint node-age posterior.

Cross-encoding model comparison uses pointwise log-likelihood by concept
(`pointwise_loglik.tsv`) so Python can compute WAIC / PSIS-LOO; **never**
compare raw AIC/BIC/WAIC across different data codings without shared
observable space (held-out pairwise cognacy / held-out cell states).

---

## 10. Strict novel-state likelihood — pseudocode

Independent re-implementation of the Rönchen et al. 2024 algorithm (Thm 10 +
Prop 7). All arithmetic in log space via helpers (`logsumexp`, `log1mexp`,
`logdiffexp`); the `p·p − q·q·ret·ret` subtractions use `logdiffexp`.

```
function T1_loglik_concept(m, tree, rates, ε_m, noise):
    # --- Step 1: data extension ---
    f = data[m]                              # partial f_m : V -> labels
    (extendable, det) = compute_determined(f, tree)   # det[v] = class or NIL
    if not extendable:
        if incompatible == error:   return ERR_INCOMPATIBLE
        if incompatible == exclude: return 0.0  (loglik -inf) with warning
        # contamination: L_tree = 0, fall through to noise only
        L_tree = 0.0
    else:
        ftilde = det                          # f̃_m : determined nodes -> class
        # --- Step 2: bottom-up recursion (iterative postorder) ---
        for v in postorder(tree):             # leaves first
            if is_leaf(v):
                p[v] = 0.0                    # log 1
                q[v] = (det[v] != NIL) ? 0.0 : -inf   # log 1 or log 0
            else:
                # assume binary here; nonbinary via Thm A.8/A.9 once appendix
                # is obtained (else resolve polytomy, documented).
                (v1, v2) = children(v)
                ch1, ret1 = branch(v, v1, rates)      # ch=1-exp(-r t), ret=exp(-r t)
                ch2, ret2 = branch(v, v2, rates)
                if det[v] == NIL:                      # Case i
                    pp   = p[v1] + p[v2]               # logs
                    sub  = q[v1] + q[v2] + log(ret1) + log(ret2)
                    p[v] = logdiffexp(pp, sub)
                    t1   = q[v1] + log(ret1) + logdiffexp(p[v2], q[v2]+log(ret2))
                    t2   = q[v2] + log(ret2) + logdiffexp(p[v1], q[v1]+log(ret1))
                    q[v] = logsumexp(t1, t2)
                else:                                  # Case ii: q[v] = p[v]
                    prob1 = (det[v1]!=NIL and det[v1]==det[v]) ? ret1 : ch1
                    prob2 = (det[v2]!=NIL and det[v2]==det[v]) ? ret2 : ch2
                    if det[v1]!=NIL and det[v2]!=NIL:           # iia
                        p[v] = p[v1] + log(prob1) + p[v2] + log(prob2)
                    elif det[v1]==NIL and det[v2]==NIL:         # iib
                        p[v] = logdiffexp(p[v1], q[v1]+log(ret1))
                             + logdiffexp(p[v2], q[v2]+log(ret2))
                    else:                                       # iic
                        if det[v1]!=NIL:   # Det(v1), ¬Det(v2)
                            p[v] = p[v1] + log(prob1)
                                 + logdiffexp(p[v2], q[v2]+log(ret2))
                        else:              # ¬Det(v1), Det(v2)
                            p[v] = p[v2] + log(prob2)
                                 + logdiffexp(p[v1], q[v1]+log(ret1))
                    q[v] = p[v]
        L_tree = exp(p[root])                  # concept tree likelihood in [0,1]
    # --- Step 3: robust contamination mixture ---
    L_noise = noise_loglik(m, data)           # Dirichlet-multinomial / empirical
    return log( (1-ε_m)*L_tree + ε_m*exp(L_noise) )    # stable via logsumexp
```

`compute_determined` finds, for each cognate class `n`, all pairs of observed
nodes labelled `n`, and marks every node on the path between them as
determined to `n`; a node determined to two different classes ⇒ not
extendable. Implemented in O(n log n) via LCA preprocessing (iterative, no
deep recursion).

**Debug mode** emits the extended-state constraints (which nodes are
determined to which classes, and any conflict).

**Simulation** (`simulate` vtable method) runs the generative process of
§3.3 forward over the tree, producing a realization `f : V → ℕ`; used for
SBC and posterior-predictive checks.

---

## 11. Independent brute-force validation — pseudocode

For trees small enough to enumerate (§13.3 of the brief). This is the
**independent** check: it does not reuse the recursion of §10, so a shared
misunderstanding cannot pass.

```
function brute_force_loglik_concept(m, tree, rates):
    f = data[m]
    # 1. enumerate all total label assignments g : V -> {0..K-1} (K = #classes
    #    observed + a few "new" classes; increase K until the likelihood
    #    converges). For tiny trees K is small.
    total = 0.0
    for g in all_total_assignments(V, K):
        if not has_path_property(g, tree):   # Prop 3: P=0 otherwise
            continue
        if not compatible(g, f):             # f_m(u)=f_m(w) <=> g(u)=g(w)
            continue
        # 2. probability of this equivalence class [g] under the generative
        #    process: product over edges of (ret if g(child)==g(parent) else ch)
        #    — this directly uses the generative definition, NOT the p/q recursion
        pr = 0.0  (log)
        for (u,v) in edges(tree):
            t   = age(u) - age(v)
            rr  = exp(-rate(v)*t)
            pr += (g[u]==g[v]) ? log(rr) : log(1-rr)
        total = logsumexp(total, pr)
    return total
```

Brute-force and the §10 recursion must agree within a documented tolerance on:
tiny binary trees, missing leaves, an observed internal node, a polytomy
(until the nonbinary recursion lands, compare against resolved-ladder),
compatible and incompatible patterns, and one vs several cognate classes.
Agreement with the authors' implementation alone is **not** sufficient
(§3.4 of the brief).

---

## 12. Provisional assumptions requiring later simulation

Each must be resolved by an SBC / recovery / posterior-predictive experiment
(M4 / M8 / M11). No published result may rely on only one sensitivity preset.

1. **Contamination noise model** (symmetric Dirichlet–multinomial vs
   empirical categorical): adequacy unverified; choose by posterior-predictive
   fit on injected-noise tests.
2. **Concept-specific contamination with hierarchical shrinkage**: deferred
   from global; identifiability of `ε_c` unverified.
3. **Relaxed-branch-clock identifiability** vs node ages: known confounded;
   needs `σ_b` posterior + calibration leave-one-out; warn when weak.
4. **Split-pulse model**: excluded from v1 unless simulation shows
   identifiability.
5. **Latent primary-synonym Gibbs updates**: mixing on real polymorphic data
   unverified; check ESS / R̂.
6. **Gauss–Hermite quadrature** for lognormal rate integration: only adopted
   if validated against dense Python grids; else grid integration.
7. **Legacy clock presets**: fallback anchors only, never universal; wide SD.
8. **Binary-CTMC ascertainment correction** adequacy under heavy polymorphism
   (King 2026 flags ascertainment × partition interactions).
9. **Cross-encoding comparability**: only via shared observable space
   (held-out pairwise cognacy / held-out cell states); raw cross-coding
   likelihoods are not comparable.
10. **Nonbinary T1 recursion**: pending the figshare appendix supplement
    (§3.3); polytomies are resolved to zero-length ladders. That resolution
    reproduces the generative process of §3.3 exactly and is verified against
    brute force on unresolved polytomies (§3.3, M6), so what is still open is
    (a) whether Theorems A.8/A.9 read a polytomy the same way, and (b) the
    direct `O(k)` recursion, which would avoid the `k − 2` inserted nodes.
    **No published polytomy result before the appendix is checked.**
11. **T1 single-primary-state assumption**: no synonym coexistence modeled;
    the paper's `k`-forms extension is a postponed mitigation
    (`docs/paper-evaluation.md`).
12. **Independence of meanings**: shared by all v1 models; posterior-predictive
    checks (§13.6) test whether real data violate it (King 2026 finds they do
    for the covarion model).

# Inference

> Milestone 0 deliverable. Pairwise posterior integration, fixed-tree MCMC,
> numerical stability, RNG, checkpointing, and outputs. Parameter
> definitions are in `docs/model-spec.md`; this document covers *how* they
> are inferred.

## 1. Pairwise posterior integration (`date-pair`, P0/P1)

Deterministic integration is preferred over MCMC for low-dimensional
pairwise models.

Implemented methods:
- **adaptive grid refinement** in `T` (and in `log r` / latent log-rates when
  not fixed);
- **log-space normalization** (`logsumexp`);
- **numerical CDF construction** from the grid (for quantiles and threshold
  probabilities `P(T > τ)`);
- **adaptive Simpson / Gauss–Kronrod** where suitable;
- **Gauss–Hermite quadrature** for lognormal rate integration **if
  validated** against dense Python grids (else fall back to grid; provisional
  assumption `docs/model-spec.md` §12.6);
- **exact marginalization over synonym candidates** (pairwise; no
  combinatorial enumeration).

The exact P0 closed form and the P1 hierarchy integral are in
`docs/model-spec.md` §8.

**Switch to MCMC** only when pairwise mode includes enough latent variables
that deterministic integration is no longer practical (e.g. many simultaneous
concept-specific `u_{f,c}`). The decision is recorded in the run manifest.

**Outputs**: the evaluated posterior grid (if requested) and the standard
summary set (mean, median, sd, mode where stable, 50/80/90/95% equal-tail,
80/90/95% HPD, threshold probabilities, prior-to-posterior comparison,
prior-dominance diagnostic).

## 2. Fixed-tree MCMC (`date-tree`, `date-pair-context`; T1/T2)

Dependency-free **Metropolis-within-Gibbs**. Target: the posterior of
`docs/model-spec.md` §9.

### 2.1 Parameter blocks

- internal node ages `a(v)`;
- global clock scale;
- family multipliers `β_f`;
- concept log-rates / deviations `μ_c, u_{f,c}`;
- relaxed branch multipliers `m_b`;
- contamination `ε` (and `ε_c` later);
- latent primary synonym states (T1);
- binary CTMC rates `λ_g, λ_l` (T2);
- uncertain calibration hyperparameters when explicitly modelled.

### 2.2 Update types

- **log-scale random walks** for positive rates (`log r`, `σ`, `m_b`);
- **logit-scale** updates for bounded proportions (`ε`, `ε_c`, weights);
- **age proposals** constrained by parent and child ages (reject invalid, or
  use the transformed-proportion parameterization §2.3);
- **block scaling** of a clade's internal ages (propose a correlated shift
  preserving order) — **implemented** (M12): once per sweep the sampler scales
  every free age by a shared `s = exp(N(0, block_step))`, rejecting order/bound
  violations, with Jacobian `s^m`. It is off by default and enabled by the CLI
  for deep trees; `test_mcmc_block_scaling_same_target` confirms it leaves the
  posterior unchanged, and it lifts ESS on the five-language IE example from
  ~13 to ~170;
- **slice sampling** for one-dimensional scale parameters where stable;
- **categorical Gibbs / Metropolis** for latent primary forms — recompute
  **only the affected concept** likelihood (cache invalidation, §4).

Adaptation **only during warm-up**; **freeze** before retained samples. Run
chains with independent deterministic seeds derived from a master seed; 4
chains by default.

### 2.3 The prior on node ages (0.3, roadmap E1)

Before 0.3 there was none: uncalibrated ages were flat on the ordered region
inside a hidden `[0, 1e15]` box. That tilts the root towards the box by
`root^(m-1)`, and with the contamination floor of T1 the unbounded posterior
is improper. On Arawan the uncalibrated root ran to about 9e14 years.

The prior is now generative and proper (`src/gd_tree_prior.h`):

- `log rho ~ N(rate_logmean, rate_logsd²)`, the clock prior;
- `log L ~ N(root_lex_logmean, root_lex_logsd²)`, where `L = rho · a_root` is
  the **lexical root**: expected replacements per concept along a root-to-tip
  path. The defaults are `-1.2` (0.30) and `1.5` (95% from 0.016 to 5.7);
- the other free ages are uniform on the ordered region below the root. This
  is the MrBayes uniform-clock construction, the λ→0 limit of a birth–death
  prior conditioned on the root. Its density is `1/V(a_root)`, where `V` is
  the volume of the region. With dated tips or pinned internal nodes, `V` is a
  polynomial in `a_root` above the oldest fixed age below the root, built
  once by a post-order recursion and checked against Monte Carlo in
  `test_tree_prior.c`.

Calibrations multiply in as factors. Because the cognates identify `L` and
not `rho`, a run with no calendar anchor has a calendar root that is exactly
the lexical posterior divided by the clock prior. `--root-lex-logmean`,
`--root-lex-logsd` and `--prior-only` are the user's handles on this. The
prior-only run reproduces the analytic marginals (`log a_root ~ N(μ_L − μ_ρ,
σ_L² + σ_ρ²)`, uniform ratios), and an SBC with the root free under this prior
passes per parameter (`test_mcmc_tree.c`). The pre-0.3 flat measure remains
as `GD_TREE_PRIOR_NONE` for tests whose truth is drawn that way. It is
refused unless every free age is bounded above.

**Moves.** A non-root age moves by `N(0, step·(hi − lo))` inside the window
that its parent, children and calibrations leave it. The root moves on the
log scale of its height above its oldest child, with the Hastings factor
`(a' − lo)/(a − lo)`. Both steps are dimensionless and adapted in warm-up.
Two scale moves act on the free internal ages; free tips (a list dated to a
narrow window) are left to their own updates:
- the block move scales them by `e^δ`;
- the **clock-preserving move** also shifts `log rho` by `−δ`, which leaves
  every lexical branch length unchanged and moves the chain along the
  rate–age ridge.

On Arawan these moves took the bulk ESS of the root from 9 to over 3,000.

### 2.3.1 Node-age parameterization (original specification)

Maintain `age(parent) > age(child)` on every positive-duration edge. Either:
- **transformed node proportions** within feasible intervals (the feasible
  bounds from `docs/calibrations.md` §2 define each node's `[lo, hi]`);
  document the Jacobian whenever a transform is used; **or**
- **constrained proposals** that reject invalid ages.

No hidden hard bounds beyond those in the configuration. The Jacobian of
whichever parameterization is used is documented in `docs/model-spec.md` and
applied so the posterior is correct on the natural (age) scale.

### 2.4 `date-pair-context`

Run the full tree MCMC, then extract the posterior of the **MRCA age** of the
target pair from the joint node-age posterior. MRCA queries use iterative
LCA (no deep recursion).

### 2.5 As implemented (M8)

`gd_mcmc.{c,h}` is the fixed-tree sampler; `gd_diagnostics.{c,h}` the
convergence diagnostics. The v1 sampler blocks are the **internal node ages**
and a **global strict-clock log-rate**; the family / relaxed-branch /
contamination / latent-synonym blocks of §2.1 extend the same driver in M9.
The likelihood is reached through a model callback
(`gd_mcmc_model.loglik(ctx, tree, rate, pointwise, out, err)`), so one sampler
serves T1, T2, or any future tree model; `gd_mcmc_t1_create` is the ready-made
strict-clock T1 adapter (every concept shares the sampled rate, no allocation
in the callback).

- **Age moves** are single-site random walks on the natural age scale; a
  proposal outside `(max child age, min(parent age, node_hi)) ∩ [node_lo,
  node_hi]` is rejected (the §2.3 constrained-proposal option, so **no
  Jacobian** is needed). A node is *free* iff `node_hi > node_lo`; fixed tips
  and fixed/anchored nodes (`node_lo == node_hi`) are held. The per-node
  calibration density enters the acceptance ratio; only the moved node's prior
  term changes.
- **The rate** moves as a log-scale random walk with a `N(rate_logmean,
  rate_logsd²)` prior on `log r`.
- **Adaptation** is Robbins-Monro on each proposal's log-step toward a 0.44
  target acceptance, applied **only while the running iteration index is below
  `n_warmup`** and frozen for every retained draw
  (`test_mcmc_no_post_warmup_adaptation` asserts the steps do not move after
  warm-up).
- **Determinism**: each chain's PCG32 stream is derived from the master seed,
  so a run is bit-reproducible (`test_mcmc_reproducible`), and the live state
  serializes (versioned + FNV-1a checksum) for chunked resume — an interrupted
  run restored mid-warm-up reproduces the uninterrupted run's draws exactly
  (`test_mcmc_checkpoint_roundtrip`).
- **Diagnostics** (Vehtari et al. 2021, ported line for line from ArviZ 0.23
  and golden-tested against it in `test_diagnostics.c`):
  - `gd_rhat_rank`, the maximum of the bulk (rank-normalised split) and tail
    (folded about the median) R̂;
  - `gd_ess_bulk` and `gd_ess_tail` (the ESS of the 5% and 95% quantiles);
  - `gd_mcse_mean`.

  A run **converges** when every sampled parameter has R̂ < 1.01 and both
  ESSs ≥ 100 per chain. `date-tree` writes its outputs either way, then exits
  with status 3 on failure unless `--allow-unconverged 1` is given.
- **Output**: retained draws carry the total log-likelihood and, optionally,
  the per-concept pointwise log-likelihood (`store_pointwise`) for WAIC /
  PSIS-LOO in Python; the pointwise vector sums to the total
  (`test_mcmc_pointwise_and_fixed_rate`).

**Validation.** `tests/simulation/test_mcmc_tree.c` runs SBC on a strict-clock
quartet (root anchored, two free internal ages, global rate): 80 replicates,
the rank of each true value among the posterior draws is uniform (χ² over 6
bins well below the 0.999 critical value), and a 4-chain run recovers the
generating rate with R̂ < 1.1.

### 2.6 Relaxed clock and robust observations (M9)

The sampler is a typed-parameter engine; the M9 blocks are switched on from the
config and are **off by default** (so the strict-clock path of §2.5 is
unchanged):

- **Relaxed branch clock.** One log-multiplier `log m_b` per branch, with the
  mean-one prior `log m_b ~ N(-σ_b²/2, σ_b²)` and a HalfNormal shrinkage
  hyperprior on `σ_b` (docs/model-spec.md §4). The effective per-branch rate a
  concept sees is `exp(logr + μ_c) · m_b`, applied through the relaxed T1
  likelihood (`gd_t1_loglik_concept_relaxed`). `σ_b` enters only the prior, so
  its Metropolis update needs no likelihood recompute; its full conditional
  includes every branch multiplier. That coupling is Neal's funnel: as `σ_b`
  shrinks the multipliers are pinned together and single-site moves stall.
  Each sweep therefore also makes a **non-centred joint move**. It shifts
  `log σ_b` by δ and maps every `log m_b` through its standardised
  deviation `z_b = (log m_b + σ_b²/2)/σ_b`, which stays fixed. The map's
  Jacobian is `n_b·δ`.
- **Contamination.** `ε ~ Beta(a,b)` sampled on the logit scale, feeding the
  robust mixture of docs/model-spec.md §6, so injected loans / miscoded
  concepts are absorbed instead of biasing the dates.
- **Per-concept rates.** Deviations `μ_c ~ N(0, s²)` around the global log-rate.

**Confounding diagnostic (docs/model-spec.md §7.2).** Branch multipliers and
node ages are confounded: only the product `r·m_b·t_b` is constrained per
branch. `gd_mcmc_run` reports `min_branch_duration_corr`, the most negative
posterior correlation between a branch's `log m_b` and its `log t_b` across
draws, and sets `branch_rate_confounded` when it passes a 0.6 magnitude. With
the ages fixed (identifiable) there is no age variance and nothing is flagged;
with the ages free the branch rate and duration trade off and the flag fires.

**Validation.** `tests/simulation/test_mcmc_relaxed.c`: branch-rate recovery
where identifiable (ages fixed → posterior-mean `log m_b` correlates ~0.98 with
the truth); the confounding warning fires when ages are free (correlation
≈ −0.77); `σ_b` shrinks toward zero on strictly-clocked data; prior-only
sampling recovers the hierarchical prior on `(σ_b, log m_b)` with `σ_b` bulk
ESS above 6,000 (175 without the joint move, which then fails the test); and the rate
estimate stays near the truth while `ε` rises as loan/miscoded concepts are
injected (robustness). Family multipliers `β_f` reduce to the global scale on a
single fixed tree (they need cross-family data, M11); the split-pulse model
stays excluded from v1 (docs/model-spec.md §12).

## 3. Numerical stability (required everywhere)

- **log likelihoods** throughout;
- `logsumexp`, `log1p`, `expm1`;
- stable `log(1 − exp(−x))` (`log1mexp`);
- `logdiffexp` for the `p·p − q·q·ret·ret` subtractions in the T1 recursion;
- scaling / normalization in pruning calculations;
- explicit **NaN and infinity** checks;
- no silent underflow to zero;
- reproducible summation order where practical.

## 4. Caches and incremental recomputation

- cache concept-local structures (the extended data `f̃_m`, determined-node
  sets);
- cache branch durations `t_b` (invalidate when ages change);
- store tree traversal orders (iterative postorder/preorder arrays);
- for **concept-rate** and **latent-primary** proposals, recompute only the
  affected concept's likelihood; for **age** proposals, recompute affected
  concepts on the touched branches;
- invalidate caches explicitly; never stale-read;
- avoid matrices sized by the (infinite) theoretical state space — T1 stores
  only `p_v`, `q_v` per node per concept;
- keep hot data contiguous; no allocation in the innermost likelihood loop.

## 5. Random-number generator

No `rand()`. Implement **PCG32** (or equivalent) with:
- independent streams (one per chain);
- published test vectors (M1 unit tests check exact sequences);
- deterministic serialization (for checkpoint/resume);
- uniform integers **without modulo bias** (rejection);
- open-interval and closed-open uniform doubles;
- normal, exponential, gamma, beta, categorical sampling.

All distributions have statistical unit tests (moments + known-vector checks)
in `tests/unit/`.

## 6. Chain checkpointing

```
gd_chain_serialize(chain, buf, &len, err)
gd_chain_restore(ctx, buf, len, err)
gd_chain_step(chain, n_iterations, progress, err)
```

The serialization format contains a **version** and a **checksum**; mismatch
on restore is an error. This is required for browser pause/resume and
responsive interfaces (chunked stepping).

## 7. Outputs

Strict TSV (schema in `docs/input-formats.md` §9):
`posterior_nodes.tsv`, `posterior_rates.tsv`,
`posterior_hyperparameters.tsv`, `posterior_concepts.tsv`, `diagnostics.tsv`,
`pointwise_loglik.tsv`, `run_manifest.tsv`, `samples.tsv` (or `samples.bin`).

### 7.1 Posterior summaries

For each node age, rate, and hyperparameter: mean, median, sd, mode (where
stable), 50/80/90/95% equal-tail intervals, 80/90/95% HPD intervals,
requested threshold probabilities (e.g. `P(age > 3000)`), prior-to-posterior
comparison, prior-dominance diagnostic.

### 7.2 MCMC diagnostics

computed in C: rank-normalized split-R̂ (bulk and folded); bulk ESS; tail
ESS; MCSE of the mean; the convergence verdict above; acceptance rate by proposal;
number of invalid proposals; chain means and quantiles; autocorrelation
summaries; Monte Carlo standard error; warnings for insufficient mixing.
**Do not** declare convergence from acceptance rate alone.

### 7.3 Identifiability diagnostics

`docs/model-spec.md` §7: absolute-scale-anchored flag; posterior
correlations (clock–root, branch-mult–duration); concepts with posterior ≈
prior; nodes determined mainly by one calibration; parameters with mass at
configured bounds; prior-to-posterior information gain.

### 7.4 Pointwise log likelihood

`pointwise_loglik.tsv` gives per-concept log-likelihood per sample, so Python
can compute **WAIC** and **PSIS-LOO** (where diagnostics permit), held-out-
concept / held-out-language / calibration prediction. Cross-encoding
comparison uses predictions in a shared observable space — **never** raw
AIC/BIC/WAIC across different data codings (`docs/model-spec.md` §3.4,
`docs/validation.md`).

### 7.5 Run manifest

program version; git commit when available; compiler; platform; native/wasm
mode; input checksums; prior-package version; random seed; configuration;
model; reference mode; warnings; start and completion metadata.

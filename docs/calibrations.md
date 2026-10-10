# Calibrations

> Milestone 0 deliverable. Calibration distributions, feasibility checking,
> and required sensitivity procedures. Follows Maurits et al. 2020
> (doi:10.1093/jole/lzz009) on best practices and `docs/model-spec.md` §7 on
> identifiability.

All times are numeric years BP relative to `time_origin_year` (default 1950).
The C engine never parses calendar strings.

## 1. Calibration types

Each calibration names a labelled node or tip (`target_id`) and has a `type`
with parameters `p1/p2/p3`, plus `is_hard` and a `source_key` for
provenance. `is_hard` matters only for `uniform_node`: `1` (or `.`) is a hard
window, `0` a soft one with 2.5% of the mass in a half-normal tail beyond
each end. Every other type is hard or soft by construction, and an `is_hard`
value that contradicts it is an error. Exact parameter meanings:

| `type` | `p1` | `p2` | `p3` | Support / shape |
|---|---|---|---|---|
| `fixed_tip` | age (years BP) | . | . | tip age fixed exactly |
| `uniform_node` | min | max | . | node age ∼ Uniform(min, max); soft tails if `is_hard=0` |
| `normal_node` | mean | sd | . | node age ∼ Normal(mean, sd), in years (truncated ≥0) |
| `lognormal_node` | log_mean | log_sd | . | node age ∼ Lognormal(log_mean, log_sd) |
| `offset_lognormal_node` | offset | log_mean | log_sd | age = offset + Lognormal(log_mean, log_sd); offset is a hard minimum |
| `hard_minimum` | min | . | . | age ≥ min (hard) |
| `hard_maximum` | max | . | . | age ≤ max (hard) |
| `soft_minimum` | min | sd | . | age ≥ min with soft Normal(sd) penalty |
| `soft_maximum` | max | sd | . | age ≤ max with soft Normal(sd) penalty |

Example:
```
calibration_id	target_id	type	p1	p2	p3	is_hard	source_key
latin_tip	latin	fixed_tip	1900	.	.	1	HIST_LATIN_01
romance_root	romance_root	offset_lognormal	1400	5.30	0.25	0	HIST_ROMANCE_02
```

**Targets.** `target_id` is a node label, or a set of tips written
`mrca:A,B,...`: the most recent common ancestor of the listed tips. A set
needs no labelled Newick and survives any relabelling of internal nodes.
Unlabelled internal nodes are reported under generated names of the same
form (`mrca:X,Y`, where X and Y are the alphabetically first tips under the
node's first two children), so every name in an output can be pasted back
into a calibration table.

**Tips.** A tip that no calibration names is an extant language, fixed at age
0. A tip with an uncertain age (an ancient or historical language) needs its
own calibration: `fixed_tip` for a known date, or any node-type distribution
for an uncertain one.

**Sampling windows.** Densities with unbounded support are sampled within
6 standard deviations (on the log scale for the lognormals), which leaves
about 1e-9 of their mass outside: below what any run can resolve.

**Widening** (`--cal-widen f`) multiplies the spread of every soft
calibration: the half-width of a uniform window about its centre, the sd of a
normal or soft bound, and the log-sd of a lognormal or offset lognormal.
Fixed tips and hard bounds are left alone.

**Never** treat a historical date as exact merely because the input gave one
number: a single attested date is `fixed_tip` only for a tip; node dates use
distributions (`normal_node`/`lognormal_node`/`offset_lognormal_node`) that
encode uncertainty.

## 2. Feasibility checking (before inference)

Before any MCMC, the engine:

1. verifies `age(parent) > age(child)` can hold on every positive-duration
   edge given the calibrations;
2. detects mutually inconsistent **hard** bounds (e.g. a hard-minimum on a
   child exceeding a hard-maximum on its parent);
3. verifies dated tips fit the tree (a `fixed_tip` age is consistent with the
   tip's ancestors being older);
4. computes **feasible lower and upper bounds** for every node by propagating
   calibrations up and down the tree (LCA / iterative propagation, no deep
   recursion);
5. **fails with an informative message** when no feasible chronology exists.

The feasible bounds double as the support for node-age proposals and for the
transformed-proportion parameterization (`docs/inference.md` §node ages).

## 3. Absolute-time identifiability

Without a `fixed_tip`, a node calibration, or a clock-scale prior, only
lexical time is identifiable (`docs/model-spec.md` §7.1). The run manifest
records whether an absolute scale is anchored; if not, outputs report
lexical-time units and state "absolute years not identifiable — add a dated
tip, node calibration, or clock-scale prior." The clock/age invariance
(double branch lengths, halve the rate ⇒ identical likelihood) is the formal
reason (Maurits et al. 2020).

## 4. Required sensitivity procedures

Every dated-tree analysis must support, and published results must report:

- **prior-only sampling** — run MCMC with the likelihood disabled to show
  what the priors alone imply (prior-dominance diagnostic);
- **calibration leave-one-out** — drop each calibration in turn, refit,
  report the change in node-age posteriors;
- **one-calibration-at-a-time** — fit with only one calibration active, for
  each calibration;
- **hard vs soft comparison** — where a hard bound exists, compare against
  the soft version;
- **posterior conflict reporting** — flag nodes whose posterior is pushed
  against a calibration bound, or where two calibrations pull a node in
  incompatible directions (calibration conflict diagnostics).

These directly address the Maurits et al. 2020 finding that revised
calibrations can shift a root age by hundreds of years and widen HPDs by
thousands.

## 5. Calibration as prior, not data

Calibrations are **priors** on named nodes/tips in the Bayesian posterior
(`docs/model-spec.md` §9). They are not hard constraints except where
`is_hard=1` (or the type is intrinsically hard, like `fixed_tip` /
`hard_minimum` / `hard_maximum` / the offset of `offset_lognormal_node`).
Uncertain calibration hyperparameters (e.g. a soft bound's `sd`) may be
sampled when explicitly modelled.

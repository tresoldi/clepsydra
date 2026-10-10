# The prior on the concept-rate spread

> 2026-10-11. Data: `concept-rate-prior.tsv`, from the 24 subtrees fitted
> in `clock-engine.md`.

`date-tree` integrates out each concept's rate multiplier, which follows a
log-normal with mean 1 and spread σ_c (roadmap E2.1). In a family of 6–9
languages σ_c is weakly identified. Its prior therefore matters, and it
shifts every lexical age in one direction.

**The problem with a vague prior.** The first default was HalfNormal(1),
which has its median at 0.67. In the misspecification simulator (`simworld`,
clean world, σ_c = 0.8, 6–9 tips), 90% intervals covered only 78% of node
lexical ages. The errors were shared within each replicate: the worst
estimated σ_c at about 0.6 and put every node about 30% too young. Too
little heterogeneity hides deep change. At the sampler level, with truth
drawn from the prior itself, coverage is nominal (the `mixrates` SBC
scenario), so the model is right and the prior is wrong for real families.

**The empirical prior.** The clock-engine runs fitted σ_c on calibrated
subtrees in six families. In the 17 where it is identified (5 or more tips,
posterior sd of log σ_c below 0.25), log σ_c has mean −0.11 and SD 0.30:
σ_c ≈ 0.9, within a factor of about 1.35. Two larger fits agree, Romance
(about 1.05) and Arawan (about 0.8). The default is now
log σ_c ~ N(−0.1, 0.35), set with `--concept-rate-prior-logmean` and
`--concept-rate-prior-logsd`. Under `--concept-rate-gamma 1` the same
numbers apply to the coefficient of variation, which is close to the
log-sd at these values.

**Caveat.** The prior comes from families that the clock experiment also
validates on. A definitive M2 must refit it inside each held-out fold.

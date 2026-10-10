# What the Arawan analysis taught us

> Compiled 2026-10-10 from `~/work/research/phyloarawan/analysis/clepsydra_timedepth/`
> and `~/work/research/arawan-paper/` (FACTS.md, DECISIONS.md, MANUSCRIPT.md, reviews/).
> Clepsydra 0.2.0, commit `18c33c4`. This is the field report that motivates
> `docs/roadmap.md`: every item below is a place where the tool failed a real
> user on a real family.

## The case

- **Data.** KAHD snapshot 2023-06-30: 9 doculects, 7 modern plus Arua
  (Chandless 1869) and a Paumari list by Natterer (1833). 253 concepts, 1,057
  binary characters. Coverage of the two old lists is thin (41 and 70
  concepts).
- **Tree.** The MCC tree of a BEAST run, with Arua's position fixed, although
  the paper calls it unresolved.
- **Calibrations.** None internal. Two dated tips, plus a root window
  Uniform(300, 3000) copied from the BEAST prior, whose only justification is
  "no reason to believe it is older".
- **Outcome.** The paper reports **no date**. BEAST was prior-dominated at
  every node under the pre-registered rule. Clepsydra was brought in as a
  cross-check, and it failed in ways that were partly its own fault.

## Failures, in the order the user hit them

| # | What happened | Root cause in clepsydra | Roadmap item |
|---|---|---|---|
| 1 | No root bound: root about 9.0e14 "years", R-hat up to 1.55, ESS 7, no warning | No tree prior; uncalibrated ages flat on `[0, 1e15]`; posterior improper under the contamination floor | R1.1–R1.4 |
| 2 | With a bound, the root sits on its upper edge (2,795 of 3,000; 7,346 of 8,000) | Flat ordered-age measure tilts ages toward the bound by `root^(m-1)` | R1.1, R1.5 |
| 3 | The rate follows the bound, not the rate prior | Same; the bound is the only real scale information | R1.1, R4 |
| 4 | R-hat 1.07, ESS 67 accepted as a result | No convergence gate; age proposal steps scaled to 1e15 | R0.9, R3.1–R3.3 |
| 5 | Strict clock only, although the paper rejects a strict clock; T2 not available | `--relaxed` not tried; T2 cannot be sampled (one-rate state) | R2.6, R7.4 |
| 6 | Unlabelled internal nodes print as `(node)`; calibrations need hand-labelled Newick | Calibrations target labels, not taxon sets | R0.11, R5.4 |
| 7 | Arua / Paumari 1833 pair dates corrected by hand | `date-pair` hard-codes tip ages to 0 | R0.6 |
| 8 | Rate given per year to `date-tree` and per kiloyear to `date-pair` | Inconsistent units | R0.7 |
| 9 | Paumari 1833 vs modern Paumari: 1,011–1,418 y against a known 165–170 y (6–9x too old) | Sparse, non-specialist old list; no per-language coding-error model | R2.5 |
| 10 | `prior_dominance` 0.01–0.28 looked reassuring while the rate prior set the whole scale | The diagnostic measures the time window, not the scale source | R1.6 |
| 11 | 208 polymorphic cells collapsed by hand to the lowest cognate ID | Synonyms resolved by argmax/row order; no latent-primary sampling | R2.3 |
| 12 | "Data-driven priors" assumed to exist | `concept_priors.tsv` is never loaded by the engine; sources are stubs | R4 |
| 13 | Topology uncertainty (Arua) ignored | Fixed topology only; no tree-set input | R5 |
| 14 | Loan and inheritance columns empty, and the engine ignores them anyway | `gd_forms` reads only four columns | R2.4 |
| 15 | Not citable ("Tresoldi, unpublished") | No archived release | R8.4 |

## The one external check, and what it says

Paumari 1833 against modern Paumari is the only pair with a known separation
(165–170 years). Clepsydra's pairwise model put it at roughly 1,000–1,400
years. The manuscript itself blames the old list's "differences of
transcription and sparse coverage", which is a data explanation, not a
method one. Both are probably true. A dating tool that cannot represent
"this tip is a poor-quality list" will turn coding noise into time, and old
wordlists by non-specialists are exactly what small Amazonian families have.
Hence R2.5: a per-language coding-error rate, which works like a
terminal-branch offset in lexical time.

## What the user asked for, in their words

> "is there any way we can get a better number? ... I basically want to have
> a confidence interval -- also consider published dates for families around
> Arawan."

What can honestly be delivered is set out in `mathematical-analysis.md` §3.
The cognates pin lexical time to about ±8%, and relative ages need no
calibration. The calendar scale is only as good as the cross-family clock
prior. With a realistic between-family rate SD of 0.3–0.4 in log units, a
2,000 BP median carries a 95% interval of roughly 1,000–4,000 BP. That is
the target: an interval of that honesty, with its sources shown, instead of
either a false ±300 years or no number at all.

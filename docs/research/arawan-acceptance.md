# Arawan acceptance log

`docs/roadmap.md` makes Arawan a continuous acceptance case: every gate
re-runs KAHD and records two things. The first is whether the root is still
glued to anything. The second is the Paumari 1833 versus modern Paumari
separation, a negative control that is never used to fit anything (true
separation 165–170 years).

Inputs are those of the original cross-check, unchanged:
`~/work/research/phyloarawan/analysis/clepsydra_timedepth/`. They are 9
doculects, 1,780 forms, the MCC topology with labelled internal nodes, and
two calibration sets: tips only, or tips plus a soft U(300, 3000) root.

## Gate E0 — 2026-10-10, commit 1463e99

Command: `date-tree --tree tree_labeled.nwk --forms forms.tsv --calibrations
<set> --warmup 4000 --samples 4000 --chains 4 --seed 42 --allow-unconverged 1`
(default Swadesh-like rate prior).

| run | root median BP | root 90% HPD | converged | notes |
|---|---|---|---|---|
| tips + U(300, 3000) root (now soft) | 2,887 | upper end 3,117 | no (R-hat 1.017, ESS 279) | still pressed against the window's top |
| tips only | 9.0e14 | to 1.0e15 | **no** (R-hat 1.68, ESS 6) | cap-defined; now flagged instead of silent |

`date-pair paumari_1833 paumari --age-a 167` (rate prior -8.8 per year,
sd 0.8) gives 68 shared concepts and 45 matches. The lineage separation
`2T − 167` has median 2,666 years and a 90% interval of 1,411–4,978, against
a true 165–170. The failure is in the data model (a sparse, non-specialist
1833 list), not the arithmetic. It is the case E2.9 must explain without
being fitted to it.

**Reading.** E0 removed the silent failures: missing cognates, free modern
tips, wrong calibration support, and no convergence check. It was not
meant to fix the dating. The two root runs show the missing tree prior
exactly as `arawan-lessons.md` items 1–3 describe it, and they are E1's
acceptance targets.

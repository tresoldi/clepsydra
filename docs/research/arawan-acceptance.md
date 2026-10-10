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

## Gate E1 (in progress) — 2026-10-10

The proper age prior (lexical root plus uniform ages below it), the new
moves, and the three-ages report. Same command, without
`--allow-unconverged`.

| run | root BP, median (90% HPD) | lexical root (90% HPD) | converged | scale source | rate contraction |
|---|---|---|---|---|---|
| tips only | 3,702 (953–8,389) | 0.539 (0.485–0.591) | yes (R-hat 1.001, ESS 3,355) | tip | 0.43, flagged prior-driven |
| tips + soft U(300, 3000) root | 2,214 (1,372–3,038) | 0.538 (0.486–0.593) | yes (R-hat 1.002, ESS 2,588) | calibration+tip | 0.88 |

**Reading.** Both of E0's failures are gone. The root is finite and
converged without any window. With the window it sits inside the window
rather than on its top edge.

The lexical root is identical in the two runs, about ±10%. Relative ages
agree between them to about 0.01 (e.g. n5 at 0.554 in both). The cognates
determine everything except the calendar scale, and that is now visible in
the output. The calendar difference between the runs is entirely the root
window. Without it, the scale comes from two dated tips 131 and 167 years
old plus the clock prior, and the run says so (`scale_from_prior`).

## Gate E2 (in progress) — 2026-10-11

The run now includes:
- integrated concept-rate heterogeneity;
- sampled ε, with incompatible concepts excluded;
- the subtree move;
- the adequacy panel.

Same command and data. A second run treats Paumari 1833 as an
`ancestor_tip` at 167 BP, i.e. as an early sample of the modern Paumari
lineage.

| run | root BP, median (90% HPD) | lexical root | σ_c | ε | incompatible |
|---|---|---|---|---|---|
| tips only | 4,197 (974–9,869) | 0.65 | 0.77 | 0.004 | 19 of 253 |
| Paumari 1833 as ancestor | 839 (502–1,278) | 0.72 | 0.85 | 0.066 | 19 of 253 |

**Concept-rate heterogeneity lengthens the lexical root by about 20%**
(0.54 at gate E1, now 0.65). This is the correction for the deep-node young
bias predicted in `mathematical-analysis.md` §2. The clock still sets almost
all of the calendar width (`root_var_share_clock` 0.998).

**The adequacy panel finds the 1833 list.** In the tips-only run, Paumari
1833 shares markedly less with every other doculect than the fitted model
predicts. Observed sharing is 0.21–0.26 against about 0.40 replicated, with
`ppp` 0.97–1.00 for each pair. The model also predicts almost no
incompatible concepts, but 7.5% are (ppp < 0.001). Both point to a list
whose coding noise the model reads as divergence.

**Treating it as an ancestor exposes that noise.** It forces all of the
1833 list's divergence from modern Paumari into 167 years. That inflates the
clock and compresses the whole tree, putting the root at 839 BP, and the
adequacy panel rejects the fit: Paumari–Paumari 1833 sharing is observed at
0.66 against 0.84 replicated (ppp 1.00). The negative control is now
diagnosed by the engine, and E2.9 (source-level observation error) is the
item that has to absorb it.

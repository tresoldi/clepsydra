# Research behind the roadmap

Compiled 2026-10-10 to decide what clepsydra needs to become a reference tool
for lexical dating. `docs/roadmap.md` is the plan these documents support.

| file | what it is |
|---|---|
| `arawan-lessons.md` | Field report: every way 0.2.0 failed on Arawan, mapped to a roadmap item |
| `clock-pilot.md` | D2: the cross-family clock pilot (pre-registered; σ_tot 0.77, a no-go for sharp uncalibrated calendar dates) |
| `clock-engine.md` | M1/M2 at pilot scale: the clock experiment with the full engine on Arawan-sized subtrees (σ_tot 0.70; still a no-go) |
| `error-prior.md` | D1.3: list-level observation error from same-language list pairs (per-list e ≈ 0.09) |
| `arawan-acceptance.md` | The Arawan re-run recorded at every roadmap gate |
| `reviews/2026-10-10-codex.md` | Independent adversarial review of the code and model (Codex) |
| `reviews/2026-10-10-claude.md` | Independent adversarial review of the code and model (Claude), with reproduced runs |
| `reviews/2026-10-10-roadmap-codex.md` | Adversarial review of the first roadmap draft (Codex) |
| `reviews/2026-10-10-roadmap-claude.md` | Adversarial review of the first roadmap draft (Claude), with power calculations |
| `literature-lexical-dating.md` | ~110 works on lexicostatistics and phylolinguistic dating; published node ages; validation evidence; ideas; open controversies |
| `literature-dating-methods.md` | ~110 works from molecular clocks, archaeological chronology, stemmatology and computational statistics, with ranked transferable recommendations |
| `mathematical-analysis.md` | The problem as a statistical one: infinite-alleles connection, rate-mixing laws, information limits (what Arawan can achieve), identifiability and scale transfer, alternative inference schemes, coverage |
| `math/` | The Python experiments behind `mathematical-analysis.md` (`python3 -I e2_tree_fisher.py`, etc.) |

The calibration database is in `reference/calibrations/`:
`calibrations.tsv` (367 rows, merged by
`tools/reference_builder/merge_calibrations.py`) from the two regional tables
`oldworld.tsv` and `newworld_pacific.tsv`, whose companion Markdown files give
sources, grades and conflicts. Grades: A strong historical, B
archaeological/inferential, C contested or a published estimate (benchmark
only, never a training calibration), U unverified.

Every DOI in these files was checked against Crossref when it was compiled;
items that could not be verified are marked in place. Numbers marked as read
from figures or secondary sources must be checked before they are quoted.

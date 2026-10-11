# Calibration-sensitivity workflows (Romance)

How much does the Romance chronology depend on its calibrations? Generated with
`clepsydra calib-sensitivity` on this example's tree/forms/calibrations (see
`docs/cli-manual.md`). Reproduce with `./run.sh`.

## `loco.tsv` — leave-one-calibration-out

Each informative calibration is dropped in turn; the tree is re-dated with the
rest, and the dropped calibration's target node is scored against its withheld
value. (Modern age-0 tips are skipped.)

| dropped | held-out BP | result |
|---|---|---|
| `proto_romance_age` | 1725 | **just missed** (inferred 1341, 90% HPD 922–1725): the data and the remaining calibrations put Proto-Romance a little younger than the window's midpoint |
| `latin_age` | 2050 | **not recovered** (inferred ~540): Latin is the only ancient tip, so dropping it removes the tree's absolute-time anchor |
| `root_age` | 2425 | **not recovered** (inferred 3081, 2523–3608): without its window the root drifts older, bounded now by the proper age prior rather than running away |

`max_other_shift_bp` shows how far the rest of the tree moves when each
calibration is removed (up to ~365 BP here, for Latin). A held-out value
falling outside its interval is the point of the exercise: it shows which
calibrations the chronology depends on.

## `widen.tsv` — loosen every calibration

Every soft calibration window is widened ×2 and the tree re-dated once. Node ages
shift modestly (root +145 BP) and interval widths roughly double, quantifying how
much the dates lean on the *tightness* of the calibrations rather than their mere
presence.

Both are `docs/model-spec.md` §4's "calibration-sensitivity diagnostics" as a
one-command workflow: always report how sensitive a published date is to the
calibrations behind it.

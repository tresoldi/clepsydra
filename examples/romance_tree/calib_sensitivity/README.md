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
| `proto_romance_age` | 1725 | **recovered** (inferred 1473, 90% HPD 1115–1884) — the data + remaining calibrations reconstruct it |
| `latin_age` | 2050 | **not recovered** (inferred ~700): Latin is the only ancient tip, so dropping it removes the tree's absolute-time anchor |
| `root_age` | 2425 | **unbounded**: the root has no other upper constraint, so its age runs away — a *load-bearing* calibration the chronology cannot do without |

`max_other_shift_bp` shows how far the rest of the tree moves when each
calibration is removed (up to ~540 BP here for the root). The `unbounded` marker
is the honest report of an identifiability limit, not a numerical failure.

## `widen.tsv` — loosen every calibration

Every soft calibration window is widened ×2 and the tree re-dated once. Node ages
shift modestly (root +299 BP) and interval widths roughly double, quantifying how
much the dates lean on the *tightness* of the calibrations rather than their mere
presence.

Both are `docs/model-spec.md` §4's "calibration-sensitivity diagnostics" as a
one-command workflow: always report how sensitive a published date is to the
calibrations behind it.

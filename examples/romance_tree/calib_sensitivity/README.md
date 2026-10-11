# Calibration-sensitivity workflows (Romance)

How much does the Romance chronology depend on its calibrations? Generated with
`clepsydra calib-sensitivity` on this example's tree/forms/calibrations (see
`docs/cli-manual.md`). Reproduce with `./run.sh`.

## `loco.tsv` — leave-one-calibration-out

Each informative calibration is dropped in turn; the tree is re-dated with the
rest, and the dropped calibration's target node is scored against its withheld
value. (Modern age-0 tips are skipped.)

Default model: relaxed clock, 4 chains × (4,000 + 4,000), seed 42.

| dropped | held-out BP | result |
|---|---|---|
| `proto_romance_age` | 1725 | **recovered** (inferred 1,576, 90% HPD 937–2,221) |
| `latin_age` | 2050 | **not recovered** (inferred ~650): Latin is the only ancient tip, so dropping it removes the tree's absolute-time anchor |
| `root_age` | 2425 | **recovered** (inferred 2,842, 2,193–3,559): without its window the root drifts older, bounded by the proper age prior |

`max_other_shift_bp` shows how far the rest of the tree moves when each
calibration is removed (up to ~290 BP here, for Latin). A held-out value
falling outside its interval is the point of the exercise: it shows which
calibrations the chronology depends on. Under the strict clock (before the
relaxed default) Proto-Romance and the root were both just missed; the relaxed
clock's wider intervals take them in.

## `widen.tsv` — loosen every calibration

Every soft calibration window is widened ×2 and the tree re-dated once. The
root moves 172 BP older and its interval widens from 540 to 951 years; the
shallow nodes barely move (Ibero-Romance +9 BP). The deep chronology leans on
the *tightness* of the root and Proto-Romance windows; the shallow one does
not.

Both are `docs/model-spec.md` §4's "calibration-sensitivity diagnostics" as a
one-command workflow: always report how sensitive a published date is to the
calibrations behind it.

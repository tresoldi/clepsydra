#!/bin/sh
# Calibration-sensitivity workflows on the Romance example (docs/cli-manual.md).
set -e
cd "$(dirname "$0")"
CLE=../../../build/clepsydra
D=../data
C=../calibrations.tsv
T=../tree.nwk

"$CLE" calib-sensitivity --tree "$T" --forms "$D/forms.tsv" --calibrations "$C" \
  --out . --mode loco --seed 42
mv calibration_sensitivity.tsv loco.tsv

"$CLE" calib-sensitivity --tree "$T" --forms "$D/forms.tsv" --calibrations "$C" \
  --out . --mode widen --factor 2.0 --seed 42
mv calibration_sensitivity.tsv widen.tsv

echo; echo "== leave-one-calibration-out =="; cat loco.tsv
echo; echo "== widen x2 =="; cat widen.tsv

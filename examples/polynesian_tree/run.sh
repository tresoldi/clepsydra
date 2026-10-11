#!/bin/sh
# Polynesian chronology from Walworth (2018) cognate data (T1, fixed topology).
# A shallow, well-dated family: the Proto-Polynesian root is anchored by a node
# calibration (archaeology), and the internal ages come out cleanly.
set -e
cd "$(dirname "$0")"
CLE=../../build/clepsydra
mkdir -p output
"$CLE" date-tree \
  --tree tree.nwk --forms data/forms.tsv --calibrations calibrations.tsv \
  --out output --seed 42 \
  --rate-logmean -8.8 --rate-logsd 1.2
echo; echo "== posterior node ages (years BP) =="; cat output/posterior_nodes.tsv

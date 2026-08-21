#!/bin/sh
# Romance chronology from IE-CoR cognate data (T1, fixed topology).
set -e
cd "$(dirname "$0")"
CLE=../../build/clepsydra
mkdir -p output
"$CLE" date-tree \
  --tree tree.nwk --forms data/forms.tsv --calibrations calibrations.tsv \
  --out output --warmup 4000 --samples 4000 --chains 4 --seed 42
echo; echo "== posterior node ages (years BP) =="; cat output/posterior_nodes.tsv

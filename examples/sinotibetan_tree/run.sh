#!/bin/sh
# Sino-Tibetan chronology from Sagart et al. (2019) cognate data (T1, fixed topology).
# A different family from the Indo-European examples, with Old Chinese, Old Tibetan
# and Old Burmese as attested tip anchors.
set -e
cd "$(dirname "$0")"
CLE=../../build/clepsydra
mkdir -p output
"$CLE" date-tree \
  --tree tree.nwk --forms data/forms.tsv --calibrations calibrations.tsv \
  --out output --warmup 5000 --samples 5000 --chains 4 --seed 42 \
  --rate-logmean -8.8 --rate-logsd 1.0
echo; echo "== posterior node ages (years BP) =="; cat output/posterior_nodes.tsv
echo; echo "The Sinitic/Tibetic/Burmish subgroup ages are well estimated; the deep"
echo "root is weakly identified under a strict clock (see README)."

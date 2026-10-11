#!/bin/sh
# Reproduction-style Indo-European tip-dating from IE-CoR (T1, fixed topology).
# Five attested ancient languages with fixed tip ages anchor the chronology.
set -e
cd "$(dirname "$0")"
CLE=../../build/clepsydra
mkdir -p output
"$CLE" date-tree \
  --tree tree.nwk --forms data/forms.tsv --calibrations calibrations.tsv \
  --out output --seed 42 \
  --rate-logmean -8.8 --rate-logsd 0.8
echo; echo "== posterior node ages (years BP) =="; cat output/posterior_nodes.tsv
echo; echo "NOTE: the IE root age is deep-time and rate-prior-sensitive; the wide"
echo "HPD spans the Steppe (~5500-6000 BP) and Anatolian (~8000-9500 BP) ranges."

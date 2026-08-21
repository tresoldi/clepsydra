#!/bin/sh
# Minimal pairwise divergence date (P1) for two Romance sisters from IE-CoR.
set -e
cd "$(dirname "$0")"
CLE=../../build/clepsydra
mkdir -p output
"$CLE" date-pair --forms data/forms.tsv --lang-a Spanish --lang-b Portuguese \
  --out output --t-lo 100 --t-hi 4000
echo; echo "== pairwise posterior =="; cat output/posterior_pair.tsv

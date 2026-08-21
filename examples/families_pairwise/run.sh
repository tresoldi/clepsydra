#!/bin/sh
# Cross-family pairwise divergence gallery (P1). Runs date-pair for a
# representative language pair in each family and tabulates the result.
# Pairwise dating needs no tree or calibration, so it works for any family.
set -e
cd "$(dirname "$0")"
CLE=../../build/clepsydra
mkdir -p output
printf '%-12s %-16s %-16s %8s %10s %14s\n' family lang_a lang_b shared retention "median_BP"
# skip the header line of pairs.tsv
tail -n +2 pairs.tsv | while IFS="$(printf '\t')" read -r fam aid bid an bn shared ret src; do
  "$CLE" date-pair --forms "$fam/forms.tsv" --lang-a "$aid" --lang-b "$bid" \
     --out "output/$fam" --t-lo 100 --t-hi 12000 >/dev/null 2>&1
  med=$(awk -F'\t' '$1=="divergence_median_bp"{print $2}' "output/$fam/posterior_pair.tsv")
  rf=$(awk -F'\t' '$1=="retention_fraction"{print $2}' "output/$fam/posterior_pair.tsv")
  ns=$(awk -F'\t' '$1=="n_shared_concepts"{print $2}' "output/$fam/posterior_pair.tsv")
  printf '%-12s %-16s %-16s %8s %10s %14s\n' "$fam" "$aid" "$bid" "$ns" "$rf" "$med"
done
echo
echo "Divergences are in lexical time scaled by the rate prior (per kiloyear);"
echo "absolute_anchored=0 in each output/<family>/posterior_pair.tsv (no tip anchor)."

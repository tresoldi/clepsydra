#!/bin/sh
# M3: the Arawan result, every layer (docs/research/arawan-result.md).
# Inputs: the phyloarawan clepsydra_timedepth directory (KAHD 2023-06-30) and
# the BEAST B1 tree sample. Each run is capped at 1 GB and runs alone.
#   sh tools/arawan_m3/run.sh OUTDIR
set -eu
OUT=${1:?usage: run.sh OUTDIR}
D=$HOME/work/research/phyloarawan/analysis/clepsydra_timedepth
B1=$HOME/work/research/phyloarawan/prod/arawan.B1.1759696249.trees
CLE=$(cd "$(dirname "$0")/../.." && pwd)/build/clepsydra
mkdir -p "$OUT"
cp "$CLE" "$OUT/clepsydra"            # pin the binary for the whole set
C="$OUT/clepsydra"
run() { systemd-run --user --scope -q -p MemoryMax=1G -p MemorySwapMax=0 "$C" "$@"; }

# Calibration sets: tip dates, plus the 1833 list's error prior (D1.3).
grep -v '^$' "$D/calibrations.tsv" > "$OUT/cal_tips.tsv"
printf 'pau1833_err\tpaumari_1833\ttip_error\t0.088\t13\t.\t.\tD1.3_EMPIRICAL\n' >> "$OUT/cal_tips.tsv"
cp "$OUT/cal_tips.tsv" "$OUT/cal_window.tsv"
printf 'root_window\troot\tuniform_node\t300\t3000\t.\t0\tBEAST_PRIOR_SCENARIO\n' >> "$OUT/cal_window.tsv"
printf 'calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n' > "$OUT/cal_modern.tsv"
echo '((((banawa,(jamamadi,jarawara)n1)n2,(deni,kulina)n4)n5,suruwaha)n6,paumari)root;' > "$OUT/tree_modern.nwk"

COMMON="--forms $D/forms.tsv --warmup 4000 --samples 4000 --chains 4 --seed 42"
run date-tree --tree "$D/tree_labeled.nwk" --calibrations "$OUT/cal_tips.tsv" --out "$OUT/strict" $COMMON || true
run date-tree --tree "$D/tree_labeled.nwk" --calibrations "$OUT/cal_tips.tsv" --out "$OUT/relaxed" $COMMON --relaxed 1 || true
run date-tree --tree "$OUT/tree_modern.nwk" --calibrations "$OUT/cal_modern.tsv" --out "$OUT/modern" $COMMON || true
run date-tree --tree "$D/tree_labeled.nwk" --calibrations "$OUT/cal_tips.tsv" --out "$OUT/clock_empirical" $COMMON \
    --rate-logmean -8.32 --rate-logsd 0.77 || true
run date-tree --tree "$D/tree_labeled.nwk" --calibrations "$OUT/cal_window.tsv" --out "$OUT/window" $COMMON || true

# Topology: the BEAST B1 sample, lower-cased to match the forms.
tr 'A-Z' 'a-z' < "$B1" | sed 's/^#nexus/#NEXUS/' > "$OUT/b1.trees"
grep -v '^root_window' "$OUT/cal_tips.tsv" > "$OUT/cal_trees.tsv"
printf 'clade_id\tmembers\nmadi\tbanawa,jamamadi,jarawara\nmadi_arawa\tarawa,banawa,jamamadi,jarawara\ndeni_kulina\tdeni,kulina\npaumari\tpaumari,paumari_1833\nnon_paumari\tarawa,banawa,jamamadi,jarawara,deni,kulina,suruwaha\n' > "$OUT/clades.tsv"
run date-trees --trees "$OUT/b1.trees" --forms "$D/forms.tsv" --calibrations "$OUT/cal_trees.tsv" \
    --clades "$OUT/clades.tsv" --out "$OUT/trees" --burnin 0.1 --ntrees 20 \
    --warmup 3000 --samples 3000 --chains 4 --seed 42 --allow-unconverged 1 || true
echo done

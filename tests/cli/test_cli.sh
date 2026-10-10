#!/bin/sh
# test_cli.sh — the command line refuses malformed invocations (roadmap E0.14)
# and reports convergence through its exit status (E0.9).
#
# Usage: tests/cli/test_cli.sh path/to/clepsydra
set -u
CLE=${1:?usage: test_cli.sh path/to/clepsydra}
EX=$(cd "$(dirname "$0")/../../examples/romance_tree" && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
fail=0; n=0

# expect STATUS DESCRIPTION -- command...
expect() {
    want=$1; what=$2; shift 3
    n=$((n + 1))
    "$@" >"$TMP/out" 2>"$TMP/err"
    got=$?
    if [ "$got" -ne "$want" ]; then
        echo "  FAIL $what: exit $got, wanted $want"; sed 's/^/    /' "$TMP/err"; fail=$((fail + 1))
    fi
}

TREE="--tree $EX/tree.nwk --forms $EX/data/forms.tsv --calibrations $EX/calibrations.tsv"

expect 1 "unknown option" -- "$CLE" date-tree $TREE --out "$TMP/o" --sampels 10
expect 1 "option without a value (last)" -- "$CLE" date-tree $TREE --out "$TMP/o" --samples
expect 1 "option followed by an option" -- "$CLE" date-tree $TREE --out "$TMP/o" --samples --chains 2
expect 1 "option given twice" -- "$CLE" date-tree $TREE --out "$TMP/o" --seed 1 --seed 2
expect 1 "stray positional argument" -- "$CLE" date-tree $TREE --out "$TMP/o" extra
expect 1 "pair: unknown language" -- "$CLE" date-pair --forms "$EX/data/forms.tsv" --lang-a Spanish --lang-b Klingon --out "$TMP/p"
expect 0 "pair: negative value is a value" -- "$CLE" date-pair --forms "$EX/data/forms.tsv" --lang-a Spanish --lang-b Portuguese --rate-logmean -8.9 --out "$TMP/p"
# A run far too short to converge writes its outputs and exits 3 ...
expect 3 "unconverged run" -- "$CLE" date-tree $TREE --out "$TMP/u" --warmup 5 --samples 12 --chains 2
[ -s "$TMP/u/posterior_nodes.tsv" ] || { echo "  FAIL unconverged run wrote no outputs"; fail=$((fail + 1)); }
grep -q '^converged	0$' "$TMP/u/run_manifest.tsv" || { echo "  FAIL manifest lacks converged 0"; fail=$((fail + 1)); }
# ... unless told it may.
expect 0 "unconverged run, allowed" -- "$CLE" date-tree $TREE --out "$TMP/u" --warmup 5 --samples 12 --chains 2 --allow-unconverged 1

# --write-draws writes every retained draw: header + samples x chains rows.
expect 0 "write draws" -- "$CLE" date-tree $TREE --out "$TMP/d" --warmup 20 --samples 30 --chains 2 --concept-rates 0 --compare-prior 0 --allow-unconverged 1 --write-draws 1
[ "$(wc -l < "$TMP/d/draws.tsv")" -eq 61 ] || { echo "  FAIL draws.tsv row count"; fail=$((fail + 1)); }

# --write-pointwise writes one row per draw and one column per concept.
expect 0 "write pointwise" -- "$CLE" date-tree $TREE --out "$TMP/pw" --warmup 20 --samples 30 --chains 2 --concept-rates 0 --compare-prior 0 --allow-unconverged 1 --write-pointwise 1
[ "$(wc -l < "$TMP/pw/pointwise.tsv")" -eq 61 ] || { echo "  FAIL pointwise.tsv row count"; fail=$((fail + 1)); }

# date-trees over a two-tree sample (the same tree twice) pools clade ages.
cat "$EX/tree.nwk" "$EX/tree.nwk" > "$TMP/two.trees"
printf 'clade_id\tmembers\nibero\tSpanish,Portuguese\n' > "$TMP/clades.tsv"
expect 0 "date-trees on a two-tree sample" -- "$CLE" date-trees --trees "$TMP/two.trees" --forms "$EX/data/forms.tsv" --calibrations "$EX/calibrations.tsv" --clades "$TMP/clades.tsv" --out "$TMP/t" --burnin 0 --ntrees 2 --warmup 200 --samples 200 --chains 2 --concept-rates 0 --allow-unconverged 1
grep -q '^ibero	Spanish,Portuguese	1.000' "$TMP/t/posterior_clades.tsv" || { echo "  FAIL date-trees clade row"; fail=$((fail + 1)); }

# average: two runs of the same model (different seeds) are combined; their
# stacking weights sum to 1. A run without pointwise output is refused.
for sd in 1 2; do
  "$CLE" date-tree $TREE --out "$TMP/av$sd" --warmup 60 --samples 60 --chains 2 --seed $sd --concept-rates 0 --compare-prior 0 --allow-unconverged 1 --write-draws 1 --write-pointwise 1 2>/dev/null
done
expect 0 "average two runs" -- "$CLE" average --runs "$TMP/av1,$TMP/av2" --labels a,b --out "$TMP/avg"
awk -F'\t' 'NR>1 { s += $NF } END { exit !(s > 0.999 && s < 1.001) }' "$TMP/avg/models.tsv" || { echo "  FAIL average weights"; fail=$((fail + 1)); }
grep -q '^root	calendar_bp	' "$TMP/avg/nodes.tsv" || { echo "  FAIL average nodes"; fail=$((fail + 1)); }
expect 1 "average refuses a run without pointwise output" -- "$CLE" average --runs "$TMP/av1,$TMP/d" --out "$TMP/avg2"

if [ "$fail" -ne 0 ]; then echo "$fail of $n CLI checks failed"; exit 1; fi
echo "$n/$n CLI checks passed"

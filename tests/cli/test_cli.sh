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

if [ "$fail" -ne 0 ]; then echo "$fail of $n CLI checks failed"; exit 1; fi
echo "$n/$n CLI checks passed"

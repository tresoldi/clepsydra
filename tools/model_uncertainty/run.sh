#!/bin/sh
# M4.3: the model-uncertainty table (docs/research/model-uncertainty.md).
# Each example family is fitted under every sensitivity variant that applies
# to it, and the variants are stacked on PSIS-LOO by `clepsydra average`.
#   strict / relaxed clock  x  log-normal / gamma concept rates
#   x  uniform-ratio / Yule age prior (Yule only where every tip is modern)
# Runs are sequential, each capped at 1 GB.
#   sh tools/model_uncertainty/run.sh OUTDIR
set -eu
OUT=${1:?usage: run.sh OUTDIR}
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$OUT"
cp "$ROOT/build/clepsydra" "$OUT/clepsydra"
C="$OUT/clepsydra"
run() { systemd-run --user --scope -q -p MemoryMax=1G -p MemorySwapMax=0 "$C" "$@"; }

family() {   # name example-dir yule(0|1) extra-options...
    name=$1; ex=$ROOT/examples/$2; yule=$3; shift 3
    runs=""; labels=""
    mkdir -p "$OUT/$name"
    for clock in strict relaxed; do
        for mix in lognormal gamma; do
            for prior in uniform yule; do
                [ "$prior" = yule ] && [ "$yule" = 0 ] && continue
                v=${clock}_${mix}_${prior}
                o="$OUT/$name/$v"
                run date-tree --tree "$ex/tree.nwk" --forms "$ex/data/forms.tsv" \
                    --calibrations "$ex/calibrations.tsv" --out "$o" \
                    --warmup 4000 --samples 4000 --chains 4 --seed 42 \
                    --compare-prior 0 --write-draws 1 --write-pointwise 1 \
                    --relaxed $([ $clock = relaxed ] && echo 1 || echo 0) \
                    --concept-rate-gamma $([ $mix = gamma ] && echo 1 || echo 0) \
                    --yule $([ $prior = yule ] && echo 1 || echo 0) \
                    --allow-unconverged 1 "$@" > "$o.log" 2>&1 || true
                runs="$runs${runs:+,}$o"; labels="$labels${labels:+,}$v"
            done
        done
    done
    "$C" average --runs "$runs" --labels "$labels" --out "$OUT/$name/average"
}

family romance romance_tree 0
family polynesian polynesian_tree 1 --rate-logmean -8.8 --rate-logsd 1.2
family five_ie five_ie_reproduction 0 --rate-logmean -8.8 --rate-logsd 0.8
echo done

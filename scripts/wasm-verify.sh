#!/usr/bin/env bash
#
# Compares the WebAssembly engine with the native one on the same seeded run.
#
#   ./scripts/wasm-verify.sh
#
# M10's acceptance gate: identical seeded native-vs-Wasm results within a
# documented floating-point tolerance. The tolerance below is a *relative*
# difference on each reported quantity. The two builds are the same C compiled
# by different toolchains at different optimisation levels, so bit-equality is
# not guaranteed by anything and is not asked for; what is asked for is that the
# page and the command line cannot disagree about a date.
#
# Needs node and a built engine (web/build-wasm.sh) plus the native CLI (make).
set -euo pipefail

TOL="${TOL:-1e-9}"

repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_dir"

if [ ! -f web/clepsydra.js ]; then
    echo "wasm-verify: web/clepsydra.js is missing; run web/build-wasm.sh" >&2
    exit 2
fi
if [ ! -x build/clepsydra ]; then
    echo "wasm-verify: build/clepsydra is missing; run make native" >&2
    exit 2
fi
if ! command -v node >/dev/null 2>&1; then
    echo "wasm-verify: node not found" >&2
    exit 2
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

WARMUP=500 SAMPLES=500 CHAINS=2 SEED=42

echo "== date-tree (Romance, ${CHAINS} x ${WARMUP}+${SAMPLES}, seed ${SEED}) =="
./build/clepsydra date-tree \
    --tree examples/romance_tree/tree.nwk \
    --forms examples/romance_tree/data/forms.tsv \
    --calibrations examples/romance_tree/calibrations.tsv \
    --out "$tmp/native" --warmup "$WARMUP" --samples "$SAMPLES" \
    --chains "$CHAINS" --seed "$SEED" --relaxed 0 >/dev/null 2>&1

node scripts/wasm_verify.js tree examples/romance_tree \
    "$WARMUP" "$SAMPLES" "$CHAINS" "$SEED" > "$tmp/wasm_tree.json"

echo "== date-pair (Spanish-Portuguese, IE-CoR) =="
./build/clepsydra date-pair \
    --forms examples/pair_simple/data/forms.tsv \
    --lang-a Spanish --lang-b Portuguese \
    --out "$tmp/native_pair" --t-lo 100 --t-hi 4000 >/dev/null 2>&1

node scripts/wasm_verify.js pair examples/pair_simple/data/forms.tsv \
    Spanish Portuguese 100 4000 > "$tmp/wasm_pair.json"

TOL="$TOL" python3 - "$tmp" <<'PY'
import json, os, sys

tmp = sys.argv[1]
tol = float(os.environ["TOL"])
worst = 0.0
worst_where = "-"
failures = []

def compare(where, native, wasm):
    global worst, worst_where
    if native == wasm:
        rel = 0.0
    else:
        scale = max(abs(native), abs(wasm), 1e-12)
        rel = abs(native - wasm) / scale
    if rel > worst:
        worst, worst_where = rel, where
    if rel > tol:
        failures.append("%-28s native=%.10g  wasm=%.10g  rel=%.3g" %
                        (where, native, wasm, rel))

# ---- date-tree: every node age summary the CLI writes ----
rows = {}
with open(os.path.join(tmp, "native", "posterior_nodes.tsv")) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows[f[0]] = dict(zip(header[1:], (float(x) for x in f[1:])))

tree = json.load(open(os.path.join(tmp, "wasm_tree.json")))
seen = 0
for node in tree["nodes"]:
    label = node["label"] or "(node)"
    if label not in rows:
        failures.append("node %s missing from the native output" % label)
        continue
    nat = rows[label]
    seen += 1
    # The CLI writes one decimal place, so compare at that resolution.
    for key, jkey in (("mean_bp", "mean"), ("median_bp", "median"),
                      ("sd", "sd"), ("hpd90_lo_bp", "hpd_lo"),
                      ("hpd90_hi_bp", "hpd_hi")):
        compare("tree %s.%s" % (label, jkey), nat[key], round(node[jkey], 1))
    compare("tree %s.rhat" % label, nat["rhat"], round(node["rhat"], 3))
    compare("tree %s.ess" % label, nat["ess"], round(node["ess"], 0))

if seen == 0:
    failures.append("no node ages were compared")
if tree.get("progress_calls", 0) <= 0:
    failures.append("the progress callback never fired")

# ---- date-pair ----
pair_native = {}
with open(os.path.join(tmp, "native_pair", "posterior_pair.tsv")) as fh:
    fh.readline()
    for line in fh:
        k, v = line.rstrip("\n").split("\t")
        pair_native[k] = v

pair = json.load(open(os.path.join(tmp, "wasm_pair.json")))
for key, jkey, places in (("divergence_mean_bp", "mean", 1),
                          ("divergence_median_bp", "median", 1),
                          ("divergence_sd_bp", "sd", 1),
                          ("et90_lo_bp", "et90_lo", 1),
                          ("et90_hi_bp", "et90_hi", 1),
                          ("retention_fraction", "retention", 4),
                          ("prior_dominance", "prior_dominance", 3)):
    compare("pair %s" % jkey, float(pair_native[key]), round(pair[jkey], places))
for key, jkey in (("n_shared_concepts", "n_shared"),
                  ("n_cognate_matches", "n_cognate")):
    compare("pair %s" % jkey, float(pair_native[key]), float(pair[jkey]))
if pair_native["model"] != pair["model"]:
    failures.append("pair model: native=%s wasm=%s" % (pair_native["model"], pair["model"]))

print("  compared %d node-age summaries and the pairwise posterior" % seen)
print("  progress callbacks from the Wasm run: %d" % tree.get("progress_calls", 0))
print("  largest relative difference: %.3g (%s), tolerance %.3g" %
      (worst, worst_where, tol))
if failures:
    print("\nFAIL — native and Wasm disagree:")
    for f in failures:
        print("  " + f)
    sys.exit(1)
print("\nOK — the browser engine and the CLI agree within tolerance.")
PY

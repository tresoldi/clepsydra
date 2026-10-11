"""D3: coverage of node lexical ages when the world breaks the model.

Each replicate draws a family of 6-9 languages (random coalescent-like
topology, all tips modern), a clock rate and per-concept rate multipliers
(log-normal, mean 1), and simulates T1 cognate data along the tree. A
scenario then breaks one assumption:

  clean        nothing (the model is right): coverage should be nominal
  borrowing    after the split, a share of concepts is copied between two
               random tips (horizontal transfer)
  coding       each tip's list replaces 8% of forms in error
  time_dep     the rate falls with depth: r(t) = r0 (1 + t / 1000)^-1
  missing_fast concepts with fast rates are missing more often (non-random)
  gamma_rates  concept multipliers are gamma (CV 1.2) instead of log-normal

The fit is date-tree with the strict clock (the default when D3 ran) and no
calendar information. For every internal node the true lexical age (rate x
age, the mean rate) is compared with the 90% HPD of lexical_* in
posterior_nodes.tsv.

    cd tools && python3 -m simworld.run --reps 20 --out ../docs/research/simworld

--tip-error 1 gives every list the D1.3 list-error prior (simworld.md §2).
Runs are sequential, each under a 1 GB cap.
"""
import argparse
import csv
import math
import os
import random
import subprocess
import sys

SCENARIOS = ["clean", "borrowing", "coding", "time_dep", "missing_fast", "gamma_rates"]
NCON = 150


def random_tree(rng, ntip, root_age):
    """Random binary topology; node ages by uniform order statistics."""
    nodes = [{"name": f"L{i}", "children": [], "age": 0.0} for i in range(ntip)]
    active = list(range(ntip))
    m = ntip - 1
    ages = sorted(rng.uniform(0, root_age) for _ in range(m - 1)) + [root_age]
    k = 0
    while len(active) > 1:
        a, b = rng.sample(active, 2)
        nodes.append({"name": f"n{k}", "children": [a, b], "age": ages[k]})
        k += 1
        active = [x for x in active if x not in (a, b)] + [len(nodes) - 1]
    return nodes, len(nodes) - 1


def newick(nodes, v, parent_age=None):
    nd = nodes[v]
    if nd["children"]:
        inner = ",".join(newick(nodes, c, nd["age"]) for c in nd["children"])
        s = f"({inner}){nd['name']}"
    else:
        s = nd["name"]
    if parent_age is not None:
        s += f":{parent_age - nd['age']:.6f}"
    return s


def simulate(rng, nodes, root, rate, mult, scenario):
    """Cognate class per tip, per concept."""
    tips = [i for i, n in enumerate(nodes) if not n["children"]]
    data = []
    nxt = [1]

    def rate_at(t):
        return rate / (1.0 + t / 1000.0) if scenario == "time_dep" else rate

    for c in range(NCON):
        lab = {root: 0}
        stack = [root]
        while stack:
            v = stack.pop()
            for ch in nodes[v]["children"]:
                t0, t1 = nodes[ch]["age"], nodes[v]["age"]
                if scenario == "time_dep":
                    # integral of r(t) over [t0, t1]
                    h = rate * 1000.0 * (math.log1p(t1 / 1000.0) - math.log1p(t0 / 1000.0))
                else:
                    h = rate * (t1 - t0)
                keep = math.exp(-h * mult[c])
                if rng.random() < keep:
                    lab[ch] = lab[v]
                else:
                    lab[ch] = nxt[0]; nxt[0] += 1
                stack.append(ch)
        row = {t: lab[t] for t in tips}
        if scenario == "borrowing" and rng.random() < 0.10:
            a, b = rng.sample(tips, 2)
            row[b] = row[a]
        if scenario == "coding":
            for t in tips:
                if rng.random() < 0.08:
                    row[t] = nxt[0]; nxt[0] += 1
        if scenario == "missing_fast":
            p = min(0.6, 0.15 * mult[c])
            for t in tips:
                if rng.random() < p:
                    row[t] = None
        data.append(row)
    return tips, data


def run_rep(scen, rep, binary, outdir, tip_error=False):
    rng = random.Random(1000 * SCENARIOS.index(scen) + rep)
    ntip = rng.randint(6, 9)
    root_age = rng.uniform(1000, 4000)
    rate = math.exp(rng.gauss(math.log(2e-4), 0.3))
    if scen == "gamma_rates":
        a = 1 / 1.2 ** 2
        mult = [rng.gammavariate(a, 1 / a) for _ in range(NCON)]
    else:
        s = 0.8
        mult = [math.exp(rng.gauss(-s * s / 2, s)) for _ in range(NCON)]
    nodes, root = random_tree(rng, ntip, root_age)
    tips, data = simulate(rng, nodes, root, rate, mult, scen)
    d = os.path.join(outdir, "runs", f"{scen}_{rep:03d}")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "tree.nwk"), "w") as f:
        f.write(newick(nodes, root) + ";\n")
    with open(os.path.join(d, "forms.tsv"), "w") as f:
        f.write("language_id\tconcept_id\tcognate_id\n")
        for c, row in enumerate(data):
            for t in tips:
                if row[t] is not None:
                    f.write(f"{nodes[t]['name']}\tc{c}\tk{row[t]}\n")
    with open(os.path.join(d, "cal.tsv"), "w") as f:
        f.write("calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n")
        if tip_error:
            # The D1.3 empirical list-error prior on every list.
            for t in tips:
                f.write(f"err_{t}\t{nodes[t]['name']}\ttip_error\t0.088\t13\t.\t.\tD1.3\n")
    cmd = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=1G", "-p", "MemorySwapMax=0",
           binary, "date-tree", "--tree", os.path.join(d, "tree.nwk"),
           "--forms", os.path.join(d, "forms.tsv"), "--calibrations", os.path.join(d, "cal.tsv"),
           "--out", os.path.join(d, "out"), "--warmup", "2000", "--samples", "2000",
           "--chains", "4", "--seed", str(rep + 1), "--compare-prior", "0",
           "--relaxed", "0",
           "--allow-unconverged", "1"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    rows = []
    pn = os.path.join(d, "out", "posterior_nodes.tsv")
    man = os.path.join(d, "out", "run_manifest.tsv")
    if not os.path.exists(pn):
        return rows
    conv = "converged\t1" in open(man).read()
    byname = {n["name"]: n for n in nodes}
    with open(pn) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            n = byname.get(r["node"])
            if n is None or not n["children"]:
                continue
            true_lex = rate * n["age"]
            lo, hi, med = float(r["lexical_hpd90_lo"]), float(r["lexical_hpd90_hi"]), float(r["lexical_median"])
            rows.append({"scenario": scen, "rep": rep, "node": r["node"], "ntip": ntip,
                         "true_lex": true_lex, "median": med, "lo": lo, "hi": hi,
                         "covered": int(lo <= true_lex <= hi),
                         "log_bias": math.log(med / true_lex), "converged": int(conv),
                         "is_root": int(n is nodes[root])})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--out", default="../docs/research/simworld")
    ap.add_argument("--binary", default=os.path.abspath("../build/clepsydra"))
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    ap.add_argument("--tip-error", type=int, default=0,
                    help="1: give every list the D1.3 tip_error prior")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    out = os.path.join(args.out, "nodes.tsv")
    cols = ["scenario", "rep", "node", "ntip", "true_lex", "median", "lo", "hi",
            "covered", "log_bias", "converged", "is_root"]
    new = not os.path.exists(out)
    with open(out, "a") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
        for scen in args.scenarios.split(","):
            for rep in range(args.reps):
                for r in run_rep(scen, rep, args.binary, args.out, bool(args.tip_error)):
                    w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
                f.flush()
                print(f"{scen} {rep}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()

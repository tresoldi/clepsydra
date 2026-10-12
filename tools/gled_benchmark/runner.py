"""Run clepsydra and ASJP dating over the benchmark targets and score them."""
import argparse
import csv
import json
import math
import os
import subprocess
import sys

from clock_pilot.glottolog import Glottolog

from . import asjp as A
from .data import Asjp, Gled
from .targets import build, roots
from .tree import newick, upgma

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLE = os.path.join(REPO, "build", "clepsydra")
PRIORS = {"default": (-8.80, 0.80), "empirical": (-8.32, 0.77)}


def implied_log_rates(out, targets):
    """{family: [log rate implied by each two-sided target]}: log(lexical age /
    known age midpoint), from the default-prior runs on disk."""
    res = {}
    for t in targets:
        if t.use not in ("interval", "point"):
            continue
        p = os.path.join(out, t.id, "default", "posterior_nodes.tsv")
        if not os.path.exists(p):
            continue
        with open(p) as f:
            root = next(r for r in csv.DictReader(f, delimiter="\t") if r["node"] == "root")
        lex = float(root["lexical_median"])
        res.setdefault(t.family, []).append(math.log(lex / math.sqrt(t.lo * t.hi)))
    return res


def lofo_prior(rates, family, floor=0.3):
    """Mean and sd of log rate from the other families (each family counts once)."""
    means = [sum(v) / len(v) for k, v in rates.items() if k != family]
    if len(means) < 2:
        return None
    m = sum(means) / len(means)
    sd = math.sqrt(sum((x - m) ** 2 for x in means) / (len(means) - 1))
    return m, max(sd, floor)
ALPHA = 0.1


def interval_score(lo, hi, y):
    s = hi - lo
    if y < lo:
        s += 2 / ALPHA * (lo - y)
    elif y > hi:
        s += 2 / ALPHA * (y - hi)
    return s


def prepare(job, gled, outdir, cals=None):
    os.makedirs(outdir, exist_ok=True)
    lid = {gc: gled.best(gc) for gc in job["gc"]}
    tree = upgma(job["gc"], lambda a, b: gled.distance(lid[a], lid[b]))
    clades = {}
    with open(os.path.join(outdir, "tree.nwk"), "w") as f:
        f.write(newick(tree, clades=clades) + "\n")
    with open(os.path.join(outdir, "forms.tsv"), "w") as f:
        f.write("language_id\tconcept_id\tcognate_id\n")
        for gc in job["gc"]:
            for c, s in gled.forms[lid[gc]].items():
                f.write(f"{gc}\t{c}\t{s}\n")
    with open(os.path.join(outdir, "cal.tsv"), "w") as f:
        f.write("calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n")
        used = []
        for cid, members, lo, hi in (cals or []):
            labs = [l for l, c in clades.items() if set(members) <= c]
            if not labs:
                continue
            lab = min(labs, key=lambda l: len(clades[l]))
            if len(clades[lab]) != len(members):      # not monophyletic in this tree
                continue
            f.write(f"{cid}\t{lab}\tuniform_node\t{lo}\t{hi}\t.\t0\t{cid}\n")
            used.append(cid)
    with open(os.path.join(outdir, "cals_used.txt"), "w") as f:
        f.write(",".join(used) if cals else "")
    return outdir


def run_clepsydra(outdir, variant, warmup, samples, threads):
    lm, ls = PRIORS[variant] if isinstance(variant, str) else variant[1]
    variant = variant if isinstance(variant, str) else variant[0]
    out = os.path.join(outdir, variant)
    cmd = [CLE, "date-tree", "--tree", f"{outdir}/tree.nwk", "--forms", f"{outdir}/forms.tsv",
           "--calibrations", f"{outdir}/cal.tsv", "--out", out, "--warmup", str(warmup),
           "--samples", str(samples), "--rate-logmean", str(lm), "--rate-logsd", str(ls),
           "--threads", str(threads), "--allow-unconverged", "1"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr[-500:])
    with open(os.path.join(out, "posterior_nodes.tsv")) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    root = next(r for r in rows if r["node"] == "root")
    man = {}
    with open(os.path.join(out, "run_manifest.tsv")) as f:
        for line in f:
            k, _, v = line.rstrip("\n").partition("\t")
            man[k] = v
    return root, man


def score(truth_use, lo, hi, med, hlo, hhi):
    """(covered, error-in-log, interval score) against the known age."""
    if truth_use in ("interval", "point"):
        mid = math.sqrt(lo * hi)
        cov = hlo <= hi and hhi >= lo          # HPD overlaps the known window
        err = math.log(med / mid)
        # score against the nearest edge of the window, zero inside it
        y = min(max(med, lo), hi)
        return cov, err, interval_score(hlo, hhi, y)
    if truth_use == "max":
        return hlo <= hi, (math.log(med / hi) if med > hi else 0.0), None
    return hhi >= lo, (math.log(med / lo) if med < lo else 0.0), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--warmup", type=int, default=4000)
    ap.add_argument("--samples", type=int, default=4000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--variants", default="default,empirical")
    ap.add_argument("--max-langs", type=int, default=24)
    ap.add_argument("--asjp-only", type=int, default=0)
    ap.add_argument("--rates-from", default="")
    ap.add_argument("--mode", default="nodes", choices=["nodes", "roots", "calroots"])
    a = ap.parse_args()
    g = Glottolog()
    gled = Gled()
    asj = Asjp(os.path.expanduser("~/lexibank_clone/gled/cldf/parameters.csv"))
    targets = (build(gled, g, max_langs=a.max_langs) if a.mode == "nodes"
               else roots(gled, g, max_langs=a.max_langs))
    calnodes = {}
    if a.mode == "calroots":
        nodes_all = build(gled, g, max_langs=8)
        keep = []
        for t in targets:
            cs = [n for n in nodes_all if n.use in ("interval", "point") and g.is_under(n.gc[0], t.family)
                  and all(g.is_under(x, t.family) for x in n.gc)]
            if not cs:
                continue
            calnodes[t.id] = cs
            extra = [x for n in cs for x in n.gc]
            t.gc = sorted(set(t.gc) | set(extra))
            for x in extra:
                t.side.setdefault(x, g.child_under(x, t.family) or x)
            keep.append(t)
        targets = keep
    node_targets = targets if a.mode == "nodes" else build(gled, g, max_langs=24)
    only = set(filter(None, a.only.split(",")))
    cols = ["id", "family", "use", "lo", "hi", "n_gled", "n_tips", "variant", "median", "hpd_lo",
            "hpd_hi", "prior_median", "prior_lo", "prior_hi", "rhat", "ess", "converged",
            "covered", "log_err", "iscore", "n_asjp", "n_pairs", "asjp_date", "asjp_covered",
            "asjp_log_err", "center", "prior"]
    os.makedirs(a.out, exist_ok=True)
    resf = os.path.join(a.out, "results.tsv")
    done = set()
    if os.path.exists(resf):
        with open(resf) as f:
            done = {(r["id"], r["variant"]) for r in csv.DictReader(f, delimiter="\t")}
    fh = open(resf, "a")
    if not done:
        fh.write("\t".join(cols) + "\n")
    for t in targets:
        if only and t.id not in only:
            continue
        job = {"gc": t.gc}
        # ASJP: same Glottolog codes where ASJP has them
        words = {gc: asj.words[asj.best(gc)] for gc in t.gc if asj.best(gc)}
        ad = A.node_date(words, {gc: t.side[gc] for gc in words})
        for v in a.variants.split(","):
            if (t.id, v) in done:
                continue
            row = {"id": t.id, "family": t.family, "use": t.use,
                   "lo": t.lo if t.lo is not None else ".", "hi": t.hi if t.hi is not None else ".",
                   "n_gled": len(t.gc_all), "n_tips": len(t.gc), "variant": v,
                   "center": f"{getattr(t, 'center', 0):.0f}" if hasattr(t, "center") else "."}
            if not a.asjp_only:
                job = {"gc": t.gc}
                cals = [(n.id, n.gc, n.lo, n.hi) for n in calnodes.get(t.id, [])]
                d = prepare(job, gled, os.path.join(a.out, t.id), cals)
                vv = v
                if v == "lofo":
                    lp = lofo_prior(implied_log_rates(a.rates_from or a.out, node_targets), t.family)
                    if lp is None:
                        continue
                    vv = ("lofo", lp)
                root, man = run_clepsydra(d, vv, a.warmup, a.samples, a.threads)
                if v == "lofo":
                    row["prior"] = f"{lp[0]:.2f}/{lp[1]:.2f}"
                med, hlo, hhi = (float(root[k]) for k in ("median_bp", "hpd90_lo_bp", "hpd90_hi_bp"))
                cov, err, isc = score(t.use, t.lo, t.hi, med, hlo, hhi)
                row.update(median=med, hpd_lo=hlo, hpd_hi=hhi, prior_median=root["prior_median_bp"],
                           prior_lo=root["prior_hpd90_lo_bp"], prior_hi=root["prior_hpd90_hi_bp"],
                           rhat=root["rhat"], ess=root["ess_bulk"],
                           converged=man.get("converged", "."), covered=int(cov),
                           log_err=f"{err:.3f}", iscore="." if isc is None else f"{isc:.0f}")
            if ad:
                lo = t.lo if t.lo is not None else 0.0
                hi = t.hi if t.hi is not None else 1e9
                if t.use in ("interval", "point"):
                    ok = lo <= ad[0] <= hi
                    ae = math.log(ad[0] / math.sqrt(lo * hi)) if ad[0] > 0 else float("nan")
                elif t.use == "max":
                    ok, ae = ad[0] <= hi, (math.log(ad[0] / hi) if ad[0] > hi else 0.0)
                else:
                    ok, ae = ad[0] >= lo, (math.log(ad[0] / lo) if ad[0] < lo else 0.0)
                row.update(n_asjp=len(words), n_pairs=ad[1], asjp_date=f"{ad[0]:.0f}",
                           asjp_covered=int(ok), asjp_log_err=f"{ae:.3f}")
            fh.write("\t".join(str(row.get(c, ".")) for c in cols) + "\n")
            fh.flush()
            print(t.id, v, {k: row.get(k) for k in ("median", "hpd_lo", "hpd_hi", "asjp_date", "converged")},
                  flush=True)


if __name__ == "__main__":
    main()

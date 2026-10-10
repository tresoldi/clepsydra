"""Run the full-engine clock experiment and write docs/research/clock-engine/.

    cd tools && python3 -m clock_engine.run [--stage all|build|date|analyse]
        [--clepsydra BINARY] [--force] [--no-scope]
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time

import numpy as np

from clock_pilot import calibrations as C
from clock_pilot import model
from clock_pilot.glottolog import DEFAULT_ZIP, Glottolog
from clock_pilot.run import aggregate, width_2000, write_tsv

from . import analysis as A
from . import dating as D
from . import subtree as S

OUT = os.path.join(C.REPO, "docs", "research", "clock-engine")
SUB = os.path.join(OUT, "subtrees")
INDEX = os.path.join(SUB, "index.tsv")
INDEX_COLS = ["node", "family", "dataset", "use", "mrca", "n_ingroup", "n_outgroup",
              "target_is_root", "n_concepts", "n_form_rows", "newick"]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def stage_build(g):
    nodes, chosen, fam_primary = S.select_nodes(g, log)
    rows = []
    for n in nodes:
        m = chosen.get(n.id)
        if m is None:
            continue
        st = S.build(n, m, g, os.path.join(SUB, n.id))
        log(f"  {n.id}: {m.dataset}, {len(st.ingroup)}+{len(st.outgroup)} doculects, "
            f"{st.n_concepts} concepts, root={int(st.target_is_root)}")
        rows.append({"node": n.id, "family": n.family, "dataset": m.dataset, "use": n.use,
                     "mrca": n.mrca, "n_ingroup": len(st.ingroup), "n_outgroup": len(st.outgroup),
                     "target_is_root": int(st.target_is_root), "n_concepts": st.n_concepts,
                     "n_form_rows": st.n_rows, "newick": st.newick})
    write_tsv(INDEX, rows, INDEX_COLS)


def read_index():
    with open(INDEX, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stage_date(binary, scope, force):
    idx = read_index()
    # two-sided nodes first: they carry the main fit
    idx.sort(key=lambda r: (r["use"] not in ("interval", "point"),))
    for i, r in enumerate(idx):
        log(f"  [{i + 1}/{len(idx)}] {r['node']}")
        D.date_node(binary, os.path.join(SUB, r["node"]), scope, force, log)


def stage_analyse(g):
    S.pin_calibrations()
    cals, _ = C.load(g)
    nodes, _ = C.dedupe(cals)
    idx = {r["node"]: r for r in read_index()}
    dated = {k: D.read(k, os.path.join(SUB, k)) for k in idx}
    est = A.as_estimates(dated)
    pilot, pilot_res = A.pilot_rows()
    by_id = {n.id: n for n in nodes}

    # ------------------------------------------------------------ nodes.tsv
    nrows = []
    total_wall, n_rerun = 0.0, 0
    for k, r in idx.items():
        n, d, p = by_id[k], dated[k], pilot.get(k, {})
        runs = D.read_runs(os.path.join(SUB, k))
        wall = sum(float(x["wall_s"]) for x in runs)
        total_wall += wall
        n_rerun += d.reruns
        with open(os.path.join(SUB, k, "doculects.tsv"), encoding="utf-8") as fh:
            docs = [f"{x['name']}[{x['role'][0]}]" for x in csv.DictReader(fh, delimiter="\t")]
        row = {"node": k, "family": n.family, "node_name": g.name.get(n.mrca, "."),
               "dataset": r["dataset"], "use": n.use, "n_ingroup": r["n_ingroup"],
               "n_outgroup": r["n_outgroup"], "target_is_root": r["target_is_root"],
               "n_concepts": r["n_concepts"], "doculects": "; ".join(docs),
               "L_median": d.L_median, "L_hpd90_lo": d.L_lo, "L_hpd90_hi": d.L_hi, "s_lnL": d.s_L,
               "rhat": d.target.get("rhat"), "ess_bulk": d.target.get("ess_bulk"),
               "incompatible": d.manifest.get("incompatible_concepts", ".").split(" (")[0],
               "epsilon": d.manifest.get("contamination_epsilon", ".").split("centre ")[-1],
               "converged": int(d.converged), "failing": ",".join(d.failing) or ".", "reruns": d.reruns, "wall_s": round(wall, 1)}
        if n.two_sided:
            c, v = n.log_moments()
            row.update({"T_years": math.exp(c), "sd_lnT": v, "beta_engine": d.ln_L - c,
                        "beta_pilot": p.get("beta_T1"), "beta_pilot_g1": p.get("beta_g1"),
                        "s_lex_pilot": p.get("s_lex_T1")})
        else:
            bd = n.bound()
            row.update({"T_years": bd, "beta_engine": d.ln_L - math.log(bd),
                        "beta_pilot": p.get("beta_T1_bound"),
                        "bound": (">=" if n.use == "max" else "<=")})
        for key in ("beta_pilot", "beta_pilot_g1", "s_lex_pilot"):
            if row.get(key) not in (None, "", "."):
                row[key] = float(row[key])
        nrows.append(row)
    ncols = ["node", "family", "node_name", "dataset", "use", "T_years", "sd_lnT", "bound",
             "L_median", "L_hpd90_lo", "L_hpd90_hi", "s_lnL", "beta_engine", "beta_pilot",
             "beta_pilot_g1", "s_lex_pilot", "converged", "failing", "reruns", "rhat", "ess_bulk",
             "n_ingroup", "n_outgroup", "target_is_root", "n_concepts", "incompatible",
             "epsilon", "wall_s", "doculects"]
    write_tsv(os.path.join(OUT, "nodes.tsv"), nrows, ncols)

    # ------------------------------------------------------------ fits
    results = []
    f, rows, (b, s2, fam, ids, logT) = A.main_fits(nodes, est)
    results += rows
    s_med = float(np.median([est[i].s_lex for i in ids]))
    fac, lo, hi = width_2000(f, s_med)
    results.append({"analysis": "main", "quantity": "width95_factor_at_2000y", "estimate": fac,
                    "lo95": lo, "hi95": hi,
                    "note": f"posterior predictive for a new family's node, median 2000 y, s_L={s_med:.3f}"})
    unconv = [k for k, d in dated.items() if not d.converged]
    if unconv:
        keep = {k for k in dated if dated[k].converged}
        results += A.main_fits(nodes, est, tag="converged_only", keep=keep)[1]
    # family means (as the pilot)
    fb, fs = [], []
    for lab in sorted(set(fam)):
        ix = [i for i, f_ in enumerate(fam) if f_ == lab]
        fb.append(float(np.mean(b[ix])))
        fs.append(float(np.mean(s2[ix]) / len(ix)))
    f1 = model.fit(np.array(fb), np.array(fs), [f"f{i}" for i in range(len(fb))])
    med, lo_, hi_ = f1.post("tot")
    results.append({"analysis": "familymeans", "quantity": "bayes_sd_between_family_means",
                    "estimate": med, "lo95": lo_, "hi95": hi_,
                    "note": "one mean b per family; sigma_f not separable from tau here"})
    results += A.censored_fit(nodes, est, f)

    # paired comparison with the pilot on the same nodes
    pb = np.array([float(pilot[i]["beta_T1"]) for i in ids])
    d_ = b - pb
    results.append({"analysis": "vs_pilot", "quantity": "mean_beta_engine_minus_pilot_T1",
                    "estimate": float(d_.mean()), "lo95": float(d_.min()), "hi95": float(d_.max()),
                    "note": "lo95/hi95 here are min and max over nodes"})
    results.append({"analysis": "vs_pilot", "quantity": "corr_beta_engine_pilot_T1",
                    "estimate": float(np.corrcoef(b, pb)[0, 1]), "lo95": None, "hi95": None,
                    "note": f"{len(b)} two-sided nodes"})
    results.append({"analysis": "vs_pilot", "quantity": "sd_beta_engine",
                    "estimate": float(np.std(b, ddof=1)), "lo95": None, "hi95": None,
                    "note": f"raw SD over nodes; pilot T1: {np.std(pb, ddof=1):.4f}"})

    # ------------------------------------------------------------ LOFO
    lrows = A.run_lofo(nodes, est, g)
    lcols = ["family", "node", "kind", "predictor", "F_train", "n_train", "log_score", "pit",
             "cov50", "cov80", "cov95", "interval_score", "width95_factor", "lo95_years",
             "hi95_years", "median_years", "z", "bound_years"]
    write_tsv(os.path.join(OUT, "lofo.tsv"), lrows, lcols)
    two = [r for r in lrows if r["kind"] in ("interval", "point")]
    for a in aggregate(two):
        for k, v in a.items():
            if k != "predictor":
                results.append({"analysis": f"lofo_{a['predictor']}", "quantity": k,
                                "estimate": v, "lo95": None, "hi95": None, "note": ""})
    results += A.skill(lrows)
    results += A.one_sided_summary(lrows)

    # ------------------------------------------------------------ pilot, for comparison
    keepq = {"bayes_sigma_f", "bayes_tau", "bayes_sigma_tot", "bayes_mu", "reml_sigma_f",
             "reml_tau", "width95_factor_at_2000y"}
    for r in pilot_res:
        if (r["analysis"] in ("main_T1", "main_g1") and r["quantity"] in keepq) or \
                r["analysis"] == "lofo_skill" or \
                (r["analysis"].startswith("lofo_") and r["quantity"] in (
                    "log_score_family_mean", "cov95_family_mean", "cov80_family_mean")):
            results.append({"analysis": "pilot_" + r["analysis"], "quantity": r["quantity"],
                            "estimate": r["estimate"], "lo95": r["lo95"], "hi95": r["hi95"],
                            "note": "copied from docs/research/clock-pilot/results.tsv"})

    results.append({"analysis": "run", "quantity": "n_nodes_dated", "estimate": len(dated),
                    "lo95": None, "hi95": None,
                    "note": f"converged: {sum(d.converged for d in dated.values())}; "
                            f"unconverged: {','.join(unconv) or 'none'}"})
    results.append({"analysis": "run", "quantity": "n_reruns", "estimate": n_rerun,
                    "lo95": None, "hi95": None, "note": "runs rerun with 2x samples"})
    results.append({"analysis": "run", "quantity": "dating_wall_hours", "estimate": total_wall / 3600,
                    "lo95": None, "hi95": None, "note": "sum of clepsydra wall times"})
    results.append({"analysis": "provenance", "quantity": "glottolog", "estimate": g.release,
                    "lo95": None, "hi95": None, "note": f"sha256={g.sha256}"})
    progs = sorted({d.manifest.get("program", ".") + " " + d.manifest.get("git_commit", ".")
                    for d in dated.values()})
    results.append({"analysis": "provenance", "quantity": "clepsydra", "estimate": "; ".join(progs),
                    "lo95": None, "hi95": None, "note": "from run_manifest.tsv"})
    write_tsv(os.path.join(OUT, "results.tsv"), results,
              ["analysis", "quantity", "estimate", "lo95", "hi95", "note"])
    log("wrote", OUT)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="all", choices=["all", "build", "date", "analyse"])
    ap.add_argument("--clepsydra", default=os.path.join(C.REPO, "build", "clepsydra"))
    ap.add_argument("--glottolog", default=DEFAULT_ZIP)
    ap.add_argument("--force", action="store_true", help="redo existing dating runs")
    ap.add_argument("--no-scope", action="store_true", help="do not wrap clepsydra in systemd-run")
    args = ap.parse_args(argv)
    t0 = time.time()
    g = None
    if args.stage in ("all", "build", "analyse"):
        g = Glottolog(args.glottolog)
    if args.stage in ("all", "build"):
        stage_build(g)
    if args.stage in ("all", "date"):
        stage_date(os.path.abspath(args.clepsydra), not args.no_scope, args.force)
    if args.stage in ("all", "analyse"):
        stage_analyse(g)
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

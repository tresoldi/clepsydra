"""Run the D2 σ_f pilot end to end and write docs/research/clock-pilot/*.tsv.

    cd tools && python3 -m clock_pilot.run [--no-lexstat] [--glottolog ZIP]
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

from . import calibrations as C
from . import lexibank as L
from . import lofo
from . import match as M
from . import model
from .glottolog import DEFAULT_ZIP, Glottolog

OUT = os.path.join(C.REPO, "docs", "research", "clock-pilot")
ALPHAS = {"T1": None, "g1": 1.0, "g0.5": 0.5, "g2": 2.0}


def fmt(x, nd=4):
    if x is None:
        return "."
    if isinstance(x, float):
        if math.isnan(x):
            return "nan"
        return f"{x:.{nd}f}"
    return str(x)


def children_covered(m: M.Match) -> int:
    if m.node.spec.kind == "crown":
        return len({m.labels[d] for p in m.pairs for d in p})
    return 2


def node_rank(m: M.Match):
    return (children_covered(m), m.median_shared, len(m.pairs), [-ord(c) for c in m.dataset])


def select(nodes, matches, rule="family"):
    """node id -> chosen Match (or None), following plan.md + deviations 2-3."""
    usable = {n.id: [m for m in matches.get(n.id, []) if m.usable] for n in nodes}
    chosen = {}
    fam_primary = {}
    if rule == "family":
        cnt = defaultdict(lambda: defaultdict(list))
        for n in nodes:
            if n.two_sided:
                for m in usable[n.id]:
                    cnt[n.family][m.dataset].append(m)
        for fam, d in cnt.items():
            best = max(d.items(), key=lambda kv: (len(kv[1]),
                                                   sum(x.median_shared for x in kv[1]) / len(kv[1]),
                                                   [-ord(c) for c in kv[0]]))
            if len(best[1]) >= 2:
                fam_primary[fam] = best[0]
    for n in nodes:
        cands = usable[n.id]
        if not cands:
            chosen[n.id] = None
            continue
        fp = fam_primary.get(n.family)
        pick = [m for m in cands if m.dataset == fp]
        chosen[n.id] = pick[0] if pick else max(cands, key=node_rank)
    return chosen, fam_primary


def estimates(nodes, chosen, alpha, pair_median=False):
    out = {}
    for n in nodes:
        m = chosen.get(n.id)
        if m is None:
            continue
        e = model.lexical(m.K, m.N, m.pair_k, m.pair_n, alpha)
        if pair_median:
            e.log_d = e.log_d_median_pairs
        out[n.id] = e
    return out


def design(nodes, est):
    b, s2, fam, ids, logT = [], [], [], [], []
    for n in nodes:
        if not n.two_sided or n.id not in est:
            continue
        e = est[n.id]
        c, v = n.log_moments()
        b.append(e.log_d - model.LN2 - c)
        s2.append(e.s_lex ** 2 + v ** 2)
        fam.append(n.family)
        ids.append(n.id)
        logT.append(c)
    return np.array(b), np.array(s2), fam, ids, np.array(logT)


def summarise_fit(tag, b, s2, fam, logT, prior=1.0):
    f = model.fit(b, s2, fam, prior_scale=prior)
    sf, tau, mu, varmu = f.reml()
    plo, phi, edge = f.profile_sf()
    tlo, thi = f.profile_tau()
    rows = []
    add = lambda q, est, lo, hi, note="": rows.append(
        {"analysis": tag, "quantity": q, "estimate": est, "lo95": lo, "hi95": hi, "note": note})
    add("n_nodes", len(b), None, None)
    add("F_families", f.F, None, None, ";".join(sorted(set(fam))))
    add("reml_sigma_f", sf, plo, phi, "profile likelihood" + (" (upper end at grid edge 3.0)" if edge else ""))
    add("reml_tau", tau, tlo, thi, "profile likelihood")
    add("reml_mu", mu, mu - 1.96 * math.sqrt(varmu), mu + 1.96 * math.sqrt(varmu), "ln rate per year")
    for q in ("sf", "tau", "tot", "mu"):
        med, lo, hi = f.post(q)
        name = {"sf": "bayes_sigma_f", "tau": "bayes_tau", "tot": "bayes_sigma_tot", "mu": "bayes_mu"}[q]
        add(name, med, lo, hi, f"posterior median, HalfNormal({prior}) on sigma_f and tau")
    if len(b) >= 3:
        sl, se = model.depth_slope(b, logT)
        add("depth_slope", sl, sl - 1.96 * se, sl + 1.96 * se, "OLS of b on ln T (unweighted)")
    return f, rows


def width_2000(f: model.Fit, s_lex: float):
    """×/÷ factor and 95% interval for a node whose median is 2000 years."""
    w, m, s = f.predictive(0.0, s_lex ** 2, tol=1e-10)
    pred = lofo.Mixture(w, m - np.sum(w * m), s)
    xs = np.linspace(-8, 8, 16001)
    cdf = pred.cdf(xs)
    med = np.interp(0.5, cdf, xs)
    lo, hi = np.interp(0.025, cdf, xs) - med, np.interp(0.975, cdf, xs) - med
    return math.exp((hi - lo) / 2), 2000 * math.exp(lo), 2000 * math.exp(hi)


def run_lofo(nodes, est_by_alpha, chosen, g):
    """Return per-node score rows and per-family one-sided rows."""
    rows, cens = [], []
    est1 = est_by_alpha["T1"]
    fams = sorted({n.family for n in nodes if n.id in est1})
    for hold in fams:
        test = [n for n in nodes if n.family == hold and n.id in est1]
        train_nodes = [n for n in nodes if n.family != hold]
        preds_per_alpha = {}
        for tag in ("T1", "g1"):
            est = est_by_alpha[tag]
            b, s2, fam, ids, logT = design(train_nodes, est)
            preds_per_alpha[tag] = (b, s2, fam, logT, ids)
        b1, _, _, logT1, ids1 = preds_per_alpha["T1"]
        train_n = [g.n_languages_under(next(n for n in nodes if n.id == i).mrca) for i in ids1]
        for n in test:
            preds = {}
            fits = {}
            for tag in ("T1", "g1"):
                b, s2, fam, logT, _ = preds_per_alpha[tag]
                e = est_by_alpha[tag][n.id]
                mix, f = lofo.clock_bayes((b, s2, fam), e)
                fits[tag] = f
                preds[f"clock_bayes_{tag}"] = mix
                preds[f"clock_reml_{tag}"] = lofo.clock_reml(f, e)
            preds["old_default_T1"] = lofo.old_default(est1[n.id])
            preds["family_size"] = lofo.family_size(logT1, train_n, g.n_languages_under(n.mrca))
            preds["climatological"] = lofo.climatological(logT1)
            if n.two_sided:
                c, _ = n.log_moments()
                for name, p in preds.items():
                    sc = lofo.score(p, n, c)
                    rows.append({"family": hold, "node": n.id, "predictor": name,
                                 "F_train": fits["T1"].F, "n_train": fits["T1"].n, **sc})
            else:
                for name, p in preds.items():
                    cens.append({"family": hold, "node": n.id, "use": n.use,
                                 "bound_years": n.bound(), "predictor": name,
                                 "censored_log_score": lofo.censored_score(p, n),
                                 "median_years": math.exp(float(np.interp(
                                     0.5, p.cdf(np.linspace(-5, 15, 4001)),
                                     np.linspace(-5, 15, 4001)))) - 50.0})
    return rows, cens


def aggregate(rows, key="predictor"):
    out = []
    preds = sorted({r[key] for r in rows})
    metrics = ["log_score", "interval_score", "cov50", "cov80", "cov95", "pit", "width95_factor"]
    for p in preds:
        rs = [r for r in rows if r[key] == p]
        fam = defaultdict(list)
        for r in rs:
            fam[r["family"]].append(r)
        agg = {"predictor": p, "n_nodes": len(rs), "n_families": len(fam)}
        for mtr in metrics:
            agg[f"{mtr}_node_mean"] = float(np.mean([r[mtr] for r in rs]))
            agg[f"{mtr}_family_mean"] = float(np.mean([np.mean([r[mtr] for r in v]) for v in fam.values()]))
        agg["sum_z2"] = float(sum(r["z"] ** 2 for r in rs))
        out.append(agg)
    return out


def write_tsv(path, rows, cols=None):
    if not rows:
        return
    cols = cols or list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([fmt(r.get(c)) for c in cols])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glottolog", default=DEFAULT_ZIP)
    ap.add_argument("--no-lexstat", action="store_true")
    ap.add_argument("--no-abvd", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    log = lambda *a: print(*a, file=sys.stderr, flush=True)

    g = Glottolog(args.glottolog)
    if g.sha256 != __import__("clock_pilot.glottolog", fromlist=["x"]).EXPECTED_SHA256:
        log("warning: Glottolog zip hash differs from the pinned v5.0 release")
    cals, _ = C.load(g)
    nodes, status = C.dedupe(cals)
    log(f"{len(cals)} rows, {len(nodes)} nodes, {sum(n.two_sided for n in nodes)} two-sided")

    names = L.expert_datasets()
    datasets, matches, prov = {}, defaultdict(list), []
    for name in names + ([] if args.no_abvd else ["abvd"]):
        ds = L.load(name, g.name)
        spc = ds.sets_per_concept
        if spc < 1.05:
            # Deviation 1 (plan.md): a coding with one cognate set per concept
            # measures nothing about replacement.
            log(f"skipping {name}: degenerate cognate coding ({spc:.3f} sets per concept)")
            prov.append({"dataset": name, "commit": ds.commit, "forms": ds.forms_total,
                         "forms_with_cognate": ds.forms_with_cognate,
                         "coverage": round(ds.coverage, 4),
                         "flag": f"DEGENERATE_COGNACY ({spc:.3f} sets/concept): excluded",
                         "role": "excluded", "ancient_dropped": ""})
            datasets[name] = {i: d.name for i, d in ds.doculects.items()}
            del ds
            continue
        anc = [d.name for d in ds.doculects.values() if d.ancient]
        prov.append({"dataset": name, "commit": ds.commit, "forms": ds.forms_total,
                     "forms_with_cognate": ds.forms_with_cognate,
                     "coverage": round(ds.coverage, 4),
                     "flag": "LOW_COVERAGE" if ds.coverage < 0.95 else "",
                     "role": "sensitivity-only" if name == "abvd" else "primary-pool",
                     "ancient_dropped": "; ".join(anc)})
        for n in nodes:
            if name == "abvd" and n.family != "aust1307":
                continue
            m = M.compute(n, ds, g)
            if m is not None:
                matches[n.id].append(m)
        # Keep only the doculect names; the forms of 59 datasets held at once
        # are several GB of Python objects (the cause of an OOM kill).
        datasets[name] = {i: d.name for i, d in ds.doculects.items()}
        del ds
    log("matching done")

    primary_matches = {k: [m for m in v if m.dataset != "abvd"] for k, v in matches.items()}
    chosen, fam_primary = select(nodes, primary_matches, "family")
    chosen_nb, _ = select(nodes, primary_matches, "node")
    chosen_abvd, _ = select(nodes, matches, "family")

    est = {tag: estimates(nodes, chosen, a) for tag, a in ALPHAS.items()}

    # ------------------------------------------------ calibration_matches.tsv
    by_id = {n.id: n for n in nodes}
    node_of = {}
    for n in nodes:
        for c in n.members:
            node_of[c.id] = n
    cm_rows = []
    for c in cals:
        n = node_of.get(c.id)
        base = {"calib_id": c.id, "family_row": c.row["family"], "grade": c.grade,
                "target": c.row["target"], "use": c.use, "node_spec": c.spec.text if c.spec else ".",
                "status": status[c.id], "node_mrca": c.mrca or ".",
                "node_name": g.name.get(c.mrca, ".") if c.mrca else ".",
                "glottolog_family": c.family or "."}
        if n is None or n.primary is not c:
            cm_rows.append({**base, "exclusion_reason": status[c.id] if n is None else "see primary row"})
            continue
        base.update({"node_use": n.use, "T_lo": n.lo, "T_hi": n.hi, "T_center": n.center,
                     "combined": n.combined or "."})
        ms = matches.get(n.id, [])
        if not ms:
            cm_rows.append({**base, "exclusion_reason": "no expert-cognate dataset has doculects on both sides"})
            continue
        for m in sorted(ms, key=lambda m: m.dataset):
            role = ("primary" if chosen.get(n.id) is m else
                    "abvd-sensitivity" if m.dataset == "abvd" else
                    "alternative" if m.usable else "unusable")
            e1 = model.lexical(m.K, m.N, m.pair_k, m.pair_n, None)
            eg = model.lexical(m.K, m.N, m.pair_k, m.pair_n, 1.0)
            ds = datasets[m.dataset]
            nm = lambda ids: ",".join(f"{ds[i]}[{m.labels[i]}]"
                                      if n.spec.kind == "crown" else ds[i] for i in ids)
            row = {**base, "dataset": m.dataset, "role": role,
                   "children_covered": children_covered(m),
                   "n_side_a": len(m.side_a), "n_side_b": len(m.side_b),
                   "side_a": nm(m.side_a), "side_b": nm(m.side_b),
                   "n_pairs": len(m.pairs), "median_shared": m.median_shared,
                   "K": e1.K, "N": e1.N, "p": e1.p, "p_corrected": e1.corrected,
                   "d_T1": e1.d, "s_lex_T1": e1.s_lex, "d_g1": eg.d, "s_lex_g1": eg.s_lex}
            if n.two_sided:
                cc, sd = n.log_moments()
                row.update({"T_geo_mid": math.exp(cc), "sd_lnT": sd,
                            "beta_T1": e1.log_d - model.LN2 - cc,
                            "beta_g1": eg.log_d - model.LN2 - cc,
                            "beta_T1_median_pairs": e1.log_d_median_pairs - model.LN2 - cc})
            else:
                bd = n.bound()
                row.update({"bound_years": bd,
                            "beta_T1_bound": e1.log_d - model.LN2 - math.log(bd),
                            "beta_bound_direction": ">=" if n.use == "max" else "<="})
            if not m.usable:
                row["exclusion_reason"] = f"median shared concepts {m.median_shared:.0f} < {M.MIN_MEDIAN_SHARED}"
            elif role == "alternative":
                row["exclusion_reason"] = "not selected (dataset rule)"
            cm_rows.append(row)
    cols = ["calib_id", "family_row", "grade", "target", "use", "node_spec", "status", "node_mrca",
            "node_name", "glottolog_family", "node_use", "T_lo", "T_hi", "T_center", "combined",
            "dataset", "role", "children_covered", "n_side_a", "n_side_b", "n_pairs", "median_shared",
            "K", "N", "p", "p_corrected", "d_T1", "s_lex_T1", "d_g1", "s_lex_g1", "T_geo_mid", "sd_lnT",
            "beta_T1", "beta_g1", "beta_T1_median_pairs", "bound_years", "beta_T1_bound",
            "beta_bound_direction", "exclusion_reason", "side_a", "side_b"]
    write_tsv(os.path.join(OUT, "calibration_matches.tsv"), cm_rows, cols)

    # ------------------------------------------------------------ main fits
    results = []
    fits = {}
    for tag in ALPHAS:
        b, s2, fam, ids, logT = design(nodes, est[tag])
        f, r = summarise_fit(f"main_{tag}", b, s2, fam, logT)
        fits[tag] = f
        results += r
        if tag in ("T1", "g1"):
            _, r = summarise_fit(f"prior0.5_{tag}", b, s2, fam, logT, prior=0.5)
            results += r
    # pair-median distances
    estm = estimates(nodes, chosen, None, pair_median=True)
    b, s2, fam, ids, logT = design(nodes, estm)
    results += summarise_fit("pairmedian_T1", b, s2, fam, logT)[1]
    # node-best dataset rule
    for tag in ("T1", "g1"):
        e = estimates(nodes, chosen_nb, ALPHAS[tag])
        b, s2, fam, ids, logT = design(nodes, e)
        results += summarise_fit(f"nodebest_{tag}", b, s2, fam, logT)[1]
        e = estimates(nodes, chosen_abvd, ALPHAS[tag])
        b, s2, fam, ids, logT = design(nodes, e)
        results += summarise_fit(f"plusabvd_{tag}", b, s2, fam, logT)[1]
    # family means only (no within-family tau): one b per family
    for tag in ("T1", "g1"):
        b, s2, fam, ids, logT = design(nodes, est[tag])
        fb, fs = [], []
        for lab in sorted(set(fam)):
            ix = [i for i, f_ in enumerate(fam) if f_ == lab]
            fb.append(float(np.mean(b[ix])))
            fs.append(float(np.mean(s2[ix]) / len(ix)))
        f1 = model.fit(np.array(fb), np.array(fs), [f"f{i}" for i in range(len(fb))])
        med, lo, hi = f1.post("tot")
        results.append({"analysis": f"familymeans_{tag}", "quantity": "bayes_sd_between_family_means",
                        "estimate": med, "lo95": lo, "hi95": hi,
                        "note": "one mean b per family; sigma_f not separable from tau here"})
    # censored fit (two-sided + one-sided)
    for tag in ("T1", "g1"):
        e = est[tag]
        rowsf = defaultdict(list)
        for n in nodes:
            if n.id not in e:
                continue
            if n.two_sided:
                c, v = n.log_moments()
                rowsf[n.family].append(("obs", e[n.id].log_d - model.LN2 - c, e[n.id].s_lex ** 2 + v ** 2))
            else:
                val = e[n.id].log_d - model.LN2 - math.log(n.bound())
                rowsf[n.family].append(("ge" if n.use == "max" else "le", val, e[n.id].s_lex ** 2))
        fam_rows = [rowsf[k] for k in sorted(rowsf)]
        sf0, tau0, mu0, _ = fits[tag].reml()
        best = model.censored_fit(fam_rows, (mu0, math.log(max(sf0, 0.05)), math.log(max(tau0, 0.05))))
        lo, hi, edge = model.censored_profile(fam_rows, best, np.linspace(0.0, 3.0, 61))
        n_one = sum(1 for k in rowsf for r in rowsf[k] if r[0] != "obs")
        results += [
            {"analysis": f"censored_{tag}", "quantity": "ml_sigma_f", "estimate": math.exp(best.x[1]),
             "lo95": lo, "hi95": hi, "note": f"ML, two-sided + {n_one} one-sided nodes; F={len(fam_rows)}"
             + ("; upper end at grid edge" if edge else "")},
            {"analysis": f"censored_{tag}", "quantity": "ml_tau", "estimate": math.exp(best.x[2]),
             "lo95": None, "hi95": None, "note": ""},
            {"analysis": f"censored_{tag}", "quantity": "ml_mu", "estimate": best.x[0],
             "lo95": None, "hi95": None, "note": ""}]
    # width for a 2000-year node, full-data fit
    s_med = float(np.median([est["T1"][i].s_lex for i in design(nodes, est["T1"])[3]]))
    for tag in ("T1", "g1"):
        fac, lo, hi = width_2000(fits[tag], s_med)
        results.append({"analysis": f"main_{tag}", "quantity": "width95_factor_at_2000y",
                        "estimate": fac, "lo95": lo, "hi95": hi,
                        "note": f"posterior predictive for a new family's node, median 2000 y, s_lex={s_med:.3f}"})

    # ------------------------------------------------------------- LOFO
    lofo_rows, cens_rows = run_lofo(nodes, est, chosen, g)
    write_tsv(os.path.join(OUT, "lofo.tsv"), lofo_rows)
    write_tsv(os.path.join(OUT, "lofo_one_sided.tsv"), cens_rows)
    agg = aggregate(lofo_rows)
    for a in agg:
        for k, v in a.items():
            if k == "predictor":
                continue
            results.append({"analysis": f"lofo_{a['predictor']}", "quantity": k, "estimate": v,
                            "lo95": None, "hi95": None, "note": ""})
    # family-clustered bootstrap of log-score differences vs baselines
    rng = np.random.default_rng(model.SEED)
    fams = sorted({r["family"] for r in lofo_rows})
    def fam_mean(pred):
        return {f: np.mean([r["log_score"] for r in lofo_rows if r["family"] == f and r["predictor"] == pred]) for f in fams}
    for clock in ("clock_bayes_T1", "clock_bayes_g1"):
        for base in ("climatological", "family_size", "old_default_T1"):
            a, bb = fam_mean(clock), fam_mean(base)
            diff = np.array([a[f] - bb[f] for f in fams])
            boots = [diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(2000)]
            results.append({"analysis": "lofo_skill", "quantity": f"{clock}_minus_{base}",
                            "estimate": float(diff.mean()), "lo95": float(np.percentile(boots, 2.5)),
                            "hi95": float(np.percentile(boots, 97.5)),
                            "note": f"family-weighted mean log-score difference; family bootstrap, F={len(fams)}"})

    # ------------------------------------------------------- uniform coding
    if not args.no_lexstat:
        try:
            from . import uniform
            results += uniform.run(nodes, chosen, datasets, g, log)
        except ImportError as exc:
            results.append({"analysis": "lexstat", "quantity": "skipped", "estimate": None,
                            "lo95": None, "hi95": None, "note": f"lingpy not available: {exc}"})

    for p in prov:
        results.append({"analysis": "provenance", "quantity": p["dataset"], "estimate": p["commit"],
                        "lo95": p["coverage"], "hi95": None,
                        "note": f"{p['role']}; forms={p['forms']}; {p['flag']}; ancient dropped: {p['ancient_dropped']}"})
    results.append({"analysis": "provenance", "quantity": "glottolog", "estimate": g.release,
                    "lo95": None, "hi95": None, "note": f"v5.0 sha256={g.sha256}"})
    for fam, d in sorted(fam_primary.items()):
        results.append({"analysis": "provenance", "quantity": f"family_primary_{fam}", "estimate": d,
                        "lo95": None, "hi95": None, "note": "dataset used for every node it matches"})
    write_tsv(os.path.join(OUT, "results.tsv"), results,
              ["analysis", "quantity", "estimate", "lo95", "hi95", "note"])
    log("wrote", OUT)


if __name__ == "__main__":
    main()

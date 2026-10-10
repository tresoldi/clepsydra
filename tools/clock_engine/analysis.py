"""Meta-analysis, censored fit and LOFO on the engine's lexical ages.

Everything statistical is the pilot's code (`clock_pilot.model`,
`clock_pilot.lofo`, and `design` / `summarise_fit` / `width_2000` /
`aggregate` from `clock_pilot.run`). The only new piece is the per-node
estimate: the pilot's estimators expose `log_d` (ln of the pair distance 2rT)
and `s_lex`; a node age L = rT corresponds to log_d = ln L + ln 2, so the
engine's estimate is passed in that form and every downstream formula
(β = log_d − ln 2 − ln T, ln T_pred = log_d − ln 2 − β) is unchanged.
"""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from clock_pilot import lofo, model
from clock_pilot.run import aggregate, design, summarise_fit, width_2000

PILOT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "docs", "research", "clock-pilot")


@dataclass
class EngineEst:
    log_d: float          # ln L + ln 2 (the pilot's interface)
    s_lex: float          # SD of ln L
    ln_L: float


def as_estimates(dated: dict) -> dict:
    return {k: EngineEst(d.ln_L + model.LN2, d.s_L, d.ln_L) for k, d in dated.items()}


def censored_fit(nodes, est, start_fit):
    rowsf = defaultdict(list)
    for n in nodes:
        if n.id not in est:
            continue
        e = est[n.id]
        if n.two_sided:
            c, v = n.log_moments()
            rowsf[n.family].append(("obs", e.ln_L - c, e.s_lex ** 2 + v ** 2))
        else:
            val = e.ln_L - math.log(n.bound())
            rowsf[n.family].append(("ge" if n.use == "max" else "le", val, e.s_lex ** 2))
    fam_rows = [rowsf[k] for k in sorted(rowsf)]
    sf0, tau0, mu0, _ = start_fit.reml()
    best = model.censored_fit(fam_rows, (mu0, math.log(max(sf0, 0.05)), math.log(max(tau0, 0.05))))
    lo, hi, edge = model.censored_profile(fam_rows, best, np.linspace(0.0, 3.0, 61))
    n_one = sum(1 for k in rowsf for r in rowsf[k] if r[0] != "obs")
    return [
        {"analysis": "censored", "quantity": "ml_sigma_f", "estimate": math.exp(best.x[1]),
         "lo95": lo, "hi95": hi, "note": f"ML, two-sided + {n_one} one-sided nodes; F={len(fam_rows)}"
         + ("; upper end at grid edge" if edge else "")},
        {"analysis": "censored", "quantity": "ml_tau", "estimate": math.exp(best.x[2]),
         "lo95": None, "hi95": None, "note": ""},
        {"analysis": "censored", "quantity": "ml_mu", "estimate": best.x[0],
         "lo95": None, "hi95": None, "note": ""}]


def run_lofo(nodes, est, g):
    """The pilot's LOFO with one estimator (the engine's)."""
    rows = []
    fams = sorted({n.family for n in nodes if n.id in est})
    for hold in fams:
        test = [n for n in nodes if n.family == hold and n.id in est]
        train = [n for n in nodes if n.family != hold]
        b, s2, fam, ids, logT = design(train, est)
        by_id = {n.id: n for n in nodes}
        train_n = [g.n_languages_under(by_id[i].mrca) for i in ids]
        for n in test:
            e = est[n.id]
            mix, f = lofo.clock_bayes((b, s2, fam), e)
            preds = {"clock_bayes": mix, "clock_reml": lofo.clock_reml(f, e),
                     "old_default": lofo.old_default(e),
                     "family_size": lofo.family_size(logT, train_n, g.n_languages_under(n.mrca)),
                     "climatological": lofo.climatological(logT)}
            for name, p in preds.items():
                base = {"family": hold, "node": n.id, "kind": n.use, "predictor": name,
                        "F_train": f.F, "n_train": f.n}
                if n.two_sided:
                    c, _ = n.log_moments()
                    rows.append({**base, **lofo.score(p, n, c)})
                else:
                    xs = np.linspace(-5, 15, 4001)
                    rows.append({**base, "bound_years": n.bound(),
                                 "log_score": lofo.censored_score(p, n),
                                 "median_years": math.exp(float(np.interp(0.5, p.cdf(xs), xs))) - 50.0})
    return rows


def skill(rows, seed=model.SEED):
    two = [r for r in rows if r["kind"] in ("interval", "point")]
    rng = np.random.default_rng(seed)
    fams = sorted({r["family"] for r in two})

    def fam_mean(pred):
        return {f: np.mean([r["log_score"] for r in two if r["family"] == f and r["predictor"] == pred])
                for f in fams}
    out = []
    for base in ("climatological", "family_size", "old_default"):
        a, bb = fam_mean("clock_bayes"), fam_mean(base)
        diff = np.array([a[f] - bb[f] for f in fams])
        boots = [diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(2000)]
        out.append({"analysis": "lofo_skill", "quantity": f"clock_bayes_minus_{base}",
                    "estimate": float(diff.mean()), "lo95": float(np.percentile(boots, 2.5)),
                    "hi95": float(np.percentile(boots, 97.5)),
                    "note": f"family-weighted mean log-score difference; family bootstrap, F={len(fams)}"})
    return out


def one_sided_summary(rows):
    out = []
    one = [r for r in rows if r["kind"] in ("max", "min")]
    for p in sorted({r["predictor"] for r in one}):
        v = [r["log_score"] for r in one if r["predictor"] == p]
        out.append({"analysis": f"lofo_one_sided_{p}", "quantity": "censored_log_score_node_mean",
                    "estimate": float(np.mean(v)), "lo95": None, "hi95": None,
                    "note": f"{len(v)} one-sided held-out nodes"})
    return out


def main_fits(nodes, est, tag="main", keep=None):
    """Bayesian/REML fit on the two-sided nodes (optionally a subset)."""
    sub = [n for n in nodes if keep is None or n.id in keep]
    b, s2, fam, ids, logT = design(sub, est)
    f, rows = summarise_fit(tag, b, s2, fam, logT)
    rows = list(rows)
    if tag == "main":
        _, r = summarise_fit("prior0.5", b, s2, fam, logT, prior=0.5)
        rows += r
    return f, rows, (b, s2, fam, ids, logT)


def pilot_rows():
    """Pilot per-node betas (primary rows) and its headline results."""
    with open(os.path.join(PILOT_DIR, "calibration_matches.tsv"), encoding="utf-8") as fh:
        per = {r["calib_id"]: r for r in csv.DictReader(fh, delimiter="\t") if r["role"] == "primary"}
    with open(os.path.join(PILOT_DIR, "results.tsv"), encoding="utf-8") as fh:
        res = list(csv.DictReader(fh, delimiter="\t"))
    return per, res

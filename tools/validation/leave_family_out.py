"""
Leave-family-out cross-family validation (M11).

When the reference package is built with ``--leave-family-out``, a consumer
dating the held-out family cannot use that family's own clock row; it must
predict the family's clock offset from the *other* families. The scientific
question is whether that out-of-sample prediction is **calibrated**: does the
leave-one-out predictive interval cover the held-out family's true clock offset
at its nominal rate?

This module answers that with a genuine Monte-Carlo coverage experiment on a
hierarchical normal model (the standard model behind a family-clock package):

    true offsets   β_g ~ N(μ, τ²)            (families vary around a grand mean)
    stored estimate y_g ~ N(β_g, ω_g²)       (each family's own fit has error ω_g)

For each held-out family f we build the leave-one-out predictive for β_f from
{(y_g, ω_g) : g ≠ f} — grand mean μ̂, between-family variance τ̂² by method of
moments, plus the uncertainty in μ̂ — and measure how often the true β_f lands
in the interval, across many independent draws and every choice of held-out
family. Coverage is an *empirical* number that can come out below or above
nominal; it is not guaranteed by construction. Three priors are compared:

  * ``pooled_loo``        — the honest hierarchical leave-one-out predictive;
  * ``unknown_fallback``  — the flat UNKNOWN row the current builder emits for a
                            left-out family (fixed sd), a baseline to beat;
  * ``family_specific``   — the held-out family's OWN estimate (in-sample,
                            cheating): a positive control that is calibrated by
                            construction and is exactly what leave-family-out
                            must not use.

That the flat fallback and the point-conditional predictive give *different*
coverage from the pooled one is the whole point: the check can fail, so a pass
is evidence. A deliberately mis-built predictive that drops τ̂² undercovers
badly (see ``run_coverage(broken=True)`` and the tests), which is what gives
this experiment teeth.

Scope. This validates the leave-one-out *procedure* on synthetic hierarchical
data. It is NOT the same as benchmarking inferred dates against real historical
calibrations for held-out families — that needs calibrated datasets run
end-to-end through the C engine and is tracked as future work in
``docs/paper-evaluation.md``. The two are complementary; only the second can
detect model misspecification against the real world.

Stdlib only (numpy/pandas not required), matching the reference builder.

Usage:
  python -m validation.leave_family_out            # full coverage table
  python -m validation.leave_family_out --tau 0.6  # a wider true dispersion
"""
import argparse
import math
import os
import sys

# Make the sibling reference_builder package importable when run directly.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from reference_builder import synthetic_sources as synth

import random


def _z_for(level):
    # inverse standard-normal via a small rational approximation (Acklam),
    # sufficient for the 0.5..0.999 range used here.
    p = 0.5 + level / 2.0
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p <= phigh:
        q = p - 0.5
        r = q*q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)


def _t_quantile(level, dof):
    """Two-sided Student-t multiplier for a central `level` interval, from the
    normal quantile via a Cornish-Fisher expansion (accurate for dof >= ~5;
    stdlib only). Accounts for estimating the predictive variance from a small
    number of families rather than knowing it."""
    z = _z_for(level)
    if dof <= 0:
        return z
    v = float(dof)
    g1 = (z**3 + z) / 4.0
    g2 = (5*z**5 + 16*z**3 + 3*z) / 96.0
    g3 = (3*z**7 + 19*z**5 + 17*z**3 - 15*z) / 384.0
    return z + g1/v + g2/(v*v) + g3/(v*v*v)


def family_omegas():
    """(family_id, ω) fit-error SDs and the UNKNOWN fallback sd, from the
    synthetic family package — the ω_g play the role of each family's own
    posterior clock SD."""
    rows = synth.generate_family_priors()
    omegas = [(r["family_id"], r["log_clock_sd"]) for r in rows
              if r["family_id"] != "UNKNOWN"]
    unknown_sd = next(r["log_clock_sd"] for r in rows
                      if r["family_id"] == "UNKNOWN")
    return omegas, unknown_sd


def pooled_loo_predictive(others, drop_tau2=False):
    """Hierarchical leave-one-out predictive N(mu, sd²) for a new family's true
    offset, from the other families' (estimate, ω) pairs.

    others: list of (y_g, omega_g). Method of moments:
      mu   = mean of y_g
      S²   = sample variance of y_g  (E[S²] = τ² + mean ω²)
      τ̂²  = max(0, S² − mean ω²)
      Var(mu) = (τ̂² + mean ω²) / m
      sd²  = τ̂² + Var(mu)           (spread of a new family's offset + mean error)
    Returns (mu, sd, tau2_hat). Needs at least 2 other families.
    """
    m = len(others)
    if m < 2:
        raise ValueError("need >= 2 other families for a pooled predictive")
    ys = [y for (y, _w) in others]
    w2 = [w * w for (_y, w) in others]
    mu = sum(ys) / m
    s2 = sum((y - mu) ** 2 for y in ys) / (m - 1)
    mean_w2 = sum(w2) / m
    tau2 = max(0.0, s2 - mean_w2)
    var_mu = (tau2 + mean_w2) / m
    sd2 = var_mu if drop_tau2 else (tau2 + var_mu)
    return mu, math.sqrt(sd2), tau2


def run_coverage(tau_true, omegas=None, unknown_sd=0.80,
                 levels=(0.5, 0.8, 0.9, 0.95), n_reps=6000, seed=20240601,
                 broken=False, n_families=None):
    """Monte-Carlo leave-one-out coverage for three priors.

    Each rep draws a whole hierarchical world (β_g ~ N(0, τ_true²),
    y_g ~ N(β_g, ω_g²)); then every family is held out in turn and its true β_f
    is checked against each prior's central interval at each level. Returns
    {arm: {level: coverage_fraction}} for arms pooled_loo, unknown_fallback,
    family_specific. Independent randomness throughout (no shared-seed trick).
    """
    if omegas is None:
        omegas, unknown_sd = family_omegas()
    om = [w for (_f, w) in omegas]
    if n_families is not None and n_families > len(om):
        # tile the fit-error SDs to a larger synthetic family count, to show the
        # leave-one-out predictive is consistent (coverage → nominal as G grows).
        om = [om[i % len(om)] for i in range(n_families)]
    G = len(om)
    rng = random.Random(seed)
    z = {lv: _z_for(lv) for lv in levels}
    # Student-t multiplier for the pooled predictive: the between-family
    # variance is estimated from m = G-1 families, so dof = m-1.
    tq = {lv: _t_quantile(lv, G - 2) for lv in levels}
    arms = ("pooled_loo", "unknown_fallback", "family_specific")
    hits = {a: {lv: 0 for lv in levels} for a in arms}
    trials = 0
    for _ in range(n_reps):
        beta = [rng.gauss(0.0, tau_true) for _ in range(G)]
        y = [rng.gauss(beta[g], om[g]) for g in range(G)]
        for f in range(G):
            others = [(y[g], om[g]) for g in range(G) if g != f]
            mu, sd, _t2 = pooled_loo_predictive(others, drop_tau2=broken)
            trials += 1
            for lv in levels:
                # pooled hierarchical leave-one-out predictive (t-tailed)
                if abs(beta[f] - mu) <= tq[lv] * sd:
                    hits["pooled_loo"][lv] += 1
                # flat UNKNOWN fallback the builder currently emits
                if abs(beta[f] - 0.0) <= z[lv] * unknown_sd:
                    hits["unknown_fallback"][lv] += 1
                # in-sample family-specific control (cheating; calibrated)
                if abs(beta[f] - y[f]) <= z[lv] * om[f]:
                    hits["family_specific"][lv] += 1
    cov = {a: {lv: hits[a][lv] / trials for lv in levels} for a in arms}
    return cov


def evaluate(tau_true=0.4, n_reps=6000, seed=20240601):
    """Report dict: coverage per arm/level plus the pooled-LOO calibration gap."""
    levels = (0.5, 0.8, 0.9, 0.95)
    omegas, unknown_sd = family_omegas()
    cov = run_coverage(tau_true, omegas, unknown_sd, levels, n_reps, seed)
    max_gap = max(abs(cov["pooled_loo"][lv] - lv) for lv in levels)
    return {
        "tau_true": tau_true,
        "n_families": len(omegas),
        "levels": list(levels),
        "coverage": cov,
        "pooled_loo_max_gap": max_gap,
        "unknown_sd": unknown_sd,
    }


def main(argv):
    ap = argparse.ArgumentParser(prog="validation.leave_family_out")
    ap.add_argument("--tau", type=float, default=0.4,
                    help="true between-family clock dispersion (log scale)")
    ap.add_argument("--reps", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=20240601)
    args = ap.parse_args(argv[1:])
    rep = evaluate(args.tau, args.reps, args.seed)
    print("Leave-one-family-out coverage (τ_true=%.2f, %d families, %d reps)"
          % (rep["tau_true"], rep["n_families"], args.reps))
    print("%-18s %s" % ("prior \\ level",
                        "  ".join("%5.2f" % lv for lv in rep["levels"])))
    for arm in ("pooled_loo", "unknown_fallback", "family_specific"):
        row = "  ".join("%5.3f" % rep["coverage"][arm][lv]
                        for lv in rep["levels"])
        print("%-18s %s" % (arm, row))
    print("pooled_loo max |coverage − nominal| = %.3f" % rep["pooled_loo_max_gap"])
    print("(pooled_loo should track nominal; unknown_fallback is typically "
          "conservative when τ_true < its fixed sd; family_specific is the "
          "in-sample control.)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

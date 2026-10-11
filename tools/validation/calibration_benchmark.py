"""
Held-out-calibration benchmark (M11, real end-to-end).

The synthetic coverage study in `leave_family_out.py` certifies that the
uncertainty *arithmetic* is honest. This module is the complementary real-world
check the model cards call for: hold out a real calibration, date the tree
**through the actual C engine** with that calibration withheld, and ask whether
the held-out age falls in the inferred interval — a genuine, non-circular
predictive test on real cognate data that can also catch model misspecification.

Leave-one-calibration-out (LOCO). A heterochronous tree (e.g. Indo-European with
attested ancient tips: Vedic, Ancient Greek, Latin, Old Irish, Old Church
Slavonic) carries several tip-age calibrations. For each one:
  1. rewrite the calibration file with that tip **freed** (its fixed age replaced
     by a broad uniform prior, so the sampler infers it — tip dating);
  2. run `clepsydra date-tree` with the remaining calibrations still fixed;
  3. read the freed node's posterior age and compare it to the withheld value.
The withheld target is the age already in the calibration file — nothing is
fabricated; the test is whether the model recovers a datum it was not given.

Aggregated over all held-out calibrations: empirical coverage (fraction whose
true age lands in the 90% HPD), signed error, RMSE, interval width, and the
interval score of the 90% interval (Gneiting & Raftery 2007), which charges the
width plus 20x the distance by which a miss falls outside. The same run's
prior-only interval for the freed node (the `prior_*` columns: every other
calibration kept, no cognate data) is scored the same way as a baseline, and
the skill is 1 - IS(posterior) / IS(prior). This is
NOT a large-sample coverage guarantee — it is a handful of held-out points from
the trees you supply, so report N and read it as a calibrated-ness spot check,
not a p-value. Deep internal nodes stay rate-prior-dominated; tips near the
leaves are far better constrained (see the per-target ess/rhat columns).

Stdlib only. The C engine does the inference; this is orchestration + scoring.

Usage:
  # real run (slow: one MCMC per held-out calibration)
  python -m validation.calibration_benchmark \\
      --binary build/clepsydra --tree examples/five_ie_reproduction/tree.nwk \\
      --forms examples/five_ie_reproduction/data/forms.tsv \\
      --calibrations examples/five_ie_reproduction/calibrations.tsv \\
      --out /tmp/loco --warmup 6000 --samples 6000 --chains 4 \\
      --rate-logmean -8.8 --rate-logsd 0.8
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

ALPHA = 0.1   # the HPD level is 90%


def interval_score(lo, hi, y, alpha=ALPHA):
    """Interval score of a central (1 - alpha) interval [lo, hi] for outcome y.
    Lower is better: the width, plus 2/alpha times any miss distance."""
    s = hi - lo
    if y < lo:
        s += 2.0 / alpha * (lo - y)
    elif y > hi:
        s += 2.0 / alpha * (y - hi)
    return s


def parse_calibrations(path):
    """Read a calibration TSV into (header, list-of-row-dicts)."""
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(dict(zip(header, line.split("\t"))))
    return header, rows


def write_calibrations(path, header, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(h, ".")) for h in header) + "\n")


def heldout_targets(rows):
    """Calibrations eligible to hold out: fixed tips with a positive age (the
    attested ancient languages). Returns [(target_id, true_age_bp), ...]."""
    out = []
    for r in rows:
        if r.get("type") == "fixed_tip":
            try:
                age = float(r["p1"])
            except (KeyError, ValueError):
                continue
            if age > 0.0:
                out.append((r["target_id"], age))
    return out


def _free_upper(rows, fallback=10000.0):
    """A generous upper bound for a freed tip's uniform prior: the largest age
    mentioned by any calibration (parent propagation constrains it further)."""
    hi = 0.0
    for r in rows:
        for k in ("p1", "p2"):
            try:
                hi = max(hi, float(r[k]))
            except (KeyError, ValueError):
                pass
    return hi if hi > 0.0 else fallback


def free_target(rows, target_id, upper=None):
    """Return a new calibration list with `target_id` freed: its fixed_tip row
    becomes a broad uniform_node prior [0, upper] so the sampler infers its age.
    """
    up = upper if upper is not None else _free_upper(rows)
    out = []
    freed = False
    for r in rows:
        if r["target_id"] == target_id and r.get("type") == "fixed_tip":
            nr = dict(r)
            nr["type"] = "uniform_node"
            nr["p1"] = "0"
            nr["p2"] = repr(up)
            nr["p3"] = "."
            nr["is_hard"] = "0"
            nr["source_key"] = "HELDOUT"
            out.append(nr)
            freed = True
        else:
            out.append(dict(r))
    if not freed:
        raise ValueError("no fixed_tip calibration for %r" % target_id)
    return out


def parse_posterior_nodes(path):
    """posterior_nodes.tsv -> {node: {median, hpd_lo, hpd_hi, rhat, ess, ...}}."""
    out = {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            node = parts[idx["node"]]

            def g(col):
                try:
                    return float(parts[idx[col]])
                except (KeyError, ValueError):
                    return float("nan")
            out[node] = {
                "median": g("median_bp"), "mean": g("mean_bp"),
                "hpd_lo": g("hpd90_lo_bp"), "hpd_hi": g("hpd90_hi_bp"),
                # "ess" before 0.3; "ess_bulk" since the tail ESS was added.
                "rhat": g("rhat"),
                "ess": g("ess_bulk") if "ess_bulk" in idx else g("ess"),
                "prior_lo": g("prior_hpd90_lo_bp"),
                "prior_hi": g("prior_hpd90_hi_bp"),
            }
    return out


def cli_runner(binary, tree, forms, settings):
    """Return a runner(cal_path, out_dir) -> posterior_nodes dict that shells out
    to the real `date-tree` command."""
    def run(cal_path, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        cmd = [binary, "date-tree", "--tree", tree, "--forms", forms,
               "--calibrations", cal_path, "--out", out_dir,
               "--warmup", str(settings["warmup"]),
               "--samples", str(settings["samples"]),
               "--chains", str(settings["chains"]),
               "--seed", str(settings["seed"]),
               "--rate-logmean", str(settings["rate_logmean"]),
               "--rate-logsd", str(settings["rate_logsd"]),
               "--allow-unconverged", "1"]
        cmd += ["--relaxed", "1" if settings.get("relaxed") else "0"]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return parse_posterior_nodes(os.path.join(out_dir, "posterior_nodes.tsv"))
    return run


def loco_one(rows, target_id, true_age, runner, out_dir):
    """Free `target_id`, run the tree, score the held-out node's posterior."""
    freed = free_target(rows, target_id)
    os.makedirs(out_dir, exist_ok=True)
    cal_path = os.path.join(out_dir, "cal_free_%s.tsv" % target_id)
    header = list(rows[0].keys())
    write_calibrations(cal_path, header, freed)
    post = runner(cal_path, os.path.join(out_dir, "run_%s" % target_id))
    node = post.get(target_id)
    if node is None:
        raise KeyError("freed node %r not in posterior output" % target_id)
    lo, hi = node["hpd_lo"], node["hpd_hi"]
    covered = 1 if (lo <= true_age <= hi) else 0
    plo = node.get("prior_lo", float("nan"))
    phi = node.get("prior_hi", float("nan"))
    is_post = interval_score(lo, hi, true_age)
    is_prior = interval_score(plo, phi, true_age)
    return {
        "target": target_id, "true_bp": true_age,
        "median_bp": node["median"], "hpd90_lo": lo, "hpd90_hi": hi,
        "covered": covered, "signed_err": node["median"] - true_age,
        "abs_err": abs(node["median"] - true_age), "width": hi - lo,
        "interval_score": is_post, "prior_lo": plo, "prior_hi": phi,
        "prior_interval_score": is_prior,
        "rhat": node["rhat"], "ess": node["ess"],
    }


def _scores(records):
    """Mean interval scores and skill (NaN-safe on the prior baseline)."""
    n = len(records)
    if not n:
        return float("nan"), float("nan"), float("nan")
    mis = sum(r["interval_score"] for r in records) / n
    mip = sum(r["prior_interval_score"] for r in records) / n
    skill = 1.0 - mis / mip if mip > 0 else float("nan")
    return mis, mip, skill


def run_benchmark(rows, runner, out_dir, targets=None):
    """LOCO over every held-out-eligible calibration. Returns (records, summary)."""
    if targets is None:
        targets = heldout_targets(rows)
    records = [loco_one(rows, tid, age, runner, out_dir) for (tid, age) in targets]
    n = len(records)
    cov = sum(r["covered"] for r in records) / n if n else float("nan")
    rmse = math.sqrt(sum(r["signed_err"]**2 for r in records) / n) if n else float("nan")
    mae = sum(r["abs_err"] for r in records) / n if n else float("nan")
    mean_w = sum(r["width"] for r in records) / n if n else float("nan")
    mis, mip, skill = _scores(records)
    summary = {"n": n, "coverage": cov, "rmse": rmse, "mae": mae,
               "mean_width": mean_w, "interval_score": mis,
               "prior_interval_score": mip, "skill": skill}
    return records, summary


REPORT_COLS = ["target", "true_bp", "median_bp", "hpd90_lo", "hpd90_hi",
               "covered", "signed_err", "abs_err", "width", "interval_score",
               "prior_lo", "prior_hi", "prior_interval_score", "rhat", "ess"]


def write_report(path, records):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(REPORT_COLS) + "\n")
        for r in records:
            f.write("\t".join(_fmt(r[c], c) for c in REPORT_COLS) + "\n")


def _fmt(x, col=None):
    if col == "rhat" and isinstance(x, float):
        return "%.3f" % x
    if isinstance(x, float):
        return "%.1f" % x if abs(x) >= 1 or x == 0 else "%.3f" % x
    return str(x)


def load_manifest(path, root=None):
    """Read a JSON benchmark manifest into (families, defaults).

    Manifest: {"defaults": {warmup,samples,chains,seed,relaxed}, "families":
    [{label, tree, forms, calibrations, rate_logmean, rate_logsd, [relaxed],
      [citation]}, ...]}. Paths are resolved relative to `root` (default: the
    manifest's own directory), so a manifest is portable.
    """
    with open(path, encoding="utf-8") as f:
        man = json.load(f)
    if root is None:
        root = os.path.dirname(os.path.abspath(path))
    defaults = man.get("defaults", {})
    fams = []
    for fam in man["families"]:
        f = dict(fam)
        for key in ("tree", "forms", "calibrations"):
            f[key] = os.path.join(root, f[key])
        fams.append(f)
    return fams, defaults


def run_family(fam, defaults, out_dir, runner_factory=cli_runner):
    """LOCO for one manifest family. Returns (label, records, summary).

    runner_factory(binary, tree, forms, settings) -> runner; injectable so tests
    can supply a fake (canned posterior) instead of shelling out to the engine.
    """
    settings = {
        "warmup": fam.get("warmup", defaults.get("warmup", 5000)),
        "samples": fam.get("samples", defaults.get("samples", 5000)),
        "chains": fam.get("chains", defaults.get("chains", 4)),
        "seed": fam.get("seed", defaults.get("seed", 42)),
        "rate_logmean": fam.get("rate_logmean", defaults.get("rate_logmean", -8.8)),
        "rate_logsd": fam.get("rate_logsd", defaults.get("rate_logsd", 0.8)),
        "relaxed": fam.get("relaxed", defaults.get("relaxed", 0)),
    }
    _hdr, rows = parse_calibrations(fam["calibrations"])
    runner = runner_factory(fam.get("binary"), fam["tree"], fam["forms"], settings)
    fam_out = os.path.join(out_dir, fam["label"])
    records, summary = run_benchmark(rows, runner, fam_out)
    for r in records:
        r["family"] = fam["label"]
    return fam["label"], records, summary


def run_multi(families, defaults, out_dir, runner_factory=cli_runner):
    """Run every manifest family; return (all_records, per_family, overall)."""
    all_records = []
    per_family = {}
    for fam in families:
        label, records, summary = run_family(fam, defaults, out_dir, runner_factory)
        per_family[label] = summary
        all_records.extend(records)
    n = len(all_records)
    overall = {
        "n": n,
        "coverage": sum(r["covered"] for r in all_records) / n if n else float("nan"),
        "rmse": math.sqrt(sum(r["signed_err"]**2 for r in all_records) / n) if n else float("nan"),
        "mae": sum(r["abs_err"] for r in all_records) / n if n else float("nan"),
    }
    overall["interval_score"], overall["prior_interval_score"], \
        overall["skill"] = _scores(all_records)
    return all_records, per_family, overall


def write_multi_report(path, records):
    cols = ["family"] + REPORT_COLS
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in records:
            f.write("\t".join(_fmt(r[c], c) for c in cols) + "\n")


def _run_single(args):
    _hdr, rows = parse_calibrations(args.calibrations)
    os.makedirs(args.out, exist_ok=True)
    settings = {"warmup": args.warmup, "samples": args.samples,
                "chains": args.chains, "seed": args.seed,
                "rate_logmean": args.rate_logmean, "rate_logsd": args.rate_logsd,
                "relaxed": args.relaxed}
    runner = cli_runner(args.binary, args.tree, args.forms, settings)
    targets = heldout_targets(rows)
    if not targets:
        print("No held-out-eligible (positive fixed_tip) calibrations found.",
              file=sys.stderr)
        return 2
    print("Leave-one-calibration-out over %d attested tips..." % len(targets))
    records, summary = run_benchmark(rows, runner, args.out)
    report_path = os.path.join(args.out, "loco_report.tsv")
    write_report(report_path, records)
    print("%-22s %8s %8s %8s %8s %6s %5s %5s"
          % ("target", "true", "median", "hpd_lo", "hpd_hi", "cov", "rhat", "ess"))
    for r in records:
        print("%-22s %8.0f %8.0f %8.0f %8.0f %6d %5.2f %5.0f"
              % (r["target"], r["true_bp"], r["median_bp"], r["hpd90_lo"],
                 r["hpd90_hi"], r["covered"], r["rhat"], r["ess"]))
    print("\nN=%d  coverage=%.2f (90%% HPD)  RMSE=%.0f BP  MAE=%.0f BP  "
          "mean width=%.0f BP" % (summary["n"], summary["coverage"],
                                  summary["rmse"], summary["mae"],
                                  summary["mean_width"]))
    print("Interval score %.0f BP (prior only %.0f BP, skill %.2f)"
          % (summary["interval_score"], summary["prior_interval_score"],
             summary["skill"]))
    print("Report: %s" % report_path)
    print("Read as a small-N calibrated-ness spot check, not a coverage "
          "guarantee; deep internal nodes stay rate-prior-dominated.")
    return 0


def _run_manifest(args):
    families, defaults = load_manifest(args.manifest)
    if args.relaxed:
        defaults["relaxed"] = args.relaxed
    for fam in families:
        fam["binary"] = args.binary
    os.makedirs(args.out, exist_ok=True)
    print("Multi-family leave-one-calibration-out over %d families%s..."
          % (len(families), " (relaxed clock)" if defaults.get("relaxed") else ""))
    records, per_family, overall = run_multi(families, defaults, args.out)
    report_path = os.path.join(args.out, "loco_multi_report.tsv")
    write_multi_report(report_path, records)
    print("%-16s %-20s %6s %6s %6s %5s %4s %5s"
          % ("family", "target", "true", "median", "cov", "rhat", "ess", ""))
    for r in records:
        print("%-16s %-20s %6.0f %6.0f %6d %5.2f %4.0f"
              % (r["family"], r["target"], r["true_bp"], r["median_bp"],
                 r["covered"], r["rhat"], r["ess"]))
    print("\nPer family (coverage @ 90% HPD):")
    for label, s in per_family.items():
        print("  %-16s N=%d  coverage=%.2f  RMSE=%.0f BP  IS=%.0f  skill=%.2f"
              % (label, s["n"], s["coverage"], s["rmse"],
                 s["interval_score"], s["skill"]))
    print("OVERALL  N=%d  coverage=%.2f  RMSE=%.0f BP  MAE=%.0f BP  "
          "IS=%.0f (prior %.0f, skill %.2f)"
          % (overall["n"], overall["coverage"], overall["rmse"], overall["mae"],
             overall["interval_score"], overall["prior_interval_score"],
             overall["skill"]))
    print("Report: %s" % report_path)
    return 0


def main(argv):
    ap = argparse.ArgumentParser(prog="validation.calibration_benchmark")
    ap.add_argument("--binary", required=True)
    ap.add_argument("--manifest", help="JSON multi-family benchmark manifest")
    ap.add_argument("--tree")
    ap.add_argument("--forms")
    ap.add_argument("--calibrations")
    ap.add_argument("--out", required=True)
    ap.add_argument("--warmup", type=int, default=6000)
    ap.add_argument("--samples", type=int, default=6000)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--rate-logmean", type=float, default=-8.8)
    ap.add_argument("--rate-logsd", type=float, default=0.8)
    ap.add_argument("--relaxed", type=int, default=0,
                    help="use the relaxed (per-branch) clock (default strict)")
    args = ap.parse_args(argv[1:])

    if args.manifest:
        return _run_manifest(args)
    if not (args.tree and args.forms and args.calibrations):
        ap.error("single-family mode needs --tree, --forms, --calibrations "
                 "(or pass --manifest)")
    return _run_single(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

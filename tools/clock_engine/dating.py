"""Run `clepsydra date-tree` on one subtree and read the target node.

Strictly one run at a time. Each run is wrapped in a transient systemd scope
with a memory cap. A run that fails the convergence check (exit status 3) is
rerun once with twice the samples (plan.md, "Dating").
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import time
from dataclasses import dataclass

WARMUP = 3000
SAMPLES = 3000
CHAINS = 4
SEED = 20261010
Z90 = 1.6448536269514722
SCOPE = ["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=1G",
         "-p", "MemorySwapMax=0"]


@dataclass
class Dated:
    node: str
    out: str                 # clepsydra output directory used
    converged: bool
    reruns: int
    L_median: float
    L_lo: float
    L_hi: float
    ln_L: float
    s_L: float
    manifest: dict
    target: dict
    failing: list


def failing_params(out):
    """Parameters that fail the engine's convergence rule (R-hat < 1.01, ESS >= 100 per chain)."""
    bad = []
    with open(os.path.join(out, "diagnostics.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                rh, eb, et = float(r["rhat"]), float(r["ess_bulk"]), float(r["ess_tail"])
            except (KeyError, ValueError):
                continue
            if rh >= 1.01 or eb < 100 * CHAINS or et < 100 * CHAINS:
                bad.append(r["param"])
    return bad


def _run(binary, sub, out, samples, scope, log):
    cmd = [binary, "date-tree", "--tree", os.path.join(sub, "tree.nwk"),
           "--forms", os.path.join(sub, "forms.tsv"),
           "--calibrations", os.path.join(sub, "calibrations.tsv"),
           "--out", out, "--warmup", str(WARMUP), "--samples", str(samples),
           "--chains", str(CHAINS), "--seed", str(SEED), "--relaxed", "0"]
    if scope:
        cmd = SCOPE + cmd
    t0 = time.time()
    with open(out + ".log", "w", encoding="utf-8") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    log(f"    {os.path.basename(out)}: exit {rc}, {dt:.0f} s")
    if rc not in (0, 3):
        raise SystemExit(f"clepsydra failed (exit {rc}); see {out}.log")
    return rc, dt


def date_node(binary, sub, scope, force, log):
    """Run (or reuse) the dating of one subtree; record wall times in run.tsv."""
    rec = os.path.join(sub, "run.tsv")
    if os.path.exists(rec) and not force:
        return
    rows = []
    rc, dt = _run(binary, sub, os.path.join(sub, "out"), SAMPLES, scope, log)
    rows.append(("out", SAMPLES, rc, dt))
    if rc == 3:
        rc, dt = _run(binary, sub, os.path.join(sub, "out_rerun"), 2 * SAMPLES, scope, log)
        rows.append(("out_rerun", 2 * SAMPLES, rc, dt))
    with open(rec, "w", encoding="utf-8") as fh:
        fh.write("run\tsamples\texit\twall_s\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.1f}\n")


def read_runs(sub):
    with open(os.path.join(sub, "run.tsv"), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read(node_id, sub) -> Dated:
    runs = read_runs(sub)
    last = runs[-1]
    out = os.path.join(sub, last["run"])
    with open(os.path.join(out, "run_manifest.tsv"), encoding="utf-8") as fh:
        man = {r[0]: r[1] for r in csv.reader(fh, delimiter="\t") if len(r) >= 2}
    with open(os.path.join(out, "posterior_nodes.tsv"), encoding="utf-8") as fh:
        tgt = next(r for r in csv.DictReader(fh, delimiter="\t") if r["node"] == "target")
    med = float(tgt["lexical_median"])
    lo, hi = float(tgt["lexical_hpd90_lo"]), float(tgt["lexical_hpd90_hi"])
    s = (math.log(hi) - math.log(lo)) / (2 * Z90)
    return Dated(node_id, out, man.get("converged") == "1", len(runs) - 1,
                 med, lo, hi, math.log(med), s, man, tgt,
                 failing_params(out))

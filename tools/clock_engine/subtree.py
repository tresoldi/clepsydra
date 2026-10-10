"""Target-like subtrees around each calibrated node (plan.md, "Subtree per node").

Node selection is the pilot's: `select_nodes` streams the expert-cognate
datasets one at a time through `clock_pilot.match` and applies the pilot's
`family` rule. `build` then reloads a node's primary dataset and writes
`tree.nwk`, `forms.tsv` and `calibrations.tsv`.
"""

from __future__ import annotations

import csv
import os
import re
import atexit
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field

from clock_pilot import calibrations as C
from clock_pilot import lexibank as L
from clock_pilot import match as M
from clock_pilot.glottolog import Glottolog
from clock_pilot.run import select

MAX_INGROUP = 6
MAX_PER_SIDE = 3
MAX_OUTGROUP = 2
# Deviation 1 (plan.md): historical languages the pilot's name rule misses.
# They would enter a subtree as modern tips (fixed at age 0).
EXTRA_ANCIENT = {"jurc1239": "Jurchen (12th-17th c.)"}
PILOT_MATCHES = os.path.join(C.REPO, "docs", "research", "clock-pilot",
                             "calibration_matches.tsv")


@dataclass
class Subtree:
    node: C.Node
    dataset: str
    ingroup: list[str]
    outgroup: list[str]
    tips: dict[str, str]                 # doculect id -> tip name
    newick: str
    target_is_root: bool
    names: dict[str, str] = field(default_factory=dict)
    n_concepts: int = 0
    n_rows: int = 0


# The calibration table as the pilot read it. D1.2 (7f38b1f) later changed the
# grade of two rows the pilot excluded, which its loader rejects.
PILOT_CALIB_COMMIT = "f92ccc1"


def pin_calibrations():
    """Point clock_pilot.calibrations at the pilot's calibrations.tsv."""
    if getattr(C, "_pinned", False):
        return
    txt = subprocess.run(["git", "-C", C.REPO, "show",
                          f"{PILOT_CALIB_COMMIT}:reference/calibrations/calibrations.tsv"],
                         capture_output=True, check=True).stdout
    fd, path = tempfile.mkstemp(suffix=".tsv", prefix="calibrations-")
    with os.fdopen(fd, "wb") as fh:
        fh.write(txt)
    atexit.register(os.unlink, path)
    C.CALIB_TSV = path
    C._pinned = True


def select_nodes(g: Glottolog, log):
    """Calibrated nodes and their primary Match, exactly as the pilot chose them."""
    pin_calibrations()
    cals, _ = C.load(g)
    nodes, _ = C.dedupe(cals)
    matches = defaultdict(list)
    for name in L.expert_datasets():
        ds = L.load(name, g.name)
        if ds.sets_per_concept < 1.05:          # pilot deviation 1
            log(f"  skip {name}: degenerate coding")
            del ds
            continue
        for n in nodes:
            m = M.compute(n, ds, g)
            if m is not None:
                matches[n.id].append(m)
        del ds
    chosen, fam_primary = select(nodes, matches, "family")
    check_against_pilot(nodes, chosen)
    return nodes, chosen, fam_primary


def check_against_pilot(nodes, chosen):
    with open(PILOT_MATCHES, encoding="utf-8") as fh:
        prim = {r["calib_id"]: r["dataset"] for r in csv.DictReader(fh, delimiter="\t")
                if r["role"] == "primary"}
    mine = {n.id: chosen[n.id].dataset for n in nodes if chosen.get(n.id) is not None}
    if mine != prim:
        raise SystemExit(f"node selection differs from the pilot:\n mine={mine}\n pilot={prim}")


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", s)


def _clean(s: str) -> str:
    return re.sub(r"\s+", "_", s.strip())


def choose_doculects(m: M.Match, ds: L.Dataset, g: Glottolog):
    node = m.node
    count = {d: len(ds.table.get(d, {})) for d in ds.doculects}
    rank = lambda d: (-count[d], d)
    groups = defaultdict(list)
    for d, lab in m.labels.items():
        if ds.doculects[d].glottocode not in EXTRA_ANCIENT:
            groups[lab].append(d)
    for v in groups.values():
        v.sort(key=rank)
    if node.spec.kind == "crown":
        order = sorted(groups, key=lambda k: (-count[groups[k][0]], k))
        ingroup, k = [], 0
        while len(ingroup) < MAX_INGROUP and any(len(groups[c]) > k for c in order):
            for c in order:
                if len(groups[c]) > k and len(ingroup) < MAX_INGROUP:
                    ingroup.append(groups[c][k])
            k += 1
    else:
        ingroup = groups.get("A", [])[:MAX_PER_SIDE] + groups.get("B", [])[:MAX_PER_SIDE]
    parent = g.parent(node.mrca)
    out = []
    if parent is not None:
        for d in ds.doculects.values():
            gc = d.glottocode
            if (d.ancient or gc in EXTRA_ANCIENT or not gc or not g.known(gc) or not ds.table.get(d.id)
                    or g.is_under(gc, node.mrca)):
                continue
            if g.mrca([gc, node.mrca]) == parent:
                out.append(d.id)
    out.sort(key=rank)
    return ingroup, out[:MAX_OUTGROUP]


def newick(docs, ds: L.Dataset, g: Glottolog, target: str, tips: dict[str, str]):
    """Glottolog classification pruned to `docs`; unary nodes collapsed."""
    trie: dict = {}
    for d in docs:
        cur = trie
        for gc in g.lineage[ds.doculects[d].glottocode] + ("\x00" + d,):
            cur = cur.setdefault(gc, {})
    found = []

    def emit(key, sub):
        if key.startswith("\x00"):
            return tips[key[1:]]
        kids = [emit(k, v) for k, v in sorted(sub.items())]
        if len(kids) == 1:
            return kids[0]
        lab = "target" if key == target else key
        if key == target:
            found.append(key)
        return "(" + ",".join(kids) + ")" + lab

    # collapse the shared stem above the MRCA
    key, sub = "", trie
    while len(sub) == 1 and not next(iter(sub)).startswith("\x00"):
        key, sub = next(iter(sub.items()))
    s = emit(key, sub) + ";"
    if found != [target]:
        raise SystemExit(f"target {target} is not a node of the pruned tree: {s}")
    return s, key == target


def build(n: C.Node, m: M.Match, g: Glottolog, outdir: str) -> Subtree:
    ds = L.load(m.dataset, g.name)
    ingroup, outgroup = choose_doculects(m, ds, g)
    gcs = [ds.doculects[d].glottocode for d in ingroup]
    if g.mrca(gcs) != n.mrca:
        raise SystemExit(f"{n.id}: ingroup MRCA {g.mrca(gcs)} != node {n.mrca}")
    docs = ingroup + outgroup
    tips, used = {}, set()
    for d in docs:
        t = base = "L_" + _safe(d)
        i = 2
        while t in used:
            t, i = f"{base}_{i}", i + 1
        used.add(t)
        tips[d] = t
    nwk, is_root = newick(docs, ds, g, n.mrca, tips)
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "tree.nwk"), "w", encoding="utf-8") as fh:
        fh.write(nwk + "\n")
    rows, concepts = 0, set()
    with open(os.path.join(outdir, "forms.tsv"), "w", encoding="utf-8") as fh:
        fh.write("language_id\tconcept_id\tcognate_id\n")
        for d in docs:
            for c, cs in sorted(ds.table[d].items()):
                fh.write(f"{tips[d]}\t{_clean(c)}\t{_clean(cs)}\n")
                rows += 1
                concepts.add(c)
    with open(os.path.join(outdir, "calibrations.tsv"), "w", encoding="utf-8") as fh:
        fh.write("calibration_id\ttarget_id\ttype\tp1\tp2\tp3\tis_hard\tsource_key\n")
        for d in docs:
            fh.write(f"tip_{tips[d]}\t{tips[d]}\tfixed_tip\t0\t.\t.\t1\tmodern_list\n")
    with open(os.path.join(outdir, "doculects.tsv"), "w", encoding="utf-8") as fh:
        fh.write("tip\tdoculect_id\tname\tglottocode\trole\tside\tn_concepts\n")
        for d in docs:
            dl = ds.doculects[d]
            role = "ingroup" if d in ingroup else "outgroup"
            side = m.labels.get(d, ".") if role == "ingroup" else "."
            fh.write(f"{tips[d]}\t{d}\t{dl.name}\t{dl.glottocode}\t{role}\t{side}\t"
                     f"{len(ds.table[d])}\n")
    names = {d: ds.doculects[d].name for d in docs}
    st = Subtree(n, m.dataset, ingroup, outgroup, tips, nwk, is_root, names,
                 len(concepts), rows)
    del ds
    return st

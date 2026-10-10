"""Calibration rows, their curation, and deduplication by node identity.

`curation.tsv` (next to this file) holds the per-row decisions of the
pre-registration (`docs/research/clock-pilot/plan.md`): how the row's age is
used (interval / point / max / min / exclude) and which Glottolog node it
dates (`crown:X`, `split:A|B`, `stem:X`).
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field

from .glottolog import Glottolog

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
CALIB_TSV = os.path.join(REPO, "reference", "calibrations", "calibrations.tsv")
CURATION_TSV = os.path.join(HERE, "curation.tsv")

KINDS = {"node", "root", "stem", "settlement_bound"}
TIP_OFFSET = 50.0          # BP1950 -> years before ~2000 CE (collection date)
POINT_LOG_SD = 0.15        # Holman's 15% for point calibrations
USE_RANK = {"interval": 0, "point": 1, "max": 2, "min": 2}


@dataclass
class Side:
    include: list[str]
    exclude: list[str] = field(default_factory=list)

    def contains(self, g: Glottolog, gc: str) -> bool:
        return (any(g.is_under(gc, a) for a in self.include)
                and not any(g.is_under(gc, x) for x in self.exclude))


def parse_side(text: str) -> Side:
    inc, exc = [], []
    for tok in text.split(","):
        parts = tok.split("-")
        inc.append(parts[0].strip())
        exc.extend(p.strip() for p in parts[1:])
    return Side(inc, exc)


@dataclass
class NodeSpec:
    kind: str                       # crown / split / stem
    text: str
    a: Side | None = None
    b: Side | None = None
    x: str | None = None

    def resolve(self, g: Glottolog) -> tuple[str | None, list[str]]:
        """(MRCA glottocode, codes defining the node) or (None, missing codes)."""
        if self.kind == "crown":
            return (self.x, []) if g.known(self.x) else (None, [self.x])
        if self.kind == "stem":
            if not g.known(self.x):
                return None, [self.x]
            return g.parent(self.x), []
        codes = self.a.include + self.b.include
        missing = [c for c in codes if not g.known(c)]
        if missing:
            return None, missing
        return g.mrca(codes), []

    def sides(self, g: Glottolog, gc: str):
        """Return 'A', 'B' or None for a doculect glottocode."""
        if not g.known(gc):
            return None
        if self.kind == "crown":
            ch = g.child_under(gc, self.x)
            return ch  # the child of X; pairs need different children
        if self.kind == "stem":
            if g.is_under(gc, self.x):
                return "A"
            p = g.parent(self.x)
            if p and g.is_under(gc, p):
                return "B"
            return None
        if self.a.contains(g, gc):
            return "A"
        if self.b.contains(g, gc):
            return "B"
        return None


def parse_node(text: str) -> NodeSpec | None:
    if text in (".", ""):
        return None
    kind, rest = text.split(":", 1)
    if kind in ("crown", "stem"):
        return NodeSpec(kind, text, x=rest)
    a, b = rest.split("|")
    return NodeSpec("split", text, a=parse_side(a), b=parse_side(b))


def _num(s: str) -> float | None:
    s = (s or "").strip()
    return None if s in ("", ".") else float(s)


@dataclass
class Calib:
    id: str
    row: dict
    use: str
    spec: NodeSpec | None
    lo: float | None
    hi: float | None
    center: float | None
    reason: str
    grade: str
    order: int
    mrca: str | None = None
    family: str | None = None
    note: str = ""


def load(g: Glottolog) -> tuple[list[Calib], list[dict]]:
    """All selected rows (with curation), and a per-row log of decisions."""
    with open(CALIB_TSV, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(CURATION_TSV, encoding="utf-8") as fh:
        cur = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    out = []
    for i, r in enumerate(rows):
        if r["grade"] not in ("A", "B") or r["target_kind"] not in KINDS:
            continue
        c = cur.get(r["id"])
        if c is None:
            raise SystemExit(f"curation.tsv lacks {r['id']}")
        cal = Calib(r["id"], r, c["use"], parse_node(c["node"]), _num(c["lo"]),
                    _num(c["hi"]), _num(c["center"]), c["reason"], r["grade"], i)
        if cal.spec is not None:
            mrca, missing = cal.spec.resolve(g)
            if mrca is None:
                cal.note = "glottocode not in Glottolog: " + ",".join(missing)
                if cal.use != "exclude":
                    cal.use = "exclude"
            else:
                cal.mrca = mrca
                cal.family = g.top(mrca)
        out.append(cal)
    if len(out) != len(cur):
        extra = set(cur) - {c.id for c in out}
        raise SystemExit(f"curation rows not in selection: {sorted(extra)}")
    return out, rows


@dataclass
class Node:
    """One calibrated node after deduplication."""
    mrca: str
    family: str
    primary: Calib
    members: list[Calib]
    use: str                  # interval / point / max / min
    lo: float | None
    hi: float | None
    center: float | None
    spec: NodeSpec
    combined: str = ""

    @property
    def id(self) -> str:
        return self.primary.id

    def log_moments(self):
        if self.use == "interval":
            a = math.log(self.lo + TIP_OFFSET)
            b = math.log(self.hi + TIP_OFFSET)
            return 0.5 * (a + b), (b - a) / math.sqrt(12.0)
        if self.use == "point":
            return math.log(self.center + TIP_OFFSET), POINT_LOG_SD
        raise ValueError(self.use)

    @property
    def two_sided(self) -> bool:
        return self.use in ("interval", "point")

    def bound(self) -> float:
        """One-sided bound in years before ~2000 CE."""
        v = self.hi if self.use == "max" else self.lo
        if v is None:
            v = self.center
        return v + TIP_OFFSET


def dedupe(cals: list[Calib]) -> tuple[list[Node], dict[str, str]]:
    """Group usable rows by MRCA; returns nodes and a row-id -> status map."""
    status: dict[str, str] = {}
    groups: dict[str, list[Calib]] = {}
    for c in cals:
        if c.use == "exclude":
            status[c.id] = "excluded: " + (c.note or c.reason)
            continue
        groups.setdefault(c.mrca, []).append(c)
    nodes = []
    for mrca, members in groups.items():
        members.sort(key=lambda c: (USE_RANK[c.use], c.grade != "A", c.order))
        p = members[0]
        use, lo, hi, center, comb = p.use, p.lo, p.hi, p.center, ""
        if use in ("max", "min"):
            mins = [c for c in members if c.use == "min"]
            maxs = [c for c in members if c.use == "max"]
            if mins and maxs and mins[0].lo < maxs[0].hi:
                use, lo, hi = "interval", mins[0].lo, maxs[0].hi
                comb = f"interval from {mins[0].id} (min) + {maxs[0].id} (max)"
        n = Node(mrca, p.family, p, members, use, lo, hi, center, p.spec, comb)
        nodes.append(n)
        for c in members:
            status[c.id] = ("primary" if c is p
                            else f"duplicate of {p.id} (same MRCA {mrca})")
            if comb and c.id in comb and c is not p:
                status[c.id] += "; combined into interval"
    nodes.sort(key=lambda n: n.primary.order)
    return nodes, status

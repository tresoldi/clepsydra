"""Benchmark targets: calibrated nodes seen through GLED and ASJP."""
from dataclasses import dataclass, field

from clock_pilot.calibrations import TIP_OFFSET, dedupe, load
from clock_pilot.glottolog import Glottolog

MIN_LANGS = 4          # tips needed to date a node
MIN_CONCEPTS = 25      # concepts a doculect must record in GLED
MAX_LANGS = 24         # tips per run


@dataclass
class Target:
    id: str
    family: str
    node: str
    use: str                    # interval / point / max / min
    lo: float | None            # years before ~2000
    hi: float | None
    gc_all: list = field(default_factory=list)   # every GLED code under the node
    gc: list = field(default_factory=list)       # the sampled codes
    side: dict = field(default_factory=dict)     # code -> 'A' / 'B' / child

    @property
    def truth_lo(self):
        return self.lo if self.lo is not None else 0.0


def build(gled, g=None, max_langs=MAX_LANGS, min_langs=MIN_LANGS):
    from .data import sample
    g = g or Glottolog()
    cals, _ = load(g)
    nodes, _ = dedupe(cals)
    out = []
    for n in nodes:
        spec = n.spec
        if spec is None or spec.kind == "stem":
            continue
        side = {}
        for gc in gled.by_gc:
            s = spec.sides(g, gc)
            if s is None:
                continue
            lid = gled.best(gc)
            if lid is None or len(gled.forms[lid]) < MIN_CONCEPTS:
                continue
            side[gc] = s
        if len(set(side.values())) < 2 or len(side) < min_langs:
            continue
        # keep both sides represented in the sample
        by = {}
        for gc, s in side.items():
            by.setdefault(s, []).append(gc)
        groups = sorted(by, key=lambda s: -len(by[s]))
        picked = []
        per = max(1, max_langs // len(groups)) if len(groups) <= max_langs else 1
        if len(groups) > max_langs:
            groups = groups[:max_langs]
        for s in groups:
            picked.extend(sample(by[s], g, per))
        picked = picked[:max_langs]
        lo = None if n.lo is None else n.lo + TIP_OFFSET
        hi = None if n.hi is None else n.hi + TIP_OFFSET
        if n.use == "point":
            lo, hi = n.center * 0.85 + TIP_OFFSET, n.center * 1.15 + TIP_OFFSET
        out.append(Target(n.id, n.family, n.primary.row.get("target", n.mrca),
                          n.use, lo, hi, sorted(side), picked,
                          {gc: side[gc] for gc in picked}))
    return out


# Published family-root estimates (docs/research/literature-lexical-dating.md,
# Part 2; status "paper" unless noted). window = union of the reported 95%
# intervals, centre = geometric mean of the central estimates, both in years
# before 2000 as published (not shifted: Bayesian studies quote "years ago").
ROOTS = [
    # id, glottocode, centre estimates, window, sources, ASJP published (Holman 2011)
    ("ROOT-IE", "indo1319", [5650, 5950, 6410, 7580, 8120, 8700], (4870, 10410),
     "Kassian 2021; Chang 2015; Bouckaert 2013; Heggarty 2023; Gray & Atkinson 2003", 4348),
    ("ROOT-URA", "ural1272", [5113], (3330, 7503), "Honkola 2013 (sec.)", 3178),
    ("ROOT-TRK", "turk1311", [2066], (1517, 2755), "Savelyev & Robbeets 2020", 3404),
    ("ROOT-SEM", "semi1276", [5750], (4400, 7400), "Kitchen 2009", None),
    ("ROOT-JPN", "japo1237", [2182], (1239, 4190), "Lee & Hasegawa 2011 (sec.)", 1564),
    ("ROOT-ST", "sino1245", [5871, 7184, 7983], (4183, 11285),
     "Zhang 2019; Sagart 2019; Zhang 2020", 5261),
    ("ROOT-DRA", "drav1251", [4433, 4650], (3000, 6500), "Kolipakam 2018", 2055),
    ("ROOT-AN", "aust1307", [5117, 5230], (4660, 5800), "Gray 2009; Greenhill 2010", 3633),
    ("ROOT-PN", "pama1250", [5671], (4455, 6966), "Bouckaert 2018 (one calibration)", 4295),
    ("ROOT-UA", "utoa1244", [4095], (3258, 5025), "Greenhill 2023", None),
    ("ROOT-BNT", "narr1281", [4846, 4940, 5110], (4500, 5770), "Grollemund 2015; Koile 2022", None),
    ("ROOT-TG", "tupi1276", [1762], (847, 2740), "Ferraz Gerardi 2023", 2043),
]


def roots(gled, g, max_langs=30):
    from .data import sample
    from math import exp, log
    out = []
    for rid, fam, cen, win, src, _asjp in ROOTS:
        side = {}
        for gc in gled.by_gc:
            if not g.is_under(gc, fam):
                continue
            lid = gled.best(gc)
            if lid is None or len(gled.forms[lid]) < MIN_CONCEPTS:
                continue
            side[gc] = g.child_under(gc, fam) or gc
        if len(side) < MIN_LANGS:
            continue
        picked = sample(list(side), g, max_langs)
        top = g.mrca(picked)               # the first branching node, not the family
        side = {gc: (g.child_under(gc, top) or gc) for gc in side}
        t = Target(rid, fam, src, "interval", float(win[0]), float(win[1]),
                   sorted(side), picked, {gc: side[gc] for gc in picked})
        t.center = exp(sum(log(x) for x in cen) / len(cen))
        out.append(t)
    return out

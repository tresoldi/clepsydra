"""GLED cognate data, ASJP word lists, and language sampling."""
import csv
import os
import re
from collections import defaultdict

GLED_DIR = os.path.expanduser("~/lexibank_clone/gled/cldf")
ASJP_DIR = os.environ.get("ASJP_CLDF", os.path.expanduser("~/scratch/asjp/cldf"))

# Doculects that are not modern: dropped, so every tip is dated 0.
ANCIENT = re.compile(
    r"(latin|sanskrit|vedic|hittite|old[ _]|middle[ _]|ancient|classical|gothic|"
    r"avest|akkad|ugarit|biblical|pali|tochar|luvian|lycian|oscan|umbrian|"
    r"ge.ez|aramaic|coptic|sumerian|etruscan|proto|reconstr)", re.I)


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


class Gled:
    """GLED: 40 concepts, automatic cognate sets, Glottolog codes."""

    def __init__(self, d=GLED_DIR):
        self.lang = {r["ID"]: r for r in read_csv(os.path.join(d, "languages.csv"))}
        self.forms = defaultdict(dict)       # lang -> concept -> cognate set
        self.segments = defaultdict(dict)
        for r in read_csv(os.path.join(d, "forms.csv")):
            if not r["Cognacy"]:
                continue
            self.forms[r["Language_ID"]].setdefault(r["Parameter_ID"], r["Cognacy"])
        self.concepts = sorted({c for f in self.forms.values() for c in f})
        self.by_gc = defaultdict(list)
        for lid, r in self.lang.items():
            if r["Glottocode"] and lid in self.forms and not ANCIENT.search(r["Name"] + " " + lid):
                self.by_gc[r["Glottocode"]].append(lid)

    def best(self, gc):
        """The doculect of a Glottolog code with the most concepts."""
        c = self.by_gc.get(gc)
        return max(c, key=lambda l: len(self.forms[l])) if c else None

    def distance(self, a, b):
        """1 - share of identical cognate sets over shared concepts (None if < 10)."""
        fa, fb = self.forms[a], self.forms[b]
        both = [c for c in fa if c in fb]
        if len(both) < 10:
            return None
        return 1.0 - sum(fa[c] == fb[c] for c in both) / len(both)


class Asjp:
    """ASJP v21 (CLDF): ASJP-coded words, 40-item subset matched by Concepticon."""

    def __init__(self, gled_param_csv, d=ASJP_DIR):
        gl = {r["Concepticon_ID"] for r in read_csv(gled_param_csv)}
        params = read_csv(os.path.join(d, "parameters.csv"))
        self.concept = {r["ID"]: r["Concepticon_ID"] for r in params
                        if r["Concepticon_ID"] in gl}
        self.lang = {r["ID"]: r for r in read_csv(os.path.join(d, "languages.csv"))}
        self.words = defaultdict(dict)       # lang -> concept id -> ASJP string
        for r in read_csv(os.path.join(d, "forms.csv")):
            p = r["Parameter_ID"]
            if p not in self.concept or r["Loan"] == "true":
                continue
            w = re.sub(r'[*"$~]', "", r["Value"].split(",")[0]).strip()
            if w and r["Language_ID"] not in ():
                self.words[r["Language_ID"]].setdefault(p, w)
        self.by_gc = defaultdict(list)
        for lid, r in self.lang.items():
            if (r["Glottocode"] and len(self.words.get(lid, ())) >= 28
                    and r["long_extinct"] != "True" and r["recently_extinct"] != "True"):
                self.by_gc[r["Glottocode"]].append(lid)

    def best(self, gc):
        c = self.by_gc.get(gc)
        return max(c, key=lambda l: len(self.words[l])) if c else None


def sample(codes, g, k):
    """Pick up to k glottocodes spread over the Glottolog subgroups: round-robin
    over the children of the deepest branching node, recursively."""
    codes = sorted(set(codes))
    if len(codes) <= k:
        return codes
    return _spread(codes, g, k)


def _spread(codes, g, k):
    if len(codes) <= k:
        return codes
    top = g.mrca(codes)
    groups = defaultdict(list)
    for c in codes:
        groups[g.child_under(c, top) or c].append(c)
    keys = sorted(groups, key=lambda x: (-len(groups[x]), x))
    quota = {x: 0 for x in keys}
    left = k
    while left > 0:
        moved = False
        for x in keys:
            if quota[x] < len(groups[x]) and left > 0:
                quota[x] += 1
                left -= 1
                moved = True
        if not moved:
            break
    out = []
    for x in keys:
        if quota[x]:
            out.extend(_spread(groups[x], g, quota[x]))
    return out

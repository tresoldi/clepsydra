"""Cross pairs across a calibrated split, and their cognate match counts."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .calibrations import Node
from .glottolog import Glottolog
from .lexibank import Dataset

MIN_MEDIAN_SHARED = 40


@dataclass
class Match:
    node: Node
    dataset: str
    side_a: list[str]            # doculect ids
    side_b: list[str]
    pairs: list[tuple[str, str]]
    concepts: list[str]
    K: np.ndarray                # per-concept matches summed over pairs
    N: np.ndarray                # per-concept comparisons summed over pairs
    pair_k: np.ndarray
    pair_n: np.ndarray
    labels: dict[str, str]       # doculect id -> side (A/B or crown child)

    @property
    def median_shared(self) -> float:
        return float(np.median(self.pair_n)) if len(self.pair_n) else 0.0

    @property
    def usable(self) -> bool:
        return len(self.pairs) > 0 and self.median_shared >= MIN_MEDIAN_SHARED


def sides(node: Node, ds: Dataset, g: Glottolog):
    """Doculect id -> side label for doculects that fall under the node."""
    lab = {}
    for d in ds.doculects.values():
        if d.ancient or not d.glottocode or d.id not in ds.table:
            continue
        s = node.spec.sides(g, d.glottocode)
        if s is not None:
            lab[d.id] = s
    return lab


def cross_pairs(node: Node, lab: dict[str, str]):
    ids = sorted(lab)
    if node.spec.kind == "crown":
        return [(a, b) for i, a in enumerate(ids) for b in ids[i + 1:]
                if lab[a] != lab[b]]
    A = [i for i in ids if lab[i] == "A"]
    B = [i for i in ids if lab[i] == "B"]
    return [(a, b) for a in A for b in B]


def compute(node: Node, ds: Dataset, g: Glottolog,
            concepts_filter: set[str] | None = None) -> Match | None:
    lab = sides(node, ds, g)
    pairs = cross_pairs(node, lab)
    if not pairs:
        return None
    docs = sorted({x for p in pairs for x in p})
    concepts = sorted({c for d in docs for c in ds.table[d]})
    if concepts_filter is not None:
        concepts = [c for c in concepts if c in concepts_filter]
    cidx = {c: i for i, c in enumerate(concepts)}
    # integer-code cognate classes per concept; -1 = missing
    M = -np.ones((len(docs), len(concepts)), dtype=np.int64)
    codes: dict[str, int] = {}
    for i, d in enumerate(docs):
        for c, cs in ds.table[d].items():
            j = cidx.get(c)
            if j is not None:
                M[i, j] = codes.setdefault(cs, len(codes))
    didx = {d: i for i, d in enumerate(docs)}
    ia = np.array([didx[a] for a, _ in pairs])
    ib = np.array([didx[b] for _, b in pairs])
    xa, xb = M[ia], M[ib]
    obs = (xa >= 0) & (xb >= 0)
    hit = obs & (xa == xb)
    if node.spec.kind == "crown":
        A = sorted({a for a, _ in pairs})
        B = sorted({b for _, b in pairs})
    else:
        A = [d for d in docs if lab[d] == "A"]
        B = [d for d in docs if lab[d] == "B"]
    return Match(node, ds.name, A, B, pairs, concepts,
                 hit.sum(0).astype(float), obs.sum(0).astype(float),
                 hit.sum(1).astype(float), obs.sum(1).astype(float), lab)

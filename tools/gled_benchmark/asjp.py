"""ASJP dating after Holman et al. (2011), with the published constants.

LDND = LDN / mean LDN of different-concept pairs, averaged over the 40
concepts; the similarity of two groups is 1 - mean LDND over the cross pairs;
t = (ln s - ln s0) / (2 ln r), s0 = 0.92, r = 0.72 per millennium. Years are
before 2000. Word strings lose the ASJP modifiers (* " ~ $) and keep the
first synonym, as in `data.Asjp`.
"""
import math
from functools import lru_cache
from itertools import combinations

S0 = 0.92
R = 0.72


@lru_cache(maxsize=None)
def lev(a, b):
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def ldn(a, b):
    return lev(a, b) / max(len(a), len(b))


def ldnd(wa, wb):
    """LDND of two languages given {concept: word}; None if < 10 shared concepts."""
    shared = [c for c in wa if c in wb]
    if len(shared) < 10:
        return None
    same = sum(ldn(wa[c], wb[c]) for c in shared) / len(shared)
    diff = [ldn(wa[i], wb[j]) for i in shared for j in shared if i != j]
    d = sum(diff) / len(diff)
    return same / d if d > 0 else None


def date_from_similarity(s):
    """Years before 2000 from a similarity s = 1 - LDND."""
    if s <= 0:
        return float("nan")
    s = min(s, 1.0)
    return (math.log(s) - math.log(S0)) / (2 * math.log(R)) * 1000.0


def node_date(words, side):
    """words: {code: {concept: str}}, side: {code: label}. Mean LDND over the
    pairs from different sides -> (date, n_pairs, mean LDND)."""
    ds = []
    for a, b in combinations(sorted(side), 2):
        if side[a] == side[b] or a not in words or b not in words:
            continue
        d = ldnd(words[a], words[b])
        if d is not None:
            ds.append(d)
    if not ds:
        return None
    m = sum(ds) / len(ds)
    return date_from_similarity(1.0 - m), len(ds), m

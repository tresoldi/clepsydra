"""D1.3: an empirical prior for list-level observation error (roadmap E2.9).

Two lists of the same language (same Glottocode) collected by different
people are separated by essentially no time, so the share of shared concepts
on which their cognate codes differ is observation error: with per-list
error probabilities e1, e2 it is about e1 + e2 - e1 e2. This script finds
every such pair in the expert-cognate Lexibank datasets, one dataset at a
time (to stay small in memory), and writes the mismatch rates.

    cd tools && python3 -m error_prior.run > ../docs/research/error-prior/pairs.tsv

Doculects flagged historical, proto or ancient are skipped; so are pairs
with fewer than 40 shared concepts, and datasets with degenerate coding
(fewer than 1.05 cognate sets per concept), as in the clock pilot.
"""
import csv
import itertools
import os
import sys

from clock_pilot import lexibank as L

MIN_SHARED = 40


def main():
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["dataset", "glottocode", "doculect_a", "doculect_b", "shared", "mismatch", "rate"])
    for name in L.expert_datasets():
        try:
            ds = L.load(name, None)
        except (OSError, KeyError) as exc:
            print(f"skip {name}: {exc}", file=sys.stderr)
            continue
        if ds.sets_per_concept < 1.05:
            continue
        by_g = {}
        for i, d in ds.doculects.items():
            if d.ancient or not d.glottocode:
                continue
            by_g.setdefault(d.glottocode, []).append(i)
        for g, ids in sorted(by_g.items()):
            for a, b in itertools.combinations(sorted(ids), 2):
                ta, tb = ds.table.get(a, {}), ds.table.get(b, {})
                shared = [c for c in ta if c in tb]
                if len(shared) < MIN_SHARED:
                    continue
                mism = sum(1 for c in shared if ta[c] != tb[c])
                out.writerow([name, g, ds.doculects[a].name, ds.doculects[b].name,
                              len(shared), mism, f"{mism / len(shared):.4f}"])
        del ds


if __name__ == "__main__":
    main()

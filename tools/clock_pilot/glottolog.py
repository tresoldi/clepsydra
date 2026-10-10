"""Glottolog classification, read straight from a Glottolog CLDF release zip.

Only two things are needed: each languoid's lineage (root first) and its level
(family / language / dialect). Both come from the CLDF release, `values.csv`
(parameter `classification`) and `languages.csv`.
"""

from __future__ import annotations

import csv
import hashlib
import io
import os
import zipfile

DEFAULT_ZIP = os.path.expanduser(
    "~/work/rgrambank/review/densify_GB/data/glottolog-cldf-v5.0.zip")
EXPECTED_SHA256 = (
    "e954bd77ac4e287f8c4588bd2216bdd662eeb86469b2c528f734319ba176a433")


class Glottolog:
    """Lineages and levels for every languoid in one Glottolog release."""

    def __init__(self, zip_path: str = DEFAULT_ZIP):
        self.path = zip_path
        with open(zip_path, "rb") as fh:
            self.sha256 = hashlib.sha256(fh.read()).hexdigest()
        self.name: dict[str, str] = {}
        self.level: dict[str, str] = {}
        self.lineage: dict[str, tuple[str, ...]] = {}
        with zipfile.ZipFile(zip_path) as zf:
            names = zf.namelist()
            lang = next(n for n in names if n.endswith("cldf/languages.csv"))
            vals = next(n for n in names if n.endswith("cldf/values.csv"))
            self.release = lang.split("/")[0]
            with zf.open(lang) as fh:
                for row in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
                    self.name[row["ID"]] = row["Name"]
                    self.level[row["ID"]] = row["Level"]
            path: dict[str, str] = {}
            with zf.open(vals) as fh:
                for row in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
                    if row["Parameter_ID"] == "classification":
                        path[row["Language_ID"]] = row["Value"]
        for gc in self.name:
            p = path.get(gc, "")
            anc = tuple(x for x in p.split("/") if x) if p else ()
            self.lineage[gc] = anc + (gc,)
        self.children: dict[str, list[str]] = {}
        for gc, lin in self.lineage.items():
            if len(lin) >= 2:
                self.children.setdefault(lin[-2], []).append(gc)

    def known(self, gc: str) -> bool:
        return gc in self.lineage

    def top(self, gc: str) -> str:
        return self.lineage[gc][0]

    def parent(self, gc: str) -> str | None:
        lin = self.lineage[gc]
        return lin[-2] if len(lin) >= 2 else None

    def is_under(self, gc: str, anc: str) -> bool:
        """True if `gc` is `anc` or one of its descendants."""
        return gc in self.lineage and anc in self.lineage[gc]

    def child_under(self, gc: str, anc: str) -> str | None:
        """The child of `anc` on the path to `gc` (None if gc == anc or not under)."""
        lin = self.lineage.get(gc, ())
        if anc not in lin:
            return None
        i = lin.index(anc)
        return lin[i + 1] if i + 1 < len(lin) else None

    def mrca(self, codes) -> str | None:
        lins = [self.lineage[c] for c in codes if c in self.lineage]
        if not lins:
            return None
        out = None
        for i in range(min(len(x) for x in lins)):
            s = {x[i] for x in lins}
            if len(s) != 1:
                break
            out = lins[0][i]
        return out

    def n_languages_under(self, gc: str) -> int:
        """Language-level descendants (dialects not counted; a language counts itself)."""
        n = 0
        stack = [gc]
        while stack:
            x = stack.pop()
            if self.level.get(x) == "language":
                n += 1
                continue  # dialects below a language are not languages
            stack.extend(self.children.get(x, []))
        return n

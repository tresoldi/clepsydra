"""Read one Lexibank CLDF dataset into (doculect, concept) -> cognate class.

Coding rule (pre-registered): per (doculect, concept) the first form in
`forms.csv` order that has a cognate judgement; per form, the first row of
`cognates.csv`. Forms without a cognate row are missing, not singletons.
Ancient / historical doculects are dropped (see `ANCIENT_RE`).
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
from dataclasses import dataclass, field

csv.field_size_limit(10 ** 9)

LEXIBANK = os.path.expanduser("~/work/chl/arcaverborum/raw/lexibank")
DATASETS_CSV = os.path.expanduser("~/work/chl/arcaverborum/datasets.csv")

ANCIENT_RE = re.compile(
    r"\b(old|middle|ancient|classical|archaic|proto|early|late|epigraphic|"
    r"medieval|vedic)\b|latin|gothic|hittite|tochar|tokhar|sanskrit|avestan|"
    r"akkadian|luvian|luwian|oscan|umbrian|gaulish|\bpali\b|prakrit|sogdian|"
    r"khotanese|parthian|bactrian|ge'ez|geez", re.I)
TRUE = {"true", "1", "yes", "t"}
# Modern doculects that the name rule would catch ("Middle" as a place name).
MODERN_BY_NAME = {"Middle Lozva Mansi", "Middle Chulym"}


def expert_datasets() -> list[str]:
    with open(DATASETS_CSV, encoding="utf-8") as fh:
        return [r["NAME"] for r in csv.DictReader(fh) if r["ExpertCognates"]]


def git_commit(name: str) -> str:
    try:
        out = subprocess.run(["git", "-C", os.path.join(LEXIBANK, name),
                              "rev-parse", "HEAD"], capture_output=True,
                             text=True, timeout=30)
        return out.stdout.strip() or "not-a-git-checkout"
    except OSError:
        return "unknown"


@dataclass
class Doculect:
    id: str
    name: str
    glottocode: str
    ancient: bool
    why: str = ""


@dataclass
class Dataset:
    name: str
    doculects: dict[str, Doculect]
    # doculect id -> concept id -> cognate class
    table: dict[str, dict[str, str]]
    concepticon: dict[str, str]           # concept id -> Concepticon ID
    forms_total: int = 0
    forms_with_cognate: int = 0
    commit: str = ""
    # (doculect, concept) -> (segments, form id) for the coded form
    segments: dict[tuple[str, str], tuple[str, str]] = field(default_factory=dict)

    @property
    def coverage(self) -> float:
        return self.forms_with_cognate / max(self.forms_total, 1)

    @property
    def sets_per_concept(self) -> float:
        """Mean number of distinct cognate classes per concept, over concepts
        with at least two coded doculects. About 1 means the coding assigns one
        set per meaning (e.g. etymological roots relative to a standard
        language), which carries no replacement signal."""
        per: dict[str, set] = {}
        cnt: dict[str, int] = {}
        for row in self.table.values():
            for par, cs in row.items():
                per.setdefault(par, set()).add(cs)
                cnt[par] = cnt.get(par, 0) + 1
        ks = [len(v) for p, v in per.items() if cnt[p] >= 2]
        return sum(ks) / len(ks) if ks else 0.0


def load(name: str, glotto_names: dict[str, str] | None = None) -> Dataset:
    base = os.path.join(LEXIBANK, name, "cldf")
    docs: dict[str, Doculect] = {}
    with open(os.path.join(base, "languages.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            gc = (r.get("Glottocode") or "").strip()
            why = ""
            if (r.get("historical") or "").strip().lower() in TRUE:
                why = "historical flag"
            elif (r.get("Is_Proto") or "").strip().lower() in TRUE:
                why = "Is_Proto flag"
            elif (r.get("Name") or "") in MODERN_BY_NAME:
                why = ""
            elif ANCIENT_RE.search(r.get("Name") or ""):
                why = "name"
            elif glotto_names and gc and ANCIENT_RE.search(glotto_names.get(gc, "")):
                why = "Glottolog name"
            docs[r["ID"]] = Doculect(r["ID"], r.get("Name", ""), gc, bool(why), why)
    conc: dict[str, str] = {}
    with open(os.path.join(base, "parameters.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            conc[r["ID"]] = (r.get("Concepticon_ID") or "").strip()
    cog: dict[str, str] = {}
    cpath = os.path.join(base, "cognates.csv")
    if os.path.exists(cpath):
        with open(cpath, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                fid = r["Form_ID"]
                if fid not in cog and (r.get("Cognateset_ID") or "").strip():
                    cog[fid] = r["Cognateset_ID"].strip()
    table: dict[str, dict[str, str]] = {}
    segs: dict[tuple[str, str], tuple[str, str]] = {}
    total = with_cog = 0
    with open(os.path.join(base, "forms.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            total += 1
            cs = cog.get(r["ID"])
            if cs is None:
                continue
            with_cog += 1
            lang, par = r["Language_ID"], r["Parameter_ID"]
            row = table.setdefault(lang, {})
            if par not in row:
                row[par] = cs
                segs[(lang, par)] = ((r.get("Segments") or "").strip(), r["ID"])
    return Dataset(name, docs, table, conc, total, with_cog, git_commit(name), segs)

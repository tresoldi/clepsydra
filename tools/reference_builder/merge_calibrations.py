"""Merge the two regional calibration tables into one common schema.

The regional tables (reference/calibrations/oldworld.tsv and
newworld_pacific.tsv) were compiled separately and use different column
names. This writes reference/calibrations/calibrations.tsv with one schema,
keeping every value verbatim; nothing is re-graded or re-dated here.

    python3 -m reference_builder.merge_calibrations   (from tools/)
"""

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2] / "reference" / "calibrations"

COLUMNS = [
    "id", "region", "family", "target_kind", "target", "node_definition",
    "glottocodes", "event", "age_min_bp", "age_max_bp", "age_point_bp",
    "bound_type", "suggested_distribution", "hard_soft", "evidence_type",
    "grade", "notes", "source_original", "source_dating_study",
]

OLD = {
    "id": "id", "family": "family", "target": "clade_or_language",
    "node_definition": "node_definition", "glottocodes": "glottocodes",
    "event": "event", "age_min_bp": "age_min_BP", "age_max_bp": "age_max_BP",
    "age_point_bp": "age_point_BP", "bound_type": "bound_type",
    "suggested_distribution": "suggested_distribution",
    "hard_soft": "hard_soft", "evidence_type": "evidence_type",
    "grade": "grade", "source_original": "source_original",
    "source_dating_study": "source_dating_study",
}

NEW = {
    "id": "id", "family": "family", "target": "target",
    "node_definition": "node_definition", "glottocodes": "glottocodes",
    "event": "event", "age_min_bp": "age_BP_min", "age_max_bp": "age_BP_max",
    "age_point_bp": "age_BP_point", "bound_type": "bound",
    "suggested_distribution": "suggested_distribution",
    "hard_soft": "hard_soft", "evidence_type": "evidence_type",
    "grade": "grade", "notes": "controversy_notes",
}


def kind_old(raw):
    raw = raw.lower()
    for k in ("stem", "tip", "root"):
        if k in raw:
            return k
    return "node"


def kind_new(raw):
    return {"node_split": "node", "tip_attestation": "tip",
            "settlement_bound": "settlement_bound"}.get(raw, raw)


def clean(v):
    v = (v or "").replace("\t", " ").replace("\n", " ").strip()
    return v if v else "."


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    out = []
    for row in read(ROOT / "oldworld.tsv"):
        rec = {c: clean(row.get(src)) for c, src in OLD.items()}
        rec["region"] = "oldworld"
        rec["target_kind"] = kind_old(row["calib_type"])
        notes = [row.get("original_spec_as_published"),
                 row.get("alt_values_other_studies"), row.get("notes")]
        rec["notes"] = clean(" | ".join(n for n in notes if n and n != "."))
        out.append(rec)
    for row in read(ROOT / "newworld_pacific.tsv"):
        rec = {c: clean(row.get(src)) for c, src in NEW.items()}
        rec["region"] = "newworld_pacific"
        rec["target_kind"] = kind_new(row["type"])
        src = [row.get("original_source"), row.get("original_doi")]
        rec["source_original"] = clean(" ".join(s for s in src if s and s != "."))
        dat = [row.get("dating_study"), row.get("dating_doi")]
        rec["source_dating_study"] = clean(" ".join(s for s in dat if s and s != "."))
        out.append(rec)
    ids = [r["id"] for r in out]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        sys.exit(f"duplicate calibration ids: {sorted(dup)}")
    for r in out:
        for c in ("age_min_bp", "age_max_bp", "age_point_bp"):
            if r[c] != "." and not re.fullmatch(r"-?\d+(\.\d+)?", r[c]):
                sys.exit(f"{r['id']}: non-numeric {c} {r[c]!r}")
    with open(ROOT / "calibrations.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows({c: r.get(c, ".") for c in COLUMNS} for r in out)
    print(f"wrote {len(out)} rows")


if __name__ == "__main__":
    main()

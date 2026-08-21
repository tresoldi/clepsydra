"""
Plain cognate-matrix -> Clepsydra TSV (M12).

For cognate data that is *not* in CLDF — the common case for a phylogenetics
project that keeps a hand-built cognate matrix. Two layouts are read:

  wide:  a table with one row per language and one column per meaning; each cell
         holds the cognate-class label(s) for that (language, meaning). The first
         column is the language id. Example (comma- or tab-separated):

             language,water,fire,hand
             Hmong,1,3,2
             Mien,1,4,2
             She,5,3,2

  long:  a tidy table with a language column, a concept column and a cognate
         column (one row per observed form), like a stripped-down forms table.

Cells that are empty or a missing marker (`?`, `-`, `.`, `NA` by default) become
no form (missing). Several cognate labels in one wide cell (synonyms / partial
cognates) may be split on a separator, each kept as a weighted candidate form.
Cognate labels only need equality *within* a meaning, so they are namespaced by
concept automatically.

This is the T1/T2 cognate-class path. A Nexus/BEAST *binary* presence-absence
matrix is a different encoding (one column per cognate set); it is not read here.

Stdlib only.

Usage:
  python -m datasets.matrix --format wide --input cognates.csv --out DIR [--family hmong_mien]
  python -m datasets.matrix --format long --input forms.csv --out DIR \\
      --lang-col Language --concept-col Concept --cognate-col Cognateset
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from datasets import lexibank   # reuse the TSV writers/headers

DEFAULT_MISSING = {"", "?", "-", ".", "na", "n/a", "nan", "none", "null"}


def _sniff_delimiter(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        first = f.readline()
    if "\t" in first:
        return "\t"
    if ";" in first and "," not in first:
        return ";"
    return ","


def _clean_id(s):
    return "".join(c if (c.isalnum() or c in "_-") else "_" for c in s.strip())


def _is_missing(v, missing):
    return v is None or v.strip().lower() in missing


def read_wide(path, delimiter=None, missing=None, synonym_sep=None, family="."):
    missing = missing or DEFAULT_MISSING
    delimiter = delimiter or _sniff_delimiter(path)
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f, delimiter=delimiter))
    if not rows:
        raise ValueError("empty matrix: %s" % path)
    header = rows[0]
    concept_labels = header[1:]
    concept_ids = ["C%s" % _clean_id(c) for c in concept_labels]

    lang_rows, form_rows = [], []
    seen_lang = set()
    for r in rows[1:]:
        if not r or _is_missing(r[0], missing):
            continue
        lang = _clean_id(r[0])
        if lang in seen_lang:
            raise ValueError("duplicate language row: %s" % r[0])
        seen_lang.add(lang)
        lang_rows.append([lang, r[0].strip(), family, "0", "."])
        for j, cid in enumerate(concept_ids):
            cell = r[j + 1] if j + 1 < len(r) else ""
            if _is_missing(cell, missing):
                continue
            vals = [v for v in (cell.split(synonym_sep) if synonym_sep else [cell])
                    if not _is_missing(v, missing)]
            if not vals:
                continue
            n = len(vals)
            for k, v in enumerate(vals):
                form_rows.append([
                    "matrix", lang, cid, "%s-%s-%d" % (lang, cid, k), ".",
                    "%s-%s" % (concept_ids[j], _clean_id(v)),
                    "%.4f" % (1.0 / n), "1.0", "0.95", "0.02", "observed",
                ])
    concept_rows = [[cid, lbl, ".", ".", "unknown", "active"]
                    for cid, lbl in zip(concept_ids, concept_labels)]
    form_rows.sort(key=lambda x: (x[1], x[2], x[3]))
    return lang_rows, concept_rows, form_rows


def read_long(path, lang_col, concept_col, cognate_col, weight_col=None,
              delimiter=None, missing=None, family="."):
    missing = missing or DEFAULT_MISSING
    delimiter = delimiter or _sniff_delimiter(path)
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = list(reader)
    langs, concepts = {}, {}
    form_rows = []
    for i, r in enumerate(rows):
        lang_raw, con_raw = r.get(lang_col), r.get(concept_col)
        cog = r.get(cognate_col)
        if _is_missing(lang_raw, missing) or _is_missing(con_raw, missing) \
                or _is_missing(cog, missing):
            continue
        lang = _clean_id(lang_raw)
        cid = "C%s" % _clean_id(con_raw)
        langs.setdefault(lang, lang_raw.strip())
        concepts.setdefault(cid, con_raw.strip())
        w = "1.0"
        if weight_col and not _is_missing(r.get(weight_col), missing):
            w = r[weight_col]
        form_rows.append([
            "matrix", lang, cid, "%s-%s-%d" % (lang, cid, i), ".",
            "%s-%s" % (cid, _clean_id(cog)), w, "1.0", "0.95", "0.02", "observed",
        ])
    lang_rows = [[lid, name, family, "0", "."] for lid, name in langs.items()]
    concept_rows = [[cid, lbl, ".", ".", "unknown", "active"]
                    for cid, lbl in concepts.items()]
    form_rows.sort(key=lambda x: (x[1], x[2], x[3]))
    return lang_rows, concept_rows, form_rows


def main(argv):
    ap = argparse.ArgumentParser(prog="datasets.matrix")
    ap.add_argument("--format", choices=["wide", "long"], required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--delimiter", default=None, help="default: sniffed")
    ap.add_argument("--family", default=".", help="family_id for all languages")
    ap.add_argument("--synonym-sep", default=None,
                    help="split multiple cognate labels in a wide cell, e.g. ';'")
    ap.add_argument("--missing", nargs="*", default=None,
                    help="extra missing markers (case-insensitive)")
    ap.add_argument("--lang-col", default="Language")
    ap.add_argument("--concept-col", default="Concept")
    ap.add_argument("--cognate-col", default="Cognateset")
    ap.add_argument("--weight-col", default=None)
    args = ap.parse_args(argv[1:])
    missing = set(DEFAULT_MISSING)
    if args.missing:
        missing |= {m.lower() for m in args.missing}
    if args.format == "wide":
        lr, cr, fr = read_wide(args.input, args.delimiter, missing,
                               args.synonym_sep, args.family)
    else:
        lr, cr, fr = read_long(args.input, args.lang_col, args.concept_col,
                               args.cognate_col, args.weight_col,
                               args.delimiter, missing, args.family)
    lexibank.write_subset(args.out, lr, cr, fr)
    print("Wrote %d languages, %d concepts, %d forms to %s"
          % (len(lr), len(cr), len(fr), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

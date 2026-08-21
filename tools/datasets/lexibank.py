"""
Generic Lexibank CLDF -> Clepsydra TSV (M12).

Converts any Lexibank/CLDF Wordlist with a `forms.csv` and a `cognates.csv`
(Form_ID -> Cognateset_ID) into Clepsydra's strict-TSV input schema
(docs/input-formats.md), restricted to a chosen set of doculects. This is the
family-agnostic converter: it drives the Indo-European (IE-CoR), Sino-Tibetan
(Sagart et al. 2019), Austronesian (ABVD), Semitic, Bantu, etc. datasets of
Lexibank's LexiCore collection uniformly. Datasets are not redistributed here;
fetch a CLDF release first, then run this converter (see examples/DATA_SOURCE.md).

LexiCore is Lexibank's core collection of expert-cognate-coded wordlists
(Blum, Barrientos, Forkel, List et al.; https://github.com/lexibank). Each
source dataset carries its own citation and (CC-BY) licence in its
`cldf-metadata.json`; record it when you publish results.

Stdlib only.

Usage:
  python -m datasets.lexibank --cldf /path/to/<dataset>/cldf \\
      --languages Sinitic_Old_Chinese Sinitic_Beijing Tibetan_Lhasa \\
      --out examples/sinotibetan_tree/data
"""
import argparse
import csv
import os
import sys

csv.field_size_limit(64 * 1024 * 1024)


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _get(row, *names, default=""):
    """First present, non-empty value among column `names`."""
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return default


def load_cldf(cldf_dir):
    langs = _read_csv(os.path.join(cldf_dir, "languages.csv"))
    params = _read_csv(os.path.join(cldf_dir, "parameters.csv"))
    forms = _read_csv(os.path.join(cldf_dir, "forms.csv"))
    form_to_cog, form_doubt = {}, {}
    cog_path = os.path.join(cldf_dir, "cognates.csv")
    if os.path.exists(cog_path):
        for c in _read_csv(cog_path):
            form_to_cog[c["Form_ID"]] = c["Cognateset_ID"]
            form_doubt[c["Form_ID"]] = (_get(c, "Doubt").lower() == "true")
    else:
        # Some datasets carry cognacy inline in forms.csv (Cognacy column).
        for f in forms:
            cog = _get(f, "Cognacy")
            if cog:
                form_to_cog[f["ID"]] = cog.split()[0]  # first partial-cognate id
    return langs, params, forms, form_to_cog, form_doubt


def build_subset(cldf_dir, language_names):
    """Return (languages, concepts, forms) rows for the Clepsydra schema.

    `language_names` matches a doculect's `Name`, `ascii_name`, or CLDF `ID`.
    """
    langs, params, forms, form_to_cog, form_doubt = load_cldf(cldf_dir)

    def lang_key(l):
        return _get(l, "ascii_name", "Name", "ID")

    index = {}
    for l in langs:
        for key in {_get(l, "ID"), _get(l, "Name"), _get(l, "ascii_name")}:
            if key:
                index[key] = l
    missing = [n for n in language_names if n not in index]
    if missing:
        raise ValueError("languages not in dataset: %s" % ", ".join(missing))
    sel = [index[n] for n in language_names]
    # A stable, ASCII-ish per-language id used as the tree label.
    lang_id = {l["ID"]: (language_names[i]) for i, l in enumerate(sel)}
    sel_ids = set(lang_id)

    param_gloss = {p["ID"]: _get(p, "Concepticon_Gloss", "Name", default=p["ID"])
                   for p in params}
    param_cnc = {p["ID"]: _get(p, "Concepticon_ID", default=".") for p in params}

    sel_forms = [f for f in forms if f["Language_ID"] in sel_ids
                 and f["ID"] in form_to_cog]
    used_concepts = sorted({f["Parameter_ID"] for f in sel_forms})

    lang_rows = []
    for l in sel:
        extinct = "1" if _get(l, "historical", "Extinct").lower() == "true" else "0"
        lang_rows.append([lang_id[l["ID"]], _get(l, "Name", default=l["ID"]),
                          _get(l, "Family", default="."), extinct,
                          _get(l, "Glottocode", default=".")])

    concept_rows = []
    for pid in used_concepts:
        cnc = param_cnc.get(pid) or "."
        ns = "Concepticon" if cnc not in (".", "") else "."
        concept_rows.append(["C%s" % pid, param_gloss.get(pid, pid), ns,
                             cnc or ".", "unknown", "active"])

    cells = {}
    for f in sel_forms:
        cells.setdefault((f["Language_ID"], f["Parameter_ID"]), []).append(f)
    form_rows = []
    for (lid, pid), cell in cells.items():
        n = len(cell)
        for f in cell:
            cog = form_to_cog[f["ID"]]
            doubt = form_doubt.get(f["ID"], False)
            loan = (_get(f, "Loan").lower() == "true")
            form_rows.append([
                "lexibank", lang_id[lid], "C%s" % pid, f["ID"],
                _get(f, "Form", "Value", default="."),
                "%s-%s" % (pid, cog),
                "%.4f" % (1.0 / n),
                "0.5" if doubt else "1.0",
                "0.2" if loan else "0.95",
                "0.8" if loan else "0.02",
                "observed",
            ])
    form_rows.sort(key=lambda r: (r[1], r[2], r[3]))
    return lang_rows, concept_rows, form_rows


LANG_HEADER = ["language_id", "name", "family_id", "is_extinct", "metadata"]
CONCEPT_HEADER = ["concept_id", "label", "external_namespace", "external_id",
                  "word_class", "status"]
FORM_HEADER = ["analysis_id", "language_id", "concept_id", "form_id", "form",
               "cognate_id", "primary_weight", "cognacy_confidence",
               "inheritance_probability", "loan_probability", "status"]


def _write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def write_subset(out_dir, lang_rows, concept_rows, form_rows):
    os.makedirs(out_dir, exist_ok=True)
    _write_tsv(os.path.join(out_dir, "languages.tsv"), LANG_HEADER, lang_rows)
    _write_tsv(os.path.join(out_dir, "concepts.tsv"), CONCEPT_HEADER, concept_rows)
    _write_tsv(os.path.join(out_dir, "forms.tsv"), FORM_HEADER, form_rows)


def main(argv):
    ap = argparse.ArgumentParser(prog="datasets.lexibank")
    ap.add_argument("--cldf", required=True, help="CLDF directory")
    ap.add_argument("--languages", nargs="+", required=True,
                    help="doculects by Name / ascii_name / CLDF ID")
    ap.add_argument("--out", required=True, help="output data directory")
    args = ap.parse_args(argv[1:])
    lr, cr, fr = build_subset(args.cldf, args.languages)
    write_subset(args.out, lr, cr, fr)
    print("Wrote %d languages, %d concepts, %d forms to %s"
          % (len(lr), len(cr), len(fr), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

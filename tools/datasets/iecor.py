"""
IE-CoR (Indo-European Cognate Relationships) -> Clepsydra TSV (M12).

Converts the CLDF release of IE-CoR (Heggarty, Anderson & Scarborough; lexibank
`iecor` v1.0, CC-BY-4.0) into Clepsydra's strict-TSV input schema
(docs/input-formats.md), restricted to a chosen set of doculects. IE-CoR is not
redistributed here; fetch it first (see `fetch.sh` in an example directory or
`--cldf` below), then run this converter to regenerate the committed subsets.

Attribution (CC-BY-4.0): Heggarty, Paul, Cormac Anderson & Matthew Scarborough.
Indo-European Cognate Relationships database (IE-CoR), version 1.0, 2019.
https://github.com/lexibank/iecor

Stdlib only.

Usage:
  python -m datasets.iecor --cldf /path/to/iecor/cldf \\
      --languages Latin Italian Spanish Portuguese French Romanian \\
      --out examples/romance_tree/data
"""
import argparse
import csv
import os
import sys

csv.field_size_limit(10 * 1024 * 1024)


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_cldf(cldf_dir):
    langs = _read_csv(os.path.join(cldf_dir, "languages.csv"))
    params = _read_csv(os.path.join(cldf_dir, "parameters.csv"))
    forms = _read_csv(os.path.join(cldf_dir, "forms.csv"))
    cognates = _read_csv(os.path.join(cldf_dir, "cognates.csv"))
    form_to_cog = {}
    form_doubt = {}
    for c in cognates:
        form_to_cog[c["Form_ID"]] = c["Cognateset_ID"]
        form_doubt[c["Form_ID"]] = (c.get("Doubt", "").lower() == "true")
    return langs, params, forms, form_to_cog, form_doubt


def build_subset(cldf_dir, language_names, concept_ids=None):
    """Return (languages, concepts, forms) rows for the Clepsydra schema.

    `language_names` matches the IE-CoR `ascii_name`; `concept_ids` (IE-CoR
    parameter IDs) defaults to all concepts present for the chosen languages.
    """
    langs, params, forms, form_to_cog, form_doubt = load_cldf(cldf_dir)

    by_ascii = {l["ascii_name"]: l for l in langs}
    missing = [n for n in language_names if n not in by_ascii]
    if missing:
        raise ValueError("languages not in IE-CoR: %s" % ", ".join(missing))
    sel_langs = [by_ascii[n] for n in language_names]
    sel_lang_ids = {l["ID"] for l in sel_langs}
    param_gloss = {p["ID"]: (p["Concepticon_Gloss"] or p["Name"]) for p in params}
    param_cnc = {p["ID"]: p["Concepticon_ID"] for p in params}

    # Forms restricted to the selected languages.
    sel_forms = [f for f in forms if f["Language_ID"] in sel_lang_ids]
    if concept_ids is not None:
        cset = set(concept_ids)
        sel_forms = [f for f in sel_forms if f["Parameter_ID"] in cset]
    used_concepts = sorted({f["Parameter_ID"] for f in sel_forms}, key=int)

    # languages.tsv rows.
    lang_rows = []
    for l in sel_langs:
        extinct = "1" if l["historical"].lower() == "true" else "0"
        lang_rows.append([l["ascii_name"], l["Name"], "indo_european",
                          extinct, l["Glottocode"] or "."])

    # concepts.tsv rows.
    concept_rows = []
    for pid in used_concepts:
        cnc = param_cnc.get(pid) or "."
        ns = "Concepticon" if cnc != "." else "."
        concept_rows.append(["C%s" % pid, param_gloss[pid], ns, cnc,
                             "unknown", "active"])

    # forms.tsv rows, grouping synonyms per (language, concept) to weight them.
    cells = {}
    for f in sel_forms:
        cells.setdefault((f["Language_ID"], f["Parameter_ID"]), []).append(f)
    form_rows = []
    for (lid, pid), cell in cells.items():
        lang_name = next(l["ascii_name"] for l in sel_langs if l["ID"] == lid)
        n = len(cell)
        for f in cell:
            cog = form_to_cog.get(f["ID"])
            if cog is None:
                continue  # form without a cognate judgement: skip
            doubt = form_doubt.get(f["ID"], False)
            loan = (f.get("Loan", "").lower() == "true")
            form_rows.append([
                "iecor",                       # analysis_id
                lang_name,                      # language_id (= tree label)
                "C%s" % pid,                    # concept_id
                f["ID"],                        # form_id
                f["Form"] or ".",               # form
                "%s-%s" % (pid, cog),           # cognate_id (namespaced by concept)
                "%.4f" % (1.0 / n),             # primary_weight (split over synonyms)
                "0.5" if doubt else "1.0",      # cognacy_confidence (y_c)
                "0.2" if loan else "0.95",      # inheritance_probability
                "0.8" if loan else "0.02",      # loan_probability
                "observed",                     # status
            ])
    # Stable ordering for reproducible output.
    form_rows.sort(key=lambda r: (r[1], int(r[2][1:]), r[3]))
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
    ap = argparse.ArgumentParser(prog="datasets.iecor")
    ap.add_argument("--cldf", required=True, help="IE-CoR CLDF directory")
    ap.add_argument("--languages", nargs="+", required=True,
                    help="IE-CoR ascii_name doculects to include")
    ap.add_argument("--out", required=True, help="output data directory")
    args = ap.parse_args(argv[1:])
    lang_rows, concept_rows, form_rows = build_subset(args.cldf, args.languages)
    write_subset(args.out, lang_rows, concept_rows, form_rows)
    print("Wrote %d languages, %d concepts, %d forms to %s"
          % (len(lang_rows), len(concept_rows), len(form_rows), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

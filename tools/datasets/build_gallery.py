"""
Cross-family pairwise gallery builder (M12).

For each family dataset, auto-selects an informative language pair (well-attested
and moderately diverged) and emits that pair's cognate data as Clepsydra TSV, so
`clepsydra date-pair` can be run across many families. Pairwise (P1) dating needs
no tree or calibration, so it scales honestly to any family in Lexibank's
LexiCore collection.

The family datasets are not redistributed; fetch each CLDF release first (see the
`sources` below and examples/DATA_SOURCE.md), then run this builder.

Stdlib only.

Usage:
  python -m datasets.build_gallery --root /path/with/<dataset>/cldf-dirs \\
      --out ../examples/families_pairwise
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from datasets import lexibank


# family label -> (lexibank dataset id, short citation for attribution)
FAMILIES = [
    ("dravidian",   "dravlex",            "Kolipakam et al. 2018 (DravLEX)"),
    ("turkic",      "savelyevturkic",     "Savelyev & List 2020 (Turkic)"),
    ("uralic",      "uralex",             "Syrjänen et al. 2018 (UraLex)"),
    ("semitic",     "felekesemitic",      "Feleke 2021 (Semitic)"),
    ("polynesian",  "walworthpolynesian", "Walworth 2018 (Polynesian)"),
    ("tupian",      "tuled",              "Gerardi et al. 2021 (TuLeD)"),
    ("arawakan",    "chaconarawakan",     "Chacon et al. (Arawakan)"),
    ("alorpantar",  "robinsonap",         "Robinson & Holton 2012 (Alor-Pantar)"),
    ("koreanic",    "leekoreanic",        "Lee 2015 (Koreanic)"),
]


def language_cognates(cldf_dir):
    """Map language ID -> {concept: primary cognate id}."""
    langs, params, forms, form_to_cog, _ = lexibank.load_cldf(cldf_dir)
    name = {l["ID"]: lexibank._get(l, "Name", "ID") for l in langs}
    per = {}
    for f in forms:
        if f["ID"] not in form_to_cog:
            continue
        lid, pid = f["Language_ID"], f["Parameter_ID"]
        cell = per.setdefault(lid, {})
        cell.setdefault(pid, form_to_cog[f["ID"]])  # first form = primary
    return name, per


def best_pair(per, target_retention=0.55, min_shared=40):
    """Pick the (a,b) language pair with the most informative divergence: enough
    shared concepts, retention nearest a moderate target. Deterministic."""
    ids = sorted(per)
    best = None
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = per[ids[i]], per[ids[j]]
            shared = [c for c in a if c in b]
            if len(shared) < min_shared:
                continue
            match = sum(1 for c in shared if a[c] == b[c])
            ret = match / len(shared)
            # favour many shared concepts and retention near the target
            score = len(shared) * (1.0 - abs(ret - target_retention))
            if best is None or score > best[0]:
                best = (score, ids[i], ids[j], len(shared), ret)
    if best is None:  # relax the shared-concept floor
        return best_pair(per, target_retention, min_shared=max(10, min_shared // 2)) \
            if min_shared > 10 else None
    return best


def build(root, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    summary = []
    for family, ds, cite in FAMILIES:
        cldf = os.path.join(root, ds, "")
        cldf = os.path.join(root, ds)
        if not os.path.isdir(cldf):
            # allow either <root>/<ds> or <root>/<ds>/cldf
            alt = os.path.join(root, ds, "cldf")
            cldf = alt if os.path.isdir(alt) else cldf
        if not os.path.isfile(os.path.join(cldf, "forms.csv")):
            print("  skip %s: no CLDF at %s" % (family, cldf))
            continue
        name, per = language_cognates(cldf)
        bp = best_pair(per)
        if bp is None:
            print("  skip %s: no suitable pair" % family)
            continue
        _, aid, bid, shared, ret = bp
        # Use the ASCII-safe CLDF IDs as language_ids (CLI-friendly labels).
        lr, cr, fr = lexibank.build_subset(cldf, [aid, bid])
        fam_dir = os.path.join(out_dir, family)
        lexibank.write_subset(fam_dir, lr, cr, fr)
        summary.append((family, aid, bid, name[aid], name[bid], shared, ret, cite))
        print("  %-11s %s / %s  (shared=%d retention=%.2f)"
              % (family, aid, bid, shared, ret))
    # manifest for the run script + README
    with open(os.path.join(out_dir, "pairs.tsv"), "w", encoding="utf-8") as f:
        f.write("family\tlang_a_id\tlang_b_id\tlang_a_name\tlang_b_name\t"
                "shared_concepts\tretention\tsource\n")
        for fam, aid, bid, an, bn, s, r, cite in summary:
            f.write("%s\t%s\t%s\t%s\t%s\t%d\t%.4f\t%s\n"
                    % (fam, aid, bid, an, bn, s, r, cite))
    return summary


def main(argv):
    ap = argparse.ArgumentParser(prog="datasets.build_gallery")
    ap.add_argument("--root", required=True, help="dir containing <dataset>/ CLDF folders")
    ap.add_argument("--out", required=True, help="examples/families_pairwise")
    args = ap.parse_args(argv[1:])
    build(args.root, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

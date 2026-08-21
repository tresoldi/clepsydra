"""
Reference-prior builder CLI.

Builds a versioned reference prior package from (currently synthetic) evidence
sources and emits strict TSV + manifest with checksums. Supports leave-family-out
packages for cross-family validation.

Usage:
  python -m reference_builder build --out reference/generated/synthetic_v1
  python -m reference_builder build --out ... --leave-family-out indo_european
  python -m reference_builder validate reference/generated/synthetic_v1
"""

import argparse
import os
import sys

from . import synthetic_sources as synth
from . import emitter
from .concepts import get_concepts


PACKAGE_VERSION = "1.0.0"


def cmd_build(args):
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    # Gather evidence (synthetic; replace with real parsers when supplied).
    pagel = synth.generate_pagel_rates(seed=args.seed)
    vejdemo = synth.generate_vejdemo_predictors(seed=args.seed + 1)
    asjp = synth.generate_asjp_stability(seed=args.seed + 2)
    families = synth.generate_family_priors()

    # Clock scale: use the legacy Swadesh-100 anchor as the global scale.
    clock_scale_log_mean = synth.LEGACY_SWADESH100["log_mean"]
    clock_scale_log_sd = synth.LEGACY_SWADESH100["log_sd"]

    warnings = [
        "SYNTHETIC DATA: Pagel 2007 Supplementary Table 2 not supplied; "
        "using synthetic rates. Replace with real data when available.",
        "SYNTHETIC DATA: Vejdemo-Hörberg 2016 Supporting Info not supplied; "
        "using synthetic predictors.",
        "SYNTHETIC DATA: ASJP database not supplied; using synthetic stability.",
    ]

    # Build concept priors.
    concept_rows, concept_header = emitter.build_concept_priors(
        pagel, vejdemo, asjp, clock_scale_log_mean, clock_scale_log_sd)
    concept_path = os.path.join(out_dir, "concept_priors.tsv")
    emitter.write_tsv(concept_path, concept_header, concept_rows)

    # Build family priors.
    family_rows, family_header = emitter.build_family_priors(families)
    family_path = os.path.join(out_dir, "family_priors.tsv")
    emitter.write_tsv(family_path, family_header, family_rows)

    # Leave-family-out: if requested, emit a separate package excluding the
    # named family's specific prior (keeping GLOBAL + other families).
    if args.leave_family_out:
        lf_families = [f for f in families
                       if f["family_id"] != args.leave_family_out]
        lf_rows, lf_header = emitter.build_family_priors(lf_families)
        lf_dir = os.path.join(out_dir, f"leave_out_{args.leave_family_out}")
        os.makedirs(lf_dir, exist_ok=True)
        emitter.write_tsv(os.path.join(lf_dir, "concept_priors.tsv"),
                          concept_header, concept_rows)
        emitter.write_tsv(os.path.join(lf_dir, "family_priors.tsv"),
                          lf_header, lf_rows)
        emitter.write_manifest(lf_dir,
                               ["concept_priors.tsv", "family_priors.tsv"],
                               {
                                   "version": PACKAGE_VERSION,
                                   "sources": ["synthetic"],
                                   "warnings": warnings + [
                                       f"Leave-family-out: {args.leave_family_out} "
                                       f"family-specific prior excluded."],
                               })
        print(f"Leave-family-out package: {lf_dir}")

    # Emit concepts.tsv for reference (the canonical mapping).
    concepts = get_concepts()
    concept_header_full = [
        "concept_id", "label", "external_namespace", "external_id",
        "word_class", "status",
    ]
    concept_full_rows = [
        [c["concept_id"], c["label"], c["external_namespace"],
         c["external_id"], c["word_class"], c["status"]]
        for c in concepts
    ]
    emitter.write_tsv(os.path.join(out_dir, "concepts.tsv"),
                      concept_header_full, concept_full_rows)

    # Manifest with checksums.
    manifest = emitter.write_manifest(
        out_dir,
        ["concept_priors.tsv", "family_priors.tsv", "concepts.tsv"],
        {
            "version": PACKAGE_VERSION,
            "sources": ["PAGEL_2007_SYNTH", "VEJDEMO_2016_SYNTH",
                         "ASJP_SYNTH"],
            "warnings": warnings,
        })

    print(f"Built reference package: {out_dir}")
    print(f"  concept_priors.tsv: {len(concept_rows)} rows")
    print(f"  family_priors.tsv: {len(family_rows)} rows")
    print(f"  concepts.tsv: {len(concepts)} rows")
    print(f"  manifest.json: {len(manifest['files'])} files with checksums")
    if warnings:
        print("  Warnings:")
        for w in warnings:
            print(f"    - {w}")
    return 0


def cmd_validate(args):
    """Validate a reference package: check TSV structure, headers, no
    duplicate IDs, manifest checksums match."""
    pkg_dir = args.package
    errors = []

    # Check required files exist.
    for fname in ["concept_priors.tsv", "family_priors.tsv", "manifest.json"]:
        if not os.path.exists(os.path.join(pkg_dir, fname)):
            errors.append(f"Missing file: {fname}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}")
        return 1

    # Validate TSV headers.
    import csv
    for fname, expected_header in [
        ("concept_priors.tsv",
         ["concept_id", "scope_id", "unit", "log_rate_mean", "log_rate_sd",
          "n_families", "n_observations", "evidence_grade", "source_key",
          "version"]),
        ("family_priors.tsv",
         ["family_id", "log_clock_mean", "log_clock_sd", "evidence_grade",
          "source_key", "version"]),
    ]:
        path = os.path.join(pkg_dir, fname)
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header != expected_header:
                errors.append(f"{fname}: bad header {header}")

    # Check for duplicate concept_ids in concept_priors.
    path = os.path.join(pkg_dir, "concept_priors.tsv")
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip header
        seen = set()
        for row in reader:
            if row and row[0]:
                if row[0] in seen:
                    errors.append(f"Duplicate concept_id: {row[0]}")
                seen.add(row[0])

    # Validate manifest checksums.
    import json
    with open(os.path.join(pkg_dir, "manifest.json"), "r") as f:
        manifest = json.load(f)
    for fname, info in manifest.get("files", {}).items():
        fpath = os.path.join(pkg_dir, fname)
        actual = emitter.sha256_file(fpath)
        if actual != info["sha256"]:
            errors.append(f"Checksum mismatch: {fname}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}")
        return 1

    print(f"Package {pkg_dir}: valid")
    print(f"  {len(seen)} unique concept priors")
    return 0


def main():
    parser = argparse.ArgumentParser(prog="reference_builder")
    sub = parser.add_subparsers(dest="command")

    p_build = sub.add_parser("build", help="Build a reference prior package")
    p_build.add_argument("--out", required=True, help="Output directory")
    p_build.add_argument("--seed", type=int, default=42, help="Random seed")
    p_build.add_argument("--leave-family-out", default=None,
                         help="Also emit a leave-family-out package")
    p_build.set_defaults(func=cmd_build)

    p_val = sub.add_parser("validate", help="Validate a reference package")
    p_val.add_argument("package", help="Package directory")
    p_val.set_defaults(func=cmd_validate)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
"""
TSV emission and manifest generation for the reference-prior builder.

Emits strict TSV files (concept_priors.tsv, family_priors.tsv) per
docs/input-formats.md, plus a manifest with checksums and provenance.
Python stdlib only.
"""

import csv
import hashlib
import json
import os
from datetime import datetime, timezone


def write_tsv(path, header, rows):
    """Write a strict TSV file with a header row and data rows.
    rows is a list of lists (or tuples); each inner list matches the header.
    Missing values are written as '.'.
    """
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            cells = []
            for val in row:
                if val is None or val == "":
                    cells.append(".")
                else:
                    cells.append(str(val))
            f.write("\t".join(cells) + "\n")


def sha256_file(path):
    """Compute the SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def center_log_rates(pagel_rates):
    """Center per-concept rates in log space per the brief §9.2:
      m̃_c = log(r_c) - weighted_mean(log r_c)
    Returns dict: concept_id -> (log_rate_centered, log_rate_sd).
    """
    log_rates = {}
    for r in pagel_rates:
        cid = r["concept_id"]
        rate = r["rate_per_10kyr"]
        if rate <= 0:
            continue
        log_rates[cid] = math.log(rate) if False else __import__("math").log(rate)
    # Weighted mean: use n_cognates as weight (Pagel uses cognate counts).
    import math
    total_w = 0.0
    total_wr = 0.0
    for r in pagel_rates:
        cid = r["concept_id"]
        if cid in log_rates:
            w = max(1, r.get("n_cognates", 1))
            total_w += w
            total_wr += w * log_rates[cid]
    weighted_mean = total_wr / total_w if total_w > 0 else 0.0
    centered = {}
    for r in pagel_rates:
        cid = r["concept_id"]
        if cid in log_rates:
            centered[cid] = log_rates[cid] - weighted_mean
    return centered, weighted_mean


def build_concept_priors(pagel_rates, vejdemo_predictors, asjp_stability,
                         clock_scale_log_mean, clock_scale_log_sd):
    """Combine evidence sources into per-concept prior rows.

    Returns (rows, header) for concept_priors.tsv.
    """
    import math
    centered, wmean = center_log_rates(pagel_rates)

    # Map predictors by concept.
    vejdemo_by_cid = {p["concept_id"]: p for p in vejdemo_predictors}
    asjp_by_cid = {s["concept_id"]: s for s in asjp_stability}

    header = [
        "concept_id", "scope_id", "unit", "log_rate_mean", "log_rate_sd",
        "n_families", "n_observations", "evidence_grade", "source_key", "version",
    ]
    rows = []

    # GLOBAL scope: combine all evidence.
    for r in pagel_rates:
        cid = r["concept_id"]
        if cid not in centered:
            continue
        # The absolute log-rate = centered + clock_scale.
        log_mean = centered[cid] + clock_scale_log_mean
        # SD: combine Pagel uncertainty with clock-scale uncertainty.
        log_sd = math.sqrt(clock_scale_log_sd ** 2 + 0.35 ** 2)
        n_fam = asjp_by_cid.get(cid, {}).get("n_families", 1)
        n_obs = r.get("n_cognates", 1)
        grade = r.get("evidence_grade", "B")
        src = r.get("source_key", "PAGEL_2007_SYNTH")
        rows.append([cid, "GLOBAL", "per_kyr", f"{log_mean:.6f}",
                     f"{log_sd:.6f}", n_fam, n_obs, grade, src, "1.0.0"])

    # Unknown concepts get a broad fallback prior.
    rows.append(["C_Custom_Unknown", "GLOBAL", "per_kyr",
                 f"{clock_scale_log_mean:.6f}",
                 f"{1.00:.6f}", 0, 0, "D", "DEFAULT", "1.0.0"])

    return rows, header


def build_family_priors(family_data):
    """Build family_priors.tsv rows from family prior data."""
    header = [
        "family_id", "log_clock_mean", "log_clock_sd", "evidence_grade",
        "source_key", "version",
    ]
    rows = []
    for f in family_data:
        rows.append([
            f["family_id"],
            f"{f['log_clock_mean']:.6f}",
            f"{f['log_clock_sd']:.6f}",
            f["evidence_grade"],
            f["source_key"],
            f["version"],
        ])
    return rows, header


def write_manifest(out_dir, files, metadata):
    """Write a manifest.json with checksums and provenance."""
    manifest = {
        "version": metadata.get("version", "1.0.0"),
        "built_at": datetime.now(timezone.utc).isoformat(),
        "sources": metadata.get("sources", []),
        "files": {},
    }
    for fname in files:
        fpath = os.path.join(out_dir, fname)
        if os.path.exists(fpath):
            manifest["files"][fname] = {
                "sha256": sha256_file(fpath),
                "size": os.path.getsize(fpath),
            }
    manifest["warnings"] = metadata.get("warnings", [])
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    return manifest
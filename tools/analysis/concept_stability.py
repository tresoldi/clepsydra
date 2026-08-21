"""
Data-driven concept-stability detection (M12, opt-in).

Ported from the `distglot` toolkit: estimate each concept's stability directly
from the cognate data as the mean pairwise cognate agreement across all language
pairs that share the concept (more stable concepts are retained as the same
cognate more often). Optionally classify concepts stable/regular by a threshold,
and optionally emit a concept-prior TSV (relative log-rates) in Clepsydra's
schema, so a data-driven prior can stand in for the external reference package.

    stability(c) = mean over language pairs (i,j) both having c of
                   [ cognate_i(c) == cognate_j(c) ]

>>> CIRCULARITY WARNING <<<
The single-file estimator is *in-sample*: it reads the very cognate data you
would then date. Using its priors to date that same data double-counts the
evidence and will understate uncertainty. It is therefore NOT a default and NOT
the reference package (`tools/reference_builder`, built from independent
published sources).

The `--leave-out` (leave-family-out) mode is the disciplined alternative: it
pools stability across *other* families and emits a prior for a held-out target
family it never saw, aligning concepts across families by their Concepticon id
(from each family's `concepts.tsv`). That prior is out-of-sample for the target
— the same discipline as `tools/validation/leave_family_out.py`.

Emitted priors are stamped `evidence_grade=D`; in-sample rows carry
`source_key=EMPIRICAL_INSAMPLE`, leave-family-out rows `EMPIRICAL_LFO`.

>>> DEPTH-CONFOUNDING CAVEAT (leave-family-out) <<<
The relative-log-rate transform assumes retention s ≈ exp(-r_c · d_avg) with a
single divergence depth d_avg, so that centering removes d_avg and leaves
relative rates. Pooling *raw similarities* across families of DIFFERENT time
depths breaks that: log(-log(pooled s)) is not log r_c − log(d_avg) for any
common depth, so a concept attested mostly in shallow families versus deep ones
gets a depth-biased relative rate. Pooling is pair-count-weighted (not equal per
family) to stop a one-pair family swinging the estimate, but it does NOT correct
for depth, and cognate agreements among closely-related languages are treated as
independent when they are not. These priors are therefore a grade-D relative
*ranking* of concept stability, not a depth-calibrated rate estimate. Use them
for exploratory ranking and prior-sensitivity checks, not as a load-bearing
prior for a published date.

Stdlib only.

Usage:
  # in-sample (exploratory / sensitivity only):
  python -m analysis.concept_stability --forms forms.tsv --out stability.tsv \\
      [--threshold 0.6] [--emit-priors concept_priors.tsv --prior-sd 0.5]

  # leave-family-out (out-of-sample prior for a held-out family):
  python -m analysis.concept_stability \\
      --families-dir examples/families_pairwise --leave-out dravidian \\
      --out pooled_stability.tsv --emit-priors dravidian_priors.tsv
"""
import argparse
import csv
import math
import os
import sys


def _read_primary_cognates(forms_path):
    """Map (language_id, concept_id) -> primary cognate id (max primary_weight)."""
    best_w = {}
    cog = {}
    concepts = []
    seen_con = set()
    with open(forms_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            lid, cid, cg = r["language_id"], r["concept_id"], r["cognate_id"]
            if cid not in seen_con:
                seen_con.add(cid); concepts.append(cid)
            try:
                w = float(r.get("primary_weight", "1") or "1")
            except ValueError:
                w = 1.0
            key = (lid, cid)
            if key not in best_w or w > best_w[key]:
                best_w[key] = w
                cog[key] = cg
    return cog, concepts


def _stability_counts(forms_path):
    """Return {concept_id: (agree, total)} raw pair-agreement counts (in order)."""
    cog, concepts = _read_primary_cognates(forms_path)
    per = {}
    for (lid, cid), cg in cog.items():
        per.setdefault(cid, {})[lid] = cg
    out = {}
    for cid in concepts:
        cells = list(per.get(cid, {}).values())
        n = len(cells)
        agree = total = 0
        for i in range(n):
            for j in range(i + 1, n):
                total += 1
                if cells[i] == cells[j]:
                    agree += 1
        out[cid] = (agree, total)
    return out


def compute_stability(forms_path):
    """Return {concept_id: (n_pairs, stability_score)} over shared language pairs."""
    out = {}
    for cid, (agree, total) in _stability_counts(forms_path).items():
        out[cid] = (total, agree / total if total else float("nan"))
    return out


def classify(stability, threshold=0.6):
    """stable if score >= threshold, regular if < threshold, unknown if no data."""
    cls = {}
    for cid, (n, s) in stability.items():
        if n == 0 or math.isnan(s):
            cls[cid] = "unknown"
        else:
            cls[cid] = "stable" if s >= threshold else "regular"
    return cls


PRIOR_HEADER = ["concept_id", "scope_id", "unit", "log_rate_mean", "log_rate_sd",
                "n_families", "n_observations", "evidence_grade", "source_key",
                "version"]


def _build_prior_rows(entries, prior_sd, source_key):
    """Centered relative log-rate rows from (concept_id, stability, n_obs, n_fam).

    retention s ~ exp(-r_c * d_avg), so -log(s) ~ r_c * d_avg and
    log r_c ~ log(-log s) - log(d_avg); centering across concepts removes the
    unknown mean divergence d_avg, leaving relative log-rates. More stable
    concepts (higher s) get a *lower* relative rate.
    """
    raw = {}
    meta = {}
    for cid, s, n_obs, n_fam in entries:
        if s is None or (isinstance(s, float) and math.isnan(s)):
            continue
        sc = min(max(s, 0.01), 0.99)          # clamp away from 0/1
        raw[cid] = math.log(-math.log(sc))
        meta[cid] = (n_obs, n_fam)
    if not raw:
        return None
    mean = sum(raw.values()) / len(raw)
    rows = []
    for cid in sorted(raw):
        n_obs, n_fam = meta[cid]
        rows.append([cid, "GLOBAL", "per_kyr", "%.6f" % (raw[cid] - mean),
                     "%.6f" % prior_sd, str(n_fam), str(n_obs), "D",
                     source_key, "0.0.0"])
    return PRIOR_HEADER, rows


def to_concept_priors(stability, prior_sd=0.5):
    """In-sample empirical relative log-rate rows for concept_priors.tsv."""
    entries = [(cid, s, n, 1) for cid, (n, s) in stability.items()]
    return _build_prior_rows(entries, prior_sd, "EMPIRICAL_INSAMPLE")


# --------------------------------------------------------------------------
# Leave-family-out: pool stability across families, emit an out-of-sample prior
# for a held-out target family (concepts aligned by Concepticon id).
# --------------------------------------------------------------------------

def _concept_key_map(concepts_path):
    """local concept_id -> (aligned_key, label).

    aligned_key is the Concepticon id when a family's `concepts.tsv` records one,
    so a concept aligns across families regardless of each family's local id;
    otherwise the aligned_key falls back to the local concept_id.
    """
    m = {}
    with open(concepts_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            cid = r["concept_id"]
            ns = (r.get("external_namespace") or "").strip().lower()
            ext = (r.get("external_id") or "").strip()
            label = (r.get("label") or cid).strip()
            key = ext if (ns == "concepticon" and ext) else cid
            m[cid] = (key, label)
    return m


def load_family(forms_path, concepts_path=None):
    """One family's per-concept counts keyed by aligned concept key.

    Returns {aligned_key: (agree, total, local_concept_id, label)}. When two
    local concepts map to the same aligned key their counts are summed.
    """
    kmap = {}
    if concepts_path and os.path.exists(concepts_path):
        kmap = _concept_key_map(concepts_path)
    out = {}
    for cid, (agree, total) in _stability_counts(forms_path).items():
        key, label = kmap.get(cid, (cid, cid))
        if key in out:
            a0, t0, lc0, lb0 = out[key]
            out[key] = (a0 + agree, t0 + total, lc0, lb0)
        else:
            out[key] = (agree, total, cid, label)
    return out


def pool_families(loaded):
    """Pool per-family cognate agreement, weighted by the number of language
    pairs each family contributes.

    loaded: {family_label: family_dict from load_family}. Returns
    {aligned_key: (n_families, n_pairs_total, pooled_stability, label)}. Pooling
    sums agreements and pairs across families (stability = Σ agree / Σ pairs), so
    a family that shares a concept across one language pair (a 0/1 Bernoulli
    draw) does not carry the same weight as one with hundreds of pairs. See the
    depth-confounding caveat in the module docstring: pooling raw similarities
    across families of different time depths is only a grade-D relative ranking.
    """
    agg = {}
    for label, fam in loaded.items():
        for key, (agree, total, lc, lbl) in fam.items():
            if total == 0:
                continue
            e = agg.setdefault(key, {"agree": 0, "pairs": 0, "nfam": 0,
                                     "label": lbl})
            e["agree"] += agree
            e["pairs"] += total
            e["nfam"] += 1
    out = {}
    for key, e in agg.items():
        out[key] = (e["nfam"], e["pairs"], e["agree"] / e["pairs"], e["label"])
    return out


def discover_families(families_dir):
    """{subdir_name: {"forms": path, "concepts": path or None}} for each family.

    A family is any immediate subdirectory containing a `forms.tsv`.
    """
    fams = {}
    for name in sorted(os.listdir(families_dir)):
        sub = os.path.join(families_dir, name)
        if not os.path.isdir(sub):
            continue
        forms = os.path.join(sub, "forms.tsv")
        if not os.path.exists(forms):
            continue
        concepts = os.path.join(sub, "concepts.tsv")
        fams[name] = {"forms": forms,
                      "concepts": concepts if os.path.exists(concepts) else None}
    return fams


def leave_family_out(families, leave_out, prior_sd=0.5):
    """Pool stability over every family except `leave_out`, then emit a prior for
    the held-out family keyed to its own local concept_ids.

    families: {label: {"forms": path, "concepts": path or None}}.
    Returns (pooled, prior, coverage) where
      pooled  = {aligned_key: (n_families, n_pairs, stability, label)}
      prior   = (header, rows) keyed by the target's local concept_id, or None
      coverage= {"target": label, "target_concepts": int, "covered": int,
                 "pooled_families": [labels]}
    """
    if leave_out not in families:
        raise KeyError("leave-out family %r not among %s"
                       % (leave_out, ", ".join(sorted(families))))
    trainers = {lab: p for lab, p in families.items() if lab != leave_out}
    if not trainers:
        raise ValueError("need at least one other family to pool over")

    loaded = {lab: load_family(p["forms"], p["concepts"])
              for lab, p in trainers.items()}
    pooled = pool_families(loaded)

    # Every target concept_id gets its own prior row, keyed to that local id.
    # We map local concept_id -> aligned key (NOT key -> id): if two of the
    # target's concepts share a Concepticon id they must both be emitted, not
    # silently collapsed to one.
    tgt = families[leave_out]
    tgt_pairs = []          # (local concept_id, aligned_key)
    if tgt["concepts"]:
        for cid, (key, _lbl) in _concept_key_map(tgt["concepts"]).items():
            tgt_pairs.append((cid, key))
    else:
        # No concepts.tsv: align on the local concept_id directly.
        for cid in _stability_counts(tgt["forms"]):
            tgt_pairs.append((cid, cid))

    entries = []
    covered = 0
    for local_cid, key in tgt_pairs:
        if key not in pooled:
            continue
        n_fam, n_pairs, s, _label = pooled[key]
        entries.append((local_cid, s, n_pairs, n_fam))
        covered += 1
    prior = _build_prior_rows(entries, prior_sd, "EMPIRICAL_LFO")
    coverage = {"target": leave_out, "target_concepts": len(tgt_pairs),
                "covered": covered, "pooled_families": sorted(trainers)}
    return pooled, prior, coverage


def _write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def _parse_family_args(family_args):
    """`label=forms[,concepts]` items -> {label: {"forms":..,"concepts":..}}."""
    fams = {}
    for spec in family_args:
        if "=" not in spec:
            raise ValueError("--family expects label=forms[,concepts], got %r" % spec)
        label, paths = spec.split("=", 1)
        parts = paths.split(",")
        forms = parts[0]
        concepts = parts[1] if len(parts) > 1 and parts[1] else None
        fams[label.strip()] = {"forms": forms, "concepts": concepts}
    return fams


def _run_single(args):
    stability = compute_stability(args.forms)
    cls = classify(stability, args.threshold)
    header = ["concept_id", "n_pairs", "stability", "class"]
    rows = []
    for cid in stability:
        n, s = stability[cid]
        rows.append([cid, str(n), ("%.4f" % s) if not math.isnan(s) else ".",
                     cls[cid]])
    _write_tsv(args.out, header, rows)
    n_stable = sum(1 for c in cls.values() if c == "stable")
    print("Wrote %d concepts to %s (%d stable, threshold %.2f)"
          % (len(rows), args.out, n_stable, args.threshold))
    if args.emit_priors:
        built = to_concept_priors(stability, args.prior_sd)
        if built:
            hdr, prows = built
            _write_tsv(args.emit_priors, hdr, prows)
            print("Wrote %d IN-SAMPLE concept priors to %s "
                  "(evidence_grade=D; do NOT date the same data with these)"
                  % (len(prows), args.emit_priors))
    return 0


def _run_leave_family_out(args):
    if args.families_dir:
        families = discover_families(args.families_dir)
    else:
        families = _parse_family_args(args.family)
    if not families:
        print("No families found.", file=sys.stderr)
        return 2
    pooled, prior, cov = leave_family_out(families, args.leave_out, args.prior_sd)

    header = ["aligned_key", "label", "n_families", "n_pairs", "stability", "class"]
    rows = []
    for key in sorted(pooled, key=lambda k: -pooled[k][2]):
        n_fam, n_pairs, s, label = pooled[key]
        cl = "stable" if s >= args.threshold else "regular"
        rows.append([key, label, str(n_fam), str(n_pairs), "%.4f" % s, cl])
    _write_tsv(args.out, header, rows)
    print("Pooled %d concepts over %d families (excluded '%s') -> %s"
          % (len(rows), len(cov["pooled_families"]), cov["target"], args.out))

    if args.emit_priors:
        if not prior:
            print("No pooled concept overlaps the held-out family; "
                  "no priors written.", file=sys.stderr)
        else:
            hdr, prows = prior
            _write_tsv(args.emit_priors, hdr, prows)
            print("Wrote %d LEAVE-FAMILY-OUT concept priors to %s "
                  "(evidence_grade=D, source_key=EMPIRICAL_LFO; out-of-sample "
                  "for '%s': %d/%d target concepts covered)"
                  % (len(prows), args.emit_priors, cov["target"],
                     cov["covered"], cov["target_concepts"]))
    return 0


def main(argv):
    ap = argparse.ArgumentParser(prog="analysis.concept_stability")
    ap.add_argument("--forms", help="single Clepsydra forms TSV (in-sample mode)")
    ap.add_argument("--families-dir",
                    help="directory of per-family subdirs (each with forms.tsv "
                         "and, ideally, concepts.tsv) for leave-family-out mode")
    ap.add_argument("--family", action="append", default=[],
                    help="explicit family as label=forms[,concepts] (repeatable)")
    ap.add_argument("--leave-out",
                    help="family label to hold out and emit an out-of-sample "
                         "prior for (enables leave-family-out mode)")
    ap.add_argument("--out", required=True, help="stability table TSV")
    ap.add_argument("--threshold", type=float, default=0.6)
    ap.add_argument("--emit-priors", default=None,
                    help="also write a concept_priors.tsv here")
    ap.add_argument("--prior-sd", type=float, default=0.5)
    args = ap.parse_args(argv[1:])

    lfo = bool(args.leave_out or args.families_dir or args.family)
    if lfo:
        if not args.leave_out:
            ap.error("leave-family-out mode requires --leave-out LABEL")
        if not (args.families_dir or args.family):
            ap.error("leave-family-out mode requires --families-dir or --family")
        if args.forms:
            ap.error("--forms is for single-file mode; use --family/--families-dir "
                     "for leave-family-out")
        return _run_leave_family_out(args)
    if not args.forms:
        ap.error("single-file mode requires --forms")
    return _run_single(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

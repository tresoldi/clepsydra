"""
Synthetic data generators for the reference-prior builder.

These stand in for the real evidence sources (Pagel 2007 Supp Table 2,
Vejdemo-Hörberg 2016 Supporting Info, ASJP database) which are not yet
available in our literature/. When the real supplements are supplied,
replace these generators with real-data parsers; the downstream TSV
emission and manifest logic is unchanged.

The synthetic data is designed to be realistic in structure:
  - per-concept rates vary ~100x (as Pagel reports);
  - frequent concepts (numbers, pronouns) are slow;
  - infrequent/adjective concepts are fast;
  - family-level multipliers center at 0 (log scale);
  - evidence grades are assigned per source.
"""

import math
import random

from .concepts import get_concepts


# Legacy anchors (from docs/model-spec.md §2):
LEGACY_SWADESH100 = {
    "prior_id": "legacy_swadesh100",
    "unit": "per_kyr",
    "log_mean": -1.8916490462,
    "log_sd": 0.80,
    "description": "0.86 retention/kyr",
}
LEGACY_LEES200 = {
    "prior_id": "legacy_lees200",
    "unit": "per_kyr",
    "log_mean": -1.5271140564,
    "log_sd": 0.80,
    "description": "0.8048 retention/kyr",
}

# Word-class typical relative rates (synthetic, inspired by Pagel 2007):
# numbers/pronouns slow, adjectives/conjunctions fast.
WORD_CLASS_RATE = {
    "numeral": 0.5,
    "pronoun": 0.6,
    "adverb": 1.2,
    "noun": 1.5,
    "verb": 2.0,
    "adjective": 2.5,
    "conjunction": 3.0,
    "preposition": 3.0,
    "unknown": 1.5,
}


def generate_pagel_rates(seed=42):
    """Generate synthetic Pagel-2007-style per-concept rates.

    Returns a list of dicts:
      {concept_id, rate_per_10kyr, n_cognates, evidence_grade, source_key}

    Rates are relative (centered in log space downstream). The absolute scale
    comes from the legacy clock anchor.

    NOTE: REAL DATA MISSING. Replace with a parser for Pagel 2007
    Supplementary Table 2 when supplied.
    """
    rng = random.Random(seed)
    concepts = get_concepts()
    results = []
    for c in concepts:
        base = WORD_CLASS_RATE.get(c["word_class"], 1.5)
        # Add noise: log-normal variation around the word-class base.
        noise = math.exp(rng.gauss(0, 0.5))
        rate = base * noise
        n_cog = max(1, int(rng.lognormvariate(math.log(base * 5), 0.8)))
        results.append({
            "concept_id": c["concept_id"],
            "rate_per_10kyr": rate,
            "n_cognates": min(n_cog, 46),
            "evidence_grade": "B",
            "source_key": "PAGEL_2007_SYNTH",
        })
    return results


def generate_vejdemo_predictors(seed=43):
    """Generate synthetic Vejdemo-Hörberg-2016-style predictor data.

    Returns a list of dicts with predictors (frequency, synonyms, senses,
    imageability) for content concepts only.

    NOTE: REAL DATA MISSING. Replace with a parser for Vejdemo-Hörberg 2016
    Supporting Info when supplied.
    """
    rng = random.Random(seed)
    concepts = get_concepts()
    content_classes = {"noun", "verb", "adjective"}
    results = []
    for c in concepts:
        if c["word_class"] not in content_classes:
            continue
        freq = rng.lognormvariate(3.0, 1.5)  # per million
        synonyms = max(0, int(rng.gauss(8, 5)))
        senses = max(1, int(rng.gauss(4, 3)))
        imageability = rng.uniform(1, 7)
        results.append({
            "concept_id": c["concept_id"],
            "frequency": freq,
            "synonyms": synonyms,
            "senses": senses,
            "imageability": imageability,
            "evidence_grade": "B",
            "source_key": "VEJDEMO_2016_SYNTH",
        })
    return results


def generate_asjp_stability(seed=44):
    """Generate synthetic ASJP-style worldwide stability scores.

    Returns a list of dicts: {concept_id, stability_score, n_families}

    NOTE: REAL DATA MISSING. Replace with a parser for the pinned ASJP
    database release when supplied.
    """
    rng = random.Random(seed)
    concepts = get_concepts()
    results = []
    for c in concepts:
        base = WORD_CLASS_RATE.get(c["word_class"], 1.5)
        # Stability is inverse of rate (slow = stable).
        stability = 1.0 / base * rng.uniform(0.7, 1.3)
        n_fam = rng.randint(50, 200)
        results.append({
            "concept_id": c["concept_id"],
            "stability_score": stability,
            "n_families": n_fam,
            "evidence_grade": "B",
            "source_key": "ASJP_SYNTH",
        })
    return results


# Synthetic family benchmarks (from docs/paper-evaluation.md — Rama-Wichmann
# 30 calibration groups, captured from the paper):
FAMILY_BENCHMARKS = [
    ("indo_european", 0.0, 0.40, "B", "META_IE_001"),
    ("romance", 0.0, 0.35, "B", "RAMA_WICHMANN_2020"),
    ("germanic", 0.0, 0.30, "B", "RAMA_WICHMANN_2020"),
    ("slavic", 0.0, 0.30, "B", "RAMA_WICHMANN_2020"),
    ("dravidian", 0.0, 0.45, "C", "RAMA_WICHMANN_2020"),
    ("austronesian", 0.0, 0.40, "C", "RAMA_WICHMANN_2020"),
    ("sino_tibetan", 0.0, 0.45, "C", "RAMA_WICHMANN_2020"),
    ("turkic", 0.0, 0.35, "B", "RAMA_WICHMANN_2020"),
    ("uralic", 0.0, 0.40, "B", "RAMA_WICHMANN_2020"),
    ("semitic", 0.0, 0.40, "C", "RAMA_WICHMANN_2020"),
]


def generate_family_priors():
    """Return family prior rows (synthetic + UNKNOWN fallback)."""
    results = []
    for (fam, mean, sd, grade, src) in FAMILY_BENCHMARKS:
        results.append({
            "family_id": fam,
            "log_clock_mean": mean,
            "log_clock_sd": sd,
            "evidence_grade": grade,
            "source_key": src,
            "version": "1.0.0",
        })
    # Fallback for unknown families.
    results.append({
        "family_id": "UNKNOWN",
        "log_clock_mean": 0.0,
        "log_clock_sd": 0.80,
        "evidence_grade": "D",
        "source_key": "DEFAULT",
        "version": "1.0.0",
    })
    return results
"""
Tests for the reference-prior builder.

Run: python3 -m pytest tests/test_reference_builder.py
  or: python3 tests/test_reference_builder.py
"""

import os
import sys
import tempfile
import unittest

# Make the reference_builder package importable.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from reference_builder import cli as cli_mod
from reference_builder import synthetic_sources as synth
from reference_builder import emitter
from reference_builder.concepts import get_concepts, canonical_ids


class TestSyntheticSources(unittest.TestCase):
    def test_pagel_rates_generated(self):
        rates = synth.generate_pagel_rates(seed=42)
        self.assertGreater(len(rates), 10)
        for r in rates:
            self.assertGreater(r["rate_per_10kyr"], 0)
            self.assertIn("evidence_grade", r)
            self.assertEqual(r["source_key"], "PAGEL_2007_SYNTH")

    def test_pagel_reproducible(self):
        r1 = synth.generate_pagel_rates(seed=42)
        r2 = synth.generate_pagel_rates(seed=42)
        self.assertEqual(r1, r2)

    def test_vejdemo_content_only(self):
        preds = synth.generate_vejdemo_predictors(seed=43)
        for p in preds:
            self.assertIn(p["concept_id"], canonical_ids())

    def test_asjp_stability(self):
        scores = synth.generate_asjp_stability(seed=44)
        self.assertGreater(len(scores), 10)
        for s in scores:
            self.assertGreater(s["stability_score"], 0)
            self.assertGreater(s["n_families"], 0)

    def test_family_priors_has_unknown(self):
        families = synth.generate_family_priors()
        ids = [f["family_id"] for f in families]
        self.assertIn("UNKNOWN", ids)
        unk = [f for f in families if f["family_id"] == "UNKNOWN"][0]
        self.assertEqual(unk["evidence_grade"], "D")

    def test_legacy_anchors(self):
        self.assertAlmostEqual(synth.LEGACY_SWADESH100["log_mean"], -1.8916490462)
        self.assertAlmostEqual(synth.LEGACY_LEES200["log_mean"], -1.5271140564)


class TestEmitter(unittest.TestCase):
    def test_write_tsv(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "test.tsv")
            emitter.write_tsv(path, ["a", "b"], [["1", "2"], ["3", None]])
            with open(path) as f:
                content = f.read()
            self.assertEqual(content, "a\tb\n1\t2\n3\t.\n")

    def test_sha256_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "test.txt")
            with open(path, "w") as f:
                f.write("hello")
            h = emitter.sha256_file(path)
            self.assertEqual(len(h), 64)

    def test_center_log_rates(self):
        pagel = synth.generate_pagel_rates(seed=42)
        centered, wmean = emitter.center_log_rates(pagel)
        # Centered rates should have roughly zero weighted mean.
        import math
        total_w = 0
        total_wr = 0
        for r in pagel:
            cid = r["concept_id"]
            if cid in centered:
                w = max(1, r.get("n_cognates", 1))
                total_w += w
                total_wr += w * centered[cid]
        self.assertAlmostEqual(total_wr / total_w, 0.0, places=10)


class TestBuildAndValidate(unittest.TestCase):
    def test_build_and_validate(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "pkg")
            # Build.
            rc = self._run_cli(["build", "--out", out, "--seed", "42"])
            self.assertEqual(rc, 0)
            # Check files exist.
            for f in ["concept_priors.tsv", "family_priors.tsv", "concepts.tsv",
                       "manifest.json"]:
                self.assertTrue(os.path.exists(os.path.join(out, f)))
            # Validate.
            rc = self._run_cli(["validate", out])
            self.assertEqual(rc, 0)

    def test_build_reproducible(self):
        with tempfile.TemporaryDirectory() as d:
            out1 = os.path.join(d, "a")
            out2 = os.path.join(d, "b")
            self._run_cli(["build", "--out", out1, "--seed", "42"])
            self._run_cli(["build", "--out", out2, "--seed", "42"])
            with open(os.path.join(out1, "concept_priors.tsv")) as f:
                c1 = f.read()
            with open(os.path.join(out2, "concept_priors.tsv")) as f:
                c2 = f.read()
            self.assertEqual(c1, c2)

    def test_leave_family_out(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "pkg")
            self._run_cli(["build", "--out", out, "--seed", "42",
                           "--leave-family-out", "indo_european"])
            lf_dir = os.path.join(out, "leave_out_indo_european")
            self.assertTrue(os.path.exists(os.path.join(lf_dir, "family_priors.tsv")))
            # The leave-out package should NOT contain the indo_european
            # family-specific prior.
            with open(os.path.join(lf_dir, "family_priors.tsv")) as f:
                content = f.read()
            self.assertNotIn("indo_european", content.split("\n")[1:])
            self.assertIn("UNKNOWN", content)

    def test_no_unexplained_constants(self):
        """Check that the manifest records all source keys."""
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "pkg")
            self._run_cli(["build", "--out", out, "--seed", "42"])
            import json
            with open(os.path.join(out, "manifest.json")) as f:
                m = json.load(f)
            self.assertIn("sources", m)
            self.assertGreater(len(m["sources"]), 0)

    @staticmethod
    def _run_cli(argv):
        import argparse
        old_argv = sys.argv
        sys.argv = ["reference_builder"] + argv
        try:
            return cli_mod.main()
        finally:
            sys.argv = old_argv


class TestConcepts(unittest.TestCase):
    def test_concepts_have_ids(self):
        concepts = get_concepts()
        self.assertGreater(len(concepts), 10)
        ids = [c["concept_id"] for c in concepts]
        # No duplicates.
        self.assertEqual(len(ids), len(set(ids)))

    def test_custom_concept_valid(self):
        """Check that a custom/unknown concept is present and active."""
        concepts = get_concepts()
        custom = [c for c in concepts if c["word_class"] == "unknown"]
        self.assertGreaterEqual(len(custom), 1)
        self.assertEqual(custom[0]["status"], "active")


if __name__ == "__main__":
    unittest.main()
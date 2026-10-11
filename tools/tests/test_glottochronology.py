"""
Tests for the classical lexicostatistical dating baselines.

Run: python3 tools/tests/test_glottochronology.py
"""
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "glottochronology"))

import dating as g


class TestFormulas(unittest.TestCase):
    def test_swadesh_half_shared(self):
        # t = ln C / (2 ln r): C = 0.5, r = 0.86 -> 2,298 years.
        self.assertAlmostEqual(g.date_swadesh(0.5, 0.86, 0, 0),
                               math.log(0.5) / (2 * math.log(0.86)) * 1000, places=6)
        self.assertAlmostEqual(g.date_swadesh(0.5, 0.86, 0, 0), 2297.9, places=1)

    def test_swadesh_dated_tip(self):
        # A tip dated tau has been separate tau fewer years: t = (ln C/ln r + tau)/2.
        t0 = g.date_swadesh(0.6, 0.805, 0, 0)
        t1 = g.date_swadesh(0.6, 0.805, 500, 0)
        self.assertAlmostEqual(t1 - t0, 250.0, places=6)

    def test_starostin_pair(self):
        # C = exp(-2 lam C t^2): C = 0.5 -> t = sqrt(ln 2 / (0.1 * 0.5)) millennia.
        t = g.date_starostin(0.5, 0.05, 0, 0, 2.0)
        self.assertAlmostEqual(t, math.sqrt(math.log(2) / 0.05) * 1000, places=6)
        C = 0.5
        self.assertAlmostEqual(math.exp(-2 * 0.05 * C * (t / 1000) ** 2), C, places=9)

    def test_identical_lists_date_to_the_tips(self):
        self.assertEqual(g.date_swadesh(1.0, 0.86, 100, 300), 200.0)
        self.assertEqual(g.date_starostin(1.0, 0.05, 100, 300), 200.0)


class TestTreeAndCounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        p = os.path.join(self.tmp.name, "forms.tsv")
        with open(p, "w", encoding="utf-8") as f:
            f.write("language_id\tconcept_id\tcognate_id\n")
            rows = [("a", "C1", "x"), ("b", "C1", "x"), ("c", "C1", "y"),
                    ("a", "C2", "p"), ("b", "C2", "q"), ("c", "C2", "q")]
            for r in rows:
                f.write("\t".join(r) + "\n")
        self.forms = g.read_forms(p)

    def tearDown(self):
        self.tmp.cleanup()

    def test_parse_newick_labels(self):
        ch, root = g.parse_newick("((a:1,b:1)n1:2,c:3)root;")
        self.assertEqual(root, "root")
        self.assertEqual(sorted(g.tips_under(ch, "n1")), ["a", "b"])

    def test_match_counts(self):
        self.assertEqual(g.match_counts(self.forms, "a", "b", ["C1", "C2"]), (1, 2))
        self.assertEqual(g.match_counts(self.forms, "b", "c", ["C1", "C2"]), (1, 2))

    def test_node_pools_cross_pairs(self):
        ch, _ = g.parse_newick("((a,b)n1,c)root;")
        fn = lambda C, x, y: g.date_swadesh(C, 0.86, x, y)
        C, n, _t = g.node_date(self.forms, ch, "root", ["C1", "C2"], fn, {})
        # a-c: 0/2, b-c: 1/2 -> pooled 1/4.
        self.assertEqual((C, n), (0.25, 4))


if __name__ == "__main__":
    unittest.main(verbosity=2)

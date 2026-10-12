"""Tests for the GLED/ASJP benchmark helpers. Run: python3 tools/tests/test_gled_benchmark.py"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gled_benchmark import asjp as A
from gled_benchmark.tree import newick, upgma


class TestAsjp(unittest.TestCase):
    def test_levenshtein(self):
        self.assertEqual(A.lev("kitten", "sitting"), 3)
        self.assertEqual(A.lev("", "abc"), 3)
        self.assertAlmostEqual(A.ldn("abc", "abd"), 1 / 3)

    def test_date_formula(self):
        # s = s0 gives age 0; s = s0 * r gives one millennium divided by two.
        self.assertAlmostEqual(A.date_from_similarity(A.S0), 0.0, places=6)
        self.assertAlmostEqual(A.date_from_similarity(A.S0 * A.R ** 2), 1000.0, places=6)

    def test_ldnd_identical_lists(self):
        w = {str(i): "abcd"[: 1 + i % 4] + "xyz"[: i % 3] for i in range(12)}
        self.assertEqual(A.ldnd(w, dict(w)), 0.0)


class TestTree(unittest.TestCase):
    def test_upgma_groups_the_close_pair(self):
        t = upgma(list("abcd"), lambda a, b: {"ab": .2, "cd": .4}.get(a + b, 1.0))
        clades = {}
        s = newick(t, clades=clades)
        self.assertTrue(s.endswith("root;"))
        self.assertIn(frozenset("ab"), clades.values())
        self.assertIn(frozenset("cd"), clades.values())


if __name__ == "__main__":
    unittest.main()

"""
Tests for the plain cognate-matrix reader (M12).

Run: python3 tools/tests/test_matrix.py
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from datasets import matrix


def _write(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


class TestWideMatrix(unittest.TestCase):
    def test_basic_wide(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.csv")
            _write(p, "language,water,fire,hand\n"
                      "Hmong,1,3,2\n"
                      "Mien,1,4,2\n"
                      "She,5,?,2\n")   # She/fire missing
            lr, cr, fr = matrix.read_wide(p)
            self.assertEqual(len(lr), 3)         # three languages
            self.assertEqual(len(cr), 3)         # three concepts
            self.assertEqual(len(fr), 8)         # 9 cells - 1 missing
            # cognate ids are namespaced by concept
            water = [row for row in fr if row[2] == "Cwater"]
            self.assertTrue(all(row[5] == "Cwater-1" for row in water
                                if row[1] in ("Hmong", "Mien")))
            # Hmong 'water' and Mien 'water' share a class; She differs
            she_water = [r for r in water if r[1] == "She"][0]
            self.assertEqual(she_water[5], "Cwater-5")

    def test_missing_and_delimiter_sniff(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.tsv")
            _write(p, "lang\twater\tfire\nA\t1\t-\nB\t.\t2\n")  # tab-sniffed
            lr, cr, fr = matrix.read_wide(p)
            self.assertEqual(len(fr), 2)        # A/water and B/fire only

    def test_synonyms_split(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.csv")
            _write(p, "language,water\nA,1;7\nB,1\n")
            lr, cr, fr = matrix.read_wide(p, synonym_sep=";")
            a = [r for r in fr if r[1] == "A"]
            self.assertEqual(len(a), 2)                     # two candidate forms
            self.assertAlmostEqual(float(a[0][6]), 0.5)     # weight split
            self.assertEqual({r[5] for r in a}, {"Cwater-1", "Cwater-7"})

    def test_duplicate_language_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.csv")
            _write(p, "language,water\nA,1\nA,2\n")
            with self.assertRaises(ValueError):
                matrix.read_wide(p)


class TestLongMatrix(unittest.TestCase):
    def test_basic_long(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "l.csv")
            _write(p, "Language,Concept,Cognateset\n"
                      "Hmong,water,1\nMien,water,1\nHmong,fire,3\n")
            lr, cr, fr = matrix.read_long(p, "Language", "Concept", "Cognateset")
            self.assertEqual(len(lr), 2)
            self.assertEqual(len(cr), 2)
            self.assertEqual(len(fr), 3)
            self.assertEqual(sorted(r[5] for r in fr),
                             ["Cfire-3", "Cwater-1", "Cwater-1"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

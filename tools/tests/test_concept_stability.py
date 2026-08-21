"""
Tests for the data-driven concept-stability tool (M12, opt-in).

Run: python3 tools/tests/test_concept_stability.py
"""
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from analysis import concept_stability as cs


FORMS_HEADER = ("analysis_id\tlanguage_id\tconcept_id\tform_id\tform\t"
                "cognate_id\tprimary_weight\tcognacy_confidence\t"
                "inheritance_probability\tloan_probability\tstatus\n")


def _write_forms(path, cells):
    """cells: list of (language, concept, cognate)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(FORMS_HEADER)
        for i, (lang, con, cog) in enumerate(cells):
            f.write("m\t%s\t%s\tf%d\t.\t%s\t1.0\t1.0\t0.95\t0.02\tobserved\n"
                    % (lang, con, i, cog))


class TestConceptStability(unittest.TestCase):
    def test_scores_and_classification(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "forms.tsv")
            # 'stable': all four languages share cognate 1 -> stability 1.0
            # 'mixed': two classes split 2/2 -> agreement among 6 pairs = 2/6
            # 'diverse': all different -> stability 0.0
            cells = []
            for lang in ("A", "B", "C", "D"):
                cells.append((lang, "stable", "1"))
            cells += [("A", "mixed", "1"), ("B", "mixed", "1"),
                      ("C", "mixed", "2"), ("D", "mixed", "2")]
            cells += [("A", "diverse", "1"), ("B", "diverse", "2"),
                      ("C", "diverse", "3"), ("D", "diverse", "4")]
            _write_forms(p, cells)
            st = cs.compute_stability(p)
            self.assertAlmostEqual(st["stable"][1], 1.0)
            self.assertAlmostEqual(st["mixed"][1], 2.0 / 6.0)
            self.assertAlmostEqual(st["diverse"][1], 0.0)
            cl = cs.classify(st, threshold=0.6)
            self.assertEqual(cl["stable"], "stable")
            self.assertEqual(cl["mixed"], "regular")
            self.assertEqual(cl["diverse"], "regular")

    def test_priors_centered_and_monotone(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "forms.tsv")
            # three concepts with distinct stabilities
            cells = []
            for lang in ("A", "B", "C", "D"):
                cells.append((lang, "hi", "1"))                 # stability 1.0
            cells += [("A", "mid", "1"), ("B", "mid", "1"),
                      ("C", "mid", "2"), ("D", "mid", "2")]     # 0.33
            cells += [("A", "lo", "1"), ("B", "lo", "2"),
                      ("C", "lo", "3"), ("D", "lo", "4")]       # 0.0
            _write_forms(p, cells)
            st = cs.compute_stability(p)
            hdr, rows = cs.to_concept_priors(st, prior_sd=0.5)
            lr = {r[0]: float(r[3]) for r in rows}
            # centered: relative log-rates average ~0
            self.assertAlmostEqual(sum(lr.values()) / len(lr), 0.0, places=6)
            # more stable => lower rate: hi < mid < lo
            self.assertLess(lr["hi"], lr["mid"])
            self.assertLess(lr["mid"], lr["lo"])
            # schema stamps the in-sample provenance
            self.assertTrue(all(r[7] == "D" and r[8] == "EMPIRICAL_INSAMPLE" for r in rows))
            self.assertEqual(hdr[0], "concept_id")

    def test_singleton_concept_is_unknown(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "forms.tsv")
            _write_forms(p, [("A", "solo", "1")])   # only one language
            st = cs.compute_stability(p)
            self.assertEqual(st["solo"][0], 0)
            self.assertTrue(math.isnan(st["solo"][1]))
            self.assertEqual(cs.classify(st)["solo"], "unknown")


CONCEPTS_HEADER = ("concept_id\tlabel\texternal_namespace\texternal_id\t"
                   "word_class\tstatus\n")


def _write_concepts(path, rows):
    """rows: list of (concept_id, label, concepticon_id)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(CONCEPTS_HEADER)
        for cid, label, ext in rows:
            f.write("%s\t%s\tConcepticon\t%s\tunknown\tactive\n"
                    % (cid, label, ext))


def _mk_family(d, name, forms_cells, concept_rows):
    sub = os.path.join(d, name)
    os.makedirs(sub)
    _write_forms(os.path.join(sub, "forms.tsv"), forms_cells)
    _write_concepts(os.path.join(sub, "concepts.tsv"), concept_rows)


class TestLeaveFamilyOut(unittest.TestCase):
    def _three_families(self, d):
        # Concepticon 100 = HAND (stable everywhere), 200 = STONE (diverse).
        # Each family uses a DIFFERENT local concept_id for the same Concepticon
        # id, so alignment must go through concepts.tsv, not the local id.
        _mk_family(
            d, "alpha",
            [("A", "C1_hand", "h"), ("B", "C1_hand", "h"), ("C", "C1_hand", "h"),
             ("A", "C2_stone", "1"), ("B", "C2_stone", "2"), ("C", "C2_stone", "3")],
            [("C1_hand", "HAND", "100"), ("C2_stone", "STONE", "200")])
        _mk_family(
            d, "beta",
            [("X", "C9_hand", "h"), ("Y", "C9_hand", "h"), ("Z", "C9_hand", "q"),
             ("X", "C8_stone", "1"), ("Y", "C8_stone", "2"), ("Z", "C8_stone", "2")],
            [("C9_hand", "HAND", "100"), ("C8_stone", "STONE", "200")])
        # target family: yet another local namespace for the same concepts.
        _mk_family(
            d, "gamma",
            [("P", "G_hand", "z"), ("Q", "G_hand", "z"),
             ("P", "G_stone", "1"), ("Q", "G_stone", "1")],
            [("G_hand", "HAND", "100"), ("G_stone", "STONE", "200")])

    def test_pool_excludes_target_and_aligns_by_concepticon(self):
        with tempfile.TemporaryDirectory() as d:
            self._three_families(d)
            fams = cs.discover_families(d)
            self.assertEqual(set(fams), {"alpha", "beta", "gamma"})
            pooled, prior, cov = cs.leave_family_out(fams, "gamma")
            # pooled over alpha+beta only, keyed by Concepticon id
            self.assertEqual(cov["pooled_families"], ["alpha", "beta"])
            self.assertIn("100", pooled)
            self.assertIn("200", pooled)
            # HAND: alpha 1.0, beta 1/3 (Z differs) -> mean 2/3; seen in 2 families
            n_fam, n_pairs, s_hand, label = pooled["100"]
            self.assertEqual(n_fam, 2)
            self.assertEqual(label, "HAND")
            self.assertAlmostEqual(s_hand, (1.0 + 1.0 / 3.0) / 2.0)
            # STONE: alpha 0.0, beta 1/3 -> mean 1/6 (more diverse than HAND)
            self.assertAlmostEqual(pooled["200"][2], (0.0 + 1.0 / 3.0) / 2.0)

    def test_prior_keyed_to_target_local_ids_and_out_of_sample(self):
        with tempfile.TemporaryDirectory() as d:
            self._three_families(d)
            fams = cs.discover_families(d)
            _pooled, prior, cov = cs.leave_family_out(fams, "gamma")
            self.assertIsNotNone(prior)
            hdr, rows = prior
            by_cid = {r[0]: r for r in rows}
            # rows are keyed by gamma's OWN local concept ids
            self.assertEqual(set(by_cid), {"G_hand", "G_stone"})
            # provenance stamped leave-family-out, evidence grade D
            self.assertTrue(all(r[7] == "D" and r[8] == "EMPIRICAL_LFO"
                                for r in rows))
            # more stable HAND -> lower relative rate than diverse STONE
            self.assertLess(float(by_cid["G_hand"][3]),
                            float(by_cid["G_stone"][3]))
            # centered across the covered concepts
            self.assertAlmostEqual(
                sum(float(r[3]) for r in rows) / len(rows), 0.0, places=6)
            # n_families column reflects the pooled family count
            self.assertEqual(by_cid["G_hand"][5], "2")
            self.assertEqual(cov["covered"], 2)
            self.assertEqual(cov["target_concepts"], 2)

    def test_explicit_family_specs(self):
        with tempfile.TemporaryDirectory() as d:
            self._three_families(d)
            specs = ["%s=%s,%s" % (n, os.path.join(d, n, "forms.tsv"),
                                   os.path.join(d, n, "concepts.tsv"))
                     for n in ("alpha", "beta", "gamma")]
            fams = cs._parse_family_args(specs)
            _pooled, prior, cov = cs.leave_family_out(fams, "alpha")
            self.assertEqual(cov["pooled_families"], ["beta", "gamma"])
            self.assertIsNotNone(prior)

    def test_unknown_target_raises(self):
        with tempfile.TemporaryDirectory() as d:
            self._three_families(d)
            fams = cs.discover_families(d)
            with self.assertRaises(KeyError):
                cs.leave_family_out(fams, "nonesuch")

    def test_target_concept_collision_not_dropped(self):
        """Two target concepts sharing a Concepticon id must BOTH get a prior
        row, not be silently collapsed to one."""
        with tempfile.TemporaryDirectory() as d:
            self._three_families(d)
            # gamma gets a second local concept (G_hand2) also mapping to HAND.
            g = os.path.join(d, "gamma")
            _write_forms(os.path.join(g, "forms.tsv"),
                         [("P", "G_hand", "z"), ("Q", "G_hand", "z"),
                          ("P", "G_hand2", "z"), ("Q", "G_hand2", "z"),
                          ("P", "G_stone", "1"), ("Q", "G_stone", "1")])
            _write_concepts(os.path.join(g, "concepts.tsv"),
                            [("G_hand", "HAND", "100"),
                             ("G_hand2", "HAND", "100"),   # same Concepticon id
                             ("G_stone", "STONE", "200")])
            fams = cs.discover_families(d)
            _pooled, prior, cov = cs.leave_family_out(fams, "gamma")
            hdr, rows = prior
            cids = {r[0] for r in rows}
            self.assertIn("G_hand", cids)
            self.assertIn("G_hand2", cids)     # not dropped
            self.assertEqual(cov["target_concepts"], 3)
            self.assertEqual(cov["covered"], 3)

    def test_pooling_is_pair_count_weighted(self):
        """A family sharing a concept over many language pairs outweighs one
        that shares it over a single pair (not equal weight per family)."""
        with tempfile.TemporaryDirectory() as d:
            # big family: 3 languages all cognate for HAND -> agree 3 / pairs 3
            _mk_family(d, "big",
                       [("A", "b_hand", "x"), ("B", "b_hand", "x"),
                        ("C", "b_hand", "x")],
                       [("b_hand", "HAND", "100")])
            # small family: 2 languages, different cognates -> agree 0 / pairs 1
            _mk_family(d, "small",
                       [("D", "s_hand", "1"), ("E", "s_hand", "2")],
                       [("s_hand", "HAND", "100")])
            _mk_family(d, "tgt",
                       [("P", "t_hand", "z"), ("Q", "t_hand", "z")],
                       [("t_hand", "HAND", "100")])
            fams = cs.discover_families(d)
            pooled, _prior, _cov = cs.leave_family_out(fams, "tgt")
            n_fam, n_pairs, s, _label = pooled["100"]
            self.assertEqual(n_fam, 2)
            self.assertEqual(n_pairs, 4)          # 3 + 1
            # pair-weighted: (3 agree + 0) / 4 pairs = 0.75, NOT equal-weight 0.5
            self.assertAlmostEqual(s, 0.75)


if __name__ == "__main__":
    unittest.main(verbosity=2)

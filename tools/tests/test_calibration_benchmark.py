"""
Tests for the held-out-calibration (LOCO) benchmark harness (M11).

The real experiment shells out to the C engine (slow); these tests exercise the
orchestration and scoring with an injected fake runner, so the coverage/error
arithmetic and calibration-file rewriting are verified deterministically and
fast.

Run: python3 tools/tests/test_calibration_benchmark.py
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from validation import calibration_benchmark as cb


CAL_HEADER = ["calibration_id", "target_id", "type", "p1", "p2", "p3",
              "is_hard", "source_key"]


def _rows(*specs):
    """specs: (target_id, type, p1) tuples -> calibration row dicts."""
    rows = []
    for i, (tid, typ, p1) in enumerate(specs):
        rows.append(dict(zip(CAL_HEADER,
                    ["c%d" % i, tid, typ, str(p1), ".", ".", "1", "SRC"])))
    return rows


class TestHeldoutTargets(unittest.TestCase):
    def test_only_positive_fixed_tips(self):
        rows = _rows(("Latin", "fixed_tip", 2050),
                     ("Spanish", "fixed_tip", 0),        # extant, excluded
                     ("root", "uniform_node", 4000))     # node, excluded
        self.assertEqual(cb.heldout_targets(rows), [("Latin", 2050.0)])

    def test_free_target_makes_uniform(self):
        rows = _rows(("Latin", "fixed_tip", 2050),
                     ("Greek", "fixed_tip", 2400))
        freed = cb.free_target(rows, "Latin", upper=6000)
        lat = next(r for r in freed if r["target_id"] == "Latin")
        self.assertEqual(lat["type"], "uniform_node")
        self.assertEqual(lat["p1"], "0")
        self.assertEqual(float(lat["p2"]), 6000.0)
        self.assertEqual(lat["source_key"], "HELDOUT")
        # the other calibration is untouched
        grk = next(r for r in freed if r["target_id"] == "Greek")
        self.assertEqual(grk["type"], "fixed_tip")

    def test_free_upper_derived_from_calibrations(self):
        rows = _rows(("Latin", "fixed_tip", 2050),
                     ("root", "uniform_node", 9000))
        freed = cb.free_target(rows, "Latin")   # upper auto
        lat = next(r for r in freed if r["target_id"] == "Latin")
        self.assertEqual(float(lat["p2"]), 9000.0)

    def test_free_unknown_raises(self):
        rows = _rows(("Latin", "fixed_tip", 2050))
        with self.assertRaises(ValueError):
            cb.free_target(rows, "Nonesuch")


def _fake_runner(posteriors):
    """posteriors: {node: (median, hpd_lo, hpd_hi)} -> runner returning them."""
    def run(cal_path, out_dir):
        return {n: {"median": m, "mean": m, "hpd_lo": lo, "hpd_hi": hi,
                    "rhat": 1.0, "ess": 500.0}
                for n, (m, lo, hi) in posteriors.items()}
    return run


class TestScoring(unittest.TestCase):
    def test_loco_one_scores_coverage_and_error(self):
        rows = _rows(("Latin", "fixed_tip", 2050))
        with tempfile.TemporaryDirectory() as d:
            runner = _fake_runner({"Latin": (2706.0, 1617.0, 3924.0)})
            rec = cb.loco_one(rows, "Latin", 2050.0, runner, d)
        self.assertEqual(rec["covered"], 1)           # 2050 in [1617, 3924]
        self.assertAlmostEqual(rec["signed_err"], 656.0)
        self.assertAlmostEqual(rec["abs_err"], 656.0)
        self.assertAlmostEqual(rec["width"], 2307.0)

    def test_loco_one_detects_miss(self):
        rows = _rows(("Latin", "fixed_tip", 2050))
        with tempfile.TemporaryDirectory() as d:
            runner = _fake_runner({"Latin": (900.0, 500.0, 1400.0)})  # 2050 outside
            rec = cb.loco_one(rows, "Latin", 2050.0, runner, d)
        self.assertEqual(rec["covered"], 0)
        self.assertAlmostEqual(rec["signed_err"], -1150.0)

    def test_run_benchmark_aggregates(self):
        rows = _rows(("Latin", "fixed_tip", 2050),
                     ("Greek", "fixed_tip", 2400))
        with tempfile.TemporaryDirectory() as d:
            runner = _fake_runner({
                "Latin": (2100.0, 1600.0, 2600.0),   # covered, err +50
                "Greek": (3000.0, 2700.0, 3300.0),   # NOT covered (2400<2700), err +600
            })
            records, summary = cb.run_benchmark(rows, runner, d)
        self.assertEqual(summary["n"], 2)
        self.assertAlmostEqual(summary["coverage"], 0.5)
        self.assertAlmostEqual(summary["mae"], (50 + 600) / 2.0)
        self.assertAlmostEqual(summary["rmse"],
                               ((50**2 + 600**2) / 2.0) ** 0.5)

    def test_missing_freed_node_raises(self):
        rows = _rows(("Latin", "fixed_tip", 2050))
        with tempfile.TemporaryDirectory() as d:
            runner = _fake_runner({"SomethingElse": (1.0, 0.0, 2.0)})
            with self.assertRaises(KeyError):
                cb.loco_one(rows, "Latin", 2050.0, runner, d)


class TestMultiFamily(unittest.TestCase):
    def _make_family(self, d, label, specs):
        sub = os.path.join(d, label)
        os.makedirs(sub)
        cal = os.path.join(sub, "cal.tsv")
        cb.write_calibrations(cal, CAL_HEADER, _rows(*specs))
        # tree/forms paths need not exist: the fake runner ignores them.
        return {"label": label, "tree": sub + "/t.nwk",
                "forms": sub + "/f.tsv", "calibrations": cal,
                "rate_logmean": -8.8, "rate_logsd": 0.8}

    def test_load_manifest_resolves_paths(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            fam = self._make_family(d, "ie", [("Latin", "fixed_tip", 2050)])
            man = {"defaults": {"warmup": 100},
                   "families": [{"label": "ie", "tree": "ie/t.nwk",
                                 "forms": "ie/f.tsv", "calibrations": "ie/cal.tsv"}]}
            mpath = os.path.join(d, "manifest.json")
            with open(mpath, "w") as f:
                json.dump(man, f)
            families, defaults = cb.load_manifest(mpath)
            self.assertEqual(defaults["warmup"], 100)
            self.assertTrue(families[0]["calibrations"].endswith("ie/cal.tsv"))
            self.assertTrue(os.path.isabs(families[0]["calibrations"]))

    def test_run_multi_aggregates_across_families(self):
        posteriors = {
            "Latin": (2100.0, 1600.0, 2600.0),          # covered
            "OldChinese": (5000.0, 4000.0, 6000.0),     # true 2800 -> NOT covered
            "OldTibetan": (1300.0, 900.0, 1700.0),      # covered
        }

        def factory(binary, tree, forms, settings):
            return _fake_runner(posteriors)

        with tempfile.TemporaryDirectory() as d:
            fams = [
                self._make_family(d, "romance", [("Latin", "fixed_tip", 2050)]),
                self._make_family(d, "st",
                                  [("OldChinese", "fixed_tip", 2800),
                                   ("OldTibetan", "fixed_tip", 1300)]),
            ]
            records, per_family, overall = cb.run_multi(
                fams, {}, d, runner_factory=factory)
        self.assertEqual(overall["n"], 3)
        self.assertAlmostEqual(overall["coverage"], 2.0 / 3.0)   # 2 of 3
        self.assertEqual(per_family["romance"]["n"], 1)
        self.assertAlmostEqual(per_family["romance"]["coverage"], 1.0)
        self.assertAlmostEqual(per_family["st"]["coverage"], 0.5)
        # every record is tagged with its family
        self.assertEqual({r["family"] for r in records}, {"romance", "st"})


class TestParsing(unittest.TestCase):
    def test_parse_posterior_nodes_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "posterior_nodes.tsv")
            with open(p, "w", encoding="utf-8") as f:
                f.write("node\tmean_bp\tmedian_bp\tsd\thpd90_lo_bp\thpd90_hi_bp\t"
                        "rhat\tess\n")
                f.write("Latin\t2702.2\t2706.3\t717.5\t1616.9\t3923.6\t1.006\t435\n")
            post = cb.parse_posterior_nodes(p)
        self.assertIn("Latin", post)
        self.assertAlmostEqual(post["Latin"]["median"], 2706.3)
        self.assertAlmostEqual(post["Latin"]["hpd_lo"], 1616.9)
        self.assertAlmostEqual(post["Latin"]["ess"], 435.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)

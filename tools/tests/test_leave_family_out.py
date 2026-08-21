"""
Tests for the leave-family-out cross-family coverage experiment (M11).

These check a genuine, falsifiable coverage experiment (not a construction that
cannot fail): the in-sample control must be calibrated, the honest pooled
leave-one-out predictive must be consistent (coverage → nominal as families
grow), and a deliberately mis-built predictive must undercover — the teeth.

Run: python3 tools/tests/test_leave_family_out.py
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from validation import leave_family_out as lfo

LEVELS = (0.5, 0.8, 0.9, 0.95)
OM, USD = lfo.family_omegas()


class TestLeaveFamilyOut(unittest.TestCase):
    def test_family_specific_control_is_calibrated(self):
        """The in-sample family-specific interval is calibrated by construction;
        if it were not, the whole harness would be suspect."""
        cov = lfo.run_coverage(0.4, OM, USD, LEVELS, n_reps=3000)
        for lv in LEVELS:
            self.assertAlmostEqual(cov["family_specific"][lv], lv, delta=0.03)

    def test_pooled_loo_is_consistent(self):
        """The honest pooled leave-one-out predictive is consistent: coverage
        climbs toward nominal as the number of families grows. At few families
        it modestly undercovers (documented small-sample effect)."""
        c10 = lfo.run_coverage(0.4, OM, USD, (0.9,), n_reps=2000, n_families=10)
        c150 = lfo.run_coverage(0.4, OM, USD, (0.9,), n_reps=2000, n_families=150)
        self.assertLess(c10["pooled_loo"][0.9], c150["pooled_loo"][0.9])
        self.assertAlmostEqual(c150["pooled_loo"][0.9], 0.9, delta=0.03)
        # small-sample undercoverage is real but bounded, not catastrophic
        self.assertGreater(c10["pooled_loo"][0.9], 0.75)

    def test_broken_predictive_undercovers(self):
        """Teeth: a predictive that drops the between-family variance τ̂²
        undercovers badly, so the experiment can tell a correct predictive from
        a wrong one."""
        good = lfo.run_coverage(0.4, OM, USD, (0.9,), n_reps=3000)
        bad = lfo.run_coverage(0.4, OM, USD, (0.9,), n_reps=3000, broken=True)
        self.assertLess(bad["pooled_loo"][0.9], good["pooled_loo"][0.9] - 0.15)
        self.assertLess(bad["pooled_loo"][0.9], 0.65)

    def test_pooled_tracks_nominal_better_than_flat_fallback(self):
        """At τ_true below the flat UNKNOWN sd, the fixed fallback over-covers
        (conservative); the pooled predictive sits closer to nominal at 0.90."""
        cov = lfo.run_coverage(0.4, OM, USD, (0.9,), n_reps=3000)
        self.assertGreater(cov["unknown_fallback"][0.9], 0.95)   # conservative
        pooled_gap = abs(cov["pooled_loo"][0.9] - 0.9)
        flat_gap = abs(cov["unknown_fallback"][0.9] - 0.9)
        self.assertLess(pooled_gap, flat_gap)

    def test_predictive_needs_two_families(self):
        with self.assertRaises(ValueError):
            lfo.pooled_loo_predictive([(0.1, 0.3)])

    def test_evaluate_report_shape(self):
        rep = lfo.evaluate(tau_true=0.4, n_reps=2000)
        self.assertEqual(rep["n_families"], len(OM))
        self.assertIn("pooled_loo", rep["coverage"])
        self.assertGreaterEqual(rep["pooled_loo_max_gap"], 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)

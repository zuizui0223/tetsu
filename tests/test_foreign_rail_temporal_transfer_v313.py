"""Behavioral guards for retrospective spatial-block transition prediction."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import analyze_foreign_rail_temporal_transfer_v313 as m
except ModuleNotFoundError as exc:
    if exc.name != "analyze_foreign_rail_temporal_transfer_v313":
        raise
    m = None


def fixture():
    rng = np.random.default_rng(313)
    species = ["Alpha alba", "Beta beta", "Gamma gamma", "Delta delta", None]
    plots, meta = {}, {}
    for i in range(15):
        s = f"site{i:02}"
        meta[s] = {"lat": 30 + i/100, "altitude": 100+i, "station": i % 2 == 0}
        for z in ("Track", "Interior Plot 1", "Interior Plot 2", "Embankment"):
            a = rng.integers(0, 3, 5).astype(float)
            a[2] = 0  # species only detected after baseline must never enter eligibility
            plots[s, 2014, z] = a
            plots[s, 2017, z] = rng.integers(0, 2, 5).astype(float)
    donors = [
        {"Alpha alba": .2, "Beta beta": .7, "Gamma gamma": .5, "Delta delta": .2},
        {"Alpha alba": .6, "Beta beta": .5, "Gamma gamma": .5},
        {"Alpha alba": .3, "Beta beta": .4},
        {"Alpha alba": .5, "Beta beta": .8},
    ]
    return species, plots, meta, donors


class TemporalTransfer(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m, "v313 implementation not yet present")

    def test_baseline_eligibility_and_whole_transect_missing_exclusion(self):
        species, plots, meta, donors = fixture()
        plots["site00", 2017, "Interior Plot 2"][:] = np.nan
        c = m.prepare_data(species, plots, meta, donors)
        self.assertEqual(c["taxa"], ["Alpha alba", "Beta beta"])
        self.assertEqual(len(c["sites"]), 14)
        self.assertNotIn("site00", c["sites"])
        self.assertEqual(np.bincount(c["blocks"]).tolist(), [3,3,3,3,2])
        altered = deepcopy(plots)
        for key in altered:
            if key[1] == 2017 and np.isfinite(altered[key]).all():
                altered[key] = 1 - (altered[key] > 0).astype(float)
        other = m.prepare_data(species, altered, meta, donors)
        self.assertEqual(c["taxa"], other["taxa"])
        np.testing.assert_array_equal(c["abundance0"], other["abundance0"])
        np.testing.assert_array_equal(c["blocks"], other["blocks"])

    def test_baseline_frequency_excludes_own_transect_and_heldout(self):
        a = np.zeros((3, 3, 1), float)
        a[0, 0, 0] = 1
        a[2, 0, 0] = 1
        f = m.baseline_frequencies(a, np.array([True, True, False]))
        np.testing.assert_allclose(f[:,0,0], [0,1,.5])
        a[2,:,:] = 100  # heldout baseline is not part of training frequency
        np.testing.assert_allclose(m.baseline_frequencies(a, np.array([True,True,False])), f)

    def test_future_outcomes_cannot_change_features_or_heldout_predictions(self):
        c = m.prepare_data(*fixture())
        changed = deepcopy(c)
        changed["presence1"][c["blocks"] == 0] = 1 - changed["presence1"][c["blocks"] == 0]
        for endpoint in ("persistence", "new_record"):
            original = m.forecast_fold(c, endpoint, 0, .01)
            altered = m.forecast_fold(changed, endpoint, 0, .01)
            for model in ("M0","M1","M1b","M2"):
                np.testing.assert_allclose(original["predictions"][model], altered["predictions"][model], atol=1e-12)
            x = m.make_rows(c, endpoint, c["blocks"] != 0)
            all_future_changed = deepcopy(c)
            all_future_changed["presence1"] = 1 - c["presence1"]
            xx = m.make_rows(all_future_changed, endpoint, c["blocks"] != 0)
            np.testing.assert_array_equal(x["x"], xx["x"])
            np.testing.assert_array_equal(x["index"], xx["index"])

    def test_weight_habitat_means_instead_of_unequal_risk_set_sizes(self):
        weights = m.habitat_weights(np.array([0,0,1,2,2,2]))
        np.testing.assert_allclose([weights[:2].sum(), weights[2], weights[3:].sum()], [.5,.25,.25])
        self.assertAlmostEqual(m.balanced_average([.2,.8,.4]), .4)
        self.assertIsNone(m.balanced_average([.2,None,.4]))
        # A block with zero positive events still has a valid loss.
        result = m.score_habitats(np.array([0,0,0]), np.array([.1,.1,.1]), np.array([0,1,2]))
        self.assertAlmostEqual(result["balanced"]["log_loss"], -np.log(.9))
        self.assertEqual(result["by_habitat"]["Track"]["events"], 0)

    def test_risk_sets_use_baseline_and_together_exhaust_eligible_rows(self):
        c = m.prepare_data(*fixture())
        a = m.make_rows(c, "persistence", c["blocks"] != 0)
        b = m.make_rows(c, "new_record", c["blocks"] != 0)
        self.assertEqual(len(a["y"]) + len(b["y"]), 15*3*2)
        self.assertFalse(set(map(tuple, a["index"])) & set(map(tuple, b["index"])))
        for row in a["index"]:
            self.assertGreater(c["abundance0"][tuple(row)], 0)
        for row in b["index"]:
            self.assertEqual(c["abundance0"][tuple(row)], 0)

    def test_unstandardized_foreign_slopes_keep_both_interior_directions(self):
        s = m.foreign_slopes(np.array([0,.4,.9,1.6]), np.array([2.,3.,4.]))
        self.assertAlmostEqual(s["Track"], .6)
        self.assertAlmostEqual(s["Interior Plot 1"], .5)
        self.assertAlmostEqual(s["Interior Plot 2"], -.1)
        self.assertAlmostEqual(s["Track_minus_interior_mean"], .4)
        self.assertAlmostEqual(s["Track_minus_Interior_Plot_1"], .1)
        self.assertAlmostEqual(s["Track_minus_Interior_Plot_2"], .7)


if __name__ == "__main__":
    unittest.main()


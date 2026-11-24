"""Synthetic regression checks for v318 original-source reanalysis logic."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_Slovenia_thermal_EIV_paired_v318 as m

class ThermalTests(unittest.TestCase):
    def test_unknown_temperature_is_missing_not_zero(self):
        traits = {"A":{"Temperature":6}, "B":{"Temperature":"x"}}
        self.assertEqual(m.temperature_cwm({"A","B"},traits,"Temperature"),
                         (6,0.5,1,2))
    def test_two_track_sides_one_independent_location(self):
        p = [dict(site="s",side=s,segment="a",shoulder=7.,toe=6.,
                  shoulder_coverage=1.,toe_coverage=1.) for s in ("L","R")]
        o=m.summarize(p,0)
        self.assertEqual(o["eligible_track_side_pairs"],2)
        self.assertEqual(o["eligible_distinct_locations"],1)
        self.assertEqual(o["mean_BS_minus_BT_EIV_T"],1)
    def test_coverage_gate_discards_ineligible_pair(self):
        p=[dict(site="s",segment="a",shoulder=7.,toe=6.,
                shoulder_coverage=.5,toe_coverage=.7)]
        self.assertEqual(m.summarize(p,.5)["eligible_distinct_locations"],1)
        self.assertEqual(m.summarize(p,.7)["eligible_distinct_locations"],0)
    def test_source_signature_identifies_only_temperature(self):
        tr={f"x{i}":{"Temperature":(5 if i<142 else ("x" if i<277 else None))}
            for i in range(312)}
        self.assertEqual(m.find_temperature_column(tr)[0],"Temperature")
    def test_unresolvable_trait_column_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"Temperature column ambiguous"):
            m.find_temperature_column({"A":{"SpeciesName":"A","LifeForm":"H"}})
if __name__ == "__main__":
    unittest.main()

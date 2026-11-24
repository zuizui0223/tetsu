"""v319 synthetic checks. No original field data are supplied here."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Slovenia_thermal_generalist_section_v319 as m

def fixture():
    traits={"A":{"EIV_T":7},"B":{"EIV_T":5},
            "X":{"EIV_T":"x"},"M":{"EIV_T":None}}
    pairs=[{"site":"s1","section":"1","side":"L",
            "BS":{"A","X"},"BT":{"B","M"}},
           {"site":"s1","section":"1","side":"R",
            "BS":{"A"},"BT":{"B"}},
           {"site":"s2","section":"2","side":"L",
            "BS":{"B"},"BT":{"A","X"}}]
    return pairs,traits

class SourceMethods(unittest.TestCase):
    def test_x_classified_generalist_not_temperature_zero(self):
        self.assertIsNone(m.eiv("x"))
        self.assertEqual(m.trait_class("x"),"generalist_x")
        self.assertEqual(m.trait_class(None),"blank")
        self.assertEqual(m.trait_class(7),"numeric")

    def test_coverage_and_generalist_fraction(self):
        pairs,traits=fixture()
        z=m.zone_measure(pairs[0]["BS"],traits)
        self.assertEqual(z["coverage"],0.5)
        self.assertEqual(z["x"],0.5)
        self.assertEqual(z["value"],7.0)

    def test_field_location_not_two_independent_sides(self):
        pairs,traits=fixture()
        s=m.field_summary(pairs,traits,0)
        self.assertEqual(s["n_sites"],2)
        self.assertEqual(s["n_sides"],3)
        # site1 mean +2; site2 -2 => equal-site mean zero
        self.assertAlmostEqual(s["mean_eiv_difference"],0)

    def test_coverage_filters_and_single_taxon_omit(self):
        pairs,traits=fixture()
        s=m.field_summary(pairs,traits,0.9)
        self.assertEqual(s["n_sites"],1)
        s2=m.field_summary(pairs,traits,0,excluded="A")
        self.assertEqual(s2["n_sites"],0)

    def test_exact_section_signflip(self):
        t=m.section_test({"s1":[1,1],"s2":[2],"s3":[3]})
        self.assertEqual(t["enumerated_section_sign_flips"],8)
        self.assertAlmostEqual(t["signed_flip_two_sided_fraction"],.25)
        self.assertEqual(t["positive_section_count"],3)

    def test_section_drop(self):
        t=m.section_drop({"a":[1,3],"b":[2],"c":[4]})
        self.assertTrue(t["all_leave_one_section_out_positive"])

if __name__=="__main__":unittest.main()

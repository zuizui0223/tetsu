"""v25.6 synthetic mathematical QA, no empirical data in this file."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_rail_microhabitat_replacement_nestedness_v256 as m


class ReplacementNestednessTests(unittest.TestCase):
    def test_pure_nestedness_and_rail_subset(self):
        x=m.components({1,2},{1,2,3,4})
        self.assertEqual(x["exclusive_rail_taxa"],0)
        self.assertEqual(x["beta_replacement_jtu"],0)
        self.assertAlmostEqual(x["beta_jac"],.5)
        self.assertAlmostEqual(x["beta_nestedness_jne"],.5)
        self.assertEqual(x["strict_rail_subset"],1.)

    def test_pure_replacement_without_shared_species(self):
        x=m.components({1,2},{3,4})
        self.assertAlmostEqual(x["beta_jac"],1.)
        self.assertAlmostEqual(x["beta_replacement_jtu"],1.)
        self.assertAlmostEqual(x["beta_nestedness_jne"],0.)
        self.assertEqual(x["replacement_fraction_of_beta"],1.)
        self.assertEqual(x["strict_rail_subset"],0.)

    def test_mixed_replacement_and_nestedness(self):
        x=m.components({1,2,3},{2,3,4,5})
        self.assertAlmostEqual(x["beta_jac"],.6)
        self.assertAlmostEqual(x["beta_replacement_jtu"],.5)
        self.assertAlmostEqual(x["beta_nestedness_jne"],.1)
        self.assertAlmostEqual(x["replacement_fraction_of_beta"],5/6)
        self.assertEqual(x["exclusive_rail_taxa"],1)

    def test_identical_sets_fraction_of_zero_beta_is_undefined(self):
        x=m.components({1,2},{1,2})
        self.assertEqual(x["beta_jac"],0)
        self.assertIsNone(x["replacement_fraction_of_beta"])
        self.assertEqual(x["strict_rail_subset"],1)
        self.assertIsNone(m.average([x])["replacement_fraction_of_beta"])

    def test_undefined_empty_zone_fails_closed(self):
        with self.assertRaises(ValueError):
            m.components(set(),{1})

    def test_combined_and_year_decomposition_under_strict_common_site(self):
        raw={
            "station":{
                "A":{
                    2014:{"Embankment":m.components({1,2},{1,2,3}),"Interior Plot 1":m.components({1,2},{1,2,3}),"Interior Plot 2":m.components({1,2},{2,3})},
                    2017:{"Embankment":m.components({1,2},{1,2,3}),"Interior Plot 1":m.components({1,2},{1,2,3}),"Interior Plot 2":m.components({1,2},{2,3})}
                },
                "B":{
                    2014:{"Embankment":m.components({1,2},{1,2,3}),"Interior Plot 1":m.components({1,2},{1,2,3}),"Interior Plot 2":m.components({1,2},{2,3})},
                    2017:{"Embankment":m.components({1,2},{1,2,3}),"Interior Plot 1":m.components({1,2},{1,2,3}),"Interior Plot 2":m.components({1,2},{2,3})}
                }
            },
            "between":{
                "C":{
                    2014:{"Embankment":m.components({1,2},{1,2}),"Interior Plot 1":m.components({1,2},{1,2}),"Interior Plot 2":m.components({1,2},{1,2})},
                    2017:{"Embankment":m.components({1,2},{1,2}),"Interior Plot 1":m.components({1,2},{1,2}),"Interior Plot 2":m.components({1,2},{1,2})}
                },
                "D":{
                    2014:{"Embankment":m.components({1,2},{1,2}),"Interior Plot 1":m.components({1,2},{1,2}),"Interior Plot 2":m.components({1,2},{1,2})},
                    2017:{"Embankment":m.components({1,2},{1,2}),"Interior Plot 1":m.components({1,2},{1,2}),"Interior Plot 2":m.components({1,2},{1,2})}
                }
            }
        }
        q=m.india_aggregate(raw,101)
        self.assertAlmostEqual(q["year_by_reference_station_minus_between"]["2014:Interior Plot 2"]["station_minus_between"]["beta_jac"],2/3)
        self.assertAlmostEqual(q["paired_site_year_change_in_station_minus_between"]["Interior Plot 2"]["observed"]["beta_jac"],0)
        self.assertEqual(q["year_by_reference_station_minus_between"]["2014:Interior Plot 2"]["bootstrap"]["beta_jac"]["draws"],101)

    def test_source_section_mean_not_side_pseudoreplication(self):
        s={"section1":[m.components({1,2},{1,2,3}) for i in range(10)],
           "section2":[m.components({1,2},{3,4})]}
        o=m.slovenia_aggregate(s,101)
        self.assertAlmostEqual(o["section_equal_observed"]["beta_jac"],(.3333333333333333+1)/2)
        self.assertEqual(o["bootstrap_source_section_count"],2)
        self.assertEqual(o["section_bootstrap"]["beta_jac"]["draws"],101)


if __name__=="__main__":
    unittest.main()

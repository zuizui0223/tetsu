"""v10.9 numerical source-taxonomy robustness diagnostics on synthetic pairs.

Real railway botanical data are downloaded separately by the source-backed
workflow. The unit-test toy assemblages do NOT constitute biological evidence.
"""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    source_guild,taxa_preference,contrast,taxon_label_permutation,run
)


def r(taxon, wind, annual, val, support=4):
    return {"taxon_code":taxon,"wind":wind,
            "stratum":"strict_annual" if annual else "other_known_lifespan",
            "n_original_locations_present_either_zone":support,
            "occupancy_preference_BS_minus_BT":val}


class OriginalRailwaySpeciesPreferenceTests(unittest.TestCase):
    def test_original_trait_wind_and_lifespan_coding_not_missing_nonwind(self):
        self.assertIsNone(source_guild({"Diaspore.strategy":None,"Lifespan.categories":"annual"}))
        self.assertIsNone(source_guild({"Diaspore.strategy":"Phragmites","Lifespan.categories":None}))
        self.assertEqual(source_guild({"Diaspore.strategy":"Epilobium",
                                        "Lifespan.categories":"annual"}),
                         {"wind":True,"stratum":"strict_annual"})
        self.assertEqual(source_guild({"Diaspore.strategy":"Lycopodium",
                                        "Lifespan.categories":"perennial"}),
                         {"wind":False,"stratum":"other_known_lifespan"})

    def test_species_occupancy_pooling_across_actual_location_not_track_side(self):
        sites={
            "L1":{"section_number":"1","BS":{"W","N"},"BT":{"N"},"n_matched_sides":2},
            "L2":{"section_number":"2","BS":set(),"BT":{"W","N"},"n_matched_sides":1},
            "L3":{"section_number":"3","BS":{"W"},"BT":set(),"n_matched_sides":1},
        }
        traits={
            "W":{"Diaspore.strategy":"Phragmites","Lifespan.categories":"annual"},
            "N":{"Diaspore.strategy":"Bidens","Lifespan.categories":"annual"},
        }
        observed,n=taxa_preference(sites,traits)
        self.assertEqual(n,3)
        wx=next(x for x in observed if x["taxon_code"]=="W")
        nx=next(x for x in observed if x["taxon_code"]=="N")
        self.assertEqual((wx["BS_only"],wx["BT_only"],wx["both"]),(2,1,0))
        self.assertAlmostEqual(wx["occupancy_preference_BS_minus_BT"],1/3)
        self.assertEqual((nx["BS_only"],nx["BT_only"],nx["both"]),(0,1,1))
        self.assertAlmostEqual(nx["occupancy_preference_BS_minus_BT"],-1/2)
        keep,n=taxa_preference(sites,traits,exclude_section="2")
        self.assertEqual(n,2)
        self.assertAlmostEqual(next(x for x in keep if x["taxon_code"]=="W")[
            "occupancy_preference_BS_minus_BT"],1)

    def test_species_equal_stratified_negative_wind_preference(self):
        d=[]
        for annual in (True,False):
            for i in range(3):
                d.append(r("W"+str(annual)+str(i),True,annual,-1))
                d.append(r("N"+str(annual)+str(i),False,annual,+1))
        c=contrast(d,1)
        self.assertTrue(c["identified"])
        self.assertEqual(c["n_qualifying_species"],12)
        self.assertAlmostEqual(c["equal_species_within_lifespan_wind_minus_other_preference"],-2)
        p=taxon_label_permutation(d,1,draws=1000,seed=42)
        self.assertLess(p["empirical_two_sided_permutation_p"],.01)
        self.assertLess(p["empirical_one_sided_wind_more_BT_p"],.01)

    def test_taxon_rare_cutoff_prevents_result_selected_on_single_occurrences(self):
        d=[]
        for annual in (True,False):
            for i in range(3):
                d.append(r("W"+str(annual)+str(i),True,annual,-1,support=1))
                d.append(r("N"+str(annual)+str(i),False,annual,+1,support=5))
        self.assertTrue(contrast(d,1)["identified"])
        self.assertFalse(contrast(d,3)["identified"])
        self.assertFalse(taxon_label_permutation(d,10,draws=100)["identified"])

    def test_species_omission_changes_result_when_one_taxon_carries_weak_signal(self):
        d=[]
        for annual in (True,False):
            for i in range(3):
                d.append(r("W"+str(annual)+str(i),True,annual,-.5 if i==0 else .5))
                d.append(r("N"+str(annual)+str(i),False,annual,0))
        orig=contrast(d,1)
        self.assertGreater(orig["equal_species_within_lifespan_wind_minus_other_preference"],0)
        omit=contrast(d,1,excluded_taxon="WTrue2")
        self.assertLess(omit["equal_species_within_lifespan_wind_minus_other_preference"],
                        orig["equal_species_within_lifespan_wind_minus_other_preference"])

    def test_requires_actual_locked_v106_receipt(self):
        with self.assertRaisesRegex(ValueError,"receipt missing"):
            run({},{"complete":False})


if __name__=="__main__":unittest.main()

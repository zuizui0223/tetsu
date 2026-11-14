"""v10.4: source original field species-level Diaspore.strategy class
math and same-location track-side source clustering tests. WIND
dispersal trait != actual train gust, observed seeds or adaptation.
"""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_wind_diaspore_vs_annual_v104 import (
    WIND_STRA,WIND_SENS,GROUPS,VALID_SOURCE_CATEGORIES,
    strategy_fraction,location_level,summarize_locations,joint_direction
)
def srcpair(lid,side,section=10,wind=0.2,annual=0.3):
    return {
        "location_id":lid,"track_side":side,
        "section_number":section,"track_type":"SINGLE","survey_year":2022,
        "source_cover_BS_minus_BT":-40.0,
        "original_source_cover_GT100":False,
        "source_strict_annual_BS_minus_BT":annual,
        "source_strict_therophyte_BS_minus_BT":annual,
        "wind_and_diaspore_group_shares_BS_minus_BT":{
            WIND_STRA:wind,
            WIND_SENS:wind,
            "source_Allium_autochory":0.05,
            "source_Bidens_epizoochory":-0.02,
            "source_Cornus_endozoochory":0.
        },
    }

class TrueRailDiasporeSourceTests(unittest.TestCase):
    def test_original_source_classes_are_explicit_and_not_train_wind(self):
        self.assertEqual(GROUPS[WIND_STRA],frozenset(("Epilobium","Phragmites")))
        self.assertIn("Lycopodium",GROUPS[WIND_SENS])
        self.assertIn("Allium",VALID_SOURCE_CATEGORIES)
    def test_missing_trait_is_not_auto_nonwind(self):
        trait={
         "A":{"Diaspore.strategy":"Epilobium"},
         "B":{"Diaspore.strategy":"Phragmites"},
         "C":{"Diaspore.strategy":"Allium"},
         "D":{"Diaspore.strategy":None},
        }
        v=strategy_fraction(set(trait),trait,GROUPS[WIND_STRA])
        self.assertAlmostEqual(v["fraction_of_source_trait_known_taxa"],2/3)
        self.assertEqual(v["n_taxa_missing_original_Diaspore_strategy"],1)
        with self.assertRaisesRegex(ValueError,"declared"):
            strategy_fraction(set(trait),trait,{"TRAIN_WIND_PROOF"})
        with self.assertRaisesRegex(ValueError,"lacking"):
            strategy_fraction({"A","MISSING"},trait,GROUPS[WIND_STRA])
    def test_same_location_two_sides_first_before_section_aggregation(self):
        pairs=[]
        for i in range(1,90):
            pairs.append(srcpair("loc"+str(i),"L",section=i%11,wind=-0.05,annual=0.1))
            if i<85:pairs.append(srcpair("loc"+str(i),"R",section=i%11,
                                       wind=-0.03,annual=0.05))
        source=location_level(pairs)
        self.assertEqual(len(source),89)
        self.assertEqual(source[0]["n_original_BS_BT_source_track_side_pairs"],2)
        wind=summarize_locations(source,WIND_STRA+"_delta")
        self.assertEqual(wind["n_independent_actual_RAIL_field_locations"],89)
        self.assertEqual(wind["n_locations_BT_larger"],89)
        self.assertEqual(wind["n_source_RAIL_segments_with_3plus_independent_locations"],11)
    def test_annual_up_wind_down_categories_not_selection_causality(self):
        rows=[srcpair("loc"+str(i),"L",section=i%15,
                      wind=-.05 if i<50 else .03,
                      annual=.12 if i<70 else -.03)
              for i in range(1,90)]
        source=location_level(rows)
        x=joint_direction(source)
        self.assertEqual(x["n_original_locations_with_valid_both_trait_source_differences"],89)
        self.assertGreater(x["contrasting_trait_directions"]["BS_annual_up_BS_wind_down"],0)
        self.assertTrue(x["source_species_lifespan_does_NOT_measure_selection_or_gust_causality"])
    def test_invalid_region_location_join_not_extra_replicates(self):
        pairs=[srcpair("S1","L"),srcpair("S1","R",section=21)]
        with self.assertRaisesRegex(ValueError,"not nested"):
            location_level(pairs)

if __name__=="__main__":unittest.main()

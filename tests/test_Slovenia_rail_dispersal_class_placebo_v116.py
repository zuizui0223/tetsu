"""v11.6 genuine source class support, no selective favorable type cherry-pick."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_dispersal_class_placebo_v116 import (
    SOURCE_CATEGORIES,COMPARISONS,label_strategy,support,run
)


class RealSourceStrategyPlaceboTests(unittest.TestCase):
    def test_every_original_diaspore_strategy_predefined(self):
        self.assertEqual(len(SOURCE_CATEGORIES),8)
        self.assertEqual(len(COMPARISONS),10)
        self.assertEqual(set(COMPARISONS["strict_wind_Epilobium_Phragmites"]),
                         {"Epilobium","Phragmites"})
        self.assertEqual(set(COMPARISONS["broad_wind_plus_Lycopodium_spore"]),
                         {"Epilobium","Phragmites","Lycopodium"})

    def test_all_known_guilds_keep_same_species_and_lifespan(self):
        codes=["E","P","A","Z"]
        rows=[{"taxon_code":c,"stratum":"strict_annual","wind":False,
               "n_original_locations_present_either_zone":2,
               "occupancy_preference_BS_minus_BT":.3,
               "site_delta":{"l1":.3}} for c in codes]
        original={
            "E":{"Diaspore.strategy":"Epilobium"},
            "P":{"Diaspore.strategy":"Phragmites"},
            "A":{"Diaspore.strategy":"Allium"},
            "Z":{"Diaspore.strategy":"Zea"}}
        mapped=label_strategy(rows,original,("Epilobium","Phragmites"))
        self.assertEqual([r["wind"] for r in mapped],[True,True,False,False])
        self.assertTrue(all(r["stratum"]=="strict_annual" for r in mapped))
        self.assertEqual(len(mapped),len(rows))

    def test_missing_strategy_fails_not_nonwind(self):
        with self.assertRaisesRegex(ValueError,"Missing or unexpected"):
            label_strategy([{"taxon_code":"X"}],{"X":{"Diaspore.strategy":None}},("Epilobium",))
        with self.assertRaises(ValueError):
            label_strategy([],{},("non_source_strategy",))

    def test_unsupported_source_lifespan_cell_is_marked_not_estimated(self):
        r=[]
        for grp in ("strict_annual","other_known_lifespan"):
            for i in range(3):
                r.append({"stratum":grp,"wind":i==0})
        counts,identified=support(r)
        self.assertFalse(identified)
        self.assertEqual(counts["strict_annual"]["focus"],1)

    def test_original_source_receipt_required(self):
        with self.assertRaisesRegex(ValueError,"authenticated original v11.5"):
            run({},{"complete":False})


if __name__=="__main__":
    unittest.main()

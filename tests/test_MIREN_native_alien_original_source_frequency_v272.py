"""v27.2 mathematical QA of regional Native/Alien and exact source margins."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_native_alien_road_local_reference_frequency_null_v272 as v


class MIRENNativeAlienTests(unittest.TestCase):
    def test_source_conflicting_regional_status_unknown_no_favorable_impute(self):
        raw={("R1","TaxonA"):{"Alien"},("R1","TaxonB"):{"Native"},
             ("R1","TaxonC"):{"Alien","Native"},
             ("R2","TaxonA"):{"Native"}}
        result,qa=v.resolve_regional_status(raw)
        self.assertEqual(result["R1","TaxonA"],"Alien")
        self.assertEqual(result["R2","TaxonA"],"Native")
        self.assertEqual(result["R1","TaxonC"],"Unknown")
        self.assertEqual(qa["conflicting_original_region_taxon_status_conservatively_UNKNOWN"],1)
        self.assertEqual(qa["unknown_region_taxon_pairs"],1)

    def test_positive_both_status_original_site_gates(self):
        original={("R1","A"):[],("R1","B"):[],
                  ("R2","A"):[],("R2","B"):[]}
        for (region,road),rows in original.items():
            for site in range(4):
                rows.append({"site":str(site),
                             "corridor":{1,site+10},
                             "reference":{1,site+20}})
        supported=v.balance_by_roads(original)
        self.assertEqual(len(supported),2)
        self.assertEqual(sum(len(q) for q in supported.values()),4)
        original["R1","B"]=original["R1","B"][:3]
        supported=v.balance_by_roads(original)
        self.assertEqual(set(supported),{"R2"})
        self.assertEqual(len(supported["R2"]),2)

    def test_native_alien_source_macro_fixed_margin_positive_support(self):
        original={}
        for reg in ("R1","R2"):
            for road in ("A","B"):
                original[(reg,road)]=[
                    {"site":str(i),"corridor":{1,i+10},
                     "reference":{1,i+20}} for i in range(5)]
        out=v.analyze_status_panel(original,draws=19,burnin=100,spacing=20,seed=19)
        self.assertTrue(out["identified"])
        self.assertEqual(out["n_source_regions"],2)
        self.assertEqual(out["n_original_roads"],4)
        self.assertEqual(out["n_original_matched_source_sites"],20)
        self.assertIn("beta_jac",out["equal_road_then_equal_region_results"])
        self.assertLessEqual(out["equal_road_then_equal_region_results"]["beta_jac"]["mean_source_region_observed"],1.)

    def test_alien_fraction_within_same_original_sites_and_roads(self):
        records=[]
        for reg in ("R1","R2"):
            for road in ("A","B"):
                for i in range(4):
                    records.append({"region":reg,"road":road,
                       "corridor_only_n":5,"reference_only_n":5,
                       "known_corridor_only_n":5,"known_reference_only_n":5,
                       "alien_share_of_known_corridor_only":.8,
                       "alien_share_of_known_reference_only":.2})
        out=v.alien_fraction_comparison(records)
        self.assertEqual(out["source_region_codes"],2)
        self.assertEqual(out["source_roads"],4)
        self.assertEqual(out["source_sites"],16)
        self.assertAlmostEqual(out["mean_region_equal_alien_fraction_diff"],.6)
        self.assertEqual(out["regions_corridor_only_more_alien"],2)
        # Unknown taxa produce undefined proportions, never zero imputations
        records[0]["alien_share_of_known_reference_only"]=None
        out2=v.alien_fraction_comparison(records)
        self.assertEqual(out2["source_sites"],8)  # R1 loses one entire road and thus Region support

    def test_source_zero_alien_plot_not_falsely_forced_into_panel(self):
        only_one_zone={("R1","road1"):[{"site":str(i),"corridor":{1},
                                              "reference":{1,2}} for i in range(4)]}
        self.assertFalse(v.analyze_status_panel(only_one_zone,draws=19,
                                                burnin=20,spacing=10)["identified"])


if __name__=="__main__":unittest.main()

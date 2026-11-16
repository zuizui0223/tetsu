"""v12.0: real source regional heterogeneity, source provenance and noncausal guards."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_rail_road_region_robustness_v120 import (
    validated_source,conditional_regions,summary,run
)


def mk(code,eq,occ,nr=10,nd=10,identified=True):
    return {"source_phytogeographic_region":code,"source_equal":eq,
            "source_occupancy":occ,"nrail":nr,"nroad":nd,"identified":identified}


class SourceCzechPhytoRegionalTests(unittest.TestCase):
    def test_original_source_receipt_is_required(self):
        with self.assertRaisesRegex(ValueError,"not present"):
            run({"version":"test","complete":False})

    def test_region_signs_and_actual_geographic_reweight_math(self):
        rr=[mk("A",.1,.2,5,20),mk("B",-.2,-.1,20,5),
            mk("C",.1,-.05,15,15)]
        x=summary(rr,draws=100,seed=123)
        self.assertEqual(x["n_original_source_eligible_phyto_regions"],3)
        self.assertEqual(x["n_original_regions_any_opposite_sign"],1)
        self.assertEqual(x["n_original_regions_positive_equal"],2)
        self.assertEqual(x["n_original_regions_positive_occupancy_weight"],1)
        self.assertAlmostEqual(x["original_region_equal_mean_species_equal"],0)
        self.assertAlmostEqual(x["original_region_equal_mean_occupancy_weighted"],
                               (.2-.1-.05)/3)
        self.assertEqual(x["n_leave_one_source_region_out"],3)
        self.assertEqual(x["n_region_resample_draws"],100)

    def test_min_source_sites_and_guild_support_gates(self):
        panel={"source_region_comparisons_not_pairwise_causal_rail_vs_road":[
            mk("A",.1,-.1,4,20),mk("B",.2,.1,8,8),
            mk("C",.2,-.2,20,15),mk("D",None,None,50,50,False)]}
        self.assertEqual(len(conditional_regions(panel,5,5)),2)
        self.assertEqual(len(conditional_regions(panel,10,10)),1)
        self.assertEqual(len(conditional_regions(panel,20,20)),0)

    def test_region_count_duplicate_fails(self):
        with self.assertRaisesRegex(ValueError,"repeated"):
            summary([mk("A",.2,.1),mk("A",.1,.2)],draws=20)

    def test_wrong_zenodo_provenance_fails_even_when_results_appear_successful(self):
        src={"complete":True,"version":"v11.9",
             "original_source_Czechveg_header_provenance":{"checksum":"FORGED"},
             "original_source_Czechveg_species_provenance":{"checksum":"FORGED"},
             "original_source_Pladias_trait_quality":{
                 "source_Pladias_original_diaspore_rows":3147,
                 "source_Pladias_original_lifeform_rows":4599,
             }}
        with self.assertRaisesRegex(ValueError,"MD5 mismatch"):
            validated_source(src)


if __name__=="__main__":unittest.main()

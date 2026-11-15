"""v16.0 synthetic-only source arithmetic and ecology-safe gate tests.

These fixtures are NOT original Czech field measurements; verified
published PBF and 39 station source fields are read only in Actions live run.
"""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_Czech_2024_centroid_geographic_identity_v154 import fake_source_sites
from test_Czech_2023_real_road_rail_heldout_identity_v156 import synthetic_pairs
from analyze_Czech_2023_true_rail_road_heldout_taxon_identity_v156 import (
    pair_records,cv_section_heldout
)
from analyze_Czech_arch_neo_true_railroad_5stretch_CV_v160 import (
    reassign_original_source_class_outcome,assert_original_reference,
    RICHNESS_MODES,N_STATIONS
)

class OriginalExplicitAlienSourceRouteTests(unittest.TestCase):
    def setUp(self):
        self.sites,self.guild=fake_source_sites()
        pairs=synthetic_pairs(self.sites)
        self.baseline,_=pair_records(pairs,self.sites,self.guild,"both_Pladias_known")
        self.early=[{t for t in s["all_original_taxa"] if t.startswith("w")} for s in self.sites]
        self.new=[{t for t in s["all_original_taxa"] if t.startswith("o")} for s in self.sites]

    def test_alien_subgroup_comparison_uses_identical_physical_road_rail_support(self):
        r1=reassign_original_source_class_outcome(self.baseline,self.sites,self.early)
        r2=reassign_original_source_class_outcome(self.baseline,self.sites,self.new)
        self.assertEqual(len(r1),len(r2))
        self.assertEqual([(r["i"],r["j"]) for r in r1],
                         [(r["i"],r["j"]) for r in r2])
        self.assertEqual([(r["source_section_i"],r["source_section_j"]) for r in r1],
                         [(r["source_section_i"],r["source_section_j"]) for r in r2])
        for a,b,original in zip(r1,r2,self.baseline):
            for f in ("road_route_ln","rail_to_road_route_ln","geo_ln",
                      "station_area_ln_mean","functional_overlap"):
                self.assertEqual(a["features"][f],original["features"][f])
                self.assertEqual(b["features"][f],original["features"][f])
            self.assertTrue(0<=a["jaccard"]<=1)
            self.assertTrue(0<=b["jaccard"]<=1)

    def test_same_original_five_source_rail_stretches_holdout_for_both_groups(self):
        for sets in (self.early,self.new):
            records=reassign_original_source_class_outcome(self.baseline,self.sites,sets)
            cv=cv_section_heldout(records)
            self.assertTrue(cv["identified"])
            self.assertEqual(cv["n_holdout_original_stretches"],5)
            self.assertEqual(cv["source_overlap_support"]["n_source_original_botanical_stations"],35)
            self.assertTrue(math.isfinite(cv["rail_added_MSE_improvement_percentage_after_road"]))

    def test_richness_sensitivity_frozen_original_vs_subgroup_preserves_source_outcome(self):
        a=reassign_original_source_class_outcome(
            self.baseline,self.sites,self.early,"subgroup_log1p_richness")
        b=reassign_original_source_class_outcome(
            self.baseline,self.sites,self.early,"whole_original_trait_known_log_richness")
        self.assertTrue(all(x["jaccard"]==y["jaccard"] for x,y in zip(a,b)))
        self.assertTrue(any(x["features"]["site_ln_species_count_mean"]!=
                            y["features"]["site_ln_species_count_mean"] for x,y in zip(a,b)))
        self.assertEqual(b[0]["features"]["site_ln_species_count_mean"],
                         self.baseline[0]["features"]["site_ln_species_count_mean"])

    def test_never_impute_double_empty_subgroup_or_missing_station_as_Jaccard_zero(self):
        absent=[set(x) for x in self.early]
        absent[0]=set()
        with self.assertRaisesRegex(ValueError,"absent"):
            reassign_original_source_class_outcome(self.baseline,self.sites,absent)
        with self.assertRaisesRegex(ValueError,"35"):
            reassign_original_source_class_outcome(self.baseline,self.sites,absent[:34])
        with self.assertRaisesRegex(ValueError,"35"):
            reassign_original_source_class_outcome(self.baseline,self.sites+[self.sites[0]],self.early)

    def test_original_v156_live_reference_is_frozen_and_requires_exact_agreement(self):
        cv=cv_section_heldout(self.baseline)
        d=cv["equal_original_stretch_weighted_mse"]
        assert_original_reference(cv,"synthetic",
          d["M1_plus_actual_2023_road_route"],
          d["M2_plus_actual_2023_road_and_rail_ratio"])
        with self.assertRaisesRegex(ValueError,"drifted"):
            assert_original_reference(cv,"synthetic",0,0)

if __name__=="__main__":unittest.main()

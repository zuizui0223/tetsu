"""v15.2 synthetic-only correctness of posthoc sensitivity grid; no fake field evidence."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_Slovenia_rail_road_heldout_identity_v151 import artificial_opportunity
from analyze_Slovenia_rail_road_heldout_species_identity_v151 import (
    model_records,section_heldout_validation,BASE)
from audit_Slovenia_OSM_rail_road_model_robustness_v152 import (
    new_outcome,heldout_variant,sensitivity_grid,SPECIFICATIONS,ALPHAS,
    EXPECTED_PBF_SHA,EXPECTED_MAIN_N_PAIRS
)

class FrozenSourceModelSensitivityTests(unittest.TestCase):
    def setUp(self):
        pairs,sites=artificial_opportunity()
        self.sites=sites
        for i,s in enumerate(sites):
            s["all_original_taxa"]=set(s["known"])|{f"unique_{i}"}
        self.rows,_=model_records(pairs,sites,True)

    def test_original_v151_model_source_prediction_reproduces_exactly(self):
        strict=section_heldout_validation(self.rows)
        same=heldout_variant(self.rows,BASE,10.)
        self.assertTrue(strict["identified"])
        self.assertEqual(same["n_heldout_source_sections"],10)
        self.assertAlmostEqual(same["source_section_equal_weight_road_MSE"],
            strict["source_section_equal_weight_MSE_by_model"]["plus_road"],places=11)
        self.assertAlmostEqual(same["source_section_equal_weight_road_plus_rail_MSE"],
            strict["source_section_equal_weight_MSE_by_model"]["plus_road_and_rail"],places=11)

    def test_observed_all_author_taxa_is_distinct_outcome(self):
        alt=new_outcome(self.rows,self.sites,True)
        self.assertEqual(len(alt),len(self.rows))
        self.assertNotEqual(alt[0]["jaccard"],self.rows[0]["jaccard"])
        self.assertEqual(new_outcome(self.rows,self.sites,False),self.rows)
        self.assertGreater(alt[0]["features"]["site_log_richness_mean"],
                           self.rows[0]["features"]["site_log_richness_mean"])

    def test_full_grid_reports_all_variants_even_negative(self):
        grid=sensitivity_grid(self.rows,self.sites)
        self.assertEqual(len(grid),2*len(SPECIFICATIONS)*len(ALPHAS))
        self.assertEqual({s["taxon_outcome"] for s in grid},{
            "261_trait_known_author_taxa","ALL_original_author_observed_BS_BT_taxa"})
        self.assertTrue(all(x["identified"] for x in grid))
        self.assertTrue(any(x["conditioning_specification"]=="sensitivity_omit_same_original_section"
                            for x in grid))
        self.assertEqual(len([x for x in grid if
           x["conditioning_specification"]=="strict_v151_source_section_and_functional"
           and x["ridge_alpha_prespecified_grid"]==10]),2)

    def test_invalid_original_all_taxa_cannot_be_synthetic_imputed(self):
        sites=[dict(x) for x in self.sites]
        sites[0].pop("all_original_taxa")
        with self.assertRaises(KeyError):
            new_outcome(self.rows,sites,True)

    def test_source_archive_SHA_and_pair_cohort_previously_frozen(self):
        self.assertEqual(len(EXPECTED_PBF_SHA),64)
        self.assertEqual(EXPECTED_MAIN_N_PAIRS,3896)

if __name__=="__main__":unittest.main()

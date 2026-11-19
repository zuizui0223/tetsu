"""Synthetic v29.6 QA only. Numerical fixtures NOT genuine biological source results."""
import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import validate_and_aggregate_four_rail_true_tree_clade_source_v296 as v

def source_fixture():
    previous={"version":"v29.5","study_backbone_source_coverage":[]}
    effects=[]
    for tree in v.BACKBONES:
        for i,study in enumerate(v.STUDIES):
            previous["study_backbone_source_coverage"].append({
                "phylogeny":tree,
                "original_study":study,
                "prespecified_exact_tip_no_grafting_inference_gate_status":
                    "DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_PASS",
                "exact_original_backbone_tip_taxon_n":str(v.SOURCE_COUNTS[tree][i]),
                "source_original_station_or_ballast_site_n":
                    str(v.EXPECTED_ORIGINAL_PHYSICAL_STUDY_SITES[i]),
                "original_field_site_tip_match_fraction_median":".83",
                "original_field_site_ge5_direct_tip_fraction":"1",
                "same_study_retained_sites_direct_tip_taxa_CWM_height_between_site_distance":"0.12",
            })
            for mode in v.MODES:
                effects.append({
                    "backbone":tree,"study":study,"null_model":mode,
                    "n_original_physical_field_sites":
                        str(v.EXPECTED_ORIGINAL_PHYSICAL_STUDY_SITES[i]),
                    "n_original_trait_known_exact_tree_tips":
                        str(v.SOURCE_COUNTS[tree][i]),
                    "n_orig_source_study_plant_geographical_groups":"3",
                    "original_observed_exact_backbone_CWM_height":"0.12",
                    "source_height_null_mean":"0.10",
                    "observed_minus_clade_null_height":"0.02",
                    "source_height_technical_upper_rank_NOT_field_p":"0.4",
                    "prespecified_support_status":("ORIGINAL_TIP_SUBSET_OCCUPANCY_NULL_REFERENCE"
                      if mode==v.MODES[0]
                      else "IDENTIFIABLE_DIRECT_BACKBONE_CLADE_LABEL_NULL"),
                    "n_technical_trait_label_shuffles_if_identifiable":str(v.N_PERM),
                    "original_site_species_geography_FIXED":"TRUE",
                    "2D_LEDA_height_and_seed_mass_paired_vector_preserved":"TRUE",
                })
    return effects,previous

class SourcePhylocladalAudit(unittest.TestCase):
    def test_4real_papers_3_alternative_nomenclatures_4_models(self):
        rows,previous=source_fixture()
        out=v.aggregate(rows,previous,"test tree DCF: no data")
        self.assertEqual(len(rows),48)
        self.assertEqual(out["true_sites_not_repeated_meta_studies"],171)
        self.assertEqual(out["original_railway_true_field_publications"],4)
        for tree in v.BACKBONES:
            for model in v.MODES:
                x=out["by_original_backbone_and_clade_null"][tree][model]
                self.assertEqual(x["original_publications_supported"],4)
                self.assertAlmostEqual(x["descriptive_height"]["original_source_pub_equal_height_excess"],.02)
                self.assertEqual(x["descriptive_height"]["original_source_publications_positive"],4)

    def test_no_pseudoreplication_from_tree_nomenclature(self):
        rows,previous=source_fixture()
        out=v.aggregate(rows,previous,"")
        self.assertEqual(out["original_railway_true_field_publications"],4)
        self.assertEqual(len(out["original_tree_versions_three_nomenclatures_not_three_independent_publications"]),3)

    def test_invalid_site_tip_denominator_fails_as_original_v295_QA(self):
        rows,previous=source_fixture()
        previous["study_backbone_source_coverage"][0]["original_field_site_tip_match_fraction_median"]="20.166666"
        with self.assertRaises(ValueError):
            v.aggregate(rows,previous,"")

    def test_unidentifiable_tree_clade_does_not_get_biological_zero(self):
        rows,previous=source_fixture()
        for x in rows:
            if x["backbone"]=="LCVP" and x["null_model"]==v.MODES[3]:
                x.update({
                    "prespecified_support_status":"NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY",
                    "n_technical_trait_label_shuffles_if_identifiable":"0",
                    "source_height_null_mean":"",
                    "observed_minus_clade_null_height":"",
                    "source_height_technical_upper_rank_NOT_field_p":"",
                })
        out=v.aggregate(rows,previous,"")
        self.assertEqual(out["by_original_backbone_and_clade_null"]["LCVP"][v.MODES[3]]
                         ["original_publications_supported"],0)
        self.assertIsNone(out["by_original_backbone_and_clade_null"]["LCVP"][v.MODES[3]]
                          ["descriptive_height"])

    def test_missing_true_paper_source_fails(self):
        rows,previous=source_fixture()
        with self.assertRaises(ValueError):
            v.aggregate(rows[:-1],previous,"")

    def test_study_tip_counts_frozen(self):
        rows,previous=source_fixture()
        rows[0]["n_original_trait_known_exact_tree_tips"]="999"
        with self.assertRaises(ValueError):
            v.aggregate(rows,previous,"")

    def test_trait_null_exact_subtraction_validated(self):
        rows,previous=source_fixture()
        rows[0]["observed_minus_clade_null_height"]="0.021"
        with self.assertRaises(ValueError):
            v.aggregate(rows,previous,"")

    def test_no_renamed_study_cohort(self):
        rows,previous=source_fixture()
        rows[0]["study"]="MadeUpRailCountry"
        with self.assertRaises(ValueError):
            v.aggregate(rows,previous,"")

    def test_nonidentified_clade_with_number_rejected(self):
        rows,previous=source_fixture()
        for x in rows:
            if x["backbone"]=="LCVP" and x["null_model"]==v.MODES[3]:
                x["prespecified_support_status"]="NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY"
                x["n_technical_trait_label_shuffles_if_identifiable"]="0"
        with self.assertRaises(ValueError):
            v.aggregate(rows,previous,"")

if __name__=="__main__":unittest.main()

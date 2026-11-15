"""v16.6 synthetic-only genuine-source audit invariants. These fake counts
are NOT the actual iNaturalist photo source observations.
"""
import copy,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_iNat_original20_SI_CZ_30m_source_eligibility_v166 import (
  run,validate_original_source,SOURCE_CANDIDATE_VERSION,PHOTO_SOURCE_VERSION,
  GEO_SHA
)

def fake_inputs():
    species=[f"Example species{i:02}" for i in range(20)]
    candidate={"version":SOURCE_CANDIDATE_VERSION,"complete":True,
        "n_distinct_original_species":20,
        "original_selected_20_species":[
           {"scientific_name_original_exact_binomial":s} for s in species]}
    photo={
      "complete":True,"version":PHOTO_SOURCE_VERSION,
      "NOT_exhaustive_global_rail_30m_count":True,
      "source_2018_2025_dates":["2018-01-01","2025-12-31"],
      "source_candidate_20_original_taxa_selection_manifest":SOURCE_CANDIDATE_VERSION,
      "n_country_species_strata_completed_without_API_cap":40,
      "n_country_species_strata_api_capped_and_NEAR30_UNESTIMATED_NOT_ZERO":0,
      "country_bbox_rail_geometry_original_source":{country:{"original_Geofabrik_PBF_SHA256":sha}
                                for country,sha in GEO_SHA.items()},
      "global_species_research_photo_totals_2018_2025_NOT_global_railway":[
            {"taxon":t,"global_iNaturalist_research_photo_total_2018_2025":10000} for t in species],
      "two_country_bbox_20taxon_real_iNat_photo_source_30m":{},
      "n_complete_strata_with_geometrically_point_within30m":1600,
      "n_complete_strata_with_photo_license_near30m":1200
    }
    for country in GEO_SHA:
        rows=[]
        for taxon in species:
            profile={
              "n_after_strict_accuracy_and_taxon_public_geoprivacy_gate":100,
              "n_point_near_30m_strict_accuracy":40,
              "n_location_certain_within_30m":20,
              "n_point_near_30m_photo_CC_licensed":30,
              "n_distinct_observers_public_within_30m":12,
              "n_distinct_100m_public_cells_near_30m":11,
              "n_point_100to1000m_strict_accuracy":40,
              "n_point_30to100m_strict_accuracy":20,
              "n_point_beyond_1000m_strict_accuracy":0,
              "near_30m_month_of_observation_counts_NOT_flowering_onset":{"05":40},
              "near_30m_explicit_flower_visibility_annotation_counts_NOT_full_phenophase":{"flower_visible_annotation":1}
            }
            rows.append({"taxon":taxon,
               "source_API_bbox_photo_total":120,
               "source_query_complete":True,
               "source_original_photo_quality_and_geo_rail_polyline_30m":profile})
        photo["two_country_bbox_20taxon_real_iNat_photo_source_30m"][country]=rows
    return candidate,photo

class RealPhotoSourceAuditSafetyTests(unittest.TestCase):
    def test_precise_source_photo_count_and_prospective_30m_gate(self):
        cand,source=fake_inputs()
        result=run(cand,source)
        self.assertEqual(result["n_original_taxa"],20)
        self.assertEqual(result["n_source_photo_rail_point_within30m_both_bboxes"],1600)
        self.assertEqual(result["n_location_within_30m_even_with_coordinate_uncertainty_circle"],800)
        self.assertEqual(result["n_country_taxon_cells_passing_predeclared_minimum_photographic_opportunity"],40)
        self.assertTrue(result["not_an_accurate_GLOBAL_near_railway_photo_count"])
        self.assertEqual(result["source_exact_original_20_species_rank_and_two_region_totals"][0][
           "source_full_global_iNat_research_photo_total_NOT_global_near30"],10000)

    def test_below_observer_gate_photo_count_does_not_imply_independence(self):
        c,s=fake_inputs()
        s["two_country_bbox_20taxon_real_iNat_photo_source_30m"]["CZ"][0][
           "source_original_photo_quality_and_geo_rail_polyline_30m"]["n_distinct_observers_public_within_30m"]=2
        a=run(c,s)
        self.assertEqual(a["n_country_taxon_cells_passing_predeclared_minimum_photographic_opportunity"],39)

    def test_near30_not_guaranteed_coordinate_accuracy_without_rechecking(self):
        c,s=fake_inputs()
        s["two_country_bbox_20taxon_real_iNat_photo_source_30m"]["SI"][0][
           "source_original_photo_quality_and_geo_rail_polyline_30m"]["n_location_certain_within_30m"]=41
        with self.assertRaisesRegex(ValueError,"guaranteed count"):
            run(c,s)

    def test_no_40_complete_source_query_becomes_false_zero(self):
        c,s=fake_inputs()
        s["n_country_species_strata_completed_without_API_cap"]=39
        with self.assertRaisesRegex(ValueError,"not completely paginated"):
            run(c,s)

    def test_country_geo_sha_or_study_row_label_tampering_fails_closed(self):
        c,s=fake_inputs()
        s["country_bbox_rail_geometry_original_source"]["CZ"]["original_Geofabrik_PBF_SHA256"]="false"
        with self.assertRaisesRegex(ValueError,"authentic historic OSM"):
            validate_original_source(c,s)
        c,s=fake_inputs()
        s["two_country_bbox_20taxon_real_iNat_photo_source_30m"]["SI"][0]["taxon"]="Taxon substitute"
        with self.assertRaisesRegex(ValueError,"candidate"):
            validate_original_source(c,s)

    def test_no_positional_accuracy_count_over_api_total(self):
        c,s=fake_inputs()
        s["two_country_bbox_20taxon_real_iNat_photo_source_30m"]["CZ"][0][
           "source_original_photo_quality_and_geo_rail_polyline_30m"]["n_after_strict_accuracy_and_taxon_public_geoprivacy_gate"]=200
        with self.assertRaisesRegex(ValueError,"source count"):
            run(c,s)

if __name__=="__main__":unittest.main()

"""v16.7 synthetic-only seasonal photo audit. No real original iNat photos.

The tests validate source metadata constraints and the phenology inference
firewall; original data are fetched from GitHub Actions study v16.6.
"""
import copy,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_iNat_original20_SI_CZ_30m_source_eligibility_v166 import fake_inputs
from audit_iNat_original20_SI_CZ_30m_source_eligibility_v166 import run as source_gate
from audit_iNat_rail30m_flower_annotation_review_gate_v167 import (
    evaluate,verified_rows,V166,MIN_ANNOTATED_FLOWER,
)

def synthetic_verified_source(*,inapplicable_horsetail=False):
    candidate,photo=fake_inputs()
    if inapplicable_horsetail:
        candidate["original_selected_20_species"][0]["scientific_name_original_exact_binomial"]="Equisetum arvense"
        photo["global_species_research_photo_totals_2018_2025_NOT_global_railway"][0]["taxon"]="Equisetum arvense"
        for region in ("SI","CZ"):
            photo["two_country_bbox_20taxon_real_iNat_photo_source_30m"][region][0]["taxon"]="Equisetum arvense"
    for region in ("SI","CZ"):
        for row in photo["two_country_bbox_20taxon_real_iNat_photo_source_30m"][region]:
            m=row["source_original_photo_quality_and_geo_rail_polyline_30m"]
            m["near_30m_month_of_observation_counts_NOT_flowering_onset"]={"04":12,"05":14,"06":14}
            m["near_30m_explicit_flower_visibility_annotation_counts_NOT_full_phenophase"]={
                "flower_visible_annotation":5,
                "no_verified_flower_visible_annotation":35
            }
    return source_gate(candidate,photo)

class OriginalSourcePhotoAnnotationFeasibilityTests(unittest.TestCase):
    def test_complete_40_strata_source_and_two_country_geo(self):
        source=synthetic_verified_source()
        self.assertEqual(len(verified_rows(source)),40)
        x=evaluate(source)
        self.assertEqual(x["source_total_gis_near30_photo_records_original"],1600)
        self.assertEqual(x["source_total_near30_tagged_FLOWERS_for_angiosperms_only"],200)
        self.assertEqual(x["n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate"],40)
        self.assertEqual(x["n_species_country_cells_with_min_5_positive_Flowers_tags"],40)
        self.assertTrue(x["none_of_new_source_counts_identifies_trait_selection"])

    def test_Equisetum_is_not_a_flowering_angiosperm(self):
        s=synthetic_verified_source(inapplicable_horsetail=True)
        x=evaluate(s)
        self.assertEqual(x["n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate"],38)
        self.assertEqual(x["source_total_near30_tagged_FLOWERS_for_angiosperms_only"],190)
        self.assertEqual(sum(v["is_spore_bearing_NON_flowering_horsetail"]
            for v in x["original_species_country_40cells_annotation_and_month_audit"]),2)
        self.assertTrue(all(not v["v167_ANGIOSPERM_positive_image_REVIEW_logistics_passed"]
            for v in x["original_species_country_40cells_annotation_and_month_audit"]
            if v["species"]=="Equisetum arvense"))

    def test_no_Flowers_annotation_does_not_mean_recorded_nonflowering(self):
        x=synthetic_verified_source()
        for row in x["source_species_country_counts"]:
            row["original_photo_explicit_flower_visibility_near30"]={
                "no_verified_flower_visible_annotation_NOT_a_negative":40,
            }
        with self.assertRaisesRegex(ValueError,"Unexpected source annotation"):
            evaluate(x)
        x=synthetic_verified_source()
        for row in x["source_species_country_counts"]:
            row["original_photo_explicit_flower_visibility_near30"]={
                "flower_visible_annotation":0,
                "no_verified_flower_visible_annotation":40
            }
        a=evaluate(x)
        self.assertEqual(a["n_species_country_cells_with_min_5_positive_Flowers_tags"],0)
        self.assertEqual(a["n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate"],40)

    def test_real_source_months_not_invented_first_flowering(self):
        x=synthetic_verified_source()
        x["source_species_country_counts"][0]["original_photo_observation_months_near30"]={"13":40}
        with self.assertRaisesRegex(ValueError,"calendar month invalid"):
            evaluate(x)
        x=synthetic_verified_source()
        x["source_species_country_counts"][0]["original_photo_observation_months_near30"]={"04":14,"05":26}
        result=evaluate(x)
        first=result["original_species_country_40cells_annotation_and_month_audit"][0]
        self.assertEqual(first["n_distinct_observed_months_with_at_least5_photos"],2)
        self.assertFalse(first["v167_ANGIOSPERM_positive_image_REVIEW_logistics_passed"])

    def test_near30_photo_total_must_match_annotated_unknown_plus_tagged(self):
        x=synthetic_verified_source()
        x["source_species_country_counts"][0]["original_photo_explicit_flower_visibility_near30"]={
            "flower_visible_annotation":3,
            "no_verified_flower_visible_annotation":35
        }
        with self.assertRaisesRegex(ValueError,"missing annotation QA"):
            evaluate(x)

    def test_reject_unverified_original_spatial_or_source(self):
        x=synthetic_verified_source()
        x["version"]="fabricated"
        with self.assertRaisesRegex(ValueError,"not authenticated"):
            evaluate(x)

if __name__=="__main__":unittest.main()

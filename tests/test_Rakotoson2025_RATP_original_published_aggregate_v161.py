"""v16.1 author-published RATP France railway 2025 source aggregate integrity.

This checks hand-transcribed PUBLIC publisher numbers and source boundary.
It is NOT a test that raw 394 survey-by-species rows were obtained.
"""
import copy,hashlib,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Rakotoson2025_RATP_original_published_aggregate_v161 import (
    assert_manifest,derive_author_published_aggregates,run
)
SOURCE=Path(__file__).resolve().parents[1]/"data/literature/Rakotoson2025_original_published_aggregate_not_raw_v161.json"

class OriginalRATPPublisherAggregateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=SOURCE.read_bytes()
        cls.source=json.loads(cls.raw)

    def test_original_394_source_field_surveys_and_170_management_reports_separate(self):
        s=self.source
        self.assertTrue(assert_manifest(s))
        site=s["field_author_samples_2021"]
        self.assertEqual(266+128,site["total_n_surveys"])
        self.assertEqual(site["total_n_surveys"],394)
        self.assertEqual(site["distinct_species_reported"],255)
        self.assertEqual(site["fixed_area_m2_per_vegetation_survey"],2)
        self.assertEqual(s["railway_management_alerts"]["source_weeding_alert_events"],170)
        self.assertFalse(s["new_plant_field_observations_collected_by_tetsu"])

    def test_published_source_mean_ratios_not_Hedges_g_or_biomass(self):
        x=derive_author_published_aggregates(self.source)
        self.assertAlmostEqual(x["mean_species_richness"]["derived_simple_unstandardized_embankment_minus_track_mean"],4.6)
        self.assertAlmostEqual(x["mean_species_richness"]["derived_descriptive_embankment_over_track_mean_RATIO"],7/2.4)
        self.assertAlmostEqual(x["mean_original_abundance_sum_scale"]["derived_descriptive_embankment_over_track_mean_RATIO"],52.6/5)
        self.assertAlmostEqual(x["mean_shannon"]["derived_simple_unstandardized_embankment_minus_track_mean"],.67)
        self.assertAlmostEqual(x["mean_pielou"]["derived_simple_unstandardized_embankment_minus_track_mean"],-.04)
        self.assertTrue(x["mean_original_abundance_sum_scale"]["NOT_meta_analysis_effect_size_or_railway_mowing_treatment_effect"])

    def test_no_plot_level_or_standardized_meta_claim_in_durable_output(self):
        result=run(self.source,self.raw)
        self.assertFalse(result["can_estimate_Hedges_g_or_valid_meta_variance_from_rounded_means"])
        self.assertFalse(result["can_estimate_mowing_causal_effect_with_170_alerts_not_treatment"])
        self.assertFalse(result["true_394_plot_x_255_taxon_author_matrix_currently_obtained"])
        self.assertTrue(result["source_manifest_is_human_transcribed_NOT_downloaded_original_raw_PLOT_DATA"])
        self.assertEqual(result["source_manifest_original_author_transcription_sha256"],hashlib.sha256(self.raw).hexdigest())

    def test_refuse_false_herbicide_before_after_or_mowed_count(self):
        x=copy.deepcopy(self.source)
        x["170_alerts_NOT_independent_mowing_treatments"]=False
        with self.assertRaisesRegex(ValueError,"limits"):
            assert_manifest(x)
        y=copy.deepcopy(self.source)
        y["railway_management_alerts"]["year_2019_herbicide_policy_change_NO_before_after_plot_cohort"]=False
        with self.assertRaisesRegex(ValueError,"glyphosate"):
            assert_manifest(y)

    def test_source_published_mean_source_count_tamper_fails_closed(self):
        x=copy.deepcopy(self.source)
        x["field_author_samples_2021"]["track_n_surveys"]=127
        with self.assertRaisesRegex(ValueError,"source count"):
            assert_manifest(x)
        x=copy.deepcopy(self.source)
        x["original_author_published_zone_average"]["mean_species_richness"]["track"]=2.5
        with self.assertRaisesRegex(ValueError,"mean_species_richness"):
            assert_manifest(x)
        x=copy.deepcopy(self.source)
        x["field_author_samples_2021"]["original_plot_by_taxon_presence_matrix_NOT_acquired"]=False
        with self.assertRaisesRegex(ValueError,"Fake claim"):
            assert_manifest(x)

if __name__=="__main__":unittest.main()

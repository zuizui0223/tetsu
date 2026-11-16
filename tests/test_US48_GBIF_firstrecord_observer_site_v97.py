"""v9.7 source-only, posthoc observer/coarse-site matched GBIF indexed
late-only GRIIS taxon category correlations; no seed immigration claims.
"""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_US48_GBIF_firstrecord_observer_site_v97 import (
    PREDICTORS,RESPONSES,CONTROLS,STATES,
    source_state_frame,rank_correlation,compute,run)

def fixture():
    observed={"complete":True,
        "nonrandom_source_probe_DO_NOT_EXTRAPOLATE_TO_ALL_INAT_USERS":True,
        "no_raw_observer_names_ids_record_ids_or_latlon_saved":True,
        "summary":{"n_actual_GBIF_iNaturalist_original_sampled_rows_before":9600,
                   "n_actual_GBIF_iNaturalist_original_sampled_rows_after":9600},
        "states":{}}
    record={"complete":True,"n_qualified_Census_US_states_for_both_year_windows":48,
            "state_FIPS_source_audit_rows":[]}
    for i in range(48):
        code=f"{i+1:02d}"
        def site(n,km):
            return {"n_late_records_in_early_sampled_geographical_cells":n,
                    "n_late_records_with_original_valid_coordinates":200,
                    "n_late_records_with_both_observer_and_coarse_location":200,
                    "n_late_records_with_same_observer_same_coarse_site_as_early":max(0,n//8),
                    "n_late_GRIIS_records_in_repeated_geographical_cells":n//4}
        n10=20+(i*13)%140
        n50=110+(i*17)%85
        observed["states"][code]={
            "source_pre_2018_21":{"n_original_rows":200},
            "source_post_2023_25":{
                "n_original_rows":200,
                "n_source_country_GRIIS_listed_plant_records_in_sample":70},
            "fraction_late_original_rows_with_repeated_prior_observer_label":
                ((i*7)%31)/100,
            "overlap":{
                "n_late_source_rows_with_observer_identity":200,
                "coarse_site_overlap":{"10km":site(n10,10),
                                       "50km":site(n50,50)}
            }
        }
        record["state_FIPS_source_audit_rows"].append({
            "Census_FIPS_state_code":code,
            "n_post_all_Plantae_indexed_records":20000+1000*i+20*(i%7),
            "post_before_all_Plantae_annual_record_ratio":1.2+(i%11)*.1,
            "post_only_fraction_of_post_recorded_taxon_categories":
                .1+((i*17)%37)/100,
            "post_only_first_recorded_taxon_categories_NOT_immigrations":
                30+(i*23)%87})
    return observed,record

class ExploratoryObserverSiteSourceAssociations(unittest.TestCase):
    def test_all_exploratory_metrics_reported_not_selected_only_if_significant(self):
        self.assertEqual(len(PREDICTORS),7)
        self.assertEqual(len(RESPONSES),2)
        self.assertEqual(STATES,48)
        self.assertEqual(len(CONTROLS),2)
        self.assertIn("sampled_observer_and_site_recurrence_50km",PREDICTORS)
    def test_all_48_same_original_state_data_sources(self):
        a,b=fixture()
        rows=source_state_frame(a,b)
        self.assertEqual(len(rows),48)
        self.assertEqual(set(rows["01"]["predictors"]),set(PREDICTORS))
        self.assertEqual(set(rows["01"]["responses"]),set(RESPONSES))
        self.assertEqual(set(rows["01"]["controls"]),set(CONTROLS))
    def test_no_imputed_observer_or_plant_recurrence_if_no_source_variance(self):
        a,b=fixture()
        r=list(source_state_frame(a,b).values())
        for row in r:row["predictors"]["sampled_observer_label_recurrence"]=0
        self.assertIsNone(rank_correlation(r,"sampled_observer_label_recurrence",RESPONSES[0]))
    def test_all_14_source_associations_include_full_state_jackknifes_no_pvalues(self):
        a,b=fixture()
        out=run(a,b)
        self.assertTrue(out["complete"])
        self.assertEqual(len(out["full_exploratory_correlation_matrix"]),7)
        for vals in out["full_exploratory_correlation_matrix"].values():
            self.assertEqual(len(vals),2)
            for result in vals.values():
                self.assertEqual(result["leave_one_Census_state_out_refits_not_independent_replicates"][
                    "n_deleted_states"],48)
                self.assertTrue(result["POSTHOC_hypothesis_no_p_value_multiple_predictors"])
        self.assertTrue(out["NO_rail_cargo_seed_causation_new_introduced_plants_mowing_wind_adaptation"])
    def test_missing_spatial_source_unit_not_assumed_botanical_absence(self):
        a,b=fixture()
        a["states"].pop("01")
        with self.assertRaisesRegex(ValueError,"source units disagree"):
            source_state_frame(a,b)
    def test_source_names_and_original_user_information_cannot_escape_to_model(self):
        a,b=fixture()
        output=run(a,b)
        content=json.dumps(output)
        self.assertNotIn("recordedBy",content)
        self.assertNotIn('"gbifID"',content)
        self.assertTrue(output["no_original_observer_labels_record_ids_or_geolocation_saved"])

if __name__=="__main__":unittest.main()

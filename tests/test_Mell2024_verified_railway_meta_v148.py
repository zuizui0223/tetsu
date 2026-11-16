"""v14.8 synthetic arithmetic/data-integrity tests, no invented source outputs."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Mell2024_verified_railway_meta_cases_v148 import rail_cases,penone_audit

def case(source,caseid,resp,lo="Railway",country="Australia",continent="Oceania"):
    return {"Study.ID":source,"Author_date":"Author_2019","Case.ID":caseid,
            "LTI":lo,"Country":country,"Continent":continent,"Bias":"Medium",
            "Response_variable":resp,"Treatment":"rail","Control":"outside",
            "N_treatment":"24","N_control":"24","Mean_treatment":"3",
            "Mean_control":"2","Var_treatment":"1","Var_control":"1",
            "Measure_var_treatment":"SD","Measure_var_control":"SD"}

def fixture():
    d=[case("none",f"road_{i}","Species_richness",lo="Road") for i in range(202)]
    d+= [case("EM_8","Agence_Gaiadomo_2012_1","Species_richness",country="France",continent="North_America"),
         case("11955","Everingham_2019_3","Species_richness"),
         case("11955","Everingham_2019_6","Abundance")]
    return d

class StrictReviewSourceGates(unittest.TestCase):
    def test_cases_not_independent_studies(self):
        d=rail_cases(fixture())
        self.assertEqual(d["all_original_meta_cases"],205)
        self.assertEqual(d["railway_cases"],3)
        self.assertEqual(d["unique_original_rail_studies"],2)
        self.assertAlmostEqual(d["rail_percent"],300/205)
        self.assertTrue(d["cases"][0]["original_geographic_metadata_conflict"])
        self.assertEqual([s["source_case_id"] for s in d["cases"]],
            ["Agence_Gaiadomo_2012_1","Everingham_2019_3","Everingham_2019_6"])
        self.assertEqual(d["cases"][1]["raw_mean_difference_rail_minus_comparator"],1.0)
    def test_negative_control_road_does_not_enter_rail(self):
        d=fixture()
        d[0]["LTI"]="Railway | Road"
        self.assertEqual(rail_cases(d)["railway_cases"],3)
    def test_wrong_site_case_source_rejected(self):
        d=fixture();d[-1]["Study.ID"]="made_up"
        with self.assertRaisesRegex(ValueError,"railway cases"):
            rail_cases(d)
    def test_205_case_pseudorep_boundary(self):
        with self.assertRaisesRegex(ValueError,"205 meta"):
            rail_cases(fixture()[:-1])
    def test_non_sd_is_rejected(self):
        d=fixture();d[-1]["Measure_var_control"]="SE"
        with self.assertRaisesRegex(ValueError,"standard deviations"):
            rail_cases(d)
    def test_penone_accepted_then_high_risk_review_is_one_cohort(self):
        base={"DOI":"10.1016/j.biocon.2012.01.041"}
        a=[base]
        c=[{**base,"Risk_of_bias":"High",
            "Other_bias":"Confounding factor: paved road and railway side by side",
            "Reasons_for_rejection":"High risk of bias (major confounding factor)"}]
        r=penone_audit(a,c)
        self.assertTrue(r["excluded_after_critical_appraisal"])
        self.assertFalse("experiment" in r["Mell_reviewers_confound"].lower())
        with self.assertRaisesRegex(ValueError,"changed"):
            penone_audit(a,[{**c[0],"Risk_of_bias":"Low"}])
if __name__=="__main__":unittest.main()

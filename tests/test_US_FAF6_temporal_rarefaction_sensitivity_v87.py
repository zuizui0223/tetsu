"""Only mathematical GBIF index source-record rarefaction + directional
FAF economic OD regression tests. No true seed movement or establishment.
"""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_macro_rail_assemblages_v22 import incidence_probabilities
from analyze_US_FAF6_temporal_rarefaction_sensitivity_v87 import (
    DEPTHS,RESPONSE,MIN_CANDIDATES,source_results,reweight_rows,
    state_probs_all_depths,COHORTS)

def mockdata():
    species={"before_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{
            str(i):10 for i in range(1001,1101)}},
        "after_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{
            str(i):10 for i in range(1020,1121)}}}
    temporal={"states":{"01":{"source_windows":species},
                         "02":{"source_windows":{
                             "before_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{
                                 str(i):10 for i in range(1070,1141)}},
                             "after_2022":{"US_country_GRIIS_positive_speciesKey_record_counts":{
                                 str(i):10 for i in range(1001,1141)}}}
                         }}}
    early1={str(i) for i in range(1001,1101)}
    early2={str(i) for i in range(1070,1141)}
    data={"01":{"early":early1},"02":{"early":early2}}
    probs=state_probs_all_depths(temporal,data)
    pair={"o":"01","d":"02","x_rail":math.log1p(100),
          "rail_2022_estimated_thousand_short_tons":100,
          "truck_2022_estimated_thousand_short_tons":20,
          "x_truck":math.log1p(20),
          "candidate_taxon_categories":70,
          "responses":{RESPONSE:.01,"raw_source_linked_first_record_fraction":.33},
          "covars":{}}
    return data,probs,pair

class KRecordSensitivityGate(unittest.TestCase):
    def test_all_depths_are_precommitted_without_selecting_ecological_direction(self):
        self.assertEqual(DEPTHS,(10,20,40,80,160))
        self.assertEqual(MIN_CANDIDATES,10)
    def test_same_candidate_cohort_all_depths_and_2022_source_not_absence(self):
        eligible,prob,row=mockdata()
        outputs={}
        for k in DEPTHS:
            result=reweight_rows([row],eligible,prob,k)[0]
            outputs[k]=result["responses"][RESPONSE]
            self.assertEqual(result["o"],"01")
            self.assertEqual(result["d"],"02")
            self.assertEqual(result["rail_2022_estimated_thousand_short_tons"],100)
            self.assertTrue(0<=outputs[k]<=1)
        self.assertEqual(len(outputs),5)
        self.assertEqual(len(eligible["01"]["early"]-eligible["02"]["early"]),69)
    def test_reverse_rail_OD_preserves_same_candidate_plant_sources(self):
        eligible,prob,row=mockdata()
        original=reweight_rows([row],eligible,prob,20)[0]
        reversed_result=reweight_rows([row],eligible,prob,20,
                                      rail_OD={("02","01"):3})[0]
        self.assertEqual(original["responses"][RESPONSE],reversed_result["responses"][RESPONSE])
        self.assertEqual(original["rail_2022_estimated_thousand_short_tons"],
                         reversed_result["rail_2022_estimated_thousand_short_tons"])
        self.assertAlmostEqual(reversed_result["x_rail"],math.log1p(3))
        self.assertAlmostEqual(original["x_rail"],math.log1p(100))
        self.assertEqual(COHORTS[1],"positive_RAIL_crop_OD")
    def test_undefined_state_pool_fails_closed_no_imputed_arrival(self):
        eligible,prob,row=mockdata()
        b=dict(row,d="01")
        with self.assertRaisesRegex(ValueError,"Candidate group"):
            reweight_rows([b],eligible,prob,20)
    def test_source_support_is_not_fake_2_state_independent_inference(self):
        eligible,prob,row=mockdata()
        x=source_results(reweight_rows([row],eligible,prob,20))
        self.assertFalse(x[COHORTS[0]]["supported"])
        self.assertIsNone(x[COHORTS[0]]["rail_rho_full_controls"])
    def test_invalid_unfrozen_K_rejected(self):
        eligible,prob,row=mockdata()
        with self.assertRaisesRegex(ValueError,"Unpredeclared"):
            reweight_rows([row],eligible,prob,25)
if __name__=="__main__":unittest.main()

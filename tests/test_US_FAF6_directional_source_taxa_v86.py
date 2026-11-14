"""v8.6 synthetic, dependency-stable source directional 2022 FAF6
freight vs 2018–21 / 2023–25 recorded taxon candidate regression.
No true colonization, transport or adaptation hypotheses validated.
"""
from __future__ import annotations
import sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_US_FAF6_directional_source_taxa_v86 import (
    K,COHORTS,PRIMARY_CONTROLS,RESPONSES,source_states,
    recorded_source_candidate,construct_directions,
    cohort,partial_rank,node_deletion)
from analyze_macro_rail_assemblages_v22 import incidence_probabilities

def state(code,before,after):
    before={str(k):15 for k in before}
    after={str(k):15 for k in after}
    return {"early":set(before),"late":set(after),
            "early_prob20":incidence_probabilities(before,k=20),
            "late_prob20":incidence_probabilities(after,k=20),
            "effort_early":1000,"effort_late":1200}

def fakefa():
    return {
        "complete":True,"source_hash_verified_vs_prior_v79":True,
        "source_trace":{"previous_verified_v79_total_and_modes_reconciled":True},
        "source_year":2022,
        "raw_official_state_OD_2022_crop_estimated_tonnage":{
            "RAIL":[{"o":"01","d":"02","estimated_thousand_short_tons":100.},
                    {"o":"02","d":"01","estimated_thousand_short_tons":1.}],
            "TRUCK":[{"o":"01","d":"02","estimated_thousand_short_tons":20.},
                     {"o":"02","d":"01","estimated_thousand_short_tons":40.}]},
    }

class RecordedNotColonized(unittest.TestCase):
    def test_2022_no_overlap_k20_and_original_train_controls(self):
        self.assertEqual(K,20)
        self.assertEqual(len(COHORTS),3)
        self.assertIn("log1p_truck_OD",PRIMARY_CONTROLS)
        self.assertEqual(RESPONSES[0],"expected_20record_source_linked_first_record_fraction")
    def test_at_risk_source_candidates_are_unrecorded_not_truly_absent(self):
        origin=state("01",range(1001,1031),range(1001,1031))
        dest=state("02",range(1020,1041),range(1010,1041))
        ans=recorded_source_candidate(origin,dest)
        self.assertEqual(ans["n_origin_before_recorded_taxa_not_recorded_destination_before"],19)
        self.assertEqual(ans["n_origin_taxa_first_recorded_destination_later"],10)
        self.assertAlmostEqual(ans["raw_source_linked_first_record_fraction"],10/19)
        self.assertTrue(ans["NOT_true_source_to_dest_propagule_arrival"])
        self.assertTrue(0<=ans[RESPONSES[0]]<=1)
    def test_no_candidate_taxa_no_invented_new_species(self):
        a=state("01",range(1001,1031),range(1001,1031))
        b=state("02",range(1001,1031),range(1001,1031))
        self.assertIsNone(recorded_source_candidate(a,b))
    def test_real_OD_direction_must_not_be_symmetrized_as_true_plant_arrival(self):
        s={"01":state("01",range(1001,1031),range(1001,1031)),
           "02":state("02",range(1020,1051),range(1010,1051))}
        g={"01":{"lat":40,"lon":-90,"land_sq_m":10000,"region":"2","division":"5"},
           "02":{"lat":42,"lon":-89,"land_sq_m":12000,"region":"2","division":"5"}}
        with patch("analyze_US_FAF6_directional_source_taxa_v86.source_states",
                   return_value=(s,{})):
            rows,matched,excluded=construct_directions({},fakefa(),g)
        self.assertEqual(len(rows),2)
        fwd=next(v for v in rows if v["o"]=="01" and v["d"]=="02")
        rev=next(v for v in rows if v["o"]=="02" and v["d"]=="01")
        self.assertEqual(fwd["rail_2022_estimated_thousand_short_tons"],100.)
        self.assertEqual(rev["rail_2022_estimated_thousand_short_tons"],1.)
        self.assertNotEqual(fwd["responses"][RESPONSES[1]],rev["responses"][RESPONSES[1]])
        self.assertEqual(len(cohort(rows,COHORTS[2])),2)
        self.assertIsNone(partial_rank(rows,RESPONSES[0]))
        self.assertTrue(node_deletion(rows,RESPONSES[0])["NOT_independent_country_replication_or_inferential_p_value"])
    def test_quality_not_eligible_state_is_not_species_absence(self):
        states={}
        for code in range(1,52):
            e={"complete":True,"source_eligible":code<=45,
               "n_all_Plantae_indexed_records":1000,
               "US_country_GRIIS_positive_speciesKey_record_counts":{
                   str(2000+i+code):15 for i in range(20)}}
            l=dict(e)
            states[f"{code:02d}"]={"source_windows":{"before_2022":e,"after_2022":l},
                "two_year_source":{"two_period_quality":code<=45}}
        x={"complete":True,"n_qualified_source_states_gate_40":True,
           "source_periods":{"before_2022":"2018,2021","after_2022":"2023,2025"},
           "source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows":True,
           "states":states}
        good,missing=source_states(x)
        self.assertEqual(len(good),45)
        self.assertEqual(len(missing),6)
    def test_no_2022_overlapping_year_window_can_enter_main_result(self):
        data={"complete":True,"n_qualified_source_states_gate_40":True,
              "source_periods":{"before_2022":"2018,2022","after_2022":"2023,2025"},
              "source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows":False,
              "states":{f"{i:02d}":{} for i in range(1,52)}}
        with self.assertRaisesRegex(ValueError,"Temporal"):
            source_states(data)
if __name__=="__main__":unittest.main()

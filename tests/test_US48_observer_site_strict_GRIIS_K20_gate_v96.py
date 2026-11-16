"""v9.6 tests of mathematically necessary, not sufficient, K20
same-observer and same 10/50km site GRIIS RECORD SOURCE support.

All synthetic only; no individual observer identifiers, original
coordinates, biological colonization or freight seed movement.
"""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US48_observer_site_strict_GRIIS_K20_gate_v96 import (
    K,MIN_MACRO_STATES,SOURCE_STATES,SOURCE_KM,support_by_state,run)


def fixture_state(obs10=12,gr10=18,obs50=30,gr50=40):
    def tile(obs,gr,site):
        return {
            "n_late_records_with_same_observer_same_coarse_site_as_early":obs,
            "n_late_GRIIS_records_in_repeated_geographical_cells":gr,
            "n_late_records_in_early_sampled_geographical_cells":site,
            "n_late_records_with_original_valid_coordinates":200,
        }
    return {
        "source_pre_2018_21":{"n_original_rows":200},
        "source_post_2023_25":{"n_original_rows":200},
        "overlap":{"coarse_site_overlap":{
            "10km":tile(obs10,gr10,70),
            "50km":tile(obs50,gr50,160)
        }},
    }


def v95_fixture():
    states={}
    for i in range(48):
        states[f"{i+1:02d}"]=fixture_state(
            obs10=25 if i<3 else 5,
            gr10=25 if i<22 else 9,
            obs50=40 if i<14 else 10,
            gr50=45 if i<36 else 12,
        )
    return {
        "complete":True,
        "bounded_deterministic_search_offsets":[0,1000],
        "no_raw_observer_names_ids_record_ids_or_latlon_saved":True,
        "nonrandom_source_probe_DO_NOT_EXTRAPOLATE_TO_ALL_INAT_USERS":True,
        "states":states
    }


class ObserverSiteMinimumSupportTests(unittest.TestCase):
    def test_frozen_thresholds_and_not_invented_independence(self):
        self.assertEqual(K,20)
        self.assertEqual(MIN_MACRO_STATES,40)
        self.assertEqual(SOURCE_STATES,48)
        self.assertEqual(SOURCE_KM,(10,50))
    def test_joint_GRIIS_support_is_only_upper_bound(self):
        a=support_by_state(fixture_state())
        self.assertEqual(a["10km"]["upper_bound_GRIIS_late_records_same_observer_AND_same_site"],12)
        self.assertFalse(a["10km"]["necessary_not_sufficient_K20_GRIIS_source_record_support"])
        self.assertEqual(a["50km"]["upper_bound_GRIIS_late_records_same_observer_AND_same_site"],30)
        self.assertTrue(a["50km"]["necessary_not_sufficient_K20_GRIIS_source_record_support"])
        self.assertIsNone(a["50km"]["n_valid_species_keys_repeat_observer_and_same_site_UNAVAILABLE"])
    def test_site_repeat_can_be_common_while_GRIIS_observer_site_not_supported(self):
        d=v95_fixture()
        v=run(d)
        self.assertEqual(v["summary_by_coarse_site_km"]["10km"][
            "n_states_for_which_K20_GRIIS_same_observer_same_site_is_EVEN_POSSIBLE_UPPER_BOUND"],3)
        self.assertEqual(v["summary_by_coarse_site_km"]["50km"][
            "n_states_for_which_K20_GRIIS_same_observer_same_site_is_EVEN_POSSIBLE_UPPER_BOUND"],14)
        self.assertFalse(v["macro_observer_site_matched_2022_FAF6_GRIIS_plant_effect_estimable"])
        self.assertTrue(v["do_not_create_freight_plant_effect_when_source_gate_fails"])
        self.assertNotIn("recordedBy",json.dumps(v))
    def test_invalid_aggregate_count_not_force_observer_species_overlap(self):
        d=fixture_state(obs10=199,gr10=18)
        with self.assertRaisesRegex(ValueError,"exceed original repeat site"):
            support_by_state(d)
        d=fixture_state()
        d["overlap"]["coarse_site_overlap"]["50km"][
            "n_late_records_with_same_observer_same_coarse_site_as_early"]=201
        with self.assertRaisesRegex(ValueError,"bounded"):
            support_by_state(d)
    def test_source_original_state_count_and_deterministic_probe_must_not_drift(self):
        d=v95_fixture()
        d["states"].pop("01")
        with self.assertRaisesRegex(ValueError,"provenance"):
            run(d)
        d=v95_fixture()
        d["bounded_deterministic_search_offsets"]=[0,50000]
        with self.assertRaisesRegex(ValueError,"provenance"):
            run(d)
    def test_stop_and_reinspect_if_more_than_40_states_could_support_K20(self):
        d=v95_fixture()
        for val in d["states"].values():
            val["overlap"]["coarse_site_overlap"]["10km"][
                "n_late_records_with_same_observer_same_coarse_site_as_early"]=30
            val["overlap"]["coarse_site_overlap"]["10km"][
                "n_late_GRIIS_records_in_repeated_geographical_cells"]=30
        with self.assertRaisesRegex(ValueError,"joint source taxon"):
            run(d)

if __name__=="__main__":unittest.main()

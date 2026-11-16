"""Synthetic four-source GBIF/GRIIS/rail/road noncausal contrast checks."""
import math
from itertools import combinations
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_griis_rail_with_road_density_v27 import (
    validate_transport_city,build_city_pairs,compare_matched_pairs,
    city_summary,aggregate_blocks,analyze,
    COUNTRIES,PERIODS
)


def sources(name,country):
    bio={"country":country,"cells":{}}
    rc={"country":country,"cell_component_cover":{},"cell_pair_connection_descriptors":[]}
    road={"country":country,"cells":{}}
    for i in range(9):
        c=f"{name}_x{i%3}_y{i//3}"
        counts={str(100+j):20 for j in range(5)}
        counts["999"]=200
        s={"full_facet_pagination_completed":True,
           "source_eligible_for_recording_effort_screen":True,
           "n_source_datasets_at_least":4,
           "indexed_occurrence_count":300,
           "n_species_key_categories":6,
           "species_key_occurrence_counts":counts}
        bio["cells"][c]={"bbox":[120+(i%3)*.1,30+(i//3)*.1,120+(i%3+1)*.1,30+(i//3+1)*.1],
                         "by_period":{p:s for p in PERIODS}}
        rc["cell_component_cover"][c]={"n_source_components":1 if i<5 else 0}
        density=.2+i*.025
        road["cells"][c]={"area_km2":50.0,"road_km":density*50,
                         "rail_km":12. if i<5 else 0.,
                         "road_length_density_km_per_km2":round(density,8)}
    ids=list(bio["cells"])
    for a,b in combinations(ids,2):
        rc["cell_pair_connection_descriptors"].append({
            "a":a,"b":b,
            "within_clipped_graph_connected":bool(
                rc["cell_component_cover"][a]["n_source_components"] and
                rc["cell_component_cover"][b]["n_source_components"])
        })
    return bio,rc,road


class RailRoadAlienSourceTests(unittest.TestCase):
    def test_real_road_cell_source_density_validation(self):
        a,b,c=sources("jp_sendai","JP")
        x=validate_transport_city(c,a,b)
        self.assertEqual(len(x),9)
        self.assertGreater(x["jp_sendai_x0_y0"]["road_density"],0)
        c["cells"]["jp_sendai_x0_y0"]["road_length_density_km_per_km2"]=8.
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            validate_transport_city(c,a,b)

    def test_alien_species_rail_graph_pairs_are_source_qualified(self):
        a,b,c=sources("jp_sendai","JP")
        r,qc=build_city_pairs(a,b,c,set(str(100+i) for i in range(100)),"2018_2025")
        self.assertEqual(len(r),36)
        self.assertEqual(qc["n_record_eligible_intro_plant_cells"],9)
        self.assertEqual(qc["n_cell_rail_presence_disagreements_connector_vs_clipped_length"],0)
        self.assertEqual(sum(v["connected"] for v in r),10)
        self.assertTrue(all(0<=v["source_rarefied20_introduced_similarity"]<=1 for v in r))

    @staticmethod
    def candidate(is_rail,idx,road_mean=0.5):
        return {
            "a":str(idx),"b":str(idx+20),"both_connector_source_rail":True,
            "connected":is_rail,"distance_km":10+idx,
            "log_record_effort_ratio":.1,
            "log_mean_record_count":5.,
            "road_mean_log_density":road_mean,
            "road_absolute_log_ratio":.2,
            "source_rarefied20_introduced_similarity":.6 if is_rail else .4
        }

    def test_road_caliper_rejects_confounded_matches(self):
        rows=[self.candidate(True,i) for i in range(3)]+[
            self.candidate(False,i,2.) for i in range(3)]
        base=compare_matched_pairs(rows,match_road=False)
        road=compare_matched_pairs(rows,match_road=True)
        self.assertEqual(base["n_covariate_matched_source_pair_contrasts"],3)
        self.assertAlmostEqual(base["median_city_pair_similarity_contrast"],.2)
        self.assertEqual(road["n_covariate_matched_source_pair_contrasts"],0)
        self.assertIsNone(road["median_city_pair_similarity_contrast"])

    def test_matched_road_noncausal_identical_source(self):
        rows=[self.candidate(True,i) for i in range(4)]+[
            self.candidate(False,i,.5) for i in range(4)]
        x=compare_matched_pairs(rows,match_road=True)
        self.assertEqual(x["n_covariate_matched_source_pair_contrasts"],4)
        self.assertAlmostEqual(x["median_city_pair_similarity_contrast"],.2)
        self.assertNotIn("p_value",x)

    def test_city_no_road_is_not_treated_as_genuine_road_free(self):
        a,b,c=sources("jp_sendai","JP")
        c["cells"]["jp_sendai_x0_y0"]["road_km"]=0
        c["cells"]["jp_sendai_x0_y0"]["road_length_density_km_per_km2"]=0.
        rows,qc=build_city_pairs(a,b,c,set(str(100+i) for i in range(100)),"2018_2025")
        self.assertEqual(qc["n_record_eligible_intro_plant_cells"],9)
        self.assertEqual(len(rows),28)

    def test_country_block_support_not_faked_from_many_city_pairs(self):
        z={"jp_sendai":{"periods":{PERIODS[1]:{"comparisons":{
            "road_intensity_all":{"source_support_gate_pass":True,
                                  "median_city_pair_similarity_contrast":.03}}}}}}
        r=aggregate_blocks(z,PERIODS[1],"road_intensity_all")
        self.assertEqual(r["n_informative_city_blocks"],1)
        self.assertFalse(r["minimum_geographic_replication_met"])

    def test_missing_prerequisite_rejected(self):
        with self.assertRaisesRegex(ValueError,"Incomplete prerequisite"):
            analyze({"source_complete":False},{"complete":True},{},{"complete":True})
        with self.assertRaisesRegex(ValueError,"Ten-country"):
            analyze({"source_complete":True,"regions":{}},
                    {"complete":True,"regions":{}},{},
                    {"complete":True,"regions":{}})


if __name__=="__main__":
    unittest.main()

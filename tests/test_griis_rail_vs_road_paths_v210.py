"""Offline five-source physical ROAD PATH × railway alien plant QA (noncausal)."""
import sys
from pathlib import Path
from itertools import combinations
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_griis_rail_vs_road_paths_v210 import (
    SENSITIVITY_CALIPERS,validate_road_graph,add_road_paths,
    compare,summarize_city_blocks,analyze
)


def source_pairs():
    cells=[f"jp_sendai_x{i%3}_y{i//3}" for i in range(9)]
    road=[]
    rail=[]
    for a,b in combinations(cells,2):
        road.append({
            "a":a,"b":b,
            "within_clipped_road_graph_physically_connected":True,
            "source_centroids_snapped_within_5km":True,
            "optimistic_road_connector_path_km":10.,
            "optimistic_road_to_euclidean_ratio":1.3,
        })
        rail.append({
            "a":a,"b":b,
            "within_clipped_graph_connected":True
        })
    return {"country":"JP","cell_pair_road_paths":road},{"country":"JP","cell_pair_connection_descriptors":rail}


def example_pair(exposed,i,path=12.,circuity=1.2):
    return {
        "a":str(i),"b":str(i+10),"distance_km":10.+i,
        "connected":exposed,"both_connector_source_rail":True,
        "log_mean_record_count":5.,"log_record_effort_ratio":.2,
        "road_mean_log_density":-.3,"road_absolute_log_ratio":.1,
        "road_network_shortest_path_km":path,
        "road_network_circuity":circuity,
        "source_rarefied20_introduced_similarity":.6 if exposed else .4
    }


class RoadPathAlienValidation(unittest.TestCase):
    def test_real_source_road_graph_matches_rail_grid(self):
        rd,rr=source_pairs()
        result=validate_road_graph("jp_sendai","JP",rd,rr)
        self.assertEqual(len(result),36)
        self.assertEqual(sum(v["road_network_source_usable"] for v in result.values()),36)
        self.assertEqual(result[("jp_sendai_x0_y0","jp_sendai_x1_y0")]["network_length_km"],10)

    def test_unconnected_road_or_unsnapped_cell_does_not_impute_path(self):
        rd,rr=source_pairs()
        rd["cell_pair_road_paths"][0]["source_centroids_snapped_within_5km"]=False
        res=validate_road_graph("jp_sendai","JP",rd,rr)
        self.assertEqual(sum(v["road_network_source_usable"] for v in res.values()),35)
        rd["cell_pair_road_paths"][0]["within_clipped_road_graph_physically_connected"]=False
        self.assertFalse(validate_road_graph("jp_sendai","JP",rd,rr)
                         [("jp_sendai_x0_y0","jp_sendai_x1_y0")]["road_network_source_usable"])

    def test_fail_closed_when_pair_frame_missing(self):
        rd,rr=source_pairs()
        rd["cell_pair_road_paths"].pop()
        with self.assertRaisesRegex(ValueError,"Incomplete source road"):
            validate_road_graph("jp_sendai","JP",rd,rr)
        rd,rr=source_pairs()
        rd["country"]="GB"
        with self.assertRaises(ValueError):
            validate_road_graph("jp_sendai","JP",rd,rr)

    def test_source_matched_positive_and_negative_example(self):
        rows=[example_pair(True,i) for i in range(4)]
        rows.extend(example_pair(False,i) for i in range(4))
        x=compare(rows,calipers=SENSITIVITY_CALIPERS["main"])
        self.assertTrue(x["city_comparison_supported"])
        self.assertEqual(x["n_matched_exposed_control_pair_contrasts"],4)
        self.assertAlmostEqual(x["city_level_median_similarity_delta"],.2)
        self.assertNotIn("p_value",x)

    def test_road_distance_circuity_controls_can_remove_support(self):
        rows=[example_pair(True,i) for i in range(4)]
        rows.extend(example_pair(False,i,path=120.,circuity=3.) for i in range(4))
        for label in SENSITIVITY_CALIPERS:
            r=compare(rows,calipers=SENSITIVITY_CALIPERS[label])
            self.assertFalse(r["city_comparison_supported"])
            self.assertIsNone(r["city_level_median_similarity_delta"])

    def test_city_block_no_pseudoreplicate_p(self):
        entry={"periods":{"2018_2025":{"comparisons":{"all":{"main":{
            "city_comparison_supported":True,
            "city_level_median_similarity_delta":.02
        }}}}}}
        z=summarize_city_blocks({"sendai":entry},"2018_2025","all","main")
        self.assertEqual(z["n_informative_independent_city_blocks"],1)
        self.assertFalse(z["minimum_five_city_blocks"])

    def test_incomplete_original_source_gates_railway_effect(self):
        with self.assertRaisesRegex(ValueError,"Incomplete fifth"):
            analyze({"source_complete":False},{"complete":True},
                    {"complete":True},{"complete":True},{"complete":True})
        with self.assertRaisesRegex(ValueError,"Ten-city"):
            analyze({"source_complete":True,"regions":{}},
                    {"complete":True,"regions":{}},
                    {"complete":True,"countries":{}},
                    {"complete":True,"regions":{}},
                    {"complete":True,"regions":{}})


if __name__=="__main__":
    unittest.main()

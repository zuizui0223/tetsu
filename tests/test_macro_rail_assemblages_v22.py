"""Source-denominator and noncausal rail connection assemblage validation tests."""
import sys
from pathlib import Path
from itertools import combinations
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_macro_rail_assemblages_v22 import (
    source_species,incidence_probabilities,expected_rarefied_overlap,jaccard,
    geographic_distance_km,calculate_pairs,matched_contrasts,analyze,
    COUNTRIES,PERIODS
)


def source(keys,records=300,datasets=3,eligible=True):
    return {
        "full_facet_pagination_completed":True,
        "source_eligible_for_recording_effort_screen":eligible,
        "n_source_datasets_at_least":datasets,
        "indexed_occurrence_count":records,
        "species_key_occurrence_counts":keys,
        "n_species_key_categories":len(keys),
    }


def fake_city(name):
    cells={}
    graph_cells={}
    for i in range(9):
        cid=f"{name}_x{i%3}_y{i//3}"
        cells[cid]={
            "bbox":[120+(i%3)*.1,30+(i//3)*.1,120+(i%3+1)*.1,30+(i//3+1)*.1],
            "by_period":{
                "2010_2017":source({"s":200,"x"+str(i%3):100}),
                "2018_2025":source({"s":150,"y"+str(i):150}),
            },
        }
        graph_cells[cid]={"n_source_components":1 if i<5 else 0}
    pairs=[]
    ids=list(cells)
    for a,b in combinations(ids,2):
        both=graph_cells[a]["n_source_components"] and graph_cells[b]["n_source_components"]
        pairs.append({"a":a,"b":b,"within_clipped_graph_connected":bool(both)})
    return {"cells":cells,"country":"JP"},{"cell_component_cover":graph_cells,
            "cell_pair_connection_descriptors":pairs,
            "n_distinct_components_in_source_window":1,
            "country":"JP"}


class RailGBIFAssemblageTests(unittest.TestCase):
    def test_source_gates_and_missingness(self):
        self.assertIsNone(source_species(source({"1":99},99)))
        self.assertIsNone(source_species(source({"1":200},300,datasets=2)))
        self.assertEqual(source_species(source({"1":200,"2":100})),{"1":200,"2":100})
        with self.assertRaises(ValueError):
            source_species(source({"1":250,"2":100},records=300))

    def test_rarefaction_if_sample_whole_source_matches_incidence(self):
        pa=incidence_probabilities({"1":50,"2":50})
        pb=incidence_probabilities({"2":70,"3":30})
        self.assertEqual(pa["1"],1.)
        self.assertEqual(pb["2"],1.)
        self.assertAlmostEqual(expected_rarefied_overlap(pa,pb),1/3)
        self.assertAlmostEqual(jaccard({"1":50,"2":50},{"2":70,"3":30}),1/3)

    def test_rarefaction_not_full_composition_when_record_count_large(self):
        p=incidence_probabilities({"common":100,"rare":1,"other":199})
        self.assertGreater(p["common"],.99)
        self.assertLess(p["rare"],.5)
        self.assertGreater(p["rare"],0.)
        with self.assertRaises(ValueError):
            incidence_probabilities({"a":25},k=100)

    def test_geographic_distance_and_null(self):
        box=[120.0,30.0,120.1,30.1]
        self.assertEqual(geographic_distance_km(box,box),0)
        self.assertGreater(geographic_distance_km(box,[120.1,30.,120.2,30.1]),9)

    def test_physical_rail_pair_counts_do_not_impute_unmapped_rails(self):
        b,g=fake_city("jp_sendai")
        pairs,cells=calculate_pairs(b,g,"2018_2025")
        self.assertEqual(cells,9)
        self.assertEqual(len(pairs),36)
        self.assertEqual(sum(x["mapped_rail_connected_within_window"] for x in pairs),10)
        self.assertEqual(sum(x["both_cells_have_source_rail"] for x in pairs),10)
        self.assertTrue(all(x["ratio_expected_rarefied_100_record_intersection_to_union"]>=0 for x in pairs))

    def test_matched_pair_deltas_are_descriptive_and_no_fake_p(self):
        rows=[]
        for i in range(4):
            for rail in (True,False):
                rows.append({"a":str(i),"b":str(i)+str(rail),
                             "both_cells_have_source_rail":True,
                             "mapped_rail_connected_within_window":rail,
                             "geographic_distance_km":10+i,
                             "log_occurrence_effort_ratio":.2,
                             "log_geomean_records":4.5,
                             "ratio_expected_rarefied_100_record_intersection_to_union":.55 if rail else .3})
        result=matched_contrasts(rows)
        self.assertEqual(result["n_distance_and_record_effort_matched_pairs"],4)
        self.assertAlmostEqual(result["median_matched_similarity_delta"],.25)
        self.assertNotIn("p_value",result)

    def test_all_city_analysis_fails_on_missing_source_not_assume_zero(self):
        with self.assertRaisesRegex(ValueError,"Incomplete source"):
            analyze({"source_complete":False},{"complete":True})
        with self.assertRaisesRegex(ValueError,"Country frame"):
            analyze({"source_complete":True,"regions":{}},{"complete":True,"regions":{}})

    def test_unequal_mapping_breaks_not_invents_species_effect(self):
        b,g=fake_city("jp_sendai")
        g["cell_pair_connection_descriptors"][0]["a"]="wrong_cell"
        with self.assertRaises(ValueError):
            calculate_pairs(b,g,"2018_2025")


if __name__=="__main__":
    unittest.main()

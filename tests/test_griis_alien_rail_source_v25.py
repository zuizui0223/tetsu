"""Synthetic triple-source GRIIS + GBIF + Overture noncausal comparison tests."""
import sys
from pathlib import Path
from itertools import combinations
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_griis_alien_rail_source_v25 import (
    valid_griis,alien_source_species,cell_observations,city_pair_analysis,analyze,
    COUNTRIES,PERIODS
)


def plant_source(valid=True):
    counts={str(100+i):20 for i in range(5)}
    counts["999"]=300
    return {
        "full_facet_pagination_completed":True,
        "source_eligible_for_recording_effort_screen":valid,
        "indexed_occurrence_count":400,
        "n_source_datasets_at_least":4,
        "species_key_occurrence_counts":counts,
        "n_species_key_categories":len(counts)
    }


def fake_city(name,country):
    cells={};graph_cells={}
    for i in range(9):
        cid=f"{name}_x{i%3}_y{i//3}"
        cells[cid]={"bbox":[120+i%3*.1,30+(i//3)*.1,120+(i%3+1)*.1,30+(i//3+1)*.1],
                    "by_period":{p:plant_source() for p in PERIODS}}
        graph_cells[cid]={"n_source_components":1 if i<5 else 0}
    edges=[]
    for a,b in combinations(cells,2):
        edges.append({"a":a,"b":b,
                      "within_clipped_graph_connected":bool(
                          graph_cells[a]["n_source_components"] and graph_cells[b]["n_source_components"])})
    return {"country":country,"cells":cells},{"country":country,"cell_component_cover":graph_cells,
            "cell_pair_connection_descriptors":edges}


def fake_griis():
    return {"complete":True,"errors":{},"countries":{
        c:{"country":c,"full_source_pages_verified":True,
           "n_directly_matched_plant_backbone_species_keys":101,
           "gbif_backbone_plant_species_key_set":[str(100+i) for i in range(101)],
           "source_dataset_reference":"https://www.gbif.org/dataset/example",
           "n_unmatched_plant_name_usages":0}
        for c in set(COUNTRIES.values())-{"TW"}
    }}


class AlienRailJoinTests(unittest.TestCase):
    def test_national_introduction_status_is_country_specific(self):
        out=valid_griis(fake_griis())
        self.assertEqual(len(out),9)
        self.assertNotIn("TW",out)
        self.assertEqual(len(out["JP"]),101)

    def test_unlisted_species_are_not_called_native(self):
        subset=alien_source_species(plant_source(),set(str(100+i) for i in range(100)))
        self.assertEqual(set(subset),set(str(100+i) for i in range(5)))
        self.assertNotIn("999",subset)
        self.assertIsNone(alien_source_species(plant_source(),None))

    def test_introduced_population_uses_twenty_record_rarefaction(self):
        city,_=fake_city("jp_sendai","JP")
        obs,union=cell_observations(city,set(str(100+i) for i in range(101)),"2018_2025")
        self.assertEqual(len(obs),9)
        self.assertEqual(len(union),5)
        self.assertEqual(len(obs["jp_sendai_x0_y0"]["probs"]),5)
        self.assertTrue(all(0<p<=1 for p in obs["jp_sendai_x0_y0"]["probs"].values()))

    def test_country_cell_pairs_are_not_independent_railway_experiments(self):
        b,g=fake_city("jp_sendai","JP")
        ans=city_pair_analysis(b,g,set(str(100+i) for i in range(101)),"2018_2025")
        self.assertEqual(ans["n_source_eligible_introduced_plant_cells"],9)
        self.assertEqual(ans["n_dependent_pairwise_source_rows"],36)
        self.assertEqual(ans["n_connected_cell_pair_sources"],10)
        self.assertEqual(ans["n_griis_introduced_species_keys_observed_in_qualified_total_plant_cells"],5)
        self.assertIn("source_comparison_adequate",ans["all_pairs_comparison"])

    def test_real_source_incomplete_cannot_be_reported(self):
        x=fake_griis();x["complete"]=False
        with self.assertRaisesRegex(ValueError,"Incomplete GRIIS"):
            valid_griis(x)
        with self.assertRaises(ValueError):
            valid_griis({"complete":True,"errors":{},"countries":{}})

    def test_complete_fake_ten_source_region_triple_join(self):
        b={"source_complete":True,"regions":{}}
        g={"complete":True,"regions":{}}
        for name,country in COUNTRIES.items():
            x,y=fake_city(name,country)
            b["regions"][name]=x;g["regions"][name]=y
        out=analyze(b,g,fake_griis())
        self.assertEqual(len(out["country_lists"]),9)
        self.assertEqual(out["regions"]["tw_hsinchu"]["alien_status"],"unassessed")
        self.assertEqual(out["rarefaction_records_per_cell"],20)
        self.assertTrue(out["no_country_nonlisted_species_treated_as_native"])
        self.assertTrue(out["not_train_seed_transport_effect"])


if __name__=="__main__":
    unittest.main()

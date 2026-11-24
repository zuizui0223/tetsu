#!/usr/bin/env python3
"""v353: original global probability sample complete inclusion and count-only QA."""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_original_iNat_world_probability_cells_v353 import selected_cells,count_source_cell,summarize
SRC=Path(__file__).resolve().parents[1]/"data/derived/v351_global_land_climate_100km_probabilistic_source_verified.json"

class V353LandCellPlantOpportunityTests(unittest.TestCase):
    def setUp(self):
        self.frame=json.loads(SRC.read_text(encoding="utf8"))
        self.all=selected_cells(self.frame)

    def test_no_original_source_cells_dropped_or_reselected(self):
        self.assertEqual(len(self.all),32)
        self.assertEqual(len({r["cell"]["source_original_fixed_grid_index"] for r in self.all}),32)
        self.assertEqual(len({r["continent"] for r in self.all}),6)
        self.assertEqual(sum(r["macrotemp_stratum"]=="warm" for r in self.all),10)

    def test_API_count_is_only_source_count(self):
        r=self.all[0]
        client=lambda url:(json.dumps({"total_results":123,
                     "results":[{"id":123,"taxon":{"name":"hidden"}}]}).encode(),200)
        out=count_source_cell(r,client=client)
        self.assertEqual(out["source_api_2025_public_plants_bbox_count"],123)
        self.assertTrue(out["source_rectangular_bbox_not_exact_landcell_geometry"])
        self.assertFalse(out["source_original_observations_QUALITY_filtered_with_per_record_GPS"])
        self.assertNotIn("id",out)

    def test_valid_zero_observation_count_not_missing(self):
        out=count_source_cell(self.all[0],
          client=lambda u:(json.dumps({"total_results":0,"results":[]}).encode(),200))
        self.assertEqual(out["source_api_2025_public_plants_bbox_count"],0)

    def test_source_strata_support_keeps_zero_and_unavailable_separate(self):
        rows=[
          {"continent":"Asia","source_climate_stratum":"warm",
           "API":{"source_api_2025_public_plants_bbox_count":0}},
          {"continent":"Asia","source_climate_stratum":"warm",
           "API":{"source_api_2025_public_plants_bbox_count":None}},
          {"continent":"Africa","source_climate_stratum":"cool",
           "API":{"source_api_2025_public_plants_bbox_count":1000}},
        ]
        out=summarize(rows)
        a=next(x for x in out if x["continent"]=="Asia" and x["macrotemp_band"]=="warm")
        self.assertEqual(a["source_probability_draws"],2)
        self.assertEqual(a["original_source_available_n"],1)
        self.assertEqual(a["source_cells_zero_observed_plants_n"],1)
        self.assertEqual(a["source_unavailable_n"],1)
        self.assertEqual(len(out),18)

if __name__=="__main__":unittest.main()

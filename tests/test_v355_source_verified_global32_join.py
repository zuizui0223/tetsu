#!/usr/bin/env python3
"""v355 original Overture/iNat source aggregate content identity checks."""
import copy,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from validate_v355_global_rail_plant_source_join import validate
ROOT=Path(__file__).resolve().parents[1]
FRAME=ROOT/"data/derived/v351_global_land_climate_100km_probabilistic_source_verified.json"
JOINED=ROOT/"data/derived/v355_source_matched_global32_rail_plant_opportunity_verified.json"

class V355SourceAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame=json.loads(FRAME.read_text())
        cls.joined=json.loads(JOINED.read_text())

    def test_original_32_true_sources_reproduce_expected_values(self):
        x=validate(self.frame,self.joined)
        self.assertEqual(x["original_source_sample_cells"],32)
        self.assertEqual(x["mapped_rail_present_source_cells"],13)
        self.assertEqual(x["original_2025_iNat_rectangular_bbox_total_source_observations_not_unique"],43025)
        self.assertTrue(x["no_source_species_temperature_trait_or_rail_biological_effect"])

    def test_no_fake_source_rail_in_plant_observation_free_cell(self):
        j=copy.deepcopy(self.joined)
        j["source_cells_aggregate_only"][0]["original_mapped_rail_km_in_exact_100km_cell"]+=10
        with self.assertRaises(ValueError):
            validate(self.frame,j)

    def test_no_missing_observation_row_without_an_original_source(self):
        j=copy.deepcopy(self.joined)
        j["source_cells_aggregate_only"].pop()
        with self.assertRaises(ValueError):
            validate(self.frame,j)

    def test_original_cell_hash_and_reference_sample_immutable(self):
        j=copy.deepcopy(self.joined)
        j["source_cells_aggregate_only"][0]["source_cell_sha256"]="falsified"
        with self.assertRaises(ValueError):
            validate(self.frame,j)

    def test_no_source_plant_bbox_as_rail_near_track_fallacy(self):
        j=copy.deepcopy(self.joined)
        j["global_sample"]["source_citizen_counts_are_rectangular_BBOX_not_exact_100km_landcell_or_near_track"]=False
        with self.assertRaises(ValueError):
            validate(self.frame,j)

if __name__=="__main__":unittest.main()

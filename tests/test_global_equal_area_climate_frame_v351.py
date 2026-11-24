#!/usr/bin/env python3
"""Mathematical/source invariants: no actual world cell/sample effect is faked."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from generate_global_equal_area_climate_land_grid_v351 import (
    CELL_M,X_START,X_END,Y_START,Y_END,CONTINENTS,TEMP_BANDS,
    sample_exact_without_replacement,temp_band,summarize,SEED
)

class V351SourceFrameTests(unittest.TestCase):
    def setUp(self):
        self.cells=[{"cell_ix_iy":f"{i}:8","source_BIO1_C":1.0,
                    "cell_center_lat":35.0,"cell_center_lon":139.0,
                    "continent":"Asia","climate_band":"cool"} for i in range(11)]

    def test_equal_area_frozen_grid_and_six_continents(self):
        self.assertEqual(CELL_M,100000)
        self.assertEqual((X_END-X_START)//CELL_M,340)
        self.assertEqual((Y_END-Y_START)//CELL_M,146)
        self.assertEqual(len(CONTINENTS),6)
        self.assertEqual(len(TEMP_BANDS),3)

    def test_macroclimate_band_boundaries(self):
        self.assertEqual(temp_band(-3),"cool")
        self.assertEqual(temp_band(7.9999),"cool")
        self.assertEqual(temp_band(8),"temperate")
        self.assertEqual(temp_band(19.999),"temperate")
        self.assertEqual(temp_band(20),"warm")
        self.assertEqual(temp_band(37),"warm")
        self.assertIsNone(temp_band(float("nan")))

    def test_exact_finite_frame_source_inclusion_probability(self):
        selected,pi=sample_exact_without_replacement(self.cells,2,SEED+"|Asia|cool")
        self.assertEqual(len(selected),2)
        self.assertEqual(pi,2/11)
        self.assertEqual(len({x["cell_ix_iy"] for x in selected}),2)
        self.assertEqual(selected,sample_exact_without_replacement(self.cells,2,SEED+"|Asia|cool")[0])

    def test_zero_size_and_complete_small_strata_are_explicit(self):
        self.assertEqual(sample_exact_without_replacement([],2,SEED),([],None))
        selected,pi=sample_exact_without_replacement(self.cells[:1],2,SEED)
        self.assertEqual(len(selected),1)
        self.assertEqual(pi,1.0)

    def test_duplicate_cells_fail_closed(self):
        with self.assertRaises(ValueError):
            sample_exact_without_replacement(self.cells[:2]*2,2,SEED)

    def test_record_count_18_strata_and_no_biological_effects(self):
        population={("Asia","cool"):self.cells}
        qa={"source_climate_supported_land_cells":11}
        x=summarize(population,qa,{"source":"synthetic ONLY"})
        self.assertEqual(x["source_strata_n"],18)
        self.assertEqual(x["n_sampled_unique_land_cells"],2)
        self.assertTrue(x["not_an_actual_global_rail_length_probability_frame"])
        self.assertTrue(x["ecological_rail_plant_thermal_effects_not_estimated"])

if __name__=="__main__":unittest.main()

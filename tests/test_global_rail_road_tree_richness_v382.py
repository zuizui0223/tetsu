import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,within_continent_direction
class TestMacro(unittest.TestCase):
    def test_matched_land_support_and_interaction(self):
        cells=[dict(key=(0,0),log_rail_km=2.)]
        rd={(0,0):{"source_GRIP4_land_fraction_of_cell":.7,"source_GRIP4_road_coverage_of_source_land":1.,"source_GRIP4_observed_road_density_m_per_km2_land":0.}}
        tr={(0,0):{"source_valid_area_fraction_of_grid_approx":1.,"source_tree_consensus_area_m2":3e9}}
        rows,qa=enrich_cells(cells,rd,tr,100000,.5)
        self.assertAlmostEqual(rows[0]["rail_tree"],.6)
        self.assertEqual(rows[0]["road_tree"],0.)
        self.assertEqual(enrich_cells(cells,rd,tr,100000,.8)[0],[])
        rd[(0,0)]["source_GRIP4_road_coverage_of_source_land"]=.8
        self.assertEqual(enrich_cells(cells,rd,tr,100000,.5)[0],[])
    def test_direction_is_descriptive_and_controls_nuisance(self):
        rows=[{"continent":"Asia","block":(i//10,0),"z":i%3,"interaction":float(i),"y":2*(i%3)-.1*i} for i in range(120)]
        self.assertAlmostEqual(within_continent_direction(rows,"y",("z","interaction"),"interaction")["estimate"],-.1)
if __name__=="__main__":unittest.main()

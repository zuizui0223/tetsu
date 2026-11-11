"""Synthetic quality/provenance tests of FRA train-count source, no seed ecology."""
import sys,unittest
from pathlib import Path
from shapely.geometry import box
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US_FRA_crossing_train_counts_v75 import (
    BLOCKS, FIELDS, numeric_count,year_of_count,source_categories,
    assign_cell,summarise_block,frozen_grid,query_params
)

class FRATrainCountSourceTests(unittest.TestCase):
    def test_original_US_blocks_are_preselected_and_have_16_cells(self):
        self.assertEqual(BLOCKS,("us_midwest","us_southeast","us_pacific"))
        cells,bounds=frozen_grid("us_midwest")
        self.assertEqual(len(cells),16)
        self.assertTrue(bounds[0]<-88 and bounds[2]>-91)
    def test_missing_and_lt_one_traffic_are_not_invented_zero(self):
        self.assertIsNone(numeric_count(None,500))
        self.assertIsNone(numeric_count(-1,500))
        self.assertIsNone(numeric_count("0",500))
        self.assertEqual(numeric_count(0,500),0)
    def test_year_is_report_year_not_train_incidence(self):
        self.assertEqual(year_of_count("2024"),2024)
        self.assertIsNone(year_of_count(None))
        self.assertIsNone(year_of_count("2017a"))
        self.assertIsNone(year_of_count("2030"))
    def test_freight_service_category_is_not_freight_train_volume(self):
        self.assertEqual(source_categories("11,12"),("11","12"))
        self.assertEqual(source_categories("1211"),())
        self.assertEqual(source_categories(None),())
    def test_point_membership_and_source_precision_gate(self):
        cells={"test":box(-90,42,-89,43)}
        self.assertEqual(assign_cell(cells,{"Longitude":-89.5,"Latitude":42.5},
                             {"x":-89.5,"y":42.5}),( "test",None))
        cell,error=assign_cell(cells,{"Longitude":-89.5,"Latitude":42.5},
                             {"x":-80,"y":42.5})
        self.assertIsNone(cell)
        self.assertEqual(error,"coordinate_source_disagreement")
    def test_aggregation_denominators_not_rail_biology(self):
        cells={"c":box(-90,42,-89,43)}
        a={"Crossing_ID":"0001","Longitude":-89.5,"Latitude":42.5,
           "Trains_Per_Week_Captured_Year":"2024",
           "Total_Daylight_Thru_Trains":9,"Total_Nighttime_Thru_Trains":2,
           "Typical_Maximum_Speed_Over_Crossing":55,
           "Type_Of_Train_Service_IDs":"11,12","Crossing_Closed":"No"}
        b=dict(a,Crossing_ID="0002",Trains_Per_Week_Captured_Year=None,
               Total_Daylight_Thru_Trains=None)
        r=summarise_block(cells,[{"attributes":a,"geometry":{"x":-89.5,"y":42.5}},
                                 {"attributes":b,"geometry":{"x":-89.5,"y":42.5}}])
        self.assertEqual(r["n_in_cell_source_crossing_records"],2)
        self.assertEqual(r["by_original_50km_cells"]["c"]["n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields"],1)
        self.assertEqual(r["by_original_50km_cells"]["c"]["median_reported_total_daily_through_trains_at_eligible_crossings"],11)
        self.assertEqual(r["by_original_50km_cells"]["c"]["n_crossings_reporting_freight_service_type_NOT_freight_train_counts"],2)
        self.assertTrue(r["DO_NOT_TREAT_CROSSINGS_AS_INDEPENDENT_RAIL_SEGMENTS"])
    def test_query_geometry_predeclared_not_selected_by_rail_outcomes(self):
        bounds=(-90,42,-89,43)
        params=query_params(bounds)
        self.assertEqual(params["where"],"1=1")
        self.assertEqual(params["geometryType"],"esriGeometryEnvelope")
        self.assertIn("Crossing_ID",FIELDS)
        self.assertNotIn("RrCont",FIELDS)
if __name__=="__main__":unittest.main()

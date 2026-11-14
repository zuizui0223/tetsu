"""Synthetic external-source QA: GBIF stateProvince name versus original
lat/lon in Census 2025 state boundary, not plant movement or ecological
treatment. Test paths never export exact record coordinates.
"""
import io,sys,unittest
from pathlib import Path
from unittest.mock import patch
from shapely.geometry import box
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_US_GBIF_stateProvince_original_coords_v90 import (
    PILOT_FIPS,PERIODS,REQUEST_OFFSETS,ROWS_PER_PAGE,
    CENSUS_SHAPE_URL,chosen_offsets,occurrence_params,
    bounded_source_points,assess_geographic_points,census_original_polygons)

def record(i,lon,lat,taxon):
    return {
      "gbifID":str(i),"key":i,
      "year":2023,"decimalLongitude":lon,"decimalLatitude":lat,
      "stateProvince":"Wyoming",
      "speciesKey":taxon,
      "datasetKey":"dataset_%s"%(i%3),
      "occurrenceStatus":"PRESENT",
    }

class GBIFOriginalStateGeoContracts(unittest.TestCase):
    def test_eight_preselected_and_disjoint_years_and_independent_Census(self):
        self.assertEqual(PILOT_FIPS,("06","12","17","19","36","37","48","53"))
        self.assertEqual(PERIODS,{"before_2022":"2018,2021","after_2022":"2023,2025"})
        self.assertEqual(REQUEST_OFFSETS,(0,1000))
        self.assertEqual(ROWS_PER_PAGE,100)
        self.assertIn("census.gov",CENSUS_SHAPE_URL)
    def test_fixed_source_bounded_page_offsets_are_NOT_uniform_global_random(self):
        self.assertEqual(chosen_offsets(0),[])
        self.assertEqual(chosen_offsets(70),[0])
        self.assertEqual(chosen_offsets(101),[0,1])
        self.assertEqual(chosen_offsets(1000000),[0,1000])
        with self.assertRaisesRegex(ValueError,"Malformed"):
            chosen_offsets(-1)
    def test_actual_GBIF_source_query_no_source_PRESENCE_enum(self):
        v=occurrence_params("New York","after_2022",50000)
        self.assertEqual(v["year"],"2023,2025")
        self.assertEqual(v["stateProvince"],"New York")
        self.assertEqual(v["limit"],100)
        self.assertNotIn("occurrenceStatus",v)
        with self.assertRaisesRegex(ValueError,"Predeclared"):
            occurrence_params("New York","2022",0)
    def test_geographical_area_conserves_all_plants_and_GRIIS_unknowns(self):
        records=[
            record(1,-85,45,"1001"),
            record(2,-95,45,"1002"),
            record(3,-90.01,45,"9999"),
            record(4,None,None,"1001"),
        ]
        s=assess_geographic_points(records,box(-90,40,-80,50),"Wyoming",{"1001","1002"})
        a=s["all_Plantae"];g=s["USA_country_GRIIS_listed_records"]
        self.assertEqual(a["n_raw_GBIF_source_looked_up"],4)
        self.assertEqual(g["n_raw_GBIF_source_looked_up"],3)
        self.assertEqual(a["source_coordinate_QA_category_counts"][
            "coordinate_strictly_within_source_state_polygon"],1)
        self.assertEqual(a["source_coordinate_QA_category_counts"][
            "coordinate_more_than_5km_outside_named_state_boundary"],1)
        self.assertEqual(a["source_coordinate_QA_category_counts"][
            "invalid_or_missing_reported_coordinates"],1)
        self.assertTrue(s["no_exact_coordinate_or_GBIF_record_ID_in_artifact"])
    def test_bounded_2023_state_probe_reports_exact_source_page_rows_not_plant_status(self):
        records=[record(i,-85,45,"1001") for i in range(1,81)]
        with patch("audit_US_GBIF_stateProvince_original_coords_v90.source_get",
                   return_value={"count":80,"results":records}):
            got,audit=bounded_source_points(None,"06","California","after_2022")
        self.assertEqual(len(got),80)
        self.assertEqual(audit["original_source_offsets_probed"],[0])
        self.assertTrue(audit["NOT_random_probability_sample_of_all_GBIF_state_records"])
    def test_country_shape_rejects_malformed_source_zip(self):
        with self.assertRaisesRegex(ValueError,"missing/truncated"):
            census_original_polygons(b"not a census file")
if __name__=="__main__":unittest.main()

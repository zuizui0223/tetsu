import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from audit_Czech_35SR70_station_local_CzechVeg_nontransport_georeference_v239 import (
    spatial_opportunity,original_reference_geo
)
class TestCzechStationSourcePrecisionOpportunity(unittest.TestCase):
    def synthetic(self):
        station=[];plots={}
        for i in range(35):
            lat=49.0+i*.032;lon=13.2+i*.05
            station.append({
              "original_locality_code":f"rail{i}",
              "source_section":str(i%5),
              "matched_czech_2026_source_station_centroid_WGS84":{"lat":lat,"lon":lon}
            })
            # 6 exact independent 2019 plots, two source papers,
            # <100m off the true registry representative point
            for k in range(6):
                plots[f"ref_{i}_{k}"]={
                  "lat":lat+0.00005*(k+1),
                  "lon":lon+0.00005*(k+1),
                  "precision_m":90,
                  "year":2019,
                  "study":("forest-paper" if k%2 else "grassland-paper"),
                  "source_region":str(i%5),"source_area_class":"small"}
        return station,plots
    def test_35_strict_source_stations_all_have_local_2015_reference(self):
        station,plots=self.synthetic()
        scores,qa=spatial_opportunity(station,plots)
        n=scores["year_since_2015_precision_250m_radius_5km"]
        self.assertEqual(n["n_original_matched_author_stations"],35)
        self.assertEqual(n["original_stations_ge5_local_reference_plots_from_ge2_publications"],35)
        self.assertEqual(n["original_stations_ge10_local_reference_plots_from_ge3_publications"],0)
    def test_missing_field_survey_GPS_never_imputed_for_4_rail_stations(self):
        station,plots=self.synthetic()
        station.append({"original_locality_code":"unmatched",
                        "source_section":"1",
                        "matched_czech_2026_source_station_centroid_WGS84":None})
        scores,qa=spatial_opportunity(station,plots)
        self.assertEqual(qa["n_original_Czech_rail_stations_matched_registry_centroid"],35)
    def test_true_source_position_precision_gate_excludes_coarse_control(self):
        station,plots=self.synthetic()
        for p in plots.values():p["precision_m"]=1500
        scores,qa=spatial_opportunity(station,plots)
        self.assertEqual(scores["year_since_2015_precision_1000m_radius_25km"][
           "original_stations_ge5_local_reference_plots_from_ge2_publications"],0)
    def test_true_suspicious_zero_metres_precision_is_excluded_only_in_sensitivity(self):
        station,plots=self.synthetic()
        for p in plots.values():p["precision_m"]=0
        raw,qa=spatial_opportunity(station,plots)
        self.assertEqual(raw["year_since_2015_precision_250m_radius_5km"][
            "original_stations_ge5_local_reference_plots_from_ge2_publications"],35)
        strict,qa=spatial_opportunity(station,{k:v for k,v in plots.items() if v["precision_m"]>0})
        self.assertEqual(strict["year_since_2015_precision_250m_radius_5km"][
            "original_stations_ge5_local_reference_plots_from_ge2_publications"],0)
    def test_true_original_source_year_gate_excludes_historic(self):
        station,plots=self.synthetic()
        for p in plots.values():p["year"]=1977
        scores,qa=spatial_opportunity(station,plots)
        self.assertEqual(scores["year_since_2000_precision_1000m_radius_25km"][
           "original_stations_ge5_local_reference_plots_from_ge2_publications"],0)
    def test_location_support_not_rail_vs_natural_effect(self):
        station,plots=self.synthetic()
        scores,qa=spatial_opportunity(station,plots)
        self.assertTrue(qa["source_2026_rail_centroids_NOT_botanical_subplot_locations"])
        self.assertTrue(qa["plot_and_station_exact_coordinates_not_published"])
if __name__=="__main__":unittest.main()

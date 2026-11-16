"""v15.3 synthetic correctness tests; no synthetic source station accepted as actual field data."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Czech_2024_S2_to_SR70_original_station_points_v153 import (
    normalized_name,match_czech_stations,MAX_SAME_PLACE_M,
)
def source_station(name="Brno-Maloměřice",code="BM"):
    return {"original_source_s2_col1_label":name,
            "original_locality_code":code,"source_section":"II",
            "source_s2_sheet_row":6}
def reg(name,lat,lon,code):
    return {"name":name,"lat":lat,"lon":lon,"sr70_code":code,
            "match_key":normalized_name(name)}

class CzechOriginalStationIdentitySourceTests(unittest.TestCase):
    def test_exact_accents_spaces_case_dash_normalized(self):
        self.assertEqual(normalized_name("Ústí nad Labem-Střekov"),
                         normalized_name("USTI NAD LABEM STREKOV"))
        self.assertEqual(normalized_name("Žďár nad Sázavou"),
                         "zdar nad sazavou")
        self.assertNotEqual(normalized_name("Praha-Bubny"),
                            normalized_name("Praha-Bubny Vltavská"))

    def test_different_stations_not_matched_by_city(self):
        ss=[source_station("Brno-Maloměřice")]
        rs=[reg("Brno-Maloměřice st.3",49.224,16.650,"380089"),
            reg("Brno dolní nádraží",49.182,16.615,"380154")]
        m=match_czech_stations(ss,rs)
        self.assertEqual(m[0]["match_status"],"no_exact_full_name_match")
        self.assertIsNone(m[0]["matched_czech_2026_source_station_centroid_WGS84"])

    def test_independent_sr70_codes_same_physical_place_collapse(self):
        m=match_czech_stations([source_station()],[
            reg("Brno-Maloměřice",49.21686,16.647632,"333559"),
            reg("Brno-Maloměřice",49.21686,16.647632,"833343")])
        self.assertEqual(m[0]["match_status"],"unique_place_by_identical_full_station_name")
        self.assertEqual(m[0]["source_registry_full_name_candidates"],2)
        self.assertEqual(m[0]["distinct_150m_place_clusters"],1)
        self.assertEqual(m[0]["matched_czech_2026_source_station_centroid_WGS84"][
           "n_sr70_registry_codes_collapsed_same_place"],2)

    def test_far_same_name_candidates_never_randomly_assigned(self):
        m=match_czech_stations([source_station("Czech Station")],[
           reg("Czech Station",49.10,16.1,"111111"),
           reg("Czech Station",50.10,16.1,"222222")])
        self.assertEqual(m[0]["match_status"],"ambiguous_different_station_locations")
        self.assertIsNone(m[0]["matched_czech_2026_source_station_centroid_WGS84"])

    def test_a_same_name_station_just_above_gate_remains_ambiguous(self):
        m=match_czech_stations([source_station("A")],[
           reg("A",49.1,16.1,"111111"),
           reg("A",49.102,16.1,"222222")])
        self.assertEqual(m[0]["match_status"],"ambiguous_different_station_locations")

    def test_empty_publisher_original_name_never_imputed_from_code(self):
        m=match_czech_stations([source_station("","BM")],[reg("BM",49,16,"123456")])
        self.assertEqual(m[0]["match_status"],"no_exact_full_name_match")
        self.assertEqual(m[0]["original_source_full_name_normalized"],"")

if __name__=="__main__":unittest.main()

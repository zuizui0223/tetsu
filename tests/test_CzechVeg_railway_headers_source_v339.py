"""v339 CzechVeg full-source habitat candidate-label strict tests."""
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_CzechVeg_railway_headers_source_v339 as m

class CzechVegHabitatGate(unittest.TestCase):
    def test_railroad_habitat_words(self):
        for text in ("železniční násep","kolejiště","vlečka","u nádraží","železnice"):
            status,kinds=m.classify(text,"")
            self.assertEqual(status,"rail_habitat_text_candidate")
            self.assertTrue(kinds)
    def test_locality_railway_landmark_not_habitat(self):
        status,kinds=m.classify("acidofilní suchá louka","u železniční trati")
        self.assertEqual(status,"rail_locality_landmark_only")
        self.assertEqual(kinds,[])
    def test_bus_station_never_assumed_railway(self):
        status,kinds=m.classify("autobusové nádraží","")
        self.assertEqual(status,"ambiguous_other_transit_habitat")
    def test_ordinary_track_word_not_assumed_rail(self):
        status,kinds=m.classify("průsek","kolej sportovního klubu")
        self.assertEqual(status,"no_rail_text_evidence")
    def test_unicode_habitat_normalization(self):
        self.assertIn("zelez",m.ascii_cz("Železniční dráha"))
    def test_decimal_metadata_missing_not_zero(self):
        self.assertIsNone(m.good_number(""))
        self.assertIsNone(m.good_number("NA"))
        self.assertIsNone(m.good_number(float("nan")))
        self.assertEqual(m.good_number("49.5"),49.5)
    def test_stable_documented_source_count(self):
        self.assertEqual(m.EXPECTED_PLOTS,117739)
        self.assertEqual(m.EXPECTED_MD5,"0bc920462d1f9eac473e7b2a231e6782")
    def test_railway_in_bibliography_not_accepted(self):
        status,kinds=m.classify("ruderal vegetation","urban fringe")
        self.assertEqual(status,"no_rail_text_evidence")

if __name__=="__main__":unittest.main()

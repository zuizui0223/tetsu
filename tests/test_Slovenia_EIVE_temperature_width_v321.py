"""v321 synthetic tests for original codebook and independent EIVE trait joins."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Slovenia_EIVE_temperature_width_v321 as m

class CrossSourceTests(unittest.TestCase):
    def test_casefold_taxon_name_no_synonym_inference(self):
        self.assertEqual(m.normal_name("  Poa   annua "),"poa annua")
        self.assertNotEqual(m.normal_name("Erigeron canadensis"),
                            m.normal_name("Conyza canadensis"))
    def test_numeric_strict(self):
        self.assertIsNone(m.float_finite("x"))
        self.assertIsNone(m.float_finite("inf"))
        self.assertEqual(m.float_finite("5.5"),5.5)
    def test_source_system_count_gate(self):
        source={"A":"poa annua","B":"galium aparine"}
        trait={"poa annua":{"position":7,"width":3,"systems":3},
               "galium aparine":{"position":5,"width":6,"systems":1}}
        self.assertEqual(m.score({"A","B"},source,trait,"position",1),(6,.9999999999999999) if False else (6,1.0))
        self.assertEqual(m.score({"A","B"},source,trait,"position",3),(7,.5))
    def test_one_location_two_track_sides(self):
        pairs=[{"site":"s1","section":"r1","side":"L","BS":{"A"},"BT":{"B"}},
               {"site":"s1","section":"r1","side":"R","BS":{"A"},"BT":{"B"}},
               {"site":"s2","section":"r2","side":"L","BS":{"A"},"BT":{"B"}}]
        names={"A":"poa annua","B":"galium aparine"}
        trait={"poa annua":{"position":7,"width":2,"systems":5},
               "galium aparine":{"position":5,"width":6,"systems":5}}
        result=m.compare(pairs,names,trait,"position",3,0)
        self.assertEqual(result["eligible_original_track_sides"],3)
        self.assertEqual(result["eligible_distinct_field_locations"],2)
        self.assertEqual(result["field_location_mean_BS_minus_BT"],2)
        result2=m.compare(pairs,names,trait,"width",3,0)
        self.assertEqual(result2["field_location_mean_BS_minus_BT"],-4)
    def test_excluding_unmapped_values_reduces_coverage(self):
        names={"A":"poa annua","B":"unmatched plant"}
        trait={"poa annua":{"position":7,"width":3,"systems":2}}
        self.assertEqual(m.score({"A","B"},names,trait,"width",1),(3,0.5))

if __name__=="__main__":unittest.main()

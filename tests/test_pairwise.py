import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.pairwise import jaccard, read_incidence, prepare_pairwise


class PairwiseTests(unittest.TestCase):
    def test_jaccard(self):
        self.assertAlmostEqual(jaccard({"a", "b"}, {"b", "c"}), 1/3)
        self.assertIsNone(jaccard(set(), {"x"}))

    def test_duplicate_taxon_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            file = Path(tmp) / "inc.csv"
            file.write_text("""region_id,taxon,status
A,x,alien
A,x,alien
""")
            with self.assertRaises(ValueError):
                read_incidence(file)

    def test_contradictory_status_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            file = Path(tmp) / "inc.csv"
            file.write_text("""region_id,taxon,status
A,x,alien
B,x,native
""")
            with self.assertRaises(ValueError):
                read_incidence(file)

    def test_connected_and_disconnected(self):
        with tempfile.TemporaryDirectory() as tmp:
            inc = Path(tmp) / "inc.csv"
            inc.write_text("""region_id,taxon,status
A,x,alien
A,y,native
B,x,alien
B,z,native
C,x,alien
C,y,native
""")
            pair = Path(tmp) / "pair.csv"
            pair.write_text("""region_a,region_b,euclid_km,road_km,rail_km,climate_distance,urban_difference,effort_a,effort_b
A,B,100,120,130,0.4,0.1,100,90
A,C,100,120,,0.4,0.1,100,90
""")
            result = prepare_pairwise(read_incidence(inc), pair)
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0]["rail_connected"], 1)
            self.assertEqual(result[1]["rail_connected"], 0)
            self.assertEqual(result[0]["alien_minus_native"], 1.0)
            self.assertEqual(result[1]["alien_minus_native"], 0.0)

    def test_duplicate_unordered_pair_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            inc = Path(tmp) / "inc.csv"
            inc.write_text("""region_id,taxon,status
A,x,alien
B,x,alien
""")
            pair = Path(tmp) / "pair.csv"
            pair.write_text("""region_a,region_b,euclid_km,road_km,rail_km,climate_distance,urban_difference,effort_a,effort_b
A,B,10,11,12,0.1,0,4,4
B,A,10,11,12,0.1,0,4,4
""")
            with self.assertRaises(ValueError):
                prepare_pairwise(read_incidence(inc), pair)


if __name__ == "__main__":
    unittest.main()

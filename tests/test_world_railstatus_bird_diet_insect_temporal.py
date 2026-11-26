import unittest
from analyze_world_railstatus_bird_diet_insect_temporal import annual_model,estimate_scores,verify_previous
from scan_world_railstatus_annual_insect_cells import candidate_cells
class TestTemporal(unittest.TestCase):
    def test_empty(self):
        z=annual_model({y:[] for y in ("2023","2024","2025")})
        self.assertEqual(z["n_common_cells"],0)
        self.assertIsNone(z["models"]["flexible"]["earlier_mean"])
    def test_replay(self):
        a={"status":"insufficient support or rank","n_cells":5}
        verify_previous(a,a)
        with self.assertRaises(ValueError):verify_previous(a,dict(a,n_cells=6))
    def test_candidates(self):
        p={"a":[0,0],"b":[1,0],"full_a":[60,10,20],"full_b":[60,10,20],"split_counts_a":[[30,5,10],[30,5,10]],"split_counts_b":[[30,5,10],[30,5,10]],"split_similarity":{}}
        src={"pairs":{"200000":{"insects":[p]}}}
        diet={"grids":{"200000":{"2023":[{"key":[0,0]}],"2024":[{"key":[1,0]}],"2025":[{"key":[0,0]},{"key":[1,0]}]}}}}
        c,old=candidate_cells(src,diet)
        self.assertEqual(c["2023"],{(0,0)});self.assertEqual(len(old),2)
if __name__=="__main__":unittest.main()

import unittest
from scan_world_railstatus_observer_annual_beta import candidate_support
from analyze_world_railstatus_observer_annual_connectivity import cells_from_pairs,subset
class TestAnnualSupport(unittest.TestCase):
    def test_frozen_candidates(self):
        a={"a":[0,0],"b":[1,0],"split_similarity":{}}
        b={"a":[1,0],"b":[2,0]}
        src={"pairs":{"200000":{"plants":[a,b],"birds":[a]}}}
        p,c=candidate_support(src)
        self.assertEqual(len(p["plants"]),1)
        self.assertEqual(c["birds"],{(0,0),(1,0)})
    def test_counts(self):
        a={"a":[0,0],"b":[1,0],"full_a":[50,5,20],"full_b":[100,8,30]}
        b={"a":[1,0],"b":[2,0],"full_a":[100,8,30],"full_b":[60,6,25]}
        self.assertEqual(len(cells_from_pairs([a,b])),3)
        with self.assertRaises(ValueError):cells_from_pairs([a,dict(b,full_a=[101,8,30])])
    def test_exact_subset(self):
        r={"a":(0,0),"b":(1,0)}
        self.assertEqual(subset([r],{((0,0),(1,0))}),[r])
        with self.assertRaises(ValueError):subset([r],{((1,0),(2,0))})
if __name__=="__main__":unittest.main()

import unittest
from analyze_world_railstatus_observer_split_connectivity import symmetric_effort,split_pair,adapt_pair
class TestSplitConnectivity(unittest.TestCase):
    def test_symmetric_effort(self):
        a=[(30,5,12),(40,8,14)];b=[(50,9,20),(25,5,10)]
        self.assertEqual(symmetric_effort(a,b),symmetric_effort(a[::-1],b))
        self.assertEqual(symmetric_effort(a,b),symmetric_effort(b,a))
    def test_paired_outcomes(self):
        p={"split_counts_a":[(30,5,12),(40,8,14)],"split_counts_b":[(50,9,20),(25,5,10)],
           "split_similarity":{k:{"jaccard":v,"simpson_similarity":v} for k,v in (("00",.4),("01",.1),("10",.3),("11",.6))},
           "observer_full":{"jaccard":.1},"same_side_observer_overlap":{"0":{"jaccard":.2},"1":{"jaccard":.4}}}
        r=split_pair({},p)
        self.assertAlmostEqual(r["same_jaccard"],.5)
        self.assertAlmostEqual(r["cross_jaccard"],.2)
        self.assertAlmostEqual(r["same_minus_cross_jaccard"],.3)
        self.assertIsNone(split_pair({},{}))
if __name__=="__main__":unittest.main()

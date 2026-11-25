import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_naturalized_quality_effort_v390 import common_support,shared_effort_comparison
class TestV390(unittest.TestCase):
    def test_common_support_requires_both_arms_in_same_geography(self):
        rows=[]
        for region,block,levels in [("A",(0,0),[0,0,1,1]),("A",(1,0),[1,1,1,1]),
                                    ("B",(0,0),[0,1,1]),("B",(1,0),[0,0,0,0])]:
            for v in levels:rows.append(dict(region=region,continent="Asia",block=block,
                                            observations=200 if v else 100,observers=20 if v else 10))
        r,k=common_support(rows)
        self.assertEqual(len(r),4);self.assertEqual(len(k),1)
        self.assertEqual(shared_effort_comparison(rows)["status"],"insufficient shared geographic support")
if __name__=="__main__":unittest.main()

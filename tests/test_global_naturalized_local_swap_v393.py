import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_naturalized_local_swap_v393 import supported_local,match,FEATURES
class TestV393(unittest.TestCase):
    def test_region_block_support_and_matching(self):
        rows=[]
        for k,n in [("A",4),("B",2),("C",4)]:
            for i in range(n):
                r={x:1. for x in FEATURES}
                r.update(region=k,continent="Asia",block=(0,0),
                         log_rail_km=float(i if k!="C" else 0))
                rows.append(r)
        s=supported_local(rows);self.assertEqual(len(s),4)
        self.assertEqual(match(s),[(0,1),(2,3)])
        self.assertTrue(all(rows[i]["region"]==rows[j]["region"] for i,j in match(rows)))
if __name__=="__main__":unittest.main()

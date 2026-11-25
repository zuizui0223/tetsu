import sys,unittest,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_world_railstatus_bird_hwi_paired_v418 import paired_average
class TestPairedHWI(unittest.TestCase):
    def test_identical_scores_cancel(self):
        def m(v):
            return {y:{"keys":[("x",1),("x",2)],"effects":{"active":{"estimate":v,"scores":np.array([.1,-.1])}}} for y in ("2023","2024")}
        x=paired_average([m(1.),m(2.),m(5.)],("2023","2024"),"active",(0.,-1.,1.))
        self.assertEqual(x["estimate"],3.)
        self.assertEqual(x["interval95_approx"],[3.,3.])
    def test_mismatch_fails(self):
        a={"2023":{"keys":[("x",1),("x",2)]}}
        b={"2023":{"keys":[("x",1),("x",3)]}}
        with self.assertRaises(ValueError):paired_average([a,b],("2023",),"active",(-1.,1.))
if __name__=="__main__":unittest.main()

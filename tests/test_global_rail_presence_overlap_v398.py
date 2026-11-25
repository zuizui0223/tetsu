import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from audit_global_rail_presence_overlap_v398 import graph,make_weights,ess,MATCH
class TestV398(unittest.TestCase):
    def test_reuse_is_not_extra_effective_sample(self):
        r=[dict(region="A",continent="Asia",block=(0,0),rail_present=p) for p in (0,0,1)]
        g=graph(r,np.zeros((3,len(MATCH))))
        self.assertEqual(g,{0:[2],1:[2]})
        a,b=make_weights(r,g)
        np.testing.assert_allclose(a,[.5,.5,0]);np.testing.assert_allclose(b,[0,0,1])
        self.assertAlmostEqual(ess(a),2);self.assertAlmostEqual(ess(b),1)
    def test_no_cross_stratum_or_caliper(self):
        r=[dict(region=k,continent="Asia",block=(0,0),rail_present=p) for k,p in (("A",0),("B",1),("A",1))]
        z=np.zeros((3,len(MATCH)));z[2,0]=.8
        self.assertEqual(graph(r,z),{})
if __name__=="__main__":unittest.main()

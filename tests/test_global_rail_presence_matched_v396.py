import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from analyze_global_rail_presence_matched_v396 import optimal_pairs,pair_weights,MATCH
class TestV396(unittest.TestCase):
    def test_max_cardinality_before_distance(self):
        # zero 0 can use either positive; zero 1 can only use positive 2.
        rows=[dict(region="A",continent="Asia",block=(0,0),key=(i,0),rail_present=int(i>=2)) for i in range(4)]
        z=np.zeros((4,len(MATCH)));z[:,0]=[.4,0.,.1,1.]
        self.assertEqual(set(optimal_pairs(rows,z)),{(0,3),(1,2)})
    def test_cross_stratum_forbidden(self):
        rows=[dict(region=k,continent="Asia",block=(0,0),key=(i,0),rail_present=i) for i,k in enumerate(("A","B"))]
        self.assertEqual(optimal_pairs(rows,np.zeros((2,len(MATCH)))),[])
    def test_equal_block_mass(self):
        np.testing.assert_allclose(pair_weights([("Asia",(0,0)),("Asia",(0,0)),("Asia",(1,0))]),[.25,.25,.5])
if __name__=="__main__":unittest.main()

import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_naturalized_joint_population_local_v394 import match,MATCH,groups,supported
class TestV394(unittest.TestCase):
    def test_population_matching_changes_support(self):
        rows=[]
        for i,p in enumerate([0.,0.,10.,10.]):
            r={k:1. for k in MATCH}
            r.update(region="A",continent="Asia",block=(0,0),log_population_density=p)
            rows.append(r)
        self.assertEqual(match(rows),[(0,1),(2,3)])
        rows[1]["region"]="B"
        self.assertEqual(match(rows),[(2,3)])
    def test_small_stratum_is_not_estimated(self):
        r=dict(region="A",continent="Asia",block=(0,0))
        self.assertFalse(supported([r]*4,("x",)))
if __name__=="__main__":unittest.main()

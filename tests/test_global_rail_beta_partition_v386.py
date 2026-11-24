import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_rail_beta_partition_v386 import partition,direction,FIELDS,RAIL
class TestV386(unittest.TestCase):
    def test_partition(self):
        self.assertEqual(partition(2,2,2,2),(0.,0.,0.))
        self.assertEqual(partition(0,2,3,5),(1.,1.,0.))
        self.assertEqual(partition(2,2,4,4),(.5,0.,.5))
        a,t,n=partition(2,3,5,6)
        self.assertAlmostEqual(a,t+n)
        self.assertAlmostEqual(t,.5)
        self.assertAlmostEqual(a,2/3)
        with self.assertRaises(ValueError):partition(3,2,4,3)
    def test_conditional_coefficient_additivity(self):
        rows=[]
        for i in range(120):
            x=(i%11)/20;y=(i%7)/20
            r={k:0. for k in FIELDS}
            r.update(continent="Asia",block=(i//5,0),log_mean_rail_km=x,absdiff_log_rail_km=y,
                     beta_jtu=.2+x*.2,beta_jne=.1+y*.1,beta_jac=.3+x*.2+y*.1)
            rows.append(r)
        v=direction(rows)
        self.assertAlmostEqual(v[RAIL[0]]["beta_jtu"],.2)
        self.assertAlmostEqual(v[RAIL[0]]["beta_jne"],0.)
        self.assertAlmostEqual(v[RAIL[1]]["beta_jac"],.1)
if __name__=="__main__":unittest.main()

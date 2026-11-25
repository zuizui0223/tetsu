import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_annual_change_v410 import make_change,RESP,YEARS
class TestV410(unittest.TestCase):
    def test_paired_change(self):
        rows={}
        for i,y in enumerate(YEARS):
            r={k:float(i) for k in RESP}
            r.update(key=(1,2),log_observations=2.+i,log_rail_km=1.)
            rows[y]=[r]
        a,f=make_change(rows,("log_rail_km","log_observations"),"2025_minus_earlier_mean")
        self.assertEqual(a[0]["log_observations"],2.5)
        self.assertEqual(a[0]["delta_log_observations"],1.5)
        self.assertEqual(a[0][RESP[0]],1.5)
        self.assertIn("delta_log_observations",f)
    def test_stable_covariate(self):
        rows={y:[dict(key=(1,2),log_rail_km=i,**{k:0 for k in RESP})] for i,y in enumerate(YEARS)}
        with self.assertRaises(ValueError):make_change(rows,("log_rail_km",),"2025_minus_earlier_mean")
if __name__=="__main__":unittest.main()

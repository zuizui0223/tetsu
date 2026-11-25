import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_turnover_climate_block_v387 import coefficient,RAIL,FULL
class TestV387(unittest.TestCase):
    def test_block_FE_removes_between_block_artifact(self):
        rows=[]
        for b in range(12):
            for j in range(12):
                x=float(b)+j/10;y=float(j%5)
                r={k:0. for k in FULL}
                r.update(continent="Asia",block=(b,0),log_mean_rail_km=x,absdiff_log_rail_km=y,
                         beta_jtu=.01*x+.02*b,beta_jne=.005*y,beta_jac=.01*x+.02*b+.005*y)
                rows.append(r)
        out=coefficient(rows,FULL,"block")
        self.assertAlmostEqual(out["direction"][RAIL[0]]["coefficients"]["beta_jtu"],.01)
        self.assertAlmostEqual(out["direction"][RAIL[0]]["coefficients"]["beta_jne"],0.)
        self.assertGreater(coefficient(rows,FULL,"continent")["direction"][RAIL[0]]["coefficients"]["beta_jtu"],.01)
if __name__=="__main__":unittest.main()

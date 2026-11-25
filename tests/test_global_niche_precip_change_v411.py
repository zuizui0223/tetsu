import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_precip_change_v411 import add_weather,MODIFIERS,PREC,DRY_INT
class TestV411(unittest.TestCase):
    def test_sign_of_precipitation_decline(self):
        r={k:2. for k in MODIFIERS};r["key"]=(0,0)
        w={y:{(0,0):n} for y,n in (("2023",100),("2024",100),("2025",50))}
        a=add_weather([r],w,"2025_minus_earlier_mean")[0]
        self.assertGreater(a[PREC],0)
        self.assertAlmostEqual(a[DRY_INT],2*a[PREC])
        del w["2024"][(0,0)]
        self.assertEqual(add_weather([r],w,"2025_minus_earlier_mean"),[])
if __name__=="__main__":unittest.main()

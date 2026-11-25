import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_observer_month_v412 import strata_arrays,balanced_draw,YEARS
class TestV412(unittest.TestCase):
    def test_exact_common_month_and_observer_month(self):
        rows={"2023":{(1,1):[1,2,3],(2,2):[4]}, "2024":{(1,1):[2,3],(1,3):[5]}, "2025":{(1,1):[6],(1,4):[7],(3,2):[8]}}
        a,n=strata_arrays(rows,"month")
        self.assertEqual(n,{1:2})
        for y in YEARS:self.assertEqual(balanced_draw(a[y],n,5)[0],2)
        a,n=strata_arrays(rows,"observer_month")
        self.assertEqual(n,{(1,1):1})
        for y in YEARS:self.assertEqual(balanced_draw(a[y],n,5)[2],{1})
if __name__=="__main__":unittest.main()

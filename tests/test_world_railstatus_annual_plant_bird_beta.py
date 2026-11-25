import unittest
from unittest.mock import patch
import analyze_world_railstatus_annual_plant_bird_beta as t

class TestAnnual(unittest.TestCase):
    def test_common_and_covariance(self):
        def row(i):return {"a":(i,0),"b":(i,1),"continent":"X","block":0}
        by={"2023":[row(0),row(1)],"2024":[row(0),row(2)],"2025":[row(0),row(3)]}
        self.assertEqual([len(r) for r in t.matched(by).values()],[1,1,1])
        scores=[(1.,{("X",0):.1,("X",1):-.1}),(3.,{("X",0):.1,("X",1):-.1}),(4.,{("X",0):.1,("X",1):-.1})]
        with patch.object(t,"annual_fit",side_effect=[({"status":"estimated"},x) for x in scores]):
            z=t.annual_panel(by)
        self.assertEqual(z["earlier_mean"]["estimate"],2.)
        self.assertEqual(z["2025_minus_earlier"]["interval95_approx"],[2.,2.])
    def test_earlier_only(self):
        by={y:[] for y in t.YEARS}
        f=(1.,{("X",0):.1,("X",1):-.1})
        with patch.object(t,"annual_fit",side_effect=[({},f),({},f),({},None)]):
            z=t.annual_panel(by)
        self.assertEqual(z["earlier_mean"]["estimate"],1.)
        self.assertNotIn("all_year_mean",z)
    def test_duplicate(self):
        r={"a":(0,0),"b":(0,1),"continent":"X","block":0}
        with self.assertRaises(ValueError):t.matched({"2023":[r,r],"2024":[r],"2025":[r]})
    def test_disjoint_direction(self):
        def fit(rows,pred,target,terms,endpoint):
            i=1-int(target.split("_")[1])
            self.assertEqual(terms,(f"plants_{i}_jaccard",))
            self.assertIn(terms[0],pred)
            self.assertNotIn(target,pred)
            return {terms[0]:(.5,{("X",0):.1,("X",1):-.1})}
        with patch.object(t,"setup",side_effect=lambda r,*a:r),patch.object(t,"score_fit",side_effect=fit):
            z,_=t.annual_fit([])
        self.assertEqual(z["effect"]["estimate"],.5)
if __name__=="__main__":unittest.main()

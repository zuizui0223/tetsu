import unittest
from unittest.mock import patch
import analyze_world_railstatus_trophic_specificity as t

class TestSpecificity(unittest.TestCase):
    def test_direct_paired_difference(self):
        calls=[]
        def fit(rows,fields,target,term):
            self.assertNotIn(target,fields)
            guild,side,_=target.split("_")
            self.assertIn(f"insects_{1-int(side)}_species",fields)
            self.assertIn(f"plants_{1-int(side)}_species",fields)
            self.assertNotIn(f"insects_{side}_species",fields)
            calls.append((target,term))
            v=(3. if guild=="d1" else 1.) if term.startswith("insects") else (2. if guild=="d0" else 1.)
            return v,{("X",0):.1,("X",1):-.1}
        rows=[{"continent":"X","block":i} for i in range(2)]
        with patch.object(t.b,"local_support",side_effect=lambda x:x),patch.object(t.b,"with_epoch",side_effect=lambda x,y:x),patch.object(t,"linear",side_effect=fit):
            r=t.matched(rows)
        self.assertEqual(r["effects"][t.TERMS[0]]["interval95_approx"],[2.,2.])
        self.assertEqual(r["effects"][t.TERMS[1]]["estimate"],1.)
        self.assertEqual(len(calls),8)
    def test_no_support(self):
        self.assertNotEqual(t.matched([])["status"],"estimated")
if __name__=="__main__":unittest.main()

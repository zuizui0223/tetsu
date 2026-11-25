import unittest,json
from scan_world_railstatus_growthform import match,family_reference,aggregate
from analyze_world_railstatus_growthform import attach
class TestGrowthform(unittest.TestCase):
    def test_missing_and_ambiguous(self):
        rows=[{"sp":"Abc def","woodiness":"H","support":"F"},{"sp":"Abc xyz","woodiness":"NA","support":"C"},{"sp":"Abc uvw","woodiness":"W","support":"F"},{"sp":"abc UVW","woodiness":"H","support":"F"}]
        j,q=match(rows,{"abc def":[1],"abc xyz":[2],"abc uvw":[3]})
        self.assertEqual(set(j),{1,2})
        a=aggregate({1,2,3},j,{}, {})
        self.assertEqual(a["known_woodiness"],1)
        self.assertEqual(a["herbaceous"],1)
        self.assertEqual(a["matched"],2)
    def test_leave_focal_species_out(self):
        j={i:{"woodiness":"W" if i<5 else "H","support":"F"} for i in range(12)}
        f={i:1 for i in j}
        a=aggregate({0,5},j,f,family_reference(j,f))
        self.assertEqual(a["family_reference_known"],2)
        self.assertAlmostEqual(a["family_expected_woody"],9/11)
    def test_coverage_gate(self):
        r={"key":(0,0),"species":50,"observations":100,"observers":10}
        t=dict(r,known_woodiness=24,woody=6)
        self.assertEqual(attach([r],{(0,0):t},.5),[])
        self.assertEqual(attach([r],{(0,0):t},.25)[0]["observed"],.25)
if __name__=="__main__":unittest.main()

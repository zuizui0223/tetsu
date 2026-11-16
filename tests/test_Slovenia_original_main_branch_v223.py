import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from audit_Slovenia_genuine_rail_main_branch_vegetation_v223 import classify_geo_nodes,class_sections
class TestSourcePilotMath(unittest.TestCase):
    def test_single_exact_nearby_tag(self):
        pts=[(46.0,14.0,('main',)),(46.01,14.0,('branch',))]
        got=classify_geo_nodes(pts,[{'lat':46.0,'lon':14.0001}])
        self.assertEqual(got[0]['class'],'main')
    def test_ambiguous_both_tags(self):
        pts=[(46.0,14.0,('main','branch'))]
        got=classify_geo_nodes(pts,[{'lat':46.0,'lon':14.0}])
        self.assertIsNone(got[0]['class'])
    def test_missing_tag_support_is_unknown(self):
        pts=[(46.0,14.0,('main',))]
        got=classify_geo_nodes(pts,[{'lat':46.3,'lon':15.0}])
        self.assertIsNone(got[0]['class'])
    def test_conflicted_section_is_not_faked_as_branch(self):
        rows=[{'section':'1','class':'main'}]*3+[{'section':'1','class':'branch'}]*2
        x=class_sections(rows)
        self.assertIsNone(x['1']['class'])
        self.assertEqual(x['1']['classified_site_n'],5)
if __name__=="__main__":unittest.main()

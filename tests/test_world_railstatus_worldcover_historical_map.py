import unittest
from shapely.geometry import LineString
from audit_world_railstatus_worldcover_historical_map import classify,query_string,site_state,summarize,DATES,a
class TestHistory(unittest.TestCase):
    def test_tags(self):
        self.assertEqual(classify({"railway":"rail"}),(a.ACTIVE,None))
        self.assertEqual(classify({"railway":"abandoned"}),(a.INACTIVE,"abandoned"))
        self.assertEqual(classify({"railway":"disused","disused:railway":"tram"}),("other",None))
        self.assertEqual(classify({"railway":"rail","tunnel":"yes"}),("other",None))
        self.assertEqual(classify({"railway":"rail","disused:railway":"rail"}),("other",None))
    def test_historical_date(self):
        self.assertIn('[date:"2021-01-01T00:00:00Z"]',query_string([(8.,49.)],DATES[0]))
        with self.assertRaises(ValueError):query_string([(8.,49.)],"latest")
    def test_spatial_confirmation(self):
        line=LineString([(-.001,0),(.001,0)])
        self.assertEqual(site_state([(a.ACTIVE,None,line)],(0.,0.),a.ACTIVE)["status"],"confirmed")
        self.assertEqual(site_state([(a.INACTIVE,"disused",line)],(0.,0.),a.ACTIVE)["status"],"competing_or_mixed")
        self.assertEqual(site_state([], (0.,0.),a.ACTIVE)["status"],"expected_line_not_confirmed")
    def test_empty(self):
        self.assertEqual(summarize([])["both_status_confirmed_at_both_dates"],0)
if __name__=="__main__":unittest.main()

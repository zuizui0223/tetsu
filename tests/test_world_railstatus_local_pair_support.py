import unittest
from audit_world_railstatus_local_pair_support import accept,window,API
class TestSupport(unittest.TestCase):
    def test_private_missing(self):
        self.assertEqual(accept({},(0,0),{},set())[1],"obscured_or_unknown")
        self.assertEqual(accept({"obscured":False},(0,0),{},set())[1],"privacy_status_not_confirmed")
    def test_disks(self):
        r={"obscured":False,"geoprivacy":None,"taxon_geoprivacy":None,"quality_grade":"research","photos":[{}],"observed_on":"2025-05-03","positional_accuracy":20,"geojson":{"type":"Point","coordinates":[0.,0.]},"taxon":{"id":1},"user":{"id":2}}
        self.assertEqual(accept(r,(0.,0.),{1:("plants",1)},{1})[0],(1,2))
        r["geojson"]["coordinates"]=[.001,0.]
        self.assertEqual(accept(r,(0.,0.),{1:("plants",1)},{1})[1],"outside_conservative_disk")
    def test_no_geographic_replacement(self):
        self.assertIsNone(window({"eligible_surface":True,"lon":179.99,"lat":0}))
        self.assertIsNotNone(window({"eligible_surface":True,"lon":10,"lat":40}))
if __name__=="__main__":unittest.main()

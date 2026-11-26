import unittest
from scan_world_railstatus_verified_plant_size import add_size,check_old
class TestSize(unittest.TestCase):
    def test_signed_scaled(self):
        l={1:{"woody":0.},2:{"woody":1.}}
        r=[{"species":"A a","ph":"-1.5","sm":"0","ssd":"NA"},{"species":"B b","ph":"2","sm":"1","ssd":"1"}]
        a=add_size(l,r,{"a a":[1],"b b":[2]},{1,2})
        self.assertEqual(l[1]["height_nonwoody"],-1.5)
        self.assertEqual(l[1]["seed_mass_nonwoody"],0.)
        self.assertNotIn("height_nonwoody",l[2])
        self.assertEqual(a["known_nonwoody"]["height_nonwoody"],1)
if __name__=="__main__":unittest.main()

import unittest
import analyze_world_railstatus_observation_support as t
class TestSupport(unittest.TestCase):
    def row(self,k,active=0):
        r={"key":k,"continent":"X","block":(0,0),**{c+"_km":0 for c in t.CATS}}
        r["active_normal_km"]=active
        return r
    def test_nested_and_missing(self):
        k=(1,1);f=(50,10,20)
        counts={"plants":{k:f},"birds":{}}
        q=t.outcomes(self.row(k,1),counts,{"plants":{k:(f,((25,5,10),(25,5,10)))},"birds":{}})
        self.assertEqual(q["plant_record"],1.)
        self.assertEqual(q["bird_record"],0.)
        self.assertEqual(q["both_full"],0.)
        self.assertEqual(q["both_split_endpoint"],0.)
        self.assertEqual(q["rail_class"],"active_positive")
    def test_split_mismatch(self):
        k=(1,1)
        with self.assertRaises(ValueError):t.outcomes(self.row(k),{"plants":{},"birds":{}},{"plants":{k:((50,10,20),())},"birds":{}})
    def test_frame_zero_and_positive(self):
        cl={(0,0):{"continent":"X","block":(0,0),"bio1":10,"abs_lat":20},(1,0):{"continent":"X","block":(0,0),"bio1":10,"abs_lat":20}}
        r=t.make_frame(cl,[{"source_equal_area_grid_ix":1,"source_equal_area_grid_iy":0,"source_real_original_Overture_normal_rail_km_clipped_into_cell":5}])
        self.assertEqual([x["rail_km"] for x in r],[0.,5.])
    def test_primary_fields(self):
        self.assertNotIn("log_observations",t.FIELDS)
        self.assertNotIn("log_observers",t.FIELDS)
        self.assertIn(t.b.ACTIVE,t.FIELDS)
if __name__=="__main__":unittest.main()

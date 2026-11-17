import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import unpack_jsonstat,summarized
class TestEurostatSegmentShape(unittest.TestCase):
    def sample(self):
        return {
            "version":"2.0","class":"dataset","label":"Synthetic, not real observed traffic",
            "id":["freq","net_seg20","train","unit","time"],"size":[1,2,2,1,1],
            "dimension":{
              "freq":{"category":{"index":{"A":0},"label":{"A":"Annual"}}},
              "net_seg20":{"category":{"index":{"SI01":0,"SI02":1},"label":{"SI01":"A - B","SI02":"C - D"}}},
              "train":{"category":{"index":{"FRT":0,"PAS":1},"label":{"FRT":"Freight train","PAS":"Passenger train"}}},
              "unit":{"category":{"index":{"NR":0},"label":{"NR":"Number"}}},
              "time":{"category":{"index":{"2020":0},"label":{"2020":"2020"}}}},
            "value":{"0":100,"1":20,"2":0,"3":10}}
    def test_indices_last_axis_fastest_and_explicit_zero(self):
        rows,qa=unpack_jsonstat(self.sample(),"2020")
        self.assertEqual(len(rows),4)
        self.assertEqual([x["value"] for x in rows],[100.,20.,0.,10.])
        self.assertEqual(rows[2]["net_seg20"]["code"],"SI02")
        self.assertEqual(rows[1]["train"]["code"],"PAS")
    def test_missing_cells_not_assumed_absent_rail_lines(self):
        fake=self.sample();del fake["value"]["2"]
        rows,qa=unpack_jsonstat(fake,"2020")
        self.assertEqual(qa["missing_2020_vs_2015_source_cells"],1)
        self.assertNotEqual(len(rows),4)
    def test_year_rejection(self):
        with self.assertRaises(ValueError):unpack_jsonstat(self.sample(),"2023")
    def test_source_stats_classification(self):
        rows,_=unpack_jsonstat(self.sample(),"2020")
        out=summarized(rows)
        self.assertEqual(out["unique_nonnull_source_segments"],2)
        self.assertEqual(out["source_train_types_counts"]["Freight train"],2)
    def test_malformed_source_rejected(self):
        d=self.sample();d["dimension"]["net_seg20"]["category"]["index"].pop("SI02")
        with self.assertRaises(ValueError):unpack_jsonstat(d,"2020")
if __name__=="__main__":unittest.main()

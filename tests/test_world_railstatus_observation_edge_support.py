import unittest
import analyze_world_railstatus_observation_edge_support as t

class TestEdges(unittest.TestCase):
    def test_independent_pairs(self):
        def row(k,active=1,obs=1):
            r={f:0. for f in t.CELL};r.update({"key":k,"continent":"X","block":(0,0),"both_record":obs,"both_effort":obs,"both_full":obs,"plant_effort":obs,"bird_effort":obs})
            r.update({c+"_km":0 for c in t.CATS});r["active_normal_km"]=active;return r
        a,b,c=row((0,0)),row((1,0)),row((2,0),0,0)
        ed={x:{} for x in t.CATS};ed["active_normal"][((0,0),(1,0))]=1
        sp={g:{((0,0),(1,0))} for g in t.GROUPS}
        p=t.make_pairs([a,b,c],ed,sp)
        self.assertEqual(len(p),2)
        self.assertEqual(p[0]["both_split_pairs"],1)
        self.assertEqual(p[0][t.ACTIVE],1)
        self.assertEqual(p[1]["both_effort"],0)
        self.assertFalse(p[1]["both_active"])
        ed["active_normal"][((1,0),(2,0))]=1
        with self.assertRaises(ValueError):t.make_pairs([a,b,c],ed,sp)
    def test_fields(self):
        self.assertEqual(len(t.FIELDS),len(set(t.FIELDS)))
        self.assertFalse(any("log_observations" in x or "log_observers" in x for x in t.FIELDS))
        self.assertIn(t.ACTIVE,t.ENDPOINT)
        self.assertIn(t.INACTIVE,t.ENDPOINT)
    def test_bad_pair(self):
        with self.assertRaises(ValueError):t.pairkey({"a":[0,0],"b":[2,0]})
if __name__=="__main__":unittest.main()

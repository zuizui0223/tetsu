import unittest
import analyze_world_railstatus_observation_stages as t
class TestStages(unittest.TestCase):
    def row(self):
        return {"key":(1,1),**{c+"_km":0 for c in t.CATS}}
    def test_effort_distinct_from_species(self):
        k=(1,1)
        for n,expected in [((20,5,15),(1,0,0,0)),((50,5,15),(1,1,0,0)),((50,10,20),(1,1,1,0))]:
            r=t.outcomes(self.row(),{"plants":{k:n},"birds":{k:n}},{"plants":{},"birds":{}})
            self.assertEqual(tuple(r[x] for x in t.STAGES),expected)
            self.assertEqual(sum(r[x] for x in t.LOSSES),1)
    def test_split_endpoint(self):
        k=(1,1);n=(50,10,20);split={k:(n,((25,5,10),(25,5,10)))}
        r=t.outcomes(self.row(),{"plants":{k:n},"birds":{k:n}},{"plants":split,"birds":split})
        self.assertEqual(tuple(r[x] for x in t.STAGES),(1,1,1,1))
        self.assertEqual(sum(r[x] for x in t.LOSSES),0)
    def test_empty(self):
        self.assertEqual(t.rates([])["conditional_pass"],{k:None for k in t.STAGES[1:]})
    def test_no_effort_leakage(self):
        self.assertNotIn("log_observations",t.FIELDS)
        self.assertNotIn("log_observers",t.FIELDS)
if __name__=="__main__":unittest.main()

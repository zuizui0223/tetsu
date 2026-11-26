import unittest
import numpy as np
import scan_world_railstatus_fixed_observer_richness as s
import analyze_world_railstatus_fixed_observer_richness as a
class TestFixed(unittest.TestCase):
    def test_all_distinct(self):
        r=s.summarize([1,2,3,4],[11,12,13,14],(4,4,4),5,reps=32,k=3)
        self.assertEqual(r["observer_balanced"]["mean_species"],3)
        self.assertEqual(r["record_balanced"]["mean_species"],3)
    def test_repeated_observer(self):
        o=[1]*30+[2,3,4];v=[11]*30+[12,13,14]
        r=s.summarize(o,v,(33,4,4),5,reps=256,k=3)
        self.assertEqual(r["observer_balanced"]["mean_species"],3)
        self.assertLess(r["record_balanced"]["mean_species"],2)
        self.assertEqual(r,s.summarize(o[::-1],v[::-1],(33,4,4),5,reps=256,k=3))
    def test_missing_observer_and_exact_counts(self):
        r=s.summarize([-1,1,2,3],[21,11,12,13],(4,3,4),8,reps=16,k=3)
        self.assertEqual(r["known_records"],3)
        self.assertEqual(r["observer_balanced"]["mean_species"],3)
        with self.assertRaises(ValueError):s.summarize([1,2],[11,12],(2,2,1),2,reps=16,k=2)
    def test_gate_no_species_threshold(self):
        z={"source_equal_area_grid_ix":0,"source_equal_area_grid_iy":0,"observations":50,"observers":20,"species":1}
        p={"grids":{str(w):{g:[z] for g in s.GROUPS} for w in s.SOURCE_REGION_M}}
        self.assertEqual(s.candidates(p)[100000]["plants"][(0,0)],(50,20,1))
    def test_base_no_effort(self):
        self.assertNotIn("log_observations",a.BASE)
        self.assertNotIn("log_observers",a.BASE)
        self.assertEqual(len(a.BASE),len(set(a.BASE)))
if __name__=="__main__":unittest.main()

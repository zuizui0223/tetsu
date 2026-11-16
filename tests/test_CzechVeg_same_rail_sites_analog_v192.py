import importlib.util,unittest,random
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/analyze_CzechVeg_same_rail_sites_analog_vs_background_v192.py"
sp=importlib.util.spec_from_file_location("v192",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
class TestMatchedRailSites(unittest.TestCase):
 def test_same_sites_and_reference_types(self):
  n=6
  plants={}
  frame={}
  for i in range(n):
   k=str(i)
   rail="rail"+k
   plants[rail]={"A a","B b","D d"}
   cs={}
   for kind in m.CATEGORIES:
    ids=[]
    for j in range(5):
     name=k+"_"+kind+str(j)
     plants[name]=({"A a","B b"} if kind=="analog" else
                   {"A a","D d"} if kind=="nonanalog" else {"A a","Y y"})
     ids.append(name)
    cs[kind]=ids
   frame[k]=[(rail,cs)]
  v=m.run_pair_rail_reference(frame,plants,draws=3)
  self.assertTrue(v["identified"])
  self.assertEqual(v["regions_count"],6)
  self.assertTrue(v["same_original_rail_plot_drawn_for_ALL_reference_choices"])
  self.assertIn("species_rail_minus_analog",v["resampling_results"])
 def test_common_rail_ids_intersection_all_refs(self):
  m0={}
  for i in range(3):
   sid="r"+str(i)
   m0[sid]={"region":"R"+str(i)}
  arms={
   "analog":{"r0":["c0"]*5,"r1":["c1"]*5,"r2":["c2"]*5},
   "nonanalog":{"r0":["d0"]*5,"r1":["d1"]*5,"r2":["d2"]*5},
   "all":{"r0":["c0"]*5,"r1":["c1"]*5}
  }
  p={k:{"A a","B b"} for k in ["r0","r1","r2","c0","c1","c2","d0","d1","d2"]}
  f,qa=m.common_frame(arms,m0,p)
  self.assertEqual(qa["source_rail_sites_with_all_three_reference_kinds"],2)
  self.assertEqual(qa["n_original_phyto_regions_with_all_reference_kinds"],2)
if __name__=="__main__":unittest.main()

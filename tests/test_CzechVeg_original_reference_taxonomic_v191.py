import importlib.util,unittest,random
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/analyze_CzechVeg_original_independent_reference_taxonomic_v191.py"
sp=importlib.util.spec_from_file_location("cv191",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
class TestCrossRegionReference(unittest.TestCase):
 def test_plant_similarity(self):
  self.assertEqual(m.jaccard({"A","B"},{"A","C"}),1/3)
  self.assertEqual(m.jaccard({"A"},{"A"}),1)
  self.assertEqual(m.jaccard({"A"},{"C"}),0)
  self.assertEqual(m.genus_set({"Acer rubrum","Acer negundo","Poa annua"}),{"Acer","Poa"})
 def test_duplicate_unobserved_background_not_artificial_natural(self):
  sites={"r1":{"study":"A","region":"region","decade":2010,"area":"1_to_4","label":"rail","strict":None,"source_class":"TH"},
         "b1":{"study":"B","region":"region","decade":2010,"area":"1_to_4","label":None,"strict":None,"source_class":"TF"},
         "b2":{"study":"A","region":"region","decade":2010,"area":"1_to_4","label":None,"strict":None,"source_class":"TH"},
         "b3":{"study":"B","region":"region","decade":2010,"area":"1_to_4","label":None,"strict":None,"source_class":"TH"}}
  blocks={("A","region",2010,"1_to_4"):{"rail":["r1"],"road":["n/a"]}}
  arms,ids,_=m.pick_source_candidates(sites,blocks)
  self.assertEqual(len(arms["analog"]),0)
  self.assertTrue("r1" not in ids)
 def test_heldout_regions_and_richness(self):
  plant={"rA":{"A aa","B bb"},"rB":{"A aa","C cc"},
         "rC":{"A aa","D dd"},"cA":{"B bb","X xx"},
         "cB":{"C cc","Y yy"},"cC":{"D dd","Z zz"}}
  src={"A":[("rA",["cA"]*5)],"B":[("rB",["cB"]*5)],"C":[("rC",["cC"]*5)]}
  result=m.run_reference_panel(src,plant,"synthetic",reps=5)
  self.assertTrue(result["identified"])
  self.assertEqual(result["n_nonindependent_cross_region_pairs_each_draw"],3)
  self.assertAlmostEqual(result["response_statistics"]["species_delta"]["mean"],1/3,places=6)
if __name__=="__main__":unittest.main()

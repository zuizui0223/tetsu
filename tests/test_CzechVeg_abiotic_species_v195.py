import importlib.util,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/analyze_CzechVeg_abiotic_matched_reference_taxonomic_v195.py"
sp=importlib.util.spec_from_file_location("v195",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

class TestAbioticMatchedNoPlantResponseLeakage(unittest.TestCase):
 def synthetic(self):
  meta={};measured={};blocks={};plant={}
  for i in range(4):
   region="r"+str(i)
   r="rail_"+region
   k=("P"+str(i),region,2010,"gt4_to_16")
   blocks[k]={"rail":[r],"road":["road_ignored"]}
   meta[r]={"study":"P"+str(i),"region":region,"decade":2010,
            "area":"gt4_to_16","label":"rail","strict":"rail","source_class":"TH"}
   measured[r]={"elev_m":200.,"slope_deg":5.,"rock_cover_pct":0.}
   plant[r]={"A aa","B bb","C cc","D dd","E ee","F ff"}
   for j in range(6):
    pid=f"ctrl_{region}_{j}"
    meta[pid]={"study":"independent"+str(j),"region":region,
               "decade":2010,"area":"gt4_to_16",
               "label":None,"strict":None,"source_class":"TH"}
    measured[pid]={"elev_m":202.+j,"slope_deg":6.,"rock_cover_pct":0.}
    plant[pid]={"A aa","B bb","G gg","H hh","I ii","J jj"}
  return meta,measured,blocks,plant
 def test_original_abiotic_source_must_not_use_syntaxon_to_select(self):
  meta,env,blocks,plant=self.synthetic()
  a,needed,qa=m.source_arms(meta,blocks,env)
  self.assertEqual(qa["source_publisher_rail_plots_with_min5_elev100"],4)
  self.assertEqual(qa["source_publisher_rail_plots_with_min5_elev100_slope10"],4)
  self.assertEqual(len(a["all"]),4)
  self.assertTrue(qa["only_original_field_elevation_slope_selected"])
  for k in a["elev100"]: self.assertEqual(len(a["elev100"][k]),6)
 def test_matching_common_rail_for_all_two_types_and_richness(self):
  meta,env,blocks,plant=self.synthetic()
  a,needed,qa=m.source_arms(meta,blocks,env)
  for scope in m.FRAME_TYPES:
   frame,receipt=m.frame_for_type(a,meta,plant,scope,richness=True)
   self.assertEqual(receipt["selected_original_rail_plots"],4)
   self.assertEqual(receipt["distinct_source_phyto_regions"],4)
   self.assertTrue(receipt["same_source_rail_sites_for_abiotic_and_all_background"])
   outcome=m.evaluate(frame,plant,seed=6,draws=7)
   self.assertTrue(outcome["identified"])
   self.assertEqual(outcome["source_region_pairs_nonindependent_per_draw"],6)
   self.assertEqual(len(outcome["summary"]),12)
   self.assertEqual(outcome["summary"]["species_rail_minus_abiotic"]["mean"],
                    outcome["summary"]["species_rail_minus_all"]["mean"])
 def test_fail_closed_when_no_phyto_regions(self):
  meta,env,blocks,plant=self.synthetic()
  a,needed,qa=m.source_arms(meta,blocks,env)
  del a["elev100"]["rail_r2"]
  del a["elev100"]["rail_r3"]
  frame,receipt=m.frame_for_type(a,meta,plant,"elev100")
  outcome=m.evaluate(frame,plant,draws=4)
  self.assertFalse(outcome["identified"])
  self.assertEqual(outcome["original_source_phyto_region_count"],2)
if __name__=="__main__":unittest.main()

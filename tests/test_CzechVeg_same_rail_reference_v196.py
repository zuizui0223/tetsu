import importlib.util, unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/analyze_CzechVeg_same_rail_abiotic_vs_syntaxon_v196.py"
sp=importlib.util.spec_from_file_location("source196",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

class TestStrictSameRailOriginalSourceFrame(unittest.TestCase):
 def source(self):
  meta={};env={};blocks={};taxa={}
  for i in range(4):
   reg="R"+str(i)
   rail="rail_"+reg
   paper="rail_paper_"+reg
   k=(paper,reg,2010,"gt4_to_16")
   blocks[k]={"rail":[rail],"road":["road_unused"]}
   meta[rail]={"study":paper,"region":reg,"decade":2010,
                "area":"gt4_to_16","label":"rail","strict":"rail","source_class":None}
   env[rail]={"elev_m":300.0,"slope_deg":5.0,"rock_cover_pct":0.0}
   taxa[rail]={"A aa","B bb","C cc","D dd","E ee","F ff"}
   for kind in ("analog","nonanalog"):
    for j in range(6):
     key=f"{kind}_{reg}_{j}"
     meta[key]={"study":"independent_"+reg+"_"+str(j),
                "region":reg,"decade":2010,"area":"gt4_to_16",
                "label":None,"strict":None,
                "source_class":"TH" if kind=="analog" else "XA"}
     env[key]={"elev_m":300+j,"slope_deg":5,
               "rock_cover_pct":0}
     taxa[key]={"A aa","B bb","C cc","D dd","G gg","H hh"}
  return meta,env,blocks,taxa
 def test_source_only_preselection(self):
  meta,env,blocks,plant=self.source()
  arm,read,qa=m.combined_arms(meta,blocks,env)
  self.assertEqual(qa["source_eligible_rail_sites_pre_species"],4)
  self.assertEqual(len(read),52)
  self.assertEqual(set(arm["rail_R1"]),set(m.CONTROL_TYPES))
  self.assertTrue(qa["all_reference_types_formed_before_species_contents"])
 def test_joint_source_site_frame_and_limited_independent_n(self):
  meta,env,blocks,plants=self.source()
  arms,_,_=m.combined_arms(meta,blocks,env)
  frame,qa=m.common_frame(arms,plants,meta,True)
  self.assertEqual(qa["n_same_original_rail_text_sites"],4)
  self.assertEqual(qa["n_covered_source_phyto_regions"],4)
  self.assertTrue(qa["same_rail_site_and_geographic_study_design_for_all_reference_choices"])
  r=m.analyze(frame,plants,draws=5)
  self.assertTrue(r["identified"])
  self.assertEqual(r["nonindependent_source_region_pairs"],6)
  self.assertIn("species_analog_minus_abiotic",r["summary"])
  self.assertTrue(r["no_field_study_p_values_or_CI"])
 def test_missing_reference_blocks_not_reported_as_biological_null(self):
  meta,env,blocks,plants=self.source()
  arms,_,_=m.combined_arms(meta,blocks,env)
  del arms["rail_R0"]
  del arms["rail_R1"]
  frame,qa=m.common_frame(arms,plants,meta)
  out=m.analyze(frame,plants,draws=4)
  self.assertFalse(out["identified"])
  self.assertEqual(out["n_geographic_regions"],2)
if __name__=="__main__":unittest.main()

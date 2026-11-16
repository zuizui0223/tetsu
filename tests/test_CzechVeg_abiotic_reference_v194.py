import importlib.util,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/audit_CzechVeg_abiotic_regional_reference_opportunity_v194.py"
spec=importlib.util.spec_from_file_location("env194",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestAbioticNonSyntaxonControl(unittest.TestCase):
 def test_numeric_missing_and_implausible_original_values(self):
  self.assertIsNone(m.number("-1",0,2500))
  self.assertIsNone(m.number("NA",0,2500))
  self.assertEqual(m.number("500",0,2500),500)
  self.assertEqual(m.number("13,5",0,90),13.5)
  self.assertIsNone(m.number("2345",0,90))
 def test_environmental_match_without_syntaxon(self):
  rail={"elev_m":200,"slope_deg":8,"rock_cover_pct":25}
  ctrl={"elev_m":240,"slope_deg":14,"rock_cover_pct":13}
  a=m.matching(rail,ctrl)
  self.assertEqual(set(a),set(m.MATCHES))
  a=m.matching(rail,{"elev_m":370,"slope_deg":14,"rock_cover_pct":13})
  self.assertEqual(a,["elev200"])
  self.assertEqual(m.matching(rail,{"elev_m":None,"slope_deg":0,"rock_cover_pct":30}),[])
 def test_independent_publication_reference_not_source_rail(self):
  def md(study,label):
   return {"study":study,"region":"R","decade":2010,"area":"gt4_to_16",
           "label":label,"strict":label,"source_class":"TH"}
  meta={"r":md("RAIL","rail"),"a":md("RAIL",None),"b":md("B",None),"c":md("C",None)}
  blocks={("RAIL","R",2010,"gt4_to_16"):{"rail":["r"],"road":["any"]}}
  options=m.source_reference_options(meta,blocks)
  self.assertEqual(set(options["r"]),{"b","c"})
  a={key:{"elev_m":150,"slope_deg":2,"rock_cover_pct":20,
          "position_uncertainty_m":10,
          "geology_recorded":True,"soil_recorded":False} for key in meta}
  b=m.audit(meta,blocks,a)
  self.assertEqual(b["source_reference_type_support"]["elev100"]["rail_plots_with_at_least_one_reference"],1)
  self.assertEqual(b["source_reference_type_support"]["elev100"]["rail_plots_with_min5_references"],0)
  self.assertFalse(b["source_reference_type_support"]["elev100"]["vegetation_source_THTF_label_used_for_filter"])
if __name__=="__main__":unittest.main()

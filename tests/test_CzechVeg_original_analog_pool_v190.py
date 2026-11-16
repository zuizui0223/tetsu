import importlib.util
from pathlib import Path
import unittest
p=Path(__file__).resolve().parents[1]/"scripts/audit_CzechVeg_original_nonrail_analog_reference_pool_v190.py"
sp=importlib.util.spec_from_file_location("pool190",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

class TestOriginalPublisherAnalogSourceFrame(unittest.TestCase):
 def test_exact_published_pladias_natural_analog_classes(self):
  self.assertEqual(set(m.CLASSES),{"TH","TF","TG"})
  self.assertEqual(m.class_two_char("THA01"),"TH")
  self.assertEqual(m.class_two_char("TF"),"TF")
  self.assertEqual(m.class_two_char("TGA01"),"TG")
  self.assertIsNone(m.class_two_char("SOMETHING UNVERIFIED"))
  self.assertIsNone(m.class_two_char("NA"))
 def test_independent_pub_sample_design(self):
  def meta(group,study="a",cl="TH",dec=2010):
   return {"study":study,"region":"B7","decade":dec,"year":2015,
           "area":"gt16_to_64","label":group,"strict":None,
           "source_class":cl}
  sites={"r1":meta("rail"),"r2":meta("rail"),
         "d1":meta("road"),"d2":meta("road"),
         "th":meta(None,"other","TH"),
         "tg":meta(None,"other","TG"),
         "tf":meta(None,"a","TF"),
         "x":meta(None,"other","XA")}
  matched,a=m.matched_rail_source(sites)
  self.assertEqual(a["n_rail_plots"],2)
  opp=m.opportunity(sites,matched)
  self.assertEqual(opp["support"]["rail_plots_analog_pool_min2"],2)
  self.assertEqual(opp["support"]["rail_plots_analog_pool_min5"],0)
  self.assertTrue(opp["nontransport_label_is_missing_NOT_confirmed_natural"])
if __name__=="__main__":unittest.main()

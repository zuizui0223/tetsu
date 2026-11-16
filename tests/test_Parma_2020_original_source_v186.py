import unittest
import importlib.util
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/audit_Parma_2020_28rail_station_flora_original_v186.py"
sp=importlib.util.spec_from_file_location("v186",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
class TestOriginalSourceAudit(unittest.TestCase):
 def test_publisher_28_codes_and_CIT_never_zero(self):
  self.assertEqual(len(m.STATION_CODES),28)
  self.assertNotIn(m.EXCLUDED_UNSURVEYED_STATION,m.STATION_CODES)
  self.assertEqual(len(set(m.STATION_CODES)),28)
 def test_frequency_is_publisher_station_number_not_imputed(self):
  self.assertEqual(m.FREQ.search("Zygophyllum fabago L., BER, PAR [2]").group(1),"2")
  self.assertFalse(m.FREQ.search("No frequency [not a number]"))
 def test_noneligible_natural_language_is_not_taxon(self):
  self.assertIsNone(m.SPECIES_START.match("Annali Museo Civico 36/20"))
  self.assertIsNotNone(m.SPECIES_START.match("Abies alba Mill."))
 def test_rail_not_equivalent_pure_ballast(self):
  self.assertTrue("Adorni" in m.ORIGINAL_URL)
if __name__=="__main__":unittest.main()

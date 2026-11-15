import importlib.util
from pathlib import Path
import unittest
p=Path(__file__).resolve().parents[1]/"scripts/audit_Parma_2020_28station_table3_table4_crosscheck_v189.py"
spec=importlib.util.spec_from_file_location("v189",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class TestSourcePublisherValidation(unittest.TestCase):
 def test_28_original_publisher_station_counts(self):
  self.assertEqual(len(m.TABLE3),28)
  self.assertEqual(sum(m.TABLE3.values()),2864)
  self.assertEqual(m.TABLE3["CHI"],32)
  self.assertEqual(m.TABLE3["SAL"],157)
 def test_27_author_table4_count_benchmarks(self):
  self.assertEqual(len(m.TABLE4),27)
  self.assertEqual(m.TABLE4["Galium aparine"],27)
  self.assertEqual(m.TABLE4["Taraxacum sect. Taraxacum"],20)
 def test_source_taxon_key_preserves_section(self):
  self.assertEqual(m.taxon_key("Taraxacum F.H.Wigg. sect. Taraxacum"),"Taraxacum sect. Taraxacum")
  self.assertEqual(m.taxon_key("Acer negundo L."),"Acer negundo")
if __name__=="__main__":unittest.main()

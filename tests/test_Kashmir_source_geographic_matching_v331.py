"""v331 source-only metadata helper unit tests."""
import unittest,sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Kashmir_source_geographic_matching_v331 as m
class GeographicMetadata(unittest.TestCase):
 def test_numeric(self):
  self.assertEqual(m.num(3.5),3.5)
  self.assertIsNone(m.num("NA"))
 def test_nonfinite_rejected(self):
  self.assertIsNone(m.num("inf"))
  self.assertIsNone(m.num(float("nan")))
 def test_original_hash_full_length(self):
  self.assertEqual(len(m.DRYAD_SHA256),64)
if __name__=="__main__":unittest.main()

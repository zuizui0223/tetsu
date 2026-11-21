"""v310 synthetic DMS geography independent integer-number and digit validation."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
 from analyze_Slovakia2016_19station_original_physical_distance_v310 import GPS_PATTERN,geodist
except ModuleNotFoundError as e:
 if e.name not in ("fitz","numpy"):raise
 GPS_PATTERN=None
@unittest.skipUnless(GPS_PATTERN is not None,"Original PDF DMS parser has optional CI dependencies")
class PublisherDMS(unittest.TestCase):
 def test_publisher_original_full_DMS_form_not_latlon_from_train_traffic(self):
  example="  1. Poprad-Tatry, 49°3´36˝, 20°17´44˝, 673 m n. m., 6987b, 26b"
  m=GPS_PATTERN.search(example)
  self.assertIsNotNone(m)
  self.assertEqual(int(m.group('id')),1)
  self.assertEqual(int(m.group('alt')),673)
 def test_true_geodesic_symmetry_and_nonzero(self):
  a=(49.05,20.25,600,"x")
  b=(49.15,20.54,500,"y")
  self.assertAlmostEqual(geodist(a,b),geodist(b,a))
  self.assertGreater(geodist(a,b),10)
  self.assertEqual(geodist(a,a),0)
if __name__=="__main__":unittest.main()

import importlib.util
from pathlib import Path
import unittest
p=Path(__file__).resolve().parents[1]/"scripts/extract_Parma_2020_28station_species_matrix_v187.py"
spec=importlib.util.spec_from_file_location("v187",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestParmaSourceParser(unittest.TestCase):
    def test_two_synthetic_records_with_wrap_and_count(self):
        raw1="Abies alba Mill. (Pinaceae), Xenof. loc., Orof. S-Europ., P scap, SAL"
        row,err=m.parse_entry(raw1,1)
        self.assertIsNone(err)
        self.assertEqual(row["taxon_original"],"Abies alba Mill.")
        self.assertEqual(row["lifeform_coarse"],"P")
        self.assertEqual(row["station_codes"],["SAL"])
        self.assertEqual(row["native_status_original"],"Xenof. loc.")
    def test_mixed_station_codes_not_zero(self):
        row,err=m.parse_entry("Acer campestre L. (Sapindaceae), Autoct., Europ.-Caucas., P scap, BER, PAR, SAL",3)
        self.assertIsNone(err)
        self.assertEqual(set(row["station_codes"]),{"BER","PAR","SAL"})
        self.assertNotIn("CIT",row["station_codes"])
        self.assertEqual(row["lifeform_coarse"],"P")
    def test_source_inconsistency_does_not_impute_station_codes(self):
        row,err=m.parse_entry("Acer negundo L. (Sapindaceae), Neof. nat., N-Americ., P scap, BUS, PAR",16)
        self.assertIsNone(err)
        self.assertEqual(row["printed_station_total"],16)
        self.assertEqual(row["parsed_station_total"],2)
        self.assertEqual(row["printed_minus_codes"],14)
    def test_source_family_required(self):
        row,err=m.parse_entry("not an identified species, SAL",1)
        self.assertIsNone(row)
        self.assertEqual(err["cause"],"no_publisher_family_pattern")
    def test_lossless_hyphen_wrap(self):
        self.assertIn("Eurimedit.",m.concatenate(["Euri-","medit., T scap"],[10,10]))
    def test_declared_source_scope(self):
        self.assertEqual(m.EXPECTED_TAXA,438)
        self.assertEqual(len(m.STATION_CODES),28)

if __name__=="__main__":unittest.main()

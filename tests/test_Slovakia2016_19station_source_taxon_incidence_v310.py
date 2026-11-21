"""Synthetic incidence conservation tests; only author PDF source run validates ecology."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
    from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import (
        clean_original_page,expand_station_codes,source_list_clauses,
        parse_clauses,EXPECTED_SOURCE_TAXA,
    )
except ModuleNotFoundError as exc:
    if exc.name not in ("fitz","numpy"):raise
    expand_station_codes=None

@unittest.skipUnless(expand_station_codes is not None,"Source CI has optional original PDF numerical libs")
class OriginalPaperIncidenceMath(unittest.TestCase):
    def test_author_1_through_19_expands_exactly(self):
        self.assertEqual(expand_station_codes("1–19"),set(range(1,20)))
    def test_mixed_ranges_and_source_singletons(self):
        self.assertEqual(expand_station_codes("1, 3–6, 12, 18–19"),{1,3,4,5,6,12,18,19})
    def test_duplicate_or_out_of_range_or_noninteger_not_accepted(self):
        for bad in ["2,2","1–4,4","0","20","1-21","1.5","1–x","","1,"]:
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):expand_station_codes(bad)
    def test_original_hyphenated_genus_epithet_reconstruction(self):
        self.assertIn("hordeaceus",clean_original_page("B. hor-\ndeaceus: 2, 4"))
        self.assertNotIn("hor-",clean_original_page("B. hor-\ndeaceus: 2, 4"))
    def test_318_original_species_strict_gate(self):
        self.assertEqual(EXPECTED_SOURCE_TAXA,318)
        with self.assertRaises(ValueError):
            parse_clauses(["Acer platanoides: 1, 2"]*317)
if __name__=="__main__":unittest.main()

import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from audit_SI_2023_PDF_vector_station_table_traffic_v234 import (
 classify_cell_row,merge_old_new
)
class TestOriginalSourcePDFFirstFourColumns(unittest.TestCase):
    def test_exact_2023_station_first_four_real_table_columns(self):
        x,line,why=classify_cell_row(
          ['Ljubljana\nZalog','557.052','561.047','3996','3','518/482'], '10')
        self.assertEqual(why,'strict_official_station_table_row')
        self.assertEqual(line,'10')
        self.assertEqual(x['name_normalized'],'ljubljana zalog')
        self.assertAlmostEqual(x['station_midpoint_km'],559.0495)
    def test_original_table_header_sets_real_national_route_code(self):
        x,line,why=classify_cell_row(['LINE 20 LJUBLJANA – JESENICE – D.M.',None,None,None], '10')
        self.assertEqual(line,'20')
        self.assertIsNone(x)
    def test_platform_length_misread_as_station_km_rejected(self):
        x,line,why=classify_cell_row(['Ljubljana Zalog','557.052','561.047','518','4','...'],'10')
        self.assertIsNone(x)
        self.assertEqual(why,'length_km_not_source_consistent')
    def test_ambiguous_source_name_refuses_silent_join(self):
        base=[{'national_line':'10','name':'Another Name','name_normalized':'another name',
               'station_start_km':557.052,'station_end_km':561.047}]
        truth=[{'national_line':'10','name':'Ljubljana Zalog','name_normalized':'ljubljana zalog',
               'station_start_km':557.052,'station_end_km':561.047}]
        got,qa=merge_old_new(base,truth)
        self.assertEqual(len(got),0)
        self.assertEqual(qa['conflicted_source_table_text_label_excluded'],1)
    def test_same_station_with_wrap_is_repaired_not_duplicated(self):
        base=[{'national_line':'10','name':'Zalog','name_normalized':'zalog',
               'station_start_km':557.052,'station_end_km':561.047}]
        truth=[{'national_line':'10','name':'Ljubljana Zalog','name_normalized':'ljubljana zalog',
               'station_start_km':557.052,'station_end_km':561.047}]
        got,qa=merge_old_new(base,truth)
        self.assertEqual(len(got),1)
        self.assertEqual(got[0]['name_normalized'],'ljubljana zalog')
        self.assertEqual(qa['verified_multiline_wrapped_station_name_repaired'],1)
if __name__=='__main__':unittest.main()

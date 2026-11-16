import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
    normalize,rows_from_text,segment_route_station_join,assign_segment
)
class TestAuthoritativeStationSourceJoin(unittest.TestCase):
    def test_parse_publisher_station_km_not_raw_track_data(self):
        x=['LINE 10 D.M. – DOBOVA – LJUBLJANA',
           'Dobova 452.685 454.306 1620',
           'Brežice 458.707 459.851 1144',
           '1 908 P,1 394,37 6,3',
           'Sevnica 485.220 486.331 1112',
           'Zidani Most 500.915 502.863 1948',
           'Ljubljana','Zalog 557.052 561.047 3996']
        ans,cur=rows_from_text(x)
        self.assertEqual(cur,'10')
        self.assertEqual(len(ans),5)
        self.assertEqual(ans[-1]['name_normalized'],'ljubljana zalog')
        self.assertAlmostEqual(ans[0]['station_start_km'],452.685)
    def test_name_accent_normalization_is_not_guessed_route_number(self):
        self.assertEqual(normalize('Ljubljana Šiška'),normalize('Ljubljana Siska'))
        self.assertNotEqual(normalize('Divača'),normalize('Ljubljana'))
    def test_exact_unique_line_station_two_endpoints(self):
        a,_=rows_from_text(['LINE 10 D.M. – DOBOVA – LJUBLJANA','Sevnica 485.220 486.331 1112','Zidani Most 500.915 502.863 1948'])
        segs=[{'source_network_segment_code':'SIS10003_TEN','source_original_label':'Sevnica - Zidani Most',
              'start_normalized':'sevnica','end_normalized':'zidani most',
              'actual_2020_freight_trains_per_year':4427,'actual_2020_passenger_trains_per_year':6540}]
        q,qa=segment_route_station_join(a,segs)
        self.assertEqual(qa['exact_same_line_unique_two_named_endpoints'],1)
        self.assertEqual(q[0]['national_line'],'10')
        self.assertEqual(q[0]['actual_2020_freight_trains_per_year'],4427)
        sites=[{'site':'siteA','section':'10','mid_source_route_chainage_km':490,'known':{'a','b'}},
               {'site':'siteB','section':'10','mid_source_route_chainage_km':470,'known':{'a','b'}}]
        assigned,coverage=assign_segment(sites,q)
        self.assertEqual(len(assigned),1)
        self.assertEqual(assigned[0]['source_site'],'siteA')
        self.assertEqual(coverage['unmapped_genuine_site_reasons']['no_segment_contains_original_site_midpoint'],1)
    def test_crosswalk_refuses_ambiguous_multi_line_same_station_names(self):
        stations=[]
        for line in ('10','20'):
            stations.extend([{'national_line':line,'name_normalized':'sevnica',
                              'station_midpoint_km':480.},
                             {'national_line':line,'name_normalized':'zidani most',
                              'station_midpoint_km':500.}])
        seg=[{'source_original_label':'Sevnica - Zidani Most','start_normalized':'sevnica',
              'end_normalized':'zidani most','source_network_segment_code':'X'}]
        res,qa=segment_route_station_join(stations,seg)
        self.assertEqual(len(res),0)
        self.assertEqual(qa['unresolved_reasons']['ambiguous_multiple_national_lines'],1)
if __name__=='__main__':unittest.main()

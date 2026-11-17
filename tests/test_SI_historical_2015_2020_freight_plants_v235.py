import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from validate_SI_historical_2015_2020_train_traffic_vs_flora_v235 import (
 matched_traffic_for_2020_site_segments
)
class TestSameEurostatSourceSegmentAcrossYears(unittest.TestCase):
    def make(self):
        names={"SIS10002_TEN":"Dobova - Sevnica",
               "SIS20002_TEN":"Sevnica - Dobova"}
        source2020=[{'source_network_segment_code':code,'source_original_label':label} for code,label in names.items()]
        old=[]
        values={
            "SIS10002_TEN":("TRN_GD",2500,"TRN_PAS",7000),
            "SIS20002_TEN":("TRN_GD",2450,"TRN_PAS",6800)}
        for code,label in names.items():
            x,y,z,w=values[code]
            for train,num in ((x,y),(z,w)):
                old.append({'net_seg15':{'code':code,'label':label},
                    'train':{'code':train,'label':train},'value':num})
        sites=[{'eurostat_segment_code':"SIS10002_TEN + SIS20002_TEN"}]
        return sites,source2020,old
    def test_exact_code_and_station_label_two_way_2015_count(self):
        sites,meta,rows=self.make()
        joins,qa=matched_traffic_for_2020_site_segments(sites,meta,rows)
        self.assertEqual(qa['2015_complete_original_both_direction_site_intervals'],1)
        self.assertEqual(joins[sites[0]['eurostat_segment_code']]['goods'],4950)
        self.assertEqual(joins[sites[0]['eurostat_segment_code']]['passenger'],13800)
    def test_mismatched_label_is_absent_not_filled_from_2020(self):
        sites,meta,rows=self.make()
        rows[1]['net_seg15']['label']='Not the same 2020 segment'
        joins,qa=matched_traffic_for_2020_site_segments(sites,meta,rows)
        self.assertEqual(len(joins),0)
        self.assertEqual(len(qa['original_2015_missing_or_drifted_exact_source_direction_rows']),1)
    def test_2015_missing_passenger_direction_rejected(self):
        sites,meta,rows=self.make()
        rows=[r for r in rows if not (r['net_seg15']['code']=='SIS20002_TEN' and r['train']['code']=='TRN_PAS')]
        joins,qa=matched_traffic_for_2020_site_segments(sites,meta,rows)
        self.assertEqual(len(joins),0)
if __name__=='__main__':unittest.main()

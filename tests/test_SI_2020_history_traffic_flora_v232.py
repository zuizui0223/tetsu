import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_SI_2020_traffic_vs_2022_2023_flora_exploratory_v232 import (
    rank,per_physical_segment,association_over_physical_segments,betas
)
class TestHistorical2020FreightNotFieldYear(unittest.TestCase):
    def test_rank_requires_four_independent_units(self):
        self.assertIsNone(rank([0,1,2],[1,2,3]))
        self.assertAlmostEqual(rank([1,2,3,4],[1,4,9,16]),1.)
        self.assertAlmostEqual(rank([1,2,3,4],[4,3,2,1]),-1.)
    def test_two_sites_on_same_interval_only_one_exposure(self):
        v=[]
        for i in range(3):
            v.append({
              'source_segment':'SIS10002 + SIS20002',
              'source_segment_endpoints':'Dobova - Sevnica / Sevnica - Dobova',
              'section':'10','year':2022,'source_original_track_type':'DOUBLE',
              'freight':8549,'passenger':11947,
              'annual_share_BS_minus_BT':.1*(i+1),
              'therophyte_share_BS_minus_BT':.1*(i+1),
              'cover_BS_minus_BT_in_original_source_units':-30,
              'species_richness_BS_minus_BT_UNEQUAL_AREAS':-3})
        seg=per_physical_segment(v)
        self.assertEqual(len(seg),1)
        self.assertEqual(seg[0]['source_n_vegetation_field_sites'],3)
        self.assertEqual(seg[0]['source_2020_goods_both_directions'],8549)
        self.assertAlmostEqual(seg[0]['observed_field_source_mean']['annual_share_BS_minus_BT'],.2)
    def test_historical_goods_vs_annual_rho_not_field_randomized(self):
        a=[]
        for i in range(7):
            a.append({
               'national_line':str(10+i),
               'source_2020_goods_both_directions':(i+1)*100,
               'source_2020_passenger_both_directions':(i+1)*200,
               'observed_field_source_mean':{
                 'annual_share_BS_minus_BT':i/10,
                 'therophyte_share_BS_minus_BT':i/12,
                 'cover_BS_minus_BT_in_original_source_units':-10*i,
                 'species_richness_BS_minus_BT_UNEQUAL_AREAS':-i}})
        result=association_over_physical_segments(a)
        self.assertAlmostEqual(result['annual_share_BS_minus_BT']['goods']['exploratory_spearman_on_segment_means'],1.)
        self.assertTrue(result['annual_share_BS_minus_BT']['goods']['NO_independent_biological_p_value'])
    def test_invalid_repeated_traffic_values_rejected(self):
        a={'source_segment':'x','source_segment_endpoints':'a-b',
           'section':'10','freight':900,'passenger':120,
           'year':2022,'source_original_track_type':'SINGLE',
           'annual_share_BS_minus_BT':0.2,'therophyte_share_BS_minus_BT':.2,
           'cover_BS_minus_BT_in_original_source_units':-2,
           'species_richness_BS_minus_BT_UNEQUAL_AREAS':1}
        with self.assertRaises(ValueError):
            per_physical_segment([a,{**a,'freight':920}])
if __name__=='__main__':unittest.main()

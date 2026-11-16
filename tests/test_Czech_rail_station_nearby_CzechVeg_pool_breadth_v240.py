import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_Czech_39station_true_local_nontransport_plant_LEDA_beta_v240 import (
    sample_local_plots,source_profile,source_beta_stats,sample_frame,
    draw_panel,source_species_to_traits
)
class TestTrueCzechLocalRailVsNontrafficExternalPools(unittest.TestCase):
    def setUp(self):
        self.v={f"Genus species{chr(97+i//26)}{chr(97+i%26)}":
                ((i%23)/7.,(i%17)/9.) for i in range(150)}
        self.station=[]
        names=sorted(self.v)
        for i in range(12):
            controls=[]
            for k in range(8):
                controls.append({'id':f'ref{i}_{k}',
                   'd':0.3+k/10., 'study':'ref_A' if k<6 else 'ref_B',
                   'taxa':set(names[(i*5+k*3)%120:(i*5+k*3)%120+20]),
                   'n_all_source_taxa':20})
            self.station.append({'id':f'station_{i}',
                'stretch':f'rail_{i%4}',
                'rail_source_trait_taxa':set(names[(i*9)%90:(i*9)%90+35]),
                'control_options':controls})
    def test_five_unique_real_control_ids_and_two_studies(self):
        s=sample_local_plots(self.station[0],random.Random(12))
        self.assertEqual(len(s),5)
        self.assertEqual(len({x['id'] for x in s}),5)
        self.assertEqual(len({x['study'] for x in s}),2)
    def test_same_exact_12_original_rail_stations_paired_with_controls(self):
        x=sample_frame(self.station,self.v,random.Random(17),rarefaction=False)
        self.assertEqual(x['n_real_rail_stations'],12)
        self.assertEqual(x['source_n_dependent_station_pairs'],66)
        self.assertEqual(x['source_reference_selected_plot_n'],60)
        self.assertTrue(math.isfinite(x['rail_minus_local_nontransport_CWM_beta']))
        self.assertTrue(-1<=x['rail_minus_local_nontransport_Jaccard']<=1)
    def test_rarefaction_matching_species_counts_both_arms(self):
        x=sample_frame(self.station,self.v,random.Random(17),rarefaction=True)
        self.assertTrue(math.isfinite(x['rail_minus_local_nontransport_CWM_beta']))
        self.assertGreaterEqual(x['absolute_rail_CWM_distance'],0)
        self.assertGreaterEqual(x['absolute_local_control_CWM_distance'],0)
    def test_same_rail_bio_sites_12_is_support_for_five_study_LOO(self):
        x=draw_panel(self.station,self.v,True,seed=7,draws=8)
        self.assertTrue(x['identified'])
        self.assertEqual(x['original_observed_Czech_author_stations_in_panel'],12)
        self.assertEqual(x['n_source_original_railway_stretches'],4)
    def test_fewer_than_ten_original_stations_cannot_estimate(self):
        x=draw_panel(self.station[:8],self.v,True,seed=7,draws=8)
        self.assertFalse(x['identified'])
    def test_CWM_source_station_self_identity(self):
        a=source_profile(self.station[0]['rail_source_trait_taxa'],self.v)
        b=source_profile(self.station[0]['rail_source_trait_taxa'],self.v)
        self.assertEqual(a['CWM'],b['CWM'])
        self.assertGreater(a['trace'],0)
if __name__=='__main__':unittest.main()

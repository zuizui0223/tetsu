import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from validate_Czech_rail_station_local_reference_no_plot_reuse_v241 import (
    globally_disjoint_controls,no_reuse_draw,evaluate
)
class TestAllOriginalReferencePlotIDsUniqueAcrossRailStations(unittest.TestCase):
    def setUp(self):
        self.names=[f'Genus spec{chr(97+i//26)}{chr(97+i%26)}' for i in range(130)]
        self.v={n:(math.sin(i/9),math.cos(i/13)) for i,n in enumerate(self.names)}
        self.stations=[]
        for i in range(16):
            unique=[{'id':f'st{i}_source_{k}',
                     'study':'published_survey_A' if k<4 else 'published_survey_B',
                     'taxa':set(self.names[(i*7+k*3)%100:(i*7+k*3)%100+15])}
                    for k in range(8)]
            common=[{'id':f'crowded_local_original_{k}',
                     'study':'published_survey_A' if k<7 else 'published_survey_B',
                     'taxa':set(self.names[(k*11)%100:(k*11)%100+15])}
                    for k in range(10)]
            self.stations.append({
              'id':f'original_rail_station_{i}',
              'stretch':f'stretch_{i%4}',
              'rail_source_trait_taxa':set(self.names[(i*7)%75:(i*7)%75+40]),
              'control_options':unique+common})
    def test_disjoint_controls_are_80_independent_plot_ids_from_min_two_papers(self):
        match,ntry=globally_disjoint_controls(self.stations,random.Random(2))
        self.assertIsNotNone(match)
        self.assertEqual(len(match),16)
        chosen=[p['id'] for arr in match.values() for p in arr]
        self.assertEqual(len(chosen),80)
        self.assertEqual(len(set(chosen)),80)
        self.assertTrue(all(len({p['study'] for p in arr})>=2 for arr in match.values()))
    def test_insufficient_unique_plots_fails_not_fake_independence(self):
        copy=[]
        for st in self.stations:
            same=[{'id':f'shared_{j}', 'study':'study1' if j<3 else 'study2',
                    'taxa':set(self.names[j:j+10])} for j in range(5)]
            copy.append({**st,'control_options':same})
        allocations,ntry=globally_disjoint_controls(copy,random.Random(1),max_attempts=8)
        self.assertIsNone(allocations)
    def test_original_16_site_unrarefied_and_rarefied_metrics(self):
        for rare in (False,True):
            x,t=no_reuse_draw(self.stations,self.v,random.Random(34),rare)
            self.assertIsNotNone(x)
            self.assertEqual(x['original_reference_plot_ids_across_all_railway_stations_reused_n'],0)
            self.assertEqual(x['source_n_dependent_station_pairs'],120)
            self.assertTrue(math.isfinite(x['rail_minus_local_nontransport_CWM_beta']))
    def test_failure_gate_with_only_six_source_stations(self):
        out=evaluate(self.stations[:6],self.v,True,7,draws=10)
        self.assertFalse(out['identified'])
if __name__=='__main__':unittest.main()

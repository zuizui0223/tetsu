import random,sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_Czech_39station_original_LEDA_height_seed_functional_beta_v236 import (
    stats,original_margins,run_chain
)
class TestTrueCzechNullMathematics(unittest.TestCase):
    def setUp(self):
        rng=random.Random(20261116)
        names=[f't{i:03d}' for i in range(100)]
        self.values={key:(rng.gauss(0,1),rng.gauss(0,1)) for key in names}
        self.sites=[]
        for i in range(39):
            zone=f'rail-stretch-{i%5}'
            keep=set(rng.sample(names,26))
            self.sites.append({'station':f'station{i:02d}','stretch':zone,
               'taxa':keep,'source_unequal_station_area_ha':None})
    def test_39_stations_have_741_dependent_pairs(self):
        x=stats(self.sites,self.values)
        self.assertGreater(x['mean_station_joint_trace'],0)
        self.assertGreater(x['all_station_pair_centroid_distance'],0)
        self.assertGreater(x['all_station_pair_taxon_Jaccard'],0)
        self.assertEqual(len(self.sites)*(len(self.sites)-1)//2,741)
    def test_within_stretch_and_station_richness_margins_fixed(self):
        before=original_margins(self.sites)
        samples,qa=run_chain(self.sites,self.values,seed=224,
                             draws=3,burnin=4500,gap=1000)
        self.assertEqual(len(samples),3)
        self.assertGreater(qa['source_original_rail_stretch_constrained_accepted_swaps'],1000)
        self.assertEqual(before,original_margins(self.sites))
        self.assertTrue(all(math.isfinite(row['within_minus_between_centroid_distance'])
                            for row in samples))
    def test_taxonomic_and_functional_null_do_not_share_units(self):
        x=stats(self.sites,self.values)
        self.assertGreaterEqual(x['all_station_pair_taxon_Jaccard'],0)
        self.assertLessEqual(x['all_station_pair_taxon_Jaccard'],1)
        self.assertGreaterEqual(x['all_station_pair_centroid_distance'],0)
    def test_original_stretch_is_not_equivalent_to_country_or_regional_natural_pool(self):
        self.assertEqual(len({x['stretch'] for x in self.sites}),5)
if __name__=='__main__':unittest.main()

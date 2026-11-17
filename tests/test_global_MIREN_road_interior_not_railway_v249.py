import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
    complete_paired_first_surveys,regional_comparison,source_profile
)
class TestMIRENRoadNotRailwayGlobalPaired(unittest.TestCase):
    def setUp(self):
        self.rng=random.Random(247)
        self.spec=[f'Genus taxa{chr(97+i//26)}{chr(97+i%26)}' for i in range(65)]
        self.traits={s:(self.rng.gauss(0,1),self.rng.gauss(0,1)) for s in self.spec}
        self.records={}
        for region in ('R1','R2','R3'):
            for road in ('RDA','RDB'):
                for i in range(6):
                    for year in (2015,2020):
                        for z in (1,3):
                            self.records[region,road,str(i),year,z]=set(
                                self.rng.sample(self.spec,15))
    def test_select_one_original_earliest_site_year_paired_across_zones(self):
        p,qa=complete_paired_first_surveys(self.records,self.traits)
        self.assertEqual(len(p),3)
        self.assertEqual(qa['qualified_unique_original_stable_road_site_first_year_only'],36)
        self.assertEqual(qa['n_main_roadside_vs_interior_paired_site_records'],36)
        self.assertTrue(all(all(site['year']==2015 for road,sites in group for site in sites) for group in p.values()))
    def test_same_original_sitepair_road_vs_interior_functional_beta_independent_metric(self):
        p,qa=complete_paired_first_surveys(self.records,self.traits)
        for region,groups in p.items():
            r=regional_comparison(groups,self.traits,region)
            self.assertEqual(r['n_original_roads'],2)
            self.assertEqual(r['n_original_spatial_sites'],12)
            self.assertTrue(math.isfinite(r['source_2D_CWM_beta_road_minus_interior']))
            self.assertTrue(-1<=r['source_species_Jaccard_road_minus_interior']<=1)
            rare=regional_comparison(groups,self.traits,region,random.Random(12),rarefy=True)
            self.assertTrue(math.isfinite(rare['source_alpha_2D_trace_road_minus_interior']))
    def test_missing_interior_is_NOT_imputed_as_absent_species(self):
        records={k:v for k,v in self.records.items() if k[-1]!=3}
        p,qa=complete_paired_first_surveys(records,self.traits)
        self.assertEqual(len(p),0)
        self.assertEqual(qa['n_main_roadside_vs_interior_paired_site_records'],0)
    def test_region_with_only_one_source_road_not_global_region_replication(self):
        records={k:v for k,v in self.records.items() if k[1]=='RDA'}
        p,qa=complete_paired_first_surveys(records,self.traits)
        self.assertEqual(len(p),0)
    def test_source_zone1_and_zone3_have_equal_technical_richness_when_rarefied(self):
        p,qa=complete_paired_first_surveys(self.records,self.traits)
        x=regional_comparison(p['R1'],self.traits,'R1',random.Random(99),True)
        self.assertEqual(x['source_site_LEDA_species_median_road'],15)
        self.assertEqual(x['source_site_LEDA_species_median_interior'],15)
if __name__=='__main__':unittest.main()

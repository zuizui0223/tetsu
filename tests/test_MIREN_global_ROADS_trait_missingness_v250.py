import random,sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import (
    valid_region_panels,describe_panel
)
class TestTrueGlobalRoadThresholdNotFalseRailCountry(unittest.TestCase):
    def setUp(self):
        rng=random.Random(12)
        sp=[f'Genus specimen{chr(97+i//26)}{chr(97+i%26)}' for i in range(125)]
        self.values={x:(rng.uniform(-2,2),rng.uniform(-3,3)) for x in sp[:85]}
        self.sample={}
        for region in ("AUT","IND","CHE"):
            for road in ("RD1","RD2"):
                for site in range(5):
                    for year in (2014,2019):
                        for zone in (1,3):
                            # 15/20 genuine original plants have LEDA traits
                            if zone==1:
                                source=set(rng.sample(sp[:85],15)+rng.sample(sp[85:],5))
                            else:
                                source=set(rng.sample(sp[:85],15)+rng.sample(sp[85:],5))
                            self.sample[region,road,str(site),year,zone]=source
    def test_minimum_baseline_fixed_two_roads_and_first_survey(self):
        groups,qa=valid_region_panels(self.sample,self.values,5,.25,2007)
        self.assertEqual(qa['source_regional_2road_4site_support_n'],3)
        self.assertEqual(qa['source_eligible_spatial_sites'],30)
        self.assertTrue(all(x['year']==2014 for group in groups.values() for road,items in group for x in items))
    def test_50percent_coverage_and_ten_trait_species_can_be_explicitly_tested(self):
        groups,qa=valid_region_panels(self.sample,self.values,10,.5,2007)
        self.assertEqual(qa['source_eligible_spatial_sites'],30)
        result=describe_panel(groups,self.values)
        self.assertTrue(result['identified'])
        self.assertEqual(result['n_source_original_region_codes'],3)
    def test_2015plus_only_source_real_second_year(self):
        groups,qa=valid_region_panels(self.sample,self.values,5,.25,2015)
        self.assertEqual(qa['source_eligible_spatial_sites'],30)
        self.assertTrue(all(x['year']==2019 for group in groups.values() for road,items in group for x in items))
    def test_no_road_interior_plot_pair_fails_source_support_not_zero_effect(self):
        data={k:v for k,v in self.sample.items() if k[-1]==1}
        groups,qa=valid_region_panels(data,self.values,5,.25,2007)
        result=describe_panel(groups,self.values)
        self.assertFalse(result['identified'])
        self.assertTrue(result['not_zero_functional_road_effect'])
    def test_single_source_road_per_region_is_not_region_replication(self):
        data={k:v for k,v in self.sample.items() if k[1]=="RD1"}
        groups,qa=valid_region_panels(data,self.values,5,.25,2007)
        self.assertEqual(len(groups),0)
if __name__=='__main__':unittest.main()

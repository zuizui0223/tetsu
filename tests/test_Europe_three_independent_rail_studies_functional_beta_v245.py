import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import (
    traits_join,source_cohort,get_metrics,regional_margin,one_null,METRICS
)
class TestIndependentThreeStudyRailSources(unittest.TestCase):
    def setup_data(self):
        names=[f"Genus spec{chr(97+i//26)}{chr(97+i%26)}" for i in range(400)]
        h={n:1+i*.08 for i,n in enumerate(names)}
        se={n:1+2*i*.3 for i,n in enumerate(names)}
        rng=random.Random(20261117)
        studies={}
        for label,n,groups in (("CZ",39,5),("IT",28,2),("HR_BA",15,3)):
            sites=[]
            for i in range(n):
                taxa=set(rng.sample(names,25))
                sites.append({'station':f'{label}{i}',
                    'source_names':taxa,'region_group':f'group{i%groups}'})
            studies[label]=sites
        return studies,h,se
    def test_common_pooled_joint_height_seed_axes(self):
        studies,h,se=self.setup_data()
        value,q=traits_join(studies,h,se)
        self.assertEqual(len(value),400)
        self.assertEqual(q['source_original_LEDA_both_canopy_height_seed_mass_species_n'],400)
        self.assertTrue(all(len(x)==2 for x in value.values()))
    def test_true_82_stations_in_three_original_studies(self):
        studies,h,se=self.setup_data()
        values,q=traits_join(studies,h,se)
        real,qc=source_cohort(studies,values)
        self.assertEqual(sum(len(x) for x in real.values()),82)
        self.assertEqual(len(real),3)
        self.assertEqual(qc['CZ']['source_total_author_station_n'],39)
        self.assertEqual(qc['IT']['source_total_author_station_n'],28)
        self.assertEqual(qc['HR_BA']['source_total_author_station_n'],15)
    def test_original_group_fixed_margin_one_source_chain(self):
        studies,h,se=self.setup_data()
        v,q=traits_join(studies,h,se)
        real,qa=source_cohort(studies,v)
        for study in real:
            before=regional_margin(real[study])
            r,mix=one_null(real[study],v,seed=2,draws=4,burnin=3500,gap=1500)
            self.assertGreater(mix['accepted_true_source_two_by_two_incidence_swaps'],1000)
            self.assertEqual(regional_margin(real[study]),before)
            self.assertEqual(len(r),4)
            self.assertEqual(set(r[0]),set(METRICS))
    def test_species_frequency_study_and_internal_CWM_different_metrics(self):
        studies,h,se=self.setup_data()
        v,q=traits_join(studies,h,se)
        real,qa=source_cohort(studies,v)
        for k,s in real.items():
            r=get_metrics(s,v)
            self.assertGreater(r['functional_CWM_centroid_distance_all_stationpairs'],0)
            self.assertGreater(r['source_station_2D_trait_trace_mean'],0)
            self.assertGreaterEqual(r['source_matched_species_Jaccard_all_stationpairs'],0)
            self.assertLessEqual(r['source_matched_species_Jaccard_all_stationpairs'],1)
if __name__=="__main__":
    unittest.main()

import math,sys,unittest,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_SI_original_82site_BS_BT_continuous_functional_beta_v229 import (
    site_zone,zone_pair_stat,pair_stat_value,summary_by_independent_sections,
    metric_aggregate,rarefy_site_once,difference_for_existing_pairs
)
class TestFunctionalBeta(unittest.TestCase):
    def setUp(self):
        self.v={str(i):(i/5.,(-1)**i*(i+1)/7.) for i in range(10)}
        self.g={str(i):int(i%2==0) for i in range(10)}
        self.p={k:site_zone(set(map(str,range(k,k+5))),self.v,self.g) for k in (0,1,2,3,4,5)}
    def test_same_input_identical_site_centroid(self):
        self.assertAlmostEqual(zone_pair_stat(self.p[0],self.p[0])['functional_centroid_euclidean'],0.)
        self.assertAlmostEqual(zone_pair_stat(self.p[0],self.p[0])['taxon_Jaccard_dissimilarity'],0.)
    def test_all_metrics_symmetry(self):
        p=zone_pair_stat(self.p[0],self.p[1]);q=zone_pair_stat(self.p[1],self.p[0])
        self.assertEqual(p,q)
    def test_label_flip_difference_is_antisymmetric(self):
        pair={'normals':{'functional_centroid_euclidean':1.5},
              'crosses':{'functional_centroid_euclidean':0.75}}
        m='functional_centroid_euclidean'
        self.assertAlmostEqual(pair_stat_value(pair,m),1.5)
        self.assertAlmostEqual(pair_stat_value(pair,m,True,True),-1.5)
        self.assertAlmostEqual(pair_stat_value(pair,m,False,True),0.75)
        self.assertAlmostEqual(pair_stat_value(pair,m,True,False),-0.75)
    def test_section_pair_aggregate_not_site_n(self):
        pairs=[]
        for i,j,a in ((0,1,.1),(1,2,-.1),(0,2,.3)):
            pairs.append({'i':i,'j':j,'sec_i':str(i),
                'sec_j':str(j),'section_pair':tuple(sorted((str(i),str(j)))),
                'same_section':False,'band':'20-40','normals':{'functional_centroid_euclidean':a},
                'crosses':{'functional_centroid_euclidean':a}})
        result=summary_by_independent_sections(pairs,'functional_centroid_euclidean')
        self.assertEqual(result['n_distinct_cross_source_section_pairs'],3)
        self.assertAlmostEqual(result['mean_equal_cross_sectionpair_delta'],.1)
        self.assertAlmostEqual(metric_aggregate(pairs,'functional_centroid_euclidean'),.1)
    def test_same_site_richness_rarefaction_keeps_pair(self):
        one={'section':'10','BS':self.p[0],'BT':self.p[1],
             'lat':46.0,'lon':14.0}
        a=rarefy_site_once([one],self.v,self.g,random.Random(4))[0]
        self.assertEqual(len(a['BS']['species']),5)
        self.assertEqual(len(a['BT']['species']),5)
    def test_convergence_direction_is_negative(self):
        # BS identical across sites; BT centroids differ
        aa=site_zone({'0','1','2','3','4'},self.v,self.g)
        bb=site_zone({'4','5','6','7','8'},self.v,self.g)
        self.assertLess(
            zone_pair_stat(aa,aa)['functional_centroid_euclidean']-
            zone_pair_stat(aa,bb)['functional_centroid_euclidean'],0)
if __name__=='__main__':unittest.main()

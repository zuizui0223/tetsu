import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_Slovenia_annual_trait_variance_decomposition_v227 import (
    component,pair_metrics,mixture_only_proxy,site_contrast,bysection_aggregate,group_lifespan
)
class TestSourceLifeHistoryTraceDecomposition(unittest.TestCase):
    def setUp(self):
        self.v={'a':(-1.,-2.),'b':(0.,0.),'c':(1.,2.),'d':(2.,-3.),
                'e':(-2.,3.),'f':(0.5,0.5),'g':(3.,-1.)}
        self.g={k:('annual' if k in ('a','c','e','g') else 'other') for k in self.v}
    def test_variance_decomposition_identity(self):
        v=component(set(self.v),self.v,self.g)
        self.assertAlmostEqual(v['total_trace'],v['within_total']+v['between_lifespan_contribution'],places=11)
        self.assertAlmostEqual(v['within_total'],v['within_annual_contribution']+v['within_nonannual_contribution'])
        self.assertAlmostEqual(v['annual_share'],4/7)
    def test_one_group_yields_between_zero(self):
        c=component({'a','c','e'},self.v,self.g)
        self.assertEqual(c['between_lifespan_contribution'],0.)
        self.assertAlmostEqual(c['within_total'],c['total_trace'])
    def test_pair_swap_sign(self):
        a={'a','b','c','d','e'}
        b={'a','b','c','f','g'}
        d=pair_metrics(a,b,self.v,self.g)
        e=pair_metrics(b,a,self.v,self.g)
        self.assertTrue(all(math.isclose(d[k],-e[k]) for k in d))
    def test_mixture_only_zero_when_group_shares_equal(self):
        a={'a','b','c','d'} # 2 annual 2 other
        b={'a','b','e','f'} # 2 annual 2 other
        self.assertAlmostEqual(mixture_only_proxy(a,b,self.v,self.g),0.)
    def test_preserve_original_source_shared_zone_and_annual_numbers(self):
        a={'a','b','c','d','e'}
        b={'a','b','c','f','g'}
        src={'z':[{'section':'10','BS':a,'BT':b,'known_joint_annual_nBS':3,
            'known_joint_annual_nBT':3,'known_joint_nBS':5,
            'known_joint_nBT':5,'shared_source_taxa':3}]}
        observed=site_contrast(src,self.v,self.g)
        vals=[site_contrast(src,self.v,self.g,random.Random(i),True)['z']['delta']['total_trace']
              for i in range(40)]
        self.assertTrue(all(math.isfinite(x) for x in vals))
        self.assertAlmostEqual(
            observed['z']['delta']['total_trace'],
            observed['z']['delta']['within_total']+observed['z']['delta']['between_lifespan_contribution'])
    def test_section_equal_vs_location_equal(self):
        a={'a','b','c'}
        b={'a','d','e'}
        src={'one':[{'section':'10','BS':a,'BT':b}],
             'two':[{'section':'20','BS':b,'BT':a}]}
        by=site_contrast(src,self.v,self.g)
        sect,score=bysection_aggregate(by)
        self.assertEqual(len(sect),2)
        self.assertAlmostEqual(score['total_trace'],0.)
    def test_mixed_lifespan_is_missing_in_strict_sensitivity(self):
        tr={'x':{'Lifespan.categories':'annual/biennial'},
            'y':{'Lifespan.categories':'annual'},
            'z':{'Lifespan.categories':None}}
        self.assertEqual(group_lifespan('x',tr),'other')
        self.assertIsNone(group_lifespan('x',tr,exclude_mixed=True))
        self.assertEqual(group_lifespan('y',tr),'annual')
        self.assertIsNone(group_lifespan('z',tr))
if __name__=='__main__':unittest.main()

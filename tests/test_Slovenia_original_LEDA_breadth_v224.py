import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    variance,components,difference,ztransform,site_effects,section_estimates,
    finite_summary,quart
)
class TestStrictSIBreadth(unittest.TestCase):
    def setUp(self):
        self.val={"a":(-1.,-2.),"b":(0.,0.),"c":(1.,2.),"d":(2.,3.),
                  "e":(-2.,-3.),"f":(0.5,0.5),"g":(3.,-1.)}
    def test_variance_population_not_fake_sample(self):
        self.assertAlmostEqual(variance([0,2]),1.)
        self.assertAlmostEqual(variance([2,2,2]),0.)
        self.assertIsNone(variance([1.]))
    def test_trace_is_sum_of_dimension_variances(self):
        x=components({"a","b","c"},self.val)
        self.assertAlmostEqual(x["joint_trace"],x["height_var"]+x["seed_var"])
        self.assertAlmostEqual(x["height_mean"],0.)
    def test_zone_order_and_equal_set(self):
        self.assertEqual(difference({"a","b","c"},{"a","b","c"},self.val)["joint_trace"],0.)
        a=difference({"a","b"},{"a","b","d"},self.val)["joint_trace"]
        self.assertAlmostEqual(a,-difference({"a","b","d"},{"a","b"},self.val)["joint_trace"])
    def test_permute_local_pool_overlaps_and_zone_sizes(self):
        # Two original same-side zones share a & b, exclusive c/d/e/f.
        pair={"section":"1","BS":{"a","b","c","d"},"BT":{"a","b","e","f"},
              "nBS":4,"nBT":4,"n_shared":2}
        sites={"s1":[pair]}
        observed=site_effects(sites,self.val)
        self.assertEqual(observed["s1"]["mean_shared_species"],2)
        rng=random.Random(42)
        vals=[site_effects(sites,self.val,draw=i,rg=rng)["s1"]["delta"]["joint_trace"] for i in range(150)]
        self.assertTrue(all(math.isfinite(x) for x in vals))
        self.assertGreater(len(set(round(x,6) for x in vals)),2)
    def test_shared_only_null_degenerate(self):
        one={"section":"1","BS":{"a","b","c"},"BT":{"a","b","c"},
              "nBS":3,"nBT":3,"n_shared":3}
        sites={"s":[one]}
        a=site_effects(sites,self.val)["s"]["delta"]
        for seed in range(10):
            b=site_effects(sites,self.val,draw=seed,rg=random.Random(seed))["s"]["delta"]
            self.assertEqual(a,b)
    def test_rarefaction_equal_counts_no_richness_artifact(self):
        one={"section":"1","BS":{"a","b","c","d","e"},"BT":{"a","b","c"},
              "nBS":5,"nBT":3,"n_shared":3}
        sites={"s":[one]}
        val=site_effects(sites,self.val,rg=random.Random(34),rarefy=True)
        self.assertTrue(math.isfinite(val["s"]["delta"]["joint_trace"]))
        self.assertEqual(val["s"]["mean_BS_species"],5)
    def test_sections_not_sides_independent(self):
        prs={"s1":[{"section":"1","BS":{"a","b"},"BT":{"a","b","c"},"nBS":2,"nBT":3,"n_shared":2},
                    {"section":"1","BS":{"a","b"},"BT":{"a","b","d"},"nBS":2,"nBT":3,"n_shared":2}],
             "s2":[{"section":"2","BS":{"a","b"},"BT":{"a","b","e"},"nBS":2,"nBT":3,"n_shared":2}]}
        loc=site_effects(prs,self.val)
        by,avg=section_estimates(loc)
        self.assertEqual(len(loc),2)
        self.assertEqual(len(by),2)
        self.assertAlmostEqual(avg["joint_trace"],sum(x["joint_trace"] for x in by.values())/2)
    def test_zscore_uses_species_not_record_frequency(self):
        with self.assertRaises(ValueError):
            ztransform({"t1":"A a","t2":"B b","t3":"C c"},
             {"A a":1,"B b":2,"C c":8},{"A a":1,"B b":4,"C c":16})
    def test_interval_order(self):
        z=finite_summary(list(range(100)),observed=0)
        self.assertLessEqual(z["q025"],z["q975"])
        self.assertEqual(z["observed_less_equal_count"],1)
        self.assertAlmostEqual(quart([1,2,3,4],.5),2.5)
if __name__=="__main__":unittest.main()

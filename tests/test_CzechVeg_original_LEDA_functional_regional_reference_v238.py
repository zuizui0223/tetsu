import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from analyze_CzechVeg_original_LEDA_rail_text_vs_abiotic_reference_functional_v238 import (
    exact_source_name,trait_db,profile_plots,candidate_arms,
    pairwise_mean,one_draw,eval_panel
)
class TestOriginalCzechVegExternalReferenceV238(unittest.TestCase):
    def setUp(self):
        self.h={f"Genus taxon{i}":float(i+1) for i in range(120)}
        self.s={f"Genus taxon{i}":float((i+1)**2) for i in range(120)}
        self.val,_=trait_db(self.h,self.s,set(self.h))
        self.plants={
          "rail_A":set(list(self.h)[:15]),
          "rail_B":set(list(self.h)[15:30]),
          "rail_C":set(list(self.h)[30:45]),
          "ctrl_1":set(list(self.h)[45:60]),
          "ctrl_2":set(list(self.h)[60:75]),
          "ctrl_3":set(list(self.h)[75:90]),
          "ctrl_4":set(list(self.h)[90:105]),
          "ctrl_5":set(list(self.h)[105:120])
        }
    def test_name_matching_preserves_species_not_aggregate(self):
        self.assertEqual(exact_source_name("Genus species"),"Genus species")
        self.assertIsNone(exact_source_name("Genus species agg."))
        self.assertIsNone(exact_source_name("Genus species cf."))
        self.assertIsNone(exact_source_name("Genus species group"))
    def test_original_all_plants_one_per_source_species(self):
        p,q=profile_plots(self.plants,self.val)
        self.assertEqual(len(p),8)
        self.assertTrue(all(abs(x["source_dual_trait_coverage"]-1.)<1e-10 for x in p.values()))
        self.assertEqual(len(p["rail_A"]["exact_species"]),15)
    def test_three_original_regions_required(self):
        p,_=profile_plots(self.plants,self.val)
        frame={"A":[("rail_A",{"altitude":["ctrl_1"]*5,"all":["ctrl_2"]*5})],
               "B":[("rail_B",{"altitude":["ctrl_3"]*5,"all":["ctrl_4"]*5})]}
        answer=eval_panel(frame,p,seed=2,draws=5)
        self.assertFalse(answer["identified"])
    def test_CWM_euclidean_nonnegative_and_Jaccard_bounded(self):
        p,_=profile_plots(self.plants,self.val)
        a=[p["rail_A"],p["rail_B"],p["rail_C"]]
        self.assertGreaterEqual(pairwise_mean(a,"CWM"),0.)
        self.assertLessEqual(pairwise_mean(a,"species"),1.)
        self.assertEqual(pairwise_mean([a[0]]*3,"CWM"),0.)
    def test_independent_publication_and_richness_control(self):
        p,_=profile_plots(self.plants,self.val)
        m={"rail_A":{"region":"A","study":"rail_paper","year":2020,"strict":"rail"},
           "ctrl_1":{"region":"A","study":"independent_nature_paper","year":2020,"strict":None},
           "ctrl_2":{"region":"A","study":"independent_nature_paper","year":2020,"strict":None}}
        controls={"all":{"rail_A":["ctrl_1"]*5},
                  "elev100":{"rail_A":["ctrl_2"]*5}}
        frame,qa=candidate_arms(m,controls,p,"all_source_rail_text")
        self.assertEqual(qa["source_rail_text_plots_matched_to_two_types"],1)
        self.assertEqual(len(frame["A"]),1)
        q=one_draw({"A":frame["A"],"B":frame["A"],"C":frame["A"]},p,["A","B","C"],random.Random(1))
        self.assertAlmostEqual(q["rail_minus_altitude_cwm_distance"],0.)
    def test_strict_label_excludes_context_only(self):
        p,_=profile_plots(self.plants,self.val)
        m={"rail_A":{"region":"A","study":"rail_paper","year":2020,"strict":None},
           "ctrl_1":{"region":"A","study":"independent_paper","year":2020,"strict":None}}
        d={"all":{"rail_A":["ctrl_1"]*5},"elev100":{"rail_A":["ctrl_1"]*5}}
        frame,qa=candidate_arms(m,d,p,"strict_original_HABITAT_rail_only")
        self.assertEqual(len(frame),0)
if __name__=="__main__":unittest.main()

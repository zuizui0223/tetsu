import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Kashmir_C4_Cynodon_specificity_v339b import (
   CY,C4_GRASS,C3_GRASS,mean_fraction_delta,exact_grass_pathway_label_shuffle
)


class PostHocGuildSpecificity(unittest.TestCase):
    def setUp(self):
        self.names=["cynodon dactylon","poa annua","setaria viridis","avena sativa"]
        self.rows=[({0:1.0,2:1.0},{1:1.0,3:1.0}),
                   ({0:9.0,1:1.0},{1:2.0,2:1.0})]

    def test_exact_grass_label_shuffle(self):
        res=exact_grass_pathway_label_shuffle(self.rows,self.names,
                  {"cynodon dactylon","setaria viridis"},
                  {"poa annua","avena sativa"})
        self.assertEqual(res["n_exact_taxon_label_assignments"],6)
        self.assertGreater(res["observed_C4_grass_track_minus_interior_fraction"],0)
        self.assertTrue(0<=res["two_sided_exact_label_shuffle_fraction_centered"]<=1)

    def test_Cynodon_effect_is_exact_difference_from_full_C4(self):
        c4={"cynodon dactylon","setaria viridis"}
        total=mean_fraction_delta(self.rows,self.names,c4)
        cy=mean_fraction_delta(self.rows,self.names,{CY})
        others=mean_fraction_delta(self.rows,self.names,c4-{CY})
        self.assertAlmostEqual(total,cy+others)

    def test_abundance_composition_differs_from_presence(self):
        group={"cynodon dactylon"}
        a=mean_fraction_delta(self.rows,self.names,group)
        b=mean_fraction_delta(self.rows,self.names,group,True)
        self.assertNotAlmostEqual(a,b)


if __name__=="__main__":unittest.main()

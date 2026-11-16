import importlib.util, unittest, random
from pathlib import Path
from collections import Counter
p=Path(__file__).resolve().parents[1]/"scripts/validate_CzechVeg_heldout_life_diaspore_trait_label_placebo_v199.py"
sp=importlib.util.spec_from_file_location("v199",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

class TestTrueSourceTraitAblationAndShuffledLabelPlacebo(unittest.TestCase):
    def synthetic(self):
        # Synthetic ONLY: no biological field data in this unit test.
        taxa=[f"T{i:02d} t{i}" for i in range(18)]
        life={s:("therophyte" if i%3==0 else "hemicryptophyte" if i%3==1 else "geophyte")
              for i,s in enumerate(taxa)}
        disp={s:("Epilobium" if i%3==0 else "Allium" if i%3==1 else "Phragmites")
              for i,s in enumerate(taxa)}
        blocks={};plots={}
        for paper_i in range(12):
            p=f"paper_{paper_i:02d}"
            block=(p,f"region_{paper_i%4}",2010,"gt4_to_16")
            arms={"rail":[],"road":[]}
            for group in ("rail","road"):
                for jj in range(2):
                    pid=f"{p}_{group}_{jj}"
                    arms[group].append(pid)
                    series=taxa[0:10] if group=="rail" else taxa[6:16]
                    plots[pid]=[(s,disp[s],life[s]) for s in series]
            blocks[block]=arms
        return blocks,plots
    def test_feature_ablation_preserves_effort_and_trait_masks(self):
        plot=[("A a","Epilobium","therophyte"),("B b","Allium","hemicryptophyte"),
              ("C c","Epilobium","therophyte"),("D d",None,"hemicryptophyte"),
              ("E e",None,None),("F f","Allium","therophyte")]
        x=m.frame_records if False else m.feature
        from validate_CzechVeg_source_publication_holdout_guild_vs_taxon_v198 import site_characteristics
        site=site_characteristics(plot)
        base=m.feature(site,"effort")
        life=m.feature(site,"lifeform")
        dispersal=m.feature(site,"dispersal")
        both=m.feature(site,"combined")
        self.assertTrue(all(k in both for k in life))
        self.assertTrue(all(k in both for k in dispersal))
        self.assertEqual(len(base),3)
        self.assertFalse(any(k.startswith("disp:") for k in life))
        self.assertFalse(any(k.startswith("life:") for k in dispersal))
    def test_real_trait_pair_and_missingness_preserved_in_null(self):
        life={"T1":"therophyte","T2":"geophyte","T3":"hemicryptophyte",
              "T4":"therophyte","T5":"hemicryptophyte","T6":"geophyte"}
        disp={"T1":"Epilobium","T2":"Allium","T3":"Cornus",
              "T4":"Epilobium","T5":"Allium","T7":"Phragmites"}
        counts=Counter({k:5 for k in set(life)|set(disp)})
        for mode in m.SHUFFLE_DESIGNS:
            x,meta=m.shuffled_taxon_mapping(life,disp,counts,random.Random(8),mode)
            self.assertEqual(set(x["life"]),set(life))
            self.assertEqual(set(x["disp"]),set(disp))
            self.assertEqual(Counter(x["life"].values()),Counter(life.values()))
            self.assertEqual(Counter(x["disp"].values()),Counter(disp.values()))
            self.assertTrue(meta["trait_missingness_mask_originally_preserved"] if "trait_missingness_mask_originally_preserved" in meta else meta["pladias_trait_coverage_unchanged"])
    def test_fixed_paper_heldout_ablation_and_placebo(self):
        blocks,plots=self.synthetic()
        folds,_=m.split_source_studies(blocks)
        panel=m.analyze_panel(blocks,plots,folds,"both_trait50_min5",permutations=2)
        self.assertTrue(panel["identified"])
        self.assertEqual(panel["n_original_source_heldout_publications"],12)
        self.assertEqual(set(panel["original_study_equal_holdout_model_auc"]),set(m.MODELS))
        for result in panel["actual_trait_mapping_vs_Pladias_label_shuffled_nulls"].values():
            self.assertEqual(result["n_synthetic_global_taxon_trait_label_reassignments"],2)
            self.assertTrue(result["trait_missingness_mask_originally_preserved"])
            self.assertTrue(result["original_trait_marginal_and_joint_distribution_preserved"])
        self.assertTrue(panel["same_source_publications_and_plot_ids_for_all_four_models"])
if __name__=="__main__":unittest.main()

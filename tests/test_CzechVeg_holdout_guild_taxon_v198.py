import importlib.util, unittest, random
from pathlib import Path
from collections import Counter
p=Path(__file__).resolve().parents[1]/"scripts/validate_CzechVeg_source_publication_holdout_guild_vs_taxon_v198.py"
sp=importlib.util.spec_from_file_location("v198",p)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

class TestOriginalFieldSourcePaperHeldoutGuildVsTaxon(unittest.TestCase):
    def toy(self):
        blocks={};plots={}
        taxa=["Echium vulgare","Poa compressa","Artemisia vulgaris",
              "Conyza canadensis","Capsella bursa-pastoris",
              "Bromus sterilis","Setaria viridis","Plantago major",
              "Taraxacum officinale","Stellaria media",
              "Bromus tectorum","Lolium perenne"]
        life={s:("therophyte" if i%2 else "hemicryptophyte")
              for i,s in enumerate(taxa)}
        disp={s:("Epilobium" if i%2 else "Allium") for i,s in enumerate(taxa)}
        for i in range(12):
            study="source_paper_"+str(i)
            key=(study,"r"+str(i%4),2010,"gt4_to_16")
            arms={"rail":[],"road":[]}
            for label in ("rail","road"):
                for j in range(3):
                    pid=f"{study}_{label}_{j}"
                    arms[label].append(pid)
                    ts=taxa[:8] if label=="rail" else taxa[4:12]
                    plots[pid]=[(t,disp[t],life[t]) for t in ts]
            blocks[key]=arms
        return blocks,plots
    def test_study_folds_and_exact_source_site_metrics(self):
        blocks,plots=self.toy()
        folds,masses=m.split_source_studies(blocks,nfold=4)
        self.assertEqual(len(folds),4)
        self.assertEqual(sum(len(s) for s in folds),12)
        for panel in m.PANEL_THRESHOLDS:
            selected,site,orig,frame=m.frame_records(blocks,plots,panel)
            self.assertEqual(frame["original_blocks_retained"],12)
            self.assertEqual(len(orig),72)
            self.assertEqual(len({orig[k]["study"] for k in orig}),12)
            for f in range(4):
                out=m.run_fold(selected,site,orig,folds,f,permutations=2)
                self.assertTrue(out["identified"])
                self.assertEqual(out["heldout_paper_n"],3)
                self.assertTrue(out["all_three_models_same_exact_heldout_source_sites"])
                for mdl in m.MODELS:
                    self.assertGreaterEqual(out["source_publication_holdout_model_AUC"][mdl]["study_balanced_auc"],0)
                    self.assertLessEqual(out["source_publication_holdout_model_AUC"][mdl]["study_balanced_auc"],1)
    def test_source_trait_coverage_no_rare_name_imputation(self):
        blocks,plots=self.toy()
        a=m.site_characteristics(plots[next(iter(plots))])
        self.assertEqual(a["n"],8)
        self.assertEqual(a["frac_life"],1)
        self.assertTrue(m.site_pass(a,"both_trait50_min5"))
        self.assertFalse(m.site_pass({"n":8,"n_life":0,"n_disp":0,
             "frac_life":0,"frac_disp":0},"both_trait50_min5"))
        base=m.features(a,"effort");guild=m.features(a,"guild")
        tax=m.features(a,"taxon")
        self.assertTrue(any(k.startswith("life:") for k in guild))
        self.assertTrue(any(k.startswith("disp:") for k in guild))
        self.assertTrue(any(k.startswith("species:") for k in tax))
        self.assertTrue(all(not x.startswith("species:") for x in guild))
        self.assertEqual(len(base),3)
    def test_source_block_auc_no_study_label_leak(self):
        self.assertEqual(m.source_block_auc([1,2],[0,0]),1)
        self.assertEqual(m.source_block_auc([1,1],[1,1]),.5)
        self.assertEqual(m.source_block_auc([0,0],[1,1]),0)
        blocks,plots=self.toy()
        folds,_=m.split_source_studies(blocks)
        train={k:v for k,v in blocks.items() if k[0] not in folds[0]}
        labels=m.permutation_of_training_labels(train,{},random.Random(12))
        for key,b in train.items():
            self.assertEqual(sum(labels[k] for k in b["rail"]+b["road"]),3)
        self.assertTrue(all(k[0] not in folds[0] for k in train))
if __name__=="__main__":unittest.main()

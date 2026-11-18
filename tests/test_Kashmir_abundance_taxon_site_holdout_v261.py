"""v26.1 synthetic QA: site-blocked cross-validation and taxon identity ablation.

Artificial data ONLY; not biological estimates. Numerical source SHA is gated
by original-source workflow before genuine Kashmir analysis.
"""
from __future__ import annotations
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import validate_Kashmir_abundance_after_species_site_holdout_v261 as v


def fake_rows(n_site=27,n_taxon=5):
    rows=[]
    for i in range(n_site):
        group="station" if i<15 else "between"
        for zone in v.ZONES:
            for t in range(n_taxon):
                position=(t+2*i+v.ZONES.index(zone))%n_taxon
                rank=2*position/(n_taxon-1)-1
                # Outcome depends on source rank, not species. This is not
                # an independent study and the values are synthetic.
                outcome=int(position>=3)
                rows.append({
                    "site":f"site{i:02d}","group":group,"zone":zone,
                    "taxon":t,"x2014_rank":rank,
                    "x2014_raw":position+1.,
                    "re_recorded_2017":outcome
                })
    return rows


class SiteHoldoutTests(unittest.TestCase):
    def test_source_within_plot_midrank_ties_and_singleton(self):
        z=v.midranks({1:5.,2:5.,3:10.,4:20.})
        self.assertAlmostEqual(z[1],-2/3)
        self.assertAlmostEqual(z[2],-2/3)
        self.assertAlmostEqual(z[3],1/3)
        self.assertAlmostEqual(z[4],1.)
        self.assertEqual(v.midranks({7:9.}),{7:0.})
        with self.assertRaises(ValueError):v.midranks({})

    def test_fixed_design_dimension_and_training_per_site(self):
        rows=fake_rows()
        y,models,ids=v.make_design(rows)
        self.assertEqual(len(ids),27)
        self.assertEqual(set(models),set(v.MODELS))
        self.assertEqual(models[v.MODELS[0]].shape,(len(rows),8))
        self.assertEqual(models[v.MODELS[2]].shape,(len(rows),8+123))
        self.assertEqual(models[v.MODELS[4]].shape,(len(rows),8+123+2))
        self.assertEqual(y.sum(),27*4*2)
        self.assertEqual(set(rows[i]["site"] for i in ids["site02"]),{"site02"})

    def test_heldout_source_rank_adds_information_after_taxon(self):
        rows=fake_rows()
        out=v.fold_evaluation(rows,1.0)
        a=out["source_predictions_by_group_and_habitat"]["station"]["Track"]
        b=out["source_predictions_by_group_and_habitat"]["between"]["Track"]
        self.assertEqual(a["n_original_source_heldout_site_folds"],15)
        self.assertEqual(b["n_original_source_heldout_site_folds"],12)
        self.assertGreater(a["M2_minus_M3_rank_added_after_taxon_logloss_gain"],.03)
        self.assertGreater(b["M2_minus_M3_rank_added_after_taxon_logloss_gain"],.03)
        self.assertGreater(a["M2_minus_M3_rank_added_after_taxon_Brier_gain"],0)
        self.assertLess(a["models"][v.MODELS[3]]["mean_heldout_site_logloss_nats"],a["models"][v.MODELS[2]]["mean_heldout_site_logloss_nats"])
        self.assertTrue(out["leave_entire_original_spatial_site_out_27_fold"])

    def test_label_values_binary_and_fixed_source_species_universe(self):
        rows=fake_rows()
        rows[0]["taxon"]=124
        with self.assertRaises(ValueError):v.make_design(rows)
        rows=fake_rows()
        rows[0]["re_recorded_2017"]=2
        with self.assertRaises(ValueError):v.make_design(rows)

    def test_loss_function_math(self):
        import numpy as np
        loss,brier=v.row_losses(np.array([1,0]),np.array([.8,.2]))
        self.assertAlmostEqual(loss[0],-math.log(.8))
        self.assertAlmostEqual(loss[1],-math.log(.8))
        self.assertAlmostEqual(brier.mean(),.04)


if __name__=="__main__":
    unittest.main()

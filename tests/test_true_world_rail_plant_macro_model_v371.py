#!/usr/bin/env python3
"""v371: synthetic maths only; no artificial original global rail/plants claims."""
from __future__ import annotations
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from model_global_railkm_plant_landscape_macro_v371 import (
  GRID_SIZES,RIDGE_LAMBDA,CONTINENTS,MIN_CONTINENT_TEST_PAIRS,
  richness_features,beta_features,block_weights,ridge_estimate,
  cv_country_holdout,RICH_BASE,BETA_BASE,RICH_RAIL,BETA_RAIL
)
class V371GenuineGlobalMacroModelTests(unittest.TestCase):
    def src(self,i,continent="Asia",railkm=12.0):
        return {"cell":(i,i+1),"continent":continent,"block":(i//10,0),
          "source_grid_width_m":100000,
          "source_original_bio1_C":14.0,"abs_lat":35.0,
          "n_photo_plant_observations":300+i%7,
          "n_public_distinct_plant_observers":35+i%6,
          "n_observed_plant_species":72+i%5,
          "rail_km":railkm}
    def test_source_original_two_truly_global_study_grid_grains(self):
        self.assertEqual(GRID_SIZES,("100000","200000"))
        self.assertEqual(RIDGE_LAMBDA,.10)
        self.assertEqual(len(CONTINENTS),6)
        self.assertEqual(MIN_CONTINENT_TEST_PAIRS,50)
    def test_beta_features_swapping_source_adjacent_cells_is_symmetric(self):
        a=self.src(0,railkm=10.)
        b=self.src(1,railkm=90.)
        p={"A":a,"B":b,"direction":0,"y":0.55,"continent":"Asia",
           "block":(0,0)}
        q={**p,"A":b,"B":a}
        self.assertEqual(beta_features(p),beta_features(q))
        self.assertTrue("log_mean_rail_km" in beta_features(p))
        self.assertTrue("abs_log_rail_contrast" in beta_features(p))
    def test_richness_baseline_excludes_direct_railkm(self):
        d=richness_features(self.src(0))
        self.assertEqual(set(RICH_BASE)|set(RICH_RAIL),set(d))
        self.assertFalse(set(RICH_BASE)&set(RICH_RAIL))
        self.assertFalse(set(BETA_BASE)&set(BETA_RAIL))
    def test_true_spatial_1000km_block_equal_weights(self):
        a=self.src(0);b=self.src(1);c=self.src(20)
        weights=block_weights([a,b,c])
        self.assertAlmostEqual(float(weights.sum()),1.0)
        self.assertAlmostEqual(float(weights[0]+weights[1]),float(weights[2]))
    def test_genuine_global_source_cell_uses_original_1000km_block_field(self):
        # Original source geo records did not expose the shorter 'block'
        # alias that synthetic units had inadvertently always supplied.
        a=self.src(0)
        a["original_equal_area_1000km_block"]=a.pop("block")
        b=self.src(1)
        b["original_equal_area_1000km_block"]=b.pop("block")
        c=self.src(20)
        c["original_equal_area_1000km_block"]=c.pop("block")
        weights=block_weights([a,b,c])
        self.assertAlmostEqual(float(weights.sum()),1.0)
        self.assertAlmostEqual(float(weights[0]+weights[1]),float(weights[2]))
    def test_global_continent_heldout_does_not_train_on_the_tested_continent(self):
        allz=[]
        for c in CONTINENTS:
            for i in range(120):
                z=self.src(i,continent=c,railkm=float(10+i%18))
                z["block"]=(i//12,list(CONTINENTS).index(c))
                allz.append(z)
        r=cv_country_holdout(allz,"richness")
        self.assertTrue(r["identified"])
        self.assertEqual(len(r["source_per_continent_heldout_reproducible_losses"]),6)
        self.assertTrue(all(f["n_training_original_macro_source_units"]==600 for f in r[
             "source_per_continent_heldout_reproducible_losses"]))
    def test_exact_railkm_incremental_math_using_synthetic_only(self):
        src=[]
        for i in range(160):
            x=self.src(i,railkm=3.+i*1.5)
            x["n_observed_plant_species"]=max(21,int(15+i*.6))
            src.append(x)
        train=src[:100];test=src[100:]
        target=lambda z: math.log1p(z["n_observed_plant_species"])
        p0=ridge_estimate(train,test,RICH_BASE,(),target,"richness")
        p1=ridge_estimate(train,test,RICH_BASE,RICH_RAIL,target,"richness")
        self.assertEqual(len(p0),60)
        self.assertTrue(math.isfinite(float(p1[0])))
if __name__=="__main__":unittest.main()

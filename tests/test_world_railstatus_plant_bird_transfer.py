import unittest
from unittest.mock import patch
import numpy as np
import analyze_world_railstatus_plant_bird_transfer as t

class TestTransfer(unittest.TestCase):
    def test_disjoint_holdout_and_increment(self):
        rows=[{"continent":c,"block":i//5,"donor":float(i)} for c in ("A","B","C","D") for i in range(110)]
        seen=[]
        def fit(train,test,fields,target,clip):
            self.assertTrue({r["continent"] for r in train}.isdisjoint({r["continent"] for r in test}))
            self.assertNotIn(target,fields)
            recipient,j,_=target.split("_")
            donor="plants" if recipient=="birds" else "birds"
            self.assertNotIn(f"{donor}_{j}_species",fields)
            seen.append((tuple(fields),target))
            return np.full(len(test),1. if f"{donor}_{1-int(j)}_species" in fields else 2.)
        with patch.object(t,"CONTINENTS",("A","B","C","D")), patch.object(t.b,"local_support",side_effect=lambda x:x), patch.object(t.b,"with_epoch",side_effect=lambda x,y:x), patch.object(t,"fit_predict",side_effect=fit), patch.object(t,"metrics",side_effect=lambda rows,pred,target:(float(np.mean(pred**2)),0.)):
            z=t.transfer(rows,flex=False)
            self.assertEqual(z["n_folds"],4)
            self.assertEqual(z["mean_delta_MSE"]["richness_given_biomass"],3.)
            self.assertEqual(z["n_positive_folds"]["richness_given_biomass"],4)
            self.assertEqual(len(seen),32)
    def test_unsupported(self):
        with patch.object(t.b,"local_support",side_effect=lambda x:x), patch.object(t.b,"with_epoch",side_effect=lambda x,y:x):
            z=t.transfer([])
            self.assertFalse(z["adequate_folds"])
            self.assertIsNone(z["relative_richness_gain_given_biomass"])
if __name__=="__main__":unittest.main()

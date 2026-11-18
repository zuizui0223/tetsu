"""v26.3 synthetic QA for species×habitat and training-only centered rank."""
import sys
import unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import validate_Kashmir_taxon_zone_centered_abundance_v263 as v


def synthetic_sites():
    rows=[]
    for i in range(27):
        for z,zone in enumerate(v.ZONES):
            for t in range(5):
                position=(t+2*i+z)%5
                rank=position/2.-1.
                rows.append({
                    "site":f"s{i:02d}","group":"station" if i<15 else "between",
                    "zone":zone,"taxon":t,
                    "x2014_rank":rank,"x2014_raw":position+1.,
                    "re_recorded_2017":int(position>=3)
                })
    return rows


class SpeciesByZoneValidation(unittest.TestCase):
    def test_training_only_taxon_zone_centering(self):
        keys=[(1,0),(1,0),(1,1),(2,0)]
        ranks=np.array([-1.,1.,.2,.8])
        tr,te,support=v.train_only_centered_rank(
            np.array([0,1,2]),np.array([3]),keys,ranks)
        self.assertEqual(tr.tolist(),[-1.,1.,0.])
        self.assertEqual(te.tolist(),[0.])
        self.assertEqual(support["n_test_events_without_train_same_taxon_zone"],1)

    def test_centering_source_zone_not_other_zone(self):
        keys=[(7,0),(7,1),(7,0)]
        ranks=np.array([-.3,.8,.5])
        tr,te,support=v.train_only_centered_rank(
            np.array([0,1]),np.array([2]),keys,ranks)
        self.assertAlmostEqual(te[0],.8)
        self.assertEqual(support["n_test_events_with_train_same_taxon_zone"],1)

    def test_species_by_zone_design_shape(self):
        rows=synthetic_sites()
        y,base,rank,meta=v.build_species_zone_design(rows)
        self.assertEqual(base.shape,(27*4*5,8+4*123))
        self.assertEqual(len(meta["by_site"]),27)
        self.assertEqual(len(meta["keys"]),len(rows))

    def test_rank_predictive_after_flexible_species_by_zone(self):
        rows=synthetic_sites()
        result=v.audit(rows,1.)
        a=result["source_group_zone_habitat_results"]["station"]["Track"]
        b=result["source_group_zone_habitat_results"]["between"]["Track"]
        self.assertEqual(a["source_site_heldout_folds"],15)
        self.assertEqual(b["source_site_heldout_folds"],12)
        self.assertGreater(a["M2Z_minus_M3Z_logloss_gain"],.02)
        self.assertGreater(b["M2Z_minus_M3Z_logloss_gain"],.02)
        self.assertGreater(a["M2Z_minus_M4Z_centered_logloss_gain"],.02)
        self.assertTrue(result["training_only_rank_centering_never_leaks_test_site"])


if __name__=="__main__":
    unittest.main()

import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
from aggregate_global_tree_consensus_EarthEnv_v380a import rollup,source_rows,NX,NY,W
class TestTree(unittest.TestCase):
    def test_2x2_conservation_and_no_missing_imputation(self):
        area=np.zeros((NY,NX));tree=area.copy()
        area[0,0]=2e9;tree[0,0]=1e9
        area[0,1]=4e9;tree[0,1]=1e9
        a=rollup(area);t=rollup(tree)
        self.assertEqual(a.sum(),area.sum())
        self.assertEqual(t.sum(),tree.sum())
        rows=source_rows(t,a,2*W)
        self.assertEqual(len(rows),1)
        self.assertAlmostEqual(rows[0]["source_tree_consensus_mean_fraction_of_valid_source_area"],1/3)
        self.assertEqual(rows[0]["source_equal_area_grid_ix"],0)
    def test_bad_source_area_rejected(self):
        with self.assertRaises(ValueError):source_rows(np.array([[2.]]),np.array([[1.]]),W)
if __name__=="__main__":unittest.main()

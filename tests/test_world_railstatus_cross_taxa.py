import unittest
import random
from collections import defaultdict
from scan_world_railstatus_cross_taxa import bucket,add,serialize
from analyze_world_railstatus_cross_taxa import paired_interval,permutation_indices,restrict,rel

class CrossTaxaTests(unittest.TestCase):
    def test_counts_not_taxon_or_observer_ids(self):
        groups={"plants":defaultdict(bucket)}
        add(groups,"plants",(1,2),10,5)
        add(groups,"plants",(1,2),10,5)
        add(groups,"plants",(1,2),11,6)
        z=serialize({100000:groups})["100000"]["plants"][0]
        self.assertEqual((z["observations"],z["species"],z["observers"]),(3,2,2))
        self.assertEqual(len(z),5)
    def test_shared_band_shuffle(self):
        rows=[dict(continent="A",bio1=10,abs_lat=20,key=(i,0)) for i in range(10)]
        rows.append(dict(continent="B",bio1=10,abs_lat=20,key=(10,0)))
        perm=permutation_indices(rows,random.Random(1))
        self.assertEqual(perm[-1],10)
        self.assertEqual(sorted(perm[:-1]),list(range(10)))
    def test_cluster_difference_cancels_identical_scores(self):
        z={("A",1):.1,("A",2):-.1}
        self.assertEqual(paired_interval(.2,z,.1,z)["block_cluster_SE"],0)
        self.assertIsNone(paired_interval(.2,z,None,{}))
    def test_support_and_relative_gain(self):
        self.assertEqual(len(restrict([{"key":(1,2)},{"key":(3,4)}],{(3,4)})),1)
        self.assertAlmostEqual(rel({"mean_dMSE":.01,"base_mean_MSE":.2}),.05)
if __name__=="__main__":unittest.main()

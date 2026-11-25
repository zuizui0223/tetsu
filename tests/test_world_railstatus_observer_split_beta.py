import unittest
from scan_world_railstatus_observer_split_beta import bucket,add,observer_side,full_counts,split_counts,pair_rows,overlap,validate_counts
class TestObserverSplit(unittest.TestCase):
    def test_global_split_and_missing(self):
        a=bucket();b=bucket()
        for x in (a,b):
            add(x,3,17);add(x,4,None)
        self.assertEqual(full_counts(a),(2,1,2))
        side=observer_side(17)
        self.assertEqual(split_counts(a)[side],(1,1,1))
        self.assertEqual(split_counts(a)[1-side],(0,0,0))
        self.assertEqual(a[4][side],b[4][side])
    def test_pairs(self):
        a=bucket();b=bucket()
        for side in (0,1):
            ids=[i for i in range(100) if observer_side(i)==side][:5]
            for i in range(25):
                add(a,i,ids[i%5]);add(b,i+2,ids[i%5])
        z=pair_rows({(0,0):a,(1,0):b})[0]
        self.assertEqual(z["species_full"]["shared"],23)
        self.assertEqual(z["observer_full"]["jaccard"],1)
        self.assertEqual(z["split_similarity"]["01"]["shared"],23)
        self.assertEqual(z["split_similarity"]["10"]["shared"],23)
        self.assertEqual(len(z["split_similarity"]),4)
    def test_small_support_no_split(self):
        a=bucket();b=bucket()
        for k in range(60):add(a,k%25,k%5);add(b,k%25,k%5)
        z=pair_rows({(0,0):a,(0,1):b})[0]
        self.assertNotIn("split_similarity",z)
if __name__=="__main__":unittest.main()

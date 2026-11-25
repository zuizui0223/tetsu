import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from scan_world_railstatus_bird_hwi_conditional_v416 import rank_groups
class TestConditionalHWI(unittest.TestCase):
    def test_equal_labels_in_each_stratum(self):
        r={i:(1,1,1,float(i)) for i in range(6)}
        r.update({10+i:(2,1,3,float(i+10)) for i in range(5)})
        lo,hi,q=rank_groups(r,True)
        self.assertEqual(lo,{0,1,2,10,11})
        self.assertEqual(hi,{3,4,5,13,14})
        self.assertEqual(q["paired_species_per_label"],5)
    def test_ties_and_incomplete(self):
        r={i:(1,1,1,10.) for i in range(6)}
        r.update({10+i:(2,None,1,float(i)) for i in range(5)})
        lo,hi,q=rank_groups(r,True)
        self.assertEqual(lo,set());self.assertEqual(hi,set())
        self.assertEqual(q["skipped_equal_HWI_pairs"],3)
if __name__=="__main__":unittest.main()

import unittest
from analyze_world_railstatus_bird_avonet_temporal_audit import EXPECTED,FIELDS,END,audits
class TestAudit(unittest.TestCase):
    def test_design_and_gates(self):
        self.assertEqual(EXPECTED[.5,"HWI"][0],537)
        self.assertEqual(len(FIELDS),len(set(FIELDS)))
        self.assertEqual(len(END),len(set(END)))
        z=audits({"2023":[],"2024":[],"2025":[]})
        self.assertNotEqual(z["block_FE"]["difference_jaccard"]["status"],"estimated")
if __name__=="__main__":unittest.main()

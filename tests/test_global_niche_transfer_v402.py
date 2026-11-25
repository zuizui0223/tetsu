import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_niche_transfer_v402 import block,fold,partitions,separated,MODELS,TARGET,COLD
class TestV402(unittest.TestCase):
    def test_same_block_and_buffer(self):
        r=[dict(continent="Asia",block=(i,j)) for i in range(15) for j in range(12)]
        for f,tr,te in partitions(r,"buffered_block"):
            self.assertTrue(all(separated(block(a),block(b)) for a in tr for b in te))
            self.assertTrue(all(fold(block(x))==int(f) for x in te))
    def test_nested_models(self):
        self.assertNotIn(TARGET,MODELS["environment"])
        self.assertNotIn(COLD,MODELS["rail_main"])
        self.assertEqual(MODELS["rail_main_and_cold"][-2:],(TARGET,COLD))
if __name__=="__main__":unittest.main()

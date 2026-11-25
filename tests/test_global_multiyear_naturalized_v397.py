import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_multiyear_naturalized_v397 import pooled,fraction_rows,support,fixed_temporal
class TestV397(unittest.TestCase):
    def test_pool_counts_species_observers_as_unions(self):
        d={"2023":{(0,0):[3,{1,2},{10,11}]},"2024":{(0,0):[2,{2,3},{11,12}]},"2025":{(0,0):[1,{3},{12}]}}
        self.assertEqual(pooled(d)[(0,0)],[6,{1,2,3},{10,11,12}])
    def test_documented_is_not_unlisted(self):
        rows=[{"key":(0,0)}]
        ids=set(range(25));names={i:str(i) for i in range(20)}
        x=fraction_rows(rows,{(0,0):[30,ids,set()]},names,{(0,0):"A"},{"A":{"1","2","23"}})[0]
        self.assertEqual(x["vascular_richness"],20);self.assertEqual(x["documented_richness"],2)
        self.assertAlmostEqual(x["known_fraction"],.1)
    def test_temporal_support_no_imputation(self):
        r={"key":(0,0),"continent":"Asia","block":(0,0),"observations":100,"observers":10}
        out=fixed_temporal([((0,0),(1,0))],{"2023":[],"2024":[],"2025":[r]},False)
        self.assertEqual(out[0]["both_years_supported_pairs"],0)
        self.assertEqual(support([r],"2023-2025","annualized_support"),[])
if __name__=="__main__":unittest.main()

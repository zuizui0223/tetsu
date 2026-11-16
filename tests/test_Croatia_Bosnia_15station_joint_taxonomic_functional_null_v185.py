import importlib.util,random,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/analyze_Croatia_Bosnia_15station_joint_taxonomic_functional_null_v185.py"
spec=importlib.util.spec_from_file_location("v185",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class TestDualEcologicalResponse(unittest.TestCase):
    def toy(self):
        # Entirely artificial and NEVER a real source estimate.
        return [{(i*7+j)%32 for j in range(16)} for i in range(15)],(["T","H","P","G","Ch"]*7)[:32]
    def test_observational_ranges_and_pair_n(self):
        ss,life=self.toy();out=m.pairwise_metrics(ss,life)
        self.assertEqual(out["number_within_climate_pairs"]+out["number_between_climate_pairs"],105)
        self.assertTrue(0<=out["tax_jaccard_mean"]<=1)
        self.assertTrue(0<=out["lifeform_overlap_mean"]<=1)
    def test_fixed_margin_and_climate_nulls(self):
        ss,life=self.toy()
        for climate in (False,True):
            a=m.generate(ss,life,climate,chains=2,draws=3,burnin=150,gap=40)
            self.assertEqual(a["draws_total_not_independent_field_replicates"],6)
            self.assertTrue(a["source_taxon_occurrence_and_site_richness_FIXED"])
            self.assertTrue(a["this_is_source_internal_not_independent_nonrail_regional_species_pool"])
            self.assertEqual(set(m.STATS),set(a["statistics"]))
            for item in a["statistics"].values():
                self.assertTrue(0<=item["observed"]<=1)
                self.assertTrue(0<=item["internal_source_null_mean"]<=1)
    def test_strict_source_dimensions(self):
        ss,life=self.toy()
        with self.assertRaises(ValueError):
            m.pairwise_metrics(ss[:14],life)
if __name__=="__main__":
    unittest.main()

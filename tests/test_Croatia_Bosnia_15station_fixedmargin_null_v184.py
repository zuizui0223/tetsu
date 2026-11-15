import importlib.util
from pathlib import Path
import unittest

p=Path(__file__).resolve().parents[1]/"scripts/analyze_Croatia_Bosnia_15station_lifeform_fixedmargin_null_v184.py"
spec=importlib.util.spec_from_file_location("v184",str(p))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
import random

class TestSourceConstrainedNull(unittest.TestCase):
    def sample(self):
        # Synthetic ONLY: 15 sites and 30 taxa, never counted as real results.
        return [{(i*3+j)%30 for j in range(12)} for i in range(15)]
    def test_unrestricted_margin_conservation(self):
        ss=self.sample()
        before=m.margins(ss,30)
        accepted=m.swaps(ss,random.Random(7),2000,False)
        self.assertGreater(accepted,0)
        self.assertEqual(before,m.margins(ss,30))
    def test_climate_stratified_conserves_within_climate_taxon_occurrences(self):
        ss=self.sample()
        before=m.margins(ss,30)
        climate=m.climate_margins(ss,30)
        accepted=m.swaps(ss,random.Random(7),2000,True)
        self.assertGreater(accepted,0)
        self.assertEqual(before,m.margins(ss,30))
        self.assertEqual(climate,m.climate_margins(ss,30))
    def test_null_tail_is_diagnostic_not_inference_on_independent_sites(self):
        ss=self.sample()
        life=["T","H","P","G","Ch"]*6
        obs=m.profile_overlap(ss,life)
        self.assertTrue(0<=obs[0]<=1)
        out=m.null_model(ss,life,False,chains=2,draws=5,burnin=100,gap=100)
        self.assertEqual(out["n_synthetic_replicates"],10)
        self.assertTrue(0<=out["one_sided_tail_ge_observed"]<=1)
        self.assertTrue(out["null_is_internal_rail_species_pool_NOT_natural_habitat_pool"])
    def test_exact_publisher_macroclimate_metadata(self):
        self.assertEqual(m.STATION_MACROCLIMATE[10],"MED")
        self.assertEqual(m.STATION_MACROCLIMATE[15],"TWS")
        self.assertEqual(sum(x=="TSS" for x in m.STATION_MACROCLIMATE.values()),5)

if __name__=="__main__": unittest.main()

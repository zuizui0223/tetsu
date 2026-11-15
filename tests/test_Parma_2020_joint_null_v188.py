import importlib.util
from pathlib import Path
import unittest
p=Path(__file__).resolve().parents[1]/"scripts/analyze_Parma_2020_joint_taxonomic_lifeform_null_v188.py"
spec=importlib.util.spec_from_file_location("v188",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestOriginalParmaNull(unittest.TestCase):
    def toy(self):
        n=28;sp=48
        sets=[{(i*7+j)%sp for j in range(18)} for i in range(n)]
        groups=["plain" if i<11 else "other" for i in range(28)]
        life=list(m.FORM6)*8
        return sets,life,groups
    def test_frozen_original_stations(self):
        self.assertEqual(len(m.PLAIN_11),11)
        self.assertEqual(len(set(m.STATION_CODES)-m.PLAIN_11),17)
        self.assertEqual(len(m.EXTRA_EFFORT_STATIONS),5)
    def test_original_site_noncrossing_group_swaps(self):
        sets,life,g=self.toy()
        before=m.margin(sets,len(life))
        group=m.group_margin(sets,len(life),g)
        num=m.swaps(sets,__import__("random").Random(41),2000,g,True)
        self.assertGreater(num,0)
        self.assertEqual(before,m.margin(sets,len(life)))
        self.assertEqual(group,m.group_margin(sets,len(life),g))
    def test_mean_tax_and_lifeforms_same_28(self):
        sets,life,g=self.toy()
        r=m.metric(sets,life,g)
        self.assertEqual(r["n_site_pairs"],378)
        self.assertEqual(r["n_within_group_pairs"]+r["n_between_group_pairs"],378)
        self.assertTrue(0<=r["tax_mean"]<=1)
        self.assertTrue(0<=r["functional_mean"]<=1)
    def test_simulation_not_independent_field_replicates(self):
        sets,life,g=self.toy()
        r=m.evaluate(sets,life,g,True,chains=2,draws=3,burnin=200,gap=60)
        self.assertTrue(r["not_a_regional_non_rail_flora_pool"])
        self.assertTrue(r["source_group_occurrences_fixed_if_stratified"])
        self.assertEqual(r["n_correlated_technical_draws"],6)
        self.assertEqual(len(r["stats"]),6)
if __name__=="__main__":unittest.main()

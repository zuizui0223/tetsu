"""v15.8 synthetic model arithmetic tests; no false actual Czech field plants."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_Czech_2024_centroid_geographic_identity_v154 import fake_source_sites
from analyze_Czech_arch_neo_5stretch_loo_contrast_v158 import paired_group_stats
from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    real_site_pairs, real_site_pairs_without_n_gate
)
from analyze_Czech_2024_real_archaeophyte_neophyte_geographic_identity_v157 import group_sets

class SameSourcePairWindowClassComparison(unittest.TestCase):
    def setUp(self):
        sites,guild=fake_source_sites()
        stat={}
        for s in sites:
            for taxon in s["all_original_taxa"]:
                stat[taxon]=("arch" if taxon.startswith("w") else
                             "neo" if taxon.startswith("o") else None)
        self.sites,self.guild,self.stat=sites,guild,stat

    def test_same_source_geography_functional_support_for_two_groups(self):
        result=paired_group_stats(self.sites,self.guild,self.stat)
        self.assertTrue(result["identified"])
        self.assertEqual(result["n_original_source_stations"],35)
        self.assertEqual(result["n_dependent_source_station_pairs"],595)
        matched=result["matched_group_contrast"]["geography_and_function"]
        self.assertTrue(matched["identified"])
        self.assertTrue(matched["identical_original_source_pair_and_distance_guild_strata"])
        self.assertAlmostEqual(matched["arch_minus_neo_source_gap_contrast"],
            matched["archaeophyte_original_matched_source_gap"]-
            matched["neophyte_original_matched_source_gap"],places=12)

    def test_source_class_prevalence_differences_are_retained_not_equalized(self):
        result=paired_group_stats(self.sites,self.guild,self.stat)
        a=result["subgroup"]["archaeophyte"]
        n=result["subgroup"]["neophyte"]
        self.assertGreater(a["source_original_taxa"],0)
        self.assertGreater(n["source_original_taxa"],0)
        self.assertTrue(result["not_a_paired_taxon_experiment"])
        self.assertEqual(a["original_station_taxon_incidences"],70)
        self.assertEqual(n["original_station_taxon_incidences"],70)

    def test_subgroup_absence_disables_false_jaccard_double_zero(self):
        sites=[dict(s) for s in self.sites]
        sites[0]["all_original_taxa"]=set(s for s in sites[0]["all_original_taxa"] if s.startswith("w"))
        r=paired_group_stats(sites,self.guild,self.stat)
        self.assertFalse(r["identified"])
        self.assertIn("lacks",r["reason"])

    def test_original_leave_one_stretch_sites_rebuilt_not_reused_pair_ids(self):
        smaller=[s for s in self.sites if s["section"]!="0"]
        self.assertEqual(len(smaller),28)
        p=real_site_pairs_without_n_gate(smaller,self.guild)
        self.assertEqual(len(p),28*27//2)
        result=paired_group_stats(smaller,self.guild,self.stat)
        self.assertTrue(result["identified"])
        self.assertEqual(result["n_original_source_stretches"],4)
        self.assertEqual(result["n_dependent_source_station_pairs"],378)

if __name__=="__main__":unittest.main()

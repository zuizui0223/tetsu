"""v14.1 genuine-source identity/functional comparison math unit tests.

Synthetic fixtures below ONLY validate the algorithm. The actual field
biology is computed solely from the SHA-audited original 2026 89-site
Slovenian study XLSX by CI, NEVER synthetic site metadata.
"""
from __future__ import annotations
import copy,random,sys,unittest
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_functional_label,jaccard,functional_overlap,make_matrices,
    row_col_margins,swap_proposals,matrix_taxon_lists,taxonomic_gap,
    pair_summary,null_diagnostics,section_leaveout
)

GUILD={"A":"strict_annual_wind","B":"strict_annual_wind",
       "C":"strict_annual_wind","D":"strict_annual_wind",
       "E":"other_known_other","F":"other_known_other"}
def site(name,section,taxa,year=2022,typ="single"):
    return {
      "site":name,"section":section,"year":year,"track_type":typ,
      "known":set(taxa),"all_original_taxa":set(taxa),
      "n_original_matched_track_sides":1,
    }
SITES=[
    site("alpha1","alpha",("A","B","E"),2022),
    site("alpha2","alpha",("A","B","E"),2023),
    site("beta1","beta",("C","D","F"),2022),
    site("beta2","beta",("C","D","F"),2023),
    site("gamma1","gamma",("A","C","E"),2022),
    site("gamma2","gamma",("B","D","F"),2023),
]


class OriginalFieldSourceFunctionalIdentityMath(unittest.TestCase):
    def test_original_strict_strategy_labels_not_taxonomic_genera(self):
        self.assertEqual(source_functional_label({
            "Diaspore.strategy":"Epilobium",
            "Lifespan.categories":"annual"}),"strict_annual_wind")
        self.assertEqual(source_functional_label({
            "Diaspore.strategy":"Phragmites",
            "Lifespan.categories":"perennial"}),"other_known_wind")
        self.assertEqual(source_functional_label({
            "Diaspore.strategy":"Allium",
            "Lifespan.categories":"annual"}),"strict_annual_other")
        self.assertIsNone(source_functional_label({
            "Diaspore.strategy":None,
            "Lifespan.categories":"annual"}))
        self.assertIsNone(source_functional_label({
            "Diaspore.strategy":"Epilobium",
            "Lifespan.categories":None}))

    def test_coarse_functional_overlap_identical_even_when_species_are_completely_different(self):
        self.assertEqual(jaccard({"A","B"},{"C","D"}),0.)
        self.assertAlmostEqual(functional_overlap({"A","B"},{"C","D"},GUILD),1.)
        self.assertAlmostEqual(functional_overlap({"A","E"},{"C","F"},GUILD),1.)
        self.assertIsNone(functional_overlap(set(),{"A"},GUILD))

    def test_original_source_site_and_species_guild_margins_exactly_preserved_by_null(self):
        matrix=make_matrices(SITES,GUILD)
        start=row_col_margins(matrix)
        n=swap_proposals(matrix,SITES,random.Random(12),6000,False)
        self.assertGreater(n,0)
        self.assertEqual(start,row_col_margins(matrix))
        original=[s["known"] for s in SITES]
        later=matrix_taxon_lists(matrix,len(SITES))
        self.assertEqual([len(a) for a in later],[len(x) for x in original])

    def test_null_stratifies_source_species_by_year_and_track_type_not_just_global(self):
        matrix=make_matrices(SITES,GUILD)
        strata=defaultdict(lambda:Counter())
        for i,s in enumerate(SITES):
            strata[(s["year"],s["track_type"])].update(s["known"])
        n=swap_proposals(matrix,SITES,random.Random(44),4000,True)
        self.assertGreater(n,0)
        after=defaultdict(lambda:Counter())
        for i,s in enumerate(SITES):
            after[(s["year"],s["track_type"])].update(matrix_taxon_lists(matrix,len(SITES))[i])
        self.assertEqual(dict(strata),dict(after))

    def test_functional_and_taxonomic_pair_metrics_are_separate_descriptive_estimands(self):
        summary=pair_summary(SITES,GUILD)
        self.assertEqual(summary["n_source_sites"],6)
        self.assertEqual(summary["n_dependent_within_section_site_pairs_NOT_independent_studies"],3)
        self.assertEqual(summary["n_dependent_between_section_site_pairs_NOT_independent_studies"],12)
        self.assertGreater(summary["taxonomic_Jaccard_similarity_dual_trait_known_same_taxa"]["same_minus_different"],0)
        self.assertGreaterEqual(summary["functional_4guild_proportional_overlap_on_dual_known_taxa"]["mean_same_rail_section"],0)
        self.assertEqual(summary["n_dependent_between_section_site_pairs_NOT_independent_studies"]+3,15)

    def test_two_source_null_chains_and_15_section_leaveout_math(self):
        d=null_diagnostics(SITES,GUILD,"all_years_both_track_types",
                           draws_per_chain=25,nchains=2,burnin=200,spacing=80,seed=15)
        self.assertEqual(d["MonteCarlo_draws_NOT_independent_field_studies"],50)
        self.assertEqual(len(d["separate_chains_and_within_chain_autocorrelation_not_exact_permutation"]),2)
        self.assertTrue(0<d["source_guild_and_species_occupancy_fixed_margin_upper_tail_fraction_plus_one"]<=1)
        c=section_leaveout(SITES,GUILD)
        self.assertEqual(c["original_source_rail_sections_left_out"],3)

    def test_missing_original_traits_are_absent_not_converted_to_other(self):
        source=copy.deepcopy(SITES)
        source[0]["all_original_taxa"].add("MISSING")
        metrics=pair_summary(source,GUILD)
        self.assertNotEqual(metrics["taxonomic_Jaccard_similarity_dual_trait_known_same_taxa"],
                            metrics["taxonomic_Jaccard_similarity_ALL_312_original_source_taxa_sensitivity"])


if __name__=="__main__":unittest.main()

"""v28.1 synthetic source matched-triad frequency/richness fixed null QA."""
from __future__ import annotations
import random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Slovenia_same_line_geodistance_matched_frequency_null_v281 as m
import audit_Slovenia_source_same_line_vs_other_geodistance_matched_v280 as design


def fake_original():
    sites=[]
    for line in range(5):
        for site in range(3):
            sites.append({
                "source_id":f"line{line}_site{site}",
                "section":str(line),
                "class":"main_G" if line%2 else "regional_R",
                "lat":46.+line*.008,"lon":14.+site*.033,
                "year":2022 if line%2 else 2023,
                "layout":"SINGLE" if line<3 else "DOUBLE",
                "BS":{1,line+10,site+30},
                "BT":{1,line+10,site+40}
            })
    return sites


class GenuineLineMatchedNullTests(unittest.TestCase):
    def setUp(self):
        self.old=(design.N_SITES,design.N_LINES,design.N_SAME_UNORDERED)
        design.N_SITES=15
        design.N_LINES=5
        design.N_SAME_UNORDERED=15

    def tearDown(self):
        design.N_SITES,design.N_LINES,design.N_SAME_UNORDERED=self.old

    def test_true_geodistance_matched_original_site_gap_reproducible(self):
        s=fake_original()
        match,qa=design.match_candidates(s,20,2.5,"geo_only")
        self.assertGreaterEqual(len(match),12)
        rows=[design.source_community(x,"BS_and_BT") for x in s]
        obs,byline=m.original_matched_gap_for_rows(rows,s,match)
        expected=design.line_equal_aggregate(match,qa)["habitat_comparisons"]["BS_and_BT"][
            "equal_source_line_same_minus_different_line_similarity"]
        self.assertAlmostEqual(obs,expected)
        self.assertEqual(len(byline),5)

    def test_exact_source_site_richness_and_species_frequency_preserved(self):
        s=fake_original()
        rows=[set(design.source_community(x,"BS_and_BT")) for x in s]
        groups=m.source_groups(s,"within_original_survey_year_and_track_layout")
        baseline=m.exact_source_margins(rows,groups)
        updated=m.swap_within_source_groups(rows,groups,random.Random(77),300)
        self.assertGreaterEqual(updated,0)
        self.assertEqual(m.exact_source_margins(rows,groups),baseline)
        self.assertEqual(sum(len(g) for g in groups),len(s))

    def test_two_chains_and_all_site_taxon_frequency_preserving_null(self):
        s=fake_original()
        match,qa=design.match_candidates(s,20,2.5,"geo_only")
        res=m.experiment(s,match,m.source_groups(s,"all_original_sites"),
                         seed=122,draws=59,burnin=200,thin=25)
        self.assertEqual(res["n_source_correlated_MCMC_draws"],118)
        self.assertEqual(len(res["two_independent_seed_chain_diagnostics"]),2)
        self.assertEqual(res["original_lines_with_distance_matched_anchors"],5)
        self.assertTrue(res["site_richness_and_taxon_frequency_preserved_per_defined_group"])
        self.assertEqual(res,m.experiment(s,match,m.source_groups(s,"all_original_sites"),
                         seed=122,draws=59,burnin=200,thin=25))

    def test_source_year_track_stratification_changes_possible_null_but_not_support(self):
        s=fake_original()
        matched,qa=design.match_candidates(s,20,2.5,"geo_only")
        res=m.experiment(s,matched,m.source_groups(s,"within_original_survey_year_and_track_layout"),
                         seed=123,draws=59,burnin=200,thin=20)
        self.assertEqual(res["source_matched_anchor_site_count"],len(matched))
        self.assertEqual(res["original_lines_with_distance_matched_anchors"],5)
        self.assertTrue(all(c["accepted_taxon_checkerboard_swaps"]>=0 for c in res["two_independent_seed_chain_diagnostics"]))
        self.assertLessEqual(abs(res["source_observed_minus_fixed_frequency_null_line_gap"]),1.)

    def test_no_unsupported_reference_null_or_confound_strata(self):
        source=fake_original()
        with self.assertRaises(ValueError):
            m.source_groups(source,"choose_best_strata_after_source_fits")
        match,_=design.match_candidates(source,20,2.5,"geo_only")
        with self.assertRaises(ValueError):
            m.experiment(source,match,m.source_groups(source,"all_original_sites"),
                         seed=1,draws=20)


if __name__=="__main__":
    unittest.main()

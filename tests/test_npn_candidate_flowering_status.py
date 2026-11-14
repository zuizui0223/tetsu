"""Synthetic matching/selection regression tests: NOT an ecological effect."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_npn_candidate_flowering_status import (
    select_candidate, match_far, site_phenometric_species,
    midpoint_doy_of_observed_onset_bracket
)


class CandidateSeriesTests(unittest.TestCase):
    def test_observed_no_yes_bracket_midpoint_not_onset(self):
        summary={
            "first_yes_prior_no_date":"2024-04-10",
            "first_yes_observed_date":"2024-04-20"
        }
        self.assertAlmostEqual(
            midpoint_doy_of_observed_onset_bracket(summary),106.0
        )
        self.assertIsNone(midpoint_doy_of_observed_onset_bracket({
            "first_yes_prior_no_date":None,
            "first_yes_observed_date":"2024-04-20"
        }))
        self.assertIsNone(midpoint_doy_of_observed_onset_bracket({
            "first_yes_prior_no_date":"2024-01-01",
            "first_yes_observed_date":"2024-04-20"
        }))

    def test_candidate_selected_by_replication_only(self):
        classes={
            7:{"near":{1,2},"far":{3,4,5,6}},
            9:{"near":{1},"far":{4,5,6,7,8}},
            31:{"near":{1,2},"far":{3,4}}
        }
        self.assertEqual(select_candidate(classes),[(7,2,4)])

    def test_geographic_pairs_prevent_distant_unmatched_rail_phenology(self):
        frame={
            1:{"xy":(0.,0.)},
            2:{"xy":(100_000.,0.)},
            3:{"xy":(1000.,0.)},
            4:{"xy":(101_000.,0.)},
            5:{"xy":(300_000.,0.)}
        }
        pairs=match_far({1,2},{3,4,5},frame)
        self.assertEqual(len(pairs),2)
        self.assertEqual(set(x[1] for x in pairs),{3,4})
        self.assertTrue(all(x[2]<75 for x in pairs))

    def test_positive_site_is_not_coded_at_unobserved_registered_site(self):
        frame={1:{"bin":"near"},2:{"bin":"near"},3:{"bin":"far"}}
        records=[
            {"site_id":1,"species_id":7,"state":"MD","phenophase_description":"Open flowers"},
            {"site_id":3,"species_id":7,"state":"MD","phenophase_description":"Open flowers"}
        ]
        x,counts=site_phenometric_species(records,frame,"MD")
        self.assertEqual(x[7]["near"],{1})
        self.assertEqual(x[7]["far"],{3})
        self.assertEqual(select_candidate(x),[])


if __name__=="__main__":
    unittest.main()

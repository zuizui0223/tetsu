"""NPN phenophase status test fixtures are synthetic, NOT biological data."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.npn_individuals import (
    normalize_open_flowers,summarize_individuals,coverage_summary
)


def r(day,status, individual="A", phase="Open flowers", site="site1"):
    return {
        "site_id":site,"individual_id":individual,
        "observation_date":f"2024-05-{day:02d}",
        "phenophase_description":phase,"phenophase_status":status
    }


class NPNIndividualTests(unittest.TestCase):
    def test_yes_no_confidence_and_no_imputation(self):
        rows=[r(1,0),r(5,0),r(12,1),r(15,1),r(16,-1),r(17,0)]
        clean,excluded=normalize_open_flowers(rows)
        self.assertEqual(len(clean),5)
        self.assertEqual(excluded["uncertain_no_imputation"],1)
        s=summarize_individuals(clean)
        self.assertEqual(len(s),1)
        self.assertEqual(s[0]["n_yes"],2)
        self.assertEqual(s[0]["n_no"],3)
        self.assertEqual(s[0]["first_yes_prior_no_date"],"2024-05-05")
        self.assertEqual(s[0]["first_yes_bracket_days"],7)

    def test_conflicting_status_same_plant_day_excluded(self):
        rows=[r(1,0),r(1,1),r(5,1)]
        clean,excluded=normalize_open_flowers(rows)
        self.assertEqual(len(clean),1)
        self.assertEqual(excluded["conflicting_same_day"],2)

    def test_duplicate_agreement_not_pseudoreplication(self):
        clean,excluded=normalize_open_flowers([r(1,0),r(1,0),r(7,1)])
        self.assertEqual(len(clean),2)
        self.assertEqual(excluded["same_day_duplicate_rows_deduplicated"],1)

    def test_2_explicit_no_confirm_monitored_interruption_only(self):
        rows=[r(1,1),r(4,0),r(8,0),r(14,1)]
        clean,_=normalize_open_flowers(rows)
        self.assertTrue(summarize_individuals(clean)[0]["has_two_no_interrupted_flowering"])

    def test_long_unmonitored_gap_is_not_reflower(self):
        rows=[r(1,1),r(4,0),r(8,0),r(31,1)]
        clean,_=normalize_open_flowers(rows)
        self.assertFalse(summarize_individuals(clean)[0]["has_two_no_interrupted_flowering"])

    def test_one_absence_cannot_show_interrupted_flowering(self):
        rows=[r(1,1),r(7,0),r(14,1)]
        clean,_=normalize_open_flowers(rows)
        self.assertFalse(summarize_individuals(clean)[0]["has_two_no_interrupted_flowering"])

    def test_no_prior_observed_absence_means_no_onset_bracket(self):
        clean,_=normalize_open_flowers([r(2,1),r(4,1)])
        s=summarize_individuals(clean)
        self.assertIsNone(s[0]["first_yes_prior_no_date"])

    def test_other_phase_is_not_open_flower(self):
        clean,excluded=normalize_open_flowers([r(1,1,phase="Leaves"),r(2,0)])
        self.assertEqual(len(clean),1)
        self.assertEqual(excluded["other_phase_not_open_flowers"],1)


if __name__=="__main__":
    unittest.main()

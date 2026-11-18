"""v26.2 mathematical QA for pre-fixed C values and frozen source score gate."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Kashmir_abundance_identity_ridge_sensitivity_v262 as v


def make(C,gain_overrides=None):
    from copy import deepcopy
    if gain_overrides is None: gain_overrides={}
    items={}
    for group in ("station","between"):
        d=v.EXPECTED_MAIN_TRACK_M2_M3[group] if C==1.0 else (.013 if group=="station" else .019)
        d=gain_overrides.get(group,d)
        items[group]={"Track":{
            "n_original_source_heldout_site_folds":15 if group=="station" else 12,
            "source_2014_at_risk_taxon_site_events":116 if group=="station" else 54,
            "source_2017_re_recorded_original_taxon_site_events":22 if group=="station" else 19,
            "M2_minus_M3_rank_added_after_taxon_logloss_gain":d,
            "M2_minus_M4_rank_and_track_interaction_added_after_taxon_logloss_gain":d*.75,
            "M2_minus_M3_rank_added_after_taxon_Brier_gain":d*.3,
            "n_source_site_folds_M3_improves_M2":12 if group=="station" else 10,
        }}
    return {"source_predictions_by_group_and_habitat":items}


class RidgeSensitivityTests(unittest.TestCase):
    def test_predeclared_C_and_same_source_cohort(self):
        xs={C:make(C) for C in v.C_VALUES}
        out=v.compare(xs)
        self.assertTrue(out["station"]["sign_consistent_in_all_three_fixed_C"]["M2_minus_M3_rank_added_after_taxon_logloss_gain"])
        self.assertEqual(out["station"]["source_independent_sites"],15)
        self.assertEqual(set(out["between"]["rank_added_after_species_logloss_gain_by_fixed_C"]),{"0.2","1.0","5.0"})

    def test_poor_sensitivity_must_not_be_misreported_positive(self):
        xs={C:make(C,{"station":-.02} if C==5.0 else {}) for C in v.C_VALUES}
        out=v.compare(xs)
        self.assertFalse(out["station"]["sign_consistent_in_all_three_fixed_C"]["M2_minus_M3_rank_added_after_taxon_logloss_gain"])

    def test_source_frozen_exact_scores_not_silently_changed(self):
        xs={C:make(C) for C in v.C_VALUES}
        xs[1.0]=make(1.0,{"station":.019})
        with self.assertRaises(ValueError):
            v.compare(xs)
        with self.assertRaises(ValueError):
            v.compare({1.:make(1.)})


if __name__=="__main__":
    unittest.main()

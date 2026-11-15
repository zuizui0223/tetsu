"""Synthetic fixtures only: no fake plot numbers are presented as Czech ecology."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_CzechVeg_rail_context_all_candidate_SIGN_balance_v177 import (
    Q,NR,NC,NS,STATION,LLO,LHI,LOO_N,
    audit_panel, eligible_without_q_sign, sign, audit
)

def row(q, rail, road, pubs, sites, lo=None, hi=None):
    return {
        Q:q,NR:rail,NC:road,NS:pubs,STATION:sites,
        LLO:lo,LHI:hi,LOO_N:3
    }

class AuditSignUnfilteredSourceNotSyntheticEcologicalInference(unittest.TestCase):
    def test_sign_and_opportunity_gates_are_independent_of_q_sign(self):
        data = {
            "positive":row(.15,5,2,2,9,.02,.25),
            "negative":row(-.10,2,6,3,8,-.25,-.03),
            "zero":row(0.,4,4,2,20,-.04,.04),
            "one_study_positive":row(.60,11,0,1,39,.50,.70),
            "low_station_negative":row(-.35,2,8,3,2,-.50,-.20)
        }
        self.assertTrue(eligible_without_q_sign(data["negative"]))
        self.assertTrue(eligible_without_q_sign(data["zero"]))
        self.assertFalse(eligible_without_q_sign(data["one_study_positive"]))
        self.assertFalse(eligible_without_q_sign(data["low_station_negative"]))
        self.assertEqual(sign(None),"unestimated")
        self.assertEqual(sign(0.),"zero")
        self.assertEqual(sign(-0.4),"negative")

    def test_all_directions_kept_before_positive_ranking(self):
        data = {
            "positive":row(.15,5,2,2,9,.02,.25),
            "negative":row(-.10,2,6,3,8,-.25,-.03),
            "zero":row(0.,4,4,2,20,-.04,.04),
            "one_study_positive":row(.60,11,0,1,39,.50,.70)
        }
        recent={"positive":{Q:-.03},"negative":{Q:.07},"zero":{Q:None}}
        strict={"positive":{Q:.01},"negative":{Q:-.07},"zero":{Q:None}}
        rows, stats = audit_panel(data,recent,strict,{"positive"})
        self.assertEqual(len(rows),3)
        self.assertEqual(stats["n_q_positive_UNSELECTED_directional_denominator"],1)
        self.assertEqual(stats["n_q_negative_UNSELECTED_directional_denominator"],1)
        self.assertEqual(stats["n_q_exactly_zero"],1)
        self.assertEqual(stats["n_eligible_source_taxa_positive_loo_for_every_single_original_study_deletion"],1)
        self.assertEqual(stats["n_eligible_source_taxa_negative_loo_for_every_single_original_study_deletion"],1)
        self.assertEqual(stats["modern_2010plus_sparse_same_original_study_sensitivity"]["n_direction_reversed"],2)
        self.assertEqual(stats["STRICT_label_only_very_small_panel_sensitivity"]["n_same_nonzero_direction_as_original_full_panel"],2)
        self.assertEqual(stats["original_top20_species_with_negative_or_zero_2010plus_score"],["positive"])
        self.assertEqual({r["species"] for r in rows},{"negative","positive","zero"})

    def test_source_audit_rejects_fake_complete_or_incompatible_grain(self):
        with self.assertRaisesRegex(ValueError,"v17.4"):
            audit({"complete":True,"version":"synthetic_fake","n_station_source_taxa":761,"n_unique_publisher_stations":39})

if __name__ == "__main__":
    unittest.main()

"""Synthetic source/eligibility checks. No factual railway effect in fixtures."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_npn_active_rail_sites import (
    validate_station, reduce_station, summarize_cohort, freeze_sites
)


def observation(site, person, taxon, day, status):
    return {
        "site_id": site, "individual_id": person, "species_id": taxon,
        "observation_date": f"2024-04-{day:02d}",
        "phenophase_description": "Open flowers", "phenophase_status": status,
    }


def four_visits(site, taxon, person="p"):
    return [
        observation(site, person, taxon, 1, 0),
        observation(site, person, taxon, 5, 0),
        observation(site, person, taxon, 8, 1),
        observation(site, person, taxon, 12, 1),
    ]


class ActiveSiteAudit(unittest.TestCase):
    def test_station_schema_accepted(self):
        self.assertEqual(
            validate_station({
                "station_id": "42", "latitude": 39.0, "longitude": -76.3
            }), (42, -76.3, 39.0)
        )
        self.assertIsNone(validate_station({
            "station_id": "x", "latitude": 39, "longitude": -76.3
        }))

    def test_four_actual_repeat_dates_required(self):
        a = reduce_station(four_visits("10", 82))
        self.assertEqual(a["n_openflower_raw"], 4)
        self.assertTrue(a["species"][82]["site_taxon_eligible"])
        b = reduce_station(four_visits("10", 82)[:3])
        self.assertFalse(b["species"][82]["site_taxon_eligible"])

    def test_photo_unknown_is_not_nonflower(self):
        rows = [observation("10", "p", 82, d, -1) for d in (1, 5, 8, 12)]
        x = reduce_station(rows)
        self.assertFalse(x["species"][82]["site_taxon_eligible"])

    def test_same_taxon_shared_in_both_exposure_groups(self):
        pairs = [(10, 20, 2500.), (30, 40, 4100.)]
        sites = {
            10: reduce_station(four_visits("10", 82)),
            20: reduce_station(four_visits("20", 82)),
            30: reduce_station(four_visits("30", 98)),
            40: reduce_station(four_visits("40", 82)),
        }
        r = summarize_cohort(sites, pairs, "MD")
        self.assertEqual(r["n_shared_eligible_species_near_far"], 1)
        self.assertEqual(r["shared_eligible_species_without_outcome_direction"][0]["npn_species_id"], 82)
        self.assertEqual(r["shared_eligible_species_without_outcome_direction"][0]["n_near_eligible_sites"], 1)
        self.assertEqual(r["shared_eligible_species_without_outcome_direction"][0]["n_far_eligible_sites"], 2)

    def test_failed_station_not_counted_as_inactive(self):
        pairs = [(10, 20, 2500.)]
        r = summarize_cohort({10: reduce_station(four_visits("10", 82))}, pairs, "MD")
        self.assertEqual(r["n_complete_station_pairs_queried"], 0)
        self.assertEqual(r["n_shared_eligible_species_near_far"], 0)

    def test_missing_species_never_claims_no_shared_species(self):
        rows = [dict(observation("10","p",82,d,1),species_id=None) for d in (1,3)]
        with self.assertRaises(ValueError):
            reduce_station(rows)


if __name__ == "__main__":
    unittest.main()

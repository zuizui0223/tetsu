"""These tests verify sampling and event-date logic on fabricated visits."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.phenology import load_events, normalize_visits, weekly_summary


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


class PhenologyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.events = Path(self.temp.name) / "events.csv"
        self.visits = Path(self.temp.name) / "visits.csv"
        self.header = (
            "visit_id", "site_id", "taxon", "obs_date", "site_class",
            "phenophase", "survey_complete", "observer_id", "effort_minutes"
        )
        write_csv(self.events, ["site_id", "event_date", "event_type", "evidence"], [
            {"site_id": "rail1", "event_date": "2026-06-01", "event_type": "cutting", "evidence": "operator_notice"},
            {"site_id": "rail1", "event_date": "2026-06-10", "event_type": "herbicide", "evidence": "operator_notice"}
        ])

    def test_nonflower_observation_is_not_missing_and_cut_date_matters(self):
        write_csv(self.visits, self.header, [
            dict(zip(self.header, ("v1", "rail1", "Species A", "2026-05-20", "rail", "vegetative", "true", "obs1", "10"))),
            dict(zip(self.header, ("v2", "rail1", "Species A", "2026-06-15", "rail", "flowering", "true", "obs1", "10"))),
            dict(zip(self.header, ("v3", "road1", "Species A", "2026-06-15", "road", "bud", "true", "obs2", "10")))
        ])
        visits, skipped = normalize_visits(self.visits, load_events(self.events))
        self.assertEqual(len(visits), 3)
        self.assertEqual(skipped["incomplete_survey"], 0)
        self.assertEqual(visits[0]["flower_observed"], 0)
        self.assertIsNone(visits[0]["days_since_verified_cut"])
        self.assertEqual(visits[1]["days_since_verified_cut"], 14)
        self.assertIsNone(visits[2]["days_since_verified_cut"])
        self.assertEqual(visits[2]["flower_observed"], 0)
        self.assertEqual(len(weekly_summary(visits)), 3)

    def test_incomplete_and_unknown_are_never_zeros(self):
        write_csv(self.visits, self.header, [
            dict(zip(self.header, ("v1", "rail1", "Species A", "2026-05-20", "rail", "vegetative", "false", "obs1", "10"))),
            dict(zip(self.header, ("v2", "rail1", "Species A", "2026-06-15", "rail", "unknown", "true", "obs1", "10"))),
            dict(zip(self.header, ("v3", "road1", "Species A", "2026-06-15", "road", "flowering", "true", "obs1", "10")))
        ])
        visits, skipped = normalize_visits(self.visits)
        self.assertEqual(len(visits), 1)
        self.assertEqual(skipped["unknown_phenophase"], 1)
        self.assertEqual(skipped["incomplete_survey"], 1)
        self.assertEqual(weekly_summary(visits)[0]["screening_only"], 1)

    def test_same_visit_species_duplicate_fails(self):
        row = dict(zip(self.header, ("v1", "rail1", "Species A", "2026-05-20", "rail", "flowering", "true", "obs1", "10")))
        write_csv(self.visits, self.header, [row, row])
        with self.assertRaises(ValueError):
            normalize_visits(self.visits)

    def test_undocumented_cut_is_not_assumed(self):
        write_csv(self.visits, self.header, [
            dict(zip(self.header, ("v1", "road1", "Species A", "2026-05-20", "road", "flowering", "true", "obs1", "10")))
        ])
        visits, _ = normalize_visits(self.visits, load_events(self.events))
        self.assertIsNone(visits[0]["days_since_verified_cut"])

    def test_weekly_aggregate_keeps_sites_count(self):
        write_csv(self.visits, self.header, [
            dict(zip(self.header, ("v1", "rail1", "Species A", "2026-06-15", "rail", "flowering", "true", "obs1", "10"))),
            dict(zip(self.header, ("v2", "rail2", "Species A", "2026-06-16", "rail", "vegetative", "true", "obs2", "10")))
        ])
        visits, _ = normalize_visits(self.visits)
        groups = weekly_summary(visits)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["n_complete_visits"], 2)
        self.assertEqual(groups[0]["n_sites"], 2)
        self.assertEqual(groups[0]["observed_flowering_proportion"], 0.5)


if __name__ == "__main__":
    unittest.main()

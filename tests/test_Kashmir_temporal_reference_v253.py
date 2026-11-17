"""Synthetic QA for descriptive source-aggregate temporal contrast v25.3."""
import importlib.util
from pathlib import Path
import unittest

MODULE = Path(__file__).resolve().parents[1] / "scripts/audit_Kashmir_temporal_reference_v253.py"
spec = importlib.util.spec_from_file_location("audit_v253", MODULE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def record(group, zones, year, n, raw, eq):
    return dict(group=group, zones=list(zones), year=year, n=n,
                all_species_distance_delta=raw, richness_equalized_delta=eq)


class TemporalAuditTests(unittest.TestCase):
    def setUp(self):
        self.rows = []
        for group in audit.GROUPS:
            for zones in audit.ZONE_PAIRS:
                self.rows.extend([
                    record(group, zones, 2014, 3, 0.10, 0.03),
                    record(group, zones, 2017, 3, 0.25, 0.06),
                ])

    def test_difference_formula(self):
        out = audit.extract({"all_species": self.rows})
        self.assertEqual(len(out["strata"]), 6)
        self.assertEqual(len(out["reference_sensitivity"]), 3)
        self.assertAlmostEqual(out["strata"][0]["raw_temporal_change"], 0.15)
        self.assertAlmostEqual(out["strata"][0]["equalized_temporal_change"], 0.03)
        self.assertAlmostEqual(out["reference_sensitivity"][0]["raw_station_minus_between_temporal_change"], 0.0)
        self.assertFalse(out["paired_site_level_temporal_change_identified"])

    def test_change_is_reference_sensitive(self):
        r = next(r for r in self.rows if r["group"] == "station" and
                 r["zones"] == ["Track", "Interior Plot 1"] and r["year"] == 2017)
        r["all_species_distance_delta"] = -0.10
        out = audit.extract({"all_species": self.rows})
        z = next(t for t in out["reference_sensitivity"] if t["reference_zone"] == "Interior Plot 1")
        self.assertAlmostEqual(z["raw_station_minus_between_temporal_change"], -0.35)

    def test_duplicates_and_missing_fail_closed(self):
        with self.assertRaises(ValueError):
            audit.extract({"all_species": self.rows + self.rows[:1]})
        with self.assertRaises(ValueError):
            audit.extract({"all_species": self.rows[:-1]})

    def test_git_blob(self):
        self.assertEqual(audit.git_blob_sha1(b"hello\n"), "ce013625030ba8dba906f756967f9e9ca394464a")


if __name__ == "__main__":
    unittest.main()

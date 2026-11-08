"""Tests validate the factual Table 4 reanalysis against its published totals.

This is not a new experiment or a test of plant adaptation.
"""
from __future__ import annotations
import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from reanalyze_abe_2012_rail_microhabitat import (
    FIELDS, REGION_TOTALS, read_and_validate, summarize, habitat_enrichment
)

INPUT = ROOT / "data/literature/abe_2012_railway_microhabitat_occurrences.csv"


class HistoricalRailMicrohabitatTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = read_and_validate(INPUT)

    def test_all_published_table4_marginals(self):
        self.assertEqual(len(self.records), 22)
        for region, expected in REGION_TOTALS.items():
            rows = [r for r in self.records if r["region"] == region]
            self.assertEqual(sum(r["total_reported"] for r in rows),
                             expected["total_reported"])
            for habitat in FIELDS:
                self.assertEqual(sum(r[habitat] for r in rows),
                                 expected[habitat])

    def test_horsetail_track_and_slope_observed_mentions(self):
        result = habitat_enrichment(
            self.records, "Equisetum arvense", habitat="trackbed"
        )
        self.assertEqual(result["taxon_occurrence_reports_in_habitat"], 24)
        self.assertEqual(result["all_taxa_occurrence_reports_in_habitat"], 70)
        self.assertEqual(result["taxon_occurrence_reports_all_habitats"], 98)
        self.assertEqual(result["all_taxa_occurrence_reports_all_habitats"], 1334)
        self.assertAlmostEqual(result["descriptive_enrichment"], 4.66706, places=4)
        slope = habitat_enrichment(
            self.records, "Equisetum arvense", habitat="embankment_slope"
        )
        self.assertEqual(slope["taxon_occurrence_reports_in_habitat"], 22)
        self.assertEqual(slope["all_taxa_occurrence_reports_in_habitat"], 642)
        self.assertAlmostEqual(
            result["taxon_share_of_habitat_reports"] /
            slope["taxon_share_of_habitat_reports"], 10.0052, places=3
        )

    def test_warm_and_cold_reproduce_horsetail_enrichment_direction(self):
        warm = habitat_enrichment(
            self.records, "Equisetum arvense", "warm", "trackbed"
        )
        cold = habitat_enrichment(
            self.records, "Equisetum arvense", "cold", "trackbed"
        )
        self.assertEqual(warm["taxon_occurrence_reports_in_habitat"], 7)
        self.assertEqual(warm["all_taxa_occurrence_reports_in_habitat"], 16)
        self.assertEqual(cold["taxon_occurrence_reports_in_habitat"], 17)
        self.assertEqual(cold["all_taxa_occurrence_reports_in_habitat"], 54)
        self.assertGreater(warm["descriptive_enrichment"], 1)
        self.assertGreater(cold["descriptive_enrichment"], 1)

    def test_solidago_has_zero_track_reports_not_entire_railway_absence(self):
        track = habitat_enrichment(
            self.records, "Solidago altissima", habitat="trackbed"
        )
        slope = habitat_enrichment(
            self.records, "Solidago altissima", habitat="embankment_slope"
        )
        self.assertEqual(track["taxon_occurrence_reports_in_habitat"], 0)
        self.assertEqual(slope["taxon_occurrence_reports_in_habitat"], 46)
        self.assertEqual(track["taxon_occurrence_reports_all_habitats"], 90)

    def test_marginal_corruption_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            rows = list(csv.DictReader(INPUT.open(encoding="utf-8")))
            rows[0]["trackbed"] = "3"  # alters both taxon and source subtotal
            corrupted = Path(td) / "table.csv"
            with corrupted.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(ValueError):
                read_and_validate(corrupted)

    def test_top10_absence_not_imputed_as_zero(self):
        out, rows = summarize(self.records)
        warm_imperata = [
            x for x in rows
            if x["taxon_group"] == "Imperata cylindrica" and
            x["region"] == "warm"
        ]
        self.assertEqual(warm_imperata, [])
        self.assertEqual(
            out["focus"]["Imperata cylindrica"].keys(), {"cold","pooled"}
        )


if __name__ == "__main__":
    unittest.main()

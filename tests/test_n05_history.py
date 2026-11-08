"""Synthetic unit tests for true N05 railway history vs GIS attribute changes.

No fabricated biological effects; historical line ID and group counters in
these tests are synthetic, not real Japanese railway records.
"""
from __future__ import annotations
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_n05_historical_railway import known_year,is_current_end,audit_records


class FakeReader:
    fields=[("DeletionFlag","C",1,0),("N05_004","C",4,0),
            ("N05_005b","C",4,0),("N05_005e","C",4,0),
            ("N05_006","C",20,0),("N05_007","C",5,0)]
    def __init__(self, records):
        self.records=records
    def iterRecords(self):
        return iter(self.records)


class N05HistoryTests(unittest.TestCase):
    def test_unknown_and_modern_year(self):
        self.assertEqual(known_year("2004"),2004)
        self.assertEqual(known_year("2024-01-01"),2024)
        self.assertIsNone(known_year("9999"))
        self.assertIsNone(known_year("999"))
        self.assertIsNone(known_year("1949"))
        self.assertIsNone(known_year("garbage"))
        self.assertTrue(is_current_end("9999"))
        self.assertFalse(is_current_end("2020"))

    def test_attribute_change_not_closed_or_new_line(self):
        rows=[
            # old GIS attribute version, same logical line group
            ["2005","2005","2019","line1","attribute-change"],
            # latest still-active version of the SAME logical line
            ["2005","2020","9999","line1","attribute-change"],
            # a truly ended *candidate*, but closure unverified
            ["2011","2011","2020","line2","abolition?"],
            ["2004","2004","9999","line3","active"]
        ]
        out=audit_records(FakeReader(rows),"RailroadSection2")
        self.assertEqual(out["n_feature_rows"],4)
        self.assertEqual(out["n_distinct_relation_groups"],3)
        self.assertEqual(out["n_relation_groups_current_2025"],2)
        self.assertEqual(out["candidate_relation_groups_newly_opened_2005_2025"]["2005"],1)
        self.assertEqual(out["candidate_relation_groups_newly_opened_2005_2025"]["2011"],1)
        self.assertEqual(out["noncurrent_group_terminal_end_years_2005_2024"]["2019"],0)
        self.assertEqual(out["noncurrent_group_terminal_end_years_2005_2024"]["2020"],1)

    def test_invalid_start_end_not_silently_hidden(self):
        out=audit_records(FakeReader([
            ["2010","2015","2012","broken",""],
        ]),"RailroadSection2")
        self.assertEqual(out["source_diagnostics"]["inconsistent_start_after_end"],1)

    def test_missing_schema_rejected(self):
        reader=FakeReader([])
        reader.fields=[("DeletionFlag","C",1,0),("N05_004","C",4,0)]
        with self.assertRaises(ValueError):
            audit_records(reader,"RailroadSection2")


if __name__=="__main__":
    unittest.main()

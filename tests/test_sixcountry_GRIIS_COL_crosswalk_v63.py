"""Synthetic GBIF official one-pass old GRIIS to COL XR cross-country source validation."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from build_sixcountry_GRIIS_COL_crosswalk_v63 import (
    parsed_selected_official_tsv, analyze,source_selected_old_keys, BLOCKS_SIX,
    REQUIRED_COLUMNS,OFFICIAL_CROSSWALK
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_macro_rail_assemblages_v22 import COUNTRIES


def fake_country_source():
    countries=set(COUNTRIES.values())-{"TW"}
    return {"complete":True,"errors":{},"countries":{
        country:{"country":country,"full_source_pages_verified":True,
          "n_directly_matched_plant_backbone_species_keys":110,
          "gbif_backbone_plant_species_key_set":[str(100+i) for i in range(110)]}
        for country in countries
    }}


def make_row(old,col,status="ACCEPTED",colkingdom="Plantae",oldrank="species",colrank="species",accepted=""):
    d={x:"" for x in REQUIRED_COLUMNS}
    d.update({"gbif:ID":str(old),"gbif:rank":oldrank,"gbif:status":"ACCEPTED",
              "col:ID":col,"col:rank":colrank,"col:scientificName":"synTaxon",
              "col:status":status,"col:acceptedID":accepted,"col:kingdom":colkingdom})
    return d


class MultiCountryGBIFCOLMigrationTests(unittest.TestCase):
    def test_all_six_country_GRIIS_lists_used_without_rail_outcome(self):
        q=source_selected_old_keys(fake_country_source())
        self.assertEqual(len(q),6)
        self.assertEqual(set(q),{BLOCKS[b][0] for b in BLOCKS_SIX})
        self.assertTrue(all(len(x)==110 for x in q.values()))
        with self.assertRaises(ValueError):
            source_selected_old_keys(fake_country_source(),blocks=["us_midwest","us_southeast"])

    def test_official_tsv_linewise_source_select_only_first_field_and_no_csv_quoted_multilines(self):
        header=list(REQUIRED_COLUMNS)+["col:classification"]
        one=make_row(100,"ColA")
        two=make_row(101,"ColB",status="SYNONYM",accepted="ColA")
        missing=make_row(9999,"NoMatch")
        stream=["\t".join(header)+"\n"]
        for row in (one,two,missing):
            stream.append("\t".join(row.get(h,'unrestricted "quoted string"') for h in header)+"\n")
        rows,n=parsed_selected_official_tsv(iter(stream),{"100","101"})
        self.assertEqual(n,3)
        self.assertEqual([x["gbif:ID"] for x in rows],["100","101"])
        self.assertEqual(rows[1]["col:acceptedID"],"ColA")

    def test_six_country_taxonomic_counts_conserve_old_keys_and_collapse(self):
        q=source_selected_old_keys(fake_country_source())
        rows=[
            make_row(100,"ColA"),
            make_row(101,"ColB",status="SYNONYM",accepted="ColA"),
            make_row(102,"ColC"),
            make_row(103,"None",colkingdom="Animalia"),
        ]
        out=analyze(q,rows,7800000,{"source":"synthetic"})
        self.assertEqual(out["summary"]["n_independent_country_taxonomic_sources"],6)
        self.assertEqual(out["summary"]["n_old_country_GRIIS_specieskey_entries_across_country_lists"],660)
        self.assertEqual(out["summary"]["n_old_country_specieskeys_direct_strict_COL_match"],12)
        self.assertEqual(out["summary"]["n_old_country_specieskeys_synonym_resolved"],6)
        self.assertEqual(out["summary"]["n_old_country_specieskeys_not_mapped_or_rejected"],642)
        for country,v in out["countries"].items():
            self.assertEqual(v["n_distinct_source_GRIIS_current_STRICT_plus_synonym_COL_species_IDs"],2)
            self.assertEqual(v["old_GRIIS_to_COL_synonym_acceptedID_map"],{"101":"ColA"})
            self.assertEqual(v["old_GRIIS_to_COL_strict_accepted_map"]["100"],"ColA")
        self.assertTrue(out["not_source_original_GBIF_occurrences_or_rail_transport_evidence"])

    def test_one_country_no_matching_crosswalk_fail_closed(self):
        q=source_selected_old_keys(fake_country_source())
        q["FR"]={"does-not-exist"}
        with self.assertRaises(ValueError):
            analyze(q,[make_row(100,"COL1")],7500000,{})

    def test_unexpected_taxonomic_TSV_source_missing_schema_fails(self):
        with self.assertRaises(ValueError):
            parsed_selected_official_tsv(iter(["gbif:ID\tcol:ID\n","100\t2\n"]),{"100"})
        self.assertTrue(OFFICIAL_CROSSWALK.startswith("https://download.checklistbank.org/"))


if __name__=="__main__":
    unittest.main()

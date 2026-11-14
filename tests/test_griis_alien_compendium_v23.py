"""Offline QA of versioned GRIIS alien-status data source gates."""
import hashlib
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_griis_alien_compendium_v23 import (
    EXPECTED_FILENAME,EXPECTED_MD5,COUNTRY_NAMES,pick_csv_source,parse_source
)


class GRIISSourceTests(unittest.TestCase):
    def test_version_and_md5_must_match_archive_metadata(self):
        meta={"files":[{"key":EXPECTED_FILENAME,
            "checksum":"md5:"+EXPECTED_MD5,
            "links":{"self":"https://zenodo.org/api/files/dummy/source.csv"}}]}
        self.assertEqual(pick_csv_source(meta),"https://zenodo.org/api/files/dummy/source.csv")
        meta["files"][0]["checksum"]="md5:00000000000000000000000000000000"
        with self.assertRaisesRegex(ValueError,"checksum"):
            pick_csv_source(meta)

    def test_country_framing_has_ten_geographically_distinct_codes(self):
        self.assertEqual(len(COUNTRY_NAMES),10)
        self.assertEqual(len(set(COUNTRY_NAMES)),10)

    def test_country_alien_checklist_csv_is_not_converted_to_native_absence(self):
        rows=["country,species,kingdom,acceptedNameUsage,gbifKey"]
        for i in range(10020):
            name="Japan" if i%2==0 else "Canada"
            rows.append(f"{name},Species ficta {i},Plantae,Species ficta {i},12345")
        raw=("\n".join(rows)+"\n").encode("utf-8")
        result=parse_source(raw,expected_md5=hashlib.md5(raw).hexdigest())
        self.assertEqual(result["n_original_country_taxon_rows"],10020)
        self.assertEqual(result["n_nonempty_species_binomial"],10020)
        self.assertEqual(result["country_label_quality"]["country"]["target_country_exact_label_counts"]["JP"],5010)
        self.assertEqual(result["country_label_quality"]["country"]["target_country_exact_label_counts"]["CA"],5010)
        self.assertEqual(result["kingdom_column_source_counts"]["kingdom"],{"Plantae":10020})
        self.assertTrue(result["no_native_status_can_be_inferred_from_list_nonmembership"])

    def test_schema_missing_country_fails(self):
        raw=("species,kingdom\n"+"Testus fictus,Plantae\n"*10001).encode()
        with self.assertRaisesRegex(ValueError,"country column"):
            parse_source(raw,expected_md5=hashlib.md5(raw).hexdigest())

    def test_too_small_or_checksum_mismatch_fails(self):
        raw=("country,species\n"+"JP,Testus fictus\n"*6).encode()
        with self.assertRaisesRegex(ValueError,"too-small"):
            parse_source(raw,expected_md5=hashlib.md5(raw).hexdigest())
        with self.assertRaisesRegex(ValueError,"MD5"):
            parse_source(raw)


if __name__=="__main__":
    unittest.main()

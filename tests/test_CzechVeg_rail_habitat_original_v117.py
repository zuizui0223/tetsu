"""v11.7 true CzechVeg source version/hash and site-level railway evidence tests."""
import csv,hashlib,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_CzechVeg_rail_habitat_original_v117 as q


class OriginalCzechVegRailHabitatTests(unittest.TestCase):
    def test_official_record_has_exact_original_header_file(self):
        payload=b"headers-real-binary-source"
        meta={"id":19231152,"files":[
            {"key":q.HEADERS_NAME,"checksum":"md5:"+hashlib.md5(payload).hexdigest(),
             "size":600001,"links":{"self":"https://zenodo.org/api/records/19231152/files/source/content"}}]}
        item=q.original_zenodo_entry(meta)
        self.assertEqual(item["checksum_algorithm"],"md5")
        self.assertEqual(item["source_key"],q.HEADERS_NAME)
        self.assertEqual(item["source_bytes"],600001)

    def test_forged_original_header_file_must_fail(self):
        b=b"SYNTHETIC_NOT_ORIGINAL_SOURCE"
        meta={"id":19231152,"files":[{"key":q.HEADERS_NAME,
            "checksum":"md5:"+hashlib.md5(b).hexdigest(),"size":600001,
            "links":{"self":"https://zenodo.org/api/records/19231152/files/source/content"}}]}
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"headers.csv";f.write_bytes(b)
            with self.assertRaisesRegex(ValueError,"size mismatch"):
                q.verify_original(meta,f)

    def test_railway_only_in_paper_title_not_accepted_as_on_rail(self):
        fields=["RELEVE_NR","SURF_AREA","HABITAT","LOCALITY",
                "REFERENCE_TITLE","YEAR","SKALICKY1988"]
        rows=[
            ["1","10","želEzničNí násep","okolí","Vegetace","2015","A"],
            ["2","25","louka","silnice Brno","Railway flora paper","2020","B"],
            ["3","16","Les","les","Railway survey","2011","A"],
        ]
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"headers.csv"
            with f.open("w",encoding="utf-8",newline="") as o:
                w=csv.writer(o);w.writerow(fields);w.writerows(rows)
            with patch.object(q,"EXPECTED_PLOTS",3):
                got=q.plot_source_scan(f)
        self.assertEqual(got["n_original_railway_site_text_strict_potential_not_independently_field_confirmed"],1)
        self.assertEqual(got["n_original_road_source_site_text"],1)
        self.assertEqual(got["n_reference_only_railway_mentioned_not_accepted_as_rail_plot"],2)
        self.assertEqual(got["strict_candidate_with_measured_positive_source_plot_area_m2"],1)
        self.assertEqual(got["strict_candidate_original_area_median_m2"],10)
        self.assertEqual(got["n_unique_source_ids"],3)

    def test_duplicate_original_releve_is_invalid(self):
        fields=["RELEVE_NR","SURF_AREA","HABITAT","LOCALITY","REFERENCE_TITLE"]
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"headers.csv"
            with f.open("w",encoding="utf-8",newline="") as o:
                w=csv.writer(o);w.writerow(fields)
                w.writerow(["1","20","železnice","",""])
                w.writerow(["1","20","železnice","",""])
            with patch.object(q,"EXPECTED_PLOTS",2):
                with self.assertRaisesRegex(ValueError,"duplicated original Czech"):
                    q.plot_source_scan(f)

if __name__=="__main__":unittest.main()

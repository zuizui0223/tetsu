"""Source-schema and quality tests for USA-NPN endpoint audit."""
import io
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_usa_npn_status import form_args, parse_stream, summarize


class FakeResponse:
    def __init__(self,source):
        self.source=source
    def iter_lines(self):
        return iter(self.source.splitlines())


class NPNSourceAuditTests(unittest.TestCase):
    def test_npn_http_form_matches_rnpn_schema(self):
        args=form_args({"year":2024,"species_id":82,"state":"MD"})
        self.assertEqual(args["request_src"],"tetsu_research_feasibility")
        self.assertEqual(args["species_id[1]"],"82")
        self.assertEqual(args["state[1]"],"MD")
        self.assertEqual(args["start_date"],"2024-01-01")

    def test_ndjson_dict_per_line(self):
        rows, sha, bytes_n=parse_stream(FakeResponse(
            b'{"observation_date":"2024-05-01","phenophase_status":1}\n'
            b'{"observation_date":"2024-05-08","phenophase_status":0}\n'))
        self.assertEqual(len(rows),2)
        self.assertGreater(bytes_n,30)
        summary=summarize(rows,{"year":2024},"abc",bytes_n)
        self.assertEqual(summary["status_raw_value_counts"]["phenophase_status"]["1"],1)
        self.assertEqual(summary["status_raw_value_counts"]["phenophase_status"]["0"],1)

    def test_json_array_one_line(self):
        rows,_,_=parse_stream(FakeResponse(b'[{"phenophase_status":0},{"phenophase_status":1}]'))
        self.assertEqual(len(rows),2)

    def test_empty_array_is_valid_no_data(self):
        rows,_,_=parse_stream(FakeResponse(b'[]'))
        self.assertEqual(rows,[])

    def test_bad_response_fails_not_fake_status(self):
        with self.assertRaises(ValueError):
            parse_stream(FakeResponse(b'["Unexpected message"]'))
        with self.assertRaises(ValueError):
            summarize([{"error":"not data"}],{"state":"MD"},"abc",30)


if __name__=="__main__":
    unittest.main()

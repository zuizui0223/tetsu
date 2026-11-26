import unittest
from audit_railstatus_disturbance_trait_sources import audit_table,csv_table,git_blob_sha
class TestAudit(unittest.TestCase):
    def test_exact_and_missing(self):
        raw=b"species,life\nA alpha,annual\nB beta,NA\nC gamma,perennial\nC gamma,annual\n"
        r=audit_table(raw,{"a alpha":[1],"b beta":[2],"c gamma":[3]}, {})
        self.assertEqual(r["name_candidates"][0]["unique_exact_species"],2)
        self.assertEqual(r["name_candidates"][0]["matched_rows_with_nonmissing_by_column"]["life"],1)
        self.assertEqual(r["column_audit"]["life"]["nonmissing"],3)
    def test_blob(self):
        self.assertEqual(git_blob_sha(b""),"e69de29bb2d1d6434b8b29ae775ad8c2e48c5391")
    def test_duplicate_columns(self):
        with self.assertRaises(ValueError):csv_table(b"a,a\n1,2\n")
if __name__=="__main__":unittest.main()

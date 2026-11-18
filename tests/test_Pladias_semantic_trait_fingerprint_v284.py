"""v28.4 synthetic canonical Pladias trait-value checksum tests."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Pladias_original_semantic_trait_fingerprint_v284 import fingerprint


class PladiasSemanticFingerprintQA(unittest.TestCase):
    def test_stable_to_source_zip_row_reordering(self):
        d1={"Plant Z":"Epilobium","Plant A":"Bidens"}
        l1={"Plant A":"therophyte","Plant Z":"hemicryptophyte"}
        d2={"Plant A":"Bidens","Plant Z":"Epilobium"}
        l2={"Plant Z":"hemicryptophyte","Plant A":"therophyte"}
        self.assertEqual(fingerprint(d1,l1),fingerprint(d2,l2))

    def test_real_category_change_not_silently_accepted(self):
        old=fingerprint({"Plant A":"Epilobium"},{"Plant A":"therophyte"})
        changed_category=fingerprint({"Plant A":"Allium"},{"Plant A":"therophyte"})
        changed_lifeform=fingerprint({"Plant A":"Epilobium"},{"Plant A":"geophyte"})
        self.assertNotEqual(old,changed_category)
        self.assertNotEqual(old,changed_lifeform)

    def test_empty_original_source_traits_fail_closed(self):
        with self.assertRaises(ValueError):fingerprint({},{"P":"therophyte"})
        with self.assertRaises(ValueError):fingerprint({"P":"Bidens"},{})


if __name__=="__main__":
    unittest.main()

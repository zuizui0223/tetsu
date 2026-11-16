"""Official GBIF Backbone 2023 → May 2026 COL26.5 XR source taxonomy crosswalk QA."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from verify_GRIIS_GBIF_to_COL_crosswalk_v62 import (
    map_country_griis_old_keys,sanitized_crosswalk_metadata,
    OFFICIAL_CROSSWALK,COL_VERSION
)


def row(old,col,status="ACCEPTED",colkingdom="Plantae",oldrank="species",
        colrank="species",accept=None):
    return {
        "gbif:ID":str(old),
        "gbif:rank":oldrank,"gbif:status":"accepted",
        "gbif:scientificName":"synthetic_species",
        "col:ID":col,"col:rank":colrank,
        "col:scientificName":"synthetic_col_name",
        "col:status":status,
        "col:acceptedID":accept or "",
        "col:kingdom":colkingdom
    }


class GRIISLegacyCOLBackboneTests(unittest.TestCase):
    def test_accepted_species_keys_from_old_GRIIS_map_to_new_COL(self):
        rows=[
            row(2979180,"8NQF"),  # official docs GBIF Acacia key -> COL
            row(1,"ABCD"),
            row(2,"SYNONYM",status="SYNONYM",accept="Q321"),
            row(3,"none",status="NONE"),
            row(999,"IGNORED")
        ]
        actual=map_country_griis_old_keys(rows,{"2979180","1","2","3"})
        self.assertEqual(actual["strict_accepted"],{"8NQF","ABCD"})
        self.assertEqual(actual["synonym_accepted_additional"],{"Q321"})
        self.assertEqual(actual["n_old_backbone_species_mapping_to_accepted_COL_species"],2)
        self.assertEqual(actual["n_old_backbone_species_mapping_to_COL_synonyms"],1)
        self.assertEqual(actual["n_old_keys_unmapped_or_no_usable_COL"],1)
        self.assertEqual(actual["n_distinct_COL_speciesIDs_with_accepted_synonym_resolution"],3)
        self.assertEqual(actual["n_mapping_rows_scanned"],5)

    def test_synonyms_never_silently_considered_independent_species_concepts(self):
        r=[row(1,"A"),row(2,"B"),row(3,"C",status="SYNONYM",accept="A")]
        got=map_country_griis_old_keys(r,{"1","2","3"})
        self.assertEqual(got["strict_accepted"],{"A","B"})
        self.assertEqual(got["strict_accepted"]|got["synonym_accepted_additional"],{"A","B"})
        self.assertEqual(got["n_distinct_COL_speciesIDs_with_accepted_synonym_resolution"],2)

    def test_many_to_one_collapse_is_explicit_not_species_count_inflation(self):
        rows=[row(1,"COL1"),row(2,"COL1")]
        v=map_country_griis_old_keys(rows,{"1","2"})
        self.assertEqual(v["n_distinct_strict_COL_speciesIDs"],1)
        self.assertEqual(v["n_many_to_one_old_species_collapsed_into_same_STRICT_COL_ID"],1)

    def test_non_species_non_plant_missing_or_duplicate_old_keys_are_rejected(self):
        rows=[row(1,"wrong",colkingdom="Animalia"),
              row(2,"x",oldrank="genus"),
              row(3,"y",colrank="genus"),
              row(4,"",status="ACCEPTED")]
        v=map_country_griis_old_keys(rows,{"1","2","3","4"})
        self.assertEqual(v["n_distinct_strict_COL_speciesIDs"],0)
        self.assertEqual(v["n_old_keys_unmapped_or_no_usable_COL"],4)
        with self.assertRaisesRegex(ValueError,"duplicate legacy"):
            map_country_griis_old_keys([row(1,"A"),row(1,"A")],{"1"})

    def test_crosswalk_raw_map_not_emitted_in_public_country_summary(self):
        v=map_country_griis_old_keys([row(1,"COL123")],{"1"})
        pub=sanitized_crosswalk_metadata(v)
        self.assertNotIn("source_old_to_strict_COL",pub)
        self.assertNotIn("source_old_to_synonym_resolution",pub)
        self.assertNotIn("strict_accepted",pub)
        self.assertEqual(pub["n_distinct_strict_COL_speciesIDs"],1)
        self.assertTrue(OFFICIAL_CROSSWALK.endswith("gbif-col265.tsv.gz"))
        self.assertIn("26.5",COL_VERSION)


if __name__=="__main__":
    unittest.main()

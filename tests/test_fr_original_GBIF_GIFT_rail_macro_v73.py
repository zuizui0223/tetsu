"""Source-only v7.3 tests. Synthetic rows NEVER represent train seed movement."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_fr_original_GBIF_GIFT_rail_macro_v73 import (
    growth_bridge,grouped,source_measure,prepared_cells,GROUPS,SOURCES,ANCHORS,K
)

class StrictOriginalAndFacetSourceTests(unittest.TestCase):
    def test_strict_same_COL_concept_deduplicates_two_legacy_keys(self):
        c={"countries":{"FR":{"old_GRIIS_to_COL_strict_accepted_map":{
            "123":"COL1","124":"COL1","125":"COL2","126":"COL3"}}}}
        g={"complete":True,"source":{"GIFT_meta_ID":"1.2.1"},
           "growth_form_by_source_gbif_species_key":{
               "123":"herbaceous","124":"herbaceous",
               "125":"woody","126":"herbaceous"}}
        known,ambiguous,mapping=growth_bridge(c,g)
        self.assertEqual(known["COL1"],"herbaceous")
        self.assertEqual(len(known),3)
        self.assertEqual(len(ambiguous),0)
        self.assertEqual(len(mapping),4)
        g["growth_form_by_source_gbif_species_key"]["124"]="woody"
        known,ambiguous,_=growth_bridge(c,g)
        self.assertNotIn("COL1",known)
        self.assertEqual(ambiguous,{"COL1"})
    def test_unmatched_taxa_not_assumed_native_or_woody(self):
        x=grouped({"COL1":3,"COL2":4,"COL3":5},{"COL1":"herbaceous","COL2":"woody"})
        self.assertEqual(x["herbaceous"],{"COL1":3})
        self.assertEqual(x["woody"],{"COL2":4})
        with self.assertRaisesRegex(ValueError,"counts invalid"):
            grouped({"COL1":0},{"COL1":"woody"})
    def test_rarefaction_can_only_use_sufficient_true_source(self):
        self.assertEqual(K,20)
        self.assertEqual(set(GROUPS),{"herbaceous","woody"})
        self.assertEqual(len(SOURCES),3)
        self.assertEqual(len(ANCHORS),2)
        self.assertIsNone(source_measure({"a":30,"b":30}))
        self.assertIsNone(source_measure({"a":2,"b":3,"c":4}))
        val=source_measure({"a":20,"b":20,"c":20})
        self.assertEqual(set(val),{"a","b","c"})
        self.assertTrue(all(0<=x<=1 for x in val.values()))
    def test_nonexhaustive_A5_warning_must_be_present(self):
        event={"source_type":"REAL interpreted original GBIF",
               "frozen_geographical_block":"fr_burgundy",
               "source_period_analyzed":"2018_2025",
               "n_original_equal_area_50km_cells":16,
               "source_independently_verified_same_all_Plantae_bbox_row_count":True,
               "source_archive_NOT_complete_regional_original_events_until_A5_partition_gap_resolved":False,
               "NO_RAIL_effect_or_train_dispersal_result":True}
        with self.assertRaisesRegex(ValueError,"provenance"):
            prepared_cells(event,{}, {}, {})
if __name__=="__main__":unittest.main()

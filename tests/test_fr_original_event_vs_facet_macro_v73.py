"""v7.3 source-independent unit checks, no simulated railway ecological result."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_fr_original_event_vs_facet_macro_v73 import (
    assign_groups,col_growth_map,source_summaries,GROUPS)

class OriginalEventVsFacetTests(unittest.TestCase):
    def test_unknown_growth_form_not_silently_zero(self):
        grouped,missing=assign_groups({"COL1":13,"COL2":7,"COL3":2},
                                      {"COL1":"herbaceous","COL2":"woody"})
        self.assertEqual(sum(grouped["herbaceous"].values()),13)
        self.assertEqual(sum(grouped["woody"].values()),7)
        self.assertEqual(missing,2)
    def test_invalid_original_counts_are_rejected(self):
        with self.assertRaisesRegex(ValueError,"Invalid original"):
            assign_groups({"COL1":0},{"COL1":"herbaceous"})
        with self.assertRaisesRegex(ValueError,"Missing source"):
            assign_groups(None,{})
    def test_conflicting_old_names_mapping_same_COL_concept_stays_unknown(self):
        cross={"countries":{"FR":{"old_GRIIS_to_COL_strict_accepted_map":{
            "11":"COL1","12":"COL1","13":"COL2"}}}}
        gift={"complete":True,"source":{"GIFT_meta_ID":"1.2.1"},
              "growth_form_by_source_gbif_species_key":{
                  "11":"woody","12":"herbaceous","13":"herbaceous"}}
        m,amb,old=col_growth_map(cross,gift)
        self.assertNotIn("COL1",m)
        self.assertIn("COL1",amb)
        self.assertEqual(m["COL2"],"herbaceous")
        self.assertEqual(set(old),{"11","12","13"})
    def test_trait_source_coverage_records_unknown_denominator(self):
        cells={"a":{"variants":{},"unknown":{}}}
        for k in ("original_event_all","original_event_precise5km","indexed_facet_2026"):
            cells["a"]["variants"][k]={"herbaceous":{"COL1":8},"woody":{"COL2":2}}
            cells["a"]["unknown"][k]=3
        result=source_summaries(cells)["original_event_all"]
        self.assertEqual(result["record_counts_by_growth_form"]["unknown"],3)
        self.assertAlmostEqual(result["percentage_classified"],76.923,places=3)
    def test_two_precommitted_growth_groups(self):
        self.assertEqual(GROUPS,("herbaceous","woody"))

if __name__=="__main__":unittest.main()

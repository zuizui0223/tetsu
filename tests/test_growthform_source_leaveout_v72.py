"""v7.2 conservative regression tests; synthetic data never imply actual dispersal."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_growthform_source_leaveout_v72 import (
    COUNT_KEYS, GROUPS, fixed_group_source_cells, equal_country_summary)

def source_fixture():
    xs={"11":20,"12":16,"13":15,"21":21,"22":20,"23":12}
    cell={"original_two_period_country_GRIIS_eligible":True,"periods":{
        "2018_2025":{"dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys":True,
                     COUNT_KEYS["full"]:dict(xs),
                     COUNT_KEYS["without_dominant_dataset"]:dict(xs)}}}
    return {"source_cells":{"a":cell,**{f"dummy_{i}":cell for i in range(15)}}}

class SourceSensitivityContractTests(unittest.TestCase):
    def test_both_source_variants_and_both_growth_forms_required(self):
        t={str(k):"herbaceous" if k<20 else "woody" for k in [11,12,13,21,22,23]}
        old={"a":{"2018_2025":{"counts":{},"all_plant_indexed_records":100}}}
        self.assertEqual(set(fixed_group_source_cells("us_midwest","2018_2025",
                         {"regions":{"us_midwest":{"cells":old}}},source_fixture(),t)),{"a"})
        f=source_fixture()
        for k in ("21","22","23"):f["source_cells"]["a"]["periods"]["2018_2025"][COUNT_KEYS["without_dominant_dataset"]].pop(k)
        self.assertFalse(fixed_group_source_cells("us_midwest","2018_2025",
                          {"regions":{"us_midwest":{"cells":old}}},f,t))
    def test_incomplete_GBIF_record_partition_fails_closed(self):
        t={"11":"herbaceous","21":"woody"}
        f=source_fixture();f["source_cells"]["a"]["periods"]["2018_2025"].pop(
            "dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys")
        with self.assertRaisesRegex(ValueError,"dataset split"):
            fixed_group_source_cells("us_midwest","2018_2025",
               {"regions":{"us_midwest":{"cells":{"a":{"2018_2025":{}}}}}},f,t)
    def test_source_category_names_and_not_assuming_native(self):
        self.assertEqual(set(COUNT_KEYS),{"full","without_dominant_dataset"})
        self.assertEqual(GROUPS,("herbaceous","woody"))
        self.assertNotIn("native",COUNT_KEYS)

if __name__=="__main__":unittest.main()

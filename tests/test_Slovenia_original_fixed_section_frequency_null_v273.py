"""v27.3 synthetic source-rail section local taxa frequency-null checks."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Slovenia_original_rail_15sections_fixed_local_frequency_null_v273 as v


def sample_original():
    source={}
    for section,size in (("s1",1),("s2",2),("s3",4),("s4",6),("s5",5)):
        for i in range(size):
            source[f"{section}_site{i}"]={
                "section":section,"matched_sides":1 if i%2 else 2,
                "BS":{1,i+10},"BT":{1,i+10,i+20}
            }
    return source


class SloveniaSourceSectionNullTests(unittest.TestCase):
    def test_source_15section_grouping_and_site_support(self):
        orig=sample_original()
        a=v.groups_from_original(orig,2)
        b=v.groups_from_original(orig,4)
        self.assertEqual(set(a),{"s2","s3","s4","s5"})
        self.assertEqual(set(b),{"s3","s4","s5"})
        self.assertEqual(sum(map(len,a.values())),17)
        self.assertEqual(sum(map(len,b.values())),15)
        with self.assertRaises(ValueError):
            v.groups_from_original(orig,1)

    def test_source_species_row_and_col_margins_and_repeatable_aggregates(self):
        grouped=v.groups_from_original(sample_original(),2)
        out=v.analyze_group(grouped,draws=19,burnin=100,spacing=20,seed=111)
        self.assertTrue(out["identified"])
        self.assertEqual(out["source_original_rail_sections"],4)
        self.assertEqual(out["source_original_spatial_locations"],17)
        q=out["source_equal_section_global_metric"]["beta_jac"]
        self.assertGreaterEqual(q["observed"],0.)
        self.assertLessEqual(q["observed"],1.)
        self.assertEqual(
            q["frequency_null_mean"],
            v.analyze_group(grouped,draws=19,burnin=100,spacing=20,seed=111)["source_equal_section_global_metric"]["beta_jac"]["frequency_null_mean"])

    def test_original_location_sides_not_pseudoreplicated(self):
        original=sample_original()
        original["s3_site0"]["matched_sides"]=2
        a=v.groups_from_original(original,2)
        self.assertEqual(len(a["s3"]),4)
        self.assertEqual(a["s3"][0]["BS"],{1,10})
        self.assertEqual(a["s3"][0]["BT"],{1,10,20})
        original["s3_site0"]["matched_sides"]=3
        with self.assertRaises(ValueError):v.groups_from_original(original,2)

    def test_no_cross_rail_section_shuffling_when_source_zero_margins(self):
        data={}
        for section in ("s1","s2"):
            data[section]=[{"original_location":f"{section}_{i}",
                             "BS":{1,i+50},"BT":{1,2}} for i in range(4)]
        out=v.analyze_group(data,draws=19,burnin=20,spacing=10)
        self.assertEqual(out["source_original_rail_sections"],2)
        self.assertEqual(out["n_source_sections_with_no_accepted_swaps"],2)


if __name__=="__main__":unittest.main()

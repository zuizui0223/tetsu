"""v17.6 synthetic-only queue math. Original 4,235 Czech iNaturalist photo
observations reloaded separately in GitHub Actions with source verification.
"""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from create_Czech_Echium_all_public_obs_label_REVIEW_queue_v176 import (
    token_id,create_queue,export_csv,verify_source,FIRST_WAVE_MAX_CONTROL,
)
from test_Czech_Echium_all_photo_annotations_no_rights_v175 import synthetic_record

class FakeOriginalRailwayPolyline:
    def distance_and_cell(self,lon,lat):
        if lon<13.2:return 10.,(1,1)
        if lon<13.9:return 200.,(3,3)
        return 5000.,(5,5)

def full_fake():
    rows=[]
    rows.extend(synthetic_record(i,lon=13.11,accuracy=5,annotation=13 if i%2 else None) for i in range(178))
    rows.extend(synthetic_record(178+i,lon=13.25,accuracy=5,annotation=21 if i%3 else None) for i in range(725))
    rows.extend(synthetic_record(903+i,lon=14,accuracy=5) for i in range(3273))
    for i in range(59):
        o=synthetic_record(4176+i,lon=14,accuracy=5)
        o["geoprivacy"]="private"
        rows.append(o)
    assert len(rows)==4235
    return rows

class AllPublicObservationLabelQueueTests(unittest.TestCase):
    def test_source_only_public_ids_and_no_licensed_image_filter(self):
        orig=full_fake()
        for r in orig:r["photos"][0]["license_code"]="all-rights-reserved"
        review,analysis,first,info=create_queue(orig,FakeOriginalRailwayPolyline(),{})
        self.assertEqual(len(review),4176)
        self.assertEqual(len(analysis),4176)
        self.assertEqual(len(first),356)
        self.assertEqual(info["n_available_original_near30_error30_source_records"],178)
        self.assertEqual(info["n_available_original_control_100to1000_error30_source_records"],725)
        self.assertEqual(info["n_FIRSTWAVE_source_control_100to1000_observations"],178)
        self.assertTrue(info["source_REVIEW_fields_left_all_blank"])
        self.assertTrue(info["photo_pixels_NEVER_downloaded_or_redistributed"])
        self.assertTrue(info["links_to_original_PUBLIC_observation_pages_can_reveal_originally_PUBLIC_localities"])
        self.assertTrue(all(x["human_REVIEW_open_flower_visible_yes_no_uncertain"]=="" for x in review))

    def test_review_source_sheet_OMITS_photo_annotations_and_rail_distance_to_reduce_assessor_bias(self):
        review,analysis,first,info=create_queue(full_fake(),FakeOriginalRailwayPolyline(),{})
        sample=review[0]
        self.assertIn("public_observation_id",sample)
        self.assertIn("public_iNaturalist_observation_page",sample)
        self.assertNotIn("source_existing_iNat_photo_annotation_EXPLICIT_only",sample)
        self.assertNotIn("source_nearest_real_2023_OSM_rail_LINE_distance_m",sample)
        context=analysis[0]
        self.assertIn("source_existing_iNat_photo_annotation_EXPLICIT_only",context)
        self.assertIn("source_nearest_real_2023_OSM_rail_LINE_distance_m",context)
        self.assertIn("source_PUBLIC_observer_pseudonym_for_cluster_adjustment",context)
        self.assertNotIn("lat",context)
        self.assertNotIn("lon",context)
        self.assertNotIn("original_observation_id",context)
        self.assertNotIn("photo_url",context)
        self.assertEqual(context["source_original_DOY"],131)
        self.assertFalse(info["not_fully_blinded_human_evaluation_by_design"] is False)

    def test_nonrandom_photo_labels_not_used_to_pick_firstwave(self):
        original=full_fake()
        review,context,first,info=create_queue(original,FakeOriginalRailwayPolyline(),{})
        chosen={r["review_token"] for r in first}
        for r in original:
            r["annotations"]=[{"controlled_attribute_id":12,"controlled_value_id":21}]
        r2,c2,first2,info2=create_queue(original,FakeOriginalRailwayPolyline(),{})
        self.assertEqual(chosen,{r["review_token"] for r in first2})

    def test_huge_public_id_duplicate_or_fake_source_count_fails_closed(self):
        src=full_fake()
        src[1]["id"]=src[0]["id"]
        with self.assertRaisesRegex(ValueError,"Repeated source public"):
            create_queue(src,FakeOriginalRailwayPolyline(),{})
        with self.assertRaisesRegex(ValueError,"4235"):
            create_queue(src[:-1],FakeOriginalRailwayPolyline(),{})
        with self.assertRaisesRegex(ValueError,"numeric identifier"):
            token_id(None)

    def test_review_csv_has_public_links_not_original_photo_pixels(self):
        import tempfile,csv
        review,context,first,info=create_queue(full_fake(),FakeOriginalRailwayPolyline(),{})
        with tempfile.TemporaryDirectory() as p:
            target=Path(p)/"public_obs.csv"
            export_csv(target,first)
            with target.open(newline="",encoding="utf8") as f:
                records=list(csv.DictReader(f))
            self.assertEqual(len(records),356)
            self.assertTrue(all(r["public_iNaturalist_observation_page"].startswith("https://www.inaturalist.org/observations/")
                                for r in records))
            self.assertTrue(all("photo" not in r["human_REVIEW_open_flower_visible_yes_no_uncertain"]
                                for r in records))

if __name__=="__main__":unittest.main()

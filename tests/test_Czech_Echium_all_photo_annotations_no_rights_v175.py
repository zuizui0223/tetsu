"""v17.5 tests use ONLY synthetic iNaturalist JSON metadata fixtures.
No real photo, botany/railway location or original plant phenology is
generated; authentic original 4,235 are re-fetched in separate source CI.
"""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Czech_Echium_all_iNat_annotations_no_photo_rights_gate_v175 import (
    TARGET,SOURCE_V166,EXPECTED_TOTAL,annotation_class,public_metadata_gate,
    process_records,phenotype_gate,validate_prior
)

class FakeOriginalRailwayPolyline:
    def distance_and_cell(self,lon,lat):
        # Source geometry mocked only for math tests: not actual Czech rails.
        if lon<13.1:return 0.,(1,1)
        if lon<13.2:return 10.,(2,2)
        if lon<13.3:return 200.,(3,3)
        return 5000.,(4,4)

def synthetic_record(index,lon=13.11,accuracy=5,annotation=None):
    r={"id":index,"taxon":{"name":TARGET},"quality_grade":"research",
       "observed_on":"2024-05-11","obscured":False,"geoprivacy":"open",
       "captive":False,"positional_accuracy":accuracy,
       "geojson":{"coordinates":[lon,49.1]},"photos":[{"license_code":None}],
       "annotations":[]}
    if annotation is not None:
        ids=annotation if isinstance(annotation,(list,tuple)) else [annotation]
        r["annotations"]=[{"controlled_attribute_id":12,"controlled_value_id":v} for v in ids]
    return r

def source_fixture():
    return {"version":SOURCE_V166,"complete":True,
      "source_original_publications_and_Geofabrik_SHA":{
          "CZ":{"original_Geofabrik_PBF_SHA256":"26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"}
      },
      "source_species_country_counts":[{
          "taxon":TARGET,"region_source_bbox":"CZ",
          "source_bbox_iNat_research_photo_total":4235,
          "n_rail_polyline_point_near_30m":178,
          "n_rail_polyline_guaranteed_circle_within_30m":153,
          "n_source_species_photo_100to1000m_controls":725,
      }]}

def synthetic_complete_source():
    # EXACT original source count math in fake geometry:
    # 153 point 10m and accuracy 5m => circle safely inside 30m
    # 25 point 10m and accuracy25m => point-near only (10+25>30)
    # 725 controls 200m away 5m accuracy
    # 3332 others >1000m, each with 5m coordinates
    source=[synthetic_record(i,lon=13.11,accuracy=5,
                             annotation=13 if i%2==0 else 21)
            for i in range(153)]
    source.extend(synthetic_record(153+i,lon=13.11,accuracy=25,annotation=None)
                  for i in range(25))
    source.extend(synthetic_record(178+i,lon=13.25,accuracy=5,
                                  annotation=13 if i%2==0 else 21)
                  for i in range(725))
    source.extend(synthetic_record(903+i,lon=14.0,accuracy=5,annotation=None)
                  for i in range(3332))
    assert len(source)==EXPECTED_TOTAL
    return source

class AllPhotoOriginalCzechEchiumLabelSafety(unittest.TestCase):
    def test_no_license_rights_filter_applied_to_public_METADATA(self):
        a=synthetic_record(1,annotation=13)
        self.assertEqual(public_metadata_gate(a)[1],"pass")
        a["photos"][0]["license_code"]="all-rights-reserved"
        self.assertEqual(public_metadata_gate(a)[1],"pass")

    def test_accurate_OFFICIAL_annotation_open_bud_fruits_no_and_unannotated(self):
        for annotation,expected in [
          (13,"photo_open_flower_visible"),
          (21,"photo_no_flower_fruit_seed_visible"),
          (15,"photo_buds_visible_but_open_flower_unknown"),
          (14,"photo_fruits_visible_but_open_flower_unknown"),
          ([13,21],"contradictory_photo_annotations_exclude"),
          (None,"unannotated_open_flower_UNKNOWN_NOT_NONFLOWERING")
        ]:
            self.assertEqual(annotation_class(synthetic_record(1,annotation=annotation)),expected)

    def test_public_accuracy_unknown_or_above30_not_silently_discarded_as_no_flower(self):
        a=synthetic_record(1,accuracy=400,annotation=13)
        v,q=public_metadata_gate(a)
        self.assertEqual(q,"pass")
        self.assertEqual(v["accuracy_class"],"accuracy_100to1000m")
        a["positional_accuracy"]=None
        v,q=public_metadata_gate(a)
        self.assertEqual(q,"pass")
        self.assertEqual(v["accuracy_class"],"accuracy_unknown")
        a["geoprivacy"]="private"
        self.assertEqual(public_metadata_gate(a)[1],"obscured_or_private")

    def test_4235_complete_synthetic_strict_original_source_geo_check_math(self):
        source=synthetic_complete_source()
        prior=validate_prior(source_fixture())
        r=process_records(source,FakeOriginalRailwayPolyline(),prior)
        self.assertTrue(r["full_original_API_photo_record_count_and_strict_geography_v166_reproduced"])
        self.assertEqual(
            sum(r["original_railway_polyline_distance_zone_anthropogenic_observation_metadata"]["strict_near30"].values()),178)
        self.assertEqual(
            sum(r["original_railway_polyline_distance_zone_anthropogenic_observation_metadata"]["guaranteed_circle_near30"].values()),153)
        self.assertEqual(
            sum(r["original_railway_polyline_distance_zone_anthropogenic_observation_metadata"]["strict_control100to1000"].values()),725)
        self.assertTrue(r["NO_photo_image_copyright_restriction_used_for_METADATA"])
        self.assertTrue(r["no_unannotated_photo_coded_negative"])
        self.assertNotIn("observation_id",str(r))
        self.assertNotIn("photo_url",str(r))

    def test_fake_4000_photos_or_bad_original_geo_audit_refused_not_imputed(self):
        prior=validate_prior(source_fixture())
        with self.assertRaisesRegex(ValueError,"not 4235"):
            process_records(synthetic_complete_source()[:-1],FakeOriginalRailwayPolyline(),prior)
        bad=source_fixture()
        bad["source_species_country_counts"][0]["n_rail_polyline_point_near_30m"]=179
        with self.assertRaisesRegex(ValueError,"strict photo\\+OSM"):
            process_records(synthetic_complete_source(),FakeOriginalRailwayPolyline(),
                            validate_prior(bad))

    def test_no_false_research_model_just_from_positive_tags_or_unannotated(self):
        from collections import Counter
        s={"strict_near30":Counter({"photo_open_flower_visible":75,
                                     "unannotated_open_flower_UNKNOWN_NOT_NONFLOWERING":200}),
           "strict_control100to1000":Counter({"photo_open_flower_visible":100,
                                      "unannotated_open_flower_UNKNOWN_NOT_NONFLOWERING":200})}
        m={"strict_near30":{"05":100},"strict_control100to1000":{"05":100}}
        f=phenotype_gate(s,m)
        self.assertFalse(f["near30_accuracy30"]["metadata_photo_label_model_gate_passed"])
        self.assertFalse(f["control100to1000_accuracy30"]["metadata_photo_label_model_gate_passed"])
        # Only explicit negative 21 changes inferability.
        s["strict_near30"]["photo_no_flower_fruit_seed_visible"]=40
        s["strict_control100to1000"]["photo_no_flower_fruit_seed_visible"]=40
        m["strict_near30"]={"01":1,"03":1,"05":1,"06":1,"08":1}
        m["strict_control100to1000"]={"01":1,"03":1,"05":1,"06":1,"08":1}
        g=phenotype_gate(s,m)
        self.assertTrue(g["near30_accuracy30"]["metadata_photo_label_model_gate_passed"])
        self.assertTrue(g["control100to1000_accuracy30"]["metadata_photo_label_model_gate_passed"])

if __name__=="__main__":unittest.main()

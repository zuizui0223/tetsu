"""v16.8 synthetic-only local matched photo opportunities. No field data in
these fixtures are actual plant photos, original user accounts or GPS.
"""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Czech_four_flower_taxa_local_near30_vs_control_photo_pairs_v168 import (
    review_priority,month_from_original_record,has_licensed_photo,
    _private_real_record,one_to_one_matched,summarize_source_species,
    PRIORITY_VERSION,MAX_LOCAL_DISTANCE_KM
)

class RailStub:
    def distance_and_cell(self,lon,lat):
        if lon<14.00001:return 0.0,(10,10)
        if lon<14.02:return 200.0,(12,12)
        return 2000.0,(300,400)

def photorecord(i,name="Echium vulgare",lon=14.0,lat=49.0,
                yearmonth="2024-05",user=None,accuracy=5):
    return {
      "id":i,"taxon":{"name":name},
      "observed_on":yearmonth+"-15",
      "obscured":False,"geoprivacy":"open","captive":False,
      "positional_accuracy":accuracy,
      "geojson":{"coordinates":[lon,lat]},
      "photos":[{"license_code":"cc-by"}],
      "annotations":[{"controlled_attribute_id":12,"controlled_value_id":13}],
      "user":{"id":user or i+10000}
    }

def synthetic_record(i,role,lat=49.0,lon=14.0,month="2024-05",uid=None):
    return {
      "internal_original_id":i,"internal_original_user_id":str(uid or i+10000),
      "lat":lat,"lon":lon,
      "internal_source_100m_cell":(int(i/2),0) if role=="near" else (int(i/2),5),
      "year_month":month,
      "group":role,
      "annotated_flowers_positive":i%2==0
    }

class FourRealTaxaLocalPhotoControlOpportunityTests(unittest.TestCase):
    def test_candidate_source_original_v167_species_selection_rejected_if_not_exact_four(self):
        src={
          "complete":True,"version":PRIORITY_VERSION,
          "source_v166_author_photo_gis_ID":"43693580409",
          "source_n_real_eligible_photo_gates_v166":11,
          "n_species_country_cells_passing_exploratory_ANGIOSPERM_image_review_gate":10,
          "exploratory_ANGIOSPERM_image_review_priority":[
            {"species":"Echium vulgare","region":"CZ","near30_original_iNat_positive_FLOWERS_annotation":21},
            {"species":"Geranium robertianum","region":"CZ","near30_original_iNat_positive_FLOWERS_annotation":7},
            {"species":"Erigeron annuus","region":"CZ","near30_original_iNat_positive_FLOWERS_annotation":5},
            {"species":"Linaria vulgaris","region":"CZ","near30_original_iNat_positive_FLOWERS_annotation":7},
            *[{"species":"Other species"+str(i),"region":"CZ",
                "near30_original_iNat_positive_FLOWERS_annotation":4} for i in range(6)]
          ],
        }
        self.assertEqual(len(review_priority(src)),4)
        src["exploratory_ANGIOSPERM_image_review_priority"][4]["near30_original_iNat_positive_FLOWERS_annotation"]=6
        with self.assertRaisesRegex(ValueError,"original four"):
            review_priority(src)

    def test_public_record_quality_and_true_rail_buffer_circle(self):
        obs=photorecord(1)
        r,status=_private_real_record(obs,RailStub(),"Echium vulgare")
        self.assertEqual(status,"pass")
        self.assertEqual(r["group"],"near30_certain")
        self.assertEqual(month_from_original_record(obs),"2024-05")
        self.assertTrue(has_licensed_photo(obs))
        a=photorecord(2,lon=14.01)
        r,status=_private_real_record(a,RailStub(),"Echium vulgare")
        self.assertEqual(r["group"],"control_100to1000")
        a["photos"][0]["license_code"]="all-rights-reserved"
        r,status=_private_real_record(a,RailStub(),"Echium vulgare")
        self.assertIsNone(r)
        self.assertEqual(status,"original_photo_no_allowed_CC_licence")

    def test_two_groups_local_same_calendar_year_month_are_1_to_1_and_no_reused_control(self):
        near=[synthetic_record(i,"near",lon=14.+i*.0001) for i in range(1,6)]
        control=[synthetic_record(100+i,"control",lon=14.005+i*.0002) for i in range(3)]
        pairs=one_to_one_matched(near,control)
        self.assertEqual(len(pairs),3)
        self.assertEqual(len({a["internal_original_id"] for a,b in pairs}),3)
        self.assertEqual(len({b["internal_original_id"] for a,b in pairs}),3)

    def test_photo_records_different_month_or_gt10km_must_not_match(self):
        n=synthetic_record(1,"near")
        wrongmonth=synthetic_record(10,"control",month="2024-06")
        far=synthetic_record(11,"control",lat=49.5,lon=14.5)
        self.assertEqual(one_to_one_matched([n],[wrongmonth,far]),[])
        close=synthetic_record(12,"control",lon=14.01)
        self.assertEqual(len(one_to_one_matched([n],[close])),1)

    def test_repeated_same_observer_can_be_flagged_as_different_observer_sensitivity(self):
        n=synthetic_record(1,"near",uid=555)
        c=synthetic_record(2,"control",lon=14.01,uid=555)
        self.assertEqual(len(one_to_one_matched([n],[c],restrict_different_observers=False)),1)
        self.assertEqual(len(one_to_one_matched([n],[c],restrict_different_observers=True)),0)

    def test_real_source_photo_aggregate_no_original_ids_or_coordinates_exported(self):
        obs=[photorecord(1),photorecord(2,lon=14.01),
             photorecord(3,lon=14.01,yearmonth="2024-06")]
        # Control band stub conserves non-source 30m metadata.
        prior={
          "source_bbox_iNat_research_photo_total":3,
          "original_public_accuracy_circle_guaranteed_near30":1,
          "original_CC_photo_near30":1,
          "original_photo_controls_100to1000m":2
        }
        summ=summarize_source_species(obs,RailStub(),"Echium vulgare",prior)
        self.assertEqual(summ["n_source_point_near30_and_errorcircle_certain_and_CC_licensed"],1)
        self.assertEqual(summ["n_same_year_month_unique_photo_pairs_within10km_different_observers"],1)
        self.assertEqual(summ["n_matched_source_distinct_months"],1)
        self.assertTrue(summ["NO_raw_original_photos_observation_ID_public_coordinates_or_usernames_exported"]
                        if "NO_raw_original_photos_observation_ID_public_coordinates_or_usernames_exported" in summ
                        else summ["image_pixels_original_photo_and_geo_or_observer_ID_NOT_exposed"])
        self.assertNotIn("original_ids",summ)
        self.assertNotIn("lat",summ)
        self.assertNotIn("lon",summ)
        self.assertNotIn("internal_original_user_id",summ)

    def test_no_strict_cc_or_exact_source_photo_total_mismatch(self):
        prior={
          "source_bbox_iNat_research_photo_total":1,
          "original_public_accuracy_circle_guaranteed_near30":0,
          "original_CC_photo_near30":0,
          "original_photo_controls_100to1000m":1
        }
        with self.assertRaisesRegex(ValueError,"source result count"):
            summarize_source_species([photorecord(1),photorecord(2)],RailStub(),"Echium vulgare",prior)
        with self.assertRaisesRegex(ValueError,"beyond original documented"):
            summarize_source_species([photorecord(1)],RailStub(),"Echium vulgare",prior)

if __name__=="__main__":unittest.main()

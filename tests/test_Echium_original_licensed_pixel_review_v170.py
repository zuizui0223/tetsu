"""v17.0 SYNTHETIC-only test pixel rights/source and blind review invariants.
Fake rectangle JPEGs and observation records below are NOT botanical field
photos, no observed railway phenotype, flowering onset or plant height.
"""
from __future__ import annotations
import csv,io,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from PIL import Image

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from build_Echium_29_railphoto_blinded_CC_contact_sheets_v170 import (
    TARGET,EXPECTED_REAL_PAIRS,EXPECTED_ORIGINAL_API_TOTAL,
    source_gate,eligible_public_photo,blind_order,
    safe_actual_image_download,choose_review_pairs,compose_blind_pair,run,
    SOURCE_V168,SOURCE_V166
)
class SourceRailStub:
    def distance_and_cell(self,lon,lat):
        return (0.0,(0,0)) if lon<=14.00001 else (200.0,(20,20))

def fake_ob(i,near=True,license_code="cc-by"):
    return {
       "id":i if near else i+1000,
       "taxon":{"name":TARGET},
       "observed_on":"2024-05-15","obscured":False,
       "geoprivacy":"open","captive":False,"positional_accuracy":5,
       "geojson":{"coordinates":[14.0 if near else 14.001,49.0]},
       "photos":[{
           "url":f"https://static.inaturalist.org/photos/{100000+i if near else 200000+i}/square.jpg",
           "license_code":license_code,"attribution":"© Synthetic test author (CC BY)"
       }],
       "user":{"id":i if near else i+10000},
       "annotations":[]
    }

def source_fixtures():
    a={
      "complete":True,"version":SOURCE_V168,
      "source_Geofabrik_original_Czech2023_PBF_SHA":
          "26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747",
      "NO_raw_original_photos_observation_ID_public_coordinates_or_usernames_exported":True,
      "source_restricted_original_4_taxa_matched_opportunity":[{
          "source_species":TARGET,
          "n_same_year_month_unique_photo_pairs_within10km_different_observers":29,
          "source_2018_2025_full_real_API_photograph_observation_count":4235,
          "n_source_point_near30_and_errorcircle_certain_and_CC_licensed":29,
      }]
    }
    b={"complete":True,"version":SOURCE_V166,
       "source_species_country_counts":[{
           "region_source_bbox":"CZ","taxon":TARGET,
           "source_bbox_iNat_research_photo_total":4235,
           "n_rail_polyline_guaranteed_circle_within_30m":29,
           "n_rail_polyline_30m_with_CC_photo_license":29,
           "n_source_species_photo_100to1000m_controls":29
       }]}
    return a,b

def fake_source_photos():
    chosen=[fake_ob(i) for i in range(1,30)]
    chosen += [fake_ob(i,near=False) for i in range(1,30)]
    # Rest of source API 4,235 fakes intentionally fail quality/taxon
    # gate rather than impute missing records as rail-site absences.
    for i in range(4235-len(chosen)):
        chosen.append({
          "id":i+20000,
          "taxon":{"name":"Other source taxon"},
          "geoprivacy":"open","photos":[{"license_code":"cc-by"}],
          "observed_on":"2024-05-15","positional_accuracy":5,
          "geojson":{"coordinates":[14.,49.]}
        })
    return chosen

class FakeResponse:
    def __init__(self,body,url):
        self.body=body;self.url=url
    def __enter__(self):return self
    def __exit__(self,*unused):return False
    def raise_for_status(self):return None
    def iter_content(self,n=65536):
        yield self.body

class FakeSession:
    def __init__(self,body):
        self.body=body;self.headers={}
    def get(self,url,**unused):
        return FakeResponse(self.body,url)

class FakeAPI:
    def __init__(self,records):
        self.records=records
    def get(self,params):
        page=int(params.get("page",1))
        per_page=int(params.get("per_page",200))
        start=(page-1)*per_page
        return {"total_results":len(self.records),
                "results":self.records[start:start+per_page]}

def imgbytes():
    raw=io.BytesIO()
    Image.new("RGB",(360,280),(35,112,77)).save(raw,format="JPEG")
    return raw.getvalue()

class OriginalEchiumLicensedPixelsTests(unittest.TestCase):
    def test_real_source_29_pair_count_and_original_CZ_2023_SHA_not_invented(self):
        a,b=source_fixtures()
        source_gate(a,b)
        a["source_restricted_original_4_taxa_matched_opportunity"][0][
            "n_same_year_month_unique_photo_pairs_within10km_different_observers"]=30
        with self.assertRaisesRegex(ValueError,"29 matched"):
            source_gate(a,b)

    def test_only_CC0_or_CC_BY_original_public_iNaturalist_hosts_and_real_attribution(self):
        ob=fake_ob(1)
        x=eligible_public_photo(ob)
        self.assertEqual(x["license_code"],"cc-by")
        self.assertIn("/medium.jpg",x["source_url"])
        self.assertIn("Synthetic test author",x["attribution"])
        ob["photos"][0]["license_code"]="cc-by-nc"
        self.assertIsNone(eligible_public_photo(ob))
        ob["photos"][0]["license_code"]="cc-by"
        ob["photos"][0]["attribution"]=""
        self.assertIsNone(eligible_public_photo(ob))
        ob["photos"][0]["attribution"]="© example"
        ob["photos"][0]["url"]="https://attacker.example/photos/fake/square.jpg"
        self.assertIsNone(eligible_public_photo(ob))

    def test_original_source_29_metadata_pairs_are_rebuilt_once_before_rights_filter(self):
        a,b=source_fixtures()
        original=fake_source_photos()
        pairs,qc=choose_review_pairs(original,SourceRailStub(),
                     a["source_restricted_original_4_taxa_matched_opportunity"][0],
                     b["source_species_country_counts"][0])
        self.assertEqual(len(pairs),29)
        self.assertEqual(qc["n_original_source_photo_pairs"],29)
        self.assertEqual(qc["n_public_both_sides_CC0_CC_BY_attributed_photo_pairs"],29)
        original[0]["photos"][0]["license_code"]="cc-by-nc"
        subset,qc=choose_review_pairs(original,SourceRailStub(),
                     a["source_restricted_original_4_taxa_matched_opportunity"][0],
                     b["source_species_country_counts"][0])
        self.assertEqual(len(subset),28)

    def test_actual_public_JPEG_bytes_are_decoded_and_hash_logged_not_fake_size(self):
        info=eligible_public_photo(fake_ob(1))
        image,metrics=safe_actual_image_download(FakeSession(imgbytes()),info)
        self.assertEqual(image.size,(360,280))
        self.assertEqual(metrics["source_original_image_dimensions_px"],[360,280])
        self.assertEqual(len(metrics["source_original_public_image_file_sha256"]),64)
        bad=dict(info);bad["source_url"]="https://evil.example/evil.jpg"
        with self.assertRaisesRegex(ValueError,"Blocked"):
            safe_actual_image_download(FakeSession(imgbytes()),bad)

    def test_contact_sheet_is_blind_and_source_attribution_is_separate(self):
        left=Image.new("RGB",(400,300),(30,144,30))
        right=Image.new("RGB",(300,400),(144,60,30))
        jpg=compose_blind_pair(left,right,1)
        self.assertEqual(jpg.size,(900,530))
        self.assertNotEqual(jpg.getpixel((100,100)),jpg.getpixel((500,100)))

    def test_full_source_requery_creates_two_BLIND_pairs_and_NO_original_private_identifiers(self):
        a,b=source_fixtures()
        records=fake_source_photos()
        with tempfile.TemporaryDirectory() as td,patch(
            "build_Echium_29_railphoto_blinded_CC_contact_sheets_v170.requests.Session",
            return_value=FakeSession(imgbytes())
        ),patch(
            "build_Echium_29_railphoto_blinded_CC_contact_sheets_v170.time.sleep",
            return_value=None
        ):
            result=run(a,b,SourceRailStub(),td,
                       client=FakeAPI(records),max_pairs=2)
            self.assertEqual(result["n_licensed_blinded_true_source_image_pairs_for_HUMAN_REVIEW"],2)
            self.assertEqual(result["n_licensed_paired_photo_assets_successfully_retrieved_and_decoded"],4)
            root=Path(td)
            self.assertTrue((root/"echium_blind_pair_01.jpg").is_file())
            self.assertTrue((root/"echium_blind_pair_02.jpg").is_file())
            with (root/"blinded_image_human_QA_BLANK_template.csv").open() as fh:
                rows=list(csv.DictReader(fh))
            self.assertEqual(len(rows),4)
            self.assertTrue(all(not r["open_flower_visible_yes_no_uncertain"] for r in rows))
            self.assertTrue(all(not r["whole_single_plant_visible_yes_no_uncertain"] for r in rows))
            attribution=(root/"image_attribution_and_verified_bytes.csv").read_text()
            self.assertIn("Synthetic test author",attribution)
            self.assertNotIn("internal_original_id",attribution)
            self.assertNotIn("internal_source_100m_cell",attribution)
            self.assertNotIn("positional_accuracy",attribution)
            self.assertNotIn("observed_on",attribution)
            self.assertTrue(result["original_photo_blinded_review_annotations_are_ALL_BLANK"])

if __name__=="__main__":unittest.main()

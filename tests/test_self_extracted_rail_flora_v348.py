#!/usr/bin/env python3
"""Synthetic spatial geometry, observation privacy and non-inference gates."""
import json,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from pilot_self_extracted_six_continent_rail_flora_v348 import (
    bounds,CITIES,NEAR_M,FAR_MIN_M,FAR_MAX_M,MIN_SUPPORT_NEAR,
    eligible_rail,overpass_rail_bytes,osm_xml_rail_bytes,project,
    point_segment_distance,rail_segments,min_rail_distance,classify,
    quality_record,plant_binomial,inat_parameters,MAX_OBSERVATIONS
)

class V348SixContinentTests(unittest.TestCase):
    def test_geography_fixed_before_species(self):
        self.assertEqual(len(CITIES),6)
        self.assertEqual(len({c[1] for c in CITIES}),6)
        self.assertEqual(MAX_OBSERVATIONS,3600)
        self.assertEqual(len(bounds(35.6,139.7)),4)

    def test_strict_rail_tags(self):
        self.assertTrue(eligible_rail({"railway":"rail"}))
        self.assertFalse(eligible_rail({"railway":"subway"}))
        self.assertFalse(eligible_rail({"railway":"rail","tunnel":"yes"}))
        self.assertFalse(eligible_rail({"railway":"rail","bridge":"yes"}))
        self.assertFalse(eligible_rail({"railway":"rail","service":"siding"}))

    def test_original_overpass_way_tag_and_geometry(self):
        payload={"elements":[
            {"type":"way","tags":{"railway":"rail"},"geometry":[
              {"lat":0.0,"lon":0.0},{"lat":0.0,"lon":0.01}]},
            {"type":"way","tags":{"railway":"rail","bridge":"yes"},"geometry":[
              {"lat":0.1,"lon":0.0},{"lat":0.1,"lon":0.01}]}
        ]}
        segments,qa=overpass_rail_bytes(json.dumps(payload).encode())
        self.assertEqual(qa["source_rail_way_count_before_tag_filters"],2)
        self.assertEqual(qa["source_rail_ways_eligible"],1)
        self.assertEqual(len(segments),1)

    def test_original_osm_primary_api_xml(self):
        raw=b'''<?xml version="1.0"?><osm version="0.6">
          <node id="1" lat="0" lon="0"/><node id="2" lat="0" lon=".01"/>
          <way id="5"><nd ref="1"/><nd ref="2"/><tag k="railway" v="rail"/></way>
          <way id="6"><nd ref="1"/><nd ref="2"/><tag k="railway" v="tram"/></way>
          </osm>'''
        ways,qa=osm_xml_rail_bytes(raw)
        self.assertEqual(len(ways),1)
        self.assertEqual(qa["source_rail_segments_eligible"],1)

    def test_point_segment_and_uncertainty_exclusion(self):
        ways=[[(0.0,0.0),(0.0,.01)]]
        seg=rail_segments(ways,0,0)
        near=min_rail_distance(.0003,.005,0,0,seg)
        mid=min_rail_distance(.002,.005,0,0,seg)
        far=min_rail_distance(.004,.005,0,0,seg)
        self.assertEqual(classify(near,20),"near")
        self.assertEqual(classify(mid,20),"ambiguous_buffer")
        self.assertEqual(classify(far,20),"background")
        self.assertEqual(classify(near,300),"ambiguous_buffer")
        self.assertIsNone(min_rail_distance(.004,.005,0,0,[]))

    def test_reject_captive_obscured_uncertain_and_genus(self):
        base={"obscured":False,"geoprivacy":None,"captive":False,
              "quality_grade":"research","positional_accuracy":12,
              "geojson":{"coordinates":[1.0,0.0]},
              "taxon":{"rank":"species","name":"Cynodon dactylon"}}
        x,err=quality_record(base)
        self.assertIsNone(err)
        self.assertEqual(x[0],"cynodon dactylon")
        for value,reason in [
            ({**base,"captive":True},"captive"),
            ({**base,"obscured":True},"obscured_or_geoprivate"),
            ({**base,"positional_accuracy":None},"accuracy_missing"),
            ({**base,"positional_accuracy":200},"accuracy_above_30m_or_invalid"),
            ({**base,"taxon":{"rank":"genus","name":"Cynodon"}},"not_exact_species_rank"),
        ]:
            self.assertEqual(quality_record(value)[1],reason)

    def test_api_parameters_source_bounded(self):
        p=inat_parameters(35.6812,139.7671)
        self.assertEqual(p["quality_grade"],"research")
        self.assertEqual(p["captive"],"false")
        self.assertEqual(p["d1"],"2025-01-01")
        self.assertEqual(p["d2"],"2025-12-31")
        self.assertEqual(p["per_page"],200)

if __name__=="__main__":
    unittest.main()

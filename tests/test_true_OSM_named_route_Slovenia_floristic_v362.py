#!/usr/bin/env python3
"""v362 synthetic true named-route membership, not station-name imputation."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_true_OSM_named_rail_route_floristic_similarity_v362 import (
    route_mode,geometry_only_matches,score_incidence,
    PANELS,MODES,GEO_MAX,GEO_CAL,NETWORK_CAL,SITE_N,
    EXPECTED_SITE_COVERAGE
)
class V362TrueNamedOSMRouteBotanyTests(unittest.TestCase):
    def test_route_is_osm_relation_membership_not_appearance_of_segment_name(self):
        self.assertTrue(route_mode({12},{12},"unique_named_route"))
        self.assertFalse(route_mode({12},{13},"unique_named_route"))
        self.assertIsNone(route_mode({12,13},{12},"unique_named_route"))
        self.assertTrue(route_mode({12,13},{13,14},"any_named_route_overlap"))
        self.assertFalse(route_mode({12,13},{14},"any_named_route_overlap"))
        self.assertIsNone(route_mode(set(),{12},"any_named_route_overlap"))

    def test_real_source_coverage_gate_is_predeclared_43_42_1_3(self):
        self.assertEqual(EXPECTED_SITE_COVERAGE,(43,42,1,3))
        self.assertEqual(PANELS,("unique_named_route","any_named_route_overlap"))
        self.assertEqual(MODES,("all_year_layout","exact_year_and_layout"))
        self.assertEqual((GEO_MAX,GEO_CAL,NETWORK_CAL),(40.,2.5,5.))

    def test_true_route_similarity_does_not_choose_plant_species_by_sign(self):
        sites=[{"lat":46.0,"lon":15.0+i*.00001,
           "year":2022,"layout":"SINGLE","section":"R1",
           "BS":{f"plant{i}"},"BT":{f"plant{i}"}} for i in range(SITE_N)]
        r=[set() for _ in sites]
        r[0]={123};r[1]={123};r[2]={789}
        network={(0,1):1.0,(0,2):1.2,(1,0):1.0,(1,2):1.1,(2,0):1.2,(2,1):1.1}
        triples,qa=geometry_only_matches(sites,network,r,"unique_named_route",False)
        self.assertEqual(len(triples),2)
        self.assertEqual(triples[0]["anchor_index"],0)
        self.assertEqual(triples[0]["within_index"],1)
        self.assertEqual(triples[0]["other_index"],2)
        self.assertLessEqual(triples[0]["abs_geo_gap_km"],GEO_CAL)
        self.assertLessEqual(triples[0]["abs_network_gap_km"],NETWORK_CAL)
        alt=[dict(s,BS={"different","plant"},BT={"different"}) for s in sites]
        a,_=geometry_only_matches(alt,network,r,"unique_named_route",False)
        self.assertEqual(a,triples)

    def test_jaccard_scoring_uses_actual_route_group_equal_weights(self):
        rows=[{"poa","trifolium"},{"poa","bromus"},{"poa"},{"rumex"}]
        t=[{"anchor_index":0,"within_index":1,"other_index":3,
            "original_anchor_route_signature_PRIVATE":(1,)}]
        score=score_incidence(rows,t,[{"section":"author"} for _ in rows],"BS_and_BT")
        self.assertEqual(score["source_n_independent_route_signatures"],1)
        self.assertGreater(score["source_line_equal_same_minus_other_Jaccard"],0)

if __name__=="__main__":unittest.main()

"""v332 synthetic validation of Kashmir geographic matching and thermal decomposition."""
import unittest,sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_geographic_station_between_EIVE_v332 as m

def s(group,lat,lon,alt,effect=0,rich=6,interior_rich=12,cover=.8):
    return {"group":group,"geo":(lat,lon,alt),
            "temp_total":effect,"temp_track":6+effect,"temp_interior":6,
            "track_richness":rich,"interior_richness":interior_rich,
            "track_trait_n":5,"interior_trait_n":10,
            "track_cov":cover,"interior_cov":cover}

class SourceMatching(unittest.TestCase):
    def test_original_dms_curly_minutes_and_spacing(self):
        self.assertAlmostEqual(m.parse_dms("34o  20’ 30.0",True),
                               34+20/60+30/3600)
        self.assertAlmostEqual(m.parse_dms("74 o 30’14.0",False),
                               74+30/60+14/3600)
    def test_impossible_DMS_rejected_not_geocoded(self):
        for x in ("","34.123","34o 60’ 10","34o 23’ 70","77o 13’ 10"):
            with self.assertRaises(ValueError):
                m.parse_dms(x,True)
    def test_haversine_zero_nearby(self):
        a=(34,74,1600);b=(34,74,1700)
        self.assertEqual(m.haversine_km(a,b),0)
        self.assertGreater(m.haversine_km(a,(34.1,74,1600)),10)
    def test_one_to_one_shortest_even_with_station_duplicates(self):
        src={
          "Between A":s("between",34.0,74.0,1600,effect=.5),
          "Between B":s("between",34.01,74.1,1600,effect=.7),
          "Station A":s("station",34.0,74.001,1601,effect=.1),
          "Station B":s("station",34.01,74.101,1600,effect=.2),
          "Station C":s("station",35,75,1770,effect=-.2)}
        pairs,q=m.source_match(src,"G1_geographic_altitude_primary")
        self.assertEqual(len(pairs),2)
        self.assertEqual(q["matched_count"],2)
        self.assertAlmostEqual(sum(p["between"]["temp_total"]-p["station"]["temp_total"] for p in pairs),.9)
        self.assertTrue(all(p["geodesic_km"]<1 for p in pairs))
    def test_tight_altitude_excludes_geographically_close(self):
        src={"Between 1":s("between",34,74,1600),
             "Station 1":s("station",34,74,1670)}
        a,_=m.source_match(src,"G1_geographic_altitude_primary")
        b,_=m.source_match(src,"G2_tight_geographic_altitude")
        self.assertEqual(len(a),1)
        self.assertEqual(len(b),0)
    def test_richness_gate_excludes_unbalanced_even_near_geo(self):
        src={"Between 1":s("between",34,74,1600,rich=5),
             "Station 1":s("station",34,74,1600,rich=12)}
        a,_=m.source_match(src,"G1_geographic_altitude_primary")
        b,_=m.source_match(src,"G3_geo_altitude_richness_coverage")
        self.assertEqual(len(a),1)
        self.assertEqual(len(b),0)
    def test_matched_track_minus_interior_component_additivity(self):
        a=s("between",34,74,1600,effect=.8)
        b=s("station",34,74,1600,effect=.3)
        pairs=[{"between":a,"station":b,"geodesic_km":0,
                "alt_abs_m":0,"alt_signed_between_minus_station_m":0}]
        result=m.summarize_pairs(pairs,{"matched_count":1},
             "G1_geographic_altitude_primary","position",2014,
             "Interior Plot 2","strict_3systems_5species_50pct",{"eligible":2})
        self.assertAlmostEqual(result["matched_between_minus_station"][
            "delta_track_minus_interior"]["mean"],.5)
    def test_signflip_not_confident_for_one_pair(self):
        self.assertIsNone(m.exact_signflip([0.5]))
        self.assertAlmostEqual(m.exact_signflip([1.,2.,3.])[
            "two_sided_exact_signflip_fraction"],.25)

if __name__=="__main__":unittest.main()

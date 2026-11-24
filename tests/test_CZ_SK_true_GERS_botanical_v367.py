#!/usr/bin/env python3
"""v367 synthetic independent source Czech station matching and Slovak decay QA."""
import importlib.util,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Czech35_true_GERS_rail_network_floristic_v367 import (
  matched_geographic_and_physical_rail_triplets,line_equal_gap,
  MAX_GEODESIC_KM,MAX_ABS_GEO_GAP_KM,MAX_ABS_GERS_NETWORK_GAP_KM,
  MIN_MATCHED_STATION_ANCHORS,MIN_AUTHOR_STRETCHES,EXPECTED_GERS_CONNECTED_STATION_PAIRS
)
from analyze_Slovakia19_true_GERS_rail_network_beta_v367 import (
  EXPECTED_GERS_CONNECTED_PAIRS as SK_EXPECT,
  SK_SOURCE_SITE_N,SK_SOURCE_TAXA_N,SK_SOURCE_SITE_TAXON_INCIDENCES,
  LAMBDA,DRAWS,BURN,THIN,station_pairs
)

class V367IndependentCzechSlovakTests(unittest.TestCase):
    def test_original_country_sample_grain_not_mixed(self):
        self.assertEqual(EXPECTED_GERS_CONNECTED_STATION_PAIRS,568)
        self.assertEqual((MAX_GEODESIC_KM,MAX_ABS_GEO_GAP_KM,
                          MAX_ABS_GERS_NETWORK_GAP_KM),(160.,10.,20.))
        self.assertEqual((MIN_MATCHED_STATION_ANCHORS,MIN_AUTHOR_STRETCHES),(12,4))
        self.assertEqual((SK_SOURCE_SITE_N,SK_SOURCE_TAXA_N,SK_SOURCE_SITE_TAXON_INCIDENCES),(19,318,1763))
        self.assertEqual(SK_EXPECT,171)
        self.assertEqual((LAMBDA,DRAWS,BURN,THIN),(.05,199,18000,650))

    def test_Slovak_real_GERS_distance_dictionary_pair_key_not_station_ordinal(self):
        # If code erroneously calls network.get(i,j), the second station
        # index j silently becomes a numeric *network distance*. Fix it
        # with an explicit regression using an obviously different source
        # physical railway route metric.
        sites=[{"lon":20.+i*.001,"lat":49.3} for i in range(19)]
        flora={i+1:{"poa",f"source{i}"} for i in range(19)}
        net={(i,j):1.0+abs(i-j)*.2 for i in range(19) for j in range(19) if i!=j}
        records=station_pairs(sites,flora,net)
        self.assertEqual(len(records),171)
        p0=next(x for x in records if x["i"]==0 and x["j"]==1)
        self.assertAlmostEqual(p0["network"],1.2)
        self.assertNotEqual(p0["network"],p0["j"])
        p1=next(x for x in records if x["i"]==0 and x["j"]==18)
        self.assertAlmostEqual(p1["network"],4.6)
        self.assertNotEqual(p1["network"],p1["j"])
        net2=dict(net);net2[0,1]=6.0
        with self.assertRaises(ValueError):
            station_pairs(sites,flora,net2)

    def test_original_Czech_must_have_real_same_and_other_routes_and_GERS_paths(self):
        sites=[]
        for i in range(35):
            sites.append({"lon":14.0+(i%7)*.08,"lat":49.8+(i//7)*.06,
                          "source_stretch":str(i//7),"taxa":{"a","b"}})
        rail={(i,j):max(1.0,abs(i-j)*3.0) for i in range(35) for j in range(35) if i!=j}
        matched,qa=matched_geographic_and_physical_rail_triplets(sites,rail)
        for t in matched:
            self.assertNotEqual(t["same"],t["other"])
            self.assertEqual(sites[t["anchor"]]["source_stretch"],sites[t["same"]]["source_stretch"])
            self.assertNotEqual(sites[t["anchor"]]["source_stretch"],sites[t["other"]]["source_stretch"])
            self.assertLessEqual(t["abs_geo_gap"],MAX_ABS_GEO_GAP_KM)
            self.assertLessEqual(t["abs_network_gap"],MAX_ABS_GERS_NETWORK_GAP_KM)
        with no_source_species_effect():
            other=[dict(x,taxa={"completely","new"}) for x in sites]
        identical,_=matched_geographic_and_physical_rail_triplets(other,rail)
        self.assertEqual(identical,matched)

    def test_no_rail_network_paths_not_zero_effect(self):
        sites=[{"lon":14.+i*.001,"lat":50.,"source_stretch":str(i//7),
                "taxa":{"poa","trifolium"}} for i in range(35)]
        triads,qa=matched_geographic_and_physical_rail_triplets(sites,{})
        self.assertEqual(triads,[])
        self.assertGreater(qa["source_true_GERS_path_missing"],0)

    def test_Czech_multiple_original_stretch_Jaccard_unweighted(self):
        sites=[{"source_stretch":str(i//2)} for i in range(8)]
        rows=[{"a","c"},{"a"},{"b","c"},{"b"},{"d","c"},{"d"},{"e","c"},{"e"}]
        matches=[
            {"anchor":i*2,"same":i*2+1,"other":((i+1)%4)*2,
             "original_source_stretch_PRIVATE":str(i)} for i in range(4)]
        score,line=line_equal_gap(rows,matches,sites)
        self.assertEqual(len(line),4)
        self.assertGreaterEqual(score,-1.)
        self.assertLessEqual(score,1.)

def no_source_species_effect():
    from contextlib import nullcontext
    return nullcontext()

if __name__=="__main__":unittest.main()

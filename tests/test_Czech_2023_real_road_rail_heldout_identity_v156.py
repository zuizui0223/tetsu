"""v15.6 SYNTHETIC source statistic tests; no synthetic outcomes presented as botanical result."""
import sys,math,unittest
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_Czech_2024_centroid_geographic_identity_v154 import fake_source_sites
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from analyze_Czech_2023_true_rail_road_heldout_taxon_identity_v156 import (
    pair_records,cv_section_heldout,rail_coefficient,ALPHA,MIN_STRETCHES
)

def synthetic_pairs(sites,regular=True):
    out=[]
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
        road=geo*1.2+0.9
        rail=road*math.exp(0.1*math.sin(i*17+j*5)) if regular else 100000
        out.append({
           "i":i,"j":j,
           "source_same_stretch":a["section"]==b["section"],
           "source_station_geodesic_km":geo,
           "network_route_comparison_eligible":True,
           "rail_graph_source_shortest_path_km":rail,
           "road_graph_source_shortest_path_km":road,
           "source_both_rail_snap_km":0.05,
           "source_both_road_snap_km":0.2,
        })
    return out

class CzechRoadConditionalModelTests(unittest.TestCase):
    def setUp(self):
        self.sites,self.guild=fake_source_sites()

    def test_original_5_stretch_leaveout_on_source_station_pairs(self):
        original=synthetic_pairs(self.sites)
        source,excluded=pair_records(original,self.sites,self.guild)
        self.assertTrue(250<=len(source)<=595)
        r=cv_section_heldout(source)
        self.assertTrue(r["identified"])
        self.assertEqual(r["n_holdout_original_stretches"],5)
        self.assertIn("M2_plus_actual_2023_road_and_rail_ratio",
                       r["equal_original_stretch_weighted_mse"])
        self.assertLessEqual(r["rail_improves_how_many_stretches"],5)
        self.assertEqual(r["ridge_alpha_fixed"],ALPHA)

    def test_all_author_observed_species_is_alternate_not_imputed_trait(self):
        original=synthetic_pairs(self.sites)
        a,_=pair_records(original,self.sites,self.guild,"both_Pladias_known")
        b,_=pair_records(original,self.sites,self.guild,"all_original_S1_current")
        self.assertEqual(len(a),len(b))
        self.assertTrue(any(x["jaccard"]!=y["jaccard"] for x,y in zip(a,b)))

    def test_unconnected_or_boundary_routes_not_treated_as_biological_zeros(self):
        original=synthetic_pairs(self.sites)
        original[0]["network_route_comparison_eligible"]=False
        a,omitted=pair_records(original,self.sites,self.guild)
        self.assertEqual(omitted.get("no_geography_and_both_true_rail_road_route_support"),1)
        self.assertLess(len(a),595)
        too_long=synthetic_pairs(self.sites,regular=False)
        b,omitted2=pair_records(too_long,self.sites,self.guild)
        self.assertEqual(len(b),0)
        self.assertGreater(omitted2.get("predeclared_extreme_route_detour",0),0)

    def test_conditional_rail_coef_descriptive_not_causal(self):
        original=synthetic_pairs(self.sites)
        source,_=pair_records(original,self.sites,self.guild)
        r=rail_coefficient(source)
        self.assertTrue(r["identified"])
        self.assertIn("no_plant_seed_carriage_or_selection_causality_from_coefficient",r)
        self.assertTrue(math.isfinite(r["source_rail_over_road_path_logratio_partial_beta"]))

    def test_no_common_source_support_honestly_not_estimable(self):
        orig=synthetic_pairs(self.sites)
        for p in orig:
            if p["i"]!=0:p["network_route_comparison_eligible"]=False
        source,_=pair_records(orig,self.sites,self.guild)
        self.assertFalse(cv_section_heldout(source)["identified"])
        self.assertEqual(MIN_STRETCHES,4)

if __name__=="__main__":unittest.main()

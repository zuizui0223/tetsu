"""v15.1 SYNTHETIC correctness tests only; actual biology is source-authenticated in CI."""
import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_road_heldout_species_identity_v151 import (
    model_records,supports_estimator,section_heldout_validation,
    rail_partial_slope,BASE,ridge_train_predict,
)

def artificial_opportunity(seed=7,change_rail=True):
    rng=random.Random(seed)
    # 60 fake sites, 10 fake route sections, 1770 fake pairs; NOT
    # experimental or genuine Slovenian field data.
    sites=[{"section":str(i%10),"year":2022+i%2,
            "track_type":"single" if i%3 else "double",
            "known":set(range((i%6)+10))} for i in range(60)]
    pairs=[]
    for i in range(len(sites)):
        for j in range(i+1,len(sites)):
            geod=2.1+.5*abs(i-j)
            road=geod*(1.10+.03*(i%7)+.06*abs(math.sin(i+j)))
            log_ratio=.15*math.sin(i*13+j*7)+.09*math.cos(i*j+7) if change_rail else 0.0
            rail=road*math.exp(log_ratio)
            sim=.7+.25*rng.random()
            source_y=.5-.18*log_ratio+.02*sim
            pairs.append({
              "i":i,"j":j,"same_original_section":sites[i]["section"]==sites[j]["section"],
              "geo_km":geod,"geographic_bin":min(6,int(math.log(geod+1,2))),
              "functional_bin":3,"functional_overlap":sim,
              "same_original_year":sites[i]["year"]==sites[j]["year"],
              "same_original_track_type":sites[i]["track_type"]==sites[j]["track_type"],
              "source_original_jaccard":source_y,
              "rail_path_km":rail,"road_path_km":road,
              "source_rail_both_site_access_km":.03+.01*(i%4+j%3),
              "source_road_both_site_access_km":.3+.02*(i%7+j%5),
              "supported_BOTH_rail_road_routes_and_geo_ge_2km":True
            })
    return pairs,sites

class OriginalRailRoadDesignMath(unittest.TestCase):
    def test_fixed_historical_rail_plus_road_features_and_outcome(self):
        pairs,sites=artificial_opportunity()
        rows,omitted=model_records(pairs,sites)
        self.assertEqual(len(rows),1770)
        self.assertFalse(omitted)
        self.assertIn("same_original_section",rows[0]["features"])
        self.assertIn("road_path_ln",rows[0]["features"])
        self.assertIn("relative_rail_road_path_ln",rows[0]["features"])
        self.assertTrue(supports_estimator(rows)["supports_main_regression"])

    def test_heldout_section_cv_can_recover_incremental_source_rail_advantage(self):
        pairs,sites=artificial_opportunity()
        rows,_=model_records(pairs,sites)
        cv=section_heldout_validation(rows)
        self.assertTrue(cv["identified"])
        self.assertEqual(cv["n_heldout_sections_estimable"],10)
        self.assertGreater(cv["n_source_sections_rail_increment_positive"],6)
        self.assertGreater(cv["source_section_mean_incremental_rail_after_road_MSE_reduction"],0)

    def test_conditional_ln_route_slope_sign_known_in_synthetic_data(self):
        pairs,sites=artificial_opportunity()
        records,_=model_records(pairs,sites)
        d=rail_partial_slope(records)
        self.assertTrue(d["identified"])
        self.assertLess(d["partial_beta_jaccard_per_1unit_ln_rail_over_road"],-.05)

    def test_collinear_rail_road_path_ratio_identified_false(self):
        pairs,sites=artificial_opportunity(change_rail=False)
        records,_=model_records(pairs,sites)
        self.assertFalse(rail_partial_slope(records)["identified"])

    def test_no_fictitious_jaccard_if_road_unreachable(self):
        pairs,sites=artificial_opportunity()
        pairs[0]["road_path_km"]=None
        pairs[0]["supported_BOTH_rail_road_routes_and_geo_ge_2km"]=False
        records,omitted=model_records(pairs,sites)
        self.assertEqual(len(records),1769)
        self.assertEqual(omitted["road_rail_missing_or_geo_under_2km"],1)

    def test_detour_trim_is_predefined_not_posthoc(self):
        pairs,sites=artificial_opportunity()
        pairs[0]["rail_path_km"]=10000.
        rec_trim,omitted=model_records(pairs,sites,True)
        rec_all,_=model_records(pairs,sites,False)
        self.assertEqual(len(rec_trim),len(rec_all)-1)
        self.assertEqual(omitted["extreme_detour_or_border_truncation"],1)

if __name__=="__main__":unittest.main()

"""v15.4 source spatial matching tests. All fixtures are SYNTHETIC math only."""
import sys,math,unittest
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    real_site_pairs,source_gap,conditional_source_gap,
    fixed_source_margin_null,source_gap,section_leaveout,
    MIN_SAME_AND_DIFF_PAIRS_EACH_STRATUM
)

def fake_source_sites():
    guild={**{f"w{n}":"th_wind" for n in range(7)},
           **{f"o{n}":"nonth_other" for n in range(7)}}
    sites=[]
    for i in range(35):
        known={f"w{(i+j)%7}" for j in (0,1)}|{
              f"o{(i*2+j)%7}" for j in (0,1)}
        sites.append({
            "site":f"fake{i}","section":str(i%5),
            "known":known,
            "all_original_taxa":known|{f"extra{i}"},
            "area_ha":float(i%7+1),
            "lat":49.+.03*(i//5),
            "lon":15.+.035*(i%5)
        })
    return sites,guild

class CzechGeoMatchingMathTests(unittest.TestCase):
    def test_no_source_595_pairs_treated_as_independent_experiments(self):
        sites,guild=fake_source_sites()
        p=real_site_pairs(sites,guild)
        self.assertEqual(len(p),595)
        self.assertEqual(sum(x["same_stretch"] for x in p),5*(7*6//2))
        self.assertEqual(sum(not x["same_stretch"] for x in p),595-105)
        self.assertTrue(source_gap(p)["same_minus_different_jaccard"] is not None)

    def test_matching_requires_actual_same_and_different_in_same_geographic_bins(self):
        pairs=[
            {"i":0,"j":1,"same_stretch":True,"geo_bin":1,"function_bin":3,
             "area_log_diff_bin":0,"taxonomic_Jaccard_original_544":.6},
            {"i":0,"j":2,"same_stretch":False,"geo_bin":1,"function_bin":3,
             "area_log_diff_bin":0,"taxonomic_Jaccard_original_544":.4},
            {"i":0,"j":3,"same_stretch":True,"geo_bin":4,"function_bin":3,
             "area_log_diff_bin":0,"taxonomic_Jaccard_original_544":.7}]
        x=conditional_source_gap(pairs,min_each=1)
        self.assertTrue(x["identified"])
        self.assertAlmostEqual(x["weighted_same_minus_different_jaccard"],.2)
        self.assertEqual(x["source_support"]["within_source_stretch_pairs_on_common_support"],1)
        self.assertEqual(x["source_support"]["different_source_stretch_pairs_on_common_support"],1)
        z=conditional_source_gap(pairs,min_each=2)
        self.assertFalse(z["identified"])
        self.assertIsNone(z["weighted_same_minus_different_jaccard"])

    def test_area_missing_is_unknown_not_equal_area(self):
        pairs=[
            {"i":0,"j":1,"same_stretch":True,"geo_bin":1,"function_bin":2,
             "area_log_diff_bin":None,"taxonomic_Jaccard_original_544":.5},
            {"i":1,"j":2,"same_stretch":False,"geo_bin":1,"function_bin":2,
             "area_log_diff_bin":0,"taxonomic_Jaccard_original_544":.4}]
        s=conditional_source_gap(pairs,mode="geography_function_and_area",min_each=1)
        self.assertFalse(s["identified"])
        self.assertEqual(s["source_support"]["excluded_missing_area_or_outside_original_declared_bins"],1)

    def test_global_source_fixed_margins_shuffle_computes_without_inventing_source_gps(self):
        sites,guild=fake_source_sites()
        pairs=[]
        for i,j in combinations(range(len(sites)),2):
            a,b=sites[i],sites[j]
            pairs.append({"i":i,"j":j,"same_stretch":a["section"]==b["section"],
                          "geo_bin":2,"function_bin":3,"area_log_diff_bin":0,
                          "taxonomic_Jaccard_original_544":len(a["known"]&b["known"])/len(a["known"]|b["known"])})
        result=fixed_source_margin_null(sites,guild,pairs,nchains=1,draws=5,
                                        burnin=100,thin=200)
        self.assertTrue(result["identified"])
        self.assertEqual(result["n_correlated_source_mcmc_draws"],5)
        self.assertEqual(result["source_chain_diagnostics"][0]["proposed"],1100)
        strict=fixed_source_margin_null(sites,guild,pairs,
              within_stretch=True,nchains=1,draws=3,burnin=100,thin=200)
        self.assertTrue(strict["within_original_stretch_species_prevalence_additionally_fixed"])

    def test_non_existing_geographic_common_support_fails_closed(self):
        sites,guild=fake_source_sites()
        p=real_site_pairs(sites,guild)
        for item in p:
            item["geo_bin"]=0 if item["same_stretch"] else 6
        m=conditional_source_gap(p)
        self.assertFalse(m["identified"])
        self.assertIsNone(m["weighted_same_minus_different_jaccard"])

if __name__=="__main__":unittest.main()

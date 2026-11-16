"""v14.2 synthetic math-only spatial matching / source common support tests.

The genuine ecological result is computed ONLY in the original actual
Slovenia 89-site author XLSX CI workflow with source SHA pin.
"""
from __future__ import annotations
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    DISTANCE_CUTS_KM,FUNCTIONAL_OVERLAP_CUTS,
    haversine_km,into_bin,spatial_common_support,spatial_null,source_geo_confounds,
)

G={"A":"strict_annual_wind","B":"strict_annual_wind",
   "C":"strict_annual_wind","D":"strict_annual_wind",
   "E":"other_known_other","F":"other_known_other"}
def site(label,sec,sp,year=2022):
    return {"site":label,"section":sec,"year":year,"track_type":"SINGLE",
            "known":set(sp),"all_original_taxa":set(sp),
            "n_original_matched_track_sides":1}
S=[
    site("one","a","ABE",2022),site("two","a","ABE",2023),
    site("three","b","CDF",2022),site("four","b","CDF",2023),
    site("five","c","ACE",2022),site("six","c","BDF",2023),
]


class RealSourceSpatialConditionalityTests(unittest.TestCase):
    def test_actual_haversine_source_WGS84_coordinate_scale(self):
        self.assertAlmostEqual(haversine_km(46,15,46,15),0.,places=9)
        self.assertTrue(100<haversine_km(46,15,47,15)<120)
        self.assertAlmostEqual(haversine_km(46,15,46,16),
                               haversine_km(46,16,46,15),places=12)

    def test_predeclared_geography_and_trait_bins_are_not_posthoc_tuned(self):
        self.assertEqual(DISTANCE_CUTS_KM,(0,5,10,20,40,80,160,500))
        self.assertEqual(FUNCTIONAL_OVERLAP_CUTS,(0,.7,.8,.9,1.000001))
        self.assertEqual(into_bin(0,DISTANCE_CUTS_KM),0)
        self.assertEqual(into_bin(5,DISTANCE_CUTS_KM),1)
        self.assertEqual(into_bin(40,DISTANCE_CUTS_KM),4)
        self.assertEqual(into_bin(.9,FUNCTIONAL_OVERLAP_CUTS),3)
        self.assertIsNone(into_bin(700,DISTANCE_CUTS_KM))

    def test_geo_matched_gaps_refuse_to_extrapolate_if_no_common_pair_bins(self):
        r=[
          {"i":0,"j":1,"same":True,"geo_bin":0,"functional_bin":3,
           "same_year":True,"same_track_type":True,
           "observed_known_species_shared_jaccard":.8},
          {"i":0,"j":2,"same":False,"geo_bin":3,"functional_bin":3,
           "same_year":True,"same_track_type":True,
           "observed_known_species_shared_jaccard":.2},
        ]
        v=spatial_common_support(r,mode="geographic_plus_functional",min_each=1)
        self.assertFalse(v["identified"])
        self.assertIsNone(v["weighted_source_identity_gap"])
        self.assertEqual(v["source_support"]["n_same_section_pairs_on_identical_source_support"],0)

    def test_geo_functional_matching_restricts_to_actual_same_bin_site_pairs_only(self):
        r=[]
        for i in range(3):
            r.append({"i":0,"j":1,"same":True,"geo_bin":2,
                      "functional_bin":3,"same_year":True,
                      "same_track_type":True,
                      "observed_known_species_shared_jaccard":.70+i*.01})
            r.append({"i":1,"j":2,"same":False,"geo_bin":2,
                      "functional_bin":3,"same_year":True,
                      "same_track_type":True,
                      "observed_known_species_shared_jaccard":.20+i*.01})
        r.append({"i":1,"j":3,"same":True,"geo_bin":0,
                  "functional_bin":3,"same_year":True,
                  "same_track_type":True,
                  "observed_known_species_shared_jaccard":1.})
        ans=spatial_common_support(r,mode="geographic_plus_functional")
        self.assertTrue(ans["identified"])
        self.assertAlmostEqual(ans["weighted_source_identity_gap"],.5)
        self.assertEqual(ans["source_support"]["n_same_section_pairs_on_identical_source_support"],3)
        self.assertEqual(ans["source_support"]["n_different_section_pairs_on_identical_source_support"],3)
        self.assertEqual(ans["source_support"]["n_total_same_section_pairs"],4)

    def test_year_and_track_strict_support_may_be_unidentified(self):
        r=[
            {"i":i,"j":i+1,"same":bool(i%2==0),
             "geo_bin":1,"functional_bin":3,
             "same_year":bool(i%2==0),
             "same_track_type":True,
             "observed_known_species_shared_jaccard":.6} for i in range(6)
        ]
        a=spatial_common_support(r,"",mode="geographic_plus_functional") if False else spatial_common_support(
            r,mode="geographic_plus_functional",min_each=2)
        self.assertTrue(a["identified"])
        b=spatial_common_support(r,mode="geographic_functional_same_year_track",min_each=2)
        self.assertFalse(b["identified"])

    def test_synthetic_guild_matched_randomization_has_two_modes_no_fake_geo_data(self):
        # Only a math fixture: not real 89-site GPS or rail distances.
        pairs=[]
        from itertools import combinations
        from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import functional_overlap,jaccard
        for i,j in combinations(range(6),2):
            pairs.append({"i":i,"j":j,"same":S[i]["section"]==S[j]["section"],
                          "geo_bin":0,"functional_bin":
                              into_bin(functional_overlap(S[i]["known"],S[j]["known"],G),
                                       FUNCTIONAL_OVERLAP_CUTS),
                          "same_year":S[i]["year"]==S[j]["year"],
                          "same_track_type":True,
                          "observed_known_species_shared_jaccard":
                              jaccard(S[i]["known"],S[j]["known"])})
        support=spatial_common_support(pairs,mode="geographic_plus_functional",min_each=1)
        self.assertTrue(support["identified"])
        # Production function uses default min_each=3: potential strata
        # insufficient is correctly reported non-estimable.
        actual=list(spatial_null(S,G,pairs,"geographic_plus_functional",
                                draws_per_chain=20,nchains=1,burnin=100,spacing=100))
        self.assertEqual({k for k,v in actual},{"source_guild_species_degree_fixed_ALL_STUDY",
                                          "source_guild_species_degree_fixed_WITHIN_YEAR_TRACK_TYPE"})
        self.assertTrue(all(v["n_correlated_chain_draws"]==20 for k,v in actual))


if __name__=="__main__":unittest.main()

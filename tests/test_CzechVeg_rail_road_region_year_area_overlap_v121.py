"""v12.1 original CzechVeg region/year/area exact common-support tests.

Synthetic fixtures test ONLY arithmetic/coverage, never botanical evidence.
Full actual 117739-plot, 2510340-row vegetation source is used in Actions.
"""
from __future__ import annotations
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    year_from_source,area_key,block_key,matched_blocks,
    block_standardized_guild,diagnostic,run,MATCHING_SPECS
)


def source(plot,cat,year,area,region="18a"):
    return {"strict":cat,"broad":cat,"source_year_int":year,
            "year":str(year),"area":area,"source_phyto_region":region,
            "region":region}


def survey_plants(codes):
    return [(c,"Epilobium" if c.startswith("W") else "Allium",
              "therophyte" if c.endswith("A") else "hemicryptophyte")
            for c in codes]


class CzechRealSourceCommonSupportMathTests(unittest.TestCase):
    def test_original_source_year_does_not_impute_year_from_database(self):
        self.assertEqual(year_from_source("2022"),2022)
        self.assertEqual(year_from_source("1960.0"),1960)
        for x in ("","NA","unknown","1000","2040","2022-2023",None):
            self.assertIsNone(year_from_source(x))

    def test_source_area_matching_bins_and_exact_areas(self):
        values=[1,4,5,16,17,64,65,256,257,1000]
        bins=[area_key(x,"class4") for x in values]
        self.assertEqual(bins,["1_to_4","1_to_4","gt4_to_16",
                      "gt4_to_16","gt16_to_64","gt16_to_64",
                      "gt64_to_256","gt64_to_256","gt256_to_1000","gt256_to_1000"])
        self.assertEqual(area_key(16.,"exact"),16.)
        with self.assertRaises(ValueError):area_key(1001,"class4")

    def test_block_sampling_only_uses_source_area_year_region(self):
        a={"a":source("a","rail",2012,16),
           "b":source("b","road",2011,16),
           "c":source("c","rail",2001,16),
           "d":source("d","road",2011,32),
           "e":source("e","rail",2011,16,"68"),
           "f":source("f","road",2011,16,"68")}
        z={x:survey_plants(["W1A"]) for x in a}
        specs=[("decade_exact",10,"exact",1),
               ("20year_class4",20,"class4",1),
               ("region_area_no_year",0,"class4",1)]
        blocks,info=matched_blocks(a,z,"strict",specs[0])
        self.assertEqual(len(blocks),2)
        self.assertEqual(info["n_source_rail_plots_on_matched_common_support"],2)
        self.assertEqual(info["n_source_road_plots_on_matched_common_support"],2)
        matched2,_=matched_blocks(a,z,"strict",specs[2])
        self.assertEqual(len(matched2),2)
        self.assertEqual(sum(len(v["rail"]) for v in matched2.values()),3)

    def test_exact_matched_common_support_standardizes_species_score(self):
        meta={}
        z={}
        for region in ("18a","68"):
            for side in ("rail","road"):
                for i in range(3):
                    plot=f"{region}-{side}-{i}"
                    meta[plot]=source(plot,side,2001,16,region)
                    codes=[]
                    for lifespan in ("A","P"):
                        for k in range(3):
                            codes.append(f"W{k}{lifespan}" if side=="rail"
                                         else f"N{k}{lifespan}")
                    z[plot]=survey_plants(codes)
        spec=("Z",10,"exact",2)
        blocks,diag=matched_blocks(meta,z,"strict",spec)
        self.assertEqual(diag["n_original_source_blocks_retaining_both_arms"],2)
        self.assertEqual(diag["n_source_rail_plots_on_matched_common_support"],6)
        m=block_standardized_guild(blocks,z)
        self.assertTrue(m["identified"])
        self.assertEqual(m["n_original_joint_trait_species_source_common_support"],12)
        self.assertAlmostEqual(m["original_source_species_equal_strict_wind_vs_other_prevalence_delta"],2.)
        self.assertAlmostEqual(m["original_source_occ_weight_strict_wind_vs_other_prevalence_delta"],2.)
        self.assertFalse(m["opposite_source_estimand_signs"])

    def test_real_source_text_groups_with_no_overlap_are_unidentifiable(self):
        x=block_standardized_guild({}, {})
        self.assertFalse(x["identified"])
        self.assertIn("No observed",x["cause"])

    def test_missing_plots_are_not_invented_zero_plant_absences(self):
        meta={str(i):source(i,"rail" if i<3 else "road",2011,16) for i in range(6)}
        source_taxa={str(i):survey_plants(["W1A"]) for i in (0,1,3,4)}
        blocks,diag=matched_blocks(meta,source_taxa,"strict",("x",10,"exact",2))
        self.assertEqual(diag["candidate_plots_with_no_source_vascular_records_excluded"],2)
        self.assertEqual(diag["n_original_source_blocks_retaining_both_arms"],1)

    def test_forged_source_v119_stops_before_independent_country_biology(self):
        with self.assertRaisesRegex(ValueError,"not present"):
            run({}, "headers","species","disp","life",{"complete":False})


if __name__=="__main__":unittest.main()

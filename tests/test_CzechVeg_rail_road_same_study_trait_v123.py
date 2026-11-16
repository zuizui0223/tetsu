"""v12.3: same-study site sample alignment and non-estimable source support."""
from __future__ import annotations
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_rail_road_same_study_trait_v123 import (
    source_study_block_key,common_original_study_blocks,comparison,run,SPECS
)


def plot(site,cl,ref,reg="18a",area=16,year=2001):
    return {
        "original_reference_title":ref,"source_year_int":year,
        "source_phyto_region":reg,"area":area,
        "strict":cl,"broad":cl,
    }


class RealStudySourceControls(unittest.TestCase):
    def test_no_study_confounding_is_silently_merged(self):
        sites={
            "a":plot("a","rail","Published vegetation study A"),
            "b":plot("b","road","Published vegetation study B"),
            "c":plot("c","rail","Published vegetation study A"),
            "d":plot("d","road","Published vegetation study B")
        }
        taxa={s:[("Abies alba","Epilobium","hemicryptophyte")] for s in sites}
        blk,stat=common_original_study_blocks(sites,taxa,"strict",
              ("source",10,"class4",1,True))
        self.assertEqual(len(blk),0)
        self.assertEqual(stat["n_source_RAIL_plots_in_common_support"],0)
        self.assertEqual(stat["n_source_ROAD_plots_in_common_support"],0)

    def test_same_study_region_decade_area_stratification_preserves_source_units(self):
        sites={
            "a":plot("a","rail","Published vegetation study A"),
            "b":plot("b","road","Published vegetation study A"),
            "c":plot("c","rail","Published vegetation study A",area=100),
            "d":plot("d","road","Published vegetation study A",year=1990),
            "e":plot("e","road","Published vegetation study A",reg="68"),
        }
        taxa={s:[("Abies alba","Epilobium","hemicryptophyte")] for s in sites}
        strict,stat=common_original_study_blocks(sites,taxa,"strict",
             ("source",10,"class4",1,True))
        self.assertEqual(len(strict),1)
        self.assertEqual(stat["n_source_RAIL_plots_in_common_support"],1)
        self.assertEqual(stat["n_source_ROAD_plots_in_common_support"],1)
        self.assertEqual(stat["effective_balanced_source_site_pair_count"],1)

    def test_unmatched_cohort_is_nonestimable_not_false_zero(self):
        sites={"r1":plot("r1","rail","Study reference name"),
               "d1":plot("d1","road","Different reference name")}
        taxa={x:[("Abies alba","Epilobium","hemicryptophyte")] for x in sites}
        obj=comparison(sites,taxa,"strict",("S",10,"class4",1,True))
        self.assertFalse(obj["source_same_publication_geo_year_area_standardized"]["identified"])
        self.assertIsNone(obj["pooled_on_same_actual_selected_study_plots"]["equal"])

    def test_min_source_sites_original_study_gate(self):
        sites={str(i):plot(str(i),"rail" if i<5 else "road",
                          "Shared paper sample 2001") for i in range(10)}
        taxa={x:[("Abies alba","Epilobium","hemicryptophyte")] for x in sites}
        for n,expected in [(1,1),(3,1),(5,1),(6,0)]:
            blocks,stats=common_original_study_blocks(
                sites,taxa,"strict",("test",10,"class4",n,True))
            self.assertEqual(stats["n_original_study_blocks_with_rail_and_road_common_support"],expected)

    def test_study_without_geography_is_distinct_from_within_region_design(self):
        a=plot("a","rail","The same report about vegetation",reg="18a")
        b=plot("b","road","The same report about vegetation",reg="18b")
        self.assertNotEqual(
            source_study_block_key(a,("x",10,"class4",1,True)),
            source_study_block_key(b,("x",10,"class4",1,True)))
        self.assertEqual(
            source_study_block_key(a,("x",0,"none",1,False)),
            source_study_block_key(b,("x",0,"none",1,False)))

    def test_original_Czechveg_v119_verified_is_mandatory(self):
        with self.assertRaisesRegex(ValueError,"not present"):
            run({},"h","s","d","l",{"complete":False})


if __name__=="__main__":unittest.main()

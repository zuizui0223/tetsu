"""v12.5 source-verified Czech vegetation same-syntaxon matching and study resample."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_syntaxon_matched_rail_road_traits_v125 import (
    source_blocks,extract_estimate,source_study_cluster_resample,
    source_leave_one_publication_out,panel,run
)


def plant_site(label,syntax="XB",study="Original paper on vegetation",region="18a",year=2001,area=16):
    return {
        "strict":label,"broad":label,
        "study":study,"source_region":region,
        "survey_year":year,"area":area,
        "SYNTAXON":syntax,"SYNTAXON_NAME":"Stellarietea mediae" if syntax else None,
        "SYNTAX_OLD":"39A" if syntax else None,"SYNTAXON_OLD_NAME":None,
    }


def demo():
    sites={}
    plants={}
    for i,paper in enumerate(("Published study A","Published study B","Published study C")):
        for arm in ("rail","road"):
            for sample in range(3):
                name=f"{paper}_{arm}_{sample}"
                sites[name]=plant_site(arm,study=paper)
                vals=[]
                for annual in (True,False):
                    for wind in (True,False):
                        for k in range(4):
                            taxon=f"taxon{annual}{wind}{k}"
                            if (sample+i+k)%4!=3:
                                vals.append((taxon,"Epilobium" if wind else "Allium",
                                    "therophyte" if annual else "hemicryptophyte"))
                plants[name]=vals
    return sites,plants


class CzechOriginalVegetationMatchedSourceTests(unittest.TestCase):
    def test_source_study_area_year_vegetation_gate(self):
        sites,plants=demo()
        blocks,stat=source_blocks(sites,plants,"strict",
          "same_paper_phyto_decade_area_SYNTAXON",2)
        self.assertEqual(stat["n_source_blocks_with_original_rail_and_road"],3)
        self.assertEqual(stat["n_published_source_studies_in_matched_blocks"],3)
        self.assertEqual(stat["n_source_rail_plots_matched"],9)
        self.assertEqual(stat["n_source_road_plots_matched"],9)
        self.assertEqual(stat["balanced_original_study_block_mass"],9)
        self.assertEqual(len(blocks),3)

    def test_original_mismatched_syntaxon_cannot_be_coded_as_same_community(self):
        sites,plants=demo()
        for id,m in sites.items():
            if m["strict"]=="road":
                m["SYNTAXON"]="XC"
        a,x=source_blocks(sites,plants,"strict",
          "same_paper_phyto_decade_area_SYNTAXON",1)
        self.assertFalse(a)
        self.assertEqual(x["n_source_rail_plots_matched"],0)
        b,y=source_blocks(sites,plants,"strict",
          "same_paper_phyto_decade_area_no_syntax",2)
        self.assertEqual(y["n_source_rail_plots_matched"],9)

    def test_absent_source_syntaxon_not_invented_or_imputed(self):
        sites,plants=demo()
        for id,m in sites.items():
            if m["strict"]=="rail":
                m["SYNTAXON"]=None
        a,x=source_blocks(sites,plants,"broad",
          "same_paper_phyto_decade_area_SYNTAXON",1)
        self.assertEqual(len(a),0)
        self.assertEqual(x["n_original_source_rail_plots_missing_syntaxon"],9)

    def test_real_source_survivor_trait_comparison_defined_when_guild_counts_sufficient(self):
        sites,plants=demo()
        result=panel(sites,plants,"strict",
            "same_paper_phyto_decade_area_SYNTAXON",2)
        self.assertEqual(result["source_original_sampling_frame"][
            "n_published_source_studies_in_matched_blocks"],3)
        self.assertTrue(result[
            "source_original_study_area_year_and_syntaxon_standardized"]["identified"])
        self.assertEqual(result[
            "source_original_study_area_year_and_syntaxon_standardized"]["n_source_taxa"],16)
        self.assertTrue("source_study_resampling_conditional" in result)

    def test_original_study_loo_never_pseudoreplicates_source_papers(self):
        sites,plants=demo()
        b,_=source_blocks(sites,plants,"strict",
             "same_paper_phyto_decade_area_SYNTAXON",2)
        loo=source_leave_one_publication_out(b,plants)
        self.assertEqual(len(loo),3)
        self.assertTrue(all(x["remaining_source_studies"]==2 for x in loo))

    def test_forged_original_source_receipts_are_rejected(self):
        with self.assertRaisesRegex(ValueError,"v12.3"):
            run({},"h","s","d","l",{},{"complete":False},{"complete":False})
        with self.assertRaisesRegex(ValueError,"v12.4"):
            run({},"h","s","d","l",{},{"complete":True,"version":"v12.3"},{"complete":False})


if __name__=="__main__":unittest.main()

"""v12.7 original Czech vegetation class same-publication area year source math."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_CzechVeg_original_vegetation_CLASS_wind_source_v127 import (
    class_key,selected_source_class_blocks,source_class_panel,
    CLASS_CONDITIONS,run
)


def plot(kind,cl="XD",paper="Published botanical field study",
         region="18a",year=2003,area=16):
    return {
        "strict":kind,"broad":kind,"study":paper,
        "source_region":region,"survey_year":year,"area":area,
        "SOURCE_CLASS_CODE":cl,
    }


def biological_demo():
    header={}
    taxa={}
    for paper in ("Original study A","Original study B","Original study C"):
        for arm in ("rail","road"):
            for i in range(3):
                pid=f"{paper}_{arm}_{i}"
                header[pid]=plot(arm,paper=paper)
                taxa[pid]=[
                    (f"plant{life}{wind}{k}",
                     "Epilobium" if wind else "Allium",
                     "therophyte" if life else "hemicryptophyte")
                    for life in (False,True)
                    for wind in (False,True)
                    for k in range(4)
                    if i%3!=k%3
                ]
    return header,taxa


class CzechvegClassSourcePlantChecks(unittest.TestCase):
    def test_vegetation_class_key_preserves_actual_publication_region_decade_area(self):
        row=plot("rail",cl="XD",year=2003,area=16)
        self.assertEqual(class_key(row,CLASS_CONDITIONS[0]),
                         ("Published botanical field study","18a",2000,
                          "gt4_to_16","XD"))
        self.assertEqual(class_key(row,CLASS_CONDITIONS[1]),
                         ("Published botanical field study","18a",2000,None,"XD"))

    def test_strict_exact_original_study_class_requires_both_arms(self):
        meta,plants=biological_demo()
        for p,v in meta.items():
            if v["strict"]=="road":v["SOURCE_CLASS_CODE"]="XB"
        groups,stats=selected_source_class_blocks(meta,plants,"strict",
                                                  CLASS_CONDITIONS[0],2)
        self.assertEqual(stats["n_source_class_matched_original_publication_study_blocks"],0)
        self.assertEqual(stats["n_original_rail_matched_plots"],0)
        self.assertEqual(stats["n_original_road_matched_plots"],0)

    def test_original_source_trait_panels_stay_same_species_after_classification(self):
        meta,plants=biological_demo()
        result=source_class_panel(meta,plants,"broad",CLASS_CONDITIONS[0],2)
        stats=result["original_source_study_and_canonical_class_support"]
        effect=result["original_matched_source_plant_diaspore_response"]
        self.assertEqual(stats["n_original_studies"],3)
        self.assertEqual(stats["n_original_rail_matched_plots"],9)
        self.assertEqual(stats["n_original_road_matched_plots"],9)
        self.assertTrue(effect["identified"])
        self.assertEqual(effect["n_source_taxa"],16)
        self.assertIn("source_publication_cluster_resampling",result)

    def test_unknown_source_class_not_matched(self):
        meta,plants=biological_demo()
        for p,r in meta.items():
            if r["strict"]=="rail":r["SOURCE_CLASS_CODE"]=None
        matched,stat=selected_source_class_blocks(meta,plants,"strict",
                                                  CLASS_CONDITIONS[0],1)
        self.assertEqual(len(matched),0)
        self.assertEqual(stat["n_original_RAIL_source_plots_with_unknown_class"],9)

    def test_verified_prior_source_receipts_are_mandatory(self):
        with self.assertRaisesRegex(ValueError,"v12.5"):
            run({},"hdr","species","d","l",{},
                {"complete":False},{"complete":False})
        with self.assertRaisesRegex(ValueError,"v12.6"):
            run({},"hdr","species","d","l",{},
                {"complete":True,"version":"v12.5"},{"complete":False})


if __name__=="__main__":unittest.main()

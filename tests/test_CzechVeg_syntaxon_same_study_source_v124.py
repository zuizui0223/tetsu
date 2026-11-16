"""v12.4 original published Czech vegetation syntaxon levels without fabricated class hierarchy."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_CzechVeg_syntaxon_same_study_source_v124 import (
    clean_syntaxon,block_key,original_source_cross_arm_support,
    MATCH_LEVELS,MIN_SIDES,run,
)

def mock_plot(kind,syntax,code="18a",year=2001,area=16,study="Published source publication"):
    return {
        "study":study,"source_region":code,"survey_year":year,
        "area":area,"strict":kind,"broad":kind,
        "SYNTAXON":syntax,"SYNTAXON_NAME":("name_"+syntax if syntax else None),
        "SYNTAX_OLD":None,"SYNTAXON_OLD_NAME":None,
    }


class CzechVegOriginalVegetationSyntaxonChecks(unittest.TestCase):
    def test_missing_original_syntaxon_cannot_be_community_category(self):
        for x in (None,"","NA","nan","N/A"," unknown "," - ","?"):
            self.assertIsNone(clean_syntaxon(x))
        self.assertEqual(clean_syntaxon("  Phragmito -   Magnocaricetea "),
                         "Phragmito - Magnocaricetea")

    def test_original_study_metadata_not_conflated_with_vegetation_rank(self):
        p=mock_plot("rail","X11")
        a=MATCH_LEVELS[0]
        self.assertEqual(block_key(p,a),("Published source publication","18a",2000,
                                            "gt4_to_16","X11"))
        p2=mock_plot("road","X11",code="18b")
        self.assertNotEqual(block_key(p,a),block_key(p2,a))
        p3=mock_plot("road","X22")
        self.assertNotEqual(block_key(p,a),block_key(p3,a))

    def test_vegen_exact_class_requires_study_area_decade_and_same_syntaxon(self):
        samples={
            "ra":mock_plot("rail","X11"),
            "rb":mock_plot("rail","X11"),
            "da":mock_plot("road","X11"),
            "db":mock_plot("road","X11"),
            "rx":mock_plot("rail","X22"),
            "dx":mock_plot("road","X11",study="Other original paper"),
            "rm":mock_plot("rail",None),
            "dm":mock_plot("road",None),
        }
        p=original_source_cross_arm_support(samples,"strict",MATCH_LEVELS[0],2)
        self.assertEqual(p["n_original_source_both_arm_blocks"],1)
        self.assertEqual(p["source_rail_sites_in_common_vegetation_support"],2)
        self.assertEqual(p["source_road_sites_in_common_vegetation_support"],2)
        self.assertEqual(p["candidate_original_rail_syntaxon_missing"],1)
        self.assertEqual(p["candidate_original_road_syntaxon_missing"],1)
        q=original_source_cross_arm_support(samples,"strict",MATCH_LEVELS[-1],2)
        self.assertEqual(q["source_rail_sites_in_common_vegetation_support"],4)
        self.assertEqual(q["source_road_sites_in_common_vegetation_support"],3)

    def test_strict_and_broad_vegetation_source_sample_units_never_invented(self):
        samples={
            "r1":mock_plot("rail","A"),
            "r2":mock_plot("rail","A"),
            "d1":mock_plot("road","A"),
            "d2":mock_plot("road","B"),
        }
        hit=original_source_cross_arm_support(samples,"broad",MATCH_LEVELS[0],2)
        self.assertEqual(hit["n_original_source_both_arm_blocks"],0)
        hit2=original_source_cross_arm_support(samples,"broad",MATCH_LEVELS[0],1)
        self.assertEqual(hit2["n_original_source_both_arm_blocks"],1)

    def test_original_source_verification_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"not present"):
            run({},"missing_header",{"complete":False})


if __name__=="__main__":unittest.main()

"""Source/QC tests with synthetic iNaturalist response objects."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_japan_rail_weed_photo_coverage import (
    classify_rail_distance, matches_focal_name,
    summarize_taxon_observations
)


def ob(taxon, accuracy=15, obscured=False, phase=None):
    return {
        "taxon":{"name":taxon}, "positional_accuracy":accuracy,
        "obscured":obscured,"geoprivacy":None,"captive":False,
        "geojson":{"type":"Point","coordinates":[139.7,35.7]},
        "observed_on":"2024-09-20", "photos":[{"id":1}],
        "user":{"id":17},
        "annotations":[
            {"controlled_attribute_id":12,"controlled_value_id":phase}
        ] if phase is not None else []
    }


class ModernRailWeedsTests(unittest.TestCase):
    def test_distance_bins_are_predefined(self):
        self.assertEqual(classify_rail_distance(100),"0-100m")
        self.assertEqual(classify_rail_distance(101),"100-250m")
        self.assertEqual(classify_rail_distance(250),"100-250m")
        self.assertEqual(classify_rail_distance(251),"250-1000m")
        self.assertEqual(classify_rail_distance(1000),"250-1000m")
        self.assertEqual(classify_rail_distance(1001),"gt_1000m")

    def test_setaria_is_genus_not_single_species(self):
        self.assertTrue(matches_focal_name(ob("Setaria viridis"),"Setaria"))
        self.assertTrue(matches_focal_name(ob("Setaria faberi"),"Setaria"))
        self.assertFalse(matches_focal_name(ob("Solidago altissima"),"Setaria"))

    def test_accepted_kudzu_variety_alias_not_broad_species(self):
        self.assertTrue(matches_focal_name(
            ob("Pueraria montana var. lobata"),"Pueraria lobata"
        ))
        self.assertTrue(matches_focal_name(
            ob("Pueraria lobata"),"Pueraria lobata"
        ))
        self.assertTrue(matches_focal_name(
            ob("Pueraria montana lobata"),"Pueraria lobata"
        ))
        self.assertFalse(matches_focal_name(
            ob("Pueraria montana"),"Pueraria lobata"
        ))
        self.assertFalse(matches_focal_name(
            ob("Pueraria montana var. thomsonii"),"Pueraria lobata"
        ))

    def test_wrong_taxon_is_rejected_not_relabelled(self):
        self.assertFalse(matches_focal_name(ob("Solidago gigantea"),"Solidago altissima"))

    def test_nonflower_annotation_is_not_missing_for_flowering_species(self):
        data=[
            ob("Solidago altissima",phase=13),
            ob("Solidago altissima",phase=21),
            ob("Solidago altissima"),
            ob("Equisetum arvense",phase=13),
            ob("Solidago altissima",accuracy=None)
        ]
        summary=summarize_taxon_observations(data,"Solidago altissima",lambda lng,lat:150)
        self.assertEqual(summary["n_valid_geo"],3)
        self.assertEqual(summary["rail_bins"]["100-250m"],3)
        self.assertEqual(summary["flower_visibility_annotation_counts"]["visible_open_flower_annotation"],1)
        self.assertEqual(summary["flower_visibility_annotation_counts"]["no_visible_flower_or_fruit_annotation"],1)
        self.assertEqual(summary["flower_visibility_annotation_counts"]["phase_not_annotated"],1)
        self.assertEqual(summary["distinct_observers_per_rail_bin"]["100-250m"],1)
        self.assertEqual(summary["quality_gate"]["returned_taxon_does_not_match_exact_group"],1)
        self.assertEqual(summary["quality_gate"]["accuracy_not_within_100m"],1)

    def test_horsetail_strobilus_never_equated_to_flowering(self):
        data=[ob("Equisetum arvense",phase=13),
              ob("Equisetum arvense",phase=21),
              ob("Equisetum arvense")]
        result=summarize_taxon_observations(
            data,"Equisetum arvense",lambda lng,lat:80
        )
        self.assertEqual(result["n_valid_geo"],3)
        self.assertEqual(
            result["flower_visibility_annotation_counts"]["not_applicable_spore_plant"],3
        )
        self.assertEqual(
            result["reproductive_scoring_mode"],"spore_cone_separate_review_required"
        )
        self.assertNotIn(
            "visible_open_flower_annotation",
            result["flower_visibility_annotation_counts"]
        )

    def test_obscured_cultivated_unconfirmed_excluded(self):
        item=ob("Setaria viridis",obscured=True)
        summary=summarize_taxon_observations([item],"Setaria",lambda lng,lat:30)
        self.assertEqual(summary["n_valid_geo"],0)
        self.assertEqual(summary["quality_gate"]["obscured_or_private"],1)


if __name__=="__main__":
    unittest.main()

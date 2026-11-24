"""v342 Czech source study/habitat and vegetation cover metadata validation."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_CzechVeg_study_habitat_confounding_v342 as m

class CzechSourceConfounds(unittest.TestCase):
    def test_rail_study_strata_cross_tab(self):
        pair=[
         {"rail":{"reference":"R","cov_trees":0,"cov_shrubs":0,"altitude":400},
          "reference":{"reference":"R","cov_trees":0,"cov_shrubs":0,"altitude":400},
          "rail_habitat":"railway_track_yard_or_structure"},
         {"rail":{"reference":"R","cov_trees":0,"cov_shrubs":0,"altitude":400},
          "reference":{"reference":"R","cov_trees":0,"cov_shrubs":0,"altitude":400},
          "rail_habitat":"railway_embankment_ditch_or_cut"}]
        out=m.group_source_pairs(pair)
        self.assertEqual(out["n_rail_studies_in_matched_pairs"],1)
        self.assertEqual(out["n_original_rail_studies_with_both_track_and_embankment"],1)
    def test_rail_habitats_without_study_overlap_no_within_study_comparison(self):
        pair=[
         {"rail":{"reference":"A"},"rail_habitat":"railway_track_yard_or_structure"},
         {"rail":{"reference":"B"},"rail_habitat":"railway_embankment_ditch_or_cut"}]
        self.assertEqual(m.group_source_pairs(pair)[
          "n_original_rail_studies_with_both_track_and_embankment"],0)
    def test_open_vegetation_both_species_not_genotype(self):
        a={"rail":{"cov_trees":0,"cov_shrubs":10},
           "reference":{"cov_trees":0,"cov_shrubs":15}}
        self.assertEqual(m.numerically_valid_veg_open_pair(a),(True,"physiognomy_open_both"))
    def test_negative_cover_not_assumed_zero(self):
        a={"rail":{"cov_trees":-1,"cov_shrubs":0},
           "reference":{"cov_trees":0,"cov_shrubs":0}}
        self.assertEqual(m.numerically_valid_veg_open_pair(a),
           (False,"missing_original_tree_or_shrub_coverage"))
    def test_forest_like_reference_label(self):
        self.assertEqual(m.land_cover_group("lesní porost"),"forest_like")
    def test_grassland_label(self):
        self.assertEqual(m.land_cover_group("travnatá louka"),"grassland_like")
    def test_source_precision_not_assumed_known_without_coordinate(self):
        self.assertIsNone(m.cover_fraction(-1))
        self.assertEqual(m.cover_fraction(10),10)
    def test_frozen_original_gates(self):
        self.assertEqual(m.EXPECTED_MATCHED["G2_same_source_extended"],39)
        self.assertEqual(m.EXPECTED_MATCHED["G4_cross_source_extended"],49)
    def test_no_source_labels_emitted_for_location_ids(self):
        a=m.small_hash("REFERENCE ORIGINAL")
        self.assertEqual(len(a),10)
        self.assertNotIn("REFERENCE",a)

if __name__=="__main__":unittest.main()

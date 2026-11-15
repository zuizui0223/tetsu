"""v14.6 conditional SOURCE margins do not identify rail-section beta."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_France_Penone2012_thesis_occupancy_v145 import load_verified_subset
from prove_France_Penone2012_marginal_nonidentification_v146 import (
 exact_bounds,construct_witness_sites,direct_site_pair_contrast
)
SOURCE=ROOT/"data/literature/France_Penone2012_thesis_six_v145.csv"

class NonidentificationProof(unittest.TestCase):
 def test_six_author_margins_feasible_opposite_signs(self):
  source=load_verified_subset(SOURCE)
  result=exact_bounds(source)
  low=result["exact_minimum_within_minus_between_mean_six_species_shared_count"]
  high=result["exact_maximum_within_minus_between_mean_six_species_shared_count"]
  self.assertAlmostEqual(low,-0.025532879818594152)
  self.assertAlmostEqual(high,1.1175510204081633)
  self.assertTrue(low<0<high)
  for label,score in (("minimum",low),("maximum",high)):
   sites=construct_witness_sites(result["witness_"+label])
   self.assertEqual(len(sites),71)
   for species,f in source.items():
    self.assertEqual(sum(species in s for s in sites),f)
   self.assertAlmostEqual(direct_site_pair_contrast(sites,35),score)
 def test_abstract_partition_must_not_be_rebranded_as_true_rail_contrast(self):
  result=exact_bounds(load_verified_subset(SOURCE))
  self.assertTrue(result["hypothetical_groups_not_actual_railway_lines"])
  self.assertTrue(result["conditional_original_observations_equal_unique_sites_UNVERIFIED"])
  self.assertFalse(result["complete_71_site_by_species_matrix_obtained"])
  self.assertFalse(result["real_railway_connectivity_or_train_seed_dispersal_identified"])
 def test_invalid_units_block(self):
  with self.assertRaises(ValueError):exact_bounds({"x":72})
  with self.assertRaises(ValueError):exact_bounds({"x":0})
  with self.assertRaises(ValueError):exact_bounds({"x":1},1,70)

if __name__=="__main__":unittest.main()

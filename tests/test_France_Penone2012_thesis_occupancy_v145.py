"""Verify original six-count transcription and no invented connectivity effect."""
import csv
import sys
import tempfile
import unittest
from math import comb
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_France_Penone2012_thesis_occupancy_v145 import (
    EXACT,load_verified_subset,marginal_moments,analyze
)
SOURCE=ROOT/"data/literature/France_Penone2012_thesis_six_v145.csv"

class SourceUnitGates(unittest.TestCase):
    def test_true_thesis_subset_golden_counts(self):
        self.assertEqual(load_verified_subset(SOURCE),
                         {name:v[0] for name,v in EXACT.items()})
    def test_exact_pair_moment(self):
        r=marginal_moments(load_verified_subset(SOURCE))
        self.assertEqual(r["source_observations_sum_six"],220)
        self.assertEqual(r["unordered_site_pairs"],2485)
        self.assertEqual(r["six_species_shared_presence_pair_total_conditional"],5345)
        self.assertAlmostEqual(r["six_species_contribution_mean_intersection_conditional"],5345/2485)
    def test_bruteforce_three_site_binary_matrix(self):
        sites=[{"a","b"},{"a","c"},{"a","b"}]
        f={s:sum(s in x for x in sites) for s in ("a","b","c")}
        brute=sum(len(a&b) for i,a in enumerate(sites) for b in sites[i+1:])/3
        self.assertEqual(brute,4/3)
        self.assertEqual(marginal_moments(f,3)["six_species_contribution_mean_intersection_conditional"],brute)
    def test_bad_transcription_and_count_fail(self):
        with self.assertRaises(ValueError):marginal_moments({"a":72},71)
        with self.assertRaises(ValueError):marginal_moments({"a":0},71)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"bad.csv"
            with p.open("w",newline="",encoding="utf-8") as f:
                w=csv.writer(f)
                w.writerow(["species","observations","page"])
                w.writerow(["Brassica napus",4,175])
            with self.assertRaises(ValueError):load_verified_subset(p)
    def test_no_fabricated_rail_or_meta_identification(self):
        r=analyze(load_verified_subset(SOURCE))
        self.assertFalse(r["distinct_site_observation_grain_verified"])
        self.assertFalse(r["complete_71_site_by_species_matrix"])
        self.assertFalse(r["can_calculate_distance_matched_Jaccard_or_route_break_effect"])
        self.assertFalse(r["can_infer_train_seed_dispersal_or_stress_adaptation"])
        self.assertFalse(r["added_independent_cohort_to_17_original_study_inventory"])

if __name__=="__main__":unittest.main()

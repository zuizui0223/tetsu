"""Synthetic contract tests only. All ecological results need original source CI."""
import importlib.util
import unittest
from pathlib import Path

p = (Path(__file__).resolve().parents[1] / "scripts" /
     "validate_CzechVeg_strict_habitat_era_holdout_v200.py")
spec = importlib.util.spec_from_file_location("v200_source_gate", p)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def synthetic():
    # Study, source-region, survey-year and plot area are originals-in-shape
    # but ALL values below are synthetic and cannot be cited as field results.
    taxa = [f"S{i:02d} species" for i in range(20)]
    life = {
        t: ("therophyte" if i < 6 else
            "hemicryptophyte" if i < 12 else "geophyte")
        for i, t in enumerate(taxa)
    }
    disp = {
        t: ("Epilobium" if i % 3 == 0 else
            "Allium" if i % 3 == 1 else "Phragmites")
        for i, t in enumerate(taxa)
    }
    meta, plants = {}, {}
    for i in range(12):
        paper = f"synthetic_author_study_{i:02d}"
        year = 1995 if i < 3 else 2005 if i < 6 else 2015
        for arm in ("rail", "road"):
            for j in range(3):
                pid = f"{paper}_{arm}_{j}"
                meta[pid] = {
                    "original_reference_title": paper,
                    "source_phyto_region": f"R{i % 4}",
                    "source_year_int": year,
                    "area": 10.0,
                    "strict": arm if i < 7 else None,
                    "broad": arm,
                }
                group = taxa[:12] if arm == "rail" else taxa[6:18]
                plants[pid] = [(t, disp[t], life[t]) for t in group]
    return meta, plants


class TestOriginalLabelAndEraSourceGate(unittest.TestCase):
    def test_source_survey_year_not_publication_year(self):
        self.assertTrue(m.year_in_era(1995, "before_2000"))
        self.assertFalse(m.year_in_era(1995, "year_2000_onward"))
        self.assertTrue(m.year_in_era(2005, "year_2000_onward"))
        self.assertFalse(m.year_in_era(2005, "year_2010_onward"))
        self.assertTrue(m.year_in_era(2015, "year_2010_onward"))
        with self.assertRaises(ValueError):
            m.year_in_era(2015, "after_2030")

    def test_strict_label_loses_locality_only_and_era_blocks(self):
        meta, plants = synthetic()
        b, _ = m.matched_label_era_blocks(meta, plants, "broad", "all_years")
        s, _ = m.matched_label_era_blocks(meta, plants, "strict", "all_years")
        modern, _ = m.matched_label_era_blocks(
            meta, plants, "broad", "year_2010_onward")
        self.assertEqual(len(b), 12)
        self.assertEqual(len(s), 7)
        self.assertEqual(len(modern), 6)
        self.assertEqual(sum(len(v["rail"]) for v in b.values()), 36)
        self.assertEqual(sum(len(v["rail"]) for v in s.values()), 21)
        self.assertEqual(sum(len(v["rail"]) for v in modern.values()), 18)

    def test_source_original_frame_identity_check(self):
        meta, plants = synthetic()
        blocks, _ = m.matched_label_era_blocks(
            meta, plants, "broad", "all_years")
        m.assert_original_broad_source_frame(blocks, blocks)
        first_key = next(iter(blocks))
        broken = {k: {arm: list(v[arm]) for arm in ("rail", "road")}
                  for k, v in blocks.items()}
        broken[first_key]["rail"].pop()
        with self.assertRaises(AssertionError):
            m.assert_original_broad_source_frame(broken, blocks)

    def test_underpowered_strict_is_not_reported_as_zero_effect(self):
        meta, plants = synthetic()
        strict, _ = m.matched_label_era_blocks(
            meta, plants, "strict", "all_years")
        result = m.run_panel(strict, plants, "all_vascular_min5")
        self.assertFalse(result["identified"])
        self.assertTrue(result["underpowered_not_zero_effect"])
        self.assertIn("under_8_original_publications", result["reasons"])

    def test_study_disjoint_out_of_source_metrics_and_chance_baseline(self):
        meta, plants = synthetic()
        broad, _ = m.matched_label_era_blocks(
            meta, plants, "broad", "all_years")
        result = m.run_panel(broad, plants, "both_trait50_min5")
        self.assertTrue(result["identified"])
        self.assertEqual(result["untrained_constant_score_baseline_AUC"], .5)
        self.assertEqual(result["source_support"][
            "n_retained_original_publications"], 12)
        self.assertEqual(len(result[
            "per_study_AUC_only_hashed_original_study_titles"]), 12)
        self.assertEqual(
            set(result["original_publication_equal_holdout_AUC"]),
            set(m.MODELS))
        self.assertTrue(result["same_original_publications_for_all_models"])
        self.assertTrue(result[
            "publication_fold_assignment_pre_model_using_only_plot_opportunity"])
        self.assertEqual(sum(result["publication_holdout_fold_size"]), 12)

    def test_descriptive_bootstrap_not_field_confidence_interval(self):
        a = m.study_bootstrap_descriptive_by_paper(
            [0.1] * 9, n=100, seed=10)
        self.assertAlmostEqual(a["mean_original_publication_equal"], .1)
        self.assertAlmostEqual(
            a["percentile_2_5_study_resampling_NOT_independent_field_CI"],
            .1)
        self.assertAlmostEqual(
            a["percentile_97_5_study_resampling_NOT_independent_field_CI"],
            .1)


if __name__ == "__main__":
    unittest.main()

import importlib.util
import unittest
from collections import Counter
from pathlib import Path

script = (
    Path(__file__).resolve().parents[1]
    / "scripts/validate_CzechVeg_publication_holdout_rail_species_v197.py"
)
spec = importlib.util.spec_from_file_location("source_v197", script)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class TestOriginalPublicationHeldoutValidation(unittest.TestCase):
    def synthetic(self):
        # Twelve independent publications, no study may enter both arms.
        blocks = {}
        original_plant = {}
        for i in range(12):
            paper = f"pub_{i:02d}"
            block = (paper, f"region_{i % 4}", 2010, "gt4_to_16")
            rids = [f"{paper}_r0", f"{paper}_r1"]
            dids = [f"{paper}_d0", f"{paper}_d1"]
            blocks[block] = {"rail": rids, "road": dids}
            for rid in rids:
                original_plant[rid] = {"Echium vulgare", "Artemisia vulgaris"}
            for did in dids:
                original_plant[did] = {"Poa compressa", "Artemisia vulgaris"}
        station = Counter({
            "Echium vulgare": 39,
            "Poa compressa": 35,
            "Artemisia vulgaris": 39,
            "Conyza canadensis": 39,
        })
        return blocks, original_plant, station

    def test_study_split_deterministic_and_disjoint(self):
        blocks, _, _ = self.synthetic()
        folds, weights = m.split_source_studies(blocks, nfold=4)
        self.assertEqual(folds, m.split_source_studies(blocks)[0])
        self.assertEqual(sum(len(x) for x in folds), 12)
        self.assertEqual(sum(weights), 24)
        self.assertTrue(all(len(x) == 3 for x in folds))
        for i in range(4):
            for j in range(i + 1, 4):
                self.assertTrue(folds[i].isdisjoint(folds[j]))

    def test_frozen_taxa_train_and_test_exclude_entire_publications(self):
        blocks, plot_taxa, station = self.synthetic()
        folds = m.fold_diagnostic(blocks, plot_taxa, station, Counter(), 100)
        self.assertEqual(len(folds), 4)
        for fold in folds:
            self.assertEqual(fold["study_n_heldout"], 3)
            self.assertEqual(fold["study_n_train"], 9)
            self.assertEqual(fold["heldout_original_RAIL_text_plot_n"], 6)
            self.assertEqual(fold["heldout_original_ROAD_text_plot_n"], 6)
            self.assertGreater(fold["locked_v174_taxa"]["Echium vulgare"]["heldout_effect"], 0)
            self.assertLess(fold["locked_v174_taxa"]["Poa compressa"]["heldout_effect"], 0)
            self.assertIsNone(fold["locked_v174_taxa"]["Conyza canadensis"]["heldout_effect"])
        summary = m.summarize(folds)
        self.assertEqual(summary["locked_prior_v174_taxa"]["Echium vulgare"]["heldout_positive_folds"], 4)
        self.assertEqual(summary["locked_prior_v174_taxa"]["Conyza canadensis"]["heldout_supported_folds"], 0)

    def test_train_positive_does_not_force_heldout_positive(self):
        def result(delta, studies=3, incidence=8):
            return {
                m.SCORE: delta,
                m.STATION_N: 39,
                m.RAIL_N: incidence,
                m.ROAD_N: 0,
                m.PAPER_N: studies,
            }
        training = {"Echium vulgare": result(0.23)}
        holdout = {"Echium vulgare": result(-0.13)}
        self.assertEqual(m.rank_train(training, "positive"), ["Echium vulgare"])
        out = m.validate_rank(training, holdout, "positive")
        self.assertEqual(out["heldout_evaluable_n"], 1)
        self.assertEqual(out["heldout_direction_replicated_n"], 0)
        self.assertEqual(out["heldout_median_weighted_plot_prevalence_difference"], -0.13)

    def test_sparse_heldout_is_missing_not_a_negative_result(self):
        def result(delta, studies, incidence):
            return {
                m.SCORE: delta, m.STATION_N: 39,
                m.RAIL_N: incidence, m.ROAD_N: 0, m.PAPER_N: studies,
            }
        train = {"Echium vulgare": result(0.23, 3, 8)}
        test = {"Echium vulgare": result(-0.23, 1, 1)}
        out = m.validate_rank(train, test, "positive")
        self.assertEqual(out["heldout_evaluable_n"], 0)
        self.assertEqual(out["heldout_missing_min_source_support_n"], 1)
        self.assertIsNone(out["heldout_direction_replicated_fraction_DESCRIPTIVE_ONLY"])


if __name__ == "__main__":
    unittest.main()

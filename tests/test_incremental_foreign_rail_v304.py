import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_incremental_foreign_rail_v304 import fit_predict, metrics, national_design, block_frequencies, crossfit_training_frequencies


class IncrementalTests(unittest.TestCase):
    def test_signal_improves_unseen_rows(self):
        x = np.linspace(-2, 2, 100)[:, None]
        y = 1 / (1 + np.exp(-2 * x[:, 0]))
        train = np.arange(100) % 2 == 0
        p, _ = fit_predict(x[train], y[train], np.ones(50), x[~train], .001)
        q, _ = fit_predict(np.empty((50, 0)), y[train], np.ones(50), np.empty((50, 0)), .001)
        self.assertGreater(metrics(y[~train], q)['log_loss'] - metrics(y[~train], p)['log_loss'], .1)

    def test_fractional_loss_equals_bernoulli_rows(self):
        a = metrics(np.array([.25]), np.array([.6]))
        b = metrics(np.array([0., 0., 0., 1.]), np.full(4, .6))
        self.assertAlmostEqual(a['log_loss'], b['log_loss'])
        self.assertAlmostEqual(a['brier'], b['brier'])

    def test_held_study_outcomes_cannot_change_design(self):
        taxa = [f'Species {i}' for i in range(25)]
        rail = {k: {t: .1 + i / 100 for i, t in enumerate(taxa)} for k in 'ABCDE'}
        refs = {k: {t: .2 for t in taxa} for k in rail}
        first = national_design(rail, refs, 'A')
        rail['A'] = {'Unseen plant': 1.}
        second = national_design(rail, refs, 'A')
        for key in ('x', 'y', 'w', 'test_x'):
            np.testing.assert_array_equal(first[key], second[key])
        self.assertEqual(first['taxa'], second['taxa'])

    def test_held_block_not_in_frequency(self):
        matrix = np.array([[1., 0.], [0., 1.], [0., 0.]])
        use = np.array([False, True, True])
        before = block_frequencies(matrix, use)
        matrix[0] = [0., 1.]
        np.testing.assert_array_equal(before, block_frequencies(matrix, use))

    def test_constant_columns_and_one_class_fit(self):
        p, _ = fit_predict(np.ones((20, 2)), np.zeros(20), np.ones(20), np.ones((3, 2)), .01)
        self.assertTrue(np.isfinite(p).all())
        self.assertTrue((p < .001).all())

    def test_target_encoding_excludes_own_training_outcome(self):
        y = np.array([[1., 0.], [0., 1.], [0., 0.]])
        train = np.array([True, True, False])
        first = crossfit_training_frequencies(y, train)
        y[0] = [0., 1.]
        second = crossfit_training_frequencies(y, train)
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[2], [.5, .5])


if __name__ == '__main__':
    unittest.main()

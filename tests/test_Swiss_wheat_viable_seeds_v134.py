"""Numerical and denominator tests; synthetic PDF layout below is not evidence."""
import importlib.util
import math
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts' / 'analyze_Swiss_wheat_viable_seeds_v134.py'
DATA = ROOT / 'data' / 'Swiss_Schulze2015_imported_wheat_germination_v134.csv'

class CohortAndSourceTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.exists(), 'viability analysis module not yet implemented')
        spec = importlib.util.spec_from_file_location('v134', SCRIPT)
        self.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.m)
        self.rows = self.m.read_data(DATA)

    def test_actual_same_cohort_fractions(self):
        result = self.m.analyze(self.rows)
        self.assertEqual(result['pooled_descriptive']['sown'], 300)
        self.assertEqual(result['pooled_descriptive']['germinated'], 75)
        self.assertEqual(result['pooled_descriptive']['GM_seedlings'], 40)
        for r, expected in zip(result['by_mill'], (16/150, 24/150)):
            self.assertAlmostEqual(r['GM_seedling_per_sown_seed'], expected)
            self.assertAlmostEqual(r['germinated_fraction'] * r['GM_fraction_of_seedlings'], expected)
        self.assertFalse(result['claims']['field_establishment_probability_measured'])

    def test_five_genotype_classes_are_mutually_exclusive(self):
        bad = [dict(r) for r in self.rows]
        bad[0]['MS8'] += bad[0]['MS8_RF3']
        with self.assertRaisesRegex(ValueError, 'genotype'):
            self.m.analyze(bad)

    def test_input_rejects_duplicate_mills_and_negative_counts(self):
        for rows in ([self.rows[0], self.rows[0]],
                     [dict(self.rows[0], nonGM=-1), self.rows[1]]):
            with self.assertRaises(ValueError):
                self.m.analyze(rows)

    def test_pool_positive_is_not_seed_prevalence(self):
        result = self.m.analyze(self.rows)
        self.assertEqual(result['pool_tests']['GT73_positive_pools'], 18)
        self.assertEqual(result['pool_tests']['total_pools'], 20)
        self.assertIsNone(result['pool_tests']['individual_seed_GT73_prevalence'])

    def test_historical_projection_units_and_joint_mill_scenarios(self):
        p = self.m.analyze(self.rows)['historical_conditional_projection']
        self.assertAlmostEqual(p['contaminant_mass_tonnes'], 3.8756)
        self.assertAlmostEqual(p['all_contaminant_seeds'], 1211125000)
        self.assertAlmostEqual(p['GM_seedling_equivalent_by_mill']['1'], 129186666.66666667)
        self.assertAlmostEqual(p['GM_seedling_equivalent_by_mill']['2'], 193780000)
        self.assertIsNone(p['rail_spill_count'])
        self.assertIsNone(p['established_plants'])
        self.assertFalse(p['two_mill_range_is_confidence_interval'])

    def test_zero_spill_or_zero_establishment_cannot_be_ignored(self):
        q = 40/300
        self.assertEqual(self.m.field_scenario(1000, q, 0, 0.5), 0)
        self.assertEqual(self.m.field_scenario(1000, q, 0.5, 0), 0)
        self.assertAlmostEqual(self.m.field_scenario(1000, q, 0.1, 0.2), 8/3)
        self.assertIsNone(self.m.field_scenario(1000, q, None, 0.2))
        with self.assertRaises(ValueError):
            self.m.field_scenario(1000, q, 1.1, 0.2)

    def test_pdf_table_verification_is_row_specific(self):
        text = '''10.1007/s11356-015-4903-y
Low level impurities in imported wheat
Table 2 a Numbers of seedlings
Grain mill No. of seeds sown GT73 MS8 RF3
1 150 24 16 3 5 4 4 8 67
2 150 51 34 7 7 5 5 27 47
b Grain mill No. of pool samples
1 10 9 10 10
2 10 9 10 10
'''
        proof = self.m.verify_table_text(text, self.rows)
        self.assertEqual(proof['rows_verified'], 4)
        with self.assertRaisesRegex(ValueError, 'Table 2a'):
            self.m.verify_table_text(text.replace('3 5 4 4 8', '4 4 4 4 8'), self.rows)

    def test_same_total_does_not_validate_changed_cell(self):
        text = '''10.1007/s11356-015-4903-y
Table 2 a Numbers of seedlings
1 150 24 16 4 4 4 4 8 67
2 150 51 34 7 7 5 5 27 47
1 10 9 10 10
2 10 9 10 10
'''
        with self.assertRaises(ValueError):
            self.m.verify_table_text(text, self.rows)

if __name__ == '__main__':
    unittest.main()

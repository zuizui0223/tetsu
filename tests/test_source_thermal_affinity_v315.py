"""Behavior checks for fixed thermal trait definitions and paired source sets."""

from collections import Counter
import importlib
import importlib.util
import os
from pathlib import Path
import re
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))


def normalize(label):
    if not isinstance(label, str):
        return None
    if re.search(r'\b(?:agg|group|cf|aff)\b', label, re.I) or '×' in label:
        return None
    match = re.match(r'^([A-Z][A-Za-z-]+)\s+([a-z][a-z-]+)(?=\s|$)', label)
    if match and match[2] not in {'sp', 'spp'}:
        return match[1] + ' ' + match[2]
    return None


class AdapterTestCase(unittest.TestCase):
    def setUp(self):
        # A missing production adapter should be an explicit RED assertion.
        self.assertIsNotNone(importlib.util.find_spec('source_thermal_affinity_v315'),
                             'The v315 thermal source adapter is not implemented')
        self.source = importlib.import_module('source_thermal_affinity_v315')


class ThermalTraits(AdapterTestCase):
    def test_tichy_indifference_and_missing_never_enter_numeric_values(self):
        trait = self.source.parse_tichy_rows([
            ('Abies alba', 4.7), ('Poa pratensis', 'x'),
            ('Pinus nigra', 'NA'), ('Acer campestre', None),
        ], normalize)
        self.assertEqual(trait['values'], {'Abies alba': 4.7})
        self.assertEqual(trait['statuses'], {
            'Abies alba': 'numeric', 'Poa pratensis': 'indifferent',
            'Pinus nigra': 'missing', 'Acer campestre': 'missing',
        })
        self.assertEqual(trait['limits'], [1, 12])

    def test_eive_zero_is_a_valid_temperature_value(self):
        trait = self.source.parse_eive_rows([
            ('Abies alba', 'Species', 0), ('Poa pratensis', 'Species', 10),
            ('Pinus nigra', 'Species', None),
        ], normalize)
        self.assertEqual(trait['values'], {'Abies alba': 0, 'Poa pratensis': 10})
        self.assertEqual(trait['statuses']['Pinus nigra'], 'missing')
        self.assertEqual(trait['limits'], [0, 10])

    def test_values_outside_each_scale_fail_closed(self):
        for value in (0, 12.01, float('nan'), float('inf'), True, '5', '?'):
            with self.subTest(source='tichy', value=value):
                with self.assertRaises(ValueError):
                    self.source.parse_tichy_rows([('Abies alba', value)], normalize)
        for value in (-.01, 10.01, float('nan'), False, 'x'):
            with self.subTest(source='eive', value=value):
                with self.assertRaises(ValueError):
                    self.source.parse_eive_rows([('Abies alba', 'Species', value)], normalize)

    def test_strict_traits_exclude_aggregates_hybrids_and_infraspecific_names(self):
        labels = ['Poa pratensis aggr.', 'Poa pratensis agg.',
                  'Poa pratensis subsp. pratensis', 'Poa pratensis ssp pratensis',
                  'Poa pratensis var. pratensis', 'Poa pratensis f. pratensis',
                  'Poa pratensis x trivialis', 'Poa pratensis × trivialis',
                  'Poa sect. Poa', 'Poa pratensis s. l.', 'Poa sp.']
        trait = self.source.parse_tichy_rows([(name, 5) for name in labels] + [
            ('Poa pratensis', 4),
        ], normalize)
        self.assertEqual(trait['values'], {'Poa pratensis': 4})
        self.assertEqual(trait['statuses'], {'Poa pratensis': 'numeric'})

    def test_mislabeled_eive_subspecies_cannot_override_species(self):
        trait = self.source.parse_eive_rows([
            ('Lupinus angustifolius', 'Species', 6.6),
            ('Lupinus angustifolius subsp angustifolius', 'Species', 5.8),
            ('Abies alba', 'Subspecies', 4),
            ('Poa pratensis', 'Species Aggregate', 3),
        ], normalize)
        self.assertEqual(trait['values'], {'Lupinus angustifolius': 6.6})
        self.assertEqual(trait['statuses'], {'Lupinus angustifolius': 'numeric'})

    def test_identical_duplicate_values_create_one_species_value(self):
        trait = self.source.parse_tichy_rows([
            ('Abies alba', 4), ('Abies alba Mill.', 4.0),
            ('Poa pratensis', 'x'), ('Poa pratensis L.', 'x'),
        ], normalize)
        self.assertEqual(trait['values'], {'Abies alba': 4})
        self.assertEqual(trait['statuses'], {'Abies alba': 'numeric',
                                            'Poa pratensis': 'indifferent'})
        self.assertEqual(trait['audit']['duplicate_normalized_keys'], 2)

    def test_conflicting_numeric_or_status_duplicates_are_unknown(self):
        rows = [('Abies alba', 4), ('Abies alba Mill.', 5),
                ('Poa pratensis', 'x'), ('Poa pratensis L.', 4),
                ('Pinus nigra', 'NA'), ('Pinus nigra J.F.Arnold', 'x')]
        for ordered in (rows, list(reversed(rows))):
            trait = self.source.parse_tichy_rows(ordered, normalize)
            self.assertEqual(trait['values'], {})
            self.assertEqual(trait['statuses'], {'Abies alba': 'conflict',
                                                'Poa pratensis': 'conflict',
                                                'Pinus nigra': 'conflict'})
            self.assertEqual(trait['audit']['conflicting_normalized_keys'], 3)

    def test_eive_duplicate_conflict_never_selects_one_habitat_value(self):
        trait = self.source.parse_eive_rows([
            ('Abies alba', 'Species', 4), ('Abies alba Mill.', 'Species', 7),
        ], normalize)
        self.assertEqual(trait['values'], {})
        self.assertEqual(trait['statuses'], {'Abies alba': 'conflict'})

    def test_bad_checksum_is_rejected_before_workbook_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, 'tichy_2023_zenodo_v2.xlsx').write_bytes(b'not a workbook')
            with self.assertRaisesRegex(ValueError, 'SHA-256'):
                self.source.load_traits(directory, normalize)


class KashmirPairs(AdapterTestCase):
    def fixtures(self):
        headers = ['Abies alba', 'Poa pratensis', 'Poa pratensis agg.',
                   'Pinus nigra', 'Acer campestre']
        plots = {
            ('Station A', 2014, 'Track'): {0, 1},
            ('Station A', 2014, 'Embankment'): {2},
            ('Station A', 2014, 'Interior Plot 1'): {1, 3},
            ('Station A', 2014, 'Interior Plot 2'): {3},
            ('Station B', 2014, 'Track'): {4},
            ('Station B', 2014, 'Embankment'): {4},
            ('Station B', 2014, 'Interior Plot 1'): {0, 1},
            ('Station B', 2014, 'Interior Plot 2'): {3},
            ('Between stations', 2014, 'Track'): {1},
            ('Station A', 2017, 'Track'): {4},
        }
        return plots, headers

    def test_pairs_use_baseline_station_unions_and_shared_ambiguity_masks(self):
        plots, headers = self.fixtures()
        pairs = self.source.pairs_from_kashmir_plots(
            plots, headers, normalize, expected_stations=2)
        self.assertEqual(pairs, [
            ({'Abies alba'}, {'Pinus nigra'}),
            ({'Acer campestre'}, {'Abies alba', 'Pinus nigra'}),
        ])

    def test_incomplete_station_pair_is_not_silently_dropped(self):
        plots, headers = self.fixtures()
        del plots['Station B', 2014, 'Interior Plot 2']
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            self.source.pairs_from_kashmir_plots(
                plots, headers, normalize, expected_stations=2)

    def test_missing_plot_or_wrong_station_count_is_rejected(self):
        plots, headers = self.fixtures()
        with self.assertRaisesRegex(ValueError, 'station count'):
            self.source.pairs_from_kashmir_plots(plots, headers, normalize)
        plots['Station B', 2014, 'Track'] = None
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            self.source.pairs_from_kashmir_plots(
                plots, headers, normalize, expected_stations=2)

    def test_species_column_indices_are_validated(self):
        plots, headers = self.fixtures()
        plots['Station A', 2014, 'Track'] = {len(headers)}
        with self.assertRaisesRegex(ValueError, 'column index'):
            self.source.pairs_from_kashmir_plots(
                plots, headers, normalize, expected_stations=2)


@unittest.skipUnless(os.environ.get('TETSU_V315_TRAITS'), 'Pinned external traits not supplied')
class PinnedTraitIntegration(AdapterTestCase):
    def test_pinned_workbooks_have_the_documented_status_and_rank_census(self):
        from analyze_matched_foreign_rail_road_v314 import normalize_name
        traits = self.source.load_traits(os.environ['TETSU_V315_TRAITS'], normalize_name)
        self.assertEqual(traits['tichy']['audit']['source_status_counts'],
                         {'numeric': 7402, 'indifferent': 1147, 'missing': 359})
        self.assertEqual(traits['eive']['audit']['source_status_counts'],
                         {'numeric': 14496, 'missing': 339})
        self.assertEqual(Counter(traits['tichy']['statuses'].values()),
                         {'numeric': 7217, 'indifferent': 1083, 'missing': 357})
        self.assertEqual(Counter(traits['eive']['statuses'].values()),
                         {'numeric': 10398, 'missing': 291})


@unittest.skipUnless(os.environ.get('TETSU_V315_DATA'), 'Pinned external inventories not supplied')
class PinnedKashmirIntegration(AdapterTestCase):
    def test_station_pairs_reproduce_pinned_v314_frequencies_after_shared_masking(self):
        from analyze_matched_foreign_rail_road_v314 import normalize_name
        from source_matched_rail_road_v314 import read_kashmir
        directory = Path(os.environ['TETSU_V315_DATA'])
        pairs = self.source.kashmir_pairs(directory, normalize_name)
        self.assertEqual(len(pairs), 16)
        rail, interior = read_kashmir(directory / 'Matrix_plot.xlsx', normalize_name)
        exclude = rail.unavailable | interior.unavailable
        for side, reference in enumerate((rail, interior)):
            counts = Counter(key for pair in pairs for key in pair[side])
            expected = {key: value for key, value in reference.frequencies.items()
                        if key not in exclude}
            self.assertEqual(set(counts), set(expected))
            for key, value in expected.items():
                self.assertAlmostEqual(counts[key] / 16, value)


if __name__ == '__main__':
    unittest.main()

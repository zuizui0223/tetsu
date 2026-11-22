"""Tests for source frequencies and their unavailable-name masks."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

import source_matched_rail_road_v314 as source


def normalize(label):
    if isinstance(label, str) and not label.endswith((' agg.', ' sp.')):
        words = label.split()
        if len(words) >= 2:
            return ' '.join(words[:2])
    return None


class SourceFrequencies(unittest.TestCase):
    def test_marginals_mask_duplicate_and_unresolved_species_keys(self):
        dataset = source.marginal_frequencies([
            ('Abies alba', 1 / 30),
            ('Abies alba subsp. alba', 2 / 30),
            ('Poa pratensis agg.', 5 / 30),
            ('Pinus nigra', 1 / 30),
        ], normalize)
        self.assertEqual(dataset.frequencies, {'Pinus nigra': 1 / 30})
        self.assertEqual(dataset.unavailable, {'Abies alba', 'Poa pratensis'})

    def test_sites_collapse_rank_and_weight_groups_equally(self):
        dataset = source.site_frequencies([
            [{'Abies alba', 'Abies alba subsp. alba'}, {'Poa pratensis agg.'}],
            [{'Abies alba', 'Pinus nigra'}],
        ], normalize, {})
        self.assertEqual(dataset.frequencies, {'Abies alba': .75, 'Pinus nigra': .5})
        self.assertEqual(dataset.unavailable, {'Poa pratensis'})

    def test_slovenian_conflicts_mask_each_candidate_only_in_matched_zones(self):
        def observation(location, zone, code, name):
            return {'location_id': location, 'zone_id': zone, 'EPPO_code': code, 'latin_name': name}

        rows = [
            observation(1, 'LBS', 'A', 'Abies alba'),
            observation(1, 'LBT', 'A', 'Abies alba var. alba'),
            observation(1, 'LBT', 'H', 'Heracleum sphondylium'),
            observation(1, 'LBS', 'H', 'Heracleum sphondylium'),
            observation(2, 'RBT', 'L', 'Lotus herbaceus'),
            observation(2, 'RBS', 'X', 'Uncertain sp.'),
            observation(2, 'IA', 'OUT', 'Pinus nigra'),
            observation(2, 'LBS', 'A', 'Abies dubia'),
        ]
        matched = {(1, 'LBS'), (1, 'LBT'), (2, 'RBS'), (2, 'RBT')}
        known = {'A': 'Abies alba', 'H': 'Heracleum sphondyllium'}
        unavailable, labels, audit = source.slovenian_conflicts(rows, matched, known, normalize)
        self.assertEqual(unavailable, {'Heracleum sphondylium', 'Heracleum sphondyllium', 'Lotus herbaceus'})
        self.assertNotIn('Pinus nigra', labels)
        self.assertNotIn('Abies dubia', labels)
        self.assertEqual(audit['field_codebook_conflicting_codes'], 1)
        self.assertEqual(audit['field_codebook_unmatched_strict_codes'], 1)
        self.assertEqual(audit['field_codebook_conflict_site_code_incidences'], 2)


if __name__ == '__main__':
    unittest.main()

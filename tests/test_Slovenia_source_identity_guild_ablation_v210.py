#!/usr/bin/env python3
"""Synthetic mathematics only, NEVER publish these synthetic counts as ecology."""
from __future__ import annotations
import sys
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_Slovenia_source_identity_guild_ablation_v210 import (
    MODELS, evaluate, occupancy_frequency_bins, probabilities,
    shuffled_guild_rows,
)


def make(site, section, taxon, guild, state):
    return {
        'internal_original_site': site,
        'author_original_section': section,
        'source_taxon_EPPO': taxon,
        'source_four_lifespan_diaspore_guild': guild,
        'observed_original_microhabitat_state': state,
    }


class FalsificationMathematics(unittest.TestCase):
    def setUp(self):
        self.rows = [
            make('1', 's1', 'A', 'strict_annual_wind', 'BS_only'),
            make('1', 's1', 'B', 'other_known_other', 'BT_only'),
            make('2', 's2', 'A', 'strict_annual_wind', 'both'),
            make('2', 's2', 'C', 'strict_annual_other', 'BT_only'),
            make('3', 's3', 'A', 'strict_annual_wind', 'BS_only'),
            make('3', 's3', 'B', 'other_known_other', 'BT_only'),
            make('4', 's4', 'A', 'strict_annual_wind', 'both'),
            make('4', 's4', 'D', 'other_known_wind', 'BT_only'),
        ]

    def test_all_models_have_proper_probabilities(self):
        fit = {
            'global': {'BS_only': .2, 'BT_only': .3, 'both': .5},
            'guild': {'g1': {'BS_only': .7, 'BT_only': .2, 'both': .1}},
            'taxon': {'A': Counter({'BS_only': 3, 'BT_only': 1})},
            'source_taxon_guild': {'A': 'g1'},
            'states': ('BS_only', 'BT_only', 'both'),
        }
        p, seen = probabilities(make('1', 's1', 'A', 'g1', 'BS_only'), fit)
        self.assertTrue(seen)
        self.assertAlmostEqual(p[MODELS[2]]['BS_only'], (3 + 10 * .2) / 14)
        self.assertAlmostEqual(p[MODELS[3]]['BS_only'], (3 + 10 * .7) / 14)
        self.assertAlmostEqual(sum(p[MODELS[2]].values()), 1)

    def test_unseen_species_ablation_reduces_to_global_vs_guild(self):
        fit = {
            'global': {'BS_only': .2, 'BT_only': .8},
            'guild': {'g1': {'BS_only': .7, 'BT_only': .3}},
            'taxon': {}, 'source_taxon_guild': {},
            'states': ('BS_only', 'BT_only'),
        }
        p, seen = probabilities(make('9', 's9', 'Z', 'g1', 'BS_only'), fit)
        self.assertFalse(seen)
        self.assertEqual(p[MODELS[2]], fit['global'])
        self.assertEqual(p[MODELS[3]], fit['guild']['g1'])

    def test_heldout_site_section_unit_and_finite_losses(self):
        r = evaluate(self.rows, 'all_observed_three_state')
        self.assertEqual(r['n_source_rail_sections'], 4)
        self.assertEqual(r['n_source_field_locations'], 4)
        self.assertEqual(r['n_source_observed_taxon_site_events'], 8)
        self.assertEqual(len(r['source_true_original_section_scores']), 4)
        self.assertTrue(all(0 < r['primary_logloss_equal_field_site_then_equal_section'][m]
                            for m in MODELS))

    def test_no_fabricated_absences_in_exclusive_secondary(self):
        r = evaluate(self.rows, 'exclusive_BS_or_BT_only')
        self.assertEqual(r['n_source_observed_taxon_site_events'], 6)
        self.assertEqual(r['n_source_field_locations'], 4)

    def test_technical_shuffles_preserve_original_site_and_taxon_states(self):
        old_bins = occupancy_frequency_bins(self.rows)
        old_source = [(x['internal_original_site'], x['source_taxon_EPPO'],
                       x['observed_original_microhabitat_state']) for x in self.rows]
        shuffled = shuffled_guild_rows(self.rows, 20261116)
        new_source = [(x['internal_original_site'], x['source_taxon_EPPO'],
                       x['observed_original_microhabitat_state']) for x in shuffled]
        self.assertEqual(old_source, new_source)
        self.assertEqual(old_bins, occupancy_frequency_bins(shuffled))
        for b in set(old_bins.values()):
            old_g = {x['source_taxon_EPPO']: x['source_four_lifespan_diaspore_guild']
                     for x in self.rows if old_bins[x['source_taxon_EPPO']] == b}
            new_g = {x['source_taxon_EPPO']: x['source_four_lifespan_diaspore_guild']
                     for x in shuffled if old_bins[x['source_taxon_EPPO']] == b}
            self.assertEqual(Counter(old_g.values()), Counter(new_g.values()))

    def test_technically_shuffled_rows_have_one_fixed_label_per_taxon(self):
        new = shuffled_guild_rows(self.rows, 42)
        per_taxon = {}
        for x in new:
            per_taxon.setdefault(x['source_taxon_EPPO'], set()).add(
                x['source_four_lifespan_diaspore_guild'])
        self.assertTrue(all(len(v) == 1 for v in per_taxon.values()))


if __name__ == '__main__':
    unittest.main()

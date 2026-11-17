#!/usr/bin/env python3
"""Synthetic-only tests for v21.1 source stratum bookkeeping."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

from analyze_Slovenia_source_identity_guild_ablation_v210 import evaluate
from analyze_Slovenia_guild_after_identity_seen_novel_v211 import summarize_strata


def rec(site, section, taxon, guild, state):
    return {
        'internal_original_site': site,
        'author_original_section': section,
        'source_taxon_EPPO': taxon,
        'source_four_lifespan_diaspore_guild': guild,
        'observed_original_microhabitat_state': state,
    }


class SourceStrataArithmetic(unittest.TestCase):
    def setUp(self):
        self.rows = [
            rec('1','s1','A','strict_annual_wind','BS_only'),
            rec('1','s1','B','other_known_other','BT_only'),
            rec('2','s2','A','strict_annual_wind','both'),
            rec('2','s2','C','strict_annual_other','BT_only'),
            rec('3','s3','A','strict_annual_wind','BS_only'),
            rec('3','s3','B','other_known_other','BT_only'),
            rec('4','s4','A','strict_annual_wind','both'),
            rec('4','s4','D','other_known_wind','BT_only'),
        ]

    def test_three_state_has_exact_source_support(self):
        result = summarize_strata(evaluate(self.rows, 'all_observed_three_state'))
        self.assertEqual(result['n_true_source_events_seen_in_other_training_sections'], 6)
        self.assertEqual(result['n_true_source_events_novel_taxon_to_training_sections'], 2)
        self.assertEqual(result['strata']['previously_seen']['n_heldout_sections_with_source_support'],4)
        self.assertEqual(result['strata']['not_seen_in_training']['n_heldout_sections_with_source_support'],2)

    def test_binary_strata_are_only_true_observed_exclusive_events(self):
        result = summarize_strata(evaluate(self.rows, 'exclusive_BS_or_BT_only'))
        self.assertEqual(result['n_true_source_events_seen_in_other_training_sections'],4)
        self.assertEqual(result['n_true_source_events_novel_taxon_to_training_sections'],2)
        self.assertAlmostEqual(result['fraction_original_heldout_events_novel_to_train'],1/3)

    def test_no_source_support_does_not_become_fake_zero(self):
        known = [r for r in self.rows if r['source_taxon_EPPO']=='A']
        output = summarize_strata(evaluate(known, 'all_observed_three_state'))
        self.assertFalse(output['strata']['not_seen_in_training']['identified'])
        self.assertTrue(output['strata']['previously_seen']['identified'])

    def test_invalid_source_counts_fail_closed(self):
        data = evaluate(self.rows,'all_observed_three_state')
        data['n_source_observed_taxon_site_events'] += 1
        with self.assertRaises(AssertionError):
            summarize_strata(data)

if __name__=='__main__':
    unittest.main()

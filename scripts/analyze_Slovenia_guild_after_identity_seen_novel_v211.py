#!/usr/bin/env python3
"""v21.1: source-verified Slovenian microhabitat heldout predictions,
separating species previously seen in other training sections from novel taxa.

Input is original author source XLSX SHA-audited via inherited v17.2 parser.
The same v21.0 M0/M1/M2/M3 predictions and original section/site weights
are used. No independent field inference, species adaptation, or new source
data are invented. Do not conflate seen and unseen strata whose source
sampling opportunities differ.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import mean

from analyze_Slovenia_source_identity_guild_ablation_v210 import (
    MODELS, evaluate
)
from analyze_Slovenia_original_guild_vs_species_microhabitat_15section_predict_v172 import (
    author_source_rows, extract_original_xlsx, MODES, N_SECTIONS, N_SITES,
    N_SITE_X_TAXON_INCIDENCES
)

EXPECTED = {
    'all_observed_three_state': {
        'M0_global': 1.0982431044178191,
        'M1_lifespan_diaspore': 1.0855888865022834,
        'M2_species_only': 1.0351061531049255,
        'M3_species_plus_guild': 1.0346230109451737,
    },
    'exclusive_BS_or_BT_only': {
        'M0_global': 0.69205610135632,
        'M1_lifespan_diaspore': 0.6777698225101966,
        'M2_species_only': 0.6497964633268455,
        'M3_species_plus_guild': 0.6482258245641337,
    }
}

STRATA = ('previously_seen', 'not_seen_in_training')


def summarize_strata(one_mode):
    sections = one_mode['source_true_original_section_scores']
    if len(sections) < 2:
        raise ValueError('No heldout section replication')
    out = {}
    for group in STRATA:
        subset = [(x['source_section'], x['seen_taxon_strata'][group])
                  for x in sections if group in x['seen_taxon_strata']]
        if not subset:
            out[group] = {'identified': False, 'reason': 'no source support'}
            continue
        models = {m: mean(x['source_equal_site_logloss'][m] for _, x in subset)
                  for m in MODELS}
        deltas = [{
            'source_section': sec,
            'n_field_sites_with_at_least_one_source_taxon_in_stratum': rec['n_sites'],
            'source_original_delta_M2_minus_M3_nats':
                rec['source_equal_site_logloss'][MODELS[2]] -
                rec['source_equal_site_logloss'][MODELS[3]],
        } for sec, rec in subset]
        out[group] = {
            'identified': True,
            'n_heldout_sections_with_source_support': len(subset),
            'n_source_field_sites_with_source_group_within_sections':
                sum(x['n_sites'] for _, x in subset),
            'original_source_equal_eligible_section_and_site_loss': models,
            'source_M2_species_only_minus_M3_species_plus_guild_nats':
                models[MODELS[2]] - models[MODELS[3]],
            'n_heldout_sections_true_guild_improves_over_species_only':
                sum(x['source_original_delta_M2_minus_M3_nats'] > 0 for x in deltas),
            'n_heldout_sections_true_guild_worsens_over_species_only':
                sum(x['source_original_delta_M2_minus_M3_nats'] < 0 for x in deltas),
            'n_heldout_sections_tied': sum(
                x['source_original_delta_M2_minus_M3_nats'] == 0 for x in deltas),
            'per_original_heldout_section_diagnostic': deltas,
        }
    prior_seen = sum(x['n_taxon_site_seen_on_other_training_section'] for x in sections)
    novel = sum(x['n_taxon_site_unseen_on_other_training_section'] for x in sections)
    if prior_seen + novel != one_mode['n_source_observed_taxon_site_events']:
        raise AssertionError('Author source observed events lost in seen/unseen partition')
    return {
        'n_true_source_events_seen_in_other_training_sections': prior_seen,
        'n_true_source_events_novel_taxon_to_training_sections': novel,
        'fraction_original_heldout_events_novel_to_train':
            novel / (prior_seen + novel),
        'stratified_section_equal_scores_note': (
            'Stratum-specific averages weight only supported original sections; '
            'cannot be arithmetically reconstructed by weighting pooled '
            'stratum means because site/section weights differ.'
        ),
        'strata': out,
        'overall_n_sections_source_M3_better_than_M2':
            one_mode['n_sections_M3_better_than_M2'],
        'overall_n_original_sections': len(sections),
    }


def run_original(original):
    rows, provenance = author_source_rows(original)
    if (len(rows) != N_SITE_X_TAXON_INCIDENCES or
        len({r['internal_original_site'] for r in rows}) != N_SITES or
        len({r['author_original_section'] for r in rows}) != N_SECTIONS):
        raise AssertionError('Original source counts drifted')
    out = {}
    for mode in MODES:
        result = evaluate(rows, mode)
        for m, expected in EXPECTED[mode].items():
            if not math.isclose(
                result['primary_logloss_equal_field_site_then_equal_section'][m],
                expected, abs_tol=1e-10
            ):
                raise AssertionError('v21.0 original source reproduction drift: ' + mode + ' ' + m)
        out[mode] = summarize_strata(result)
    return {
        'version': 'v21.1_original_Slovenia_rail_seen_vs_novel_taxon_guild_after_identity',
        'complete': True,
        'source_provenance': provenance,
        'genuine_original_author_XLSX_SHA_and_v21_0_scores_reproduced': True,
        'modes': out,
        'microhabitat_strata_are_conditional_on_original_taxon_detection': True,
        'not_a_causal_mowing_wind_freight_or_adaptation_result': True,
        'not_independent_countries_or_nonrail_controls': True,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--original-xlsx', required=True)
    p.add_argument('--source-audit102', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    original = extract_original_xlsx(a.original_xlsx, a.source_audit102)
    result = run_original(original)
    dst = Path(a.output)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf8')
    print(json.dumps({
        'version': result['version'],
        'real_source_reproduced': result['genuine_original_author_XLSX_SHA_and_v21_0_scores_reproduced'],
        'modes': {m: {
            'n_seen_site_x_taxon': v['n_true_source_events_seen_in_other_training_sections'],
            'n_novel_site_x_taxon': v['n_true_source_events_novel_taxon_to_training_sections'],
            'novel_fraction': v['fraction_original_heldout_events_novel_to_train'],
            'n_original_sections_M3_better_overall': v['overall_n_sections_source_M3_better_than_M2'],
            'strata': {k: {
                'n_supported_sections': sv.get('n_heldout_sections_with_source_support', 0),
                'source_per_model_loss': sv.get('original_source_equal_eligible_section_and_site_loss'),
                'delta_M2_minus_M3': sv.get('source_M2_species_only_minus_M3_species_plus_guild_nats'),
                'n_sections_M3_better': sv.get('n_heldout_sections_true_guild_improves_over_species_only'),
                'n_sections_M3_worse': sv.get('n_heldout_sections_true_guild_worsens_over_species_only'),
            } for k, sv in v['strata'].items()}
        } for m, v in result['modes'].items()}
    }, ensure_ascii=False))


if __name__ == '__main__':
    main()

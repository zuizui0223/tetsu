#!/usr/bin/env python3
"""v21.0: original Slovenian railway microhabitats, guild-after-identity falsification.

Pre-frozen four models, evaluated by original railway section leaveout:
 M0: global BS_only/BT_only/both rates;
 M1: original author four lifespan x diaspore guilds, shrunk to M0;
 M2: source species identity, shrunk to M0 (NO trait information);
 M3: source species identity, shrunk to M1 (the published v17.2 model).
Original publisher XLSX is SHA checked by the existing v17.2 parser.
The 63 frequency-bin matched taxon-trait label shuffles preserve every
original species-site occurrence and the species occupancy distribution.
They are technical mapping placebos, NOT randomized railway field trials.
Source zone width and surveyed areas differ. No evidence of evolved traits,
train-mediated seed delivery, actual mowing, worldwide convergence, or
causal habitat selection can be inferred from this conditional prediction.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_original_guild_vs_species_microhabitat_15section_predict_v172 import (
    author_source_rows, evaluate_study, extract_original_xlsx, mode_data,
    predict, source_smoothed_model, MODELS as OLD_MODELS, MODES,
    N_SECTIONS, N_SITES, N_JOINT_TRAIT_OBSERVED_TAXA,
    N_SITE_X_TAXON_INCIDENCES, SPECIES_PSEUDO,
)

MODELS = ('M0_global', 'M1_lifespan_diaspore', 'M2_species_only',
          'M3_species_plus_guild')
N_PLACEBO = 63
SEED = 20261116


def probabilities(row, fit):
    """Species-only prior is global M0; it never reads the species guild."""
    old = predict(row, fit)
    count = fit['taxon'].get(row['source_taxon_EPPO'], Counter())
    n = sum(count.values())
    global_p = fit['global']
    species_only = {
        s: (count[s] + SPECIES_PSEUDO * global_p[s]) / (n + SPECIES_PSEUDO)
        for s in fit['states']
    }
    out = dict(zip(MODELS, (
        old[OLD_MODELS[0]], old[OLD_MODELS[1]], species_only,
        old[OLD_MODELS[2]]
    )))
    if any(not math.isclose(sum(p.values()), 1, abs_tol=1e-10) or
           min(p.values()) <= 0 for p in out.values()):
        raise AssertionError('Invalid source probabilities')
    return out, bool(old['original_heldout_species_has_prior_original_train_observations'])


def evaluate(rows, mode):
    """Equal weight to field sites within heldout sections, and to sections."""
    observed, states = mode_data(rows, mode)
    sections = sorted({r['author_original_section'] for r in rows})
    if len(sections) < 2:
        raise ValueError('Need multiple genuine source railway sections')
    by_section = []
    for sec in sections:
        train = [r for r in observed if r['author_original_section'] != sec]
        test = [r for r in observed if r['author_original_section'] == sec]
        if not train or not test:
            raise ValueError('No common original source section support')
        fit = source_smoothed_model(train, states)
        sites = defaultdict(list)
        seen = {True: defaultdict(list), False: defaultdict(list)}
        for r in test:
            p, prior = probabilities(r, fit)
            state = r['observed_original_microhabitat_state']
            vals = {m: -math.log(p[m][state]) for m in MODELS}
            site = r['internal_original_site']
            sites[site].append(vals)
            seen[prior][site].append(vals)
        def site_mean(group):
            return {m: mean(mean(x[m] for x in v) for v in group.values())
                    for m in MODELS}
        section = {
            'source_section': sec, 'n_original_sites': len(sites),
            'n_source_site_taxon_events': sum(map(len, sites.values())),
            'source_equal_site_logloss': site_mean(sites),
            'n_taxon_site_seen_on_other_training_section': sum(
                len(v) for v in seen[True].values()),
            'n_taxon_site_unseen_on_other_training_section': sum(
                len(v) for v in seen[False].values()),
        }
        section['seen_taxon_strata'] = {
            ('previously_seen' if key else 'not_seen_in_training'):
                {'n_sites': len(v), 'source_equal_site_logloss': site_mean(v)}
            for key, v in seen.items() if v
        }
        by_section.append(section)
    avg = {m: mean(x['source_equal_site_logloss'][m] for x in by_section)
           for m in MODELS}
    return {
        'source_original_response_mode': mode,
        'n_source_rail_sections': len(sections),
        'n_source_field_locations': len({r['internal_original_site'] for r in observed}),
        'n_source_observed_taxon_site_events': len(observed),
        'primary_logloss_equal_field_site_then_equal_section': avg,
        'M2_species_only_minus_M3_species_plus_guild_nats':
            avg[MODELS[2]] - avg[MODELS[3]],
        'M0_global_minus_M1_guild_nats': avg[MODELS[0]] - avg[MODELS[1]],
        'M0_global_minus_M2_species_only_nats': avg[MODELS[0]] - avg[MODELS[2]],
        'n_sections_M3_better_than_M2': sum(
            x['source_equal_site_logloss'][MODELS[3]] <
            x['source_equal_site_logloss'][MODELS[2]] for x in by_section),
        'source_true_original_section_scores': by_section,
    }


def occupancy_frequency_bins(rows):
    count = Counter(r['source_taxon_EPPO'] for r in rows)
    def bin_of(n):
        return ('01' if n == 1 else '02_03' if n <= 3 else
                '04_07' if n <= 7 else '08_15' if n <= 15 else '16_plus')
    return {taxon: bin_of(n) for taxon, n in count.items()}


def shuffled_guild_rows(rows, seed):
    """Exactly one shuffled taxon label everywhere; preserve occurrence masks."""
    gb = {}
    for r in rows:
        t, g = r['source_taxon_EPPO'], r['source_four_lifespan_diaspore_guild']
        if t in gb and gb[t] != g:
            raise ValueError('Original taxon changed author trait guild')
        gb[t] = g
    bins = occupancy_frequency_bins(rows)
    grouped = defaultdict(list)
    for t in sorted(gb):
        grouped[bins[t]].append(t)
    rng = random.Random(seed)
    perm = {}
    for b, taxa in sorted(grouped.items()):
        traits = [gb[t] for t in taxa]
        rng.shuffle(traits)
        perm.update(zip(taxa, traits))
        if Counter(perm[t] for t in taxa) != Counter(gb[t] for t in taxa):
            raise AssertionError('Trait inventory changed in occupancy bin')
    out = []
    for r in rows:
        x = r.copy()
        x['source_four_lifespan_diaspore_guild'] = perm[r['source_taxon_EPPO']]
        out.append(x)
    return out


def source_full_run(original):
    rows, provenance = author_source_rows(original)
    if (len(rows) != N_SITE_X_TAXON_INCIDENCES or
        len({r['internal_original_site'] for r in rows}) != N_SITES or
        len({r['author_original_section'] for r in rows}) != N_SECTIONS or
        len({r['source_taxon_EPPO'] for r in rows}) != N_JOINT_TRAIT_OBSERVED_TAXA):
        raise AssertionError('Original Slovenia source data support drift')
    modes = {}
    for mode in MODES:
        actual = evaluate(rows, mode)
        reproduced = evaluate_study(rows, mode)
        old = reproduced['original_source_study_holdout_scores'][
            'out_of_original_source_section_log_loss'][
            'original_source_section_equal_site_equal_section']
        current = actual['primary_logloss_equal_field_site_then_equal_section']
        for m, prior in zip((MODELS[0], MODELS[1], MODELS[3]), OLD_MODELS):
            if not math.isclose(current[m], old[prior], abs_tol=1e-10):
                raise AssertionError('Previously source-verified v17.2 score drift: ' + m)
        placebos = []
        for index in range(N_PLACEBO):
            shuffled = shuffled_guild_rows(rows, SEED + 1009 * (index + 1))
            placebos.append(evaluate(shuffled, mode)[
                'M2_species_only_minus_M3_species_plus_guild_nats'])
        actual_gain = actual['M2_species_only_minus_M3_species_plus_guild_nats']
        sorted_null = sorted(placebos)
        actual['source_real_vs_occupancy_bin_matched_guild_label_placebo'] = {
            'n_predeclared_technical_trait_remappings': N_PLACEBO,
            'true_original_incremental_guild_gain_after_identity_nats': actual_gain,
            'placebo_min': min(sorted_null),
            'placebo_median': sorted_null[len(sorted_null) // 2],
            'placebo_max': max(sorted_null),
            'n_placebo_gains_ge_original': sum(x >= actual_gain for x in placebos),
            'empirical_rank_fraction_for_mapping_diagnostic_ONLY_not_field_p_value':
                (1 + sum(x >= actual_gain for x in placebos)) / (N_PLACEBO + 1),
            'original_species_presence_and_source_state_data_fixed_in_placebos': True,
            'original_lifespan_diaspore_four_category_counts_fixed_within_occupancy_bins': True,
            'unmeasured_trait_phylogenetic_structure_not_preserved': True,
            'source_taxon_present_elsewhere_is_not_randomized_field_treatment': True,
        }
        modes[mode] = actual
    return {
        'version': 'v21.0_original_Slovenia_species_only_vs_traits_after_species_ablation',
        'complete': True,
        'genuine_original_publisher_source': provenance,
        'v17_2_genuine_source_leave_section_out_scores_reproduced': True,
        'predeclared_source_placebo_count': N_PLACEBO,
        'modes': modes,
        'no_independent_nonrail_or_road_controls': True,
        'no_transport_seed_mowing_wind_genetic_adaptation_or_global_homogenization_estimated': True,
        'site_species_events_not_independent_field_replicates': True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--original-xlsx', required=True)
    parser.add_argument('--source-audit102', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    original = extract_original_xlsx(args.original_xlsx, args.source_audit102)
    output = source_full_run(original)
    p = Path(args.output)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding='utf8')
    print(json.dumps({
        'version': output['version'],
        'original_source_reproduced': output['v17_2_genuine_source_leave_section_out_scores_reproduced'],
        'main': {name: {
            'source_original_eligible_events': x['n_source_observed_taxon_site_events'],
            'source_original_sections': x['n_source_rail_sections'],
            'source_new_M0_M1_M2_M3_equal_section_site_logloss':
                x['primary_logloss_equal_field_site_then_equal_section'],
            'source_real_guide_after_identity_nats':
                x['M2_species_only_minus_M3_species_plus_guild_nats'],
            'source_technical_placebo':
                x['source_real_vs_occupancy_bin_matched_guild_label_placebo'],
        } for name, x in output['modes'].items()}
    }, ensure_ascii=False))


if __name__ == '__main__':
    main()

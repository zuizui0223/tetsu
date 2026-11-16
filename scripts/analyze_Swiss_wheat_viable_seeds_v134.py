#!/usr/bin/env python3
"""v13.4: imported-wheat contaminant seeds -> germination -> genotype.

Source: Schulze, Brodmann, Oehen & Bagutti (2015), Table 2,
DOI 10.1007/s11356-015-4903-y. Secondary analysis of published data.
Two mills, three pooled collection dates per mill, one germination
tray/group per mill. Seeds/seedlings are not independent mill replicates.
No field establishment or train-specific loss probability is estimated.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

DOI = '10.1007/s11356-015-4903-y'
SOURCE_PDF = ('https://www.cb.iee.unibe.ch/e58878/e337393/e337414/'
              'e483217/e615152/Schultzeetal2015_eng.pdf')
FIELDS = ['mill_id', 'seeds_sown', 'seedlings', 'GT73', 'MS8_RF3',
          'MS8', 'RF3', 'nonGM', 'seed_pools', 'pool_size',
          'GT73_positive_pools', 'MS8_positive_pools', 'RF3_positive_pools']
GENOTYPES = ['GT73', 'MS8_RF3', 'MS8', 'RF3', 'nonGM']


def validate(rows: list[dict]) -> None:
    if len(rows) != 2 or {r['mill_id'] for r in rows} != {1, 2}:
        raise ValueError('Expected two distinct published mill groups')
    for r in rows:
        if any(type(r[k]) is not int or r[k] < 0 for k in FIELDS):
            raise ValueError('Counts must be nonnegative integers')
        if r['seeds_sown'] != 150 or r['seedlings'] > r['seeds_sown']:
            raise ValueError('Invalid original sown/germinated denominator')
        if sum(r[k] for k in GENOTYPES) != r['seedlings']:
            raise ValueError('Mutually exclusive genotype counts must sum to seedlings')
        if r['seedlings'] == 0 or r['pool_size'] != 10 or r['seed_pools'] != 10:
            raise ValueError('Original germination or pooled-test sample frame differs')
        if any(r[k] > r['seed_pools'] for k in FIELDS if k.endswith('_positive_pools')):
            raise ValueError('Positive pools exceed sampled pools')


def read_data(path: str | Path) -> list[dict]:
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != FIELDS:
            raise ValueError('Source Table 2 CSV schema differs')
        rows = [{k: int(v) for k, v in r.items()} for r in reader]
    validate(rows)
    return sorted(rows, key=lambda r: r['mill_id'])


def verify_table_text(text: str, rows: list[dict]) -> dict:
    """Match every numerical Table 2 row, not merely grand totals."""
    if DOI not in text:
        raise ValueError('Source DOI absent from PDF text')
    caption = re.search(r'Table\s+2\s+a\s+Numbers\s+of\s+seedlings', text, re.I)
    if caption is None:
        raise ValueError('Original Table 2 caption not found')
    table = text[caption.start():caption.start()+11000]
    for r in rows:
        gm = sum(r[k] for k in GENOTYPES[:-1])
        a = [r['mill_id'], r['seeds_sown'], r['seedlings'],
             round(100*r['seedlings']/r['seeds_sown'])]
        a += [r[k] for k in GENOTYPES]
        a += [round(100*gm/r['seedlings'])]
        b = [r['mill_id'], r['seed_pools'], r['GT73_positive_pools'],
             r['MS8_positive_pools'], r['RF3_positive_pools']]
        for label, values in [('2a', a), ('2b', b)]:
            pattern = r'^\s*' + r'\s+'.join(map(str, values)) + r'\s*$'
            if len(re.findall(pattern, table, flags=re.M)) != 1:
                raise ValueError(f'Table {label}: mill {r["mill_id"]} exact row not verified')
    return {'rows_verified': 4, 'table_2a_numeric_cells_verified': 20,
            'table_2b_numeric_cells_verified': 10, 'original_DOI': DOI}


def verify_pdf(path: str | Path, rows: list[dict]) -> dict:
    path = Path(path)
    raw = path.read_bytes()
    if not raw.startswith(b'%PDF-') or len(raw) < 10000:
        raise ValueError('Downloaded source is not a substantive PDF')
    proc = subprocess.run(['pdftotext', '-layout', str(path), '-'],
                          check=True, capture_output=True, timeout=90)
    proof = verify_table_text(proc.stdout.decode('utf-8'), rows)
    return {**proof, 'pdf_sha256': hashlib.sha256(raw).hexdigest(),
            'pdf_bytes': len(raw), 'native_text_not_OCR': True,
            'source_table_checked_against_PDF': True}


def field_scenario(seeds: float, laboratory_gm_yield: float,
                   spill_fraction: float | None,
                   establishment_given_spill: float | None) -> float | None:
    """User-specified scenario only; an unknown transition stays unknown."""
    for p in (laboratory_gm_yield, spill_fraction, establishment_given_spill):
        if p is not None and (not math.isfinite(p) or not 0 <= p <= 1):
            raise ValueError('Scenario probabilities must be finite and within [0,1]')
    if not math.isfinite(seeds) or seeds < 0:
        raise ValueError('Seed input must be finite and nonnegative')
    if spill_fraction is None or establishment_given_spill is None:
        return None
    return seeds * laboratory_gm_yield * spill_fraction * establishment_given_spill


def analyze(rows: list[dict]) -> dict:
    validate(rows)
    mills = []
    for r in sorted(rows, key=lambda x: x['mill_id']):
        gm = r['seedlings'] - r['nonGM']
        mills.append({'mill_id': r['mill_id'], 'sown': r['seeds_sown'],
                      'germinated': r['seedlings'], 'GM_seedlings': gm,
                      'germinated_fraction': r['seedlings']/r['seeds_sown'],
                      'GM_fraction_of_seedlings': gm/r['seedlings'],
                      'GM_seedling_per_sown_seed': gm/r['seeds_sown'],
                      'exclusive_genotype_counts': {k: r[k] for k in GENOTYPES}})
    n = sum(r['seeds_sown'] for r in rows)
    germinated = sum(r['seedlings'] for r in rows)
    gm = sum(r['seedlings']-r['nonGM'] for r in rows)
    # Historical inputs and extrapolation assumptions from the 2015 paper;
    # the mass impurity estimate was available for mill 1, not both mills.
    tonnes, impurity_fraction, grams_per_1000 = 77512, 0.00005, 3.2
    contaminant_mass = tonnes * impurity_fraction
    seeds = contaminant_mass * 1_000_000 / (grams_per_1000/1000)
    projected = {str(r['mill_id']): seeds*r['GM_seedling_per_sown_seed'] for r in mills}
    g = [r['germinated_fraction'] for r in mills]
    h = [r['GM_fraction_of_seedlings'] for r in mills]
    return {
        'complete': True, 'version': 'v13.4', 'source_DOI': DOI,
        'by_mill': mills,
        'pooled_descriptive': {'sown': n, 'germinated': germinated,
                              'GM_seedlings': gm,
                              'germinated_fraction': germinated/n,
                              'GM_fraction_of_seedlings': gm/germinated,
                              'GM_seedling_per_sown_seed': gm/n,
                              'independent_mills': 2,
                              'pooled_collection_dates_per_mill': 3,
                              'not_a_multi_study_meta_analysis': True},
        'pool_tests': {'total_pools': sum(r['seed_pools'] for r in rows),
                       'seeds_per_pool': 10,
                       'GT73_positive_pools': sum(r['GT73_positive_pools'] for r in rows),
                       'individual_seed_GT73_prevalence': None,
                       'joint_MS8_RF3_genotype_inferred_from_pools': False},
        'historical_conditional_projection': {
            'historical_import_window': '2010-2013', 'sampling_year': 2014,
            'wheat_import_tonnes_per_year': tonnes,
            'OSR_impurity_mass_fraction': impurity_fraction,
            'impurity_estimated_from_only_mill_1': True,
            'assume_same_impurity_in_common_and_durum_wheat': True,
            'thousand_seed_mass_g': grams_per_1000,
            'contaminant_mass_tonnes': contaminant_mass,
            'all_contaminant_seeds': seeds,
            'GM_seedling_equivalent_by_mill': projected,
            'pooled_descriptive_GM_seedling_equivalent': seeds*gm/n,
            'cross_mill_marginal_extrema_scenario': [seeds*min(g)*min(h), seeds*max(g)*max(h)],
            'two_mill_range_is_confidence_interval': False,
            'not_actual_track_deposition_or_current_2026_flow': True,
            'rail_spill_count': None, 'established_plants': None},
        'claims': {'cargo_contaminant_germination_observed': True,
                   'GM_genotype_measured_in_germinated_seedlings': True,
                   'field_establishment_probability_measured': False,
                   'individual_cargo_to_field_parentage_proven': False,
                   'train_specific_loss_fraction_measured': False,
                   'GM_relative_germination_advantage_identified': False,
                   'native_flora_mowing_or_wind_adaptation_tested': False},
        'limits': [
            'Two pooled mill groups with different sowing dates, not 300 independent shipments.',
            'Seeds were selected from an enriched cleaning-waste fraction; mass scaling is essential.',
            'Laboratory dormancy-breaking and four-week emergence do not give ballast establishment rates.',
            'Non-germinated seed genotypes were not individually measured; genotype-specific viability cannot be estimated.',
            'Commercial event-type overlap with feral plants supports a route but is not parentage identification.',
            '2013 field plants and 2014 shipment samples are not one prospectively followed cohort.',
            'National input scenarios are not spill measurements, confidence intervals or multi-country estimates.']}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--csv', required=True)
    ap.add_argument('--pdf')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    rows = read_data(args.csv)
    result = analyze(rows)
    result['source_verification'] = verify_pdf(args.pdf, rows) if args.pdf else {
        'source_table_checked_against_PDF': False,
        'status': 'CSV transcribed from retrieved primary Table 2; PDF-byte audit not performed in this run'}
    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False))

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""v13.5: original Basel 2015 report, route monitoring and cargo viability.

Sections 2.5.6 and 2.5.7 are different observation frames. The cargo
experiment repeats the Schulze 2015 study; it is not another experiment.
No route-specific arrival-to-establishment probability is estimated.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

SOURCE_URL = ('https://media.bs.ch/original_file/'
              '0d096108b159f34916e7b5203adfab20e885ac1a/jahresbericht2015.pdf')


def section(text: str, number: str, following: str) -> str:
    # Use last body heading, not the earlier table of contents.
    starts = list(re.finditer(re.escape(number) + r'\s+', text))
    if not starts:
        raise ValueError('Original report section missing: ' + number)
    start = starts[-1].end()
    end = re.search(re.escape(following) + r'\s+', text[start:])
    if not end:
        raise ValueError('Original next report section missing: ' + following)
    return text[start:start + end.start()]


def analyze_text(text: str) -> dict:
    text = re.sub(r'\s+', ' ', text.replace('\u00ad', '').replace('\x02', '-'))
    route = section(text, '2.5.6', '2.5.7')
    cargo = section(text, '2.5.7', '2.5.8')
    anchors = {
        'route_sample_count': (route, r'\b1314\s+Rapsproben\b'),
        'targeted_negatives': (route, r'Alle\s+Proben\s+wurden\s+GV-negativ\s+getestet'),
        'other_sample_count': (route, r'\b391\s+Rapsproben\b'),
        'other_survey_distance': (route, r'\b30\s*km\b'),
        'positive_cluster': (route, r'sieben\s+Rapsproben\s+an\s+einem\s+Standort\s+in\s+Giornico'),
        'marker': (route, r'\bGT73\b'),
        'mills': (route, r'Jowa.*Swissmill'),
        'mill_sown_seeds': (cargo, r'\b150\s+Rapssamen\b'),
        'germination_aggregate': (cargo, r'\b25\s*%'),
        'wheat_impurity_mass_fraction': (cargo, r'\b0[.,]005\s*%'),
        'GM_event': (cargo, r'\bMS8xRF3\b'),
        'gene_source': (cargo, r'\bGT73\b'),
        'approximate_input': (cargo, r'zwei\s+keimf[aä]hige\s+transgene\s+Rapssamen\s+pro\s+kg'),
    }
    failed = [k for k, (s, pat) in anchors.items() if re.search(pat, s, re.I) is None]
    if failed:
        raise ValueError('Primary government source anchors differ: ' + ', '.join(failed))
    return {
        'complete': True, 'version': 'v13.5', 'government_report_year': 2015,
        'source_url': SOURCE_URL, 'source_sections': ['2.5.6', '2.5.7'],
        'source_numeric_and_context_anchors_verified': sorted(anchors),
        'route_monitoring': {
            'targeted_route': {
                'frame': 'Kleinhuningen port to Jowa Wildegg and Swissmill Zurich',
                'tested_samples': 1314, 'positive_samples': 0,
                'positive_fraction': 0.0, 'total_locations': None,
                'number_of_independent_plants': None,
                'interpretation': 'No GM detected in tested rape samples; not zero seed arrival.'},
            'other_survey': {
                'frame': 'Non-risk-based survey of 30 km of the SBB network',
                'surveyed_length_km': 30, 'tested_samples': 391,
                'positive_samples': 7, 'positive_fraction': 7/391,
                'positive_locations': 1, 'positive_locality': 'Giornico (TI)',
                'marker': 'GT73', 'total_locations': None}},
        'cargo_experiment_corroboration': {
            'sampling_year': 2014, 'mills': 2, 'sown_seeds_per_mill': 150,
            'reported_mean_germination_fraction': 0.25,
            'reported_impurity_mass_fraction_mill_1': 0.00005,
            'reported_approximate_germinable_GM_seeds_per_kg_uncleaned_wheat': 2,
            'same_journal_study_DOI': '10.1007/s11356-015-4903-y'},
        'same_cargo_study_republished_not_independent_replication': True,
        'journal_Table2_cell_audit_performed': False,
        'field_establishment_rate_identified': False,
        'rail_spill_rate_identified': False,
        'randomized_route_comparison': False,
        'limits': [
            'Cargo collection in 2014 and 2015 route monitoring are not the same seed cohort.',
            '7 positive samples cluster at one location; no independent-sample causal test is justified.',
            'Total field sites, unique plant count, detection probability and individual exposure are unreported here.',
            'The laboratory study is corroborated by the same institution, not an independent replication.',
            'No-gene detection does not prove no cargo spills, no non-GM recruitment or zero GM population.',
            'This annual report does not reproduce every mill-specific Table 2 cell of the journal article.']}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--government-pdf', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    pdf = Path(args.government_pdf)
    raw = pdf.read_bytes()
    if not raw.startswith(b'%PDF-') or len(raw) < 100000:
        raise ValueError('Government source is not a substantive PDF')
    proc = subprocess.run(['pdftotext', '-layout', str(pdf), '-'],
                          check=True, capture_output=True, timeout=120)
    result = analyze_text(proc.stdout.decode('utf-8'))
    result['source_verification'] = {
        'actual_government_pdf_bytes': len(raw),
        'actual_government_pdf_sha256': hashlib.sha256(raw).hexdigest(),
        'native_text_not_OCR': True,
        'primary_report_context_and_numeric_anchors_checked': True}
    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False))

if __name__ == '__main__':
    main()

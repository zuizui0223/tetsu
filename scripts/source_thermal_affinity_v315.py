"""Read fixed external temperature indicators and the baseline Kashmir pairs.

This adapter returns species-level maps and in-memory occurrence sets only.
It does not calculate community means, habitat contrasts, or climate effects.
Original trait workbooks and occurrence matrices remain outside the repository.
"""

from collections import Counter, defaultdict
import hashlib
import math
from pathlib import Path
import re


PLAN_COMMIT = '12688cc91aab2a3844c47d512f3ce2b90e9086b8'
TRAIT_SOURCES = {
    'tichy': {
        'filename': 'tichy_2023_zenodo_v2.xlsx',
        'published_filename': 'Indicator.values-tables-2022-11-07-Zenodo.v2.xlsx',
        'sha256': 'b3fa7c3a6e9d356e52a45472549f807967bd0c3189e58e7e8a22e0389c29cc4e',
        'url': 'https://zenodo.org/api/records/7427088/files/Indicator.values-tables-2022-11-07-Zenodo.v2.xlsx/content',
        'dataset_doi': '10.5281/zenodo.7427088',
        'publication_doi': '10.1111/jvs.13168',
        'license': 'CC BY 4.0',
        'sheet': 'Tab-IVs-Tichy-et-al2023',
        'limits': [1, 12],
        'source_rows': 8908,
        'source_status_counts': {'numeric': 7402, 'indifferent': 1147, 'missing': 359},
    },
    'eive': {
        'filename': 'EIVE_Paper_1.0_SM_08.xlsx',
        'published_filename': 'EIVE_Paper_1.0_SM_08.xlsx',
        'sha256': '9551b585ba13dfe3fddbd2287f65fc98ac9f4ec5f80b9ae6f6b985bed73fb4a6',
        'url': 'https://zenodo.org/api/records/7534792/files/EIVE_Paper_1.0_SM_08.xlsx/content',
        'dataset_doi': '10.5281/zenodo.7534792',
        'publication_doi': '10.3897/VCS.98324',
        'license': 'CC BY 4.0',
        'sheet': 'mainTable',
        'limits': [0, 10],
        'source_rows': 14835,
        'source_status_counts': {'numeric': 14496, 'missing': 339},
    },
}
TICHY_HEADER_ROWS = (
    (None, None, 'LIGHT', 'TEMPERATURE', 'MOISTURE', 'REACTION', 'NUTRIENTS', 'SALINITY'),
    ('SeqID', 'Taxon', 'Average', 'Average', 'Average', 'Average', 'Average', 'Average'),
)
EIVE_COLUMNS = (
    'TaxonConcept', 'UUID', 'TaxonRank', 'AccordingTo',
    'EIVEres-M', 'EIVEres-M.nw3', 'EIVEres-M.n',
    'EIVEres-N', 'EIVEres-N.nw3', 'EIVEres-N.n',
    'EIVEres-R', 'EIVEres-R.nw3', 'EIVEres-R.n',
    'EIVEres-L', 'EIVEres-L.nw3', 'EIVEres-L.n',
    'EIVEres-T', 'EIVEres-T.nw3', 'EIVEres-T.n',
)
EIVE_RANK_COUNTS = {
    'Species': 11148, 'Subspecies': 2899, 'Species Aggregate': 664,
    'Species Aggregate s. l.': 60, 'Variety': 42, 'Section': 22,
}
KASHMIR_ZONES = ('Track', 'Embankment', 'Interior Plot 1', 'Interior Plot 2')


def _require_normalizer(normalize):
    if not callable(normalize):
        raise TypeError('A species-name normalizer is required')


def _strict_trait_key(label, normalize):
    """Exclude explicit non-species concepts before exact name normalization."""
    if not isinstance(label, str):
        return None
    label = ' '.join(label.split())
    if re.search(r'\b(?:agg|aggr|aggregate|group|complex|sect|section|cf|aff|'
                 r'subsp|ssp|var|forma)\b|(?:^|\s)f\.(?=\s|$)', label, re.I):
        return None
    if (re.search(r'\bs\.\s*(?:l|lat)\.', label, re.I)
            or 'sensu lato' in label.lower() or '×' in label
            or re.search(r'(?:^|\s)x(?:\s|$)', label, re.I)):
        return None
    key = normalize(label)
    if not isinstance(key, str) or not re.fullmatch(r'[A-Z][A-Za-z-]+ [a-z][a-z-]+', key):
        return None
    return key


def _temperature_status(value, limits, allow_indifferent):
    if value is None or value == 'NA':
        return 'missing', None
    if allow_indifferent and value == 'x':
        return 'indifferent', None
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or not limits[0] <= value <= limits[1]):
        raise ValueError('Invalid temperature indicator encoding or value outside its scale')
    return 'numeric', float(value)


def _parse_trait_rows(rows, normalize, definition):
    _require_normalizer(normalize)
    limits = TRAIT_SOURCES[definition]['limits']
    candidates = defaultdict(list)
    source_counts, rank_counts = Counter(), Counter()
    source_rows = excluded_rows = 0
    for label, raw_temperature, rank in rows:
        source_rows += 1
        status = _temperature_status(raw_temperature, limits, definition == 'tichy')
        source_counts[status[0]] += 1
        if definition == 'eive':
            rank_counts[rank] += 1
        key = _strict_trait_key(label, normalize)
        if key is None or (definition == 'eive' and rank != 'Species'):
            excluded_rows += 1
            continue
        candidates[key].append(status)

    values, statuses = {}, {}
    for key, entries in sorted(candidates.items()):
        unique = set(entries)
        if len(unique) != 1:
            statuses[key] = 'conflict'
            continue
        status, value = next(iter(unique))
        statuses[key] = status
        if status == 'numeric':
            values[key] = value
    audit = {
        'source_rows': source_rows,
        'source_status_counts': dict(sorted(source_counts.items())),
        'excluded_taxon_rows': excluded_rows,
        'eligible_species_rows': source_rows - excluded_rows,
        'normalized_keys': len(candidates),
        'duplicate_normalized_keys': sum(len(entries) > 1 for entries in candidates.values()),
        'conflicting_normalized_keys': sum(status == 'conflict' for status in statuses.values()),
        'retained_status_counts': dict(sorted(Counter(statuses.values()).items())),
    }
    if definition == 'eive':
        audit['source_rank_counts'] = dict(rank_counts)
    return {'values': values, 'statuses': statuses, 'limits': list(limits), 'audit': audit}


def parse_tichy_rows(rows, normalize):
    """Read (taxon, temperature) pairs; retain x separately from missing NA."""
    return _parse_trait_rows(((label, value, None) for label, value in rows),
                             normalize, 'tichy')


def parse_eive_rows(rows, normalize):
    """Read (taxon, rank, temperature) triples without collapsing lower ranks."""
    return _parse_trait_rows(((label, value, rank) for label, rank, value in rows),
                             normalize, 'eive')


def _read_pinned_trait(directory, definition, normalize):
    descriptor = TRAIT_SOURCES[definition]
    path = Path(directory) / descriptor['filename']
    with path.open('rb') as stream:
        checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
    if checksum != descriptor['sha256']:
        raise ValueError(f'Original {definition} trait workbook SHA-256 mismatch')

    from openpyxl import load_workbook

    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        if descriptor['sheet'] not in workbook.sheetnames:
            raise ValueError(f'Original {definition} temperature sheet missing')
        sheet = workbook[descriptor['sheet']]
        expected_shape = (8910, 8) if definition == 'tichy' else (14836, 19)
        if (sheet.max_row, sheet.max_column) != expected_shape:
            raise ValueError(f'Original {definition} temperature table shape changed')
        iterator = sheet.iter_rows(values_only=True)
        if definition == 'tichy':
            if tuple(next(iterator) for _ in range(2)) != TICHY_HEADER_ROWS:
                raise ValueError('Original Tichy temperature column headers changed')
            rows = list(iterator)
            if [row[0] for row in rows] != list(range(1, 8909)):
                raise ValueError('Original Tichy sequential taxon rows changed')
            trait = parse_tichy_rows(((row[1], row[3]) for row in rows), normalize)
        else:
            if next(iterator) != EIVE_COLUMNS:
                raise ValueError('Original EIVE temperature column headers changed')
            trait = parse_eive_rows(((row[0], row[2], row[16]) for row in iterator), normalize)
    finally:
        workbook.close()
    audit = trait['audit']
    if (audit['source_rows'] != descriptor['source_rows']
            or audit['source_status_counts'] != descriptor['source_status_counts']):
        raise ValueError(f'Original {definition} temperature status census changed')
    if definition == 'eive' and audit['source_rank_counts'] != EIVE_RANK_COUNTS:
        raise ValueError('Original EIVE taxonomic rank census changed')
    metadata = {
        **descriptor,
        'bytes': path.stat().st_size,
        'sha256_verified': True,
        'source_rows': audit['source_rows'],
        'source_status_counts': dict(audit['source_status_counts']),
    }
    return trait, metadata


def load_traits(directory, normalize):
    """Return two separate, fixed-scale thermal maps with aggregate provenance."""
    _require_normalizer(normalize)
    result, metadata = {}, {}
    for definition in ('tichy', 'eive'):
        result[definition], metadata[definition] = _read_pinned_trait(
            directory, definition, normalize)
    result['metadata'] = {
        'design_commit': PLAN_COMMIT,
        'sources': metadata,
        'matching': 'Strict species names; no inferred synonyms or infraspecific averaging',
        'duplicate_policy': 'Identical values/statuses collapse; conflicts remain unknown',
        'scale_policy': 'One fixed value per species and definition, shared across habitats',
        'tichy_x_policy': 'Indifferent; never converted to zero, missing, or a midpoint',
        'eive_policy': 'Separate consensus index retaining the original regional-indifference treatment',
    }
    return result


def pairs_from_kashmir_plots(plots, headers, normalize, expected_stations=16):
    """Normalize complete 2014 station unions and apply shared ambiguity masks.

The output contains sets only. Original station identifiers remain internal.
"""
    _require_normalizer(normalize)
    if any(not isinstance(label, str) or not label.strip() for label in headers):
        raise ValueError('Invalid original Kashmir species column label')
    stations = sorted({site for site, year, _ in plots
                       if year == 2014 and not site.startswith('Between ')})
    if len(stations) != expected_stations:
        raise ValueError('Unexpected baseline Kashmir station count')
    raw_pairs = []
    for station in stations:
        observed = [plots.get((station, 2014, zone)) for zone in KASHMIR_ZONES]
        if any(indices is None for indices in observed):
            raise ValueError('Incomplete baseline Kashmir station habitat pair')
        if any(isinstance(index, bool) or not isinstance(index, int)
               or not 0 <= index < len(headers)
               for indices in observed for index in indices):
            raise ValueError('Invalid original Kashmir species column index')
        raw_pairs.append(({headers[index] for index in observed[0] | observed[1]},
                          {headers[index] for index in observed[2] | observed[3]}))

    from source_matched_rail_road_v314 import audit_names, unresolved_keys

    labels = set().union(*(side for pair in raw_pairs for side in pair))
    audit = audit_names(labels, normalize)
    unavailable = unresolved_keys(audit, labels, normalize)
    return [(set(audit[label] for label in rail) - {None} - unavailable,
             set(audit[label] for label in interior) - {None} - unavailable)
            for rail, interior in raw_pairs]


def kashmir_pairs(data, normalize):
    """Read the same 16 baseline stations as v314, with no exported IDs."""
    _require_normalizer(normalize)
    from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
    from openpyxl import load_workbook

    path = Path(data) / 'Matrix_plot.xlsx'
    plots, _ = parse_xlsx_matrix(path)
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = workbook.worksheets[0]
        if (sheet.max_row, sheet.max_column) != (249, 129):
            raise ValueError('Original Kashmir workbook shape changed')
        rows = sheet.iter_rows(values_only=True)
        headers = next(rows)[6:]
        selected = 0
        for row in rows:
            site = ' '.join(str(row[0]).split())
            if row[4] != 2014 or site.startswith('Between '):
                continue
            selected += 1
            key = (site, 2014, ' '.join(str(row[5]).split()))
            values = row[6:]
            if (any(isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or value < 0 for value in values)
                    or plots.get(key) != {index for index, value in enumerate(values) if value > 0}):
                raise ValueError('Kashmir XML and worksheet baseline observations disagree')
        if selected != 64:
            raise ValueError('Expected four original quadrats at each of 16 baseline stations')
    finally:
        workbook.close()
    return pairs_from_kashmir_plots(plots, headers, normalize)

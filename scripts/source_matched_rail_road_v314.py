"""Read the pinned v314 sources and aggregate identifiable name records in memory."""

import hashlib
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


SI_XLSX = 'si/source_original/SLO_2026_rail_vegetation_microhabitats_original.xlsx'
SI_JSON = 'si/results/audit/Slovenia_2026_RAIL_92sites_plant_supplement_source_v102.json'
SOURCE_SHA256 = {
    'cz1.xlsx': 'df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06',
    'cz2.xlsx': '00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2',
    'parma.pdf': 'c09bd96ff7307a936030195cc4510bf2591576a713de54019625c04aa9333787',
    'hrba.pdf': '90d489ca8cbc29464d3e6a2b1d8b5d9568c9940e953b1ce720964991366ab956',
    SI_XLSX: '26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f',
    SI_JSON: 'c3973ab0bc55740de019b470e1a2768562267d5df3a54c1916962ca112b603d1',
    'sk30_2.pdf': 'ff8271a396ecda01d92bba9817721577bcdd44fdac1309cc4884240c9e3c9838',
    'Matrix_plot.xlsx': '529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474',
    'miren.csv': 'e966419c14879861198c6c1e01c49ee0eca00c3905d784102d6e2c229ac56169',
    'tree1.xlsx': '6ce8e4ea068f644e9ad08eab97a4352e88c071fb46bea038d151f792eb8a21de',
    'sk_urban.bin': '5c70e4b02a17b2886cafb0dde5221a6e7ad440780f55a51249b2fb95ff6d2f22',
    'wa_rail_reserve.pdf': '5e392a5ff5f705fa01e728d19e8a0b8257dfb62bfb5235388ccf07952c4f50bd',
}
WA_COHORT_SHA256 = 'd0931d3377df1f08b3a4b1bb5393d8198bf4dcdde0bf1dd9d0e97290a0c51a0f'


@dataclass
class SourceFrequency:
    frequencies: dict
    unavailable: set
    name_audit: dict
    metadata: dict


def checked_sources(data):
    digests = {}
    for name, expected in SOURCE_SHA256.items():
        with (data / name).open('rb') as source:
            actual = hashlib.file_digest(source, 'sha256').hexdigest()
        if actual != expected:
            raise ValueError(f'Original source checksum mismatch: {name}')
        digests[name] = actual
    return digests


def audit_names(labels, normalize):
    return {label: normalize(label) for label in sorted(set(labels))}


def unresolved_keys(audit, labels, normalize):
    # A rejected aggregate can obscure detection of its namesake species.
    # Genus-only and other unresolvable labels imply no specific species key.
    bases = {' '.join(label.split()[:2]) for label in labels if audit[label] is None}
    return {normalize(base) for base in bases} - {None}


def make_frequency(frequencies, unavailable, audit, metadata):
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in frequencies.values()):
        raise ValueError('Source frequencies must be finite fractions')
    frequencies = {key: value for key, value in frequencies.items()
                   if key not in unavailable and value > 0}
    return SourceFrequency(frequencies, unavailable, audit, {
        **metadata,
        'raw_labels': len(audit),
        'raw_labels_with_key': sum(key is not None for key in audit.values()),
        'raw_labels_without_key': sum(key is None for key in audit.values()),
        'unavailable_keys': len(unavailable),
        'frequency_keys': len(frequencies),
    })


def site_frequencies(groups, normalize, metadata):
    """Each group has equal weight. Each site within a group has equal weight."""
    if not groups or any(not group for group in groups):
        raise ValueError('Source frequency needs at least one site per group')
    labels = set().union(*(site for group in groups for site in group))
    audit = audit_names(labels, normalize)
    unavailable = unresolved_keys(audit, labels, normalize)
    group_frequencies = []
    for group in groups:
        counts = Counter(key for site in group
                         for key in {audit[label] for label in site} - {None})
        group_frequencies.append({key: count / len(group) for key, count in counts.items()})
    frequencies = {
        key: sum(group.get(key, 0) for group in group_frequencies) / len(groups)
        for key in set().union(*group_frequencies)
    }
    return make_frequency(frequencies, unavailable, audit, {
        **metadata,
        'sites': sum(map(len, groups)),
        'raw_site_label_incidences': sum(len(site) for group in groups for site in group),
    })


def marginal_frequencies(rows, normalize):
    """Preserve identifiable single marginals and mask unresolved or duplicate keys."""
    audit = audit_names([label for label, _ in rows], normalize)
    unavailable = unresolved_keys(audit, [label for label, v in rows if v > 0], normalize)
    values = defaultdict(list)
    for label, value in rows:
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('Invalid original marginal frequency')
        if audit[label] is not None:
            values[audit[label]].append(value)
    collisions = {key for key, entries in values.items() if len(entries) > 1}
    return make_frequency({key: entries[0] for key, entries in values.items()},
                          unavailable | collisions, audit, {
                              'original_marginal_rows': len(rows),
                              'duplicate_key_marginals': len(collisions),
                          })


def slovenian_conflicts(rows, matched_zones, codebook, normalize):
    """Compare occurrence and codebook labels only in the selected BS/BT zones."""
    codebook_keys = {code: normalize(label) for code, label in codebook.items()}
    conflicting, unmatched, unavailable = set(), set(), set()
    field_labels, affected_site_codes = set(), set()
    included_rows = affected_rows = 0
    for row in rows:
        location = row['location_id']
        if (location, row['zone_id']) not in matched_zones:
            continue
        included_rows += 1
        label = row['latin_name']
        if isinstance(label, str) and label.strip():
            field_labels.add(label)
        code = row['EPPO_code']
        field_key, codebook_key = normalize(label), codebook_keys.get(code)
        if field_key == codebook_key:
            continue
        if code in codebook_keys:
            conflicting.add(code)
        else:
            unmatched.add(code)
        unavailable.update({field_key, codebook_key} - {None})
        affected_site_codes.add((location, code))
        affected_rows += 1
    return unavailable, field_labels, {
        'field_codebook_scope': 'occurrence labels only in matched same-side BS and BT zones',
        'field_codebook_matched_zones': len(matched_zones),
        'field_codebook_occurrence_rows': included_rows,
        'field_codebook_conflicting_codes': len(conflicting),
        'field_codebook_unmatched_strict_codes': len(unmatched),
        'field_codebook_conflict_occurrence_rows': affected_rows,
        'field_codebook_conflict_site_code_incidences': len(affected_site_codes),
        'field_codebook_conflict_sites': len({location for location, _ in affected_site_codes}),
        'field_codebook_unavailable_keys': len(unavailable),
    }


def read_europe_rail(data, normalize):
    import pymupdf

    import validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 as eu
    from analyze_Slovenia_rail_microhabitats_v103 import (
        extract_original_xlsx, pair_source, validate_original_keys,
    )
    from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
    from select_two_original_rail_20_species_for_inat_v164 import (
        original_slovenian_EPPO_codebook_binomial_mapping,
    )

    studies = {}
    for key, count, source in (
        ('CZ', 39, eu.publisher_czech(data / 'cz1.xlsx', data / 'cz2.xlsx')),
        ('IT', 28, eu.publisher_italy(data / 'parma.pdf')),
        ('HRBA', 15, eu.publisher_hr_ba(data / 'hrba.pdf')),
    ):
        sites, audit = source
        if len(sites) != count:
            raise ValueError(f'Unexpected original station count: {key}')
        studies[key] = site_frequencies([[site['source_names'] for site in sites]], normalize, {
            'grain': 'railway station',
            'original_table_taxa': audit['original_n_taxa'],
        })

    workbook = extract_original_xlsx(data / SI_XLSX, data / SI_JSON)
    sites, _, _ = source_sites_and_guilds(workbook)
    legacy, codebook_audit = original_slovenian_EPPO_codebook_binomial_mapping(workbook)
    latin = {str(row['EPPO_code']): row['latin_name'] for row in workbook['species_codebook']}
    observed_codes = set().union(*(site['all_original_taxa'] for site in sites))
    code_incidences = sum(len(site['all_original_taxa']) for site in sites)
    if (len(sites), len(observed_codes), code_incidences) != (89, 285, 2814):
        raise ValueError('Unexpected original Slovenian sites or code incidences')
    if any(eu.binomial(latin.get(code)) != legacy.get(code) for code in observed_codes):
        raise ValueError('Original Slovenian raw name no longer agrees with source codebook')
    raw_sites = [{latin[code] for code in site['all_original_taxa']
                  if isinstance(latin.get(code), str) and latin[code].strip()} for site in sites]
    geos, zones, species, traits, _, _ = validate_original_keys(workbook)
    pairs = pair_source(geos, zones, species, traits)
    matched_zones = {(pair['location_id'], pair['track_side'] + zone)
                     for pair in pairs for zone in ('BS', 'BT')}
    matched_site_codes = defaultdict(set)
    for location, zone in matched_zones:
        matched_site_codes[str(location)].update(species[location, zone])
    if (len(matched_zones) != 346 or matched_site_codes !=
            {site['site']: site['all_original_taxa'] for site in sites}):
        raise ValueError('Original Slovenian matched zone and site cohorts disagree')
    conflicts, field_labels, conflict_audit = slovenian_conflicts(
        workbook['species'], matched_zones, latin, normalize)
    slovenia = site_frequencies([raw_sites], normalize, {
        'grain': 'matched ballast shoulder/toe site union',
        'original_site_code_incidences': code_incidences,
        'original_observed_codes': len(observed_codes),
        'observed_codes_without_codebook_row': len(observed_codes - latin.keys()),
        'publisher_codebook_rows': codebook_audit['source_codebook_total_original_EPPO'],
    })
    studies['SI'] = make_frequency(slovenia.frequencies, slovenia.unavailable | conflicts,
                                   slovenia.name_audit | audit_names(field_labels, normalize), {
                                       **slovenia.metadata, **conflict_audit,
                                   })

    with pymupdf.open(data / 'sk30_2.pdf') as doc:
        if len(doc) != 9:
            raise ValueError('Unexpected Slovak Table S1 page count')
        table = '\n'.join(page.get_text() for page in doc)
    rows = [(' '.join(label.split()), int(percent)) for label, percent in
            re.findall(r'(?m)^([^\n]+)\n(\d+)[ \t]*%[ \t]*$', table)]
    if len(rows) != 413 or len({label for label, _ in rows}) != 413:
        raise ValueError('Unexpected Slovak Table S1 taxon row count')
    counts = {(100 * occupied + 15) // 30: occupied for occupied in range(31)}
    if len(counts) != 31 or any(not counts.get(pct) for _, pct in rows):
        raise ValueError('Printed Slovak percentages do not identify counts out of 30')
    studies['SK'] = marginal_frequencies([(label, counts[pct] / 30) for label, pct in rows], normalize)
    studies['SK'].metadata.update(sites=30, grain='railway station marginals',
                                  inferred_occupied_station_incidence=sum(counts[pct] for _, pct in rows))
    return studies


def read_kashmir(path, normalize):
    import pandas as pd

    from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix

    plots, workbook_audit = parse_xlsx_matrix(path)
    frame = pd.read_excel(path)
    headers = frame.columns[6:]
    if frame.shape != (248, 129):
        raise ValueError('Unexpected original Kashmir workbook shape')
    zones = ('Track', 'Embankment', 'Interior Plot 1', 'Interior Plot 2')
    stations = sorted({site for site, year, _ in plots if year == 2014 and not site.startswith('Between ')})
    if len(stations) != 16:
        raise ValueError('Expected the same 16 Kashmir stations for both 2014 habitat sets')
    selected = frame[frame.iloc[:, 4].eq(2014) & ~frame.iloc[:, 0].str.startswith('Between')]
    if len(selected) != 64:
        raise ValueError('Expected four original plots each at 16 Kashmir stations')
    for row in selected.itertuples(index=False, name=None):
        key = (' '.join(str(row[0]).split()), int(row[4]), ' '.join(str(row[5]).split()))
        values = row[6:]
        if (any(pd.isna(value) for value in values) or
                plots[key] != {i for i, value in enumerate(values) if value > 0}):
            raise ValueError('Kashmir XML and worksheet plot readings disagree')
    rail, interior = [], []
    for station in stations:
        observed = [plots.get((station, 2014, zone)) for zone in zones]
        if any(record is None for record in observed):
            raise ValueError('Incomplete habitat pair at a 2014 Kashmir station')
        rail.append({headers[i] for i in observed[0] | observed[1]})
        interior.append({headers[i] for i in observed[2] | observed[3]})
    meta = {'year': 2014, 'original_species_columns': len(headers),
            'source_missing_species_cells_all_years': workbook_audit['source_missing_species_cells']}
    return (site_frequencies([rail], normalize, {**meta, 'grain': 'station Track + Embankment union'}),
            site_frequencies([interior], normalize, {**meta, 'grain': 'same-station Interior 1 + 2 union'}))


def read_miren(path, normalize):
    import pandas as pd

    columns = ['Region', 'Road', 'Site', 'Plot', 'Year', 'Protocol', 'Accepted.Name.MIREN']
    records = pd.read_csv(path, usecols=columns, dtype={'Plot': str}, low_memory=False)
    if len(records) != 167804:
        raise ValueError('Unexpected MIREN v3 source row count')
    eligible = records[records.Plot.eq('1') & records.Protocol.isin(['v.2017', 'v.2022'])]
    earliest = eligible.groupby(['Region', 'Road', 'Site']).Year.transform('min')
    eligible = eligible[eligible.Year.eq(earliest)]
    vectors = {}
    for country, region, roads, sites in (('CZ', 'CZE', 1, 21), ('IT', 'ITA', 3, 60), ('IN', 'IND', 3, 28)):
        sub = eligible[eligible.Region.eq(region)]
        groups = [[set(site['Accepted.Name.MIREN'].dropna()) for _, site in road.groupby('Site')]
                  for _, road in sub.groupby(['Region', 'Road'])]
        if len(groups) != roads or sum(map(len, groups)) != sites:
            raise ValueError(f'Unexpected qualifying MIREN road/site counts: {country}')
        vectors[country] = site_frequencies(groups, normalize, {
            'grain': 'roadside plot 1, first eligible year per site',
            'region': region, 'roads': roads, 'years': sorted(map(int, sub.Year.unique())),
            'weighting': 'equal sites within roads, then equal roads',
        })
    return vectors


def read_tree_plots(path, normalize):
    import pandas as pd

    tree = pd.read_excel(path, header=None).iloc[4:]
    if tree.shape != (337, 9) or tree.iloc[:, [0, 1, 2, 3, 7]].isna().any().any():
        raise ValueError('Unexpected original urban tree-plot table')
    labels = list(tree.iloc[:, 0])
    vectors = {}
    for country, columns in {'SI': [1], 'HR': [2, 7], 'BA': [3]}.items():
        percentages = tree[columns].mean(axis=1)
        vectors[country] = marginal_frequencies(list(zip(labels, percentages / 100)), normalize)
        vectors[country].metadata.update(plots=50 * len(columns), cities=len(columns),
                                         grain='urban tree-pit printed frequencies')
    hr, ba = vectors['HR'], vectors['BA']
    frequencies = {
        key: (13 * hr.frequencies.get(key, 0) + 2 * ba.frequencies.get(key, 0)) / 15
        for key in hr.frequencies.keys() | ba.frequencies.keys()
    }
    vectors['HRBA'] = make_frequency(frequencies, hr.unavailable | ba.unavailable, hr.name_audit, {
        'grain': 'urban tree-pit printed frequencies', 'plots': 150, 'cities': 3,
        'country_weights': {'HR': 13, 'BA': 2},
        'original_marginal_rows': len(tree),
        'duplicate_key_marginals': hr.metadata['duplicate_key_marginals'],
    })
    return vectors['SI'], vectors['HRBA']


def read_malacky(path, normalize):
    import pymupdf

    with pymupdf.open(path) as doc:
        passage = doc[35].get_text().split('Tab. 1. Relev')[1]
    lines = [line.strip() for line in passage.splitlines() if line.strip()]
    plots = [set() for _ in range(18)]
    parsed = 0
    for index, label in enumerate(lines[:-19]):
        cells = lines[index + 1:index + 20]
        observations = cells[:17] + cells[18:]
        if not (re.fullmatch(r'\d+', cells[17]) and
                all(re.fullmatch(r'[.+rab12345]', value) for value in observations)):
            continue
        if abs(100 * sum(value != '.' for value in cells[:17]) / 17 - int(cells[17])) > 0.51:
            raise ValueError('Malacky Table 1 entries disagree with printed frequency')
        parsed += 1
        for plot, value in zip(plots, observations):
            if value != '.':
                plot.add(label)
    singleton_text = passage.split('Taxa found only in one relev')[1].split('E1:')[1]
    singletons = re.findall(r'([A-Z][a-z]+\s+[a-z]+)\s*\((\d+)\)', ' '.join(singleton_text.split()))
    if parsed != 24 or len(singletons) != 8 or len(re.findall(r'\(\d+\)', singleton_text)) != 8:
        raise ValueError('Unexpected published Malacky Table 1 taxon rows')
    for label, plot in singletons:
        plots[int(plot) - 1].add(label)
    return site_frequencies([plots[:14] + plots[15:]], normalize, {
        'grain': 'Malacky Table 1 urban relevés, excluding the railway relevé',
        'original_releves': 18, 'original_table_taxa': parsed + len(singletons),
    })


def read_wa(path, normalize):
    from extract_australian_rail_reserve_v303 import extract as extract_wa

    source = extract_wa(path)
    matrix_sites = {row['site'] for row in source['matrix']}
    if len(matrix_sites) != 137 or len(source['matrix']) != 3819:
        raise ValueError('Western Australia original Appendix I extraction mismatch')
    sheet_sites = defaultdict(set)
    for row in source['sheet_records']:
        if not row['out_of_plot']:
            sheet_sites[row['site']].add(row['taxon'])
    eligible_sites = {row['site'] for row in source['metadata'] if row['area_m2'] == 2500}
    cohort = sorted(eligible_sites & matrix_sites & sheet_sites.keys())
    cohort_digest = hashlib.sha256(('\n'.join(cohort) + '\n').encode()).hexdigest()
    sites = [sheet_sites[site] for site in cohort]
    strata = dict(Counter(site[:2] for site in cohort))
    if (cohort_digest != WA_COHORT_SHA256 or strata != {'YA': 63, 'YB': 57} or
            len(set().union(*sites)) != 320 or sum(map(len, sites)) != 3412):
        raise ValueError('Western Australia original Appendix G cohort mismatch')
    return site_frequencies([sites], normalize, {
        'grain': '2500 m² railway-reserve quadrats, Appendix G, OP excluded',
        'year': 2008, 'original_area_m2': 2500, 'cohort_sha256': cohort_digest, 'strata': strata,
    })


def load(data, normalize):
    """Provide frequency mappings, aggregate QA and in-memory name audit/masks.

    Rail and road maps use short study keys. Background maps use target and
    habitat keys. Unavailability uses short rail keys plus qualified names
    such as ``road:CZ`` and ``background:SK:urban``. False zeros are excluded
    from all maps whenever an ambiguous original label identifies a base key.
    """
    if not callable(normalize):
        raise TypeError('A source-name normalizer is required')
    data = Path(data)
    digests = checked_sources(data)
    rail = read_europe_rail(data, normalize)
    rail['IN'], interior = read_kashmir(data / 'Matrix_plot.xlsx', normalize)
    rail['AU_WA'] = read_wa(data / 'wa_rail_reserve.pdf', normalize)
    road = read_miren(data / 'miren.csv', normalize)
    si, hrba = read_tree_plots(data / 'tree1.xlsx', normalize)
    background = {
        'SI': {'urban': si},
        'HRBA': {'urban': hrba},
        'SK': {'urban': read_malacky(data / 'sk_urban.bin', normalize)},
        'IN': {'interior': interior, 'miren_road': road['IN']},
    }
    flat_road = {f'road:{country}': value for country, value in road.items()}
    flat_background = {f'background:{country}:{label}': value
                       for country, choices in background.items() for label, value in choices.items()}
    flat_rail = {f'rail:{country}': value for country, value in rail.items()}
    datasets = {**flat_rail, **flat_road, **flat_background}
    return {
        'rail': {key: value.frequencies for key, value in rail.items()},
        'road': {key: value.frequencies for key, value in road.items()},
        'background': {key: {label: value.frequencies for label, value in choices.items()}
                       for key, choices in background.items()},
        'metadata': {
            'source_sha256': digests,
            'rail': {key: value.metadata for key, value in rail.items()},
            'road': {key: value.metadata for key, value in road.items()},
            'background': {key: {label: value.metadata for label, value in choices.items()}
                           for key, choices in background.items()},
        },
        'raw_name_audit': {key: value.name_audit for key, value in datasets.items()},
        'unavailable': {key: value.unavailable
                        for key, value in {**rail, **flat_road, **flat_background}.items()},
    }

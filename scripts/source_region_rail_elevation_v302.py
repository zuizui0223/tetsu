"""Source construction for region-heldout and independent Kashmir rail transfer."""
from collections import defaultdict
from pathlib import Path
import hashlib
import math
import numpy as np
from openpyxl import load_workbook
import evaluate_MIREN_native_elevation_heldout_stature_predictability_v300 as src
from analyze_MIREN_original_road_out_transfer_v301 import road_out_calibrated_transfer
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix, DRYAD_SHA256
from validate_Kashmir_out_of_Europe_rail_nonrail_height_beta_v288 import species_headers
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial
from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import (
    FROZEN_3_STUDY_LOG_MEANS, FROZEN_3_STUDY_LOG_SDS)


def read_sources(folder):
    folder = Path(folder)
    for filename, key in [('canopy_height.txt', 'height'), ('seed_mass.txt', 'seed')]:
        assert hashlib.sha256((folder/filename).read_bytes()).hexdigest() == src.LEDASHA[key]
    flora, flora_qa = src.read_source_miren(folder/'miren.csv')
    height, _ = src.read_leda(folder/'canopy_height.txt', 'height_m')
    seed, _ = src.read_leda(folder/'seed_mass.txt', 'seed_mass_mg')
    local, traits, _ = src.trait_pool_original_matched_regions(flora, height, seed)
    elevations, elevation_qa = src.read_original_elevation_csv(folder/'miren.csv', flora)
    assert elevation_qa['n_source_plots_with_qualified_original_Elevation'] == 6446
    selected = {}
    for panel, (minimum, share, year, nr, nd, ns) in src.PANEL_DEFINITIONS.items():
        regions, _ = src.valid_region_panels(flora, local, minimum, share, year)
        assert (len(regions), sum(len(r) for r in regions.values()),
                sum(len(s) for r in regions.values() for _, s in r)) == (nr, nd, ns)
        selected[panel] = regions
    return dict(flora=flora, height=height, seed=seed, traits=traits,
                elevations=elevations, selected=selected, source_qa=flora_qa)


def tree_context(folder, tree):
    obj, qa = src.inspect_tree_rdata(Path(folder)/'trees'/f'GBOTB.extended.{tree}.rda', tree)
    tips = [str(t).replace('_', ' ') for t in np.asarray(obj['tip.label']).reshape(-1)]
    lookup = {t: i+1 for i, t in enumerate(tips)}
    parent, depth, age = src.original_root_depths(
        np.asarray(obj['edge']), np.asarray(obj['edge.length']), len(tips),
        int(np.asarray(obj['Nnode']).reshape(-1)[0]))
    return dict(lookup=lookup, parent=parent, depth=depth,
                cut=.75*age['original_total_root_to_tip_distance_median'], qa=qa)


def region_panels(source, tree, panel_name):
    minimum = src.PANEL_DEFINITIONS[panel_name][0]
    regions, audit = {}, {}
    for region, original in sorted(source['selected'][panel_name].items()):
        panel, qa = src.construct_original_region_tree_panel(
            original, source['traits'], tree['lookup'], minimum)
        audit[region] = {'phylogeny': qa, 'status': 'insufficient_phylogeny'}
        if panel is None:
            continue
        incidence, geos, habitats, roads = src.source_study_incidence(panel)
        clades = src.original_true_clade_groups(
            [tree['lookup'][t] for t in panel['species']], tree['parent'], tree['depth'], tree['cut'])
        bins = src.original_global_taxon_occupancy_quintiles(incidence)
        groups = src.joint_true_clade_and_occupancy_groups(clades, bins, len(panel['species']))
        exchange = src.original_source_joint_exchangeability(incidence, groups, geos, habitats, roads)
        audit[region]['exchangeability'] = exchange
        if not exchange['joint_phylogenetic_clade_AND_source_occupancy_null_identifiable']:
            audit[region]['status'] = 'insufficient_null_exchangeability'
            continue
        native, env = src.construct_native_environment_cohort(
            original, region, panel, source['elevations'], minimum, set(tree['lookup']))
        audit[region]['elevation'] = env
        if not env['original_environmental_altitude_support']:
            audit[region]['status'] = 'insufficient_elevation'
            continue
        retained, cv = src.original_environment_roads_with_cv_support(native, env)
        audit[region]['cv'] = cv
        if not cv['original_region_out_of_sample_heldout_elevation_support']:
            audit[region]['status'] = 'insufficient_prediction_support'
            continue
        regions[region] = dict(roads=retained, traits=panel['source_trait_vectors'], groups=groups)
        audit[region]['status'] = 'eligible'
    return regions, audit


def verify_v301(regions, reference, tree, panel):
    strict = panel.startswith('strict')
    expected = reference['strict_panels' if strict else 'primary_panels'][tree]
    gains = [road_out_calibrated_transfer(r['roads'], r['traits'], 'lowest')
             ['roadside']['gain_zlog_height_RMSE'] for r in regions.values()]
    wanted = expected['lowest_real_gain'] if strict else expected['lowest']['actual']
    observed = float(np.mean(gains))
    if abs(observed-wanted) > 1e-10:
        raise ValueError(f'v301 observed source reproduction failed {tree}/{panel}: {observed} vs {wanted}')
    return observed


def rail_records(folder, source, tree, minimum, share, height_only=False):
    """Return original complete plot observations, retaining missing-cell exclusions."""
    path = Path(folder)/'Matrix_plot.xlsx'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == DRYAD_SHA256
    plots, qa = parse_xlsx_matrix(path)
    raw_names = species_headers(path)
    names = [binomial(t) for t in raw_names]
    workbook = load_workbook(path, read_only=True, data_only=True)
    meta = {}
    for row in list(workbook.active.values)[1:]:
        site = ' '.join(str(row[0]).split())
        altitude = float(row[3])
        if not math.isfinite(altitude):
            raise ValueError('Missing original railway altitude')
        if site in meta and altitude != meta[site]:
            raise ValueError('Inconsistent original railway altitude across zones or years')
        meta[site] = altitude
    workbook.close()
    assert len(meta) == 31
    height = source['height']
    trait_names = set(height) if height_only else set(height) & set(source['seed']) & set(tree['lookup'])
    valid_names = {t for t in trait_names if height[t] > 0 and math.isfinite(height[t])}
    zheight = {t: (math.log(height[t])-FROZEN_3_STUDY_LOG_MEANS[0])/FROZEN_3_STUDY_LOG_SDS[0]
               for t in valid_names}
    cohorts = defaultdict(list)
    audits = defaultdict(lambda: {'available': 0, 'missing': 0, 'low_coverage': 0,
                                 'known_species': [], 'coverage': []})
    zones = ['Track', 'Embankment', 'Interior Plot 1', 'Interior Plot 2', 'Track+Embankment']
    for site, altitude in sorted(meta.items()):
        kind = 'between' if site.startswith('Between') else 'station'
        for year in (2014, 2017):
            for zone in zones:
                key = f'{kind}/{year}/{zone}'
                audit = audits[key]
                audit['available'] += 1
                original = [plots[site, year, z] for z in zone.split('+')]
                if any(v is None for v in original):
                    audit['missing'] += 1
                    continue
                taxa = {names[i] for v in original for i in v} - {None}
                known = taxa & valid_names
                coverage = len(known)/len(taxa) if taxa else 0
                audit['known_species'].append(len(known))
                audit['coverage'].append(coverage)
                if len(known) < minimum or coverage < share:
                    audit['low_coverage'] += 1
                    continue
                cohorts[key].append(dict(site=site, x=altitude/1000,
                                         y=float(np.mean([zheight[t] for t in sorted(known)]))))
    # Only aggregate metadata and estimates are exported by the analysis entry point.
    return cohorts, dict(audits), qa

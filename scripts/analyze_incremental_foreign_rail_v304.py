"""Source-replayed study holdout and Kashmir spatial-block incremental predictions.

Only aggregate results are written. Source species/site rows remain in memory.
See docs/tetsu_v304_incremental_foreign_rail_prediction_plan_ja.md.
"""
from pathlib import Path
import argparse
from collections import Counter
import hashlib
import json
import re
import runpy
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
PENALTIES = (.01, .001, .1)


def logit(x):
    x = np.clip(np.asarray(x, float), .01, .99)
    return np.log(x / (1 - x))


def metrics(y, p):
    y, p = np.asarray(y, float), np.clip(np.asarray(p, float), 1e-12, 1 - 1e-12)
    return {'log_loss': float(np.mean(-y * np.log(p) - (1-y) * np.log1p(-p))),
            'brier': float(np.mean(y * (1-p)**2 + (1-y) * p**2))}


def fit_predict(x, y, weights, test_x, penalty):
    x, test_x, y = np.asarray(x, float), np.asarray(test_x, float), np.asarray(y, float)
    assert np.isfinite(x).all() and np.isfinite(test_x).all() and np.isfinite(y).all()
    assert np.all((y >= 0) & (y <= 1)) and x.shape[1] == test_x.shape[1]
    weights = np.asarray(weights, float)
    weights = weights / weights.sum()
    mean = np.sum(x * weights[:, None], axis=0)
    sd = np.sqrt(np.sum((x-mean)**2 * weights[:, None], axis=0))
    sd[sd < 1e-10] = 1.
    design = np.column_stack([np.ones(len(x)), (x-mean)/sd])
    new = np.column_stack([np.ones(len(test_x)), (test_x-mean)/sd])

    def objective(b):
        eta = design @ b
        val = weights @ (np.logaddexp(0, eta) - y * eta) + penalty/2 * np.sum(b[1:]**2)
        grad = design.T @ (weights * (expit(eta)-y))
        grad[1:] += penalty*b[1:]
        return val, grad

    start = np.zeros(design.shape[1])
    start[0] = float(logit(weights @ y))
    opt = minimize(objective, start, jac=True, method='L-BFGS-B',
                   options={'maxiter': 1000, 'ftol': 1e-12, 'gtol': 1e-8})
    if not opt.success:
        raise RuntimeError(opt.message)
    return expit(new @ opt.x), opt.x.tolist()


def mean_frequencies(groups):
    taxa = set().union(*(set(v) for v in groups))
    return {t: float(np.mean([v.get(t, 0.) for v in groups])) for t in taxa}


def national_design(rail, refs, held):
    train_ids = [s for s in rail if s != held]
    samples = {}
    for s in [*train_ids, held]:
        donors = [rail[d] for d in train_ids if d != s]
        support = Counter(t for donor in donors for t in donor)
        taxa = sorted(t for t, n in support.items() if n >= 2 and refs[s].get(t, 0.) > 0)
        pred = mean_frequencies(donors)
        x = np.array([[float(logit(refs[s][t])), float(logit(pred[t]))] for t in taxa]).reshape(-1, 2)
        samples[s] = {'x': x, 'y': np.array([rail[s].get(t, 0.) for t in taxa]), 'taxa': taxa}
    eligible = [s for s in train_ids if len(samples[s]['taxa']) >= 20]
    if len(eligible) < 3 or len(samples[held]['taxa']) < 20:
        return {'not_identifiable': True, 'training_sizes': {s: len(samples[s]['taxa']) for s in train_ids},
                'target_size': len(samples[held]['taxa'])}
    return {'x': np.concatenate([samples[s]['x'] for s in eligible]),
            'y': np.concatenate([samples[s]['y'] for s in eligible]),
            'w': np.concatenate([np.full(len(samples[s]['y']), 1/len(samples[s]['y'])) for s in eligible]),
            'test_x': samples[held]['x'], 'test_y': samples[held]['y'], 'taxa': samples[held]['taxa'],
            'training': eligible, 'training_sizes': {s: len(samples[s]['taxa']) for s in eligible}}


def source_replay(data):
    sys.path.insert(0, str(ROOT/'research/validation_20261117'))
    private = data/'v304_country_source_replay.json'
    previous = sys.argv[:]
    try:
        sys.argv = ['country_validation.py', str(data), str(private)]
        source = runpy.run_path(str(ROOT/'research/validation_20261117/country_validation.py'))
    finally:
        sys.argv = previous
    fresh = json.loads(private.read_text(encoding='utf-8'))
    frozen = json.loads((ROOT/'research/validation_20261117/country_validation.json').read_text(encoding='utf-8'))
    if fresh['sha256'] != frozen['sha256']:
        raise ValueError('Source bytes differ from existing frozen country analysis')
    for key in frozen['comparisons']:
        for panel in ('common', 'extended'):
            for measure in ('foreign', 'domestic', 'delta', 'mae_foreign', 'mae_domestic'):
                a, b = fresh['comparisons'][key][panel][measure], frozen['comparisons'][key][panel][measure]
                if a is None or b is None:
                    assert a == b
                else:
                    assert abs(a-b) < 1e-9, (key, panel, measure, a, b)
    return source, {'input_sha256': fresh['sha256'], 'replayed_comparisons': len(frozen['comparisons']),
                    'numeric_checks': len(frozen['comparisons'])*2*5,
                    'maximum_allowed_error': 1e-9}


def national(source):
    original, controls = source['target'], source['control']
    rail = {k: original[k] for k in ('CZ', 'IT', 'SI', 'IN')}
    rail['AD'] = source['studies']['HR_BA15']
    urban = {t: (13*controls['HR', 'urban'].get(t, 0) + 2*controls['BA', 'urban'].get(t, 0))/15
             for t in set(controls['HR', 'urban']) | set(controls['BA', 'urban'])}
    answer = []
    for cz_kind in ('urban', 'road'):
        refs = {'CZ': controls['CZ', cz_kind], 'IT': controls['IT', 'road'],
                'SI': controls['SI', 'urban'], 'IN': controls['IN', 'road'], 'AD': urban}
        for penalty in PENALTIES:
            for held in [*rail, 'AU', 'CN']:
                target_rail, target_refs = rail, refs
                if held in ('AU', 'CN'):
                    # Both cities of the same publication are absent from all training.
                    target_rail = dict(rail, **{held: original[held]})
                    target_refs = dict(refs, **{held: controls[held, 'road']})
                d = national_design(target_rail, target_refs, held)
                row = {'held_rail_study': held, 'cz_reference': cz_kind, 'penalty': penalty,
                       'small_sample_external_diagnostic': held in ('AU', 'CN')}
                if d.get('not_identifiable'):
                    row.update(d)
                    answer.append(row)
                    continue
                row.update({'target_taxa': len(d['taxa']), 'target_recorded_taxa': int(sum(d['test_y'] > 0)),
                            'target_incidence_coverage': float(sum(d['test_y'])/sum(target_rail[held].values())),
                            'training_studies': d['training'], 'training_sizes': d['training_sizes']})
                models = {}
                for name, cols in {'intercept': [], 'domestic': [0], 'foreign': [1], 'domestic_foreign': [0, 1]}.items():
                    p, beta = fit_predict(d['x'][:, cols], d['y'], d['w'], d['test_x'][:, cols], penalty)
                    models[name] = dict(metrics(d['test_y'], p), coefficients=beta)
                row['models'] = models
                row['foreign_increment'] = {k: models['domestic'][k]-models['domestic_foreign'][k]
                                            for k in ('log_loss', 'brier')}
                answer.append(row)
    return answer


def block_frequencies(matrix, train_sites):
    selected = np.asarray(matrix)[train_sites]
    if not len(selected) or not np.isfinite(selected).all():
        raise ValueError('Incomplete/empty training reference')
    return selected.mean(axis=0)


def crossfit_training_frequencies(matrix, train_sites):
    """Leave each training site out of its own outcome-derived frequency feature."""
    matrix, train_sites = np.asarray(matrix), np.asarray(train_sites, bool)
    n = int(train_sites.sum())
    if n < 2:
        raise ValueError('At least two training sites required')
    total = matrix[train_sites].sum(axis=0)
    return np.array([(total-matrix[j])/(n-1) if train_sites[j] else total/n
                     for j in range(len(matrix))])


def latitude(raw):
    values = re.findall(r'\d+(?:\.\d+)?', str(raw))
    if len(values) != 3:
        raise ValueError('Unrecognized original latitude')
    d, m, s = map(float, values)
    if not (0 <= m < 60 and 0 <= s < 60):
        raise ValueError('Invalid original latitude')
    return d + m/60 + s/3600


def local_source(data, source):
    import pandas as pd
    from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
    original, audit = parse_xlsx_matrix(data/'Matrix_plot.xlsx')
    frame = pd.read_excel(data/'Matrix_plot.xlsx')
    species = [source['eu'].binomial(t) for t in frame.columns[6:]]
    known = [t for t in species if t]
    if len(set(known)) != len(known):
        raise ValueError('Duplicate binomial mapping: cannot manufacture disjoint species')
    plots = {}
    meta = {}
    for _, r in frame.iterrows():
        site = ' '.join(str(r.iloc[0]).split())
        key = site, int(r.iloc[4]), str(r.iloc[5])
        values = r.iloc[6:].to_numpy(dtype=float)
        independently_parsed = None if np.isnan(values).any() else set(np.flatnonzero(values > 0))
        assert independently_parsed == original[key], key
        plots[key] = values
        meta[site] = {'lat': latitude(r.iloc[1]), 'altitude': float(r.iloc[3]),
                      'station': not site.startswith('Between')}
    return species, plots, meta, audit


def local(data, source):
    species, plots, meta, audit = local_source(data, source)
    foreign_studies = [source['studies'][s] for s in ('CZ39', 'IT28', 'HR_BA15', 'SI89')]
    support = Counter(t for f in foreign_studies for t in f)
    foreign = mean_frequencies(foreign_studies)
    rail2 = mean_frequencies([source['studies'][s] for s in ('CZ39', 'IT28')])
    road2 = mean_frequencies([source['control'][c, 'road'] for c in ('CZ', 'IT')])
    results = []
    for year in (2014, 2017):
        zones = ('Track', 'Interior Plot 1', 'Interior Plot 2')
        sites = sorted((s for s in meta if all((s, year, z) in plots and
                                              np.isfinite(plots[s, year, z]).all() for z in zones)),
                       key=lambda s: (meta[s]['lat'], s))
        all_y = np.array([plots[s, year, 'Track'] > 0 for s in sites], float)
        blocks = np.empty(len(sites), int)
        for b, idx in enumerate(np.array_split(np.arange(len(sites)), 5)):
            blocks[idx] = b
        for panel in ('foreign_two_study_support', 'all_original_binomials'):
            indices = [i for i, t in enumerate(species) if t and
                       (panel == 'all_original_binomials' or support[t] >= 2)]
            taxa = [species[i] for i in indices]
            if len(taxa) < 10 or len(sites) < 15:
                results.append({'year': year, 'panel': panel, 'not_identifiable': True,
                                'sites': len(sites), 'taxa': len(taxa)})
                continue
            y = all_y[:, indices]
            i1 = np.array([plots[s, year, zones[1]][indices] > 0 for s in sites], float)
            i2 = np.array([plots[s, year, zones[2]][indices] > 0 for s in sites], float)
            f = logit([foreign.get(t, 0) for t in taxa])
            fr = logit([rail2.get(t, 0) for t in taxa])
            fd = logit([road2.get(t, 0) for t in taxa])
            for penalty in PENALTIES:
                predictions = {n: np.full(y.shape, np.nan) for n in
                               ('local', 'local_foreign', 'local_track', 'local_track_foreign',
                                'local_foreign_rail2', 'local_foreign_road2')}
                foreign_beta = []
                for b in range(5):
                    train = blocks != b
                    local_frequency = logit(block_frequencies((i1+i2)/2, train))
                    target_frequency = logit(crossfit_training_frequencies(y, train))
                    x = np.array([np.column_stack([
                        local_frequency, i1[j], i2[j],
                        np.full(len(taxa), meta[s]['altitude']),
                        np.full(len(taxa), meta[s]['station']),
                        f, target_frequency[j], fr, fd]) for j, s in enumerate(sites)])
                    if len(np.unique(y[train])) != 2:
                        raise ValueError('Training fold lacks both outcomes')
                    for name, cols in {
                        'local': [0, 1, 2, 3, 4], 'local_foreign': [0, 1, 2, 3, 4, 5],
                        'local_track': [0, 1, 2, 3, 4, 6], 'local_track_foreign': [0, 1, 2, 3, 4, 6, 5],
                        'local_foreign_rail2': [0, 1, 2, 3, 4, 7],
                        'local_foreign_road2': [0, 1, 2, 3, 4, 8],
                    }.items():
                        train_x = x[train][:, :, cols].reshape(-1, len(cols))
                        test_x = x[~train][:, :, cols].reshape(-1, len(cols))
                        p, beta = fit_predict(train_x, y[train].ravel(), np.ones(train.sum()*len(taxa)), test_x, penalty)
                        predictions[name][~train] = p.reshape(-1, len(taxa))
                        if name == 'local_foreign':
                            foreign_beta.append(beta[-1])
                model_metrics = {n: metrics(y, p) for n, p in predictions.items()}
                diagnostic = {}
                if panel == 'all_original_binomials':
                    # Post-primary diagnostic, no refitting or new selection of primary models.
                    for label, mask in (
                        ('at_least_two_European_rail_studies', np.array([support[t] >= 2 for t in taxa])),
                        ('fewer_than_two_European_rail_studies', np.array([support[t] < 2 for t in taxa]))):
                        diagnostic[label] = {'taxa': int(mask.sum()), 'recorded_incidences': int(y[:, mask].sum()),
                                             'models': {n: metrics(y[:, mask], p[:, mask]) for n, p in predictions.items()}}
                comparisons = {}
                for name, baseline, extended in (
                    ('foreign_beyond_local', 'local', 'local_foreign'),
                    ('foreign_beyond_local_and_track', 'local_track', 'local_track_foreign'),
                    ('matched_foreign_rail_vs_road', 'local_foreign_road2', 'local_foreign_rail2')):
                    differences = {k: [] for k in ('log_loss', 'brier')}
                    for j in range(len(sites)):
                        left, right = metrics(y[j], predictions[baseline][j]), metrics(y[j], predictions[extended][j])
                        for k in differences:
                            differences[k].append(left[k]-right[k])
                    comparison = {}
                    for k, values in differences.items():
                        values = np.asarray(values)
                        rng = np.random.default_rng(3042026)
                        draws = []
                        for _ in range(999):
                            sampled = rng.integers(0, 5, 5)
                            draws.append(float(np.concatenate([values[blocks == b] for b in sampled]).mean()))
                        comparison[k] = {
                            'gain': float(values.mean()),
                            'block_gains': [float(values[blocks == b].mean()) for b in range(5)],
                            'leave_block_mean_range': [float(min(values[blocks != b].mean() for b in range(5))),
                                                      float(max(values[blocks != b].mean() for b in range(5)))],
                            'block_bootstrap_95_percentile': np.quantile(draws, [.025, .975]).tolist(),
                            'station_gain': float(values[[meta[s]['station'] for s in sites]].mean()),
                            'between_gain': float(values[[not meta[s]['station'] for s in sites]].mean())}
                    comparisons[name] = comparison
                results.append({'year': year, 'panel': panel, 'penalty': penalty,
                                'sites': len(sites), 'station_sites': sum(meta[s]['station'] for s in sites),
                                'taxa': len(taxa), 'block_sizes': [int(sum(blocks == b)) for b in range(5)],
                                'target_recorded_incidences': int(y.sum()),
                                'target_all_taxon_incidences': int(all_y.sum()),
                                'incidence_coverage': float(y.sum()/all_y.sum()),
                                'foreign_beta_range': [min(foreign_beta), max(foreign_beta)],
                                'models': model_metrics, 'comparisons': comparisons,
                                'posthoc_fixed_prediction_support_diagnostic': diagnostic})
    return {'source': audit, 'panels': results}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, replay = source_replay(args.data.resolve())
    print('Original country analysis and hashes reproduced', flush=True)
    result = {'design': 'v304 exploratory; national study holdout and within-Kashmir spatial holdout are distinct',
              'replay': replay, 'national': national(source)}
    print('National study-holdout models completed', flush=True)
    result['local_kashmir'] = local(args.data.resolve(), source)
    result['rail_donor_sha256'] = {name: hashlib.sha256((args.data/name).read_bytes()).hexdigest()
                                  for name in ('cz1.xlsx', 'cz2.xlsx', 'parma.pdf', 'hrba.pdf')}
    result['slovenia_source_sha256'] = hashlib.sha256(next((args.data/'si').rglob('*.xlsx')).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('All planned panels saved', flush=True)


if __name__ == '__main__':
    main()

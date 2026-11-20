"""v30.2 original-data region-heldout and non-European railway prediction.
Run: python scripts/analyze_region_rail_elevation_transfer_v302.py --data DATA --output RESULT
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import source_region_rail_elevation_v302 as sources

ANCHORS = ('lowest', 'median', 'none')
ZONES = ('roadside', 'interior')
DRAWS = 398


def slope(x, y):
    dx = x-x.mean()
    if dx @ dx <= 1e-12:
        raise ValueError('No elevation variation')
    return float(dx @ (y-y.mean())/(dx @ dx))


def xy_regions(regions, zone, randomized=None):
    out = {}
    for name, region in regions.items():
        traits = region['traits'] if randomized is None else randomized[name]
        out[name] = [(np.asarray(r['elev_m'])/1000,
                      r['matrices'][zone] @ traits[:, 0]) for r in region['roads']]
    return out


def fit_foreign(training, minimum_sources=2):
    """Every region, then road, then site has equal weight within its parent."""
    if len(training) < minimum_sources:
        raise ValueError('At least two training regions required')
    slopes = [np.mean([slope(x, y) for x, y in roads]) for roads in training.values()]
    x, y, w = [], [], []
    for roads in training.values():
        for xx, yy in roads:
            x.extend(xx)
            y.extend(yy)
            w.extend([1/(len(training)*len(roads)*len(xx))]*len(xx))
    x, y, w = map(np.asarray, (x, y, w))
    mx, my = float(w @ x), float(w @ y)
    b = float((w*(x-mx)) @ (y-my)/(w @ ((x-mx)**2)))
    return dict(slope=float(np.mean(slopes)), cold_slope=b, cold_intercept=my-b*mx,
                baseline=my, training_elevation_km=[float(x.min()), float(x.max())])


def predict_score(x, y, model, anchor):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) != len(y) or len(x) < 3 or not np.isfinite(x+y).all():
        raise ValueError('Invalid evaluation arrays')
    if anchor == 'none':
        observed = y
        pred = model['cold_intercept']+model['cold_slope']*x
        base = np.full(len(y), model['baseline'])
        xx = x
    else:
        if anchor == 'lowest':
            index = int(np.argmin(x))
        elif anchor == 'median':
            index = int(np.argmin(abs(x-np.median(x))))
        else:
            raise ValueError('Unknown anchor rule')
        take = np.arange(len(x)) != index
        xx, observed = x[take], y[take]
        pred = y[index]+model['slope']*(xx-x[index])
        base = np.full(sum(take), y[index])
    rmse = float(np.sqrt(np.mean((observed-pred)**2)))
    baseline = float(np.sqrt(np.mean((observed-base)**2)))
    lo, hi = model['training_elevation_km']
    return dict(gain=baseline-rmse, rmse=rmse, baseline_rmse=baseline,
                evaluated_sites=len(observed), outside_training_elevation=int(sum((xx<lo)|(xx>hi))))


def average_scores(scores):
    return {k:float(np.mean([s[k] for s in scores])) for k in ('gain','rmse','baseline_rmse')}


def region_holdout(arrays):
    results = {}
    for excluded, roads in arrays.items():
        train = {r:v for r,v in arrays.items() if r != excluded}
        model = fit_foreign(train)
        scores = {}
        for anchor in ANCHORS:
            road_scores = [predict_score(x, y, model, anchor) for x, y in roads]
            scores[anchor] = dict(**average_scores(road_scores),
                                  positive_roads=sum(v['gain']>0 for v in road_scores),
                                  n_roads=len(roads), n_sites=sum(len(x) for x,y in roads),
                                  tested_sites=sum(v['evaluated_sites'] for v in road_scores),
                                  outside_training_elevation=sum(v['outside_training_elevation'] for v in road_scores))
        results[excluded] = dict(training_regions=sorted(train), model=model,
                                observed_slope=float(np.mean([slope(x,y) for x,y in roads])), scores=scores)
    return results


def region_analysis(regions, seed):
    observed = {zone:region_holdout(xy_regions(regions, zone)) for zone in ZONES}
    null = {(z, r, a):[] for z in ZONES for r in regions for a in ANCHORS}
    rng = np.random.default_rng(seed)
    for _ in range(DRAWS):
        shuffled = {name:sources.src.permute_paired_traits(r['traits'],r['groups'],rng)
                    for name,r in sorted(regions.items())}
        for zone in ZONES:
            current = region_holdout(xy_regions(regions,zone,shuffled))
            for region,row in current.items():
                for anchor in ANCHORS:
                    null[zone,region,anchor].append(row['scores'][anchor]['gain'])
    macro = {}
    for zone in ZONES:
        arrays = xy_regions(regions,zone)
        for region,row in observed[zone].items():
            for anchor in ANCHORS:
                mean_null=float(np.mean(null[zone,region,anchor]))
                row['scores'][anchor]['null_mean_gain']=mean_null
                row['scores'][anchor]['excess_gain']=row['scores'][anchor]['gain']-mean_null
                removed={}
                if len(regions)>=4:
                    for omitted in regions:
                        if omitted==region:
                            continue
                        model=fit_foreign({r:v for r,v in arrays.items() if r not in (region,omitted)})
                        removed[omitted]=average_scores([predict_score(x,y,model,anchor) for x,y in arrays[region]])['gain']
                row['scores'][anchor]['extra_training_region_removed']=removed
        macro[zone]={}
        for anchor in ANCHORS:
            values=[r['scores'][anchor] for r in observed[zone].values()]
            gains=[r['gain'] for r in values]
            macro[zone][anchor]=dict(gain=float(np.mean(gains)),excess=float(np.mean([r['excess_gain'] for r in values])),
                                     positive_regions=sum(v>0 for v in gains),n_regions=len(values),
                                     baseline_rmse=float(np.mean([r['baseline_rmse'] for r in values])),
                                     rmse=float(np.mean([r['rmse'] for r in values])),
                                     delete_one_region_mean_range=[min(float(np.mean(np.delete(gains,i))) for i in range(len(gains))),
                                                                   max(float(np.mean(np.delete(gains,i))) for i in range(len(gains)))])
    stable={name:dict(r, roads=sources.src.original_year_consistent_roads(r['roads'])) for name,r in regions.items()}
    stable={name:r for name,r in stable.items() if len(r['roads'])>=2 and
            sum(len(x['elev_m']) for x in r['roads'])>=.6*sum(len(x['elev_m']) for x in regions[name]['roads'])}
    sensitivity={z:region_holdout(xy_regions(stable,z)) for z in ZONES} if len(stable)>=3 else None
    return dict(regions=observed,macro=macro,single_year_roads=sensitivity,null_draws=DRAWS)


def rail_score(rows, training):
    x=np.array([r['x'] for r in rows]);y=np.array([r['y'] for r in rows])
    out=dict(n=len(rows),elevation_span_m=float(np.ptp(x)*1000) if len(x) else 0,
             status='insufficient_sites_or_elevation')
    if len(x)<6 or np.ptp(x)<.15-1e-12:
        return out
    model=fit_foreign(training)
    out.update(status='eligible',observed_slope=slope(x,y),model=model,
               scores={a:predict_score(x,y,model,a) for a in ANCHORS})
    delete=[slope(np.delete(x,i),np.delete(y,i)) for i in range(len(x))]
    out['delete_one_site_slope_range']=[min(delete),max(delete)]
    out['negative_delete_one_site_slopes']=sum(v<0 for v in delete)
    nohigh=np.arange(len(x))!=np.argmax(x)
    out['highest_site_removed']=dict(n=int(sum(nohigh)),elevation_span_m=float(np.ptp(x[nohigh])*1000),
                                    slope=slope(x[nohigh],y[nohigh]),
                                    scores={a:predict_score(x[nohigh],y[nohigh],model,a) for a in ANCHORS})
    out['extra_training_region_removed']={}
    if len(training)>=3:
        for omitted in training:
            alternative=fit_foreign({r:v for r,v in training.items() if r!=omitted})
            out['extra_training_region_removed'][omitted]={a:predict_score(x,y,alternative,a)['gain'] for a in ANCHORS}
    return out


def rail_analysis(cohorts,audit,training):
    results={}
    for key, qa in sorted(audit.items()):
        rows=cohorts.get(key,[])
        other=key.replace('/2014/','/2017/') if '/2014/' in key else key.replace('/2017/','/2014/')
        shared={r['site'] for r in cohorts.get(other,[])}
        known=qa['known_species'];coverage=qa['coverage']
        results[key]=dict(result=rail_score(rows,training),
                         same_sites_both_years=rail_score([r for r in rows if r['site'] in shared],training),
                         source_available=qa['available'],missing=qa['missing'],low_coverage=qa['low_coverage'],
                         known_species_min_median_max=[min(known),float(np.median(known)),max(known)] if known else None,
                         coverage_min_median_max=[min(coverage),float(np.median(coverage)),max(coverage)] if coverage else None)
    return results


def run(data,output):
    base=Path(__file__).resolve().parents[1]
    frozen=json.loads((base/'data/derived/tetsu_v301_original_MIREN_ROAD_out_calibrated_source_macro_summary.json').read_text())
    source=sources.read_sources(data)
    result=dict(version='v30.2',reference_commit='94024e126b21724f229538f2e1c0f4a21c8313e7',
                primary_rail='station/2014/Track+Embankment',panels={})
    for ti,tree_name in enumerate(sources.src.TREES):
        tree=sources.tree_context(data,tree_name)
        for pi,panel in enumerate(sources.src.PANEL_DEFINITIONS):
            print('Running',tree_name,panel,flush=True)
            regions,audit=sources.region_panels(source,tree,panel)
            observed=sources.verify_v301(regions,frozen,tree_name,panel)
            new=region_analysis(regions,302200+ti*100+pi)
            training={r:arr for r,arr in xy_regions(regions,'roadside').items() if r!='IND'}
            minimum,share=sources.src.PANEL_DEFINITIONS[panel][:2]
            cohorts,rail_audit,_=sources.rail_records(data,source,tree,minimum,share)
            rail=rail_analysis(cohorts,rail_audit,training)
            item=dict(v301_observed_gain_reproduced=observed,source_support=audit,
                      road=new,rail=rail,tree_source=tree['qa'])
            if tree_name=='LCVP':
                hc,ha,_=sources.rail_records(data,source,tree,minimum,share,height_only=True)
                item['rail_height_only_all_LEDA_sensitivity']=rail_analysis(hc,ha,training)
            result['panels'][tree_name+'/'+panel]=item
            output.parent.mkdir(parents=True,exist_ok=True)
            output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.data,args.output)

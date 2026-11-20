"""Exploratory v305; raw sources remain local, exports contain aggregates only.

See docs/tetsu_v305_traffic_urbanization_analysis_plan_ja.md (committed before join).
OLS coefficients describe associations; neither site pairs nor MCMC draws
are treated as independent biological replicates. No naive pairwise p-values.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import numpy as np
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites, jaccard_similarity
from analyze_Slovenia_historical_train_frequency_local_species_coupling_v275 import site_specific_null, mapped_segment_data
from analyze_rail_local_site_coupling_v257 import slovenia_location_sets
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from audit_Slovenia_2015_2020_real_train_local_coupling_temporal_v276 import historical_site_rematch, EUROSTAT_SHA
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import eurostat_2020_segments
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import query_series, unpack_jsonstat

BUFFERS=('r1000_inner0','r500_inner0','r1000_inner100')

def fit(rows, response, predictors, weight_group=None):
    names=['intercept']+list(predictors)
    if not rows:return {'status':'no_support','n':0}
    x=np.array([[1.]+[float(r[p]) for p in predictors] for r in rows])
    y=np.array([r[response] for r in rows],float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Missing covariates must not become zero')
    counts=Counter(r[weight_group] for r in rows) if weight_group else None
    w=np.array([1/counts[r[weight_group]] if counts else 1 for r in rows])
    xw=x*np.sqrt(w[:,None]);yw=y*np.sqrt(w)
    rank=int(np.linalg.matrix_rank(xw));df=len(rows)-rank
    out={'n':len(rows),'rank':rank,'residual_df':df,
         'condition_number':float(np.linalg.cond(xw)),
         'weighting':weight_group or 'equal_units'}
    if rank<len(names) or df<3:
        return dict(out,status='not_identifiable')
    coef=np.linalg.lstsq(xw,yw,rcond=None)[0]
    return dict(out,status='descriptive',coefficients=dict(zip(names,map(float,coef))),
                weighted_rmse=float(np.sqrt(np.average((y-x@coef)**2,weights=w))))

def model(rows,response,predictors,target,weight_group=None,delete_group='section'):
    out=fit(rows,response,predictors,weight_group)
    vals=[]
    for g in sorted({r[delete_group] for r in rows}):
        f=fit([r for r in rows if r[delete_group]!=g],response,predictors,weight_group)
        if f['status']=='descriptive':vals.append(f['coefficients'][target])
    out['leave_one_group_out']={'group':delete_group,'n_estimable':len(vals),
        'positive':sum(x>0 for x in vals),'negative':sum(x<0 for x in vals),
        'min':min(vals) if vals else None,'max':max(vals) if vals else None}
    return out

def urban_rows(sites,urban,buffer):
    ix={str(r['source_id']):r for r in urban['sites']}
    if len(ix)!=len(urban['sites']):raise ValueError('Duplicate urban source IDs')
    rows=[]
    for s in sites:
        u=ix[s['source_id']]['buffers'][buffer]
        if u['valid_area_fraction']<.95:continue
        rows.append(dict(s,urban=u['built_up_fraction'],tree=u['tree_cover_fraction'],
                         crop=u['cropland_fraction'],main=float(s['class']=='main_G'),
                         year2023=float(s['year']==2023),double=float(s['layout']=='DOUBLE')))
    return rows

def class_models(rows):
    out={}
    ranges={c:(min(r['urban'] for r in rows if r['class']==c),
               max(r['urban'] for r in rows if r['class']==c)) for c in ('main_G','regional_R')}
    lo=max(v[0] for v in ranges.values());hi=min(v[1] for v in ranges.values())
    out['urban_ranges']=ranges;out['urban_common_range']=[lo,hi]
    out['class_covariates']={c:{'n_sites':len(rs),'n_lines':len({r['section'] for r in rs}),
        'site_mean_built_up':float(np.mean([r['urban'] for r in rs])),
        'site_mean_tree':float(np.mean([r['tree'] for r in rs])),
        'site_mean_crop':float(np.mean([r['crop'] for r in rs]))}
        for c in ranges for rs in [[r for r in rows if r['class']==c]]}
    for subset,rs in [('all',rows),('single',[r for r in rows if not r['double']]),
                      ('common_urban',[r for r in rows if lo<=r['urban']<=hi])]:
        for response in ('local_excess_beta','local_similarity'):
            for adjust in (False,True):
                predictors=['main']+(['urban','year2023'] if adjust else [])
                if adjust and subset!='single':predictors+=['double']
                out[f'{subset}/{response}/adjusted={adjust}']=model(rs,response,predictors,'main','section')
    return out

def traffic_models(rows,traffic):
    ix={r['source_id']:r for r in rows};groups=defaultdict(list)
    for t in traffic:
        if str(t['sid']) in ix:groups[t['source_segment']].append((ix[str(t['sid'])],t))
    segments=[]
    for key,rs in sorted(groups.items()):
        if len({(t['freight'],t['passenger'],t['section']) for r,t in rs})!=1:
            raise ValueError('Inconsistent traffic within physical segment')
        segments.append({'segment':key,'section':rs[0][1]['section'],
            'freight':float(np.log1p(rs[0][1]['freight'])),
            'passenger':float(np.log1p(rs[0][1]['passenger'])),
            **{k:float(np.mean([r[k] for r,t in rs])) for k in
               ('urban','local_excess_beta','local_similarity')}})
    out={'n_segments':len(segments),'n_sites':sum(map(len,groups.values())),
         'n_lines':len({r['section'] for r in segments}),'models':{}}
    lines=sorted({r['section'] for r in segments})
    for s in segments:
        for line in lines[1:]:s['line_'+line]=float(s['section']==line)
    for response in ('local_excess_beta','local_similarity'):
        for exposure in ('freight','passenger'):
            for adjust in ('none','urban','urban_line'):
                predictors=[exposure]+(['urban'] if adjust!='none' else [])
                if adjust=='urban_line':predictors+=['line_'+x for x in lines[1:]]
                v=model(segments,response,predictors,exposure,delete_group='segment')
                v['leave_one_line_out']=model(segments,response,predictors,exposure,delete_group='section')['leave_one_group_out']
                out['models'][f'{response}/{exposure}/{adjust}']=v
    return out

def traffic_pair_models(rows,traffic):
    """Nine segment vertices, NOT 36 independent traffic replicates."""
    ix={r['source_id']:r for r in rows};groups=defaultdict(list)
    for t in traffic:
        if str(t['sid']) in ix:groups[t['source_segment']].append((ix[str(t['sid'])],t))
    keys=sorted(groups);pairs=[]
    for i,ka in enumerate(keys):
        for kb in keys[i+1:]:
            values=[]
            for a,ta in groups[ka]:
                for b,tb in groups[kb]:
                    values.append({'distance':haversine_km(a['lat'],a['lon'],b['lat'],b['lon']),
                        'urban_mean':(a['urban']+b['urban'])/2,
                        'urban_diff':abs(a['urban']-b['urban']),
                        **{h:jaccard_similarity(a['BS']|a['BT'] if h=='union' else a[h],
                                               b['BS']|b['BT'] if h=='union' else b[h]) for h in ('BS','BT','union')}})
            t1=groups[ka][0][1];t2=groups[kb][0][1]
            r={k:float(np.mean([v[k] for v in values])) for k in values[0]}
            r.update(left=ka,right=kb,log_distance=float(np.log1p(r['distance'])),
                     freight=(np.log1p(t1['freight'])+np.log1p(t2['freight']))/2,
                     passenger=(np.log1p(t1['passenger'])+np.log1p(t2['passenger']))/2)
            pairs.append(r)
    out={'n_segment_vertices':len(keys),'n_dependent_pairs':len(pairs),
         'warning':'Descriptive dyadic regression; residual df is algebraic, not biological; no p-values',
         'models':{}}
    for habitat in ('BS','BT','union'):
        for exposure in ('freight','passenger'):
            for adjusted in (False,True):
                predictors=[exposure,'log_distance']+(['urban_mean','urban_diff'] if adjusted else [])
                result=fit(pairs,habitat,predictors)
                slopes=[]
                for key in keys:
                    f=fit([r for r in pairs if key not in (r['left'],r['right'])],habitat,predictors)
                    if f['status']=='descriptive':slopes.append(f['coefficients'][exposure])
                result['leave_one_segment_vertex_out']={'n':len(slopes),
                    'positive':sum(v>0 for v in slopes),'negative':sum(v<0 for v in slopes),
                    'min':min(slopes) if slopes else None,'max':max(slopes) if slopes else None}
                out['models'][f'{habitat}/{exposure}/urban={adjusted}']=result
    return out

def spatial(rows):
    bins=[(0,5),(5,20),(20,50),(50,float('inf'))]
    groups=defaultdict(lambda:defaultdict(list))
    for i,a in enumerate(rows):
        for b in rows[i+1:]:
            d=haversine_km(a['lat'],a['lon'],b['lat'],b['lon'])
            band=next(f'{lo}-{hi}' for lo,hi in bins if lo<=d<hi)
            cl='/'.join(sorted([a['class'],b['class']]))
            linepair=tuple(sorted([a['section'],b['section']]))
            for habitat in ('BS','BT','union'):
                aa=a['BS']|a['BT'] if habitat=='union' else a[habitat]
                bb=b['BS']|b['BT'] if habitat=='union' else b[habitat]
                groups[(band,cl,habitat)][linepair].append((jaccard_similarity(aa,bb),a['source_id'],b['source_id']))
    summaries=[]
    for (band,cl,habitat),ls in sorted(groups.items()):
        summaries.append({'distance_km':band,'classes':cl,'habitat':habitat,
            'n_line_pairs':len(ls),'n_lines':len({s for pair in ls for s in pair}),
            'n_sites':len({s for vs in ls.values() for r in vs for s in r[1:]}),
            'n_pairs_not_independent':sum(map(len,ls.values())),
            'line_pair_equal_similarity':float(np.mean([np.mean([r[0] for r in v]) for v in ls.values()]))})
    matched={}
    for mode in ('geo_only','year_layout','year_layout_urban'):
        triples=[]
        for a in rows:
            candidates=[]
            for b in rows:
                if a['source_id']==b['source_id']:continue
                if mode!='geo_only' and (a['year'],a['layout'])!=(b['year'],b['layout']):continue
                d=haversine_km(a['lat'],a['lon'],b['lat'],b['lon'])
                if d<=40:candidates.append((d,b))
            allowed=[]
            for d,b in candidates:
                if a['section']!=b['section']:continue
                for e,c in candidates:
                    if a['section']==c['section'] or abs(d-e)>2.5:continue
                    if mode=='year_layout_urban' and abs(b['urban']-c['urban'])>.05:continue
                    allowed.append((abs(d-e),d+e,b['source_id'],c['source_id'],b,c))
            if allowed:
                best=min(allowed,key=lambda v:v[:4]);triples.append((a,best[-2],best[-1]))
        lines={a['section'] for a,b,c in triples}
        record={'n_anchors':len(triples),'n_lines':len(lines),
                'status':'supported' if len(triples)>=12 and len(lines)>=4 else 'not_identifiable'}
        if record['status']=='supported':
            for habitat in ('BS','BT','union'):
                by=defaultdict(list)
                for a,b,c in triples:
                    sets=[r['BS']|r['BT'] if habitat=='union' else r[habitat] for r in (a,b,c)]
                    by[a['section']].append(jaccard_similarity(sets[0],sets[1])-jaccard_similarity(sets[0],sets[2]))
                means=[float(np.mean(x)) for x in by.values()]
                record[habitat]={'same_minus_other':float(np.mean(means)),
                                 'positive_lines':sum(x>0 for x in means),
                                 'line_effects':means}
        matched[mode]=record
    return {'distance_bins_descriptive':summaries,'matched':matched}

def main():
    p=argparse.ArgumentParser()
    for name in ('source_xlsx','source_audit','urban','station_pdf','output'):
        p.add_argument('--'+name.replace('_','-'),required=True)
    p.add_argument('--draws',type=int,default=99)
    a=p.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.source_audit)
    grouped,qa=to_original_sites(original);sites=[r for rs in grouped.values() for r in rs]
    source,_=slovenia_location_sets(original)
    local,lqa=site_specific_null(source,draws=a.draws)
    traffic,tqa=mapped_segment_data(original,Path(a.station_pdf))
    raw20,qa20=eurostat_2020_segments()
    raw15,qa15=query_series('rail_tf_ns15_si')
    if qa15['sha256_original_json_response']!=EUROSTAT_SHA['2015'] or qa20['source']['sha256_original_json_response']!=EUROSTAT_SHA['2020']:
        raise ValueError('Historical traffic source drift')
    direction15,_=unpack_jsonstat(raw15,'2015')
    traffic15,historical_qa=historical_site_rematch(traffic,raw20,direction15)
    urban=json.loads(Path(a.urban).read_text(encoding='utf8'))
    if urban.get('complete') is not True or len(urban['sites'])!=89:
        raise ValueError('Incomplete urban acquisition')
    if urban['source_plant_sha256']!=hashlib.sha256(Path(a.source_xlsx).read_bytes()).hexdigest():
        raise ValueError('Urban and botanical source differ')
    out={'version':305,'plan':'f997281 + 13b0617','source_qa':qa,'local_qa':lqa,
         'traffic_qa':tqa,'draws':a.draws,'urban_sha256':hashlib.sha256(Path(a.urban).read_bytes()).hexdigest(),
         'historical_2015_qa':historical_qa,
         'scope':'One original railway study; BS and BT both railway habitats; observational associations',
         'buffers':{}}
    for buffer in BUFFERS:
        allrows=urban_rows(sites,urban,buffer)
        rows=[dict(r,**{k:v for k,v in local[r['source_id']].items() if k!='section'},
                   local_similarity=1-local[r['source_id']]['source_same_site_jaccard_beta'])
              for r in allrows if r['source_id'] in local]
        out['buffers'][buffer]={'valid_urban_sites':len(allrows),'local_support_sites':len(rows),
            'class':class_models(rows),'traffic':traffic_models(rows,traffic),
            'traffic_pairs':traffic_pair_models(allrows,traffic),
            'traffic_2015':traffic_models(rows,traffic15),
            'traffic_pairs_2015':traffic_pair_models(allrows,traffic15),
            'spatial':spatial(allrows)}
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
    print(json.dumps({'output':str(dest),'buffers':list(out['buffers']),'sites':len(sites),'draws':a.draws}))

if __name__=='__main__':main()

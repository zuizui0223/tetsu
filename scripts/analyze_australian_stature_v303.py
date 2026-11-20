"""Exploratory within-study elevation associations, with plot-area and habitat controls."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


def norm(s):
    return ' '.join(s.split())


def habitat(text):
    t=(text or '').lower()
    if any(s in t for s in ['drain','flood','creek','river']): return 'drainage_floodplain'
    if 'sandplain' in t: return 'sandplain'
    if 'plain' in t: return 'other_plain'
    return 'other'


def estimate(rows, adjusted=False):
    if len(rows)<10: return dict(n=len(rows),status='fewer_than_10_plots')
    y=np.array([r['log_height'] for r in rows]);x=np.array([r['elevation_m']/100 for r in rows])
    groups=sorted({r['habitat'] for r in rows})
    X=np.column_stack([np.ones(len(rows)),x]+([
        np.array([r['habitat']==g for r in rows],float) for g in groups[1:]] if adjusted else []))
    if np.linalg.matrix_rank(X)<X.shape[1]:return dict(n=len(rows),status='rank_deficient')
    coef=np.linalg.lstsq(X,y,rcond=None)[0][1]
    blocks=np.array([r['block'] for r in rows]);unique=sorted(set(blocks));loo=[]
    for b in unique:
        take=blocks!=b
        if np.linalg.matrix_rank(X[take])==X.shape[1]:loo.append(float(np.linalg.lstsq(X[take],y[take],rcond=None)[0][1]))
    rng=np.random.default_rng(303);boot=[]
    for _ in range(999):
        ids=np.concatenate([np.flatnonzero(blocks==b) for b in rng.choice(unique,len(unique),replace=True)])
        if np.linalg.matrix_rank(X[ids])==X.shape[1]:boot.append(float(np.linalg.lstsq(X[ids],y[ids],rcond=None)[0][1]))
    return dict(n=len(rows),blocks=len(unique),elevation_range_m=[float(x.min()*100),float(x.max()*100)],
        slope_log_height_per_100m=float(coef),leave_block_out_range=[min(loo),max(loo)] if loo else None,
        block_bootstrap_percentile_95=np.quantile(boot,[.025,.975]).tolist() if boot else None,
        valid_bootstraps=len(boot),status='exploratory_one_study',habitat_adjusted=adjusted)


def run(folder):
    folder=Path(folder);data=json.loads((folder/'wa_extracted.json').read_text(encoding='utf-8'))
    dem={r['site']:r['elevation_m'] for r in json.loads((folder/'wa_dem.json').read_text())['points']}
    h=pd.read_csv(folder/'austraits_height.csv',dtype=str)
    h=h[(h.trait_name=='plant_height')&(h.unit=='m')&(h.life_stage=='adult')&
        (h.entity_type=='species')&(h.value_type=='maximum')].copy()
    h['value']=pd.to_numeric(h.value,errors='coerce');h=h[h.value>0]
    h['taxon_name']=h.taxon_name.map(norm)
    # Prevent repeated observations from making a contributor dominate the taxon estimate.
    t=h.groupby(['taxon_name','dataset_id']).value.median().groupby(level=0).median().to_dict()
    sheet=defaultdict(set);all_sheet=defaultdict(set);matrix=defaultdict(set)
    for r in data['sheet_records']:
        all_sheet[r['site']].add(norm(r['taxon']))
        if not r['out_of_plot']:sheet[r['site']].add(norm(r['taxon']))
    for r in data['matrix']:matrix[r['site']].add(norm(r['taxon']))
    rows=[];excluded=[]
    for m in data['metadata']:
        site=m['site'];species=sheet[site];known=species&set(t)
        coverage=len(known)/len(species) if species else 0
        # Fixed area is read from the individual sheet, never inferred from richness.
        if m['area_m2']!=2500 or len(known)<5 or coverage<.5 or site not in matrix:
            excluded.append(dict(site=site,area=m['area_m2'],known=len(known),coverage=coverage));continue
        rows.append(dict(site=site,n_species=len(species),n_known=len(known),coverage=coverage,
            elevation_m=dem[site],log_height=float(np.mean([np.log(t[s]) for s in known])),
            habitat=habitat(m['habitat']),block=int(m['northing']//10000),
            matrix_sheet_complete=matrix[site]<=all_sheet[site]))
    panels={'all_2500m2':rows,'matrix_sheet_fully_reconciled':[r for r in rows if r['matrix_sheet_complete']],
            'coverage_75pct':[r for r in rows if r['coverage']>=.75]}
    for g in sorted({r['habitat'] for r in rows}):panels['habitat_'+g]=[r for r in rows if r['habitat']==g]
    out={name:estimate(rr) for name,rr in panels.items()};out['all_habitat_adjusted']=estimate(rows,True)
    # Fixed 50 m bands, retained regardless of estimated sign; no outcome-based cutpoints.
    for lower in range(400,550,50):
        band=[r for r in rows if lower<=r['elevation_m']<lower+50]
        out[f'elevation_band_{lower}_{lower+50}']=estimate(band)
        out[f'elevation_band_{lower}_{lower+50}_habitat_adjusted']=estimate(band,True)
    return dict(estimator='equal-presence mean log adult maximum height; AusTraits, contributor-balanced',
        area_m2=2500,unit='natural log metres per 100m elevation',source_studies=1,
        warning='railway reserve vegetation, not confirmed ballast; DEM elevation; exploratory associations',
        qa=dict(extraction=data['summary'],included=len(rows),excluded=len(excluded),
            coverage_range=[min(r['coverage'] for r in rows),max(r['coverage'] for r in rows)],
            fully_reconciled=sum(r['matrix_sheet_complete'] for r in rows)),results=out),rows


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--output',required=True)
    p.add_argument('--local-rows');a=p.parse_args();result,rows=run(a.data)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if a.local_rows:Path(a.local_rows).write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['results'],indent=2))

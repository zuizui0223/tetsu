#!/usr/bin/env python3
"""v25.5: source-frozen MIREN ROAD corridor reference-choice audit.

This reads only a previously source-verified derived JSON. MIREN roads are
NOT railways. All metrics are descriptive, not a transport intervention.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

SHA = '895fd24d2176ef1f19df9ba9eaab34f8326b7f7f'
PANELS = {'all_binomials': (0, 'jaccard'),
          'traits': (1, 'cwm'),
          'strict_traits': (1, 'cwm')}


def blob_sha(raw: bytes) -> str:
    return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()


def average(numbers):
    if not numbers: raise ValueError('Empty source group')
    return sum(numbers)/len(numbers)


def make_panel(item: dict, axis: int, metric: str) -> dict:
    rows=[]
    for region, roads in sorted(item['regions'].items()):
        if len(roads)<2 or any(r['sites']<4 for r in roads):
            raise ValueError('Unsupported original region/road site count')
        info={'region':region,'source_roads':len(roads),
              'source_sites':sum(r['sites'] for r in roads)}
        for mode in ('raw','equal'):
            means=[average([r[mode][zone][axis] for r in roads]) for zone in (0,1,2)]
            if not all(math.isfinite(x) for x in means):
                raise ValueError('Nonfinite original region average')
            info[mode]={'roadside_minus_middle':means[0]-means[1],
                        'roadside_minus_interior':means[0]-means[2],
                        'middle_minus_interior':means[1]-means[2],
                        'strictly_monotone':means[0]<means[1]<means[2]}
            if abs(info[mode]['roadside_minus_interior'] - info[mode]['roadside_minus_middle'] - info[mode]['middle_minus_interior'])>1e-12:
                raise AssertionError('Reference cancellation identity broken')
        rows.append(info)
    label={'raw':'raw','equal':'equal'}
    def sign(x):return -1 if x<0 else (1 if x>0 else 0)
    out={}
    for mode in label:
        by={k:[r[mode][k] for r in rows] for k in ('roadside_minus_middle','roadside_minus_interior','middle_minus_interior')}
        for key,values in by.items():
            original_key=mode+'_'+metric+'_'+('1minus2' if key=='roadside_minus_middle' else '1minus3' if key=='roadside_minus_interior' else '2minus3')
            ref=item['summary'].get(original_key)
            if ref is None or ref['n']!=len(rows) or abs(average(values)-ref['mean'])>1e-11:
                raise ValueError(f'Exact original MIREN source mean not replicated: {original_key}')
        out[mode]={'regional_mean_roadside_minus_middle':average(by['roadside_minus_middle']),
                   'regional_mean_roadside_minus_interior':average(by['roadside_minus_interior']),
                   'regional_mean_reference_only_middle_minus_interior':average(by['middle_minus_interior']),
                   'negative_roadside_minus_middle_regions':sum(x<0 for x in by['roadside_minus_middle']),
                   'negative_roadside_minus_interior_regions':sum(x<0 for x in by['roadside_minus_interior']),
                   'reference_sign_flip_regions':[r['region'] for r in rows if sign(r[mode]['roadside_minus_middle'])!=sign(r[mode]['roadside_minus_interior'])],
                   'strictly_monotone_regions':sum(r[mode]['strictly_monotone'] for r in rows)}
    out['n_regions']=len(rows)
    out['n_roads']=sum(r['source_roads'] for r in rows)
    out['n_original_sites']=sum(r['source_sites'] for r in rows)
    out['region_sign_change_rarefaction_roadside_vs_interior']=[r['region'] for r in rows if sign(r['raw']['roadside_minus_interior'])!=sign(r['equal']['roadside_minus_interior'])]
    out['by_region']=rows
    return out


def audit(j: dict) -> dict:
    original=j['results']
    outputs={}
    for name,(axis,metric) in PANELS.items():
        outputs[name]=make_panel(original[name],axis,metric)
    if (outputs['all_binomials']['n_regions']!=17 or outputs['all_binomials']['n_roads']!=50 or
        outputs['all_binomials']['n_original_sites']!=822):
        raise ValueError('Full original 17-region 50-road 822-site MIREN common support changed')
    if (outputs['traits']['n_regions']!=6 or outputs['strict_traits']['n_regions']!=4):
        raise ValueError('Original trait panel region support changed')
    return {'version':'v25.5','source_git_blob_sha1':SHA,'source_kind':'original-source-verified MIREN ROAD aggregated regional pairwise distances',
            'railway_direct_effect_identified':False,'distinguish_road_from_railway':True,
            'not_original_source_reingestion':True,'results':outputs}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raw=a.source.read_bytes()
    if blob_sha(raw)!=SHA:raise SystemExit('SHA-mismatched source-derived MIREN reference data')
    out=audit(json.loads(raw));a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:{'n_regions':v['n_regions'],'raw':v['raw'],'equal':v['equal'],
                            'reference_sign_flips_equal':v['equal']['reference_sign_flip_regions']}
                      for k,v in out['results'].items()},ensure_ascii=False))

if __name__=='__main__': main()

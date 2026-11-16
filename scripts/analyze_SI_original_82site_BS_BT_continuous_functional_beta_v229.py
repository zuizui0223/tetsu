#!/usr/bin/env python3
"""v22.9: authentic same-site paired Slovenian railway functional beta.

Hypothesis: BS has higher annual proportion but need not have less within-site
trait breadth; are different BS communities more functionally similar ACROSS
sites than the corresponding BT communities? The SAME unordered site pairs,
source geographic distances, survey years and true rail sections are used
for BOTH BS-BS and BT-BT comparisons.

CWM: incidence-equal centroid of species-level median z(log LEDA canopy
height), z(log LEDA seed mass); pairwise Euclidean centroid distance.
Auxiliary: absolute difference of author-defined strict annual proportions,
and Jaccard dissimilarity of EXACT SAME trait-complete species universe.
Negative (BS−BT) = greater convergence in BS for that measure.
Original microhabitat BS/BT survey widths differ; this is observational
composition, not selection, adaptation, or seed-transport flux.

Site cohort: v22.4 exactly 82 genuine localities from 89, 139 source
same-side matched 5+trait-species pairs; within-location sides unioned
for each zone. Real GPS midpoint from published 100m transects.

Pairwise pseudo-replication prevented by aggregation into 15 original
same-section means and 105 potentially independent-looking but dependent
original between-section pair means. No pairwise t tests or p-values.

Technical null: exchange BS and BT entire assemblage labels randomly
within each original site (999), preserving zone composition, source
positions, original sections and site-paired ecological differences. This
is not a randomized field treatment.
Richness equalization: sample min(nBS,nBT) source species in EACH zone of
each site for each of 199 technical draws and compute CWM distance again.
Source original date / genuine 2022–23 survey year and original GPS only.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx,validate_original_keys
)
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_sites_and_guilds,jaccard
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    original_gps_coordinates,haversine_km
)
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    HEIGHT_SHA,SEED_SHA,all_original_matched_taxa,ztransform,create_sites,
    components,quart
)

DRAWS=999
RARE_DRAWS=199
SEED=20261116+9
DIST_CUTS=(0,5,10,20,40,80,160,500)
METRICS=('functional_centroid_euclidean','annual_share_absolute_difference',
         'taxon_Jaccard_dissimilarity')

def site_zone(zset,values,group):
    if len(zset)<5:raise ValueError('Original source zone does not meet 5+ trait species')
    x=mean(values[t][0] for t in zset)
    y=mean(values[t][1] for t in zset)
    pa=mean(group[t] for t in zset)
    return {'species':set(zset),'centroid':(x,y),'annual_share':pa}

def zone_pair_stat(one,two):
    x1,y1=one['centroid'];x2,y2=two['centroid']
    return {
       'functional_centroid_euclidean':math.hypot(x1-x2,y1-y2),
       'annual_share_absolute_difference':abs(one['annual_share']-two['annual_share']),
       'taxon_Jaccard_dissimilarity':1-jaccard(one['species'],two['species'])}

def paired_site_groups(original,height,seed,minimum=5):
    locs,zones,sp,traits,cb,qc=validate_original_keys(original)
    names,nameqa=original_slovenian_EPPO_codebook_binomial_mapping(original)
    orig_universe=all_original_matched_taxa(zones,sp)
    v,zq=ztransform({t:binomial(names.get(t)) for t in orig_universe},height,seed)
    if len(v)!=224:raise ValueError('v22.4 original real 224 joined LEDA species changed')
    original_life={r['EPPO_code']:r for r in original['traits']}
    group={}
    for t in v:
        raw=original_life[t].get('Lifespan.categories')
        if raw is None or str(raw).strip()=='':
            raise ValueError('v22.7 confirmed every joint-trait species had source lifespan')
        group[t]=int(str(raw).strip().lower()=='annual')
    source_pairs,original_all_sides=create_sites(original,v,min_taxa=5)
    if minimum>5:
        source_pairs={site:[p for p in rows if p['nBS']>=minimum and p['nBT']>=minimum]
                      for site,rows in source_pairs.items()}
        source_pairs={site:rows for site,rows in source_pairs.items() if rows}
    qa={'all_publisher_matched_same_side_pairs':len(original_all_sides),
        'minimum_joint_trait_species_per_original_same_side_zone':minimum,
        'qualifying_joint_trait_matched_sides':sum(len(v) for v in source_pairs.values()),
        'qualifying_original_field_sites':len(source_pairs),
        'qualifying_original_line_sections':len({p['section'] for v in source_pairs.values() for p in v})}
    source_sites,guild,originalqa=source_sites_and_guilds(original)
    gps,gpsqa=original_gps_coordinates(original,source_sites)
    g={str(s['site']):s for s in gps}
    arr=[]
    for sid,prs in sorted(source_pairs.items()):
        A=set().union(*(r['BS'] for r in prs))
        B=set().union(*(r['BT'] for r in prs))
        if len(A)<5 or len(B)<5:raise ValueError('Source same-site union disobeyed 5+ support')
        field=g[sid]
        if any(p['section']!=str(field['section']) for p in prs):
            raise ValueError('Source site section inconsistent with original independent GPS')
        arr.append({'sid':sid,'section':str(field['section']),
            'year':field['year'],'track_type':field['track_type'],
            'lat':field['lat'],'lon':field['lon'],
            'BS':site_zone(A,v,group),'BT':site_zone(B,v,group),
            'source_same_side_pairs':len(prs)})
    expect={5:(82,15,139),10:(40,14,56)}
    if minimum in expect and (len(arr),len({s['section'] for s in arr}),
                              qa['qualifying_joint_trait_matched_sides'])!=expect[minimum]:
        raise ValueError('Author original matched sites/sections/side count drift under trait support threshold')
    return arr,qa,v,group

def all_site_pairs(sites):
    out=[]
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        dist=haversine_km(a['lat'],a['lon'],b['lat'],b['lon'])
        band=next((f'{lo}-{hi}' for lo,hi in zip(DIST_CUTS[:-1],DIST_CUTS[1:])
                   if lo<=dist<hi),None)
        if band is None:raise ValueError('Original geographical pair out of expected range')
        secpair=tuple(sorted((a['section'],b['section']),key=int))
        bb=zone_pair_stat(a['BS'],b['BS'])
        tt=zone_pair_stat(a['BT'],b['BT'])
        bt=zone_pair_stat(a['BS'],b['BT'])
        tb=zone_pair_stat(a['BT'],b['BS'])
        out.append({'i':i,'j':j,'sec_i':a['section'],'sec_j':b['section'],
                    'section_pair':secpair,'same_section':a['section']==b['section'],
                    'km':dist,'band':band,'normals':{m:bb[m]-tt[m] for m in METRICS},
                    'crosses':{m:bt[m]-tb[m] for m in METRICS}})
    if len(out)!=len(sites)*(len(sites)-1)//2:
        raise ValueError('Original independent locations pair count mismatch')
    return out

def pair_stat_value(pair,metric,flip_i=False,flip_j=False):
    a=pair['normals'][metric];x=pair['crosses'][metric]
    if not flip_i and not flip_j:return a
    if flip_i and flip_j:return -a
    if not flip_i and flip_j:return x
    return -x

def summary_by_independent_sections(pairs,metric,flips=None):
    per_sectionpair=defaultdict(list)
    per_band_sectionpair=defaultdict(list)
    raw=[]
    for p in pairs:
        x=pair_stat_value(p,metric,flips[p['i']],flips[p['j']]) if flips is not None else p['normals'][metric]
        raw.append(x)
        per_sectionpair[p['section_pair']].append(x)
        per_band_sectionpair[p['band'],p['section_pair']].append(x)
    equal_pair={s:mean(x) for s,x in per_sectionpair.items()}
    within=[v for (s,t),v in equal_pair.items() if s==t]
    cross=[v for (s,t),v in equal_pair.items() if s!=t]
    grouped=defaultdict(list)
    for (band,sectionpair),d in per_band_sectionpair.items():
        grouped[band].append(mean(d))
    return {'n_dependent_site_pairs':len(raw),
        'n_unique_same_source_section_pairs':len(within),
        'n_distinct_cross_source_section_pairs':len(cross),
        'mean_dependent_site_pair_delta':mean(raw),
        'mean_equal_within_section_delta':mean(within) if within else None,
        'mean_equal_cross_sectionpair_delta':mean(cross) if cross else None,
        'distance_strata_equal_section_pair':{
          k:{'n_original_section_pairs':len(v),'delta':mean(v)}
          for k,v in sorted(grouped.items())},
        'original_cross_sectionpair_leave_one_section_out':{
          'min':min(mean(v for (s,t),v in equal_pair.items()
              if s!=t and dropped not in (s,t))
              for dropped in sorted({z for st in equal_pair for z in st})),
          'max':max(mean(v for (s,t),v in equal_pair.items()
              if s!=t and dropped not in (s,t))
              for dropped in sorted({z for st in equal_pair for z in st}))} if len(cross)>10 else None}

def metric_aggregate(pairs,metric,flipvec=None):
    sec=defaultdict(list)
    for p in pairs:
        if p['same_section']:continue
        x=pair_stat_value(p,metric,flipvec[p['i']],flipvec[p['j']]) if flipvec is not None else p['normals'][metric]
        sec[p['section_pair']].append(x)
    return mean(mean(x) for x in sec.values())

def rarefy_site_once(source,vals,groups,rng):
    new=[]
    for r in source:
        k=min(len(r['BS']['species']),len(r['BT']['species']))
        zones={name:site_zone(set(rng.sample(sorted(r[name]['species']),k)),vals,groups)
               for name in ('BS','BT')}
        new.append({**r,**zones})
    return new

def difference_for_existing_pairs(rarefied,original_pairs):
    # All three beta metrics must be recomputed from exactly the same
    # rarefied original plant sets, not just continuous functional CWM.
    by={met:defaultdict(list) for met in METRICS}
    for p in original_pairs:
        if p['same_section']:continue
        a=rarefied[p['i']];b=rarefied[p['j']]
        bb=zone_pair_stat(a['BS'],b['BS'])
        tt=zone_pair_stat(a['BT'],b['BT'])
        for met in METRICS:
            by[met][p['section_pair']].append(bb[met]-tt[met])
    return {met:mean(mean(values) for values in by[met].values()) for met in METRICS}


def analyze(original,height,seed,draws=DRAWS,raredraws=RARE_DRAWS):
    sites,source_support,vals,groups=paired_site_groups(original,height,seed)
    pairs=all_site_pairs(sites)
    res={'version':'v22.9','status':'ORIGINAL_SI_CONTINUOUS_FUNCTIONAL_BETA',
         'author_source_DOI':'10.1007/s44473-026-00155-z',
         'source_joint_2d_trait_species':len(vals),
         'source_support':source_support,
         'source_locations_with_mapped_trait_beta':len(sites),
         'original_route_sections':len({s['section'] for s in sites}),
         'n_original_dependent_site_pairs':len(pairs),
         'original_species_annual_proportions_zone_means':{
           z:mean(s[z]['annual_share'] for s in sites) for z in ('BS','BT')},
         'source_zone_richness_mean_over_sites':{
           z:mean(len(s[z]['species']) for s in sites) for z in ('BS','BT')},
         'responses':{}}
    fliprng=random.Random(SEED)
    for met in METRICS:
        obs=summary_by_independent_sections(pairs,met)
        vals_null=[]
        for it in range(draws):
            flips=[bool(fliprng.getrandbits(1)) for _ in sites]
            vals_null.append(metric_aggregate(pairs,met,flips))
        score=obs['mean_equal_cross_sectionpair_delta']
        obs['site_pair_whole_zone_label_shuffle_null']={
           'n_draws':draws,'null_mean':mean(vals_null),
           'null_q025':quart(vals_null,.025),'null_q975':quart(vals_null,.975),
           'observed_minus_null_mean':score-mean(vals_null),
           'lower_tail_technical_rank_plusone':
               (1+sum(v<=score for v in vals_null))/(draws+1)}
        res['responses'][met]=obs
    rg=random.Random(SEED+33)
    rare_d={m:[] for m in METRICS}
    for _ in range(raredraws):
        rar=rarefy_site_once(sites,vals,groups,rg)
        scores=difference_for_existing_pairs(rar,pairs)
        for met in METRICS:rare_d[met].append(scores[met])
    res['sitewise_species_richness_equalized_all_beta']={
       met:{'technical_draws':raredraws,'mean_cross_sectionpair_delta':mean(v),
            'q025':quart(v,.025),'q975':quart(v,.975),
            'same_site_BS_BT_richness_equalized':True}
       for met,v in rare_d.items()}
    stricter,support10,unused1,unused2=paired_site_groups(original,height,seed,minimum=10)
    strictpairs=all_site_pairs(stricter)
    res['strict_both_zones_ge10_trait_species_source_sensitivity']={
       'source_support':support10,
       'mean_equal_cross_sectionpair_delta_by_metric':{
          met:summary_by_independent_sections(strictpairs,met)['mean_equal_cross_sectionpair_delta']
          for met in METRICS},
       'n_qualifying_locations':len(stricter),
       'n_source_sections':len({r['section'] for r in stricter})}
    res['interpretation_limits']=[
        'CWM centroid distance is not total functional beta hypervolume; within-site FDis/variance is separate',
        'Each BS-BT pair uses the same site GPS and rail section, so physical separation is exactly paired, but habitat sampling area remains potentially unequal',
        '15 original route sections and 105 overlapping section pairs are not 105 independent studies or field experiments',
        'Full-site zone label flips are a technical null and do not randomize field treatments',
        'Stricter >=10-species panels select a different 40-site community subset and are a sensitivity, not independent replication',
        'Rarefaction draws are technical site species sample resamplings, not sample-level confidence intervals',
        'LEDA traits are source species medians, not local individual acclimation or evolved phenotypes',
        'No freight or passenger train rate, departure count, historical operating service, direct dispersal or management pressure measured',
        'No independent nonrail regional natural species-pool reference exists at these plots']
    return res

def main():
    p=argparse.ArgumentParser()
    for v in ('original_xlsx','source_audit','height','seed','output'):
        p.add_argument('--'+v.replace('_','-'),required=True)
    p.add_argument('--draws',type=int,default=DRAWS)
    p.add_argument('--raredraws',type=int,default=RARE_DRAWS)
    a=p.parse_args()
    org=extract_original_xlsx(a.original_xlsx,a.source_audit)
    if hashlib.sha256(Path(a.height).read_bytes()).hexdigest()!=HEIGHT_SHA:
        raise ValueError('Original LEDA height SHA changed')
    if hashlib.sha256(Path(a.seed).read_bytes()).hexdigest()!=SEED_SHA:
        raise ValueError('Original LEDA seed mass SHA changed')
    h,hqa=read_leda(a.height,'height_m');s,sqa=read_leda(a.seed,'seed_mass_mg')
    ans=analyze(org,h,s,draws=a.draws,raredraws=a.raredraws)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ans,ensure_ascii=False,indent=2),encoding='utf8')
    summary={'version':ans['version'],'status':ans['status'],'source_support':ans['source_support'],
             'source_locations_with_mapped_trait_beta':ans['source_locations_with_mapped_trait_beta'],
             'original_species_annual_proportions_zone_means':ans['original_species_annual_proportions_zone_means'],
             'responses':ans['responses'],
             'sitewise_species_richness_equalized_all_beta':ans['sitewise_species_richness_equalized_all_beta'],
             'strict_both_zones_ge10_trait_species_source_sensitivity':ans['strict_both_zones_ge10_trait_species_source_sensitivity']}
    print('V229_REAL_CONTINUOUS_FUNCTIONAL_BETA='+json.dumps(summary,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

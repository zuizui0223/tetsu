#!/usr/bin/env python3
"""v22.7: decompose authentic Slovenian BS–BT functional breadth into lifespan mixture.

Reuses immutable 2026 publisher-author original species-zone XLSX, 2022-23
field sites, and independently source-verified public LEDA species-median
canopy height and seed-mass tables. 89 locations, 173 same-side BS/BT pairs,
15 original railway sections. This is observational assembly, NOT evolved
railway plant size or train-wind selection.

An exact population-variance identity separates joint (standardized log
height + log seed mass) trait variance into:
 pA*varA + (1-pA)*varOther + pA*(1-pA)*||meanA-meanOther||².
The first two are the within-life-history contribution; the last is
between-life-history separation. 'Other' includes any source-known
non-exact-annual lifespan; missing labels are excluded.
Primary considers ONLY source species with both LEDA traits and a
source-reported lifespan; min 5 per BS and BT zone. Retain an explicit
sensitivity requiring >=2 annual and >=2 other species per each zone,
and another excluding mixed lifespan descriptions.

MIXTURE-ONLY proxy fixes within-lifespan distributions to the UNION of
actually observed source species from a same-side paired BS∪BT pool,
allowing ONLY observed annual proportions to vary. This is an algebraic
decomposition, not the causal counterfactual of changing actual railway
selection.

Technical null: within each original same-side BS∪BT pool, hold shared
taxa fixed in both zones, preserve each zone's numbers of ANNUAL and
OTHER source taxa, shuffle ONLY the taxon assignments within classes.
It cannot replicate an independent nonrail/region species pool.
Equal sides within location, then equal locations within original railway
section, then equal 15 source sections. 499 Monte Carlo shuffles are
correlated technical draws, NOT independent biological experiments.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx,validate_original_keys
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    HEIGHT_SHA,SEED_SHA,all_original_matched_taxa,ztransform,variance,quart,
    site_effects,section_estimates
)

DRAWS=499
SEED=20261116+7
METRICS=('total_trace','within_annual_contribution',
         'within_nonannual_contribution','between_lifespan_contribution',
         'within_total','annual_share','height_mean','seed_mean')

def norm(v):
    return str(v).strip().lower() if v is not None else None

def group_lifespan(code, traits, exclude_mixed=False):
    original=norm(traits[code].get('Lifespan.categories'))
    if original is None or not original:return None
    if original=='annual':return 'annual'
    if exclude_mixed and ('annual' in original or 'biennial' in original):
        return None
    return 'other'

def component(codes, values, groups):
    if len(codes)<2:raise ValueError('Too few source species for full trait variance')
    a=[s for s in codes if groups[s]=='annual']
    o=[s for s in codes if groups[s]=='other']
    if len(a)+len(o)!=len(codes):raise ValueError('Missing original lifespan category')
    n=len(codes);p=len(a)/n
    def mu(v,k):
        return mean(values[s][k] for s in v) if v else 0.
    def vtrace(v):
        return sum(variance([values[s][k] for s in v]) or 0. for k in (0,1))
    va=vtrace(a) if a else 0.
    vo=vtrace(o) if o else 0.
    annual=p*va;nonannual=(1-p)*vo
    between=(p*(1-p)*sum((mu(a,k)-mu(o,k))**2 for k in (0,1))) if a and o else 0.
    joint=vtrace(codes)
    if not math.isclose(joint,annual+nonannual+between,rel_tol=1e-10,abs_tol=1e-9):
        raise AssertionError('Within/between law of total variance violated')
    return {
       'total_trace':joint,
       'within_annual_contribution':annual,
       'within_nonannual_contribution':nonannual,
       'between_lifespan_contribution':between,
       'within_total':annual+nonannual,
       'annual_share':p,
       'height_mean':mu(codes,0),'seed_mean':mu(codes,1),
       'raw_annual_n':len(a),'raw_other_n':len(o)}

def pair_metrics(a,b,values,groups):
    x=component(a,values,groups);y=component(b,values,groups)
    return {k:x[k]-y[k] for k in METRICS}

def mixture_only_proxy(a,b,values,groups):
    # Not observed difference: only p annual changes and the class-specific
    # means/variances are held at the same-site/side union of source plants.
    union=a|b
    x={g:{s for s in union if groups[s]==g} for g in ('annual','other')}
    nA=len(x['annual']);nO=len(x['other'])
    if nA+nO!=len(union):raise ValueError('Missing true source lifespan')
    if not nA or not nO:return 0.
    mu={g:[mean(values[s][d] for s in x[g]) for d in (0,1)] for g in x}
    vv={g:sum(variance([values[s][d] for s in x[g]]) or 0. for d in (0,1)) for g in x}
    bmu=sum((mu['annual'][d]-mu['other'][d])**2 for d in (0,1))
    def pred(p):return p*vv['annual']+(1-p)*vv['other']+p*(1-p)*bmu
    return pred(sum(groups[s]=='annual' for s in a)/len(a))-pred(sum(groups[s]=='annual' for s in b)/len(b))

def make_pairs(original,values,source_groups,minimum=5,within_stratum_ge2=False):
    loc,zones,species,traits,cb,qc=validate_original_keys(original)
    valid=set(values)&set(source_groups)
    paired=defaultdict(list);n_all=0;n_valid_species=set()
    for (lid,zid),record in zones.items():
        if record['area']!='BS' or zid not in ('LBS','RBS'):continue
        back=(lid,zid[0]+'BT')
        if back not in zones or zones[back]['area']!='BT':continue
        n_all+=1
        a=species.get((lid,zid),set())&valid
        b=species.get(back,set())&valid
        if len(a)<minimum or len(b)<minimum:continue
        if within_stratum_ge2:
            if any(sum(source_groups[t]==c for t in zone)<2
                   for zone in (a,b) for c in ('annual','other')):continue
        section=str(loc[lid]['track_section_number'])
        paired[str(lid)].append({
            'section':section,'BS':set(a),'BT':set(b),
            'known_joint_annual_nBS':sum(source_groups[t]=='annual' for t in a),
            'known_joint_annual_nBT':sum(source_groups[t]=='annual' for t in b),
            'known_joint_nBS':len(a),'known_joint_nBT':len(b),
            'shared_source_taxa':len(a&b)})
        n_valid_species|=a|b
    if n_all!=173:raise ValueError('Publisher original BS-BT 173 paired sides drifted')
    return paired,{'real_original_pairs_all':n_all,
         'qualified_pairs':sum(map(len,paired.values())),
         'qualified_sites':len(paired),
         'qualified_sections':len({p['section'] for pairs in paired.values() for p in pairs}),
         'qualified_species':len(n_valid_species),
         'mean_qualified_BS_known_species_per_side':mean(p['known_joint_nBS'] for v in paired.values() for p in v) if paired else None,
         'mean_qualified_BT_known_species_per_side':mean(p['known_joint_nBT'] for v in paired.values() for p in v) if paired else None,
         'mean_shared_species_per_side':mean(p['shared_source_taxa'] for v in paired.values() for p in v) if paired else None}

def site_contrast(pairs,values,groups,rng=None,permuted=False):
    bysite={}
    for site,rows in sorted(pairs.items()):
        diffs=[];proxies=[]
        for r in rows:
            a,b=r['BS'],r['BT']
            if permuted:
                common=a&b
                exclusive=a^b
                outA=set(common);outB=set(common)
                for grp in ('annual','other'):
                    avail=sorted(s for s in exclusive if groups[s]==grp)
                    n=sum(groups[t]==grp for t in a-common)
                    chosen=set(rng.sample(avail,n))
                    outA|=chosen;outB|=set(avail)-chosen
                a,b=outA,outB
                if len(a)!=r['known_joint_nBS'] or len(b)!=r['known_joint_nBT']:
                    raise AssertionError('Conditional fixed-margin lifespan shuffle changed zone richness')
                if len(a&b)!=r['shared_source_taxa']:raise AssertionError('Shared source taxa changed')
            diffs.append(pair_metrics(a,b,values,groups))
            if not permuted:proxies.append(mixture_only_proxy(a,b,values,groups))
        bysite[site]={
            'section':rows[0]['section'],
            'delta':{k:mean(d[k] for d in diffs) for k in METRICS},
            'mix_only':mean(proxies) if proxies else None,
            'n_original_matched_sides':len(rows)}
    return bysite

def bysection_aggregate(locations):
    d=defaultdict(list)
    for x in locations.values():d[x['section']].append(x)
    per={}
    for s,xs in d.items():
        per[s]={k:mean(x['delta'][k] for x in xs) for k in METRICS}
        if all(x['mix_only'] is not None for x in xs):
            per[s]['mix_only']=mean(x['mix_only'] for x in xs)
    global_vals={k:mean(x[k] for x in per.values()) for k in METRICS}
    global_vals['mix_only']=mean(x['mix_only'] for x in per.values() if 'mix_only' in x) if all('mix_only' in x for x in per.values()) else None
    return per,global_vals

def quantile_summary(values,obs=None):
    if not values:return None
    ret={'null_draws':len(values),'null_average':mean(values),
         'null_q025':quart(values,.025),'null_q975':quart(values,.975)}
    if obs is not None:
        ret['observed_minus_null_mean']=obs-ret['null_average']
        ret['technical_lower_tail_rank_plus_one']=(1+sum(x<=obs for x in values))/(len(values)+1)
    return ret

def analysis(original,height,seed,draws=DRAWS):
    loc,z,sp,original_lifespans,cb,base=validate_original_keys(original)
    observed_codes=all_original_matched_taxa(z,sp)
    names,nameqc=original_slovenian_EPPO_codebook_binomial_mapping(original)
    values,zqc=ztransform({t:binomial(names.get(t)) for t in observed_codes},height,seed)
    if len(values)!=224:raise ValueError('v22.0 true plant joint-LEDA names failed exact 224 support')
    original_trait_map={r['EPPO_code']:r for r in original['traits']}
    all_group={t:group_lifespan(t,original_trait_map) for t in observed_codes}
    groups={t:g for t,g in all_group.items() if t in values and g is not None}
    out={'version':'v22.7','source original BS/BT sites':89}
    out['source_joint_trait_species_n']=len(values)
    out['source_joint_trait_AND_lifespan_known_species_n']=len(groups)
    out['source_joint_trait_but_lifespan_missing_n']=len(values)-len(groups)
    result={}
    for panel in ('all_known_lifespans','both_lifespan_groups_ge2','mixed_lifespan_excluded'):
        gg=groups.copy()
        if panel=='mixed_lifespan_excluded':
            gg={t:group_lifespan(t,original_trait_map,exclude_mixed=True) for t in groups}
            gg={t:g for t,g in gg.items() if g is not None}
        pairs,qc=make_pairs(original,values,gg,minimum=5,within_stratum_ge2=(panel=='both_lifespan_groups_ge2'))
        if qc['qualified_sites']<20 or qc['qualified_sections']<8:
            result[panel]={'status':'INSUFFICIENT_SUPPORT','source_support':qc};continue
        obs=site_contrast(pairs,values,gg)
        section,score=bysection_aggregate(obs)
        q={
           'status':'OBSERVED','source_support':qc,
           'original_section_equal_delta_BS_minus_BT':score,
           'n_source_sections_positive_total_trace':sum(x['total_trace']>0 for x in section.values()),
           'n_source_sections_positive_annual_share':sum(x['annual_share']>0 for x in section.values()),
           'mean_effect_not_mix_only_residual':score['total_trace']-score['mix_only'] if score['mix_only'] is not None else None,
           'sections_pos_within_total':sum(x['within_total']>0 for x in section.values()),
           'sections_pos_between':sum(x['between_lifespan_contribution']>0 for x in section.values()),
        }
        # Only primary uses technical shuffle and all original source sections,
        # restrictive panels are explicit sensitivity not second main result.
        if panel=='all_known_lifespans' and qc['qualified_sections']==15:
            rg=random.Random(SEED)
            vals={k:[] for k in METRICS}
            for ix in range(draws):
                _,sample=bysection_aggregate(site_contrast(pairs,values,gg,rg,permuted=True))
                for k in METRICS:vals[k].append(sample[k])
            q['conditional_lifespan_fixed_localpool_shuffles']={
                k:quantile_summary(vals[k],score[k]) for k in METRICS}
            q['sampled_technical_shuffles_not_independent_field_p_values']=True
        result[panel]=q
    out['panels']=result
    out['method']='trace(Var of source species medians in z log height/log seed) = abundance-free annual/nonannual within trace + between group dispersion; equal sides->sites->route sections'
    out['inference_limits']=[
      'No direct measurements of individuals, genotype, mowing, ballast selection severity or train-count exposure',
      'Annual is source exact lifespan label; other comprises all source known non-exact annual categories',
      'Zone sampled widths and detectability differ; rarity thresholds change study sites',
      'Local BS union BT is a railway-internal observation, never an independent nonrail species pool',
      'Mixture-only proxy is algebraic standardization and does not identify causal mechanism',
      'The national set of 15 route sections are nested in a single original survey',
      'Technical shuffle ranks are not randomized field-trial p-values']
    return out

def main():
    p=argparse.ArgumentParser()
    for field in ('original_xlsx','source_audit','height','seed','output'):
        p.add_argument('--'+field.replace('_','-'),required=True)
    p.add_argument('--draws',type=int,default=DRAWS)
    a=p.parse_args()
    orig=extract_original_xlsx(a.original_xlsx,a.source_audit)
    if hashlib.sha256(Path(a.height).read_bytes()).hexdigest()!=HEIGHT_SHA:
        raise ValueError('LEDA source height exact SHA changed')
    if hashlib.sha256(Path(a.seed).read_bytes()).hexdigest()!=SEED_SHA:
        raise ValueError('LEDA source seed exact SHA changed')
    h,qh=read_leda(a.height,'height_m')
    se,qs=read_leda(a.seed,'seed_mass_mg')
    data=analysis(orig,h,se,a.draws)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf8')
    print('V227_REAL_ANNUAL_TRACE_DECOMPOSITION='+json.dumps(data,separators=(',',':'),ensure_ascii=False))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""v28.8 source-authentic India Kashmir rail/nonrail functional beta.
Year and original station sites are NOT independent publication replicates.
No estimates if a prespecified original-trait coverage gate fails.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import random
from pathlib import Path
from statistics import mean, median
from openpyxl import load_workbook
from analyze_Kashmir_same_sites_reference_year_v254 import (
    DRYAD_SHA256, YEARS, GROUPS, parse_xlsx_matrix, matching_sites,
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda, binomial

LEDA_HEIGHT_SHA = 'cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92'
LEDA_SEED_SHA = '34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef'
# v24.5 verified reference from three European studies; frozen BEFORE looking at India.
ZLOG_MU = (-0.5802365754696401, 0.2936113186962433)
ZLOG_SD = (1.3801000545693773, 2.1657660559995824)
ZONES = ('Track', 'Embankment', 'Interior Plot 1', 'Interior Plot 2')
REFERENCES = ('Interior Plot 1', 'Interior Plot 2')
MIN_DUAL_TRAIT_PER_ZONE = 5
MIN_COVERAGE_PER_ZONE = .30
MIN_SITES = {'station': 10, 'between': 8}
MIN_DISTINCT_TAXA = 20
N_RAREFACTION = 100
SEED = 288119

def strict_sha(path: Path, expected: str) -> None:
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError('Original author/LEDA source hash mismatch: '+str(path))

def species_headers(path: Path) -> tuple[str, ...]:
    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        hdr = next(wb.worksheets[0].values)
    finally:
        wb.close()
    if len(hdr) != 129 or len(hdr[6:]) != 123 or any(not isinstance(x,str) for x in hdr[6:]):
        raise ValueError('Author original Dryad 123 species columns do not match source schema')
    return tuple(hdr[6:])

def cov(plots:dict, trait:dict, names:tuple, site:str, year:int, zone:str) -> tuple[set[str],float,int]:
    source = plots[site,year,zone]
    if source is None:
        raise ValueError('Missing original source floristic observation')
    original = {binomial(names[i]) for i in source}
    original.discard(None)
    known = original & set(trait)
    return known, (len(known)/len(original) if original else 0.), len(original)

def metric(profiles:list[dict], field:str)->float:
    if len(profiles)<2: raise ValueError('Fewer than two independent original spatial sites')
    vals = []
    for i,a in enumerate(profiles):
        for b in profiles[i+1:]:
            if field=='joint':v=math.dist(a['joint'],b['joint'])
            elif field=='height':v=abs(a['joint'][0]-b['joint'][0])
            elif field=='seed':v=abs(a['joint'][1]-b['joint'][1])
            else:raise ValueError(field)
            vals.append(v)
    return mean(vals)

def site_profile(species:set[str],traits:dict) -> dict:
    return {'joint':(mean(traits[t][0] for t in species),mean(traits[t][1] for t in species))}

def contrasts(paired: list[dict], year:int, rarefy:bool=False, seed:int=SEED)->dict:
    rng=random.Random(seed)
    a=[]; b=[]
    for rec in paired:
        r,c=rec[year]['Track'],rec[year]['reference']
        if rarefy:
            target=min(len(r),len(c))
            r=set(rng.sample(sorted(r),target))
            c=set(rng.sample(sorted(c),target))
        a.append(site_profile(r, rec['traits']))
        b.append(site_profile(c, rec['traits']))
    return {field: {'rail':metric(a,field),'inland':metric(b,field),
                    'rail_minus_inland':metric(a,field)-metric(b,field),
                    'log_rail_over_inland': math.log(metric(a,field)/metric(b,field))
                    if metric(b,field)>0 and metric(a,field)>0 else None}
            for field in ('joint','height','seed')}

def eligible_pairs(plots, names, traits, group, reference):
    # SAME original matched sites for both years. Source-prespecified gate only.
    matched=matching_sites(plots,group,ZONES)
    paired=[]; audit=[]
    upper_bound_by_original_richness=0
    potential_if_no_30pct_coverage_gate=0
    for site in matched:
        rec={'traits':traits}
        allpass=True
        raw_possible=True
        known_possible=True
        for year in YEARS:
            r,cr,nr=cov(plots,traits,names,site,year,'Track')
            b,cb,nb=cov(plots,traits,names,site,year,reference)
            audit.append({'year':year,'rail_known':len(r),'inland_known':len(b),
                          'rail_coverage':cr,'inland_coverage':cb,
                          'rail_richness':nr,'inland_richness':nb})
            rec[year]={'Track':r,'reference':b}
            raw_possible &= min(nr,nb) >= MIN_DUAL_TRAIT_PER_ZONE
            known_possible &= min(len(r),len(b)) >= MIN_DUAL_TRAIT_PER_ZONE
            allpass &= (min(len(r),len(b)) >= MIN_DUAL_TRAIT_PER_ZONE and
                        min(cr,cb) >= MIN_COVERAGE_PER_ZONE)
        if raw_possible: upper_bound_by_original_richness+=1
        if known_possible: potential_if_no_30pct_coverage_gate+=1
        if allpass:paired.append(rec)
    n_taxa=(len(set().union(*(rec[y][z] for rec in paired for y in YEARS
                              for z in ('Track','reference')))) if paired else 0)
    status=('ESTIMABLE' if len(paired)>=MIN_SITES[group] and n_taxa>=MIN_DISTINCT_TAXA
            else 'NOT_IDENTIFIABLE_TRAIT_SUPPORT')
    aggregate={
        'original_sites_with_both_field_years_and_four_zones':len(matched),
        'n_original_sites_meeting_predeclared_trait_gate_both_years':len(paired),
        'absolute_max_eligible_sites_even_with_100pct_trait_coverage_given_observed_plot_richness':upper_bound_by_original_richness,
        'sites_with_observed_5plus_joint_trait_species_both_years_ignoring_coverage_fraction':potential_if_no_30pct_coverage_gate,
        'original_rail_plot_raw_richness_min_median_max':[min(x['rail_richness'] for x in audit),median(x['rail_richness'] for x in audit),max(x['rail_richness'] for x in audit)],
        'original_rail_plot_matched_joint_trait_species_min_median_max':[min(x['rail_known'] for x in audit),median(x['rail_known'] for x in audit),max(x['rail_known'] for x in audit)],
        'original_inland_plot_raw_richness_min_median_max':[min(x['inland_richness'] for x in audit),median(x['inland_richness'] for x in audit),max(x['inland_richness'] for x in audit)],
        'n_distinct_original_binomials_with_both_traits_in_eligible_sites':n_taxa,
        'min_site_n':MIN_SITES[group],
        'min_taxa_n':MIN_DISTINCT_TAXA,
        'status':status,
        'site_count_source_QA_only':{
            'low_coverage_site_years':sum(min(x['rail_coverage'],x['inland_coverage'])
                                          <MIN_COVERAGE_PER_ZONE for x in audit),
            'low_taxon_site_years':sum(min(x['rail_known'],x['inland_known'])
                                       <MIN_DUAL_TRAIT_PER_ZONE for x in audit)}}
    if status=='ESTIMABLE':
        vals={str(year): contrasts(paired,year) for year in YEARS}
        rare={str(year):[contrasts(paired,year,True,SEED+year*1000+i)
                         for i in range(N_RAREFACTION)] for year in YEARS}
        aggregate['observed_rail_minus_inland_CWM_between_site_beta']=vals
        aggregate['equal_trait_known_species_per_paired_site_zone_technical_draws']={
            str(year):{field:{'mean_delta':mean(v[field]['rail_minus_inland']
                                                  for v in rare[str(year)]),
              'negative_draw_fraction':sum(v[field]['rail_minus_inland']<0
                                             for v in rare[str(year)])/N_RAREFACTION}
             for field in ('joint','height','seed')} for year in YEARS}
        aggregate['leave_one_original_site_out_signs']={
            str(year):{field:{'negative':sum(
                contrasts([v for j,v in enumerate(paired) if j!=i],year)[field]['rail_minus_inland']<0
                for i in range(len(paired))),
                'n':len(paired)} for field in ('joint','height','seed')} for year in YEARS}
    return aggregate

def run(dryad:Path,height:Path,seed:Path)->dict:
    for p,h in ((dryad,DRYAD_SHA256),(height,LEDA_HEIGHT_SHA),(seed,LEDA_SEED_SHA)):
        strict_sha(p,h)
    plots,provenance=parse_xlsx_matrix(dryad)
    names=species_headers(dryad)
    if (provenance['all_spatial_sites']!=31 or
        provenance['plot_year_zone_rows']!=248 or
        provenance['source_missing_species_cells']!=4):
        raise ValueError('Original 31 sites/248 zone-by-year rows/4 missing-source cells changed')
    h,hqa=read_leda(height,'height_m')
    s,sqa=read_leda(seed,'seed_mass_mg')
    full={binomial(x) for x in names};full.discard(None)
    joined=full & set(h) & set(s)
    traits={t:((math.log(h[t])-ZLOG_MU[0])/ZLOG_SD[0],
               (math.log(s[t])-ZLOG_MU[1])/ZLOG_SD[1]) for t in joined}
    out={'version':'v28.8','study':'Rashid et al 2021 original 2014/2017 Kashmir',
         'original_Dryad':provenance,
         'source_original_123_taxa_strict_binomial_n':len(full),
         'source_2D_LEDA_both_trait_species_n':len(joined),
         'source_2D_LEDA_trait_species_fraction':len(joined)/len(full),
         'original_LEDA_source_QA':{'height_sha':hqa['sha256_actual_download'],
             'seed_sha':sqa['sha256_actual_download'],
             'height_trait_species':hqa['unique_positive_trait_species'],
             'seed_trait_species':sqa['unique_positive_trait_species']},
         'predeclared_site_min_dual_trait_n':MIN_DUAL_TRAIT_PER_ZONE,
         'predeclared_zone_min_trait_coverage':MIN_COVERAGE_PER_ZONE,
         'two_years_are_repeated_survey_of_same_original_sites_not_independent_studies':True,
         'author_source_plot_area_m2_per_zone':100,
         'reference_is_actual_adjacent_inland_field_zone_not_confirmed_gravel_habitat':True,
         'traits_are_original_species_medians_not_local_phenotypes':True,
         'source_four_zone_repeated_field_site_panels':{}}
    for group in GROUPS:
        out['source_four_zone_repeated_field_site_panels'][group]={}
        for ref in REFERENCES:
            out['source_four_zone_repeated_field_site_panels'][group][ref]=eligible_pairs(
                plots,names,traits,group,ref)
    out['n_estimable_rail_nonrail_panels']=sum(
        x['status']=='ESTIMABLE'
        for grp in out['source_four_zone_repeated_field_site_panels'].values()
        for x in grp.values())
    out['note']='no study-level meta-analysis, no causal adaptation, seed dispersal, traffic or herbicide effect'
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--height',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=run(a.source,a.height,a.seed)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps({'original_n_taxa':result['source_original_123_taxa_strict_binomial_n'],
                      'both_LEDA_n_taxa':result['source_2D_LEDA_both_trait_species_n'],
                      'n_estimable_panels':result['n_estimable_rail_nonrail_panels'],
                      'statuses':{g:{r:x['status'] for r,x in panels.items()}
                                  for g,panels in result['source_four_zone_repeated_field_site_panels'].items()}},indent=2))

if __name__=='__main__': main()

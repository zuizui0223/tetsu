#!/usr/bin/env python3
"""v29.2: does MIREN road/interior height-sorting survive a taxonomic genus null?

Original site×taxon incidence FROZEN; no species records shuffled. The
unrestricted/occupancy-conditioned trait-label null in v29.1 may detect genus
turnover and phylogenetic trait conservation rather than within-genus signal.
Swap each genuine species' PAIRED height/seed traits with other observed
species of the SAME original genus. A genus with only one eligible taxon is
fully fixed. Report coverage and fail closed when too few multispecies genera
or source incidences are exchangeable.

One true MIREN observational field programme, not independent countries or
a fifth railway field dataset. Technical draws are not ecological p-values.
"""
from __future__ import annotations
import argparse, hashlib, json, math, random, re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
    read_source_miren, traits_from_source, regional_comparison, MIREN_DOI,
)
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 import (
    source_mean_distances, selected_species_with_frequency, PANEL_DEFINITIONS,
    PINNED_PRIOR_ROAD_RESULT, METRICS, FROZEN_3_STUDY_LOG_MEANS,
    FROZEN_3_STUDY_LOG_SDS,
)
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA

DRAWS_PER_CHAIN=199
SEEDS=(292119,292120)
MIN_EXCHANGEABLE_GENERA_PER_REGION=10
MIN_EXCHANGEABLE_INCIDENCE_BY_HABITAT=.25
MIN_REGION_REPLICATES_FOR_SUMMARY=3
FROZEN_V291_PRIMARY_RAW_HEIGHT_ROAD_MINUS_INTERIOR=-0.0931493207367427
GENUS_RX=re.compile(r'^[A-Z][A-Za-z-]+$')

def genus(taxon: str)->str:
    p=taxon.split()
    if len(p)!=2 or not GENUS_RX.fullmatch(p[0]) or not p[1].islower():
        raise ValueError('Original binomial taxonomy is ambiguous: '+repr(taxon))
    return p[0]

def genus_groups(taxa:set[str])->list[list[str]]:
    by=defaultdict(list)
    for taxon in sorted(taxa):
        by[genus(taxon)].append(taxon)
    return [sorted(v) for _,v in sorted(by.items())]

def genus_support(roads,traits)->dict:
    occ=selected_species_with_frequency(roads)
    if not set(occ).issubset(traits):raise ValueError('Original taxon trait missing')
    groups=genus_groups(set(occ))
    movable={t for g in groups if len(g)>=2 for t in g}
    incidences={}
    for zone in ('roadside','interior'):
        original_total=sum(len(record[zone]) for _,sites in roads for record in sites)
        mobile=sum(len(record[zone]&movable) for _,sites in roads for record in sites)
        if original_total<=0:raise ValueError('No original source field observations in habitat')
        incidences[zone]={'n_original_trait_known_site_taxon_incidence':original_total,
            'n_original_trait_known_incidence_in_multispecies_genus':mobile,
            'exchangeable_incidence_fraction':mobile/original_total}
    n_gen=sum(len(g)>=2 for g in groups)
    eligible=(n_gen>=MIN_EXCHANGEABLE_GENERA_PER_REGION and
              min(x['exchangeable_incidence_fraction'] for x in incidences.values())
              >=MIN_EXCHANGEABLE_INCIDENCE_BY_HABITAT)
    # This gate depends on source taxonomy/trait support only, never observed CWM effects.
    return {
        'n_original_joint_trait_species_in_region':len(occ),
        'n_original_genera':len(groups),
        'n_multispecies_genera':n_gen,
        'n_exchangeable_joint_trait_species':len(movable),
        'n_immobile_singleton_genus_joint_trait_species':len(occ)-len(movable),
        'n_distinct_trait_vector_multispecies_genera':sum(
            len(g)>=2 and len({traits[t] for t in g})>=2 for g in groups),
        'original_site_habitat_trait_known_taxon_incidence':incidences,
        'predeclared_min_multispecies_genera':MIN_EXCHANGEABLE_GENERA_PER_REGION,
        'predeclared_min_exchangeable_incidence_fraction_each_habitat':
             MIN_EXCHANGEABLE_INCIDENCE_BY_HABITAT,
        'genus_conditioned_null_identifiable':eligible,
        'status':'ELIGIBLE_GENUS_CONDITIONED_NULL' if eligible
                  else 'NOT_IDENTIFIABLE_GENUS_CONDITIONED_TRAIT_SHUFFLE',
    }

def permute_within_genus(traits:dict,groups:list[list[str]],rng:random.Random)->dict:
    answer=dict(traits)
    for names in groups:
        if len(names)<2:continue
        donors=names.copy()
        rng.shuffle(donors)
        for target,origin in zip(names,donors):
            answer[target]=traits[origin]
    if set(answer)!=set(traits) or sorted(answer.values())!=sorted(traits.values()):
        raise ValueError('Original genus null did not conserve all bivariate trait values')
    for names in groups:
        if sorted(answer[t] for t in names)!=sorted(traits[t] for t in names):
            raise ValueError('Paired traits moved across genus boundaries')
    return answer

def genus_placebo(roads,traits,draws:int=DRAWS_PER_CHAIN)->dict:
    if draws<5:raise ValueError('Too few technical draws')
    support=genus_support(roads,traits)
    if not support['genus_conditioned_null_identifiable']:
        return {'source_support':support,'trait_beta_effects':None,
                'status':'NOT_IDENTIFIABLE_GENUS_CONDITIONED_TRAIT_SHUFFLE'}
    occ=selected_species_with_frequency(roads)
    observed_traits={t:traits[t] for t in occ}
    groups=genus_groups(set(occ))
    observed=source_mean_distances(roads,observed_traits)
    draws_by_zone={z:{m:[] for m in METRICS}
                   for z in ('roadside','interior','road_minus_interior')}
    for seed in SEEDS:
        rng=random.Random(seed+len(occ))
        for _ in range(draws):
            shuffled=permute_within_genus(observed_traits,groups,rng)
            trial=source_mean_distances(roads,shuffled)
            for z in draws_by_zone:
                for m in METRICS:
                    draws_by_zone[z][m].append(trial[z][m])
    effects={}
    for z in draws_by_zone:
        effects[z]={}
        for m,drawn in draws_by_zone[z].items():
            actual=observed[z][m]
            estimate=mean(drawn)
            effects[z][m]={
              'original_observed':actual,
              'within_genus_trait_label_null_mean':estimate,
              'observed_minus_within_genus_null':actual-estimate,
              'technical_upper_rank_not_field_p':
                    (1+sum(v>=actual for v in drawn))/(len(drawn)+1),
              'n_technical_label_assignments':len(drawn)}
    return {
        'source_support':support,
        'trait_beta_effects':effects,
        'status':'SOURCE_GENUS_CONDITIONED_EMPIRICAL_ESTIMABLE',
        'original_site_taxon_presence_frequency_and_community_fixed':True,
        'both_traits_remain_paired_within_same_genus':True,
    }

def region_meta(regions:dict,values:dict,draws:int)->dict:
    by_reg={}
    for reg,roads in sorted(regions.items()):
        s=genus_placebo(roads,values,draws)
        by_reg[reg]=s
    ok={r:v for r,v in by_reg.items()
        if v['status']=='SOURCE_GENUS_CONDITIONED_EMPIRICAL_ESTIMABLE'}
    stats={'original_region_codes_all':len(regions),
           'source_genus_support_eligible_n':len(ok),
           'source_region_codes_eligible':sorted(ok),
           'original_source_region_fails_support_n':len(regions)-len(ok),
           'no_pooled_genus_effect_when_fewer_than_3_eligible_regions':
                 len(ok)<MIN_REGION_REPLICATES_FOR_SUMMARY,
           'descriptive_effects':None}
    if len(ok)>=MIN_REGION_REPLICATES_FOR_SUMMARY:
        out={}
        for metric in METRICS:
            for zone in ('roadside','interior','road_minus_interior'):
                key=f'{zone}_{metric}'
                eff={r:v['trait_beta_effects'][zone][metric][
                    'observed_minus_within_genus_null'] for r,v in ok.items()}
                loo=[mean(v for r,v in eff.items() if r!=removed)
                     for removed in eff] if len(eff)>=4 else []
                out[key]={
                    'source_region_equal_mean_genus_conditioned_excess':
                        mean(eff.values()),
                    'n_original_region_codes_genus_conditional_positive':
                        sum(x>0 for x in eff.values()),
                    'n_original_region_codes_genus_conditional_negative':
                        sum(x<0 for x in eff.values()),
                    'original_regions':len(eff),
                    'source_loo_region_equal_mean_range':
                        [min(loo),max(loo)] if loo else None,
                    'all_original_region_code_signs':eff}
        stats['descriptive_effects']=out
    return {'source_region_taxonomy_support':stats,
            'source_region_aggregates_only':by_reg}

def run(miren:Path,height:Path,seed:Path,draws:int=DRAWS_PER_CHAIN)->dict:
    for p,key in ((height,'height'),(seed,'seed')):
        if hashlib.sha256(p.read_bytes()).hexdigest()!=LEDASHA[key]:
            raise ValueError('Original LEDA trait source digest changed')
    source,src_qa=read_source_miren(miren)
    h,hqa=read_leda(height,'height_m')
    s,sqa=read_leda(seed,'seed_mass_mg')
    prior_value,prior_qa=traits_from_source(source,h,s)
    if len(prior_value)!=1096 or src_qa['original_source_record_QA'][
        'original_source_species_observation_rows']!=167804:
        raise ValueError('Original global MIREN source joint trait support changed')
    traits={taxon:((math.log(h[taxon])-FROZEN_3_STUDY_LOG_MEANS[0])/
                   FROZEN_3_STUDY_LOG_SDS[0],
                  (math.log(s[taxon])-FROZEN_3_STUDY_LOG_MEANS[1])/
                   FROZEN_3_STUDY_LOG_SDS[1])
            for taxon in prior_value}
    panels={}
    for label,(threshold,share,minyear,nregions,nroads,nsites) in PANEL_DEFINITIONS.items():
        regions,opportunity=valid_region_panels(source,prior_value,threshold,share,minyear)
        actual_support=(len(regions),sum(len(v) for v in regions.values()),
            sum(len(sites) for group in regions.values() for _,sites in group))
        if actual_support!=(nregions,nroads,nsites):
            raise ValueError('Previously frozen MIREN v29.1 original panel support drift')
        if label=='primary_5species_25percent':
            m=mean(regional_comparison(roads,prior_value,reg)[
                'source_2D_CWM_beta_road_minus_interior'] for reg,roads in regions.items())
            if abs(m-PINNED_PRIOR_ROAD_RESULT)>1e-9:
                raise ValueError('v24.9 genuine original MIREN 2D road beta drift')
            height_obs=mean(source_mean_distances(roads,traits)[
                'road_minus_interior']['height'] for roads in regions.values())
            if abs(height_obs-FROZEN_V291_PRIMARY_RAW_HEIGHT_ROAD_MINUS_INTERIOR)>1e-9:
                raise ValueError('v29.1 frozen Europe-axis true roadside height beta drift')
        panels[label]={
            'original_source_site_habitat_support':opportunity,
            'genus_conditioned_placebo':region_meta(regions,traits,draws)}
    return {
        'version':'v29.2','source':'Genuine MIREN global mountain ROADS',
        'source_DOI':MIREN_DOI,'source_count_QA':src_qa,
        'source_original_LEDA_QA':{'joint_taxa':len(traits),
                'height_sha':hqa['sha256_actual_download'],
                'seed_sha':sqa['sha256_actual_download']},
        'predeclared_original_source_species_site_layout_AND_Europe_Z_axis':
            True,
        'panels':panels,
        'original_genus_phylogeny_proxy_not_true_phylogeny':True,
        'null_retains_taxon_identity_and_species_locations':True,
        'not_independent_replicate_of_v29.1':True,
        'no_field_adaptation_or_transport_inferred':True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--miren',required=True,type=Path)
    p.add_argument('--height',required=True,type=Path)
    p.add_argument('--seed',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--draws',default=DRAWS_PER_CHAIN,type=int)
    a=p.parse_args()
    out=run(a.miren,a.height,a.seed,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({'version':out['version'],
      'panels':{k:{
        'support':v['genus_conditioned_placebo']['source_region_taxonomy_support'],
        'region_genus_support':{r:z['source_support'] for r,z in
            v['genus_conditioned_placebo']['source_region_aggregates_only'].items()}
        } for k,v in out['panels'].items()}
        },indent=2,ensure_ascii=False))
if __name__=='__main__': main()

#!/usr/bin/env python3
"""v29.1: genuine MIREN ROAD negative control for Europe railway v28.5.

Freeze the original species x site x habitat incidence matrix. Reassign only
paired canopy-height/seed-mass species-trait labels *within source regions*.
Do roadside zones show a trait-spatial signature beyond a random label
assignment and, crucially, is the same seen in interior vegetation?

Not a fifth railway study. ROAD/INTERIOR regions are related MIREN field
programs, not independent randomized ecosystem experiments.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
from statistics import mean

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
    read_source_miren, traits_from_source, regional_comparison,
    MIREN_MD5, MIREN_DOI, MIREN_STUDY_DOI,
)
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import (
    FROZEN_3_STUDY_LOG_MEANS, FROZEN_3_STUDY_LOG_SDS,
)

PANEL_DEFINITIONS={
    'primary_5species_25percent':(5,.25,2007,8,22,272),
    'strict_10species_50percent':(10,.50,2007,5,13,192),
}
METRICS=('joint','height','seed')
MODES=('region_trait_labels_unrestricted','region_trait_labels_occupancy_quintile')
SEEDS=(291119,291120)
DRAWS_PER_CHAIN=99
PINNED_PRIOR_ROAD_RESULT=-0.10381093795211047

def source_mean_distances(roads: list, vectors:dict) -> dict:
    """Each original road gets one vote, then each original MIREN region."""
    by={'roadside':{k:[] for k in METRICS},
        'interior':{k:[] for k in METRICS}}
    for road,sites in roads:
        if len(sites)<2: raise ValueError('Need independent field sites within original road')
        current={'roadside':{k:[] for k in METRICS},
                 'interior':{k:[] for k in METRICS}}
        for zone in ('roadside','interior'):
            profiles=[]
            for record in sites:
                if not record[zone] or not set(record[zone]).issubset(vectors):
                    raise ValueError('Missing original joint trait pool for source observation')
                profiles.append((
                    mean(vectors[t][0] for t in record[zone]),
                    mean(vectors[t][1] for t in record[zone])))
            for i,(a,b) in enumerate(profiles):
                for c,d in profiles[i+1:]:
                    current[zone]['joint'].append(math.hypot(a-c,b-d))
                    current[zone]['height'].append(abs(a-c))
                    current[zone]['seed'].append(abs(b-d))
        for zone in current:
            for metric in METRICS:
                by[zone][metric].append(mean(current[zone][metric]))
    if not roads:raise ValueError('No roads in field region')
    result={zone:{metric:mean(vals) for metric,vals in part.items()}
            for zone,part in by.items()}
    result['road_minus_interior']={
        metric:result['roadside'][metric]-result['interior'][metric]
        for metric in METRICS}
    return result

def selected_species_with_frequency(roads:list)->Counter:
    occ=Counter()
    for road,sites in roads:
        for rec in sites:
            for zone in ('roadside','interior'):
                occ.update(rec[zone])
    return occ

def occupancy_quintile_groups(occ:Counter)->list[list[str]]:
    """Disjoint frequency-ordered bins, no post-hoc traits in grouping."""
    names=sorted(occ,key=lambda t:(occ[t],t))
    if len(names)<10:raise ValueError('Region has too few observed trait species')
    n=len(names)
    groups=[names[(n*i)//5:(n*(i+1))//5] for i in range(5)]
    if min(map(len,groups))<2:raise ValueError('Too few original taxa per occupancy bin')
    if len(set().union(*(set(x) for x in groups)))!=n:
        raise ValueError('Taxa not partitioned uniquely')
    return groups

def permute_region_labels(vectors:dict,groups:list[list[str]],rng:random.Random)->dict:
    """Source species/site incidence unchanged; paired 2D vectors never split."""
    keys=set(vectors)
    result=dict(vectors)
    for group in groups:
        if not set(group).issubset(keys):raise ValueError('Stratum contains unknown source species')
        labels=list(group)
        rng.shuffle(labels)
        for target, donor in zip(group,labels):
            result[target]=vectors[donor]
    if sorted(result.values())!=sorted(vectors.values()):
        raise ValueError('Paired original traits were not preserved in multiset')
    return result

def region_null(roads:list,traits:dict,draws:int=DRAWS_PER_CHAIN)->dict:
    if draws<4:raise ValueError('Too few technical label assignments')
    occ=selected_species_with_frequency(roads)
    if not set(occ).issubset(traits):raise ValueError('Original trait missingness drift')
    local={taxon:traits[taxon] for taxon in occ}
    observed=source_mean_distances(roads,local)
    groupdefs={
        MODES[0]:[sorted(occ)],
        MODES[1]:occupancy_quintile_groups(occ),
    }
    null={}
    for mode,groups in groupdefs.items():
        chains=[]
        for chain_seed in SEEDS:
            rng=random.Random(chain_seed+sum(ord(c) for c in mode)+len(occ))
            draws_values={z:{k:[] for k in METRICS} for z in
                          ('roadside','interior','road_minus_interior')}
            for _ in range(draws):
                trial=source_mean_distances(roads,permute_region_labels(local,groups,rng))
                for zone in draws_values:
                    for key in METRICS:
                        draws_values[zone][key].append(trial[zone][key])
            chains.append(draws_values)
        mode_out={}
        for zone in ('roadside','interior','road_minus_interior'):
            mode_out[zone]={}
            for key in METRICS:
                values=chains[0][zone][key]+chains[1][zone][key]
                mu=mean(values);real=observed[zone][key]
                mode_out[zone][key]={
                    'original_observed':real,
                    'trait_label_null_mean':mu,
                    'observed_minus_trait_label_null':real-mu,
                    'n_technical_assignments_not_field_studies':len(values),
                    'technical_upper_rank_not_biological_p':
                        (1+sum(v>=real for v in values))/(len(values)+1),
                }
        null[mode]={'effects':mode_out,'n_original_bins':len(groups)}
    return {
        'n_source_original_roads':len(roads),
        'n_source_field_sites':sum(len(x) for _,x in roads),
        'n_joint_trait_species_original_region_pool':len(local),
        'source_original_species_site_habitat_incidence_frozen':True,
        'source_bivariate_trait_correlation_preserved_per_assignment':True,
        'source_null_by_mode':null,
    }

def run(miren:Path,height:Path,seed:Path,draws:int=DRAWS_PER_CHAIN)->dict:
    for src,exp in ((height,LEDASHA['height']),(seed,LEDASHA['seed'])):
        if hashlib.sha256(src.read_bytes()).hexdigest()!=exp:
            raise ValueError('Original Oldenburg LEDA source SHA mismatch')
    rec,miren_source=read_source_miren(miren)
    h,hqa=read_leda(height,'height_m')
    s,sqa=read_leda(seed,'seed_mass_mg')
    original_v249_values,coverage=traits_from_source(rec,h,s)
    if len(original_v249_values)!=1096:
        raise ValueError('Original MIREN source to LEDA species source count changed')
    if miren_source['original_source_record_QA']['original_source_species_observation_rows']!=167804:
        raise ValueError('Genuine original MIREN observation row count changed')
    calibrated={
        taxon:(
            (math.log(h[taxon])-FROZEN_3_STUDY_LOG_MEANS[0])/FROZEN_3_STUDY_LOG_SDS[0],
            (math.log(s[taxon])-FROZEN_3_STUDY_LOG_MEANS[1])/FROZEN_3_STUDY_LOG_SDS[1],
        )
        for taxon in original_v249_values}
    result={
       'version':'v29.1',
       'original_DOI':MIREN_DOI,'original_publication_DOI':MIREN_STUDY_DOI,
       'scope':'MIREN 2007-2023 global mountain ROADS; not railway botanical sample',
       'original_source_QA':miren_source,
       'joint_trait_original_coverage_QA':coverage,
       'source_original_LEDA_sha':{'height':hqa['sha256_actual_download'],
                                  'seed':sqa['sha256_actual_download']},
       'frozen_Europe_v245_zlog_axis_for_out_of_source_transfer':{
          'mu':FROZEN_3_STUDY_LOG_MEANS,'sd':FROZEN_3_STUDY_LOG_SDS},
       'strict_source_rail_vs_road_causal_effect_not_identified':True,
       'source_plant_taxon_sites_not_exported':True,
       'panels':{}
    }
    for panel,(minimum,share,year,nr,nroads,nsites) in PANEL_DEFINITIONS.items():
        regions,support=valid_region_panels(rec,original_v249_values,minimum,share,year)
        if (len(regions),sum(len(roads) for roads in regions.values()),
           sum(len(x) for group in regions.values() for _,x in group))!=(nr,nroads,nsites):
            raise ValueError('Previously verified v25.0 MIREN region/road/site support drift: '+panel)
        if panel=='primary_5species_25percent':
            m=mean(regional_comparison(roads,original_v249_values,reg)[
                'source_2D_CWM_beta_road_minus_interior']
                for reg,roads in sorted(regions.items()))
            if abs(m-PINNED_PRIOR_ROAD_RESULT)>1e-9:
                raise ValueError('Previously verified original MIREN beta source baseline drift')
        estimates={reg:region_null(roads,calibrated,draws)
                   for reg,roads in sorted(regions.items())}
        summary={}
        for mode in MODES:
            summary[mode]={}
            for zone in ('roadside','interior','road_minus_interior'):
                summary[mode][zone]={}
                for metric in METRICS:
                    vals={reg:est['source_null_by_mode'][mode]['effects'][zone][metric][
                        'observed_minus_trait_label_null'] for reg,est in estimates.items()}
                    loo=[mean(v for k,v in vals.items() if k!=skip) for skip in vals] if len(vals)>=3 else []
                    summary[mode][zone][metric]={
                        'source_region_equal_mean_observed_minus_trait_random_null':mean(vals.values()),
                        'source_region_count_positive_excess':sum(v>0 for v in vals.values()),
                        'source_region_n':len(vals),
                        'source_region_LOO_mean_minmax':[min(loo),max(loo)] if loo else None,
                        'original_region_codes_signed_values':vals,
                        'no_region_codes_are_random_country_or_publication_replicates':True
                    }
        result['panels'][panel]={
            'source_original_matched_study_support_QA':support,
            'by_region_only_source_aggregated_effects':estimates,
            'region_equal_descriptive_results':summary,
            'not_two_independent_papers_the_same_MIREN_source_for_both_panels':True,
        }
    result['provenance_warning']='Original regional source plant composition, climate, phylogeny and LEDA missingness are NOT controlled by trait-label shuffling'
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--miren',type=Path,required=True)
    p.add_argument('--height',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--draws',type=int,default=DRAWS_PER_CHAIN)
    a=p.parse_args()
    d=run(a.miren,a.height,a.seed,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({
      'status':'SOURCE_AUTHENTICATED_MIREN_ROADS_ONLY',
      'n_panels':len(d['panels']),
      'summary':{panel:{mode:{
         zone:{met:{
           'mean':stats['source_region_equal_mean_observed_minus_trait_random_null'],
           'positive_regions':stats['source_region_count_positive_excess'],
           'n_regions':stats['source_region_n']}
           for met,stats in metrics.items()}
         for zone,metrics in zones.items()}
         for mode,zones in p['region_equal_descriptive_results'].items()}
         for panel,p in d['panels'].items()}
    },ensure_ascii=False,indent=2))

if __name__=='__main__':main()

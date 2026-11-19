#!/usr/bin/env python3
"""v29.3 four REAL independent railway studies, genus-conditioned trait null.

Original author site×taxon matrix and original geography fully FROZEN, paired
LEDA height + seed mass values permuted only between species of same GENUS.
This is a phylogenetic-proxy falsification of v28.5 positive height sorting,
not evidence of within-species adaptation, ballast selection or train dispersal.

The v28.5 published original-source result is re-derived at the original
draw count before any genus-only null and anchored against frozen outcomes.
HR+BA is ONE original field publication, not two. Slovenia's 89 observations
are 100-m ballast locations, not railway stations; other 3 survey station-wide.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean
import numpy as np

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import (
    STUDIES, N, METRICS, MODES, LEDASHA,
    publisher_czech, publisher_italy, publisher_hr_ba,
    extract_original_xlsx, frozen_reference, build_matrix, metric_values,
    source_study_null,
)
from validate_MIREN_genus_conditioned_trait_label_null_v292 import genus

DRAWS_PER_CHAIN=299
GENUS_SEEDS=(293119,293120)
# Exact source original two-seed ×299 technical v28.5 occupancy-band null
# reported to six decimals in prior v28.5 original study report. 5e-6
# tolerance is rounding of display, NOT post hoc new biological selection.
FROZEN_V285_QUINTILE_HEIGHT={
 'CZ_Czechia':.054819,'IT_Italy_Parma':.038810,
 'HR_BA_2017_one_study':.099189,'SI_Slovenia_100m_track':.078781,
}
FROZEN_V285_QUINTILE_JOINT={
 'CZ_Czechia':.058534,'IT_Italy_Parma':.048722,
 'HR_BA_2017_one_study':.070119,'SI_Slovenia_100m_track':.080191,
}
MIN_SOURCE_MULTISPECIES_GENERA=10
MIN_SOURCE_TRAIT_KNOWN_SITE_INCIDENCE_MOVABLE=.25
MIN_ORIGINAL_GEOGRAPHY_GROUPS_WITH_MOVABLE_SOURCE=2
MIN_STUDIES_FOR_DESCRIPTIVE_SUMMARY=3
METRIC_KEYS=('cwm2D_mean_within_original_geography',
             'cwm_height_absolute_mean_within_original_geography',
             'cwm_seed_mass_absolute_mean_within_original_geography')

def original_species_and_genus_index(sites, traits)->tuple[list[str],dict[str,list[int]]]:
    taxa=sorted(set().union(*(site['taxa'] for site in sites)))
    if not taxa or any(t not in traits for t in taxa):
        raise ValueError('Original source study trait-known species were lost')
    groups=defaultdict(list)
    for i,t in enumerate(taxa):
        groups[genus(t)].append(i)
    return taxa,dict(sorted(groups.items()))

def genus_source_support(sites,traits)->dict:
    species,groups=original_species_and_genus_index(sites,traits)
    move={species[i] for idx in groups.values() if len(idx)>=2 for i in idx}
    total=sum(len(site['taxa']) for site in sites)
    movable=sum(len(site['taxa']&move) for site in sites)
    if not total:raise ValueError('No original site taxa in field source')
    by_geo=defaultdict(lambda:[0,0])
    for site in sites:
        n=len(site['taxa'])
        if n<=0:raise ValueError('Empty original site')
        by_geo[site['group']][0]+=len(site['taxa']&move)
        by_geo[site['group']][1]+=n
    geo_fractions={g:mob/count for g,(mob,count) in by_geo.items()}
    eligible=(sum(len(x)>=2 for x in groups.values())>=MIN_SOURCE_MULTISPECIES_GENERA and
              movable/total>=MIN_SOURCE_TRAIT_KNOWN_SITE_INCIDENCE_MOVABLE and
              sum(f>=MIN_SOURCE_TRAIT_KNOWN_SITE_INCIDENCE_MOVABLE
                  for f in geo_fractions.values())>=MIN_ORIGINAL_GEOGRAPHY_GROUPS_WITH_MOVABLE_SOURCE)
    return {
        'n_original_true_field_sites':len(sites),
        'n_original_joint_trait_taxa':len(species),
        'n_original_genus':len(groups),
        'n_multispecies_trait_known_genus':sum(len(x)>=2 for x in groups.values()),
        'n_original_trait_known_species_in_exchangeable_genus':len(move),
        'n_original_singleton_genus_species_immobile':len(species)-len(move),
        'original_source_trait_known_plot_taxon_incidence_total':total,
        'source_site_incidence_in_exchangeable_genus':movable,
        'source_incidence_exchangeable_fraction':movable/total,
        'original_geography_group_count':len(by_geo),
        'original_geography_group_exchangeable_incidence_fraction':geo_fractions,
        'original_geography_groups_above_25pct':sum(
            f>=MIN_SOURCE_TRAIT_KNOWN_SITE_INCIDENCE_MOVABLE for f in geo_fractions.values()),
        'min_multispecies_genus':MIN_SOURCE_MULTISPECIES_GENERA,
        'min_exchangeable_incidence_fraction':MIN_SOURCE_TRAIT_KNOWN_SITE_INCIDENCE_MOVABLE,
        'min_supported_source_geography_groups':MIN_ORIGINAL_GEOGRAPHY_GROUPS_WITH_MOVABLE_SOURCE,
        'source_genus_trait_shuffle_eligible':eligible,
        'status':'ELIGIBLE_SOURCE_GENUS_TRAIT_NULL' if eligible
                 else 'NOT_IDENTIFIABLE_GENUS_STRATIFIED_TRAIT_SHUFFLE',
    }

def genus_numpy_labels(T:np.ndarray,groups:dict[str,list[int]],rng:np.random.Generator)->np.ndarray:
    if T.ndim!=2 or T.shape[1]!=2:raise ValueError('Need true paired continuous traits')
    if sum(len(x) for x in groups.values())!=len(T):
        raise ValueError('Original genus binning does not partition source trait taxa')
    out=T.copy()
    for idx in groups.values():
        if len(idx)<2:continue
        p=rng.permutation(idx)
        out[idx,:]=T[p,:]
        if not np.array_equal(np.sort(out[idx,:],axis=0),np.sort(T[idx,:],axis=0)):
            raise ValueError('Height or seed mass escaped original genus stratum')
    if not np.array_equal(np.sort(out,axis=0),np.sort(T,axis=0)):
        raise ValueError('Genuine original height or seed mass changed by placebo')
    return out

def source_genus_null(sites,traits,study:str,draws:int=DRAWS_PER_CHAIN)->dict:
    if draws<5:raise ValueError('At least five technical relabels per source chain')
    support=genus_source_support(sites,traits)
    if not support['source_genus_trait_shuffle_eligible']:
        return {'support':support,'source_effects':None,
                'status':'NOT_IDENTIFIABLE_GENUS_STRATIFIED_TRAIT_SHUFFLE'}
    W,T,geo,qa=build_matrix(sites,traits)
    names,groups=original_species_and_genus_index(sites,traits)
    if len(names)!=T.shape[0]:raise ValueError('Original W columns and taxa changed')
    ii,jj,within=(qa['i'],qa['j'],qa['within'])
    actual=metric_values(W,T,ii,jj,within)
    samples={key:[] for key in METRIC_KEYS}
    for seed in GENUS_SEEDS:
        rng=np.random.default_rng(seed+17911*STUDIES.index(study))
        for _ in range(draws):
            sh=genus_numpy_labels(T,groups,rng)
            metrics=metric_values(W,sh,ii,jj,within)
            for key in samples:samples[key].append(metrics[key])
    return {'support':support,
            'n_technical_label_draws':2*draws,
            'field_source_CWM_original_within_geography':{
                key:actual[key] for key in METRIC_KEYS},
            'source_effects':{
                key:{
                    'original_observed':actual[key],
                    'within_genus_label_null_mean':mean(vals),
                    'observed_minus_genus_label_null':actual[key]-mean(vals),
                    'technical_upper_rank_not_field_p':
                        (1+sum(x>=actual[key] for x in vals))/(len(vals)+1)}
                for key,vals in samples.items()},
            'source_taxon_site_matrix_is_frozen':True,
            'source_bivariate_genus_traits_preserved':True,
            'status':'SOURCE_GENUS_TRAIT_NULL_ESTIMABLE'}

def run(cz_s1:Path,cz_s2:Path,it_pdf:Path,hr_pdf:Path,si_xlsx:Path,
        si_audit:Path,height:Path,seed:Path,draws:int=DRAWS_PER_CHAIN)->dict:
    for p,k in ((height,'height'),(seed,'seed')):
        if hashlib.sha256(p.read_bytes()).hexdigest()!=LEDASHA[k]:
            raise ValueError('Oldenburg LEDA source changed')
    h,_=read_leda(height,'height_m')
    s,_=read_leda(seed,'seed_mass_mg')
    published={'CZ_Czechia':publisher_czech(cz_s1,cz_s2),
               'IT_Italy_Parma':publisher_italy(it_pdf),
               'HR_BA_2017_one_study':publisher_hr_ba(hr_pdf)}
    studies={k:v[0] for k,v in published.items()}
    origqa={k:v[1] for k,v in published.items()}
    original_sl=extract_original_xlsx(si_xlsx,si_audit)
    cohorts,trait,baseline=frozen_reference(
        (studies,{**origqa,'slovenia':original_sl}),{'height':h,'seed':s})
    if set(cohorts)!=set(STUDIES) or sum(len(cohorts[k]) for k in STUDIES)!=171:
        raise ValueError('4 independent study field units drifted')
    result={}
    for key in STUDIES:
        if len(cohorts[key])!=N[key]:
            raise ValueError('Original source author field sample count drifted')
        previous=source_study_null(cohorts[key],trait[key],key,draws=299)
        old=previous['source_null_results'][MODES[1]]['metrics']
        hres=old['cwm_height_absolute_mean_within_original_geography'][
            'original_observed_minus_trait_randomized_null']
        jres=old['cwm2D_mean_within_original_geography'][
            'original_observed_minus_trait_randomized_null']
        if (abs(hres-FROZEN_V285_QUINTILE_HEIGHT[key])>5e-6 or
            abs(jres-FROZEN_V285_QUINTILE_JOINT[key])>5e-6):
            raise ValueError('Prior source v28.5 original cohort + occupancy null drift: '+key)
        newer=source_genus_null(cohorts[key],trait[key],key,draws=draws)
        result[key]={'original_field_study_count':1,
                     'source_original_cohort_qa':previous['source_QA'],
                     'v285_exact_original_species_and_geography_occupancy_quintile_null':{
                         m:old[m]['original_observed_minus_trait_randomized_null']
                         for m in METRIC_KEYS},
                     'genus_conditioned_source':newer}
    eligible={k:v['genus_conditioned_source'] for k,v in result.items()
              if v['genus_conditioned_source']['status']=='SOURCE_GENUS_TRAIT_NULL_ESTIMABLE'}
    summary=None
    if len(eligible)>=MIN_STUDIES_FOR_DESCRIPTIVE_SUMMARY:
        summary={}
        for metric in METRIC_KEYS:
            delta={k:v['source_effects'][metric]['observed_minus_genus_label_null']
                   for k,v in eligible.items()}
            loo=[mean(x for k,x in delta.items() if k!=name) for name in delta]
            summary[metric]={
                'publication_equal_mean_genus_conditioned_excess':mean(delta.values()),
                'n_positive_original_field_publications':sum(x>0 for x in delta.values()),
                'original_independent_field_publications':len(delta),
                'study_leave_one_out_mean_min_max':[min(loo),max(loo)],
                'study_results_not_station_level_pseudoreplications':delta}
    return {'version':'v29.3',
        'source_FOUR_independent_original_EU_railway_publications':True,
        'original_study_names':list(STUDIES),
        'original_3study_european_trait_axis_and_slovenia_holdout_QA':baseline,
        'exact_v285_original_2x299_occupancy_label_null_recreated':True,
        'sources':result,
        'source_genus_matched_publication_equal_summary_if_supported':summary,
        'no_global_phylogenetic_fraction_explained':True,
        'genus_not_full_molecular_phylogeny':True,
        'original_paired_LEDA_source_species_trait_matched_only':True,
        'not_direct_evidence_of_rail_specific_mowing_wind_or_train_dispersal':True,
        'does_not_supply_physically_matched_nonrail_original_field_controls':True}
def main():
    p=argparse.ArgumentParser()
    for arg in ('cz-s1','cz-s2','it-pdf','hr-pdf','si-xlsx','si-audit',
                'height','seed','output'):
        p.add_argument('--'+arg,required=True,type=Path)
    p.add_argument('--draws',type=int,default=DRAWS_PER_CHAIN)
    a=p.parse_args()
    r=run(a.cz_s1,a.cz_s2,a.it_pdf,a.hr_pdf,a.si_xlsx,a.si_audit,
          a.height,a.seed,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps({'version':'v29.3',
        'source_study_statuses':{
            k:{'original_field_sites':v['original_field_study_count'],
               'original_taxa':v['source_original_cohort_qa']['source_actual_plant_species_with_both_traits'],
               'source_genus_support':v['genus_conditioned_source']['support'],
               'source_genus_effect':v['genus_conditioned_source']['source_effects'],
               'original_v285_occupancy_null_excess':v[
                   'v285_exact_original_species_and_geography_occupancy_quintile_null']}
            for k,v in r['sources'].items()},
        'source_publication_equal_summary_if_supported':
              r['source_genus_matched_publication_equal_summary_if_supported']},
        indent=2,ensure_ascii=False))
if __name__=='__main__':main()

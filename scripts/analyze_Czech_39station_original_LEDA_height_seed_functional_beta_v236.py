#!/usr/bin/env python3
"""v23.6: original Czech 39-station 2020–22 plant flora continuous functional beta.

Source: Kutlvašr et al. 2024 Preslia 96:247-266 genuine S1 recent "rec"
survey columns 39 stations, 761 taxa, 6941 real species×station incidences.
Source: publisher's separate corrected erratum S2 verified 5 original
author railway stretches. Source: original University Oldenburg LEDA
species-median canopy height [m] and seed mass [mg]. 525 source taxa
have BOTH matched without fuzzy/genus guessing, per v22.0 SOURCE PASS.

Response:
 station-specific mean z(log height), mean z(log seed), trace covariance,
 and across station pair continuous functional centroid distance and
 exact matched-source taxonomic Jaccard. A single broad rail-corridor
 source 39 stations, not 39 independent studies and NOT the Slovenian
 matched same-location ballast shoulder/toe contrast.

Null:
 Explicit two-chain checkerboard swaps ONLY within each of the FIVE
 original author railway stretches, preserving every source station's
 exact joint-trait-known taxon richness and every source taxon's source
 station occupancy frequency WITHIN stretch. Taxon identity/median
 measured trait values NEVER assigned to another species. This tests
 within-study trait arrangement conditional on actual regional station
 source species pools, not convergence against independent natural
 regional nonrail flora. No active freight/wind/mowing/genotype measured.

Run source-authenticated publisher + LEDA data on GitHub Actions,
archive only aggregate estimates and fingerprints; no raw S1/S2/LEDA.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from analyze_Czech_2024_39station_functional_identity_v143 import read_original_czech_publisher
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import ztransform,components,quart

LEDA_HASH={
  'height_m':'cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92',
  'seed_mass_mg':'34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef',
}
DRAW_PER_CHAIN=120
BURNIN=30000
GAP=3000
SEED=20261116+236
N_SOURCE_STATIONS=39
N_SOURCE_STRETCHES=5
N_SOURCE_DUAL_TRAITS=525

def source_cohort(s1,s2,h,seed):
    floras,meta,source=read_original_czech_publisher(s1,s2)
    recent=set().union(*floras.values())
    if len(recent)!=761 or len(floras)!=39 or sum(len(s) for s in floras.values())!=6941:
        raise ValueError('Original genuine 2024 39station 761taxa 6941incidence drift')
    joined,zqc=ztransform({x:binomial(x) for x in recent},h,seed)
    if len(joined)!=N_SOURCE_DUAL_TRAITS:
        raise ValueError('v22.0 original Czech 525 exact height+mass taxa mismatch')
    src=[]
    for station in sorted(floras):
        codes=floras[station]&set(joined)
        if len(codes)<10:raise ValueError('Original Czech station <10 dual trait plants')
        info=meta[station]
        src.append({'station':station,'stretch':info['stretch'],
                    'taxa':set(codes),
                    'source_unequal_station_area_ha':info['original_source_station_total_area_ha']})
    groups=defaultdict(list)
    for i,x in enumerate(src):groups[x['stretch']].append(i)
    if len(groups)!=N_SOURCE_STRETCHES or min(len(q) for q in groups.values())<3:
        raise ValueError('Original 5 railway stretch source sampling support insufficient')
    return src,joined,{
      'author_Preslia_source':source,'source_exact_trait_zlog_species_scaling':zqc,
      'n_original_stations':len(src),
      'n_author_railway_stretches':len(groups),
      'stations_per_publisher_stretch':{s:len(v) for s,v in groups.items()},
      'n_exact_joint_LEDA_taxa':len(joined),
      'original_station_taxon_incidence_total':6941,
      'joint_trait_station_incidence':sum(len(x['taxa']) for x in src),
      'source_station_joint_trait_richness_min_med_max':[
        min(len(x['taxa']) for x in src),median(len(x['taxa']) for x in src),
        max(len(x['taxa']) for x in src)],
      'source_area_known_station_n':sum(x['source_unequal_station_area_ha'] is not None for x in src),
      'not_equal_original_station_area_or_true_rail_BS_BT_areas':True}

def original_margins(sites):
    group=defaultdict(list)
    for s in sites:group[s['stretch']].append(s['taxa'])
    return {
      'site_by_source_order_known_trait_richness':[len(s['taxa']) for s in sites],
      'original_each_stretch_species_occupancy':{
         g:dict(Counter(t for flora in fs for t in flora)) for g,fs in group.items()}}

def stats(sites,traits):
    profiles=[]
    for s in sites:
        d=components(s['taxa'],traits)
        profiles.append(d)
    within={'tax_Jaccard':[],'functional_centroid_distance':[],
            'source_trait_width_trace':[]}
    between={'tax_Jaccard':[],'functional_centroid_distance':[]}
    allpairs={'tax_Jaccard':[],'functional_centroid_distance':[]}
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i]['taxa'],sites[j]['taxa']
        jac=len(a&b)/len(a|b)
        distance=math.hypot(
            profiles[i]['height_mean']-profiles[j]['height_mean'],
            profiles[i]['seed_mean']-profiles[j]['seed_mean'])
        typ=within if sites[i]['stretch']==sites[j]['stretch'] else between
        typ['tax_Jaccard'].append(jac)
        typ['functional_centroid_distance'].append(distance)
        allpairs['tax_Jaccard'].append(jac)
        allpairs['functional_centroid_distance'].append(distance)
    if len(allpairs['tax_Jaccard'])!=741:
        raise ValueError('Czech 39 unique stations must have 741 *dependent* pairs')
    trmean=mean(q['joint_trace'] for q in profiles)
    return {
      'mean_station_joint_trace':trmean,
      'mean_station_height_variance':mean(q['height_var'] for q in profiles),
      'mean_station_seed_mass_variance':mean(q['seed_var'] for q in profiles),
      'all_station_pair_taxon_Jaccard':mean(allpairs['tax_Jaccard']),
      'all_station_pair_centroid_distance':mean(allpairs['functional_centroid_distance']),
      'within_original_stretch_pair_taxon_Jaccard':mean(within['tax_Jaccard']),
      'between_original_stretch_pair_taxon_Jaccard':mean(between['tax_Jaccard']),
      'within_original_stretch_pair_centroid_distance':mean(within['functional_centroid_distance']),
      'between_original_stretch_pair_centroid_distance':mean(between['functional_centroid_distance']),
      'within_minus_between_Jaccard':mean(within['tax_Jaccard'])-mean(between['tax_Jaccard']),
      'within_minus_between_centroid_distance':mean(within['functional_centroid_distance'])-mean(between['functional_centroid_distance']),
    }

def run_chain(sites,traits,seed,draws,burnin,gap):
    rng=random.Random(seed)
    group=defaultdict(list)
    for i,s in enumerate(sites):group[s['stretch']].append(i)
    labels=sorted(group)
    current=[set(s['taxa']) for s in sites]
    original=original_margins(sites)
    acc=0
    def advance(count):
        nonlocal acc
        for _ in range(count):
            g=rng.choice(labels)
            i,j=rng.sample(group[g],2)
            a=current[i]-current[j];b=current[j]-current[i]
            if not a or not b:continue
            # Strict one-for-one source checkerboard shuffle in same rail stretch.
            x=rng.choice(sorted(a));y=rng.choice(sorted(b))
            current[i].remove(x);current[i].add(y)
            current[j].remove(y);current[j].add(x)
            acc+=1
    advance(burnin)
    vals=[]
    for z in range(draws):
        advance(gap)
        pert=[{**src,'taxa':current[i]} for i,src in enumerate(sites)]
        if original_margins(pert)!=original:
            raise ValueError('Original 39-station richness or source stretch-specific taxon occurrence drifted under null')
        vals.append(stats(pert,traits))
    if acc<1000:raise ValueError('Original within-stretch checkerboard swaps did not mix adequately')
    return vals,{'source_original_rail_stretch_constrained_accepted_swaps':acc,
                 'proposals':burnin+draws*gap,'draws':draws}

def analyze(s1,s2,h,seed,draws=DRAW_PER_CHAIN,burnin=BURNIN,gap=GAP):
    sites,values,qa=source_cohort(s1,s2,h,seed)
    real=stats(sites,values)
    metrics=tuple(real.keys())
    nulls=[];bychain=[];chainqc=[]
    for chain in range(2):
        arr,mix=run_chain(sites,values,SEED+7193*chain,draws,burnin,gap)
        nulls.extend(arr);bychain.append(arr);chainqc.append(mix)
    result={}
    for k in metrics:
        observed=real[k];ns=[r[k] for r in nulls]
        avg=mean(ns)
        result[k]={
          'observed_original_author_stations':observed,
          'source_stretch_fixed_margin_null_mean':avg,
          'observed_minus_null':observed-avg,
          'technical_null_quantile_025':quart(ns,.025),
          'technical_null_quantile_975':quart(ns,.975),
          'technical_lower_tail_rank_plusone_NOT_field_p':
             (1+sum(v<=observed for v in ns))/(1+len(ns)),
          'difference_in_two_null_chain_means':
             abs(mean(x[k] for x in bychain[0])-mean(x[k] for x in bychain[1])),
          'direction_of_convergence_for_metric':(
               'negative_is_more_functionally_similar' if 'centroid_distance' in k else
               'positive_is_more_taxonomically_similar' if 'Jaccard' in k else
               'negative_is_narrower_station_functional_breadth')}
    return {'version':'v23.6',
       'status':'SOURCE_VERIFIED_CZECH_39RAIL_STATION_CONTINUOUS_TRAIT_SOURCE_INTERNAL_NULL',
       'source_genuine_Preslia_39stations_5stretches_2020_22':qa,
       'n_source_dependent_station_pairs':741,
       'n_source_null_chains':len(bychain),
       'draws_per_chain':draws,
       'source_original_train_GBIF_freight_count_not_observed':True,
       'null_chain_QA':chainqc,
       'source_observed_and_within_stretch_trait_preserving_margins_comparisons':result,
       'study_source_limits':[
          'Source 39 Czech station presences pool 3 different habitats, NOT 39×3 zone resolved source or same 100m railway transect scale as Slovenian BS-BT',
          'LEDA species median height+seed mass is not local phenotypic size or railway adaptation',
          'Fixed-margin within-original-railway-stretch null holds the already-railway-selected regional floristic pool; no independent natural nonrail plot control',
          'Original 5 stretches are not 5 independent countries or independent ecosystem natural experiments',
          'Source station areas and sampling efforts unequal; 39 observed stations have no precise GPS in S1/S2 for matching climate and land use',
          'Station pair counts are dependent, 240 technical null draws are not independent field experiments or causal p-values',
          'No original actual goods/passenger train counts, mowing dates, wind speed, freight OD or seed intercept observations',
          'No phylogeny and no genetic/phenotypic railway-selected alleles measured',
          'Do not equate within-station functional variance with between-station functional centroid beta, or pool with Slovenia as identical study design']}
def main():
    p=argparse.ArgumentParser()
    for arg in ('cz_s1','cz_s2','height','seed','output'):
        p.add_argument('--'+arg.replace('_','-'),required=True)
    p.add_argument('--draws',type=int,default=DRAW_PER_CHAIN)
    a=p.parse_args()
    for file,digest in ((a.height,LEDA_HASH['height_m']),(a.seed,LEDA_HASH['seed_mass_mg'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=digest:
            raise ValueError('Not the exact original source-verified LEDA height / seed-mass bytes')
    h,qh=read_leda(a.height,'height_m');se,qs=read_leda(a.seed,'seed_mass_mg')
    report=analyze(a.cz_s1,a.cz_s2,h,se,draws=a.draws)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    summary={k:report[k] for k in ('status','source_genuine_Preslia_39stations_5stretches_2020_22',
       'n_source_dependent_station_pairs','null_chain_QA',
       'source_observed_and_within_stretch_trait_preserving_margins_comparisons')}
    print('V236_REAL_CZECH_HEIGHT_SEED_FUNCTIONAL_BETA='+json.dumps(summary,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

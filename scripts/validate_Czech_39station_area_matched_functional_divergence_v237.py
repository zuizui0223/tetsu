#!/usr/bin/env python3
"""v23.7: 39 authentic Czech rail stations functional beta station-area guard.

Follow v23.6 same genuine Preslia 2024 publisher S1 39 station rec taxa +
corrected erratum S2 five original stretches, original 525 exact LEDA
species-median CANOPY height and seed mass. NO synthetic field data.

Sensitivity:
 Original source station survey area (all 39 known) is a strong potential
 confounder for functional centroid differences. Form within each of the
 five original railway stretches two station-area-rank strata (small vs
 large), with source n >=2 per stratum. Preserve EACH original station's
 exact two-trait-known richness and EACH original species' occupancy
 FREQUENCY WITHIN source stretch×area stratum. Shuffle only checkerboard
 occurrence pairs within that stratum, never move a taxon across
 geographic stretch or station-area bin. This is strictly more
 conditioning than v23.6 (within-stretch only). No natural nonrail
 reference is created; area is an exposure/survey opportunity proxy.

 Repeated 2 independent random walks, 100 draws each. Every draw
 computes all 39 site centroid beta and five leave-one-source-stretch
 aggregated station subsets. Do not treat 741 dependent station pairs or
 200 autocorrelated null draws as independent field replicates.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from analyze_Czech_39station_original_LEDA_height_seed_functional_beta_v236 import (
    source_cohort,stats,LEDA_HASH
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import components,quart
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import spearman_ties_pairs

SEED=20261116+237
DRAWS=100
BURNIN=25000
GAP=3000

def area_bands(sites):
    s=defaultdict(list)
    for i,r in enumerate(sites):
        a=r['source_unequal_station_area_ha']
        if a is None or not isinstance(a,(float,int)) or a<=0 or not math.isfinite(a):
            raise ValueError('Publisher area data missing, nonpositive, or not numeric; no area adjustment')
        s[r['stretch']].append((float(a),r['station'],i))
    grouped=defaultdict(list);qa={}
    for stretch,rows in sorted(s.items()):
        order=sorted(rows)
        n=len(order);split=n//2
        if min(split,n-split)<2:
            raise ValueError('Within source stretch not enough low/high original survey area stations')
        if order[split-1][0]==order[split][0]:
            raise ValueError('Original source area bin boundary splits identical hectares')
        for area,name,index in order[:split]:grouped[(stretch,'small_area')].append(index)
        for area,name,index in order[split:]:grouped[(stretch,'large_area')].append(index)
        qa[stretch]={'n_original_station_areas':n,'n_smaller':split,'n_larger':n-split,
             'original_area_ha_min':order[0][0],
             'original_area_ha_lower_bin_max':order[split-1][0],
             'original_area_ha_upper_bin_min':order[split][0],
             'original_area_ha_max':order[-1][0]}
    if len(grouped)!=10 or len(s)!=5:raise ValueError('Author original 5 stretches×2 station survey-area bins changed')
    return grouped,qa

def source_margins(sites,groups):
    by={}
    size=[len(s['taxa']) for s in sites]
    for key,idxs in groups.items():
        by[str(key)]=dict(Counter(t for i in idxs for t in sites[i]['taxa']))
    return {'station_39_richness_unchanged':size,'stretch_by_area_source_species_incidence':by}

def means_and_LOO(sites,values):
    center={}
    for i,r in enumerate(sites):
        d=components(r['taxa'],values)
        center[i]=(d['height_mean'],d['seed_mean'],d['joint_trace'])
    dists=[];within=[];jacs=[];area=[(math.log(r['source_unequal_station_area_ha']),
                                       *center[i]) for i,r in enumerate(sites)]
    influence=defaultdict(list);stretch_ids=sorted(set(r['stretch'] for r in sites))
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        v=math.hypot(center[i][0]-center[j][0],center[i][1]-center[j][1])
        dists.append(v)
        jacs.append(len(a['taxa']&b['taxa'])/len(a['taxa']|b['taxa']))
        if a['stretch']==b['stretch']:within.append(v)
        for excluded in stretch_ids:
            if a['stretch']!=excluded and b['stretch']!=excluded:
                influence[excluded].append(v)
    if len(dists)!=741 or len(stretch_ids)!=5 or not all(len(x)>100 for x in influence.values()):
        raise ValueError('Source 39 stations 5 original stretches and 741 dependent pairs drift')
    return {
       'all_39station_pair_functional_centroid_distance':mean(dists),
       'within_original_stretch_functional_centroid_distance':mean(within),
       'all_39station_pair_taxon_Jaccard':mean(jacs),
       'mean_39station_two_trait_variance':mean(x[2] for x in center.values()),
       'leave_one_source_stretch_out_beta':{k:mean(v) for k,v in sorted(influence.items())},
       'original_area_descriptive_source_correlations':{
          'log_survey_area_ha_vs_species_CWM_zlog_height':
               spearman_ties_pairs([(x[0],x[1]) for x in area]),
          'log_survey_area_ha_vs_species_CWM_zlog_seed_mass':
               spearman_ties_pairs([(x[0],x[2]) for x in area]),
          'log_survey_area_ha_vs_station_joint_trace':
               spearman_ties_pairs([(x[0],x[3]) for x in area])}}

def shuffle_null(sites,values,groups,draws=DRAWS,burnin=BURNIN,gap=GAP):
    origin=source_margins(sites,groups)
    samples=[];bychain=[];qcs=[]
    for ch in (0,1):
        rng=random.Random(SEED+ch*10007)
        current=[set(site['taxa']) for site in sites]
        labels=sorted(groups)
        accepts=0
        def step(n):
            nonlocal accepts
            for _ in range(n):
                key=rng.choice(labels)
                i,j=rng.sample(groups[key],2)
                a=current[i]-current[j];b=current[j]-current[i]
                if not a or not b:continue
                u=rng.choice(sorted(a));v=rng.choice(sorted(b))
                current[i].remove(u);current[i].add(v)
                current[j].remove(v);current[j].add(u)
                accepts+=1
        step(burnin)
        obs=[]
        for i in range(draws):
            step(gap)
            scenario=[{**s,'taxa':current[ix]} for ix,s in enumerate(sites)]
            if source_margins(scenario,groups)!=origin:
                raise ValueError('Original station richness or regional×area source species occupancy not conserved')
            v=means_and_LOO(scenario,values)
            obs.append(v)
            samples.append(v)
        if accepts<1000:raise ValueError('Original source within area bin checkerboard null failed to mix')
        bychain.append(obs)
        qcs.append({'technical_chain':ch+1,'source_swap_proposals':burnin+draws*gap,
              'original_source_swaps_accepted':accepts,'autocorrelated_diagnostic_draws':draws})
    return samples,bychain,qcs

def report(original,values,draws=DRAWS):
    groups,areaqc=area_bands(original)
    actual=means_and_LOO(original,values)
    prior=stats(original,values)
    expected=prior['all_station_pair_centroid_distance']
    if abs(expected-.16520612655692127)>1e-9 or abs(actual[
         'all_39station_pair_functional_centroid_distance']-expected)>1e-9:
        raise ValueError('v23.6 authentic Czech 39 station continuum effect changed')
    samples,chains,cqa=shuffle_null(original,values,groups,draws=draws)
    metrics={}
    for k in ('all_39station_pair_functional_centroid_distance',
              'within_original_stretch_functional_centroid_distance',
              'all_39station_pair_taxon_Jaccard',
              'mean_39station_two_trait_variance'):
        o=actual[k]
        nul=[q[k] for q in samples]
        metrics[k]={'true_original_observed':o,
          'source_stretch_AND_station_area_ranks_fixed_expected':mean(nul),
          'observed_minus_area_stratified_fixed_margin_null':o-mean(nul),
          'null_q025':quart(nul,.025),'null_q975':quart(nul,.975),
          'technical_lower_tail_rank':(1+sum(x<=o for x in nul))/(1+len(nul)),
          'two_chains_mean_difference':abs(mean(x[k] for x in chains[0])-
                                           mean(x[k] for x in chains[1]))}
    left={}
    for stretch,v in actual['leave_one_source_stretch_out_beta'].items():
        nul=[q['leave_one_source_stretch_out_beta'][stretch] for q in samples]
        left[stretch]={'source_observed_centroid_distance_after_dropping_stretch':v,
          'null_mean_after_dropping_stretch':mean(nul),
          'observed_minus_null':v-mean(nul),
          'technical_lower_tail_rank':(1+sum(x<=v for x in nul))/(1+len(nul))}
    return {'version':'v23.7',
      'status':'SOURCE_VERIFIED_CZECH_39RAIL_STATIONS_CONTINUOUS_TRAIT_DIVERGENCE_AREA_STRATIFIED',
      'true_source_39_stations_525_height_and_seed_median_taxa':True,
      'original_five_stretch_survey_area_strata':areaqc,
      'source_station_area_rank_strata_count':len(groups),
      'station_areas_all_positive_and_known':True,
      'source_2020_22_Czech_original_station_count':len(original),
      'fixed_margin_area_matched_2chain_technical_QA':cqa,
      'observed_and_original_area_stratified_null':metrics,
      'source_leave_one_stretch_out_centroid_discrepancy':left,
      'original_area_covariate_descriptive_rank_correlations':
          actual['original_area_descriptive_source_correlations'],
      'original_stretch_only_v236_centroid_observed_minus_null_reference':0.04457802597382583,
      'limitations':[
         'No actual natural nonrail plots; preserve observed within-railway species source pools only',
         'Station survey areas vary by station and are only coarsely conditioned via within-original-stretch median area rank',
         'Survey habitat zones pooled; this is not Slovenia same-side BS vs BT paired trial',
         'Source 5 original stretches from one paper, no independent regional-country replicates',
         'Species-median CANOPY height is top photosynthetic tissue, not total flowering stem height',
         'Higher source functional centroid difference does not imply genetic or ecological adaptation',
         'Site pairs, technical permutations and leave-one-stretch recalculations do not form independent field experiments',
         'No measured 2020–22 local train movements, freight goods, mowing, herbicide, wind or regionally matched natural species pool']}
def main():
    p=argparse.ArgumentParser()
    for key in ('cz_s1','cz_s2','height','seed','output'):
        p.add_argument('--'+key.replace('_','-'),required=True)
    p.add_argument('--draws',type=int,default=DRAWS)
    args=p.parse_args()
    for file,digest in ((args.height,LEDA_HASH['height_m']),
                        (args.seed,LEDA_HASH['seed_mass_mg'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=digest:
            raise ValueError('Original LEDA height and seed mass source SHA drift')
    h,hqa=read_leda(args.height,'height_m')
    se,sqa=read_leda(args.seed,'seed_mass_mg')
    src,values,provenance=source_cohort(args.cz_s1,args.cz_s2,h,se)
    r=report(src,values,draws=args.draws)
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf8')
    print('V237_REAL_CZECH_STATION_AREA_NULL='+json.dumps({
       'status':r['status'],
       'source_2020_22_Czech_original_station_count':r['source_2020_22_Czech_original_station_count'],
       'original_five_stretch_survey_area_strata':r['original_five_stretch_survey_area_strata'],
       'observed_and_original_area_stratified_null':r['observed_and_original_area_stratified_null'],
       'source_leave_one_stretch_out_centroid_discrepancy':r['source_leave_one_stretch_out_centroid_discrepancy'],
       'original_area_covariate_descriptive_rank_correlations':r['original_area_covariate_descriptive_rank_correlations']
       },ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

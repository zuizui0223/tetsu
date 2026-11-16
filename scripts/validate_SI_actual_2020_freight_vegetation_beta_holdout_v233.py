#!/usr/bin/env python3
"""v23.3: historical 2020 GOODS train intensity vs 2022/23 original flora beta,
geography / passenger frequency / same national line adjusted, held-out
physical segment validation.

Read true 2023 official railway station PDF SHA, 2026 author-released XLSX
(2022–23 field sites) SHA and real Eurostat 2020 directional G/P counts.
Use EXACT v23.1 unique 27 real field sites -> 9 true two-way physical
intervals on 5 national lines. Never treat a directional pair as two
independent railway stations. All site pair Jaccards are aggregated
FIRST to 36 *dependent* physical segment pair means.

Compare baseline y ~ log geodesic separation + log passenger-count
difference + same national source line + fraction same field-sampling year
against baseline + log freight-count difference, evaluating full vs
baseline in leave-one-physical-segment-out cross validation (9 folds).
A 9-segment CV prediction score is a feasibility diagnostic not a
generalizable causal or independent station-trial effect.

Within-original-national-line exact label reassignment of historical
freight counts (48 finite permutations for this cohort) reveals whether
the coefficient depends only on differences between the five source
lines. These counts are still not 2022/2023 service exposure or rail
transport origin-destination. No inferential p-values are reported.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    jaccard,functional_overlap
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
    extract_authoritative_stations,source_vegetation_sites,
    eurostat_2020_segments,segment_route_station_join,
    collapse_exact_reverse_direction_segments,assign_segment
)
from analyze_SI_2020_traffic_vs_2022_2023_flora_exploratory_v232 import (
    SOURCE_PDF_SHA,segment_traffic_sites,per_physical_segment,betas
)

MIN_GEO_KM=2.
METRICS=("taxonomic_Jaccard","functional_four_guild_overlap")
BASE_COLS=("log_sitepair_geographic_distance","passenger_log_traffic_difference",
           "same_original_national_rail_line","same_original_survey_year_pair_share")

def build_segment_dyads(sites,guild):
    import numpy as np
    grouped=defaultdict(list)
    seg_info={}
    for site in sites:
        key=site['source_segment']
        vals=seg_info.get(key)
        now=(site['section'],site['freight'],site['passenger'])
        if vals is not None and vals!=now:
            raise ValueError('One original true physical Eurostat segment has inconsistent source traffic')
        seg_info[key]=now
    keys=sorted(seg_info)
    if len(keys)!=9:raise ValueError('Original v23.1 mapped segment universe drifted')
    index={k:i for i,k in enumerate(keys)}
    for a,b in itertools.combinations(sites,2):
        if a['source_segment']==b['source_segment']:continue
        dist=haversine_km(a['lat'],a['lon'],b['lat'],b['lon'])
        if dist<MIN_GEO_KM:continue
        i,j=sorted((index[a['source_segment']],index[b['source_segment']]))
        grouped[i,j].append({
            'taxonomic_Jaccard':jaccard(a['source_species_all'],b['source_species_all']),
            'functional_four_guild_overlap':
                functional_overlap(a['known_original_lifespan_diaspore_species'],
                                   b['known_original_lifespan_diaspore_species'],guild),
            'log_sitepair_geographic_distance':math.log1p(dist),
            'same_original_survey_year_pair_share':float(a['year']==b['year'])})
    if len(grouped)!=36:
        raise ValueError('Expected 36 distinct source physical segment pairs in the genuine 27-site cohort')
    frame=[]
    for (i,j),rs in sorted(grouped.items()):
        line1,freight1,passenger1=seg_info[keys[i]]
        line2,freight2,passenger2=seg_info[keys[j]]
        v={'index_i':i,'index_j':j,'original_source_site_pair_n':len(rs),
           'same_original_national_rail_line':float(line1==line2),
           'passenger_log_traffic_difference':abs(math.log1p(passenger1)-math.log1p(passenger2)),
           'log_freight_traffic_difference':abs(math.log1p(freight1)-math.log1p(freight2))}
        for k in METRICS+('log_sitepair_geographic_distance','same_original_survey_year_pair_share'):
            v[k]=mean(q[k] for q in rs)
        frame.append(v)
    return frame,{i:{"national_line":seg_info[k][0],
                     "freight_2020":seg_info[k][1],
                     "passenger_2020":seg_info[k][2]}
                  for i,k in enumerate(keys)}

def design(rows,include_goods=True,goods_override=None):
    import numpy as np
    means={}
    # No scaling based on the holdout outcome; raw source predictors are
    # fixed original observed values. Stable goods scaling globally.
    x=[]
    for r in rows:
        vals=[1.]+[r[k] for k in BASE_COLS]
        if include_goods:
            if goods_override is None:
                val=r['log_freight_traffic_difference']
            else:
                val=abs(goods_override[r['index_i']]-goods_override[r['index_j']])
            vals.append(val)
        x.append(vals)
    return np.asarray(x,dtype=float)

def fit(X,y):
    import numpy as np
    if X.shape[0]<X.shape[1]+3 or np.linalg.matrix_rank(X)<X.shape[1]:
        return None
    return np.linalg.lstsq(X,y,rcond=None)[0]

def score(y,true):
    import numpy as np
    return float(np.mean((y-true)**2))

def leave_one_segment_cv(rows,outcome,goods_override=None):
    import numpy as np
    base=design(rows,False)
    full=design(rows,True,goods_override)
    y=np.asarray([r[outcome] for r in rows],dtype=float)
    r=[]
    for k in range(9):
        tr=np.asarray([i for i,row in enumerate(rows)
                       if k not in (row['index_i'],row['index_j'])],dtype=int)
        te=np.asarray([i for i,row in enumerate(rows)
                       if k in (row['index_i'],row['index_j'])],dtype=int)
        if len(te)!=8:raise ValueError('Each original physical segment should have eight neighbors')
        b=fit(base[tr],y[tr]);f=fit(full[tr],y[tr])
        if b is None or f is None:return None
        r.append({'segment_index':k,'n_source_training_traffic_pairs':len(tr),
                  'n_source_heldout_traffic_pairs':len(te),
                  'baseline_mean_squared_error':score(y[te],base[te]@b),
                  'freight_augmented_mean_squared_error':score(y[te],full[te]@f)})
    return {'source_physical_segment_holdout_folds':len(r),
            'source_test_pairs_per_fold':8,
            'baseline_mean_fold_mse':mean(x['baseline_mean_squared_error'] for x in r),
            'freight_augmented_mean_fold_mse':mean(x['freight_augmented_mean_squared_error'] for x in r),
            'gain_baseline_minus_freight_MSE':mean(x['baseline_mean_squared_error']-
                                                  x['freight_augmented_mean_squared_error'] for x in r),
            'n_folds_improved_by_historical_freight':sum(x['freight_augmented_mean_squared_error']<
                                                          x['baseline_mean_squared_error'] for x in r),
            'fold_source_pairs_NOT_independent_field_experiments':True}

def observed_adjustment(rows,outcome,goods_override=None):
    import numpy as np
    X0=design(rows,False);X1=design(rows,True,goods_override)
    Y=np.asarray([r[outcome] for r in rows],dtype=float)
    b0=fit(X0,Y);b1=fit(X1,Y)
    if b0 is None or b1 is None:return {'status':'UNIDENTIFIABLE_RANK_DEFICIENT'}
    mse0=score(Y,X0@b0);mse1=score(Y,X1@b1)
    return {'status':'DESCRIPTIVE_ONLY',
        'n_source_dependent_physical_segment_pairs':len(rows),
        'baseline_in_sample_MSE':mse0,'add_historical_freight_in_sample_MSE':mse1,
        'in_sample_MSE_gain_baseline_minus_freight':mse0-mse1,
        'historical_logfreight_difference_conditional_coefficient':float(b1[-1]),
        'baseline_coefficients_named':dict(zip(('intercept',)+BASE_COLS,map(float,b0))),
        'goodness_of_fit_is_not_a_causal_rail_displacement_effect':True}

def exact_within_line_labels(meta):
    original={k:math.log1p(v['freight_2020']) for k,v in meta.items()}
    group=defaultdict(list)
    for i,item in meta.items():group[item['national_line']].append(i)
    sets=[]
    for line in sorted(group):
        seq=sorted(group[line])
        sets.append((seq,list(itertools.permutations([original[j] for j in seq]))))
    total=math.prod(len(p) for ids,p in sets)
    if total>1000:raise ValueError('Unexpected >1000 source-label permutations')
    yield total,None # first item declares size for source-QA
    for choices in itertools.product(*(items for _,items in sets)):
        x={}
        for (indices,_),vals in zip(sets,choices):
            for i,v in zip(indices,vals):x[i]=v
        yield total,x

def run(original,pdfpath):
    stat,sqa=extract_authoritative_stations(pdfpath)
    if sqa['sha256_original_pdf_source']!=SOURCE_PDF_SHA:
        raise ValueError('Authentic 2023 operator station PDF SHA changed')
    sites,locqa=source_vegetation_sites(original)
    net,nqa=eurostat_2020_segments()
    resolved,eq=segment_route_station_join(stat,net)
    complete,physicalqa=collapse_exact_reverse_direction_segments(resolved)
    mapped,crossqa=assign_segment(sites,complete)
    real,guild=segment_traffic_sites(original,mapped)
    pairs,seginfo=build_segment_dyads(real,guild)
    results={}
    for outcome in METRICS:
        observed=observed_adjustment(pairs,outcome)
        cv=leave_one_segment_cv(pairs,outcome)
        if cv is None:raise ValueError('v23.3 original controls rank deficient under leave-one-segment-out')
        observed['heldout_actual_source_traffic_segment_prediction']=cv
        n=0;betas=[];gains=[]
        for num,swap in exact_within_line_labels(seginfo):
            if swap is None:
                n=num
                continue
            test=observed_adjustment(pairs,outcome,swap)
            cvp=leave_one_segment_cv(pairs,outcome,swap)
            if test['status']=='DESCRIPTIVE_ONLY' and cvp:
                betas.append(test['historical_logfreight_difference_conditional_coefficient'])
                gains.append(cvp['gain_baseline_minus_freight_MSE'])
        if not betas or len(betas)!=n:
            raise ValueError('Historical freight source within-line comparison not fully evaluated')
        target=observed['historical_logfreight_difference_conditional_coefficient']
        actualgain=cv['gain_baseline_minus_freight_MSE']
        observed['source_exact_within_original_national_line_relabelled_freight_QA']={
           'all_unique_source_linewise_arrangements_including_original':n,
           'actual_2020_observed_freight_coefficient':target,
           'technical_line_constrained_relative_signed_coefficient_rank':
                sum(x<=target for x in betas)/n,
           'technical_permutation_absolute_coefficient_rank':
                sum(abs(x)<=abs(target) for x in betas)/n,
           'actual_LOSO_MSE_gain_baseline_minus_freight':actualgain,
           'technical_permutations_with_LOSO_gain_at_least_actual':
                sum(x>=actualgain for x in gains),
           'min_median_max_LOSO_MSE_gain_under_technical_linewise_reassignment':
                [min(gains),sorted(gains)[len(gains)//2],max(gains)],
           'NOT_a_randomized_field_design_or_formal_p_value':True}
        results[outcome]=observed
    return {
       'version':'v23.3','status':'ORIGINAL_SOURCE_GEO_AND_PASSENGER_ADJUSTED_FREIGHT_BETA_HOLDOUT',
       'n_original_2022_2023_flora_sites_with_unique_2020_segments':len(mapped),
       'n_unique_physical_2020_traffic_segments':len(seginfo),
       'n_national_source_railway_lines':len({v['national_line'] for v in seginfo.values()}),
       'n_original_dependent_physical_segment_pairs':len(pairs),
       'real_2023_station_PDF_SHA256':SOURCE_PDF_SHA,
       'real_2020_EUROSTAT_original_json_SHA256':nqa['source']['sha256_original_json_response'],
       'source_outcomes':results,
       'controls':list(BASE_COLS),
       'original_unbiased_traffic_segment_coverage':crossqa,
       'major_limitations':[
         'historical 2020 pandemic-affected train traffic is NOT 2022/23 field-year concurrent exposure',
         'only 9 real freight intervals, nested in 5 original national lines; dyads are dependent',
         'held-out folds 9 are not 9 independent countries, but original segment opportunity groups',
         'author 2022/23 original sites were nonrandom vegetated railway plots, NOT station-wide standardized flora',
         'Eurostat physical traffic segment codes carry no freight origin–destination; freight similarity is not shared direct train connection',
         'within-line permutations do not establish ecological exchangeability or independent observational p-values',
         'road/urbanization/botanical source regions and rail maintenance are not causally identified',
         'cannot infer railway seed dispersal, genetic adaptation, or rail-overrides-biogeography']}
def main():
    p=argparse.ArgumentParser()
    for field in ('original_xlsx','source_audit','official_station_pdf','output'):
        p.add_argument('--'+field.replace('_','-'),required=True)
    arg=p.parse_args()
    original=extract_original_xlsx(arg.original_xlsx,arg.source_audit)
    ans=run(original,arg.official_station_pdf)
    dst=Path(arg.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(ans,ensure_ascii=False,indent=2),encoding='utf8')
    print('V233_ADJUSTED_REAL_TRAIN_SOURCE='+json.dumps(ans,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""v23.5 authentic 2015/2020 freight movement historical sensitivity.

Paired historical Eurostat observed rail network direction datasets:
  rail_tf_ns15_si (2015) and rail_tf_ns20_si (2020).
Only strictly SAME original Eurostat source segment code AND endpoint label
are permitted for both years (never map 2015 by order, station geography
alone, or assume route counts were stable after 2020).

Independent source A: 2023 Slovenia operator network statement Annex2B
authenticated 16-page station km PDF SHA pinned.
Source B: genuine Šabić et al 2026 author source, 2022/23 railway vegetation.
Source C: Eurostat true historical direction-specific goods/passenger numbers.

2020 certified 27 field sites -> 9 complete two-way physical railway
segments on 5 actual national routes. Exact source traffic code/label
2015/2020 match at same 9 segments is required. 2015 traffic counts
are 7–8 years BEFORE plant field surveys and do NOT constitute
contemporaneous exposure or a second independent vegetation test.

For each historical year, same 36 dependent physical-segment dyads,
same original botany, same leave-one-PHYSICAL-SEGMENT-out 9 folds,
same geographic/passenger/same-national-line/year-matched baseline,
and all source 48 exact line-restricted historical freight relabellings.
This is a negative control for traffic-year dependence, not causal
evidence of train-mediated seed transport. No p-values.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import (
   query_series,unpack_jsonstat,source_segment_count_rows,spearman_ties_pairs
)
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
   extract_authoritative_stations,source_vegetation_sites,eurostat_2020_segments,
   segment_route_station_join,collapse_exact_reverse_direction_segments,assign_segment
)
from analyze_SI_2020_traffic_vs_2022_2023_flora_exploratory_v232 import (
   SOURCE_PDF_SHA,segment_traffic_sites
)
from validate_SI_actual_2020_freight_vegetation_beta_holdout_v233 import (
   METRICS,build_segment_dyads,observed_adjustment,leave_one_segment_cv,
   exact_within_line_labels
)

EUROSTAT_SHA={
    '2015':'44edbd7a2fc857fe363612ce44c25de6c77d4f83c9d915b3455846c548efe2c3',
    '2020':'05bc92c3555a47947f4acb8cf5aa1ebfe72607487e81713177570cc0b9d5114d',
}

def paired_2015_2020_traffic(segs2020,source2015):
    """The only authorized crossyear key is exact original code + label."""
    all2015=source_segment_count_rows(source2015)
    by2020={x['source_network_segment_code']:x for x in segs2020}
    out={}
    for physical in sorted({m['eurostat_segment_code'] for m in segs2020 if
                            'eurostat_segment_code' in m}):
        raise ValueError("Cannot pass site mappings into original directional 2020 metadata")
    return out

def matched_traffic_for_2020_site_segments(mapped,raw_2020_direction,original_2015_direction):
    old=source_segment_count_rows(original_2015_direction)
    index={x['source_network_segment_code']:x for x in raw_2020_direction}
    result={};unmatched=[]
    for key in sorted({x['eurostat_segment_code'] for x in mapped}):
        if ' + ' not in key:raise ValueError("Source 2020 physical interval is not exact two-direction id pair")
        ids=key.split(' + ')
        if len(ids)!=2 or ids[0]==ids[1]:
            raise ValueError("Not exactly two unique 2020 direction records")
        pieces=[]
        for code in ids:
            row=index.get(code)
            if row is None:raise ValueError("Missing genuine 2020 origin direction code")
            previous=old.get((code,row['source_original_label']))
            if previous is None:
                unmatched.append({'source_segment_code':code,'reason':'2015 code+label not identical'})
                continue
            if 'TRN_GD' not in previous or 'TRN_PAS' not in previous:
                unmatched.append({'source_segment_code':code,'reason':'2015 goods/passenger incomplete'})
                continue
            pieces.append(previous)
        if len(pieces)!=2:continue
        result[key]={'goods':sum(q['TRN_GD'] for q in pieces),
                     'passenger':sum(q['TRN_PAS'] for q in pieces),
                     'directional_source_records':2}
    return result,{'n_original_2020_unique_physical_source_intervals':len({x['eurostat_segment_code'] for x in mapped}),
      '2015_complete_original_both_direction_site_intervals':len(result),
      'original_2015_missing_or_drifted_exact_source_direction_rows':unmatched,
      'not_impute_missing_2015_as_same_as_2020':True}

def source_outcome_comparison(sites,guild):
    dyads,segmeta=build_segment_dyads(sites,guild)
    if len(dyads)!=36 or len(segmeta)!=9:
        raise ValueError("Original v23.3 36 between 9 source physical intervals changed")
    out={}
    for target in METRICS:
        observed=observed_adjustment(dyads,target)
        cv=leave_one_segment_cv(dyads,target)
        if observed['status']!='DESCRIPTIVE_ONLY' or cv is None:
            raise ValueError('Genuine original physical segment inference is rank deficient')
        perms=[]
        for n,shuffle in exact_within_line_labels(segmeta):
            if shuffle is None:
                if n!=48:raise ValueError('Genuine source nine segments did not produce original 48 line-restricted arrangements')
                continue
            c=leave_one_segment_cv(dyads,target,shuffle)
            if c is None:raise ValueError('Historical freight assignment makes source model unidentifiable')
            perms.append(c['gain_baseline_minus_freight_MSE'])
        if len(perms)!=48:raise ValueError('Not all exact 48 original line-restricted labels were tested')
        improvement=cv['gain_baseline_minus_freight_MSE']
        out[target]={
          'unseen_true_physical_segment_fold_count':9,
          'n_original_dependent_physical_pairs':36,
          'baseline_mean_fold_MSE':cv['baseline_mean_fold_mse'],
          'augmented_freight_mean_fold_MSE':cv['freight_augmented_mean_fold_mse'],
          'MSE_gain_positive_is_better':improvement,
          'relative_baseline_reduction_pct':100*improvement/cv['baseline_mean_fold_mse'],
          'n_holdout_folds_improved_by_historical_freight':cv['n_folds_improved_by_historical_freight'],
          'historical_logfreight_difference_adjusted_in_sample_coef':
              observed['historical_logfreight_difference_conditional_coefficient'],
          'line_constrained_exact_48_historical_label_QA':{
              'n_all_original_national_line_label_permutations_including_true':48,
              'n_permutations_at_least_observed_CV_improvement':sum(v>=improvement-1e-12 for v in perms),
              'min_median_max_permuted_gain':[min(perms),median(perms),max(perms)],
              'not_randomized_natural_trial_p_value':True}
          }
    return out,segmeta

def run(original,pdfpath):
    stat,qa=extract_authoritative_stations(pdfpath)
    if qa['sha256_original_pdf_source']!=SOURCE_PDF_SHA:
        raise ValueError('True 2023 station official annex SHA differs')
    sites,meta=source_vegetation_sites(original)
    raw2020,qa20=eurostat_2020_segments()
    if qa20['source']['sha256_original_json_response']!=EUROSTAT_SHA['2020']:
        raise ValueError('Eurostat 2020 true original dataset changed')
    one,joinqa=segment_route_station_join(stat,raw2020)
    two,qa_phys=collapse_exact_reverse_direction_segments(one)
    mapped,site_cov=assign_segment(sites,two)
    if len(mapped)!=27 or site_cov['n_unique_2020_traffic_segments_mapped']!=9:
        raise ValueError('v23.3 original genuine source 27 site × 9 physical interval scope changed')
    # The 2015 source dataset requires identical codes AND exact endpoint
    # text; no fuzzy year-to-year segment matching.
    raw2015,qa15=query_series('rail_tf_ns15_si')
    if qa15['sha256_original_json_response']!=EUROSTAT_SHA['2015']:
        raise ValueError('Eurostat 2015 true original dataset changed')
    rows15,sourceqa15=unpack_jsonstat(raw2015,'2015')
    pairs2015,qamap=matched_traffic_for_2020_site_segments(mapped,raw2020,rows15)
    if len(pairs2015)!=9:
        return {'version':'v23.5','status':'STOP_INCOMPLETE_IDENTICAL_2015_2020_PHYSICAL_SEGMENT_COEFFICIENTS',
          '2015_2020_source_join_QA':qamap,'source_27_sites_9_original_physical_intervals':True,
          'NO_FAKE_MISSING_HISTORICAL_TRAFFIC':True}
    real20,guild=segment_traffic_sites(original,mapped)
    real15=[]
    for site in real20:
        known=pairs2015.get(site['source_segment'])
        if known is None:raise ValueError("Site 2015 historical source incompletely mapped")
        real15.append({**site,'freight':known['goods'],'passenger':known['passenger']})
    meta20={}
    for site in real20:
        meta20[site['source_segment']]=(site['freight'],site['passenger'])
    meta15={k:(v['goods'],v['passenger']) for k,v in pairs2015.items()}
    if set(meta15)!=set(meta20):
        raise ValueError('Physical intervals differ across years')
    goods=[(meta15[k][0],meta20[k][0]) for k in sorted(meta20)]
    passg=[(meta15[k][1],meta20[k][1]) for k in sorted(meta20)]
    results20,m20=source_outcome_comparison(real20,guild)
    results15,m15=source_outcome_comparison(real15,guild)
    # Reproduce v23.3 source known 2020 taxon result before claiming
    # any robustness; abort on source drift even if the new 2015 matches.
    tax20=results20['taxonomic_Jaccard']
    if abs(tax20['baseline_mean_fold_MSE']-.0028816162866543826)>1e-9 or abs(
           tax20['augmented_freight_mean_fold_MSE']-.002566217858035404)>1e-9:
        raise ValueError('2020 v23.3 reference real source heldout coefficients drifted')
    return {'version':'v23.5',
      'status':'REAL_2015_2020_SAME_PHYSICAL_RAIL_SEGMENT_HISTORICAL_FREIGHT_TEMPORAL_ROBUSTNESS',
      'botanical_2022_2023_true_sites_mapped':len(real20),
      'original_national_rail_line_count':len({v['national_line'] for v in m20.values()}),
      'exact_same_2015_2020_two_direction_physical_interval_count':len(meta20),
      'dependently_reused_source_plant_dyads':36,
      'original_2023_official_station_16page_pdf_SHA256':SOURCE_PDF_SHA,
      'original_eurostat_source_sha':EUROSTAT_SHA,
      'traffic_2015_2020_exact_original_code_label_join_QA':qamap,
      'true_2015_2020_source_segment_goods_spearman_rank':spearman_ties_pairs(goods),
      'true_2015_2020_source_segment_passenger_spearman_rank':spearman_ties_pairs(passg),
      'true_2015_2020_goods_annual_9_segment_range':{
          '2015':[min(x for x,y in goods),max(x for x,y in goods)],
          '2020':[min(y for x,y in goods),max(y for x,y in goods)]},
      'same_flora_9segment_GEO_passenger_adjusted_holdout':{
          '2020_historical_traffic':results20,
          '2015_historical_traffic':results15},
      '2020_2015_single_plant_response_temporal_predictor_sensitivity_only':True,
      'hard_limitations':[
        '2015 and 2020 railway movements are historic and 2022–23 field vegetation observed later',
        '2015 and 2020 are two versions of the same predictor on identical botanical response, NOT independent ecological replication',
        'only 9 physical traffic intervals on five national rail lines; 36 segment dyads share endpoints',
        'actual contemporary 2022/23 goods trains and freight origin–destination remain unmeasured',
        'exact national-line constrained traffic permutations are an observational diagnostic, NOT field causal p-values',
        'no independent natural regional flora reference and missing unmapped 62 botanical sites remain unobserved traffic exposure',
        'cannot infer railway dispersal, genotypic adaptation, selection and global biogeographic override']}
def main():
    p=argparse.ArgumentParser()
    for arg in ('official_station_pdf','original_xlsx','source_audit','output'):
        p.add_argument('--'+arg.replace('_','-'),required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.original_xlsx,a.source_audit)
    result=run(original,a.official_station_pdf)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print('V235_EUROSTAT_TIME_SENSITIVITY='+json.dumps(result,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

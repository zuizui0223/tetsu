#!/usr/bin/env python3
"""v23.2 exploratory source-authenticated historical 2020 freight intensity × 2022/23 rail flora.

Nested original sources:
- author Šabić et al. 2026 authentic botanical locations field surveyed 2022/23
- Slovenia government railway operator 2023 official 16-page station km PDF
- Eurostat observed 2020 goods and passenger trains, original opposite
  directions explicitly summed only for complete two-way physical intervals

Original 89 vegetation localities; v23.1 fail-closed unique joins yield
27 eligible real vegetation sites, 9 physical network traffic intervals,
5 national author railway line sections. This is a nonrandom selected
subset, NOT 27 independent station flora samples. Measured 2020 freight
is historical only, not 2022/23 contemporary causal exposure.

Unit for a univariate check is EUROSTAT PHYSICAL SEGMENT (not source
author site) and route/section is an even higher non-independent cluster.
For between-site floristic Jaccard, take averages across each distinct
physical SEGMENT pair (not hundreds of independent site pairs).
Geography and source original same-year, local flora/landuse remain
uncontrolled; never produce causal p-values. Traffic service composition
and 2020 pandemic preclude an evolutionary adaptation mechanism claim.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from analyze_Slovenia_rail_microhabitats_v103 import (
   extract_original_xlsx,validate_original_keys,pair_source,locations_independent
)
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_sites_and_guilds,jaccard,functional_overlap
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    original_gps_coordinates,haversine_km
)
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import spearman_ties_pairs
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
    extract_authoritative_stations,source_vegetation_sites,
    eurostat_2020_segments,segment_route_station_join,
    collapse_exact_reverse_direction_segments,assign_segment
)

SOURCE_PDF_SHA="83f12dfdb91784a653b9a39a8d927a027e7eeb5d8af18ce17b850b5318796425"
CUTOFFS=(0,5,10,20,40,80,160,500)
N_MAPPED=27
N_TRAFFIC=9
N_ORIGINAL_LINES=5

def rank(x,y):
    if len(x)!=len(y) or len(x)<4:return None
    return spearman_ties_pairs(list(zip(x,y)))

def segment_traffic_sites(original,mapped):
    base,groups,meta=source_sites_and_guilds(original)
    gps,qa=original_gps_coordinates(original,base)
    gps_ix={s['site']:s for s in gps}
    loc,zone,sp,tr,cb,qc=validate_original_keys(original)
    effects=locations_independent(pair_source(loc,zone,sp,tr))
    life_ix={str(s['rail_location_id']):s for s in effects}
    sites=[]
    for m in mapped:
        sid=str(m['source_site'])
        if sid not in gps_ix or sid not in life_ix:raise ValueError('Publisher source original site cannot be joined to true eco response')
        s=gps_ix[sid];l=life_ix[sid]
        if str(s['section'])!=m['national_line'] or str(l['section_number'])!=m['national_line']:
            raise ValueError('Source sections inconsistent after traffic physical route match')
        if m['historical_2020_passenger_trains_per_year'] is None:
            raise ValueError('No complete two-way passenger source traffic in accepted mapped sample')
        sites.append({'sid':sid,'section':m['national_line'],
           'year':s['year'],
           'source_original_track_type':s['track_type'],
           'source_segment':m['eurostat_segment_code'],
           'source_segment_endpoints':m['eurostat_source_segment_name'],
           'freight':m['historical_2020_goods_trains_per_year'],
           'passenger':m['historical_2020_passenger_trains_per_year'],
           'annual_share_BS_minus_BT':l['BS_minus_BT_annual_strict'],
           'therophyte_share_BS_minus_BT':l['BS_minus_BT_therophyte_strict'],
           'cover_BS_minus_BT_in_original_source_units':l['BS_minus_BT_cover'],
           'species_richness_BS_minus_BT_UNEQUAL_AREAS':l['BS_minus_BT_richness'],
           'source_species_all':s['all_original_taxa'],
           'known_original_lifespan_diaspore_species':s['known'],
           'lat':s['lat'],'lon':s['lon']})
    if len(sites)!=N_MAPPED or len({x['section'] for x in sites})!=N_ORIGINAL_LINES or len({x['source_segment'] for x in sites})!=N_TRAFFIC:
        raise ValueError('v23.1 genuine traffic-to-source species site universe drifted')
    return sites,groups

def per_physical_segment(sites):
    d=defaultdict(list)
    for x in sites:d[x['source_segment']].append(x)
    rows=[]
    for seg,source in sorted(d.items()):
        lines={x['section'] for x in source}
        if len(lines)!=1:raise ValueError('One Eurostat physical segment mapped to more than one Slovenian official national line')
        freights={x['freight'] for x in source}
        passeng={x['passenger'] for x in source}
        if len(freights)!=1 or len(passeng)!=1:raise ValueError('One physical Eurostat interval has contradictory original traffic value')
        metrics=('annual_share_BS_minus_BT','therophyte_share_BS_minus_BT','cover_BS_minus_BT_in_original_source_units','species_richness_BS_minus_BT_UNEQUAL_AREAS')
        rows.append({'source_physical_segment_codes':seg,
           'source_endpoint_name_both_directions':source[0]['source_segment_endpoints'],
           'national_line':source[0]['section'],
           'source_n_vegetation_field_sites':len(source),
           'source_original_survey_years':dict(Counter(str(x['year']) for x in source)),
           'source_original_track_type':dict(Counter(str(x['source_original_track_type']) for x in source)),
           'source_2020_goods_both_directions':next(iter(freights)),
           'source_2020_passenger_both_directions':next(iter(passeng)),
           'observed_field_source_mean':{
               met:mean(x[met] for x in source if x[met] is not None)
               if any(x[met] is not None for x in source) else None
               for met in metrics}})
    return rows

def association_over_physical_segments(rows):
    out={}
    if len(rows)<6:return {'status':'LESS_THAN_6_COMPLETE_EUROSTAT_PHYSICAL_SEGMENTS_NO_DIAGNOSTIC'}
    for outcome in ('annual_share_BS_minus_BT','therophyte_share_BS_minus_BT',
                    'cover_BS_minus_BT_in_original_source_units',
                    'species_richness_BS_minus_BT_UNEQUAL_AREAS'):
        sub=[r for r in rows if r['observed_field_source_mean'][outcome] is not None]
        y=[r['observed_field_source_mean'][outcome] for r in sub]
        metrics={}
        for source in ('goods','passenger'):
            x=[math.log1p(r['source_2020_'+source+'_both_directions']) for r in sub]
            actual=rank(x,y)
            leave=[rank([x[i] for i in range(len(x)) if i!=z],
                        [y[i] for i in range(len(y)) if i!=z])
                   for z in range(len(x))]
            leave=[v for v in leave if v is not None]
            # Route-level support guard, no unreported pseudo p-values.
            blocked=[]
            for line in sorted({r['national_line'] for r in sub},key=int):
                xi=[math.log1p(r['source_2020_'+source+'_both_directions']) for r in sub if r['national_line']!=line]
                yi=[r['observed_field_source_mean'][outcome] for r in sub if r['national_line']!=line]
                r=rank(xi,yi)
                if r is not None:blocked.append(r)
            metrics[source]={'source_n_distinct_physical_segments':len(sub),
                'n_source_national_route_lines':len({r['national_line'] for r in sub}),
                'exploratory_spearman_on_segment_means':actual,
                'leave_one_segment_out_minmax':[min(leave),max(leave)] if leave else None,
                'leave_one_original_line_out_minmax':[min(blocked),max(blocked)] if blocked else None,
                'NO_independent_biological_p_value':True}
        out[outcome]=metrics
    return out

def betas(sites,guild):
    data=defaultdict(list);bybin=defaultdict(list);same=Counter()
    for i,j in combinations(sites,2):
        if i['source_segment']==j['source_segment']:
            same['same_source_segment_site_pairs_excluded_from_between_segment_beta']+=1
            continue
        geo=haversine_km(i['lat'],i['lon'],j['lat'],j['lon'])
        if geo<2:continue
        freightgap=abs(math.log1p(i['freight'])-math.log1p(j['freight']))
        passengap=abs(math.log1p(i['passenger'])-math.log1p(j['passenger']))
        jk=jaccard(i['source_species_all'],j['source_species_all'])
        fun=functional_overlap(i['known_original_lifespan_diaspore_species'],
                               j['known_original_lifespan_diaspore_species'],guild)
        pair=tuple(sorted((i['source_segment'],j['source_segment'])))
        band=next((f"{lo}-{hi}" for lo,hi in zip(CUTOFFS[:-1],CUTOFFS[1:]) if lo<=geo<hi),None)
        x={'jaccard':jk,'functional_overlap':fun,'geo_km':geo,
           'log_goods_movement_count_difference':freightgap,
           'log_passenger_movement_count_difference':passengap,
           'same_original_national_line':i['section']==j['section'],
           'same_original_field_year':i['year']==j['year']}
        data[pair].append(x)
        if band is not None:bybin[(pair,band)].append(x)
    pairs=[]
    for key,vals in sorted(data.items()):
        pairs.append({'true_Eurostat_physical_segment_pair':list(key),
           'source_dependent_original_site_pairs':len(vals),
           'mean_species_Jaccard':mean(x['jaccard'] for x in vals),
           'mean_original_four_guild_overlap':mean(x['functional_overlap'] for x in vals),
           'mean_geodesic_distance_km':mean(x['geo_km'] for x in vals),
           'mean_abs_log_goods_traffic_difference':mean(x['log_goods_movement_count_difference'] for x in vals),
           'mean_abs_log_passenger_traffic_difference':mean(x['log_passenger_movement_count_difference'] for x in vals),
           'same_source_national_line':all(x['same_original_national_line'] for x in vals)})
    n=len(pairs)
    if n<8:return {'source_physical_traffic_segment_pairs':n,'status':'UNSUPPORTED_FEWER_THAN_8_SEGMENT_PAIR_UNITS'}
    bx=[r['mean_species_Jaccard'] for r in pairs]
    fg=[r['mean_abs_log_goods_traffic_difference'] for r in pairs]
    pg=[r['mean_abs_log_passenger_traffic_difference'] for r in pairs]
    geo=[r['mean_geodesic_distance_km'] for r in pairs]
    fun=[r['mean_original_four_guild_overlap'] for r in pairs]
    out={'source_n_traffic_physical_segment_pairs':n,
         'source_n_all_dependent_original_site_pair_comparisons':sum(x['source_dependent_original_site_pairs'] for x in pairs),
         'site_pairs_excluded_same_physical_segment':dict(same),
         'spearman_species_Jaccard_vs_abs_loggoods_difference':rank(fg,bx),
         'spearman_species_Jaccard_vs_abs_logpassenger_difference':rank(pg,bx),
         'spearman_species_Jaccard_vs_geodesic_distance':rank(geo,bx),
         'spearman_functional_four_guild_overlap_vs_abs_loggoods_difference':rank(fg,fun),
         'n_traffic_pairs_same_national_route':sum(x['same_source_national_line'] for x in pairs)}
    # Geographical bin poststratification, descriptive only. Each original
    # physical pair may enter different geographic bins when multiple sites
    # are mapped to each physical segment, and is dependent across bins.
    bins={}
    for (key,band),val in bybin.items():
        bins.setdefault(band,[]).append({
           'source_physical_segment_pair':key,
           'jaccard':mean(x['jaccard'] for x in val),
           'guild':mean(x['functional_overlap'] for x in val),
           'abs_log_goods_diff':mean(x['log_goods_movement_count_difference'] for x in val),
           'geographic_km':mean(x['geo_km'] for x in val)})
    out['geo_bin_QA_descriptive']={
        band:{'source_n_segment_pairs_in_bin':len(v),
              'rho_Jaccard_vs_abs_loggoods_diff':
               rank([p['abs_log_goods_diff'] for p in v],[p['jaccard'] for p in v])
               if len(v)>=6 else None}
        for band,v in sorted(bins.items())}
    # Without physical freight OD, pairwise "similar freight intensity"
    # does NOT represent trains moving from one source segment to another.
    out['source_intersegment_train_freight_OD_UNMEASURED']=True
    out['no_pseudo_independent_sitepair_pvalue']=True
    return out

def analysis(original,pdf):
    stat,qa=extract_authoritative_stations(pdf)
    if qa['sha256_original_pdf_source']!=SOURCE_PDF_SHA:
        raise ValueError('Genuine 2023 station PDF changed')
    original_sites,gps=source_vegetation_sites(original)
    traffic,eq=eurostat_2020_segments()
    matched,mq=segment_route_station_join(stat,traffic)
    phys,pqa=collapse_exact_reverse_direction_segments(matched)
    mapped,cov=assign_segment(original_sites,phys)
    if len(mapped)!=27 or cov['n_unique_2020_traffic_segments_mapped']!=9:
        raise ValueError('Original v23.1 source joined 27 sites /9 physical segments changed')
    paired,groups=segment_traffic_sites(original,mapped)
    segment=per_physical_segment(paired)
    if len(segment)!=9:raise ValueError('Original independent physical segment count changed')
    return {'version':'v23.2',
      'status':'REAL_HISTORICAL_2020_TRAIN_INTENSITY_EXPLORATORY_NOT_2023_CAUSAL',
      '2023_rail_operator_official_station_source_sha':SOURCE_PDF_SHA,
      '2020_eurostat_actual_traffic_response_sha':eq['source']['sha256_original_json_response'],
      'source_actual_plant_survey_years':'2022–2023',
      'source_true_original_site_mapping_and_network_coverage':cov,
      'source_observed_2020_two_way_traffic_segment_level_2022_23_vegetation_effects':segment,
      'descriptive_segment_unit_2020_traffic_vs_BS_BT_original_field_traits':association_over_physical_segments(segment),
      'historical_2020_train_intensity_differences_vs_between_segment_plant_species_beta':betas(paired,groups),
      'no_2022_23_actual_service_data':True,
      'no_freight_origination_destination':True,
      'limitations':[
          'This is a year-mismatched historical transport proxy, not 2022–2023 contemporaneous traffic',
          'Site selection into 27/89 mapped sites nonrandom and depends on legibility/availability of station endpoints',
          '9 source physical segments nested in 5 original railway lines within only one national observational survey',
          'Each original author site and location has repeated species/zone observations, not station-wide area-standardized independent botanical sampling',
          'Route-level spatial, substrate, mowing and land-use effects cannot be controlled without comparable nonrail plants',
          'Actual train seed cargo, freight OD and propagule establishment are not observed',
          'Train number vs railway microhabitat is not a measured genetic adaptation or plasticity test',
          'Taxonomic pair dependence not resolved by averaging over the same source traffic segment multiple times',
          'No causal estimates, no bootstrap-significance claims, no inference to unmapped source sections']}
def main():
    p=argparse.ArgumentParser()
    for k in ('official_station_pdf','original_xlsx','source_audit','output'):
        p.add_argument('--'+k.replace('_','-'),required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.original_xlsx,a.source_audit)
    out=analysis(original,a.official_station_pdf)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    print('V232_REAL_2020_TRAIN_2022_23_PLANT_SOURCE='+json.dumps(out,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

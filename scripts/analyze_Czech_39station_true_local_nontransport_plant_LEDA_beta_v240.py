#!/usr/bin/env python3
"""v24.0 genuine Czech 39 published railway station floras vs independent
nearby nontransport CzechVeg botanical plots (same local geography).

Actual separately sourced data:
 A: Preslia 2024 39 station pooled flora 2020–22, original S1 rec + S2,
    publisher SHA locked. B: 2026 Czech SR70 official-derived station
    CENTROID points, 35 exact full-name unique stations matched to A.
 C: CzechVeg-Open 2026 v1.1 original 117739 geo headers (MD5) + 2.51M
    vascular taxa plot rows (MD5). D: original LEDA canopy height +
    seed mass species medians, SHA anchored.

Main external local reference: 2015–2023 observation year, WGS84
 coordinate source precision STRICTLY POSITIVE <=250m (do NOT trust 0
 as exact), source plot nontransport HABITAT+LOCALITY text, <=10km of
 2026 railway station registered representative point. Requires >=5
 source plots from >=2 original studies. Origin 35/39 station geocodes;
 original v23.9b showed 19 stations satisfying metadata support.

This is NOT a paired railway-site vs equal-area natural-plot survey:
 publisher railway station source pools 3 habitats over 2–90 hectares,
 while "nontransport" source relevés are small, heterogeneously sampled
 and NOT independently proven truly natural. Never infer causal rail
 filtering, biological convergence or adaptive evolution.

For each station, aggregate 5 control vegetation PLOTS selected within
 10km from >=2 different original publication studies, unioning their
 real vascular species names. Original station source taxa are fixed.
 Compute species-equal CWM in z log canopy height & z log seed mass.

Across EXACT same matched source railway station pairs:
 - CWM centroid pairwise Euclidean beta rail minus local control.
 - Same station pair mean source Jaccard rail vs controls.
 - mean station/community 2D trait alpha trace.

Draw 299 technical matched reference pools (5 plots) and sitewise
 source-rarefaction to min(station original known taxa, five-plot control
 unique known taxa) in BOTH arms to mitigate huge sampling grain.
 This does NOT fully correct unequal area, vegetation type or effort.
 Report source region and original railway-stretch leave-one-out
 descriptive support, not falsely independent 171 station-pair p-values.

Publish only aggregate SHA/support and no original source named
 location IDs, plot coordinates, or raw species×site data.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from audit_Czech_2024_S2_to_SR70_original_station_points_v153 import (
    run as original_station_geolocation
)
from analyze_Czech_2024_39station_functional_identity_v143 import (
    read_original_czech_publisher
)
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import get_meta
from audit_CzechVeg_abiotic_regional_reference_opportunity_v194 import number
from analyze_CzechVeg_original_independent_reference_taxonomic_v191 import (
    source_species,jaccard
)
from audit_Czech_35SR70_station_local_CzechVeg_nontransport_georeference_v239 import (
    REQUIRED_GEOGRAPHIC_SOURCE_COLUMNS,CZECH_BOUNDS
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_CzechVeg_original_LEDA_rail_text_vs_abiotic_reference_functional_v238 import (
    exact_source_name
)
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    components,quart
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    haversine_km
)

LEDA_SHA={
 "height_m":"cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92",
 "seed_mass_mg":"34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef"}
MAX_RADIUS_KM=10
MIN_SOURCE_YEAR=2015
MAX_SOURCE_YEAR=2023
MIN_SOURCE_PRECISION_M=0
MAX_SOURCE_PRECISION_M=250
MIN_CTRL_PLOTS=5
MIN_CONTROL_ORIGINAL_PAPERS=2
MIN_TRAIT_CODES_PER_PLOT=5
MIN_TRAIT_COVERAGE=0.55
DRAWS=299
SEED=20261117+240

def source_station_and_reference_coords(s1,s2,sr70,headers):
    geo=original_station_geolocation(s1,s2,sr70)
    source_stations,src_stretch,srcqa=read_original_czech_publisher(s1,s2)
    active=[]
    for s in geo['all_39_station_original_source_identity_and_candidate_provenance']:
        xy=s['matched_czech_2026_source_station_centroid_WGS84']
        if xy is None:continue
        code=s['original_locality_code']
        if code not in source_stations:raise ValueError('Original Preslia 39 station source code missing')
        active.append({'id':code,'stretch':src_stretch[code]['stretch'],
                       'original_taxa':source_stations[code],
                       'lat':xy['lat'],'lon':xy['lon']})
    if len(active)!=35:raise ValueError('Original v15.3 exactly 35 official SR70 station geo locations drifted')
    meta,hq=get_meta(headers)
    good={}
    counts=Counter()
    with open(headers,encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f)
        if set(REQUIRED_GEOGRAPHIC_SOURCE_COLUMNS)-set(r.fieldnames or []):
            raise ValueError('True CzechVeg WGS84 decimal degree/precision fields absent')
        for row in r:
            counts['original_CzechVeg_total_plot_headers']+=1
            id=str(row.get('RELEVE_NR') or '').strip()
            m=meta.get(id)
            if m is None:continue
            if m['label'] is not None or m['strict'] is not None:continue
            if not MIN_SOURCE_YEAR<=m['year']<=MAX_SOURCE_YEAR:continue
            lat=number(row.get('DEG_LAT'),-90,90)
            lon=number(row.get('DEG_LON'),-180,180)
            precision=number(row.get('PRECISION'),0,1000000)
            if lat is None or lon is None or precision is None:continue
            if not(MIN_SOURCE_PRECISION_M<precision<=MAX_SOURCE_PRECISION_M):continue
            if not(CZECH_BOUNDS[0]<=lat<=CZECH_BOUNDS[1] and
                   CZECH_BOUNDS[2]<=lon<=CZECH_BOUNDS[3]):continue
            good[id]={'lat':lat,'lon':lon,'year':m['year'],
                      'study':m['study'],'source_phyto_region':m['region']}
    if counts['original_CzechVeg_total_plot_headers']!=117739:
        raise ValueError('Original published CzechVeg header row count differs')
    control=defaultdict(list)
    for station in active:
        for id,p in good.items():
            if abs(station['lat']-p['lat'])>MAX_RADIUS_KM/110+.01:continue
            km=haversine_km(station['lat'],station['lon'],p['lat'],p['lon'])
            if km<=MAX_RADIUS_KM:
                control[station['id']].append((id,km))
    accepted=[]
    for st in active:
        ids=control[st['id']]
        if len(ids)<MIN_CTRL_PLOTS:continue
        papers={good[x]['study'] for x,d in ids}
        if len(papers)<MIN_CONTROL_ORIGINAL_PAPERS:continue
        accepted.append({**st,'metadata_candidate_plot_ids':ids})
    if len(accepted)!=19:
        raise ValueError('v23.9b original 19/35 true source 2015+ 10km positive precision station support drifted; found '+str(len(accepted)))
    toread={id for st in accepted for id,d in st['metadata_candidate_plot_ids']}
    return accepted,good,toread,{
        'original_rail_station_source_publisher_39':len(source_stations),
        'original_verified_true_station_centroids':len(active),
        'n_station_centroids_with_5_original_nontransport_candidate_plots_2_studies':len(accepted),
        'n_unique_reference_original_plant_plot_ids_preselected_before_species_content':len(toread),
        'original_source_railway_stretches_supported':dict(Counter(st['stretch'] for st in accepted)),
        'original_source_phyto_region_reference_counts':dict(Counter(good[id]['source_phyto_region'] for id in toread)),
        'original_CzechVeg_headers_MD5':hq['publisher_original_headers_md5'],
        'original_source_plots_coordinate_precision_in_m_positive_and_le250':True,
        '2026_station_registry_representative_point_not_2020_original_botany_GPS':True}

def source_species_to_traits(stations,original_plants,height,seed):
    station_names=set().union(*(s['original_taxa'] for s in stations))
    ref_names=set().union(*(a for x in original_plants.values() for a in x)) if original_plants else set()
    universe=station_names|ref_names
    joined={}
    for name in universe:
        b=exact_source_name(name)
        if b and b in height and b in seed:
            joined[name]=(math.log(height[b]),math.log(seed[b]))
    if len(joined)<100:raise ValueError('Original station and ref plot joint exact LEDA trait universe too small')
    mu=[mean(v[j] for v in joined.values()) for j in (0,1)]
    sd=[math.sqrt(mean((v[j]-mu[j])**2 for v in joined.values())) for j in (0,1)]
    if min(sd)<=0:raise ValueError('One original source standardized trait has zero variance')
    traits={key:tuple((raw[j]-mu[j])/sd[j] for j in (0,1)) for key,raw in joined.items()}
    return traits,{
       'unique_original_source_taxa_rail_station_and_local_ref_combined':len(universe),
       'unique_unambiguous_taxa_joined_source_LEDA_height_seed':len(traits),
       'n_original_study_station_taxa_matched_trait':len(station_names&set(traits)),
       'n_original_independent_nontransport_reference_taxa_matched_trait':len(ref_names&set(traits)),
       'trait_standardization_z_from_union_unique_taxa_not_group_specific':True}

def source_valid_local_control_plots(stations,sourceplants,geo,traits):
    supported=[];qc=Counter()
    for s in stations:
        options=[]
        for id,d in s['metadata_candidate_plot_ids']:
            spp=sourceplants.get(id,set())
            taxa=spp&set(traits)
            if len(taxa)<MIN_TRAIT_CODES_PER_PLOT:
                qc['source_reference_plot_rejected_joint_trait_species_lt5']+=1;continue
            if len(taxa)/len(spp)<MIN_TRAIT_COVERAGE:
                qc['source_reference_plot_rejected_trait_coverage_below55pct']+=1;continue
            options.append({'id':id,'d':d,'study':geo[id]['study'],'taxa':set(taxa),
                            'n_all_source_taxa':len(spp)})
        stat=s['original_taxa']&set(traits)
        if len(stat)<30:raise ValueError('Original 2020/22 source Preslia railway station unexpectedly few valid LEDA trait species')
        if len(options)>=5 and len({p['study'] for p in options})>=2:
            supported.append({**s,'control_options':options,'rail_source_trait_taxa':stat})
        else:
            qc['station_lost_after_actual_source_species_and_LEDA_QA']+=1
    return supported,{
       'original_19_station_metadata_eligible':len(stations),
       'n_stations_with_at_least5_real_plots_and_2_original_studies_after_traitQC':len(supported),
       'station_count_by_original_rail_stretch':dict(Counter(s['stretch'] for s in supported)),
       'actual_source_reference_plot_rejections':dict(qc),
       'n_selected_source_station_known_joint_traits_min_median_max':(
           [min(len(s['rail_source_trait_taxa']) for s in supported),
            median(len(s['rail_source_trait_taxa']) for s in supported),
            max(len(s['rail_source_trait_taxa']) for s in supported)] if supported else None),
       'no_species_imputed_at_genus_level':True}

def sample_local_plots(station,rng):
    candidates=station['control_options']
    papers=defaultdict(list)
    for row in candidates:papers[row['study']].append(row)
    if len(candidates)<5 or len(papers)<2:
        raise ValueError('Source local plot reference must contain 5 distinct field plots and at least 2 original publications')
    # Force two different source publications BEFORE inspecting selected
    # source species; choose 3 others from the remaining true plot IDs.
    selected_papers=rng.sample(sorted(papers),2)
    chosen=[rng.choice(papers[paper]) for paper in selected_papers]
    first_ids={x['id'] for x in chosen}
    remainder=[x for x in candidates if x['id'] not in first_ids]
    if len(remainder)<3:
        raise ValueError('Original five-plot control without reuse impossible')
    selected=chosen+rng.sample(remainder,3)
    if len({x['id'] for x in selected})!=5 or len({x['study'] for x in selected})<2:
        raise AssertionError('Fixed original publication or unique control plot selection failed')
    return selected

def source_profile(codes,traits):
    d=components(codes,traits)
    return {'CWM':(d['height_mean'],d['seed_mean']),
            'trace':d['joint_trace'],'taxa':set(codes)}

def source_beta_stats(original,reference):
    if set(original)!=set(reference):raise ValueError('Source station cohort not matched within same exact registry centroids')
    ids=sorted(original)
    diffs={'rail_cwm':[],'ref_cwm':[],
           'rail_jaccard':[],'ref_jaccard':[],
           'rail_trace':[],'ref_trace':[]}
    for id in ids:
        diffs['rail_trace'].append(original[id]['trace'])
        diffs['ref_trace'].append(reference[id]['trace'])
    for i,j in combinations(ids,2):
        x,y=original[i],original[j]
        a,b=reference[i],reference[j]
        diffs['rail_cwm'].append(math.dist(x['CWM'],y['CWM']))
        diffs['ref_cwm'].append(math.dist(a['CWM'],b['CWM']))
        diffs['rail_jaccard'].append(jaccard(x['taxa'],y['taxa']))
        diffs['ref_jaccard'].append(jaccard(a['taxa'],b['taxa']))
    return {
       'rail_minus_local_nontransport_CWM_beta':mean(diffs['rail_cwm'])-mean(diffs['ref_cwm']),
       'rail_minus_local_nontransport_Jaccard':mean(diffs['rail_jaccard'])-mean(diffs['ref_jaccard']),
       'rail_minus_local_nontransport_2D_alpha_trace':mean(diffs['rail_trace'])-mean(diffs['ref_trace']),
       'absolute_rail_CWM_distance':mean(diffs['rail_cwm']),
       'absolute_local_control_CWM_distance':mean(diffs['ref_cwm']),
       'source_n_dependent_station_pairs':len(diffs['rail_cwm']),
       'n_real_rail_stations':len(ids)}

def sample_frame(stations,traits,rng,rarefaction):
    rail={};ref={};selected_control_plot_ids=[]
    for s in stations:
        controls=sample_local_plots(s,rng)
        selected_control_plot_ids.extend(x['id'] for x in controls)
        pool=set().union(*(p['taxa'] for p in controls))
        staxa=s['rail_source_trait_taxa']
        if len(pool)<5 or len(staxa)<5:
            raise ValueError('Source reference 5-plot union not enough actual trait-known original species')
        if rarefaction:
            k=min(len(staxa),len(pool))
            staxa=set(rng.sample(sorted(staxa),k))
            pool=set(rng.sample(sorted(pool),k))
        rail[s['id']]=source_profile(staxa,traits)
        ref[s['id']]=source_profile(pool,traits)
    out=source_beta_stats(rail,ref)
    out['source_reference_selected_plot_n']=len(selected_control_plot_ids)
    out['source_reference_selected_unique_plot_n']=len(set(selected_control_plot_ids))
    out['n_actual_plot_reuses_across_station_references']=len(selected_control_plot_ids)-len(set(selected_control_plot_ids))
    return out

def draw_panel(stations,traits,rarefaction,seed,draws=DRAWS):
    if len(stations)<10 or len({s['stretch'] for s in stations})<3:
        return {'identified':False,'n_sites':len(stations),
            'n_true_original_railway_stretches':len({s['stretch'] for s in stations}),
            'reason':'Fewer than 10 original railway stations or fewer than three true source stretches after original reference species quality gates'}
    rng=random.Random(seed)
    draw=[sample_frame(stations,traits,rng,rarefaction) for i in range(draws)]
    focus=('rail_minus_local_nontransport_CWM_beta',
           'rail_minus_local_nontransport_Jaccard',
           'rail_minus_local_nontransport_2D_alpha_trace')
    res={}
    for met in focus+('absolute_rail_CWM_distance','absolute_local_control_CWM_distance'):
        vals=[x[met] for x in draw]
        res[met]={'mean':mean(vals),
          'q025_technical_not_field_CI':quart(vals,.025),
          'q975_technical_not_field_CI':quart(vals,.975),
          'technical_draws':draws,
          'fraction_negative':sum(v<0 for v in vals)/len(vals)}
    reuse=[d['n_actual_plot_reuses_across_station_references'] for d in draw]
    loo={}
    original_stretches=sorted({s['stretch'] for s in stations})
    for i,stretch in enumerate(original_stretches):
        rest=[s for s in stations if s['stretch']!=stretch]
        if len(rest)<8 or len({s['stretch'] for s in rest})<3:continue
        rng2=random.Random(seed+9917+i)
        x=[sample_frame(rest,traits,rng2,rarefaction) for z in range(59)]
        loo[stretch]={m:mean(q[m] for q in x) for m in focus}
    return {'identified':True,
      'original_observed_Czech_author_stations_in_panel':len(stations),
      'n_source_original_railway_stretches':len(original_stretches),
      'original_source_geographical_reference_radius_km':MAX_RADIUS_KM,
      'source_local_reference_plots_per_station':5,
      'original_publication_study_types_min_per_reference':2,
      'rarefaction_to_equal_source_LEDA_known_richness_in_both_arms':rarefaction,
      'station_pair_beta_observed_vs_local_reference':res,
      'source_reference_plot_reuse_QA':{
         'min':min(reuse),'median':median(reuse),'max':max(reuse),
         'draws_with_any_plot_reused_between_two_stations':sum(v>0 for v in reuse)},
      'source_leave_one_stretch_out_comparison':loo,
      'no_rail_causal_effect_or_independent_field_CI':True}

def run(s1,s2,sr70,header,species,height,seed,draws=DRAWS):
    src,geoplots,toread,meta=source_station_and_reference_coords(s1,s2,sr70,header)
    plants,speciesqa=source_species(species,toread)
    trait,zq=source_species_to_traits(src,plants,height,seed)
    station,quality=source_valid_local_control_plots(src,plants,geoplots,trait)
    if len(station)<10 or len({s['stretch'] for s in station})<3:
        panel={'status':'INSUFFICIENT_FIELD_REFERENCE_SUPPORT','identified':False}
        return {'version':'v24.0','status':'SOURCE_HARD_FIELD_SUPPORT_GATE_NOT_PASSED',
            'original_Czech_Preslia_CzechVeg_local_species_reference_metadata':meta,
            'source_actual_species_by_reference_QA':speciesqa,
            'source_LEDA_exact_height_seed_join_QA':zq,
            'source_local_reference_source_quality':quality,
            'biological_analysis':panel,
            'source_reference_missing_NOT_ZERO_RAIL_EFFECT':True}
    result={'version':'v24.0',
       'status':'AUTHENTIC_CZECH39_RAIL_STATIONS_VS_LOCAL_VEGETATION_REFERENCE_TRAIT_BETA_DESCRIPTIVE',
       'original_Czech_Preslia_CzechVeg_local_species_reference_metadata':meta,
       'source_actual_species_by_reference_QA':speciesqa,
       'source_LEDA_exact_height_seed_join_QA':zq,
       'source_local_reference_source_quality':quality,
       'biological_analysis':{
          'source_5_plot_union_fixed_railway_actual_taxa':draw_panel(
             station,trait,False,SEED,draws),
          'source_5_plot_union_species_count_rarefied_RAIL_AND_CONTROL':draw_panel(
             station,trait,True,SEED+51,draws),
       },
       'strict_inference_limits':[
          'Preslia 39 stations are 3-zone entire station floras, not same-area plots with CzechVeg nontransport source relevés',
          'Author Preslia plant field time 2020/22 vs source CzechVeg 2015/23, sometimes earlier or later',
          '2026 SR70 station point is not original Preslia botanical plot GPS',
          'CzechVeg missing transport words are not proof of genuine natural or rail-free habitat',
          'Five reference plots' + ' pooled species lists are a different sampling process from single author railway station-wide census',
          'Matching the number of LEDA-known taxa cannot eliminate habitat selection, sample area, species frequency or study protocol',
          'Control relevés may be reused across nearby original railway station matching radii; the recorded reuse is a confounding diagnostic',
          'Dependent station pairs, shared author source stretches, reused reference plots and technical draws are not independent field studies',
          'Physical soil/gravel/maintenance, contemporary freight OD and plant hereditary adaptation remain unobserved',
          'No claim of global railway species or trait convergence, railway adaptive effect or seed train transport supported']}
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ('s1','s2','sr70','headers','species','height','seed','output'):
        p.add_argument('--'+name,required=True)
    p.add_argument('--draws',type=int,default=DRAWS)
    a=p.parse_args()
    for file,sha in ((a.height,LEDA_SHA['height_m']),(a.seed,LEDA_SHA['seed_mass_mg'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=sha:
            raise ValueError('Public LEDA source original height/seed SHA changed')
    h,hq=read_leda(a.height,'height_m')
    seed,sqa=read_leda(a.seed,'seed_mass_mg')
    out=run(a.s1,a.s2,a.sr70,a.headers,a.species,h,seed,a.draws)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    small={'version':'v24.0','status':out['status'],
      'original_source_support':out['original_Czech_Preslia_CzechVeg_local_species_reference_metadata'],
      'source_actual_species_plant_rows':out['source_actual_species_by_reference_QA'],
      'source_traits_QA':out['source_LEDA_exact_height_seed_join_QA'],
      'source_reference_QA':out['source_local_reference_source_quality'],
      'biological_analysis':out['biological_analysis']}
    print('V240_CZECH_REAL_STATION_LOCAL_REFERENCE_BETA='+json.dumps(small,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

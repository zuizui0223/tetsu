#!/usr/bin/env python3
"""v24.3: maximum-source-support sensitivity to independent feasible control assignments.

SOURCE: actual Preslia 2020-22 station flora (39), 35 conservatively
georegistered 2026 station representative points, CzechVeg 2015-23
independent original vegetation plots <=10km (original PRECISION>0
<=250m), independently original 2+ papers / station and exactly five
unique original plot IDs per selected station, and source-authenticated
LEDA height+seed species medians.

Reference v24.2 exact mixed-integer MILP optimum is 14 of 16 eligible
stations; 70 distinct actual local vegetation plots across four rail
stretches. Its 199 source species rarefaction draws fixed ONE choice.
This run seeks multiple different OPTIMAL-SUPPORT 14-station assignments
using randomized, BOTANY-BLIND costs only based on source true distance,
source study/plot opportunity IDs and seeded tie-breaking. Each 14-station
solution is certified feasible by binary MILP, although the 14 maximum
is certified separately by v24.2. No plant species/trait outcome enters
the matching objective.

Each assignment gives mean BS?? No: this is NOT Slovenia paired BS/BT.
It compares species-equal 2D CWM centroid distance across true Czech
RAILWAY STATIONS to distance across nearby independently sampled
NONTRANSPORT-TEXT VEGETATION-POOL UNIONS (five source plots/station).
Both samples are rarefied to the SAME species count within each
station. 20 technical rarefaction draws per independent source assignment.
Leave each original Preslia railway stretch out and repeat cross-station
functional beta in that same true source matched subset.

Additional support:
 - number of unique original assignment solutions
 - n originally distinct plot IDs per assignment (must equal 70)
 - exact source 4-decimal duplicate geocoordinate cells counts (same
   plot may be distinct ID at same original source point)
 - n distinct original CzechVeg publisher studies per assignment
 - reference source sample community richness relative to original
   station source trait-known species richness (sampling-grain QA).

None of these data establish causality, physical exclusion of true rail
plot from CzechVeg nontransport source, equal sampling area, identical
habitat/management, or railway seed train-mediated connectivity.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_Czech_39station_true_local_nontransport_plant_LEDA_beta_v240 import (
    LEDA_SHA, source_station_and_reference_coords,source_species_to_traits,
    source_valid_local_control_plots,source_species,
    source_profile,source_beta_stats,quart
)
from audit_Czech_local_vegetation_exact_MILP_nonoverlap_capacity_v242 import (
    solve_full_source_global_exclusivity
)

SOURCE_MAX_SUPPORT=14
MIN_SOURCE_LINE_STRETCHES=3
ALTERNATIVE_OBJECTIVE_SEEDS=20
RAREFIED_DRAWS=16
ORIGINAL_SEED=20261117+243

def solve_alternate_maximum(stations,source_max=SOURCE_MAX_SUPPORT,seed=ORIGINAL_SEED,
                            time_limit_s=25):
    """Exact integer solution with unchanged biological/metadata constraints,
    fixed number of max-supported stations and independent plot selection
    objective randomization using reference opportunity metadata only.
    """
    import numpy as np
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    n=len(stations)
    if n!=16 or source_max<10:raise ValueError("Original source station availability changed")
    rng=random.Random(seed)
    idx={}
    flat=[]
    for i,station in enumerate(stations):
        for p in sorted(station['control_options'],key=lambda x:x['id']):
            ix=len(flat);flat.append((i,p))
            idx[i,p['id']]=ix
    studykeys={}
    for i,station in enumerate(stations):
        for study in sorted({p['study'] for p in station['control_options']}):
            studykeys[i,study]=len(flat)+n+len(studykeys)
    zs={i:len(flat)+i for i in range(n)}
    stretches=sorted({s['stretch'] for s in stations})
    ws={g:len(flat)+n+len(studykeys)+j for j,g in enumerate(stretches)}
    total=len(flat)+n+len(studykeys)+len(ws)
    cost=np.zeros(total,dtype=float)
    # Primary fixed sum(z)=14, objective only chooses among valid
    # full-support source allocations. No original plant names or
    # original functional beta can enter this cost.
    for ix,(i,p) in enumerate(flat):
        cost[ix]=rng.uniform(-.45,.45)+p['d']*.03
    for i,ix in zs.items():
        cost[ix]=rng.uniform(-.01,.01)
    groupplot=defaultdict(list)
    for (i,p),ix in zip(flat,range(len(flat))):
        groupplot[p['id']].append(ix)
    rows=[];cols=[];data=[];low=[];high=[]
    def con(seq,lo=-math.inf,hi=math.inf):
        row=len(low)
        for ix,coef in seq:
            rows.append(row);cols.append(ix);data.append(coef)
        low.append(lo);high.append(hi)
    for i in range(n):
        con([(idx[i,p['id']],1) for p in stations[i]['control_options']]+[(zs[i],-5)],lo=0,hi=0)
    for plot,arr in groupplot.items():con([(ix,1) for ix in arr],hi=1)
    for (i,study),iy in studykeys.items():
        cols_for_study=[idx[i,p['id']] for p in stations[i]['control_options'] if p['study']==study]
        con([(x,1) for x in cols_for_study]+[(iy,-1)],lo=0)
        con([(x,1) for x in cols_for_study]+[(iy,-5)],hi=0)
    for i in range(n):
        con([(iy,1) for (j,paper),iy in studykeys.items() if j==i]+[(zs[i],-2)],lo=0)
    for stretch,iw in ws.items():
        con([(zs[i],1) for i,s in enumerate(stations) if s['stretch']==stretch]+[(iw,-1)],lo=0)
    con([(ix,1) for ix in ws.values()],lo=MIN_SOURCE_LINE_STRETCHES)
    # Certificate from original v24.2: maximal station number 14.
    # Here force 14 EXACTLY, never silently choose 12/13 easier stations.
    con([(ix,1) for ix in zs.values()],lo=source_max,hi=source_max)
    mat=coo_matrix((data,(rows,cols)),shape=(len(low),total)).tocsr()
    result=milp(c=cost,integrality=np.ones(total),
        bounds=Bounds(np.zeros(total),np.ones(total)),
        constraints=LinearConstraint(mat,np.array(low),np.array(high)),
        options={'time_limit':time_limit_s,'presolve':True,'mip_rel_gap':0.0})
    if result.x is None or int(result.status)!=0:
        return None,{"status":"NO_OPTIMAL_ALTERNATIVE_IDENTIFIED",
            "status_code":int(result.status),"solver_message":str(result.message)[:220],
            "not_a_biological_null_effect":True}
    selected={i:[p for j,p in flat if j==i and result.x[idx[j,p['id']]]>.5]
              for i,j in zs.items() if result.x[j]>.5}
    if len(selected)!=source_max:
        raise AssertionError("Exact max-supported 14 station preservation failed")
    allplot=[p['id'] for arr in selected.values() for p in arr]
    if len(allplot)!=5*source_max or len(set(allplot))!=len(allplot):
        raise AssertionError("Exact source global original reference plot disjointness violated")
    if any(len(arr)!=5 or len({p['study'] for p in arr})<2 for arr in selected.values()):
        raise AssertionError("Original source five plot and at least two papers constraint violated")
    if len({stations[i]['stretch'] for i in selected})<MIN_SOURCE_LINE_STRETCHES:
        raise AssertionError("Original national railway stretch support lost")
    allocation_digest=hashlib.sha256("|".join(
        f'{stations[i]["id"]}:{",".join(sorted(p["id"] for p in arr))}'
        for i,arr in sorted(selected.items())).encode()).hexdigest()
    return selected,{"status":"FULL_SUPPORT_ALTERNATIVE_EXACT_MILP_FEASIBLE",
        "original_source_stations":len(selected),
        "true_source_distinct_releve_ids":len(set(allplot)),
        "original_preslia_rail_stretch_counts":dict(Counter(stations[i]['stretch'] for i in selected)),
        "original_reference_distinct_publisher_studies":len({p['study'] for arr in selected.values() for p in arr}),
        "botany_blind_allocation_sha256":allocation_digest,
        "source_objective_uses_only_plot_distance_random_tie_breaks":True,
        "MIP_optimality_gap":float(getattr(result,'mip_gap',0) or 0),
        "source_binary_decision_variables":total}

def station_reference_profiles(stations,assignments,geo,traits,rarefy_seed):
    rng=random.Random(rarefy_seed)
    rail={};ref={}
    original_rich=[];local_rich=[];physical_grid=Counter()
    for i,plots in sorted(assignments.items()):
        station=stations[i]
        codes=station['rail_source_trait_taxa']
        control=set().union(*(p['taxa'] for p in plots))
        original_rich.append(len(codes))
        local_rich.append(len(control))
        size=min(len(codes),len(control))
        if size<5:raise ValueError("True source field-matched plot too small")
        codes=set(rng.sample(sorted(codes),size))
        control=set(rng.sample(sorted(control),size))
        rail[station['id']]=source_profile(codes,traits)
        ref[station['id']]=source_profile(control,traits)
        for p in plots:
            coord=geo[p['id']]
            physical_grid[(round(coord['lat'],4),round(coord['lon'],4))]+=1
    statistics=source_beta_stats(rail,ref)
    return rail,ref,statistics,{
       "median_source_original_railway_station_LEDA_known_richness":median(original_rich),
       "median_source_five_nontransport_plots_joint_LEDA_richness":median(local_rich),
       "n_source_original_plot_ids_same_exact_4decimal_coordinate_cell":
          sum(n-1 for n in physical_grid.values() if n>1),
       "n_true_unique_original_plot_grid_cells":len(physical_grid)}

def stratum_loo(rail,ref,stations,alloc):
    out={}
    for drop in sorted({stations[i]['stretch'] for i in alloc}):
        keepids={stations[i]['id'] for i in alloc if stations[i]['stretch']!=drop}
        if len(keepids)<8 or len({stations[i]['stretch'] for i in alloc
                if stations[i]['stretch']!=drop})<3:continue
        x=source_beta_stats({k:v for k,v in rail.items() if k in keepids},
                            {k:v for k,v in ref.items() if k in keepids})
        out[drop]=x['rail_minus_local_nontransport_CWM_beta']
    return out

def evaluate(stations,geo,traits,assignments,seed):
    result=[];loo=defaultdict(list);qa=[]
    for d in range(RAREFIED_DRAWS):
        rail,ref,stats,q=station_reference_profiles(stations,assignments,geo,traits,seed+149*d)
        result.append(stats)
        qa.append(q)
        split=stratum_loo(rail,ref,stations,assignments)
        for line,val in split.items():loo[line].append(val)
    measures=("rail_minus_local_nontransport_CWM_beta",
              "rail_minus_local_nontransport_Jaccard",
              "rail_minus_local_nontransport_2D_alpha_trace",
              "absolute_rail_CWM_distance","absolute_local_control_CWM_distance")
    return {
      "source_n_preslia_original_stations":len(assignments),
      "original_rail_stretch_counts":dict(Counter(stations[i]['stretch'] for i in assignments)),
      "source_original_5localplot_unique_ID_assignments":5*len(assignments),
      "source_rarefaction_reps_not_independent_field_data":len(result),
      "metrics":{key:mean(q[key] for q in result) for key in measures},
      "source_leave_one_preslia_stretch_out_CWM_mean":{
          k:mean(v) for k,v in loo.items()},
      "source_plot_geolocation_collision_QA":{
          key:median(v[key] for v in qa)
          for key in qa[0]}}

def run(s1,s2,sr70,headers,species,height,seed):
    rail,geo,read,source=source_station_and_reference_coords(s1,s2,sr70,headers)
    plants,plantqa=source_species(species,read)
    traits,traitqa=source_species_to_traits(rail,plants,height,seed)
    stations,elig=source_valid_local_control_plots(rail,plants,geo,traits)
    if len(stations)!=16 or len({s['stretch'] for s in stations})!=4:
        raise ValueError("Source v24.2 original 16 stations/4 original line stretches changed")
    original,optimum=solve_full_source_global_exclusivity(stations)
    if original is None or optimum['status']!='EXACT_SOURCE_BINARY_MILP_OPTIMAL' or optimum[
         'max_exact_unique_reference_original_rail_stations']!=SOURCE_MAX_SUPPORT:
        raise ValueError("v24.2 real exact maximum 14 source station support failed")
    scenarios=[];digests=set();fails=Counter()
    variants=[("v242_original_min_distance",original,{
         "status":"EXACT_SOURCE_BINARY_MILP_OPTIMAL","botany_blind_allocation_sha256":None})]
    for index in range(ALTERNATIVE_OBJECTIVE_SEEDS):
        result,q=solve_alternate_maximum(stations,seed=ORIGINAL_SEED+index*1459)
        if result is None:
            fails[q['status']]+=1
            continue
        if q['botany_blind_allocation_sha256'] in digests:
            fails['duplicate_exact_station_and_source_plot_assignment']+=1
            continue
        digests.add(q['botany_blind_allocation_sha256'])
        variants.append((f'botany_blind_source_alternative_{index:02d}',result,q))
    # Do not report a source variation analysis if all alternatives
    # point to the same source selected flora.
    if len(variants)<8:
        return {'version':'v24.3','status':'INSUFFICIENT_DISTINCT_MAX_SUPPORT_SOURCE_MATCHINGS',
            'n_distinct_optimal_source_assignments_found':len(variants),
            'n_source_alternate_solver_failures':dict(fails),
            'NO_source_robust_beta_inference':True}
    for j,(name,choice,qa) in enumerate(variants):
        e=evaluate(stations,geo,traits,choice,seed=ORIGINAL_SEED+20000*j)
        scenarios.append({'label':name,'source_n_sites':e['source_n_preslia_original_stations'],
          'source_rail_stretch_counts':e['original_rail_stretch_counts'],
          'source_original_LEDA_2D_rarefied_beta_means':e['metrics'],
          'source_leave_original_stretch_out_CWM_beta':e['source_leave_one_preslia_stretch_out_CWM_mean'],
          'source_plot_precision_and_exact_coord_collisions':e['source_plot_geolocation_collision_QA'],
          'matching_QA':qa})
    focus=('rail_minus_local_nontransport_CWM_beta',
           'rail_minus_local_nontransport_Jaccard',
           'rail_minus_local_nontransport_2D_alpha_trace')
    summaries={}
    for met in focus:
        x=[v['source_original_LEDA_2D_rarefied_beta_means'][met] for v in scenarios]
        summaries[met]={
          'n_distinct_max_support_original_source_matching_arms':len(x),
          'min_mean_over_alternative_true_source_14station_assignments':min(x),
          'median_mean_over_alternative_true_source_14station_assignments':median(x),
          'max_mean_over_alternative_true_source_14station_assignments':max(x),
          'n_source_matchings_negative':sum(a<0 for a in x),
          'n_source_matchings_positive':sum(a>0 for a in x),
          'original_source_rarefaction_draws_per_assignment':RAREFIED_DRAWS}
    loo=defaultdict(list)
    for s in scenarios:
        for stretch,v in s['source_leave_original_stretch_out_CWM_beta'].items():
            loo[stretch].append(v)
    return {
      'version':'v24.3',
      'status':'Genuine_original_railway_vs_local_vegetation_multiple_full_support_disjoint_reference_assignments',
      'source_true_railway_station_cohort_max_certified':SOURCE_MAX_SUPPORT,
      'source_original_true_16_station_4_stretch_reference_cohort':len(stations),
      'n_distinct_full_support_true_source_plot_ID_solutions':len(scenarios),
      'source_20_random_botany_blind_MILP_alternative_draws':ALTERNATIVE_OBJECTIVE_SEEDS,
      'source_solver_failure_or_duplicate_solutions':dict(fails),
      'source_Preslia_and_CzechVeg_original_species_rows_MD5':plantqa['publisher_species_csv_md5'],
      'source_all_2D_LEDA_cross_arm_species_coverage_QA':traitqa,
      'source_original_rail_censored_spatial_reference_QA':source,
      'trait_cwm_comparison_conditional_on_source_14station_assignment':summaries,
      'source_four_rail_stretch_leave_one_out_CWM_direction_in_alternative_matchings':{
          s:{'n_assignment_loo':len(v),'negative':sum(x<0 for x in v),
             'min':min(v),'median':median(v),'max':max(v)}
          for s,v in sorted(loo.items())},
      'source_multiple_source_assignment_plot_duplicate_rounded_geopoints':{
         'n_solutions_with_at_least_one_rounded_4decimal_duplicate_coordinate_cell':
           sum(v['source_plot_precision_and_exact_coord_collisions'][
               'n_source_original_plot_ids_same_exact_4decimal_coordinate_cell']>0
               for v in scenarios),
         'max_duplicate_4decimal_coord_cells':max(v['source_plot_precision_and_exact_coord_collisions'][
               'n_source_original_plot_ids_same_exact_4decimal_coordinate_cell']
               for v in scenarios)},
      'original_all_14_station_rarefied_results_across_exact_source_assignments':scenarios,
      'stop_limitations':[
        'True v24.2 14 maximum sites certified; alternate MILPs objective constraints preserve station count=14 and original published plot uniqueness',
        'Different optimal source plot assignments and possible alternative 14-station sets are SOURCE sensitivity, not independent field replicates',
        'All candidate source reference plots within 10km, but nontransport HABITAT word absence not actual verified natural/rail-free substrate',
        'Preslia station flora area huge (2020-22 author 3 habitat zones), different from five small CzechVeg 2015-23 original study plots',
        'Source station registry points from 2026 are not actual original plant sample coordinates',
        'Nearby source plots may have same approximate coordinates despite distinct source IDs; counts in QA',
        'Data are only Czech station regional beta, not a global random sample of railroads',
        'Species median height and seed mass are not local measured plant size/plasticity/evolution',
        'No freight OD, trains per 2020-22, direct seed interception, herbicide/mowing measured',
        'No regression or biological field p-values from 14 selected source stations or technical species draws']}

def main():
    p=argparse.ArgumentParser()
    for name in ('s1','s2','sr70','headers','species','height','seed','output'):
        p.add_argument('--'+name,required=True)
    a=p.parse_args()
    for file,sha in ((a.height,LEDA_SHA['height_m']),
                     (a.seed,LEDA_SHA['seed_mass_mg'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=sha:
            raise ValueError("University original public source LEDA traits SHA differs")
    height,hq=read_leda(a.height,'height_m')
    seed,sq=read_leda(a.seed,'seed_mass_mg')
    result=run(a.s1,a.s2,a.sr70,a.headers,a.species,height,seed)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    summary={'version':result['version'],'status':result['status'],
        'true_source_railway_station_cohort_max_certified':result.get('source_true_railway_station_cohort_max_certified'),
        'n_distinct_full_support_true_source_plot_ID_solutions':result.get('n_distinct_full_support_true_source_plot_ID_solutions'),
        'source_solver_failure_or_duplicate_solutions':result.get('source_solver_failure_or_duplicate_solutions'),
        'trait_cwm_comparison_conditional_on_source_14station_assignment':
            result.get('trait_cwm_comparison_conditional_on_source_14station_assignment'),
        'source_four_rail_stretch_leave_one_out_CWM_direction_in_alternative_matchings':
            result.get('source_four_rail_stretch_leave_one_out_CWM_direction_in_alternative_matchings'),
        'source_multiple_source_assignment_plot_duplicate_rounded_geopoints':
            result.get('source_multiple_source_assignment_plot_duplicate_rounded_geopoints')}
    print('V243_14SOURCE_STATIONS_ALT_EXACT_REFERENCES='+json.dumps(summary,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

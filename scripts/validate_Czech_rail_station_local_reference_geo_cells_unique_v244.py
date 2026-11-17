#!/usr/bin/env python3
"""v24.4 exact original CZECH railway stations vs local independent vegetation
reference support with NO duplicate true field coordinates, not just plot IDs.

Original Preslia 2024 39 stations field-surveyed 2020-22, verified 35
exact full-name 2026 SR70 representative railway station points.
CzechVeg 2026 117739 original headers + 2,510,340 species incidence,
source WGS84 and position precision>0<=250m, observation 2015-23,
<=10km of station and 2+ original studies across 5 source relevés per
selected station. True Oldenburg LEDA species-median canopy height
and seed mass as in v24.0-24.3.

Previous v24.2 MAX 14 stations with 70 distinct original plot IDs.
But v24.3 found even these can include coincident coordinate cells
at ~10m. v24.4 adds GLOBAL GIS-coordinate uniqueness:
 - All source unique reference relevé IDs (no source plot used twice)
 - additionally, round ORIGINAL source WGS84 latitude and longitude
   each to four decimals (~7-11m latitude/longitude), then one source
   plot per grid cell across ALL selected railway stations, including
   within one station.
 - secondly 3 decimal digits (~70-111m longitude/latitude):
   separate stricter sensitivity. These are source coordinate precision
   guardrails; they are NOT exact plot spatial non-overlap because true
   source positional precision may be as coarse as 250m.
 - Exactly 5 separate source relevés per station, >=2 true original
   distinct SOURCE publications per station, preserve >=3 original
   Preslia railway stretches, source plant traits do not enter objective.

MILP objective lexicographically maximizes retained ORIGINAL station
count, then minimizes original station-to-plot geographic distance
in each source cell-constrained optimum. HiGHS optimality certificate
required to call max support; on failure no fake ecological zero.
Report accepted station counts and original study/section coverage,
compare rarefied source 2D continuous functional centroid beta to
local nontransport vegetation only if >=10 exact stations, >=3
original railway stretches. 100 tech species rarefaction draws are
NOT confidence intervals or independent repeated field experiments.

Even when original reference *IDs* and rounded coordinates differ,
the station-wide field survey area and 5 nontransport local small
plots remain incomparable protocol/habitat-scale, and original
CzechVeg nontransport-text absence is not proof of natural habitat.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_Czech_39station_true_local_nontransport_plant_LEDA_beta_v240 import (
    LEDA_SHA,source_station_and_reference_coords,source_species_to_traits,
    source_valid_local_control_plots,source_species,source_beta_stats,quart
)
from audit_Czech_local_vegetation_exact_MILP_nonoverlap_capacity_v242 import (
    solve_full_source_global_exclusivity
)
from validate_Czech_rail_local_reference_alternative_14station_assignments_v243 import (
    station_reference_profiles,stratum_loo
)

MIN_ORIGINAL_ROUTE_STRETCHES=3
SOURCE_KNOWN_MAX_STATIONS=14
SOURCE_REQUIRED_16_STATIONS=16
SOURCE_REPS=100
SEED=20261117+244
COORD_PRECISION_ROUNDINGS=(4,3)

def approximate_coord_cell(plot,sourcegeo,decimals):
    data=sourcegeo.get(plot['id'])
    if data is None:raise ValueError('Original source WGS84 not found for verified CzechVeg relevé')
    lat,lon=data.get('lat'),data.get('lon')
    if not isinstance(lat,(int,float)) or not isinstance(lon,(int,float)):
        raise ValueError('Original source location invalid')
    if not (48<=lat<=52 and 12<=lon<=20):
        raise ValueError('Original true control coordinate outside Czech study region')
    return round(lat,decimals),round(lon,decimals)

def exact_max_stations_with_cell_uniqueness(stations,sourcegeo,
                                           decimals=4,timeout_s=120):
    import numpy as np
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    if len(stations)!=SOURCE_REQUIRED_16_STATIONS:
        raise ValueError('v24.2 source qualified 16 Preslia rail stations changed')
    n=len(stations)
    options={}
    seen_cells=defaultdict(list);seen_plots=defaultdict(list)
    for i,station in enumerate(stations):
        for p in station['control_options']:
            key=(i,p['id'])
            if key in options:raise ValueError('Duplicate source original plot ID inside one station')
            options[key]=p
    x={(i,pid):k for k,(i,pid) in enumerate(sorted(options))}
    for key,index in x.items():
        cell=approximate_coord_cell(options[key],sourcegeo,decimals)
        seen_cells[cell].append(index)
        seen_plots[key[1]].append(index)
    z={i:len(x)+i for i in range(n)}
    y={}
    for i,st in enumerate(stations):
        for study in sorted({p['study'] for p in st['control_options']}):
            y[i,study]=len(x)+n+len(y)
    stretches=sorted({s['stretch'] for s in stations})
    w={s:len(x)+n+len(y)+k for k,s in enumerate(stretches)}
    total=len(x)+n+len(y)+len(w)
    c=np.zeros(total,dtype=float)
    for key,index in x.items():
        # Source distance is an observational matching opportunity;
        # no botanical trait/identity response is used.
        c[index]=.0001*options[key]['d']/10.
    for i,index in z.items():c[index]=-100.
    rows=[];cols=[];vals=[];lb=[];ub=[]
    def cst(q,lo=-math.inf,hi=math.inf):
        r=len(lb)
        for ix,value in q:
            rows.append(r);cols.append(ix);vals.append(value)
        lb.append(lo);ub.append(hi)
    for i in range(n):
        cst([(x[i,p['id']],1) for p in stations[i]['control_options']]+[(z[i],-5)],
            lo=0,hi=0)
    # Global source plot ID disjoint, even if same original source plot
    # erroneously has geocoding in more than one source metadata record.
    for key,indices in seen_plots.items():
        cst([(ix,1) for ix in indices],hi=1)
    # Stronger source original near-coordinate cell disjointness.
    for cell,indices in seen_cells.items():
        cst([(ix,1) for ix in indices],hi=1)
    for (i,study),index in y.items():
        entries=[x[i,p['id']] for p in stations[i]['control_options'] if p['study']==study]
        cst([(ix,1) for ix in entries]+[(index,-1)],lo=0)
        cst([(ix,1) for ix in entries]+[(index,-5)],hi=0)
    for i in range(n):
        cst([(ix,1) for (j,st),ix in y.items() if j==i]+[(z[i],-2)],lo=0)
    for stretch,index in w.items():
        cst([(z[i],1) for i,s in enumerate(stations) if s['stretch']==stretch]
            +[(index,-1)],lo=0)
    cst([(ix,1) for ix in w.values()],lo=MIN_ORIGINAL_ROUTE_STRETCHES)
    A=coo_matrix((vals,(rows,cols)),shape=(len(lb),total)).tocsr()
    result=milp(c=c,integrality=np.ones(total),
        bounds=Bounds(np.zeros(total),np.ones(total)),
        constraints=LinearConstraint(A,np.asarray(lb),np.asarray(ub)),
        options={'time_limit':timeout_s,'mip_rel_gap':0.0,'presolve':True})
    if result.x is None or int(result.status)!=0:
        return None,{'status':'MIP_NOT_CERTIFIED_STOP',
            'source_original_WGS84_coord_cell_digits':decimals,
            'solver_status':int(result.status),'solver_message':str(result.message)[:500],
            'cannot_claim_maximum_Czech_station_geographic_reference_support':True}
    matched={}
    for i,iz in z.items():
        if result.x[iz]>.5:
            chosen=[options[(j,pid)] for (j,pid),index in x.items()
                    if j==i and result.x[index]>.5]
            if len(chosen)!=5 or len({p['id'] for p in chosen})!=5 or len({p['study'] for p in chosen})<2:
                raise ValueError('Original source true station matching integrity broken')
            matched[i]=chosen
    flat=[p for xs in matched.values() for p in xs]
    cellkeys=[approximate_coord_cell(p,sourcegeo,decimals) for p in flat]
    if len(flat)!=5*len(matched) or len(set(p['id'] for p in flat))!=len(flat):
        raise ValueError('Original source source plot IDs reused across railway stations')
    if len(cellkeys)!=len(set(cellkeys)):
        raise ValueError('Original source control plot coordinate rounded-cells repeated')
    if len({stations[i]['stretch'] for i in matched})<MIN_ORIGINAL_ROUTE_STRETCHES:
        raise ValueError('Source rail stretch coverage below three original railway stretches')
    if len(matched)>SOURCE_KNOWN_MAX_STATIONS:
        raise ValueError('This coordinate-constrained maximum cannot exceed v24.2 original maximum 14')
    return matched,{
        'status':'EXACT_SOURCE_MILP_GEO_UNIQUE_OPTIMAL',
        'coord_source_WGS84_lat_lon_digits':decimals,
        'max_original_railway_stations_with_five_100pct_unique_control_geo_cells':len(matched),
        'n_source_independent_releve_plot_IDs':len(flat),
        'n_source_unique_rounded_geocoordinate_cells':len(set(cellkeys)),
        'n_original_publications_in_reference_set':len({p['study'] for p in flat}),
        'source_original_preslia_rail_stretches':dict(Counter(stations[i]['stretch'] for i in matched)),
        'source_regional_GIS_precision_not_individual_field_polygon_nonoverlap':True,
        'source_mip_result_certified_optimal_status':int(result.status),
        'source_mip_optimization_gap':float(getattr(result,'mip_gap',0) or 0),
        'source_active_mip_variables':total,
        'source_active_mip_constraints':len(lb),
        'source_unique_actual_CzechVeg_control_coordinate_groups_considered':len(seen_cells),
        'source_botany_blind_distance_cost_only':True}

def source_rarefied_within_stretch_check(stations,choice,geo,traits,draws=SOURCE_REPS):
    if choice is None or len(choice)<10 or len({stations[i]['stretch'] for i in choice})<3:
        return {'identified':False,
           'reason':'Original source strict unique geographic reference cell plot coverage <10 stations or <3 stretches',
           'not_a_null_or_zero_biological_effect':True}
    rows=[];loo=defaultdict(list);coordqa=[]
    for seed in range(draws):
        rail,ref,result,qa=station_reference_profiles(
            stations,choice,geo,traits,SEED+291*seed)
        rows.append(result);coordqa.append(qa)
        subset=stratum_loo(rail,ref,stations,choice)
        for k,v in subset.items():loo[k].append(v)
    keys=('rail_minus_local_nontransport_CWM_beta',
          'rail_minus_local_nontransport_Jaccard',
          'rail_minus_local_nontransport_2D_alpha_trace',
          'absolute_rail_CWM_distance','absolute_local_control_CWM_distance')
    metrics={}
    for k in keys:
        vals=[x[k] for x in rows]
        metrics[k]={'mean':mean(vals),
           'q025_technical_NOT_independent_field_CI':quart(vals,.025),
           'q975_technical_NOT_independent_field_CI':quart(vals,.975),
           'positive_technical_draws':sum(x>0 for x in vals),
           'negative_technical_draws':sum(x<0 for x in vals)}
    return {'identified':True,
         'n_source_rail_stations':len(choice),
         'n_source_dependent_rail_station_pairs':len(choice)*(len(choice)-1)//2,
         'n_source_technical_species_draws':len(rows),
         'source_railway_minus_nearby_CzechVeg_beta_comparisons':metrics,
         'original_railway_stretch_leave_one_out_CWM_difference':{
            k:{'n_draws':len(v),'mean':mean(v),'negative_technical_draws':sum(x<0 for x in v)}
            for k,v in loo.items()},
         'source_reference_rounded_4decimal_collision_not_necessarily_0_for_stricter_3decimal':
             median(q['n_source_original_plot_ids_same_exact_4decimal_coordinate_cell'] for q in coordqa),
         'source_species_rarefaction_tech_draws_NOT_original_field_replicates':True}

def run(s1,s2,sr70,headers,species,height,seed):
    from audit_Czech_local_vegetation_exact_MILP_nonoverlap_capacity_v242 import solve_full_source_global_exclusivity
    orig,geo,read,source=source_station_and_reference_coords(s1,s2,sr70,headers)
    plants,sourceplants=source_species(species,read)
    traits,q=source_species_to_traits(orig,plants,height,seed)
    stations,elig=source_valid_local_control_plots(orig,plants,geo,traits)
    if len(stations)!=SOURCE_REQUIRED_16_STATIONS:
        raise ValueError('Original authenticated 16 railway stations with LEDA coverage changed')
    reference,prior=solve_full_source_global_exclusivity(stations)
    if reference is None or prior['max_exact_unique_reference_original_rail_stations']!=SOURCE_KNOWN_MAX_STATIONS:
        raise ValueError('v24.2 original exact source max14 failed')
    outcomes={}
    for digits in COORD_PRECISION_ROUNDINGS:
        result,qa=exact_max_stations_with_cell_uniqueness(stations,geo,digits)
        outcomes[str(digits)+'_digit_decimal_WGS84']={
          'source_global_source_station_plot_optimum_QA':qa,
          'source_LEDA_station_vs_local_reference_if_identified':
             source_rarefied_within_stretch_check(stations,result,geo,traits)
             if result is not None else {'identified':False,
                 'reason':'Coordinate unique solution not certified by exact MILP'}}
    return {'version':'v24.4',
       'status':'AUTHENTIC_CZECH_16RAILWAY_STATION_CZECHVEG_SOURCE_CELL_INDEPENDENCE_SENSITIVITY',
       'source_original_botanical_39railways_35geo_16trait_qualified':elig,
       'reference_original_no_plot_reuse_14_station_optimum':prior,
       'source_original_WGS84_reference_coord_rounding_digits':list(COORD_PRECISION_ROUNDINGS),
       'source_inference_for_each_geocoordinate_uniqueness_constraint':outcomes,
       'original_species_source_MD5':sourceplants['publisher_species_csv_md5'],
       'source_true_union_rail_and_reference_LEDA_2trait_QA':q,
       'interpretation_limits':[
        'Original source field positional precision <=250m; 4/3-decimal coordinate-cell uniqueness does NOT prove 250m spatial independence',
        'Original 2015-23 nontransport absence in HABITAT and LOCALITY does not prove physically nonrail or natural vegetation',
        'Czech Preslia whole railway station taxon lists combine 3 microhabitats over hectares, whereas reference flora pools five small CzechVeg vegetation plots',
        'Source matched 2026 national rail station point is not the actual 2020-22 botanical survey location',
        'The source Bayesian or random site/study bias is not estimated by technical rarefaction or exact optimization',
        'Source original publication titles and exact plot IDs/coordinates are not exposed in the aggregate result',
        'Train-propagule flow, station train frequency, wind/mowing pressure, railway adaptation evolution not identified']}
def main():
    p=argparse.ArgumentParser()
    for k in ('s1','s2','sr70','headers','species','height','seed','output'):
        p.add_argument('--'+k,required=True)
    a=p.parse_args()
    for file,digest in ((a.height,LEDA_SHA['height_m']),
                        (a.seed,LEDA_SHA['seed_mass_mg'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=digest:
            raise ValueError('Original Oldenburg public source LEDA 2D trait bytes changed')
    h,hqa=read_leda(a.height,'height_m')
    se,sqa=read_leda(a.seed,'seed_mass_mg')
    out=run(a.s1,a.s2,a.sr70,a.headers,a.species,h,se)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf8')
    slim={'status':out['status'],
          'original_source_exact_no_plot_ID_reuse_max':out['reference_original_no_plot_reuse_14_station_optimum']['max_exact_unique_reference_original_rail_stations'],
          'coord_uniqueness':out['source_inference_for_each_geocoordinate_uniqueness_constraint']}
    print('V244_EXACT_CZECH_CONTROL_GEO_COORD_CELL_UNIQUENESS='+json.dumps(slim,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

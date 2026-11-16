#!/usr/bin/env python3
"""v24.2: EXACT BINARY mixed-integer max-support railway station vs independent
local source plots with ZERO shared plot IDs across any two stations.

Original v24.0 source 2020–22 Czech 39 genuine Preslia rail station floras
(35 official SR70 effective-2026 coordinates), CzechVeg 2026 original
2015–23 vascular source relevés within 10 km, location precision>0 <=250m,
true source >=5 jointly trait-known species / plot and >=55% coverage,
original botanist publications >=2 per station, source original LEDA
height+seed species medians.

The previous heuristic v24.1 found ZERO 16-station disjoint allocations
in 299 batches of 120 random attempts. That is a technical design failure,
not proof global matching is mathematically impossible. Therefore solve
the exact finite binary mixed-integer program:

  x[i,p] binary choose original CzechVeg plot p for rail station i.
  z[i] binary retain authentic rail station i.
  y[i,study] binary a published study contributes to retained station i.
  w[stretch] binary an original Preslia rail stretch is represented.
  sum_p x[i,p] = 5*z[i]
  sum_i x[i,p] <= 1 for every source original plot ID
  y[i,study] <= sum_{p in study} x[i,p] <= 5*y[i,study]
  sum_study y[i,study] >= 2*z[i]
  w[stretch] <= sum_{stations in stretch} z[i]
  sum_stretch w[stretch] >= 3.
  Maximize sum_i z[i] (tiny positive distance costs break ties, only
  original geographic opportunity data and source study; no trait effects).

There are only 16 source eligible railway stations over 4 stretches.
Return *optimal* supported number and original-stretch group counts.
If <10 sites/3 stretches, report no external beta test identified,
never a zero biology effect. Otherwise, conditional on ONE exact
source-optimized five-plot physical matching, evaluate rarefied
two-trait CWM beta vs original rail station source and Jaccard using
199 technical species draws; NOT independent sampling uncertainty.

Caveats: source station whole-yard floras and CzechVeg 5 small plots
remain different botanical grain, source study site habitat and
true rail/station sample locations are not equal even after matching.
2026 station representative point is not 2020–22 botanical plot GPS.
No independent natural habitat proof, moving seed or genetic adaptation.
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

MIN_STRETCHES=3
MIN_SITES_FOR_BETA=10
REPS=199
SEED=20261117+242

def solve_full_source_global_exclusivity(stations,min_stretches=MIN_STRETCHES):
    import numpy as np
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    if len(stations)==0:raise ValueError("No source vegetation stations")
    # Each station had already passed >=5 original source plots with
    # source species LEDA coverage and two distinct ORIGINAL studies.
    ids=sorted({row['id'] for st in stations for row in st['control_options']})
    idplot={id:j for j,id in enumerate(ids)}
    idmap={i:{p['id']:p for p in st['control_options']} for i,st in enumerate(stations)}
    stretches=sorted({s['stretch'] for s in stations})
    if len(stretches)<min_stretches:raise ValueError("Insufficient original source rail stretches")
    x={}
    k=0
    for i in range(len(stations)):
        for plot in sorted(idmap[i]):
            x[i,plot]=k;k+=1
    z={i:k+i for i in range(len(stations))}
    k+=len(stations)
    y={}
    for i in range(len(stations)):
        for study in sorted({p['study'] for p in idmap[i].values()}):
            y[i,study]=k;k+=1
    w={s:k+i for i,s in enumerate(stretches)}
    k+=len(stretches)
    c=np.zeros(k,dtype=float)
    for (i,plot),index in x.items():
        # cost depends only on genuine original geographic DISTANCE,
        # not on botanical taxa/trait outcome. Source station matches
        # maximized even when they are far from generic town areas.
        c[index]=1e-5*idmap[i][plot]['d']/10
    for i,index in z.items():c[index]=-100.
    rows=[];cols=[];data=[];lower=[];upper=[]
    def constrain(entries,lo=-math.inf,hi=math.inf):
        row=len(lower)
        for ix,coef in entries:
            rows.append(row);cols.append(ix);data.append(coef)
        lower.append(lo);upper.append(hi)
    # 5 original controls per retained railway station
    for i in range(len(stations)):
        constrain([(ix,1.) for (j,p),ix in x.items() if j==i]+[(z[i],-5.)],
                  lo=0.,hi=0.)
    # one original CzechVeg relevé can be used once globally
    for plot in ids:
        constrain([(index,1.) for (i,p),index in x.items() if p==plot],
                  hi=1.)
    # Study y[i,paper] iff at least one original field plot from paper
    for (i,study),ix in y.items():
        source=[x[i,plot] for plot,p in idmap[i].items() if p['study']==study]
        constrain([(j,1.) for j in source]+[(ix,-1.)],lo=0.)
        constrain([(j,1.) for j in source]+[(ix,-5.)],hi=0.)
    # each retained station must have >=2 genuinely distinct source papers
    for i in range(len(stations)):
        constrain([(ix,1.) for (j,p),ix in y.items() if j==i]+[(z[i],-2.)],lo=0.)
    for st,ix in w.items():
        constrain([(z[i],1.) for i,s in enumerate(stations) if s['stretch']==st]+[(ix,-1.)],lo=0.)
    constrain([(ix,1.) for ix in w.values()],lo=float(min_stretches))
    mat=coo_matrix((data,(rows,cols)),shape=(len(lower),k)).tocsr()
    result=milp(c=c,integrality=np.ones(k),bounds=Bounds(np.zeros(k),np.ones(k)),
                constraints=LinearConstraint(mat,np.asarray(lower),np.asarray(upper)),
                options={'time_limit':105,'mip_rel_gap':0.0,'presolve':True})
    if result.x is None or result.status!=0:
        return None,{'status':'EXACT_SOLVER_NOT_OPTIMAL_NO_MAX_SUPPORT_CLAIM',
           'solver_status':int(result.status),'solver_message':str(result.message),
           'n_original_source_stations':len(stations),'n_original_control_plot_candidates':len(ids)}
    active=sorted(i for i,j in z.items() if result.x[j]>.5)
    allocated={}
    for i in active:
        chosen=[idmap[i][plot] for (j,plot),index in x.items()
                if j==i and result.x[index]>.5]
        if len(chosen)!=5 or len({p['id'] for p in chosen})!=5 or len({p['study'] for p in chosen})<2:
            raise ValueError('Source exact MILP solution violated five unique original plots or two papers')
        allocated[i]=chosen
    allplots=[p['id'] for choices in allocated.values() for p in choices]
    if len(set(allplots))!=len(allplots):
        raise ValueError('Exact MILP source final plot allocation reused original plot ID')
    if len({stations[i]['stretch'] for i in active})<min_stretches:
        raise ValueError('Exact source solver violated 3 original rail stretches')
    return allocated,{
       'status':'EXACT_SOURCE_BINARY_MILP_OPTIMAL',
       'original_source_available_rail_stations':len(stations),
       'max_exact_unique_reference_original_rail_stations':len(active),
       'rail_stretches_with_optimal_supported_stations':
          dict(Counter(stations[i]['stretch'] for i in active)),
       'source_total_reference_plots_in_optimal_assignment':len(allplots),
       'source_total_distinct_reference_plot_IDs':len(set(allplots)),
       'source_original_reference_studies_at_least_two_per_rail_station':True,
       'source_paired_plot_study_matching_is_STRICT_not_study_year_exposure':True,
       'n_unique_original_eligible_control_plot_ids':len(ids),
       'n_source_binary_allocation_candidate_edges':len(x),
       'n_true_original_publication_station_study_variables':len(y),
       'n_original_rail_stretch_variables':len(w),
       'source_ilp_constraints':len(lower),
       'source_exact_mip_solver_fun':float(result.fun),
       'source_exact_mip_solver_gap':float(getattr(result,'mip_gap',0) or 0),
       'source_no_plant_trait_effect_used_in_station_or_reference_plot_selection':True}

def evaluate_matched(stations,selected,traits,reps=REPS):
    active={stations[i]['id']:stations[i] for i in selected}
    ref={stations[i]['id']:set().union(*(p['taxa'] for p in controls))
         for i,controls in selected.items()}
    if len(active)<MIN_SITES_FOR_BETA or len({x['stretch'] for x in active.values()})<MIN_STRETCHES:
        return {'identified':False,
           'reason':'Fewer than 10 true railway station and 3 source rail stretch support after exact disjoint source controls',
           'original_source_rail_station_n':len(active),
           'not_an_estimated_zero_biological_effect':True}
    rng=random.Random(SEED)
    out=[]
    for r in range(reps):
        rail_profile={};local_profile={}
        for id,s in active.items():
            rail_taxa=s['rail_source_trait_taxa']
            source_taxa=ref[id]
            size=min(len(rail_taxa),len(source_taxa))
            if size<5:raise ValueError('True source station local reference has no valid species')
            a=set(rng.sample(sorted(rail_taxa),size))
            b=set(rng.sample(sorted(source_taxa),size))
            rail_profile[id]=source_profile(a,traits)
            local_profile[id]=source_profile(b,traits)
        out.append(source_beta_stats(rail_profile,local_profile))
    metrics={}
    for k in ('rail_minus_local_nontransport_CWM_beta',
              'rail_minus_local_nontransport_Jaccard',
              'rail_minus_local_nontransport_2D_alpha_trace',
              'absolute_rail_CWM_distance',
              'absolute_local_control_CWM_distance'):
        v=[x[k] for x in out]
        metrics[k]={'source_true_fixed_optimal_plot_assignment_species_rarefied_mean':mean(v),
              'q025_technical_NOT_independent_field_CI':quart(v,.025),
              'q975_technical_NOT_independent_field_CI':quart(v,.975),
              'tech_draws':len(v),'fraction_negative_technical':sum(x<0 for x in v)/len(v)}
    return {'identified':True,'source_rail_stations':len(active),
      'source_national_original_rail_stretches':len({x['stretch'] for x in active.values()}),
      'source_unique_control_original_releve_count':5*len(active),
      'source_technical_species_rarefaction_draws':len(out),
      'not_resampled_original_control_plot_assignments_only_species_rarefaction':True,
      'same_original_stations_both_rail_and_local_reference':True,
      'source_original_rail_and_local_LEDA_known_species_rarefied_to_min_per_station':True,
      'metrics':metrics}

def run(s1,s2,sr70,headers,species,height,seed):
    src,geo,ids,source=source_station_and_reference_coords(
        s1,s2,sr70,headers)
    plants,plantqa=source_species(species,ids)
    traits,traitqa=source_species_to_traits(src,plants,height,seed)
    stations,sampleqa=source_valid_local_control_plots(
        src,plants,geo,traits)
    if len(stations)!=16 or len({s['stretch'] for s in stations})!=4:
        raise ValueError('v24.0 original field-certified 16 station 4 stretch cohort changed')
    choice,optim=solve_full_source_global_exclusivity(stations)
    return {'version':'v24.2',
       'status':'ORIGINAL_SOURCE_BINARY_MILP_MAX_FEASIBLE_DISJOINT_RAIL_LOCAL_PLOT_SUPPORT',
       'source_full_2020_22_Preslia_rail_station_2015_23_CzechVeg_reference_metadata':source,
       'source_original_2510340_plant_rows_checked':plantqa,
       'source_original_LEDA_dual_trait_cross_arm_QA':traitqa,
       'source_16station_trait_qualified_references_before_full_independence':sampleqa,
       'original_exact_source_global_max_unique_reference':optim,
       'fixed_optimal_source_reference_assignment_exploratory_species_beta':
           evaluate_matched(stations,choice,traits) if choice is not None else {
             'identified':False,'reason':'Original source exact MILP solver failed to establish optimum',
             'not_an_estimated_zero_biological_effect':True},
       'no_source_plot_ids_or_precise_field_coordinates_published':True,
       'interpretation_limits':[
         'Exact binary solver establishes max feasible source original station sample given source five plots two original studies and no plot reuse',
         'One optimal source plot allocation chosen using only physical distance and study metadata, not taxa CWM response',
         'Population-level uncertainty in alternative plot allocations is NOT evaluated, so tech species-rarefaction quantiles are not independent field CI',
         'Railway station flora spans original large three-habitat survey, local reference is union of five heterogeneous small plots',
         'Different historical field years and 2026 station centroid vs 2020–22 original source botanic sampling locations remain unmatched',
         'No proof of natural land habitat from absence of transport words, no equal soil climate mowing or geographic source vegetation class',
         'Original Czech railway 16 station field research and CzechVeg vegetation publication are independent surveys but not randomized rail treatment',
         'No train-seed direct transport, station passenger/freight timetables or adaptive genetics measured']}
def main():
    p=argparse.ArgumentParser()
    for name in ('s1','s2','sr70','headers','species','height','seed','output'):
        p.add_argument('--'+name,required=True)
    a=p.parse_args()
    for file,hash in ((a.height,LEDA_SHA['height_m']),(a.seed,LEDA_SHA['seed_mass_mg'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=hash:
            raise ValueError('Original Oldenburg LEDA public trait source SHA differs')
    h,hqa=read_leda(a.height,'height_m')
    s,sqa=read_leda(a.seed,'seed_mass_mg')
    result=run(a.s1,a.s2,a.sr70,a.headers,a.species,h,s)
    pth=Path(a.output);pth.parent.mkdir(parents=True,exist_ok=True)
    pth.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf8')
    summary={'status':result['status'],
       'exact_source_optimum':result['original_exact_source_global_max_unique_reference'],
       'field_species_reference_comparison':result['fixed_optimal_source_reference_assignment_exploratory_species_beta']}
    print('V242_EXACT_SOURCE_MAX_NONOVERLAP_RAIL_REFERENCE='+json.dumps(summary,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

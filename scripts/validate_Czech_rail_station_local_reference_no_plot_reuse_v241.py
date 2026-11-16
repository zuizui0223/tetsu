#!/usr/bin/env python3
"""v24.1: do true Czech original railway station vs independent local
nontransport vegetation 2D functional beta contrasts persist when NO
CzechVeg reference plot is reused across any two railway stations?

Original source join is v24.0 with STRICT 10km, source vegetation year
2015–2023, source coordinate PRECISION>0 <=250m, 5 original source
control plots from >=2 independently published studies, actual Preslia
2020–22 39 railway stations, 2026 derived SR70 station representative
locations (35/39), 2026 CzechVeg 117739 metadata and 2510340 plant rows,
same Oldenburg LEDA source height+seed medians all SHA/MD5 gated.

Critical v24.0 vulnerability: on 299 technical 16-station reference
draws, controls were repeated 8–18 times across stations, violating
reference pseudo-independence. All v24.0 conclusions are descriptive,
and source 2020 railway station FLORA vs pooled five different small
nonrail field plots still differs strongly in area and botanic sampling.

This v24.1 tests stricter GLOBAL MATCHING: exactly 5 different true
CzechVeg original field plots per original railway station, source
original publication count >=2, and **one reference plot belongs to no
more than one railway station in each draw**. A real bipartite augmenting
matching algorithm assigns scarce source plots after reserving two
source publication anchors for each station. If no such global matching
can be constructed after 120 attempts, declare NOT IDENTIFIED, never
allow hidden plot reuse or fake independent reference observations.

Same source 16 stations/4 original railway stretches; 299 technical
resamplings of original plot opportunities, original species matched to
LEDA and source TAXON richness-rarefaction per station in both arms.
No actual natural habitat verified, local field plot area different,
historical 2026 station point NOT 2020–22 field GPS. Technical draws
are neither field repetitions nor causal selection experiments.
"""
from __future__ import annotations
import argparse,hashlib,json,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_Czech_39station_true_local_nontransport_plant_LEDA_beta_v240 import (
    LEDA_SHA,DRAWS,source_station_and_reference_coords,source_species_to_traits,
    source_valid_local_control_plots,source_profile,source_beta_stats,
    source_species,sample_local_plots,quart
)

SEED=20261117+241
MAX_MATCH_TRIES=120
MIN_STATIONS=10
MIN_RAIL_STRETCHES=3

def globally_disjoint_controls(stations,rng,max_attempts=MAX_MATCH_TRIES):
    """Two independent source publication anchors/station first, then
    augmenting bipartite matching to fill 3 additional original
    field-plot slots; never repeat one source relevé across stations.
    """
    options=[]
    for x in stations:
        byid={p['id']:p for p in x['control_options']}
        if len(byid)<5 or len({p['study'] for p in byid.values()})<2:
            raise ValueError('Source original local station failed five distinct plots/two study constraints')
        options.append(byid)
    n=len(stations)
    multiplicity=Counter(pid for item in options for pid in item)
    for iteration in range(max_attempts):
        used=set();anchored=defaultdict(list)
        ordering=list(range(n))
        rng.shuffle(ordering)
        ordering.sort(key=lambda i:(len(options[i]),len({v['study'] for v in options[i].values()})))
        failed=False
        for i in ordering:
            st=defaultdict(list)
            for id,p in options[i].items():
                if id not in used:st[p['study']].append(id)
            if len(st)<2:
                failed=True;break
            studies=list(st)
            rng.shuffle(studies)
            studies.sort(key=lambda study:len(st[study]))
            # Vary the first two paper sources over technical draws,
            # prioritize scarce legal source opportunities when present.
            first,second=studies[:2]
            for paper in (first,second):
                picks=list(st[paper])
                rng.shuffle(picks)
                picks.sort(key=lambda pid:multiplicity[pid])
                chosen=picks[0]
                used.add(chosen);anchored[i].append(chosen)
        if failed:continue
        # Each station still needs THREE more ORIGINAL source plots, and
        # no plot may be reused by another rail station's local reference.
        tasks=[(i,k) for i in range(n) for k in range(3)]
        rng.shuffle(tasks)
        tasks.sort(key=lambda slot:sum(pid not in used for pid in options[slot[0]]))
        match={}
        def augment(slot,seen):
            i,_=slot
            choices=[pid for pid in options[i] if pid not in used]
            rng.shuffle(choices)
            choices.sort(key=lambda pid:multiplicity[pid])
            for pid in choices:
                if pid in seen:continue
                seen.add(pid)
                old=match.get(pid)
                if old is None or augment(old,seen):
                    match[pid]=slot
                    return True
            return False
        if not all(augment(slot,set()) for slot in tasks):
            continue
        result=defaultdict(list)
        for i,anchors in anchored.items():
            result[i].extend(options[i][pid] for pid in anchors)
        for pid,slot in match.items():
            i,_=slot
            result[i].append(options[i][pid])
        if len(result)!=n or any(len(result[i])!=5 for i in range(n)):
            raise AssertionError('Source five reference plot slots not all assigned uniquely')
        flat=[p['id'] for vals in result.values() for p in vals]
        if len(flat)!=5*n or len(set(flat))!=5*n:
            raise AssertionError('CzechVeg physical source plot reuse forbidden across all railway station references')
        if any(len({p['study'] for p in result[i]})<2 for i in range(n)):
            raise AssertionError('Source two independent original publications per station missing')
        return dict(result),iteration+1
    return None,max_attempts

def no_reuse_draw(stations,trait,rng,rarefied):
    alloc,attempts=globally_disjoint_controls(stations,rng)
    if alloc is None:return None,attempts
    rail={};ref={}
    for i,s in enumerate(stations):
        plots=alloc[i]
        x=set().union(*(p['taxa'] for p in plots))
        y=s['rail_source_trait_taxa']
        if min(len(x),len(y))<5:raise ValueError('True original source control/reference species missing')
        if rarefied:
            k=min(len(x),len(y))
            x=set(rng.sample(sorted(x),k))
            y=set(rng.sample(sorted(y),k))
        ref[s['id']]=source_profile(x,trait)
        rail[s['id']]=source_profile(y,trait)
    val=source_beta_stats(rail,ref)
    val['original_reference_plot_ids_across_all_railway_stations_reused_n']=0
    val['original_source_matching_attempts']=attempts
    return val,attempts

def evaluate(stations,trait,rarefied,seed,draws=DRAWS):
    if len(stations)<MIN_STATIONS or len({s['stretch'] for s in stations})<MIN_RAIL_STRETCHES:
        return {'identified':False,'reason':'Original source <10 stations or <3 stretches after true functional trait support'}
    rng=random.Random(seed)
    draws_data=[];fails=0
    for z in range(draws):
        value,trycount=no_reuse_draw(stations,trait,rng,rarefied)
        if value is None:
            fails+=1
            continue
        draws_data.append(value)
    if len(draws_data)<max(30,int(.8*draws)):
        return {'identified':False,'reason':'Source globally disjoint five-control plot matching does not support at least 80% technical replicates',
                'n_successful_source_disjoint_matching_attempts':len(draws_data),
                'n_source_disjoint_matching_failures':fails,
                'underpowered_plot_overlap_not_a_biological_zero':True}
    metrics={}
    for k in ('rail_minus_local_nontransport_CWM_beta',
              'rail_minus_local_nontransport_Jaccard',
              'rail_minus_local_nontransport_2D_alpha_trace',
              'absolute_rail_CWM_distance',
              'absolute_local_control_CWM_distance'):
        vals=[v[k] for v in draws_data]
        metrics[k]={'mean':mean(vals),'q025_technical_NOT_field_CI':quart(vals,.025),
           'q975_technical_NOT_field_CI':quart(vals,.975),'n_source_technical_draws':len(vals),
           'technical_draws_direction_negative':sum(x<0 for x in vals)/len(vals)}
    return {'identified':True,
       'true_source_original_railway_stations':len(stations),
       'true_source_rail_stretches':len({s['stretch'] for s in stations}),
       'original_control_plots_per_rail_station':5,
       'n_original_control_plots_unique_per_draw':5*len(stations),
       'n_source_matching_failures':fails,
       'source_matching_tries_median':median(v['original_source_matching_attempts'] for v in draws_data),
       'n_source_technical_resamples_NOT_field_experiments':len(draws_data),
       'source_plots_never_reused_between_original_rail_stations':True,
       'source_reference_publications_per_station_at_least_two':True,
       'species_richness_rarefaction_same_size_per_actual_station':rarefied,
       'metric_estimates':metrics}

def run(s1,s2,sr70,headers,species,h,seed,draws=DRAWS):
    original,source_ref_geo,read,sourceqa=source_station_and_reference_coords(
        s1,s2,sr70,headers)
    plants,plantqa=source_species(species,read)
    traits,traitqa=source_species_to_traits(original,plants,h,seed)
    stations,eligible=source_valid_local_control_plots(
       original,plants,source_ref_geo,traits)
    if len(stations)!=16 or len({s['stretch'] for s in stations})!=4:
        raise ValueError('Authenticated corrected v24.0 original 16 stations / 4 stretches source cohort changed')
    outcome={
      'version':'v24.1','status':'SOURCE_AUTHENTICATED_NO_DUPLICATE_LOCAL_REFERENCE_PLOTS',
      'same_original_2020_22_railway_sites_2015_23_nontransport_reference_source_qa':sourceqa,
      'genuine_original_2510340_species_rows_qa':plantqa,
      'genuine_both_study_plant_trait_comparability':traitqa,
      'source_rail_station_field_vs_5_local_CzechVeg_plots_qa':eligible,
      'primary_full_source_16stations_4stretches_global_unique_reference':{
         'as_original_collected_species_richness':evaluate(stations,traits,False,SEED,draws),
         'station_and_control_known_trait_species_equalized':evaluate(stations,traits,True,SEED+171,draws)},
      'source_overlap_and_study_limits':[
        'True 2015–23 independent CzechVeg original nontransport source plots are disjoint across all railway stations in each technical draw',
        'Preslia source railway 2020–22 station-wide three managed habitats and 2–90ha different botanical area than small vegetation plot controls',
        '2026 SR70 matched registry point is NOT exact Preslia original field plot GPS; radius 10km',
        'Nontransport-text source background is not proof that controls are nonrail, natural or same gravel/soil/substrate',
        'Control plots may have unequal representation of vegetation types, plant guilds and landscape mosaics',
        'CzechVeg original plot compilation may have overlapping sites despite different recorded plot IDs, no hard field-coordinate uniqueness verification',
        '299 plot assignments are technical replications of same fixed published field vegetation, NOT independent experiments',
        'Source original railway station 16 across four Preslia stretches not 16 independent regional ecosystems',
        'Railway train frequency, freight source/OD, mowing, plasticity and hereditary railway selection NOT measured']}
    return outcome

def main():
    p=argparse.ArgumentParser()
    for name in ('s1','s2','sr70','headers','species','height','seed','output'):
        p.add_argument('--'+name,required=True)
    p.add_argument('--draws',type=int,default=DRAWS)
    a=p.parse_args()
    for f,sha in ((a.height,LEDA_SHA['height_m']),(a.seed,LEDA_SHA['seed_mass_mg'])):
        if hashlib.sha256(Path(f).read_bytes()).hexdigest()!=sha:
            raise ValueError('True 2026 Oldenburg LEDA original SHA differs')
    height,hq=read_leda(a.height,'height_m')
    seedmass,sq=read_leda(a.seed,'seed_mass_mg')
    result=run(a.s1,a.s2,a.sr70,a.headers,a.species,height,seedmass,a.draws)
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    summary={'status':result['status'],'source_true_reference_support':result['source_rail_station_field_vs_5_local_CzechVeg_plots_qa'],
       'joined_taxa':result['genuine_both_study_plant_trait_comparability'],
       'results':result['primary_full_source_16stations_4stretches_global_unique_reference']}
    print('V241_DISJOINT_TRUE_LOCAL_PLOTS='+json.dumps(summary,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()

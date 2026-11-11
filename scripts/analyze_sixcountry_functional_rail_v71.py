#!/usr/bin/env python3
"""Six-country GRIIS source flora × real rail/road paths × source-vetted GIFT growth form.

Primary result is only an observational trait-stratified SOURCE correlation.
True cutting, wind, active freight, adaptation and dispersal are unmeasured.
Group comparisons use identical geographic cell pairs, never a favorable
separate cohort. Pair count is NOT the number of independent replicates.
"""
from __future__ import annotations
import argparse, json, math
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from statistics import median

from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,PERIODS,STRATEGIES,combine_source_regions,MIN_SOURCE_CELLS,MIN_VALID_PAIRS
from analyze_macro_50km_rail_road_aliens_v34 import (
    _eligible_both_period_cells,_routes,partial_rank_r
)
from analyze_macro_rail_assemblages_v22 import incidence_probabilities,expected_rarefied_overlap

GROUPS=('herbaceous','woody')
MIN_GROUP_TAXA=3
MIN_GROUP_RECORDS=20
MIN_TRAIT_COUNTRY_SOURCE_COVERAGE=.35


def validate_traits(t):
    if (t.get('complete') is not True or t.get('source',{}).get('GIFT_meta_ID')!='1.2.1'
        or t.get('source',{}).get('GIFT_biasref') is not False
        or t.get('source',{}).get('GIFT_biasderiv') is not False
        or t.get('source',{}).get('GIFT_categorical_min_agreement',0)<.66):
        raise ValueError('Strict independently sourced GIFT 1.2.1 trait provenance missing')
    if not set(BLOCKS_SIX).issubset(t.get('by_block',{})):
        raise ValueError('Six-country trait source coverage missing')
    vals=t.get('growth_form_by_source_gbif_species_key')
    if not isinstance(vals,dict) or any(not str(k).isdigit() or v not in GROUPS for k,v in vals.items()):
        raise ValueError('Invalid trait accepted-GBIF species key and known category mapping')
    return vals


def eligible_group_counts(p,traits,group):
    values={k:v for k,v in p['counts'].items() if traits.get(k)==group}
    if len(values)<MIN_GROUP_TAXA or sum(values.values())<MIN_GROUP_RECORDS:
        return None
    return values


def make_pair_rows(block,period,strategy,cells,rail,road,traits):
    ar=_routes(rail,block,strategy,'pair_routes')
    rd=_routes(road,block,strategy,'pair_paths')
    if set(ar)!=set(rd):raise ValueError('RAIL and ROAD source geographic pair universes disagree')
    group_cells={g:{} for g in GROUPS}
    for cid,bytime in cells.items():
        for g in GROUPS:
            counts=eligible_group_counts(bytime[period],traits,g)
            if counts:
                group_cells[g][cid]={
                    'counts':counts,
                    'probs':incidence_probabilities(counts,k=MIN_GROUP_RECORDS),
                }
    rows=[]
    for a,b in combinations(sorted(cells),2):
        p=(a,b);r=ar[p];d=rd[p]
        rratio=r.get('rail_to_geodesic_ratio');dratio=d.get('road_to_geodesic_ratio')
        rp=r.get('cell_centroid_snap_augmented_rail_path_km');dp=d.get('centroid_snap_augmented_road_path_km')
        geo=r.get('source_cell_centroid_geodesic_km')
        if not all(isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) and x>0
                   for x in (rratio,dratio,rp,dp,geo)):
            continue
        if rratio<.95 or dratio<.95: continue
        n1=cells[a][period]['all_plant_indexed_records']
        n2=cells[b][period]['all_plant_indexed_records']
        if min(n1,n2)<=0:raise ValueError('Missing index-record source opportunity covariates')
        x={'a':a,'b':b,'x':math.log(rp/dp),
           'covars':[math.log(geo),math.log(dp),abs(math.log(n1/n2)),(math.log(n1)+math.log(n2))/2]}
        for g in GROUPS:
            pa=group_cells[g].get(a);pb=group_cells[g].get(b)
            x[g]=expected_rarefied_overlap(pa['probs'],pb['probs']) if pa and pb else None
        rows.append(x)
    return rows,group_cells


def rho_for_rows(rows,key):
    recs=[r for r in rows if isinstance(r.get(key),(float,int)) and math.isfinite(r[key])]
    n_cells=len({c for r in recs for c in (r['a'],r['b'])})
    if n_cells<MIN_SOURCE_CELLS or len(recs)<MIN_VALID_PAIRS:
        return {'partial_rho':None,'n_source_pairs_sharing_cells':len(recs),
                'n_distinct_50km_cells':n_cells,'source_support':False}
    coef=partial_rank_r([r['x'] for r in recs],[r[key] for r in recs], [r['covars'] for r in recs])
    return {'partial_rho':coef,'n_source_pairs_sharing_cells':len(recs),
            'n_distinct_50km_cells':n_cells,'source_support':coef is not None}


def analyze_one(block,period,strategy,cells,rail,road,traits,block_trait_coverage):
    rows,gc=make_pair_rows(block,period,strategy,cells,rail,road,traits)
    group_report={g:rho_for_rows(rows,g) for g in GROUPS}
    common=[r for r in rows if all(r.get(g) is not None for g in GROUPS)]
    common_report={g:rho_for_rows(common,g) for g in GROUPS}
    # With same exact physical pairs, DIFFERENCE in signed descriptive rhos
    # is a sensitivity contrast, NOT a causal interaction coefficient/p-value.
    rho=[common_report[g]['partial_rho'] for g in GROUPS]
    delta=round(rho[0]-rho[1],6) if all(x is not None for x in rho) else None
    if block_trait_coverage is None or block_trait_coverage<MIN_TRAIT_COUNTRY_SOURCE_COVERAGE:
        # Do not selectively discard low-coverage regions while claiming
        # country-level trait differences, retain diagnostics but block inference.
        for g in GROUPS:
            group_report[g]['partial_rho']=None
            group_report[g]['source_support']=False
            common_report[g]['partial_rho']=None
            common_report[g]['source_support']=False
        delta=None
    return {'n_original_source_cells_both_periods':len(cells),
        'source_country_trait_coverage_fraction':block_trait_coverage,
        'coverage_minimum_precommitted':MIN_TRAIT_COUNTRY_SOURCE_COVERAGE,
        'growth_group_eligible_cell_counts':{g:len(gc[g]) for g in GROUPS},
        'n_all_original_source_plant_railroad_pairs':len(rows),
        'trait_group_rarefied20_source_similarity':group_report,
        'same_cell_pair_cohort_trait_sensitivity':{
            'n_same_cell_pair_rows':len(common),
            'herbaceous':common_report['herbaceous'],
            'woody':common_report['woody'],
            'descriptive_delta_rho_herbaceous_minus_woody':delta},
        'NOT_measured_wind_cutting_or_evolution':True}


def country_equal_summary(regions,period,strategy):
    results={}
    for g in GROUPS:
        rhos={k:rec['periods'][period][strategy]['trait_group_rarefied20_source_similarity'][g]['partial_rho']
            for k,rec in regions.items()}
        x={k:v for k,v in rhos.items() if v is not None}
        results[g]={
            'country_partial_rhos':x,
            'n_independent_country_blocks_supported':len(x),
            'country_equal_median_partial_rho':round(median(x.values()),6) if x else None,
            'n_negative_rail_advantage_hypothesis_compatible':sum(v<0 for v in x.values()),
            'n_positive_contradictory':sum(v>0 for v in x.values()),
            'five_or_more_countries_with_support':len(x)>=5,
            'leave_country_out_medians':{c:round(median(v for k,v in x.items() if k!=c),6) for c in x if len(x)>1},
        }
    delta={k:r['periods'][period][strategy]['same_cell_pair_cohort_trait_sensitivity']['descriptive_delta_rho_herbaceous_minus_woody']
           for k,r in regions.items()}
    d={k:v for k,v in delta.items() if v is not None}
    results['SAME_COHORT_GROUP_DIFFERENCE']={
        'country_descriptive_herb_vs_woody_rho_differences':d,
        'n_independent_country_blocks_supported':len(d),
        'country_equal_median_difference':round(median(d.values()),6) if d else None,
        'no_causal_stress_or_seeding_inference':True}
    return results


def run(flora,rail,road,trait_src):
    traits=validate_traits(trait_src)
    if any(x.get('complete') is not True for x in (flora,rail,road)):
        raise ValueError('Original source artifact incomplete')
    if set(rail.get('regions',{}))!=set(BLOCKS_SIX) or set(road.get('regions',{}))!=set(BLOCKS_SIX):
        raise ValueError('Rail and road macro comparison must use exact frozen six blocks')
    out={'complete':True,'utc':datetime.now(timezone.utc).isoformat(),
        'study':'Source-vetted GIFT growth form × country GRIIS indexed flora × physical rail/road shortest paths',
        'source_growth_form_trait_id':'1.2.1',
        'no_interpretation_as_measured_mowing_or_wind_tolerance':True,
        'NO_train_mediated_seed_disperal_mechanism_established':True,
        'NO_plant_local_adaptation_or_plasticity_established':True,
        'six_original_independent_country_blocks':list(BLOCKS_SIX),
        'regions':{},'country_equal_summary':{}}
    for block in BLOCKS_SIX:
        cells=_eligible_both_period_cells(block,flora)
        v=trait_src['by_block'].get(block,{})
        cov=v.get('source_GIFT_coverage_species_fraction')
        if cov is not None and (not isinstance(cov,(int,float)) or cov<0 or cov>1):
            raise ValueError('Invalid source trait coverage fraction')
        out['regions'][block]={'n_qualified_two_period_source_cells':len(cells),'periods':{}}
        for period in PERIODS:
            out['regions'][block]['periods'][period]={}
            for strat in STRATEGIES:
                out['regions'][block]['periods'][period][strat]=analyze_one(
                    block,period,strat,cells,rail,road,traits,cov)
    for period in PERIODS:
        out['country_equal_summary'][period]={}
        for strat in STRATEGIES:
            out['country_equal_summary'][period][strat]=country_equal_summary(out['regions'],period,strat)
    out['inference_limitations']=[
        'GIFT growth form is species-level baseline morphology, NOT independently measured tolerance to railway cutting or slipstream wind.',
        '2026 mapped rail routes do not prove operated/freight train service during 2010–2025 species records.',
        'GRIIS-positive GBIF facets do not establish biological absences or seed carriage; count rarefaction does not remove dataset bias.',
        'Same 50km cells recur in pairwise comparisons; only country blocks are independent descriptive summaries.',
        'Trait availability is nonrandom. Low coverage (<35% of recorded species) suppresses ecological coefficients.',
        'Herbaceous versus woody is NOT a direct annuality/regrowth height/seed aerodynamic trait.',
        'Rho contrasts use exact same geographic pairs, not a statistical stress-by-trait interaction.',
        'No measured cutting frequency, gust speed, active train traffic, or field genotype/common-garden evidence available.'
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--flora',required=True)
    p.add_argument('--growth-form',required=True)
    p.add_argument('--rail-sources',required=True,nargs='+')
    p.add_argument('--road-sources',required=True,nargs='+')
    p.add_argument('--output',required=True)
    x=p.parse_args()
    load=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
    flora=load(x.flora); trait=load(x.growth_form)
    rail=combine_source_regions([load(p) for p in x.rail_sources]);road=combine_source_regions([load(p) for p in x.road_sources])
    result=run(flora,rail,road,trait)
    d=Path(x.output);d.parent.mkdir(parents=True,exist_ok=True)
    d.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'complete':True,'source_regions':len(result['regions']),
        'late_nearest_independent_countries':{
          key:result['country_equal_summary']['2018_2025']['nearest_in_cell'][key]['n_independent_country_blocks_supported']
          for key in ('herbaceous','woody','SAME_COHORT_GROUP_DIFFERENCE')},
        'NO_stress_adaptation_or_seed_transport_claim':True}),flush=True)

if __name__=='__main__':main()
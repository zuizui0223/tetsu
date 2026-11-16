#!/usr/bin/env python3
"""v7.4: six-country GRIIS alien flora, GIFT growth forms, Baselga 2012
Jaccard beta-diversity decomposition versus physical RAIL/ROAD routes.

TRUE measured source facets, not observed railway seed transport. FULL vs
WITHOUT dominant GBIF dataset and herbaceous vs woody use identical cell
pairs within each period/anchor. 50km pairs are dependent. No p-values.
Rarefied20 beta is an algebraic ratio-of-expected-inclusion extension of
Baselga's incidence decomposition, NOT E[beta] or standardized surveys.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from itertools import combinations
from statistics import median
from datetime import datetime,timezone

from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,PERIODS,STRATEGIES,combine_source_regions,MIN_SOURCE_CELLS,MIN_VALID_PAIRS
from analyze_GRIIS_topdataset_leaveout_rail_macro_v57 import merge_species_sources
from analyze_growthform_source_leaveout_v72 import fixed_group_source_cells,COUNT_KEYS
from analyze_sixcountry_functional_rail_v71 import GROUPS,validate_traits,eligible_group_counts
from analyze_macro_50km_rail_road_aliens_v34 import _routes,partial_rank_r
from analyze_macro_rail_assemblages_v22 import incidence_probabilities,expected_rarefied_overlap

METRICS=("beta_jac","beta_jtu","beta_jne")
FAMILIES=("raw","rarefied20")
K=20


def baselga_jaccard_components(p,q):
    """On incidence p,q: Jaccard dissimilarity = turnover + nestedness.
    
    For probabilities, decompose the ratio of expected shared/exclusive
    occurrences; not the expected value of beta diversity.
    """
    names=set(p)|set(q)
    if not names:return None
    if any(type(v) not in (int,float) or not math.isfinite(v) or not 0<=v<=1
           for v in [*p.values(),*q.values()]):
        raise ValueError("Species inclusion must be finite in [0,1]")
    a=sum(p.get(k,0)*q.get(k,0) for k in names)
    b=sum(p.get(k,0)*(1-q.get(k,0)) for k in names)
    c=sum(q.get(k,0)*(1-p.get(k,0)) for k in names)
    if a+b+c<=0:return None
    total=(b+c)/(a+b+c)
    m=min(b,c)
    if a+2*m<=0:return None
    turn=2*m/(a+2*m)
    nest=total-turn
    if not math.isclose(turn+nest,total,abs_tol=1e-8) or nest< -1e-8:
        raise ValueError("Baselga decomposition identity failed")
    return {"beta_jac":total,"beta_jtu":turn,"beta_jne":max(0.,nest)}


def pair_beta(counts_a,counts_b):
    raw=baselga_jaccard_components(
        {k:1 for k in counts_a},{k:1 for k in counts_b})
    pa=incidence_probabilities(counts_a,k=K)
    pb=incidence_probabilities(counts_b,k=K)
    rare=baselga_jaccard_components(pa,pb)
    if rare is not None:
        similarity=expected_rarefied_overlap(pa,pb)
        if similarity is None or not math.isclose(rare["beta_jac"],1-similarity,abs_tol=1e-8):
            raise ValueError("Source-v7.2 rarefied Jaccard invariant failed")
    return {"raw":raw,"rarefied20":rare}


def make_source_common_pairs(block,period,anchor,flora,split,rail,road,traits):
    if split.get("complete") is not True or split.get("source_errors") or split.get("source_geo_block")!=block:
        raise ValueError("Original GBIF source split failed provenance gate")
    selected=fixed_group_source_cells(block,period,flora,split,traits)
    ar=_routes(rail,block,anchor,"pair_routes")
    rd=_routes(road,block,anchor,"pair_paths")
    if set(ar)!=set(rd):raise ValueError("Physical RAIL/ROAD source pairs disagree")
    rows=[]
    for a,b in combinations(sorted(selected),2):
        r=ar[a,b];d=rd[a,b]
        rp=r.get("cell_centroid_snap_augmented_rail_path_km")
        dp=d.get("centroid_snap_augmented_road_path_km")
        geo=r.get("source_cell_centroid_geodesic_km")
        ratios=(r.get("rail_to_geodesic_ratio"),d.get("road_to_geodesic_ratio"))
        if not all(type(v) in (int,float) and math.isfinite(v) and v>0
                   for v in (rp,dp,geo,*ratios)):continue
        if min(ratios)<.95:continue
        c1=selected[a][period]["all_plant_indexed_records"]
        c2=selected[b][period]["all_plant_indexed_records"]
        if min(c1,c2)<=0:raise ValueError("Missing all-Plantae source record effort")
        response={}
        for source,field in COUNT_KEYS.items():
            response[source]={}
            for group in GROUPS:
                counts=[]
                for cell in (a,b):
                    p=split["source_cells"][cell]["periods"][period][field]
                    x=eligible_group_counts({"counts":p},traits,group)
                    if x is None:raise ValueError("Eligible provider/group source vanished")
                    counts.append(x)
                response[source][group]=pair_beta(*counts)
        if any(response[s][g][family] is None
               for s in COUNT_KEYS for g in GROUPS for family in FAMILIES):continue
        rows.append({"a":a,"b":b,"x":math.log(rp/dp),
                     "z":[math.log(geo),math.log(dp),abs(math.log(c1/c2)),
                          (math.log(c1)+math.log(c2))/2],
                     "response":response})
    return rows,len(selected)


def one_country(rows,eligible_cells):
    n_cells=len({c for row in rows for c in (row["a"],row["b"])})
    enough=n_cells>=MIN_SOURCE_CELLS and len(rows)>=MIN_VALID_PAIRS
    output={"n_common_cells_with_both_forms_and_providers":eligible_cells,
            "n_cells_in_physical_pairs":n_cells,
            "n_dependent_common_50km_cell_pairs":len(rows),
            "source_support":enough,"coefficients":{}}
    for source in COUNT_KEYS:
        output["coefficients"][source]={}
        for group in GROUPS:
            output["coefficients"][source][group]={}
            for family in FAMILIES:
                v={}
                for metric in METRICS:
                    y=[row["response"][source][group][family][metric] for row in rows]
                    rho=partial_rank_r([row["x"] for row in rows],y,
                                       [row["z"] for row in rows]) if enough else None
                    v[metric]=rho
                med={metric:round(median(row["response"][source][group][family][metric]
                                          for row in rows),6)
                     for metric in METRICS} if rows else {}
                output["coefficients"][source][group][family]={
                    "partial_rho_rail_ratio_vs_dissimilarity":v,
                    "median_dependent_pair_beta":med}
    return output


def equal_country_summary(regions,period,anchor):
    result={}
    for source in COUNT_KEYS:
        result[source]={}
        for group in GROUPS:
            result[source][group]={}
            for family in FAMILIES:
                result[source][group][family]={}
                for metric in METRICS:
                    x={name:rec["periods"][period][anchor]["coefficients"][source][group][family]
                       ["partial_rho_rail_ratio_vs_dissimilarity"][metric]
                       for name,rec in regions.items()}
                    vals={name:v for name,v in x.items() if v is not None}
                    result[source][group][family][metric]={
                        "country_partial_rhos":vals,
                        "n_independent_country_blocks":len(vals),
                        "positive_dissimilarity_sign_hypothesis_compatible":sum(v>0 for v in vals.values()),
                        "negative_dissimilarity_sign":sum(v<0 for v in vals.values()),
                        "country_equal_median_rho":round(median(vals.values()),6) if vals else None,
                        "leave_one_country_out_medians":{
                            c:round(median(v for name,v in vals.items() if name!=c),6)
                            for c in vals if len(vals)>1},
                        "five_country_support":len(vals)>=5}
    return result


def run(flora,gift,splits,rail,road):
    traits=validate_traits(gift)
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Original flora source is incomplete")
    if any(set(x.get("regions",{}))!=set(BLOCKS_SIX) for x in (splits,rail,road)):
        raise ValueError("Original six-country geographic frame mismatch")
    out={"complete":True,"utc":datetime.now(timezone.utc).isoformat(),
         "posthoc_metric_decomposition_not_pre_registered_ecological_hypothesis":True,
         "method":"Baselga 2012 Jaccard turnover and nestedness components, raw incidence and algebraic rarefied inclusion ratio",
         "repetition_unit":"one selected country block, NOT 50km cell pairs",
         "metric_sign":"positive partial rho between RAIL/ROAD route ratio and beta-dissimilarity is railway-hypothesis-compatible",
         "source_GIFT_growth_form_ID":"1.2.1","source_record_rarefaction_K":K,
         "regions":{},"country_equal_summary":{},
         "NO_observed_freight_seed_flow_mowing_wind_or_adaptation":True}
    for block in BLOCKS_SIX:
        out["regions"][block]={"periods":{}}
        for period in PERIODS:
            out["regions"][block]["periods"][period]={}
            for anchor in STRATEGIES:
                rows,n=make_source_common_pairs(block,period,anchor,flora,splits["regions"][block],rail,road,traits)
                out["regions"][block]["periods"][period][anchor]=one_country(rows,n)
    for period in PERIODS:
        out["country_equal_summary"][period]={
            anchor:equal_country_summary(out["regions"],period,anchor)
            for anchor in STRATEGIES}
    out["limits"]=[
        "Physical 2026 rail networks are not observed historical active train services or seed movements.",
        "Only 6 selected geographic country blocks; common pairs share cells and yield no pairwise p-values.",
        "Source GBIF facets are indexed recorded species, not independent standardized plots or biological absences.",
        "GIFT growth forms are not measured plant stress tolerance or true adaptation.",
        "Dominant GBIF dataset UUID removal is a source sensitivity check, not observer-level independence.",
        "Raw observed taxa beta is vulnerable to source record depth; rarefied beta partitions the ratio of EXPECTED inclusion and is NOT expectation of Baselga beta.",
        "Turnover or nestedness associations can result from unmeasured climate/land cover/recording effort, not measured railway delivery/filtering.",
        "Both source definitions and growth forms use identical cell pair cohort within period and anchor."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for k in ("flora","gift","output"):p.add_argument("--"+k,required=True)
    for k in ("split-sources","rail-sources","road-sources"):
        p.add_argument("--"+k,nargs="+",required=True)
    a=p.parse_args();load=lambda x:json.loads(Path(x).read_text(encoding="utf-8"))
    flora=load(a.flora);gift=load(a.gift)
    splits=merge_species_sources([load(x) for x in a.split_sources])
    rail=combine_source_regions([load(x) for x in a.rail_sources])
    road=combine_source_regions([load(x) for x in a.road_sources])
    result=run(flora,gift,splits,rail,road)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    main=result["country_equal_summary"]["2018_2025"]["nearest_in_cell"]
    print(json.dumps({"complete":True,"period":"2018_2025","K":K,
        "primary":main,"no_transport_or_stress_causality":True},ensure_ascii=False),flush=True)

if __name__=="__main__":main()

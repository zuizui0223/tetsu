#!/usr/bin/env python3
"""v11.5: Does the REAL Slovenia RAIL source 261-taxon weighting reversal
survive keeping original left/right sides separate?

Original Šabić, Grašič & Rovanšek (2026), DOI 10.1007/s44473-026-00155-z;
89 true sampled locations, 173 matched same-side BS/BT pairs, 15 source
rail sections. Secondary analysis of authors' source observations, NOT
new original tetsu field surveys.

Three diagnostics:
1) Original union(L/R) per site versus within-site average of paired
   track-side BS/BT presence preference, identical taxon universe.
2) Explicit frequency exponent 0=equal taxon, 1=taxon source site
   opportunity weighting with root (if any).
3) Exact ALL 2^15 section-level BS/BT label-flip distribution holding
   ALL original within-site co-occurrences, taxon source occurrence
   frequencies, lifespan and genus groups FIXED. This is a descriptive
   exchangeability stress test, NOT causal inference because BS/BT
   sampling areas/conditions differ and sites were not randomized.

Cannot measure actual train gust, freight seed delivery, mowing,
phenotypic plasticity or adaptive genetic evolution.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_DOI, SOURCE_SHA256, extract_original_xlsx,
    validate_original_keys, pair_source,
)
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    original_location_matched_sets, taxa_preference, STRATA,
)
from analyze_Slovenia_rail_species_occupancy_weighting_v110 import group_estimator

N_SOURCE_SECTIONS=15
N_SOURCE_SITES=89
N_SOURCE_SIDES=173
N_ORIGINAL_TAXA=261


def extract_side_pairs(original):
    """Leave true source same-side field-pair structure unmodified."""
    loc,zones,occ,traits,_,quality = validate_original_keys(original)
    paired = pair_source(loc,zones,occ,traits)
    result=defaultdict(list)
    for x in paired:
        lid=x["location_id"]; side=x["track_side"]
        if side not in ("L","R"):
            raise ValueError("Invalid real source track side")
        result[lid].append({
            "section":str(x["section_number"]),
            "side":side,
            "bs":occ[(lid,side+"BS")],
            "bt":occ[(lid,side+"BT")],
        })
    if len(result)!=N_SOURCE_SITES or sum(map(len,result.values()))!=N_SOURCE_SIDES:
        raise ValueError("Real original 89 sites and 173 sides changed")
    if sum(len(v)==2 for v in result.values())!=84:
        raise ValueError("Original real 84 double-side sites changed")
    if any(len({z["side"] for z in v})!=len(v) or
           len({z["section"] for z in v})!=1 for v in result.values()):
        raise ValueError("Duplicate real source sides or conflicting section")
    return result,traits,quality


def preference_site_matrix(sides, traits, *, mode="union", both_sides_only=False):
    """One row per taxon: each *site* gets ONE q contribution.

    'union': Site BS = union of its own eligible matched-side BS plots,
    same for BT.  A taxon present both BS/BT gets zero.
    'side_average': If a taxon occurs in >0 sides of a given site,
    first average (BS-present - BT-present) over ONLY the matched
    sides where that taxon was detected. Then average these per-site
    values over all occupied source sites for that taxon.
    Both estimands count one site as one occupancy opportunity.
    """
    if mode not in ("union","side_average"):
        raise ValueError("Undefined source side aggregation")
    source_matrix=defaultdict(dict)
    sections={}
    eligible_sites=sorted(lid for lid,ss in sides.items()
                         if not both_sides_only or len(ss)==2)
    for lid in eligible_sites:
        ss=sides[lid]
        sections[lid]=ss[0]["section"]
        observed=set().union(*(z["bs"]|z["bt"] for z in ss))
        for taxon in sorted(observed):
            if mode=="union":
                in_bs=any(taxon in z["bs"] for z in ss)
                in_bt=any(taxon in z["bt"] for z in ss)
                delta=float(in_bs)-float(in_bt)
            else:
                within=[float(taxon in z["bs"])-float(taxon in z["bt"])
                        for z in ss if taxon in z["bs"] or taxon in z["bt"]]
                if not within:
                    raise ValueError("Taxon appears on no matched side")
                delta=mean(within)
            source_matrix[taxon][lid]=delta
    rows=[]
    for taxon,values in sorted(source_matrix.items()):
        info=traits[taxon]
        wind=info.get("Diaspore.strategy")
        lifespan=info.get("Lifespan.categories")
        if wind is None or wind=="" or lifespan is None or lifespan=="":
            continue
        if not set(values).issubset(set(eligible_sites)):
            raise ValueError("Site matrix has invented field sites")
        rows.append({
            "taxon_code":taxon,
            "wind":wind in ("Epilobium","Phragmites"),
            "stratum":"strict_annual" if lifespan=="annual" else "other_known_lifespan",
            "n_original_locations_present_either_zone":len(values),
            "occupancy_preference_BS_minus_BT":mean(values.values()),
            "site_delta":values
        })
    return rows,sections


def continuum(rows, power):
    """At exponent 0 => v11.0 equal-species, 1 => v11.0 opportunity."""
    if not math.isfinite(power) or power<0:
        raise ValueError("Source weighting exponent must be nonnegative")
    sub={s:{"wind":[],"other":[]} for s in STRATA}
    for r in rows:
        n=r["n_original_locations_present_either_zone"]
        if n<1: raise ValueError("No original site opportunity")
        sub[r["stratum"]]["wind" if r["wind"] else "other"].append((r,n**power))
    if any(min(len(sub[s]["wind"]),len(sub[s]["other"]))<3 for s in STRATA):
        return None
    total=0.0
    num=0.0
    for s in STRATA:
        z=sub[s]
        dw=sum(w for _,w in z["wind"])
        dn=sum(w for _,w in z["other"])
        mw=sum(r["occupancy_preference_BS_minus_BT"]*w for r,w in z["wind"])/dw
        mn=sum(r["occupancy_preference_BS_minus_BT"]*w for r,w in z["other"])/dn
        total+=dw+dn
        num+=(dw+dn)*(mw-mn)
    return num/total


def crossing(rows):
    v0=continuum(rows,0)
    v1=continuum(rows,1)
    curve=[{"frequency_power":round(i/20,2),
            "wind_minus_other_BS_preference":round(continuum(rows,i/20),9)}
           for i in range(21)]
    if v0 is None or v1 is None:
        return {"identified":False,"curve":curve}
    root=None
    if v0==0:root=0.0
    elif v1==0:root=1.0
    elif v0*v1<0:
        lo,hi=0.0,1.0
        for _ in range(48):
            mid=(lo+hi)/2
            if continuum(rows,mid)*v0>0:lo=mid
            else:hi=mid
        root=(lo+hi)/2
    return {"identified":True,"n_species":len(rows),
            "species_equal":v0,"site_opportunity_weighted":v1,
            "source_frequency_exponent_zero_crossing_0_to_1":root,
            "curve":curve,"root_is_descriptive_not_a_selection_threshold":True}


def section_contributions(rows,sections,power):
    """Exact linear additive site contributions, frequencies held fixed."""
    grouped={s:{"wind":[],"other":[]} for s in STRATA}
    for r in rows:grouped[r["stratum"]]["wind" if r["wind"] else "other"].append(r)
    totals={}
    for stratum in STRATA:
        for guild in ("wind","other"):
            arr=grouped[stratum][guild]
            if len(arr)<3:return None
            totals[(stratum,guild)]=sum(
                r["n_original_locations_present_either_zone"]**power for r in arr)
    whole=sum(totals.values())
    coefs={}
    for stratum in STRATA:
        weight=(totals[(stratum,"wind")]+totals[(stratum,"other")])/whole
        for guild in ("wind","other"):
            denom=totals[(stratum,guild)]
            sign=1.0 if guild=="wind" else -1.0
            for r in grouped[stratum][guild]:
                n=r["n_original_locations_present_either_zone"]
                coefs[r["taxon_code"]]=sign*weight*(n**power)/denom/n
    out={sec:0.0 for sec in sorted(set(sections.values()))}
    for r in rows:
        a=coefs[r["taxon_code"]]
        for lid,delta in r["site_delta"].items():
            out[sections[lid]]+=a*delta
    if len(out)!=N_SOURCE_SECTIONS and len(sections)==N_SOURCE_SITES:
        raise ValueError("Real source railway sections changed")
    if not math.isclose(sum(out.values()),continuum(rows,power),abs_tol=1e-10):
        raise ValueError("Mathematically incorrect additive source site statistic")
    return out


def exact_section_label_swaps(rows,sections):
    """All 32768 Rademacher sign patterns: swap BS and BT for entire section."""
    equal=section_contributions(rows,sections,0)
    weighted=section_contributions(rows,sections,1)
    if equal is None or weighted is None or set(equal)!=set(weighted):
        raise ValueError("Wind/lifespan group source metrics not identifiable")
    keys=sorted(equal)
    if len(keys)!=N_SOURCE_SECTIONS:
        raise ValueError("Cannot enumerate exactly 15 original rail sections")
    base_e=sum(equal.values())
    base_w=sum(weighted.values())
    vals_e=[]
    vals_w=[]
    n_extreme_e=n_extreme_w=n_extreme_gap=n_opposite=0
    gap=base_e-base_w
    for pattern in itertools.product((-1.,1.),repeat=len(keys)):
        e=sum(sign*equal[k] for sign,k in zip(pattern,keys))
        w=sum(sign*weighted[k] for sign,k in zip(pattern,keys))
        vals_e.append(e);vals_w.append(w)
        if abs(e)+1e-12>=abs(base_e):n_extreme_e+=1
        if abs(w)+1e-12>=abs(base_w):n_extreme_w+=1
        if abs(e-w)+1e-12>=abs(gap):n_extreme_gap+=1
        if e>0 and w<0:n_opposite+=1
    n=2**len(keys)
    if len(vals_e)!=n:
        raise ValueError("Exact sign-flip enumeration incomplete")
    return {
        "original_15_rail_section_species_equal":base_e,
        "original_15_rail_section_occupancy_weighted":base_w,
        "original_gap_equal_minus_occupancy":gap,
        "n_possible_section_exchangeability_flips":n,
        "two_sided_full_15section_flip_p_species_equal":n_extreme_e/n,
        "two_sided_full_15section_flip_p_occ_weighted":n_extreme_w/n,
        "two_sided_full_15section_flip_p_weighting_gap":n_extreme_gap/n,
        "fraction_null_flips_with_equal_positive_and_occ_negative":n_opposite/n,
        "label_swaps_preserve_original_taxon_frequency_and_within_section_plant_composition":True,
        "original_BS_BT_zone_area_unequal_exchangeability_unverified":True,
        "NOT_randomized_causal_or_global_RAIL_inference":True
    }


def run(original,receipt):
    if (receipt.get("complete") is not True or receipt.get("version")!="v11.0"
        or receipt.get("original_source_XLSX_sha256")!=SOURCE_SHA256):
        raise ValueError("Locked original v11.0 actual source data required")
    sites,traits,integrity=original_location_matched_sets(original)
    actual,sections=preference_site_matrix(*extract_side_pairs(original)[:2],mode="union")
    sides,traits2,_=extract_side_pairs(original)
    if traits!=traits2:
        raise ValueError("Original field species trait source drift")
    if len(actual)!=N_ORIGINAL_TAXA or len(sections)!=N_SOURCE_SITES:
        raise ValueError("True original 261 species/89 locations not reproduced")
    v0=continuum(actual,0)
    v1=continuum(actual,1)
    if not math.isclose(v0,receipt["full_1plus_species_equal_wind_minus_other_preference"],abs_tol=1e-12):
        raise ValueError("Original v11.0 species-equal contrast not reproduced")
    if not math.isclose(v1,receipt["full_1plus_occupancy_weighted_wind_minus_other_preference"],abs_tol=1e-12):
        raise ValueError("Original v11.0 occurrence-weight contrast not reproduced")
    side,side_sections=preference_site_matrix(sides,traits,mode="side_average")
    if {r["taxon_code"] for r in side}!={r["taxon_code"] for r in actual}:
        raise ValueError("Source original taxa changed after within-site paired-side corrections")
    double_union,_=preference_site_matrix(sides,traits,mode="union",both_sides_only=True)
    double_side,_=preference_site_matrix(sides,traits,mode="side_average",both_sides_only=True)
    info={
        "complete":True,"version":"v11.5",
        "original_paper_DOI":SOURCE_DOI,
        "original_source_SHA256":SOURCE_SHA256,
        "n_original_true_RAIL_field_sites":len(sides),
        "n_original_true_matched_track_sides":sum(len(v) for v in sides.values()),
        "n_original_sites_with_two_matched_sides":sum(len(v)==2 for v in sides.values()),
        "n_original_source_rail_sections":len(set(sections.values())),
        "n_original_trait_known_detected_source_taxa":len(actual),
        "source_original_LR_pooled_reproduces_v11":crossing(actual),
        "same_261_original_taxa_average_each_matched_side_within_site":crossing(side),
        "only_84_both-side_sites_union_original_taxa_subset":crossing(double_union),
        "only_84_both-side_sites_side_average_original_taxa_subset":crossing(double_side),
        "exact_32768_rail_section_BS_BT_label_swaps_original_union":exact_section_label_swaps(actual,sections),
        "exact_32768_rail_section_BS_BT_label_swaps_side_average":exact_section_label_swaps(side,side_sections),
        "cautions":[
            "BS/BT are original source sampling zones; their sampled footprint areas may differ, so rail section label exchangeability is not a valid randomized intervention.",
            "Source track-side averaging and source-location union are different defensible but non-equivalent presence estimands.",
            "The observed site preference score is not frequency, individual abundance, seed carriage or plant fitness; true wind dispersal under ambient vs passing train gusts was not measured.",
            "One original country with non-random field sampling cannot show global rail adaptation or train seed transport.",
        ],
    }
    return info


def main():
    p=argparse.ArgumentParser()
    for n in ("original-xlsx","source-audit102","verified110","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    data=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    verified=json.loads(Path(a.verified110).read_text(encoding="utf-8"))
    result=run(data,verified)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v11.5",
        "actual_sources":{"sites":result["n_original_true_RAIL_field_sites"],
                          "matched_sides":result["n_original_true_matched_track_sides"],
                          "taxa":result["n_original_trait_known_detected_source_taxa"]},
        "original":{k:v for k,v in result["source_original_LR_pooled_reproduces_v11"].items() if k!="curve"},
        "within_site_side_average":{k:v for k,v in result["same_261_original_taxa_average_each_matched_side_within_site"].items() if k!="curve"},
        "84_both_sides_union":{k:v for k,v in result["only_84_both-side_sites_union_original_taxa_subset"].items() if k!="curve"},
        "84_both_sides_side_average":{k:v for k,v in result["only_84_both-side_sites_side_average_original_taxa_subset"].items() if k!="curve"},
        "exact_15section_swap":result["exact_32768_rail_section_BS_BT_label_swaps_original_union"],
        "exact_15section_side_swap":result["exact_32768_rail_section_BS_BT_label_swaps_side_average"],
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

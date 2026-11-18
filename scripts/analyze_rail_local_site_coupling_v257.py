#!/usr/bin/env python3
"""v25.7 spatial site-coupling null for two authentic railway flora studies.

Original same-site rail/reference overlap is contrasted with overlap under
reference SITE-identity derangements within a defensible source stratum:
Kashmir station-vs-between groups (same derangement applied to both years);
Slovenia original line section (14 with >=2 observed site locations).

This keeps every *observed* plant species assemblage and plot richness fixed;
only local pairing labels are changed. The permutation diagnostic CANNOT
separate direct train dispersal, physical microenvironment, geographic
distance, landscape species availability, management, plot areas or detection.
It is not a field treatment p-value, genetic effect, or a rail/road comparison.

No source GPS or plant-by-site raw tables are written to output.
"""
from __future__ import annotations
import argparse
import json
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix, matching_sites
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx, validate_original_keys

KZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
KREFS=("Embankment","Interior Plot 1","Interior Plot 2")
GROUPS=("station","between")
YEARS=(2014,2017)
SEED=257118


def jac_similarity(a: set, b: set) -> float:
    if not a or not b:raise ValueError("Only nonempty source occurrence sets")
    return len(a&b)/len(a|b)


def derangement(ids: list, rng: random.Random) -> dict:
    """Bijective source-site reassignment that never retains the same site."""
    if len(ids)<2:
        raise ValueError("A singleton stratum cannot identify spatial pairing null")
    # Fisher-Yates shuffles, rejected if any original site remains paired:
    # conditional on no fixed points, this is uniform over all derangements.
    for attempt in range(500):
        shuffled=list(ids)
        rng.shuffle(shuffled)
        if all(a!=b for a,b in zip(ids,shuffled)):
            if len(set(shuffled))!=len(ids):
                raise AssertionError("Non-bijective source reference reassignment")
            return dict(zip(ids,shuffled))
    raise RuntimeError("Could not draw original-site derangement")


def quantile(a: list[float],p: float) -> float:
    x=sorted(a);k=(len(x)-1)*p;i=int(k)
    return x[i] if i==len(x)-1 else x[i]*(1-k+i)+x[i+1]*(k-i)


def diagnostic(observed: float, null: list[float]) -> dict:
    if len(null)<99:
        raise ValueError("Insufficient fixed-assemblage label reassignment draws")
    return {
        "same_site_observed_jaccard_similarity":observed,
        "different_site_label_null_mean":mean(null),
        "same_site_minus_reassigned_mean":observed-mean(null),
        "technical_q025":quantile(null,.025),
        "technical_q975":quantile(null,.975),
        "n_derangement_draws":len(null),
        "n_reassignments_equal_or_more_similar":sum(z>=observed-1e-12 for z in null),
        "upper_tail_mapping_rank_plusone_NOT_a_field_p_value":(1+sum(z>=observed-1e-12 for z in null))/(len(null)+1),
        "original_source_plant_sets_and_plot_richness_unchanged":True,
    }


def analyze_kashmir(xlsx: Path,draws:int=999) -> dict:
    plots,source=parse_xlsx_matrix(xlsx)
    rng=random.Random(SEED)
    result={}
    for group in GROUPS:
        sites=matching_sites(plots,group,KZONES)
        if len(sites)<10:raise ValueError("Strict original Kashmir common site cohort too small")
        actual={
            str(year):{
                ref:mean(jac_similarity(plots[s,year,"Track"],plots[s,year,ref]) for s in sites)
                for ref in KREFS
            } for year in YEARS
        }
        reference_null={str(year):{ref:[] for ref in KREFS} for year in YEARS}
        for _ in range(draws):
            mapping=derangement(sites,rng)
            for year in YEARS:
                for ref in KREFS:
                    score=mean(jac_similarity(plots[s,year,"Track"],plots[mapping[s],year,ref]) for s in sites)
                    reference_null[str(year)][ref].append(score)
        result[group]={
            "n_original_spatial_sites":len(sites),
            "same_4_zone_repeated_sites_both_years":True,
            "source_group_not_real_train_frequency":True,
            "by_year":{
                str(year):{
                    ref:diagnostic(actual[str(year)][ref],reference_null[str(year)][ref])
                    for ref in KREFS
                } for year in YEARS
            }
        }
    return {
        "original_Dryad_source_qc":source,
        "reference_sites_permuted_only_within_station_between_source_stratum":True,
        "same_reference_site_derangement_applied_both_years_in_every_draw":True,
        "source_stratum_can_have_spatial_or_climate_gradients_not_measured_here":True,
        "results":result,
    }


def slovenia_location_sets(original: dict) -> tuple[dict,dict]:
    loc,zone,sp,tr,cb,audit=validate_original_keys(original)
    items=defaultdict(lambda:{"BS":set(),"BT":set(),"section":None,"matched_sides":0})
    matched=0
    for (lid,zid),one in zone.items():
        if one["area"]!="BS":continue
        if zid not in ("LBS","RBS"):raise ValueError("Not published source shoulder side")
        counterpart=zid[0]+"BT"
        other=zone.get((lid,counterpart))
        if not other:continue
        if other["area"]!="BT" or one.get("area_side")!=other.get("area_side"):
            raise ValueError("Original matched left/right track zone mismatch")
        q=items[lid]
        section=str(loc[lid]["track_section_number"])
        if q["section"] is not None and q["section"]!=section:
            raise ValueError("A true source field location cannot span sections")
        q["section"]=section
        q["matched_sides"]+=1
        q["BS"].update(sp[(lid,zid)])
        q["BT"].update(sp[(lid,counterpart)])
        matched+=1
    if matched!=173 or len(items)!=89 or len({v["section"] for v in items.values()})!=15:
        raise ValueError("True source Slovenia 173 side pairs / 89 sites / 15 original sections drift")
    if any(not q["BS"] or not q["BT"] for q in items.values()):
        raise ValueError("Empty actual original same-site matched zone incidence")
    return items,{
        "original_source_qc":audit,
        "n_paired_source_sides":matched,
        "n_original_locations":len(items),
        "n_original_sections":15,
        "both_sides_pooled_within_same_original_location_to_preserve_site_unit":True,
        "source_BS_BT_zone_area_equivalence_unverified":True,
    }


def analyze_slovenia(xlsx: Path,audit_path:Path,draws:int=999) -> dict:
    original=extract_original_xlsx(xlsx,audit_path)
    sites,source=slovenia_location_sets(original)
    sections=defaultdict(list)
    for lid,q in sites.items():
        sections[q["section"]].append(lid)
    retained={sec:sorted(locations,key=str) for sec,locations in sections.items() if len(locations)>=2}
    excluded={sec:len(locs) for sec,locs in sections.items() if len(locs)<2}
    if len(retained)<12:
        raise ValueError("Not enough independent source rail sections for within-section location mapping null")
    raw_by_section={
        sec:mean(jac_similarity(sites[lid]["BS"],sites[lid]["BT"]) for lid in ids)
        for sec,ids in retained.items()
    }
    actual=mean(raw_by_section.values())
    rng=random.Random(SEED+289)
    null=[]
    for _ in range(draws):
        by=[]
        for sec,ids in sorted(retained.items()):
            perm=derangement(ids,rng)
            by.append(mean(jac_similarity(sites[lid]["BS"],sites[perm[lid]]["BT"]) for lid in ids))
        null.append(mean(by))
    return {
        "source":source,
        "original_source_multi_location_sections":len(retained),
        "single_location_original_sections_excluded_no_fabricated_match":len(excluded),
        "original_locations_in_pairable_sections":sum(map(len,retained.values())),
        "observation_unit":"site union across observed L/R same-side matched BS/BT; section-equal",
        "within_section_reference_site_derangement_not_randomized_field_treatment":True,
        "diagnostic":diagnostic(actual,null),
        "leave_one_source_section_out_observed_min_max":{
            "min":min(mean(v for s,v in raw_by_section.items() if s!=omit) for omit in retained),
            "max":max(mean(v for s,v in raw_by_section.items() if s!=omit) for omit in retained),
        },
    }


def run(kashmir:Path,slovenia:Path,audit:Path,draws:int=999) -> dict:
    if draws<99:raise ValueError("Require >=99 synthetic label assignments")
    india=analyze_kashmir(kashmir,draws)
    si=analyze_slovenia(slovenia,audit,draws)
    return {
        "version":"v25.7",
        "status":"TWO_GENUINE_RAILWAY_STUDIES_SITE_PAIRING_NULL_IF_SOURCE_RUN_SUCCEEDS",
        "interpretation":"Constrained source site-identity reassignment diagnostics on observed communities",
        "Kashmir":india,
        "Slovenia":si,
        "studies_analyzed_separately":2,
        "no_nonrail_experimental_control_no_passenger_or_freight_ORIGIN_DESTINATION":True,
        "same_site_vs_other_source_site_does_not_identify_specific_environmental_or_dispersal_mechanism":True,
        "site_reassignment_null_unadjusted_actual_geographic_distance_and_microclimate":True,
        "derangement_draws_not_independent_field_replicates":True,
        "derived_aggregate_outputs_only_no_source_taxon_or_GPS_export":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-source",type=Path,required=True)
    p.add_argument("--slovenia-xlsx",type=Path,required=True)
    p.add_argument("--slovenia-audit",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=999)
    a=p.parse_args()
    out=run(a.kashmir_source,a.slovenia_xlsx,a.slovenia_audit,a.draws)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"version":"v25.7","sources":2,
        "Kashmir_actual_station_2014_track_interior2":
        out["Kashmir"]["results"]["station"]["by_year"]["2014"]["Interior Plot 2"]["same_site_observed_jaccard_similarity"],
        "Kashmir_matched_minus_rematched_station_2014_track_interior2":
        out["Kashmir"]["results"]["station"]["by_year"]["2014"]["Interior Plot 2"]["same_site_minus_reassigned_mean"],
        "Slovenia_section_equal_observed":
        out["Slovenia"]["diagnostic"]["same_site_observed_jaccard_similarity"],
        "Slovenia_matched_minus_rematched":
        out["Slovenia"]["diagnostic"]["same_site_minus_reassigned_mean"]
    },ensure_ascii=False))


if __name__=="__main__":
    main()

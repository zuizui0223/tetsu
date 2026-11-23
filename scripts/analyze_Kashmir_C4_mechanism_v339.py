#!/usr/bin/env python3
"""v339: strict named C4 pathway contribution to Kashmir railway thermal EIVE.

Original Dryad SHA and published independent EIVE trait checks inherited.
Original sites, positive source abundance and geographical matches never exported.
C4 is an EXPLORATORY published-pathway candidate: unknown != C3.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import (
    species_names_from_publisher, site_effect
)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_geographic_station_between_EIVE_v332 import (
    site_metadata, TRAIT_GATES
)
from analyze_Kashmir_cover_weighted_matched_temperature_v334 import (
    read_original_plant_values, cohort, group_comparison
)

C4 = frozenset((
    "amaranthus caudatus", "bothriochloa ischaemum", "cynodon dactylon",
    "cyperus rotundus", "echinochloa colona", "euphorbia prostrata",
    "fimbristylis dichotoma", "saccharum spontaneum", "setaria viridis",
    "sorghum halepense",
))
C4_GRASS = frozenset((
    "bothriochloa ischaemum", "cynodon dactylon", "echinochloa colona",
    "saccharum spontaneum", "setaria viridis", "sorghum halepense"
))
C4_AMBIGUOUS = frozenset(("euphorbia hispida",))
C3_GRASS = frozenset((
    "arundo donax", "avena sativa", "bromus japonicus", "phleum pratense",
    "poa annua", "poa bulbosa",
))
PANELS = (
    (2014, "strict_3systems_5species_50pct", "G1_geographic_altitude_primary"),
    (2014, "relaxed_1system_3species_30pct", "G1_geographic_altitude_primary"),
    (2014, "relaxed_1system_3species_30pct", "G3_geo_altitude_richness_coverage"),
    (2017, "sparse_3systems_1species", "G1_geographic_altitude_primary")
)
EXCLUSIONS = (
    ("none", frozenset()),
    ("strict_C4", C4),
    ("C4_grasses_only", C4_GRASS),
    ("C4_plus_ambiguous_sensitivity", C4|C4_AMBIGUOUS),
)
ZONES = ("Track", "Interior Plot 2")
WEIGHTS = ("presence", "original_species_abundance")


def mean_safe(x):
    return mean(x) if x else None


def mean_traits(ab, names, eive, systems, weighting="presence", exclude=frozenset()):
    """Site EIVE thermal position, with zero/missing abundance explicitly excluded."""
    if ab is None:
        return None
    values = []
    for i,a in ab.items():
        name=names[i]
        if (name is None or name in exclude or name not in eive
                or eive[name]["systems"] < systems):
            continue
        if not math.isfinite(a) or a<=0:
            raise ValueError("Nonpositive observed source abundance entered")
        values.append((float(eive[name]["position"]), float(a)))
    if not values:
        return None
    w = [1.0 if weighting=="presence" else x[1] for x in values]
    return sum(x[0]*w0 for x,w0 in zip(values,w))/sum(w)


def source_fraction(ab, names, selected, weighted=False):
    if not ab:
        return None
    denom=sum(ab.values()) if weighted else len(ab)
    numerator=sum((val if weighted else 1) for idx,val in ab.items()
                  if names[idx] in selected)
    if denom<=0:
        raise ValueError("Source nonpositive abundance denominator")
    return numerator/denom


def site_c4_fractions(source, abund, names, year):
    """True same-location Track minus Interior2 across original 31 sites."""
    panels=[]
    for site in sorted({k[0] for k in source}):
        p=(site,year,ZONES[0])
        q=(site,year,ZONES[1])
        if source.get(p) is None or source.get(q) is None:
            continue
        if abund.get(p) is None or abund.get(q) is None:
            continue
        dat={"group":"between" if site.startswith("Between ") else "station"}
        for group_name,group_set in (("C4",C4),("C4_grass",C4_GRASS),("C3_grass",C3_GRASS)):
            for mode in (False,True):
                key=group_name+("_source_abundance" if mode else "_species_occurrence")
                dat[key]=source_fraction(abund[p],names,group_set,mode) - source_fraction(abund[q],names,group_set,mode)
        dat["n_C4_track"]=sum(names[i] in C4 for i in source[p])
        dat["n_C4_interior"]=sum(names[i] in C4 for i in source[q])
        panels.append(dat)
    summary={"year":year,"same_source_site_pairs":len(panels),
             "same_site_paired_fraction_results":{}}
    for group in ("all","station","between"):
        sub=[p for p in panels if group=="all" or p["group"]==group]
        keys=sorted(k for k in panels[0] if k not in ("group","n_C4_track","n_C4_interior")) if panels else []
        summary["same_site_paired_fraction_results"][group]={
            "n":len(sub),
            "paired_delta":{k:{"mean":mean_safe([x[k] for x in sub]),
                             "n_track_higher":sum(x[k]>0 for x in sub),
                             "n_interior_higher":sum(x[k]<0 for x in sub),
                             "n_equal":sum(x[k]==0 for x in sub)}
                            for k in keys},
            "mean_C4_species_track":mean_safe([x["n_C4_track"] for x in sub]),
            "mean_C4_species_interior":mean_safe([x["n_C4_interior"] for x in sub])
        }
    return summary


def attenuation(full,remaining):
    if full is None or remaining is None or full<=0:
        return None
    return {"thermal_delta_full":full,
            "thermal_delta_after_C4_removal":remaining,
            "C4_removal_subtraction":full-remaining,
            "fraction_thermal_excess_removed":(full-remaining)/full,
            "fraction_thermal_excess_remaining":remaining/full}


def true_rail_site_cohorts(presence,abund,names,eive,year,traitgate):
    """Cohort fixed by source eligibility with ALL species, not by C4 effects."""
    systems,min_n,coverage = TRAIT_GATES[traitgate]
    eligible={}
    for site in sorted({k[0] for k in presence}):
        delta,meta=site_effect(presence,names,eive,site,year,ZONES,
                               "position",systems,min_n,coverage)
        if delta is not None:
            eligible[site]=delta
    if not eligible:
        return {"n_source_eligible_sites":0}
    byscheme={}
    for weighting in WEIGHTS:
        complete=[]
        skipped=0
        for site,full_source in eligible.items():
            plot_a=abund[(site,year,ZONES[0])]
            plot_b=abund[(site,year,ZONES[1])]
            a=mean_traits(plot_a,names,eive,systems,weighting)
            b=mean_traits(plot_b,names,eive,systems,weighting)
            if a is None or b is None:
                raise AssertionError("Frozen v327 cohort has no source measured baseline")
            if weighting=="presence" and not math.isclose(a-b,full_source,rel_tol=0,abs_tol=1e-9):
                raise AssertionError("Original v327 full EIVE temperature cohort drift")
            results={}
            for lab,excluded in EXCLUSIONS:
                aa=mean_traits(plot_a,names,eive,systems,weighting,excluded)
                bb=mean_traits(plot_b,names,eive,systems,weighting,excluded)
                results[lab]=None if aa is None or bb is None else aa-bb
            complete.append({"source_baseline":a-b,**results})
        base=mean([x["none"] for x in complete])
        filtered={}
        for lab,excluded in EXCLUSIONS:
            sub=[x for x in complete if x[lab] is not None]
            if not sub:
                filtered[lab]={"n_complete_same_original_eligible_sites":0}
                continue
            comparable_full=mean(x["none"] for x in sub)
            current=mean(x[lab] for x in sub)
            filtered[lab]={
              "n_complete_same_original_eligible_sites":len(sub),
              "mean_track_minus_interior_EIVE":current,
              "same_cohort_full_EIVE":comparable_full,
              "mean_C4_removal_effect_within_same_sites":comparable_full-current,
              "attenuation_on_comparable_cohort":attenuation(comparable_full,current)
            }
        one_taxon={}
        for sp in sorted(C4):
            changes=[]
            for site in eligible:
                a,b=abund[(site,year,ZONES[0])],abund[(site,year,ZONES[1])]
                aa=mean_traits(a,names,eive,systems,weighting,frozenset((sp,)))
                bb=mean_traits(b,names,eive,systems,weighting,frozenset((sp,)))
                if aa is not None and bb is not None:
                    changes.append(aa-bb)
            one_taxon[sp]={"n_retained":len(changes),
                          "thermal_delta":mean_safe(changes)}
        byscheme[weighting]={
            "n_frozen_source_sites":len(eligible),
            "original_full_thermal_delta":base,
            "exclusion_sensitivities":filtered,
            "one_C4_taxon_removed":one_taxon,
        }
    group=Counter("between" if site.startswith("Between ") else "station"
                  for site in eligible)
    return {"year":year,"source_trait_gate":traitgate,
            "n_source_eligible_sites":len(eligible),
            "eligible_original_station_sites":group["station"],
            "eligible_original_between_sites":group["between"],
            "weightings":byscheme}


def four_zone_delta(record_values,siteB,siteS,year,names,eive,systems,weighting,exclude):
    x={}
    for key,site,zone in (("BT",siteB,"Track"),("BI",siteB,"Interior Plot 2"),
                          ("ST",siteS,"Track"),("SI",siteS,"Interior Plot 2")):
        x[key]=mean_traits(record_values.get((site,year,zone)),names,eive,
                           systems,weighting,exclude)
        if x[key] is None:
            return None
    return x["BT"]-x["BI"]-x["ST"]+x["SI"]


def four_zone_matched_cohort(presence,abund,names,eive,meta,year,traitgate,matchgate):
    frozen,_,qa = cohort(presence,names,eive,meta,year,traitgate,matchgate)
    if not frozen:
        return {"year":year,"trait_gate":traitgate,"match_gate":matchgate,
                "n_matched_original_geography_pairs":0}
    systems=TRAIT_GATES[traitgate][0]
    scheme_results={}
    for scheme in WEIGHTS:
        original_scheme="presence" if scheme=="presence" else "original_positive_source_value"
        before=group_comparison(abund,frozen,names,eive,year,traitgate,
                                original_scheme,"position")
        all_pairs=[]
        for m in frozen:
            b,s=m["between_source_site"],m["station_source_site"]
            v={}
            for lab,exc in EXCLUSIONS:
                v[lab]=four_zone_delta(abund,b,s,year,names,eive,systems,scheme,exc)
            all_pairs.append(v)
        full=mean(x["none"] for x in all_pairs)
        if not math.isclose(full,before["mean_between_minus_station_track_relative_interior"],
                            abs_tol=1e-9,rel_tol=0):
            raise AssertionError("Original v334 matched 4-zone baseline drift")
        sens={}
        for lab,exc in EXCLUSIONS:
            same=[x for x in all_pairs if x[lab] is not None]
            if same:
                full_cmp=mean(x["none"] for x in same)
                remaining=mean(x[lab] for x in same)
                sens[lab]={"n_original_pairs_with_full_and_excluded":len(same),
                            "thermal_full_on_identical_pairs":full_cmp,
                            "thermal_after_exclusion":remaining,
                            "thermal_excess_change":full_cmp-remaining,
                            "attenuation_on_comparable_cohort":attenuation(full_cmp,remaining)}
            else:
                sens[lab]={"n_original_pairs_with_full_and_excluded":0}
        scheme_results[scheme]={
          "n_original_geographic_pairs":len(frozen),
          "full_source_confirmed_four_plot_EIVE":full,
          "exclusion_sensitivities":sens,
        }
    return {"year":year,"trait_gate":traitgate,"match_gate":matchgate,
            "n_matched_original_geography_pairs":len(frozen),
            "source_matching_QA":qa,"weightings":scheme_results}


def analyze(source_path,eive_path):
    source,sourceqa=parse_xlsx_matrix(Path(source_path))
    names=species_names_from_publisher(source_path)
    eive,eiveqa=eive_table(eive_path)
    meta,_=site_metadata(source_path)
    abund,abundqa=read_original_plant_values(source_path,source)
    matched_taxa=set(x for x in names if x)
    if not C4 <= matched_taxa:
        raise ValueError("All strict C4 names must occur in original source")
    if not C4_GRASS <= C4:
        raise ValueError("C4 grass subset exceeds C4 universe")
    if not C4_AMBIGUOUS <= matched_taxa:
        raise ValueError("Expected C4-sensitive taxon absent")
    out={"version":"v339-C4-railway-thermal-attribution",
        "source_qa":sourceqa,"source_abundance_qa":abundqa,
        "EIVE_md5":eiveqa["EIVE_original_md5"],
        "candidate_taxonomy":{
           "strict_source_named_C4_species":sorted(C4),
           "strict_C4_species_count":len(C4),
           "strict_C4_grass_count":len(C4_GRASS),
           "ambiguous_source_C4_sensitivity_only":sorted(C4_AMBIGUOUS),
           "known_C3_grass_only":sorted(C3_GRASS),
           "source_exact_binomials":len(matched_taxa),
           "C4_group_is_candidate_trait_label_not_lab_confirmation":True,
           "other_source_species_UNCLASSIFIED_not_C3":True,
           "strict_C4_with_independent_EIVE_systems_at_least_1":sum(s in eive for s in C4),
           "strict_C4_with_independent_EIVE_systems_at_least_3":sum(s in eive and eive[s]["systems"]>=3 for s in C4)
        },
        "same_original_site_C4_fractions":[
            site_c4_fractions(source,abund,names,year) for year in (2014,2017)],
        "track_vs_interior_fixed_2014_2017_site_cohorts":[],
        "station_vs_between_frozen_geographic_matched_pairs":[],
        "limitations":[
            "EXPLORATORY class annotation after previous warm EIVE result, not prospective physiological experiment",
            "One Kashmir railway study, spatially dependent sites, same source repeated 2014/2017",
            "Europe-specific EIVE transferred to Kashmir; not physical temperature or thermal tolerance",
            "Named C4 candidate subset only: unlabelled species not automatically C3",
            "Original abundance data from quadrats have source numerical units not fully resolved",
            "No causal control for differences in disturbance, ballast heat, abundance and station/between landscape",
            "No original source plot by species, original GPS or site ID is output"
        ]}
    for year,gate,match in PANELS:
        if year==2014 or gate.startswith("sparse"):
            out["track_vs_interior_fixed_2014_2017_site_cohorts"].append(
                true_rail_site_cohorts(source,abund,names,eive,year,gate))
            out["station_vs_between_frozen_geographic_matched_pairs"].append(
                four_zone_matched_cohort(source,abund,names,eive,meta,year,gate,match))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--kashmir-original",required=True)
    ap.add_argument("--eive-original",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    out=analyze(args.kashmir_original,args.eive_original)
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    compact={"version":out["version"],"taxonomy":out["candidate_taxonomy"],
             "site_C4_fractions":out["same_original_site_C4_fractions"],
             "source_site_EIVE":out["track_vs_interior_fixed_2014_2017_site_cohorts"],
             "four_site_matched_EIVE":out["station_vs_between_frozen_geographic_matched_pairs"]}
    print(json.dumps(compact,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

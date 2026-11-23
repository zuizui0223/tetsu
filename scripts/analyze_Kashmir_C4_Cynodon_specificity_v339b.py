#!/usr/bin/env python3
"""Post-result exploratory check v339b: C4-guild vs one dominant grass species.

This is deliberately post hoc after v339 showed strong C4 attenuation.
Tests one true original railway same-site cohort (2014 31 sites),
a strict EIVE 2014 20-site cohort and frozen original 14 matched pairs.
Within-Poaceae 6-of-12 exact label permutations preserve all field records.
No site IDs or original species-by-plot rows are exported.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
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
from analyze_Kashmir_C4_mechanism_v339 import (
    C4, C4_GRASS, C3_GRASS, source_fraction, mean_traits,
    mean_safe, four_zone_delta
)

CY="cynodon dactylon"
GRASS_UNIVERSE=tuple(sorted(C4_GRASS|C3_GRASS))
SITE_ZONES=("Track","Interior Plot 2")
YEAR=2014


def comparable_original_sites(presence,abund,names,year):
    sites=[]
    for site in sorted({k[0] for k in presence}):
        p=abund.get((site,year,SITE_ZONES[0]))
        q=abund.get((site,year,SITE_ZONES[1]))
        if p and q:
            sites.append((p,q))
    return sites


def mean_fraction_delta(pairs,names,group,weighted=False):
    return mean(source_fraction(p,names,group,weighted)
                -source_fraction(q,names,group,weighted) for p,q in pairs)


def exact_grass_pathway_label_shuffle(pairs,names,c4group,c3group,weighted=False):
    """Randomly assign 6/12 known Poaceae labels to C4; not a causal null."""
    universe=tuple(sorted(c4group|c3group))
    if set(c4group)&set(c3group) or len(set(universe))!=len(universe):
        raise ValueError("Overlapping source C3 and C4 grass labels")
    obs=mean_fraction_delta(pairs,names,c4group,weighted)
    null=[mean_fraction_delta(pairs,names,set(group),weighted)
          for group in itertools.combinations(universe,len(c4group))]
    mu=mean(null)
    centered=sum(abs(a-mu)>=abs(obs-mu)-1e-12 for a in null)/len(null)
    upper=sum(a>=obs-1e-12 for a in null)/len(null)
    rank=sum(a<=obs+1e-12 for a in null)
    return {"n_source_pairs":len(pairs),
        "n_C4_grass_labels":len(c4group),
        "n_C3_grass_labels":len(c3group),
        "n_exact_taxon_label_assignments":len(null),
        "observed_C4_grass_track_minus_interior_fraction":obs,
        "null_mean":mu,"null_min":min(null),"null_max":max(null),
        "two_sided_exact_label_shuffle_fraction_centered":centered,
        "positive_tail_fraction":upper,
        "number_null_no_more_than_observed":rank,
        "no_field_randomization_and_no_species_prevalence_matching":True}


def qualifying_sites(presence,names,eive,year,systems,min_n,coverage):
    selected=[]
    for site in sorted({k[0] for k in presence}):
        d,qa=site_effect(presence,names,eive,site,year,SITE_ZONES,
                         "position",systems,min_n,coverage)
        if d is not None:
            selected.append((site,d))
    return selected


def single_taxon_specificity(presence,abund,names,eive,year,trait_gate,match_gate):
    systems,n_min,cov=TRAIT_GATES[trait_gate]
    elig=qualifying_sites(presence,names,eive,year,systems,n_min,cov)
    geo,_=site_metadata(global_source_path)
    fixed,_,qa=cohort(presence,names,eive,geo,year,trait_gate,match_gate)
    exclusions=(("none",frozenset()),("Cynodon_only",frozenset((CY,))),
                ("other_C4_except_Cynodon",C4-{CY}),("all_strict_C4",C4))
    site_stats={}
    for weighting in ("presence","original_species_abundance"):
        rows=[]
        for site,original_effect in elig:
            p,q=abund.get((site,year,"Track")),abund.get((site,year,"Interior Plot 2"))
            pairs={}
            for lab,exc in exclusions:
                a=mean_traits(p,names,eive,systems,weighting,exc)
                b=mean_traits(q,names,eive,systems,weighting,exc)
                pairs[lab]=None if a is None or b is None else a-b
            if not math.isclose(pairs["none"],original_effect,rel_tol=0,abs_tol=1e-9) and weighting=="presence":
                raise AssertionError("Original source paired CWM mismatch")
            rows.append(pairs)
        site_stats[weighting]={
            lab:{"n_complete":sum(x[lab] is not None for x in rows),
                 "mean_effect":mean_safe([x[lab] for x in rows if x[lab] is not None])}
             for lab,_ in exclusions}
    pair_stats={}
    for weighting in ("presence","original_species_abundance"):
        out={}
        for lab,exc in exclusions:
            values=[four_zone_delta(abund,m["between_source_site"],
                       m["station_source_site"],year,names,eive,systems,weighting,exc)
                    for m in fixed]
            valid=[x for x in values if x is not None]
            out[lab]={"n_source_matched_pairs":len(valid),
                      "mean_4zone_EIVE":mean_safe(valid)}
        original_scheme=("presence" if weighting=="presence" else
                         "original_positive_source_value")
        old=group_comparison(abund,fixed,names,eive,year,trait_gate,original_scheme,"position")
        if not math.isclose(out["none"]["mean_4zone_EIVE"],
                            old["mean_between_minus_station_track_relative_interior"],
                            abs_tol=1e-9,rel_tol=0):
            raise AssertionError("Original v334 4-plot matched EIVE mismatch")
        pair_stats[weighting]=out
    return {"year":year,"trait_gate":trait_gate,"match_gate":match_gate,
           "source_eligible_sites":len(elig),
           "source_matched_station_between_pairs":len(fixed),
           "site_CWM":site_stats,"matched_4zone_CWM":pair_stats}


def analyze(original,eivepath):
    global global_source_path
    global_source_path=original
    presence,qa=parse_xlsx_matrix(Path(original))
    names=species_names_from_publisher(original)
    eive,eiveqa=eive_table(eivepath)
    abund,abqa=read_original_plant_values(original,presence)
    if CY not in C4_GRASS or len(GRASS_UNIVERSE)!=12:
        raise AssertionError("Original C3/C4 Poaceae taxonomy unexpected")
    cohorts=[(y,comparable_original_sites(presence,abund,names,y)) for y in (2014,2017)]
    result={"version":"v339b-post-outcome-Cynodon-specificity",
            "original_source_hash":qa["original_file_sha256"],
            "EIVE_source_md5":eiveqa["EIVE_original_md5"],
            "source_taxon_labels":list(GRASS_UNIVERSE),
            "source_C4_specificity":[]}
    for year,pairs in cohorts:
        a={"source_year":year,"n_true_field_site_track_interior_pairs":len(pairs),
           "species_presence_fraction":{
               "strict_C4_all":mean_fraction_delta(pairs,names,C4),
               "C4_excluding_Cynodon":mean_fraction_delta(pairs,names,C4-{CY}),
               "Cynodon_only":mean_fraction_delta(pairs,names,frozenset((CY,))),
               "C4_grasses_only":mean_fraction_delta(pairs,names,C4_GRASS),
               "C4_grasses_excluding_Cynodon":mean_fraction_delta(pairs,names,C4_GRASS-{CY}),
               "known_C3_grasses":mean_fraction_delta(pairs,names,C3_GRASS)
           },
           "source_abundance_fraction":{
               "strict_C4_all":mean_fraction_delta(pairs,names,C4,True),
               "C4_excluding_Cynodon":mean_fraction_delta(pairs,names,C4-{CY},True),
               "Cynodon_only":mean_fraction_delta(pairs,names,frozenset((CY,)),True),
               "C4_grasses_only":mean_fraction_delta(pairs,names,C4_GRASS,True),
               "C4_grasses_excluding_Cynodon":mean_fraction_delta(pairs,names,C4_GRASS-{CY},True),
               "known_C3_grasses":mean_fraction_delta(pairs,names,C3_GRASS,True)
           },
           "within_grass_exact_6_of_12_shuffle":exact_grass_pathway_label_shuffle(
               pairs,names,C4_GRASS,C3_GRASS,False),
           "within_grass_exact_5_of_11_excluding_Cynodon":exact_grass_pathway_label_shuffle(
               pairs,names,C4_GRASS-{CY},C3_GRASS,False),
           "within_grass_abundance_label_shuffle":exact_grass_pathway_label_shuffle(
               pairs,names,C4_GRASS,C3_GRASS,True)
        }
        result["source_C4_specificity"].append(a)
    result["Cynodon_C4_attribution_EIVE"]=[
       single_taxon_specificity(presence,abund,names,eive,2014,
          "strict_3systems_5species_50pct","G1_geographic_altitude_primary"),
       single_taxon_specificity(presence,abund,names,eive,2014,
          "relaxed_1system_3species_30pct","G1_geographic_altitude_primary")
    ]
    result["caveats"]=[
       "This particular taxon-specificity extension was designed AFTER seeing v339; exploratory",
       "Randomly relabeling 12 source Poaceae is a taxon-level technical null, not a field randomization",
       "Species prevalence and warm-indicator value are not held constant across grass label shuffles",
       "All original field sites are from one observational railway corridor",
       "Source abundance units uncertain and do not establish biomass",
       "Cynodon dactylon taxon identity and its source distribution can explain most guild effect without demonstrating C4 physiology as causal factor",
       "No GPS, source original site species membership, abundance original plot rows are exported"]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    x=p.parse_args()
    r=analyze(x.kashmir_original,x.eive_original)
    path=Path(x.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

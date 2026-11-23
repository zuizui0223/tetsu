#!/usr/bin/env python3
"""v338: Kashmir authentic species-conditioned quantitative-abundance null.

Hold 2014/2017 species-by-original-plot incidence and EIVE values fixed.
Permute genuine source *species abundance* numbers only among plots in
which each SAME species actually occurs, either globally over the same
year or more tightly *within the same original railway habitat zone*.
Those two nulls preserve species-specific source abundance distributions
and are stricter than v337's arbitrary within-plot species-weight shuffle.
They do not preserve each plot's total abundance and are not causal field
randomization tests. No field species×site rows or coordinates exported.
"""
import argparse
import json
import math
import random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_geographic_station_between_EIVE_v332 import site_metadata,TRAIT_GATES
from analyze_Kashmir_cover_weighted_matched_temperature_v334 import (
    read_original_plant_values,cohort,group_comparison
)

PANEL_DESIGN=(
  (2014,"relaxed_1system_3species_30pct","G1_geographic_altitude_primary"),
  (2014,"strict_3systems_5species_50pct","G1_geographic_altitude_primary"),
  (2014,"relaxed_1system_3species_30pct","G3_geo_altitude_richness_coverage"),
  (2017,"sparse_3systems_1species","G1_geographic_altitude_primary")
)
ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
METRIC="position"
REFERENCE="Interior Plot 2"
MODES=("species_all_zones","species_within_original_zone")
N_DRAWS=1999
SEED=338

def site_trait_abundances(source_abundances,names,eive,year,systems):
    world={}
    for (site,y,zone),ab in sorted(source_abundances.items()):
        if y!=year or ab is None:continue
        v={}
        for code,amount in ab.items():
            name=names[code]
            if name is None or name not in eive or eive[name]["systems"]<systems:
                continue
            if amount<=0 or not math.isfinite(amount):raise ValueError("Invalid source abundance")
            if name in v:raise ValueError("Source same accepted name duplicated")
            v[name]=float(amount)
        if v:world[site,y,zone]=v
    return world

def matched_keys(matched,year):
    out=[]
    for row in matched:
        bs=row["between_source_site"]
        st=row["station_source_site"]
        groups={"BT":(bs,year,"Track"),"BI":(bs,year,REFERENCE),
                "ST":(st,year,"Track"),"SI":(st,year,REFERENCE)}
        out.append(groups)
    return out

def cwm(plot,eive):
    if not plot:raise ValueError("Empty source plot in matched cohort")
    temp={sp:float(eive[sp][METRIC]) for sp in plot}
    unweighted=mean(temp.values())
    weighted=sum(a*temp[sp] for sp,a in plot.items())/sum(plot.values())
    return unweighted,weighted,weighted-unweighted

def matched_effect(world,keys,eive):
    if not keys:raise ValueError("Original matched site pairs absent")
    deltas=[]
    for record in keys:
        c={}
        for lab,key in record.items():
            if key not in world:raise ValueError("Missing source matched location")
            c[lab]=cwm(world[key],eive)
        sign={"BT":1,"BI":-1,"ST":-1,"SI":1}
        oo=sum(sign[s]*c[s][0] for s in sign)
        aa=sum(sign[s]*c[s][1] for s in sign)
        dd=sum(sign[s]*c[s][2] for s in sign)
        if not math.isclose(aa-oo,dd,rel_tol=0,abs_tol=1e-10):
            raise AssertionError("Source abundance incremental CWM identity failed")
        deltas.append((oo,aa,dd))
    return {"matched_pairs":len(keys),
        "presence_four_zone_CWM":mean(r[0] for r in deltas),
        "abundance_four_zone_CWM":mean(r[1] for r in deltas),
        "abundance_minus_presence_increment":mean(r[2] for r in deltas),
        "number_positive_incremental_pairs":sum(r[2]>0 for r in deltas)}

def shuffle_setup(world,mode):
    groups=defaultdict(list)
    for key,values in sorted(world.items()):
        for species,ab in sorted(values.items()):
            grp=(species,) if mode=="species_all_zones" else (species,key[2])
            groups[grp].append((key,species,ab))
    movable=[x for x in groups.values() if len(x)>=2 and len({r[2] for r in x})>=2]
    species={g[0] for g in groups}
    movable_species={r[1] for v in movable for r in v}
    qa={"original_source_species_groups":len(groups),
        "source_matched_eive_species":len(species),
        "groups_with_two_or_more_occurrences":sum(len(v)>=2 for v in groups.values()),
        "groups_with_two_or_more_different_abundance_values":len(movable),
        "source_species_in_nontrivially_permutable_groups":len(movable_species),
        "original_source_species_abundance_records":sum(len(v) for v in groups.values()),
        "source_abundance_records_in_nontrivially_permutable_groups":sum(len(v) for v in movable),
        "unmoved_group_values_retained_as_source_original":True}
    return groups,movable,qa

def reshuffle_world(world,movable,rng):
    output={key:dict(v) for key,v in world.items()}
    for records in movable:
        vals=[r[2] for r in records]
        rng.shuffle(vals)
        for (key,species,_),num in zip(records,vals):
            output[key][species]=num
    return output

def technical_null(world,keys,eive,mode,draws=N_DRAWS,seed=SEED):
    groups,movable,qa=shuffle_setup(world,mode)
    obs=matched_effect(world,keys,eive)
    rng=random.Random(seed)
    sample=[]
    for _ in range(draws):
        assigned=reshuffle_world(world,movable,rng)
        stat=matched_effect(assigned,keys,eive)
        # Species membership, source plots and site-year-zone identity fixed.
        if not math.isclose(stat["presence_four_zone_CWM"],
                            obs["presence_four_zone_CWM"],
                            rel_tol=0,abs_tol=1e-10):
            raise AssertionError("Species-conditioned source null changed plot taxon occurrence")
        sample.append(stat["abundance_minus_presence_increment"])
    sample.sort()
    null_mu=mean(sample)
    val=obs["abundance_minus_presence_increment"]
    count=sum(abs(x-null_mu)>=abs(val-null_mu)-1e-12 for x in sample)
    low=sum(x<=val+1e-12 for x in sample)
    return {"permutation_mode":mode,
        "source_group_permutability_QA":qa,
        "draws":draws,"seed":seed,
        "observed":obs,
        "source_species_conditional_null_mean":null_mu,
        "technical_null_95pct":[sample[int(.025*draws)],sample[int(.975*draws)]],
        "technical_two_sided_fraction_centered_on_null_mean_plus_one":
            (count+1)/(draws+1),
        "technical_lower_tail_fraction_plus_one":(low+1)/(draws+1),
        "technical_abundance_label_permutation_not_field_randomization":True}

def analyze(original,eivepath):
    plotdata,plotqa=parse_xlsx_matrix(Path(original))
    names=species_names_from_publisher(original)
    eive,eiveqa=eive_table(eivepath)
    geo,geoqa=site_metadata(original)
    abundance,abundqa=read_original_plant_values(original,plotdata)
    out={"version":"v338",
        "status":"original-source-species-conditioned-abundance-null",
        "original_kashmir_sha256":plotqa["original_file_sha256"],
        "EIVE_source_original_md5":eiveqa["EIVE_original_md5"],
        "species_abundance_origin":"Rashid et al. 2021 section 2.2 sampled twelve 1x1m quadrats per source site for species abundance; individual Dryad cell arithmetic units not fully specified",
        "source_abundance_QA":abundqa,"panels":[]}
    for idx,(year,trait_gate,match_gate) in enumerate(PANEL_DESIGN):
        matched,srcqa,matchqa=cohort(plotdata,names,eive,geo,year,
                                    trait_gate,match_gate)
        panel={"original_year":year,"source_EIVE_coverage_gate":trait_gate,
               "source_matching_gate":match_gate,
               "source_eligibility":srcqa,"source_original_match_QA":matchqa}
        if not matched:
            panel["estimable"]=False
            out["panels"].append(panel);continue
        systems=TRAIT_GATES[trait_gate][0]
        world=site_trait_abundances(abundance,names,eive,year,systems)
        keys=matched_keys(matched,year)
        observed=matched_effect(world,keys,eive)
        rawcheck=group_comparison(abundance,matched,names,eive,year,trait_gate,
                                   "original_positive_source_value",METRIC)
        prescheck=group_comparison(abundance,matched,names,eive,year,trait_gate,
                                   "presence",METRIC)
        if not math.isclose(observed["abundance_four_zone_CWM"],
                            rawcheck["mean_between_minus_station_track_relative_interior"],
                            abs_tol=1e-10,rel_tol=0):
            raise AssertionError("v334 original-source numeric abundance CWM differs in v338")
        if not math.isclose(observed["presence_four_zone_CWM"],
                            prescheck["mean_between_minus_station_track_relative_interior"],
                            abs_tol=1e-10,rel_tol=0):
            raise AssertionError("v334 original-source species occurrence CWM differs in v338")
        panel.update({"estimable":True,"n_original_matched_station_between_pairs":len(keys),
                      "original_source_observed":observed,
                      "conditional_nulls":[
                         technical_null(world,keys,eive,mode,N_DRAWS,
                                        SEED+100*idx+(0 if mode==MODES[0] else 19))
                         for mode in MODES]})
        out["panels"].append(panel)
    out["limitations"]=[
       "v338 was designed after v337 observed low technical conditional fraction and is exploratory, not registered-confirmatory",
       "Species-conditioned null conserves each species original year-specific distribution of abundance values, but NOT each plot total abundance",
       "Tight species×zone null additionally conserves source Track/interior habitat-specific species abundances; station/between labels are observational",
       "Species with single occurrences or invariant source abundance cannot permute; counts of movable source records are reported",
       "EIVE is a Europe-derived thermal ecology indicator not a measured species physiological threshold in Kashmir",
       "Source numbers refer to quadrat species abundance but original workbook cell units still not fully identified from sampling methods",
       "No temperature sensors, randomized railway controls, replicated climatic gradients or direct train dispersal data",
       "No source site coordinates or species-by-site field plot records included in output"]
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze(a.kashmir_original,a.eive_original)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v338","results":[
        {"year":p["original_year"],"source_trait_gate":p["source_EIVE_coverage_gate"],
         "match_gate":p["source_matching_gate"],
         "n":p["source_original_match_QA"]["matched_count"],
         **({"observed":p["original_source_observed"],
             "conditional_nulls":p["conditional_nulls"]} if p.get("estimable") else {})}
        for p in out["panels"]]},ensure_ascii=False),flush=True)
if __name__=="__main__":main()

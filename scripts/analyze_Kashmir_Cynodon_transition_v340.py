#!/usr/bin/env python3
"""v340: genuine same-site 2014->2017 Cynodon recorded persistence by railway habitat.

Source-validated Kashmir Dryad 99735 original 31 sites, 4 habitat zones, 2 years.
A recorded transition is not individual survival, demographic establishment,
or an experiment in train dispersal. No original per-site fields are emitted.
Protocol: docs/tetsu_v340_Cynodon_2014_2017_transition_frozen_design_ja.md
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import Counter
from pathlib import Path
from statistics import mean, median

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Kashmir_cover_weighted_matched_temperature_v334 import read_original_plant_values
from analyze_Kashmir_C4_mechanism_v339 import C4, C3_GRASS

CYNODON="cynodon dactylon"
PRIMARY=("Track","Interior Plot 2")
ALTERNATE=("Track","Interior Plot 1")
ZONES=PRIMARY
YEARS=(2014,2017)
GROUPS={
    "Cynodon_dactylon": frozenset((CYNODON,)),
    "other_C4_candidates_ex_Cynodon": C4-frozenset((CYNODON,)),
    "all_strict_C4_candidates": C4,
    "C3_Poaceae_reference_6species": C3_GRASS,
}
N_FLIPS=19999
N_BOOTSTRAP=4999


def quantile(values,p):
    if not values:
        return None
    v=sorted(values)
    q=(len(v)-1)*p
    a=int(math.floor(q))
    b=int(math.ceil(q))
    return v[a]+(v[b]-v[a])*(q-a)


def source_sites(abund,zones):
    names=sorted({site for site,year,zone in abund})
    eligible=[]
    exclusions=Counter()
    for site in names:
        if any(abund.get((site,y,z)) is None for y in YEARS for z in zones):
            exclusions["missing_at_least_one_original_site_year_zone"]+=1
        else:
            eligible.append(site)
    if len(names)!=31:
        raise ValueError("Original Kashmir study no longer has 31 spatial sites")
    return eligible,dict(exclusions)


def recorded_species(positive_record,species_names,nameset):
    if positive_record is None:
        raise ValueError("Do not infer absence from source-incomplete record")
    return any(species_names[i] in nameset for i in positive_record)


def transition(a,b):
    return str(int(bool(a)))+str(int(bool(b)))


def zone_counts(rows,zone,species_names,nameset):
    t=Counter()
    for site,field in rows:
        t[transition(recorded_species(field[2014,zone],species_names,nameset),
                     recorded_species(field[2017,zone],species_names,nameset))]+=1
    for k in ("00","01","10","11"):
        t.setdefault(k,0)
    n=sum(t.values())
    return {"source_complete_sites":n,
            "transition_2014_to_2017":{k:t[k] for k in ("00","01","10","11")},
            "2014_occupied_sites":t["10"]+t["11"],
            "2017_occupied_sites":t["01"]+t["11"],
            "2014_presence_to_2017_reobserved_fraction":
               t["11"]/(t["10"]+t["11"]) if (t["10"]+t["11"]) else None,
            "2014_absence_to_2017_new_record_fraction":
               t["01"]/(t["00"]+t["01"]) if (t["00"]+t["01"]) else None,
            "2014_to_2017_occupancy_change_fraction":
               (t["01"]-t["10"])/n if n else None}


def exact_mcnemar_double_baseline(rows,zones,species_names,nameset):
    both=[]
    for site,field in rows:
        if (recorded_species(field[2014,zones[0]],species_names,nameset)
              and recorded_species(field[2014,zones[1]],species_names,nameset)):
            both.append((
                recorded_species(field[2017,zones[0]],species_names,nameset),
                recorded_species(field[2017,zones[1]],species_names,nameset)
            ))
    n10=sum(a and not b for a,b in both)
    n01=sum(b and not a for a,b in both)
    n11=sum(a and b for a,b in both)
    n00=sum(not a and not b for a,b in both)
    n=n10+n01
    p=min(1.0,2*sum(math.comb(n,k) for k in range(min(n10,n01)+1))/(2**n)) if n else None
    return {"n_original_sites_occupied_both_habitats_in_2014":len(both),
            "2017_both_reobserved":n11,
            "2017_track_only_reobserved":n10,
            "2017_interior_only_reobserved":n01,
            "2017_neither_reobserved":n00,
            "exact_McNemar_technical_two_sided":p,
            "not_estimable_direct_zone_difference_if_both2014_positive_sites_under5":len(both)<5,
            "not_independent_individual_survival":True}


def technical_site_did(rows,zones,species_names,nameset,draws=N_FLIPS,boot=N_BOOTSTRAP):
    observed_d=[]
    for site,field in rows:
        vals={(y,z):int(recorded_species(field[y,z],species_names,nameset))
              for y in YEARS for z in zones}
        difference=(vals[2017,zones[0]]-vals[2014,zones[0]]
                    -vals[2017,zones[1]]+vals[2014,zones[1]])
        observed_d.append(difference)
    if not observed_d:
        return None
    d=mean(observed_d)
    rng=random.Random(340)
    flips=[]
    for _ in range(draws):
        flips.append(sum(x*(1 if rng.getrandbits(1) else -1)
                         for x in observed_d)/len(observed_d))
    p=(1+sum(abs(x)>=abs(d)-1e-12 for x in flips))/(draws+1)
    rngboot=random.Random(341)
    bootvals=[mean([observed_d[rngboot.randrange(len(observed_d))]
                    for _ in observed_d]) for _ in range(boot)]
    return {"n_same_original_sites":len(observed_d),
            "observed_change_in_track_minus_interior_occupancy_fraction":d,
            "n_positive_site_changes":sum(x>0 for x in observed_d),
            "n_negative_site_changes":sum(x<0 for x in observed_d),
            "n_no_relative_change":sum(x==0 for x in observed_d),
            "site_swap_signflip_draws":draws,
            "technical_signflip_two_sided_fraction_plus_one":p,
            "site_cluster_resampling_draws":boot,
            "site_resampling_95pct_percentiles_not_independent_field_CI":
               [quantile(bootvals,.025),quantile(bootvals,.975)],
            "inference_is_observational_spatially_dependent":True}


def source_species_amounts(rows,zones,species_names):
    idx=[i for i,sp in enumerate(species_names) if sp==CYNODON]
    if len(idx)!=1:raise ValueError("Cynodon must be one exact original source species column")
    taxon_index=idx[0]
    result={}
    for zone in zones:
        pos_by_year={y:[] for y in YEARS}
        share_by_year={y:[] for y in YEARS}
        richness_by_year={y:[] for y in YEARS}
        retained_changes=[]
        retained_log=[]
        present_both=0
        for _,fields in rows:
            original2014=fields[2014,zone]
            original2017=fields[2017,zone]
            v0=original2014.get(taxon_index)
            v1=original2017.get(taxon_index)
            if v0 is not None and v1 is not None:
                retained_changes.append(v1-v0)
                retained_log.append(math.log(v1/v0))
                present_both+=1
            for yr,rec,v in ((2014,original2014,v0),(2017,original2017,v1)):
                richness_by_year[yr].append(len(rec))
                if v is not None:
                    pos_by_year[yr].append(v)
                share_by_year[yr].append((v or 0.0)/sum(rec.values()))
        result[zone]={
           "n_source_complete_sites":len(rows),
           "year_stats":{
              str(y):{
                  "n_positive_Cynodon_source_values":len(pos_by_year[y]),
                  "mean_Cynodon_source_value_among_recorded_positive":
                      mean(pos_by_year[y]) if pos_by_year[y] else None,
                  "median_Cynodon_source_value_among_recorded_positive":
                      median(pos_by_year[y]) if pos_by_year[y] else None,
                  "mean_Cynodon_share_of_plot_recorded_source_abundance":
                      mean(share_by_year[y]) if share_by_year[y] else None,
                  "mean_recorded_community_richness":mean(richness_by_year[y])}
              for y in YEARS
           },
           "same_source_plot_Cynodon_positive_both_years":{
              "n":present_both,
              "mean_log_2017_div_2014":mean(retained_log) if retained_log else None,
              "median_log_2017_div_2014":median(retained_log) if retained_log else None,
              "geometric_mean_ratio_2017_div_2014":(
                   math.exp(mean(retained_log)) if retained_log else None),
              "median_raw_numeric_difference":(
                   median(retained_changes) if retained_changes else None),
              "n_source_values_decreased":sum(x<0 for x in retained_changes),
              "n_source_values_increased":sum(x>0 for x in retained_changes),
              "n_source_values_unchanged":sum(x==0 for x in retained_changes),
              "absence_cells_not_used_as_zero_abundance":True,
           }
        }
    return result


def panel(abund,names,zones):
    cohort,exclusions=source_sites(abund,zones)
    rows=[(site,{(y,z):abund[site,y,z] for y in YEARS for z in zones})
          for site in cohort]
    groups={}
    for class_name,species_set in GROUPS.items():
        stats={}
        for layer,sub in (("all",rows),
                          ("station",[x for x in rows if not x[0].startswith("Between ")]),
                          ("between",[x for x in rows if x[0].startswith("Between ")])):
            zone0=zone_counts(sub,zones[0],names,species_set)
            zone1=zone_counts(sub,zones[1],names,species_set)
            stats[layer]={
               "n_original_sites":len(sub),
               "zones":{zones[0]:zone0,zones[1]:zone1},
               "both_2014_positive_habitats":
                   exact_mcnemar_double_baseline(sub,zones,names,species_set),
               "change_in_zone_occupancy_difference":
                   technical_site_did(sub,zones,names,species_set,
                      draws=N_FLIPS if layer=="all" else 1999,boot=N_BOOTSTRAP if layer=="all" else 999)
            }
        groups[class_name]=stats
    return {"source_zone_comparison":list(zones),
            "number_source_original_sites":31,
            "number_common_complete_source_sites":len(rows),
            "excluded_sites_for_missing_in_any_four_source_zone_year_cells":exclusions,
            "complete_original_stations":sum(not x[0].startswith("Between ") for x in rows),
            "complete_original_between_sites":sum(x[0].startswith("Between ") for x in rows),
            "species_class_has_at_least_one_observed_taxon_not_individual_survival":True,
            "group_transitions":groups,
            "Cynodon_positive_source_abundance":source_species_amounts(rows,zones,names)}


def analyze(source):
    x,source_qa=parse_xlsx_matrix(Path(source))
    names=species_names_from_publisher(source)
    abund,abundance_qa=read_original_plant_values(source,x)
    if CYNODON not in names:
        raise ValueError("Exact genuine source Cynodon dactylon taxon missing")
    if not C4<=set(n for n in names if n) or not C3_GRASS<=set(n for n in names if n):
        raise ValueError("Previously verified C4/C3 candidate taxonomy changed")
    out={
      "version":"v340-genuine-Cynodon-temporal-reobservation",
      "source_Dryad_SHA256":source_qa["original_file_sha256"],
      "original_31_sites_and_2_source_years":source_qa,
      "species_abundance_source_QA":abundance_qa,
      "original_spatial_site_repeated_not_two_independent_studies":True,
      "primary":panel(abund,names,PRIMARY),
      "nearby_interior_plot1_distance_sensitivity":panel(abund,names,ALTERNATE),
      "limitations":[
        "Recorded 2014 to 2017 presence transition is neither individual plant survival nor verified population establishment",
        "Missing original field records not coerced to absence or numeric zero",
        "Positive original source species-abundance units not fully defined; do not claim percent cover, biomass or plant density",
        "Study is one Kashmir railway, 31 originally surveyed locations; no independent route replication",
        "All significance fractions are technical label and site-resampling diagnostics not randomized assignment",
        "Station and between are not true observed train frequency measures",
        "No microclimate, ballast temperature, train dispersal, seed rain or flowering measured",
        "Plot-level source species abundance, site names and raw GPS never exported",
      ]
    }
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--kashmir-original",required=True)
    a.add_argument("--output",required=True)
    args=a.parse_args()
    result=analyze(args.kashmir_original)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

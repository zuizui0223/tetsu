#!/usr/bin/env python3
"""v335: genuine Kashmir source-value CWM sign-reversal robustness.

Same 1-to-1 station/between pairs as source-verified v332/v334. Source positive
numeric values remain unnamed units (no cover percent/biomass assumption).
Result-selected exponent/cap search is NOT allowed: panels fixed in v335 frozen
design and labelled post hoc diagnostics.
"""
import argparse,json,math
from collections import Counter
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_geographic_station_between_EIVE_v332 import site_metadata,TRAIT_GATES,exact_signflip
from analyze_Kashmir_cover_weighted_matched_temperature_v334 import (
    read_original_plant_values,cohort,group_comparison
)

PANELS=(
 (2014,"relaxed_1system_3species_30pct","G1_geographic_altitude_primary"),
 (2014,"strict_3systems_5species_50pct","G1_geographic_altitude_primary"),
 (2014,"relaxed_1system_3species_30pct","G3_geo_altitude_richness_coverage"),
 (2017,"sparse_3systems_1species","G1_geographic_altitude_primary")
)
EXPONENTS=(0.0,0.25,0.5,0.75,1.0)
CAPS_AT_ONE=(None,100,20,10,5)
REFERENCE="Interior Plot 2"

def accepted_plot(record,names,eive,systems):
    if record is None:raise ValueError("Matched original site-year plot was missing")
    taxa={}
    for code,v in record.items():
        name=names[code]
        if name in eive and eive[name]["systems"]>=systems:
            if name in taxa:raise ValueError("Duplicate source EIVE name within original plot")
            taxa[name]=v
    if not taxa:raise ValueError("No original source numeric EIVE known species")
    return taxa

def raw_matched_plot_values(cover,matched,names,eive,year,systems):
    output=[]
    for m in matched:
        groups={}
        b=m["between_source_site"];s=m["station_source_site"]
        for key,site,zone in (
                ("BT",b,"Track"),("BI",b,REFERENCE),
                ("ST",s,"Track"),("SI",s,REFERENCE)):
            groups[key]=accepted_plot(cover.get((site,year,zone)),names,eive,systems)
        output.append(groups)
    return output

def community_mean(src,eive,metric,power,cap,excluded=None):
    selected={s:w for s,w in src.items() if s!=excluded}
    if not selected:return None
    weight={}
    for s,v in selected.items():
        if v<=0:raise ValueError("Source original recorded plant weight not positive")
        w=min(v,cap) if cap is not None else v
        w=w**power
        if not w>0 or not math.isfinite(w):
            raise ValueError("Source weighting not finite")
        weight[s]=w
    denom=sum(weight.values())
    return sum((w/denom)*eive[s][metric] for s,w in weight.items())

def effect(pairs,eive,metric,power,cap,excluded=None):
    values=[]
    for p in pairs:
        cwm={k:community_mean(arr,eive,metric,power,cap,excluded)
             for k,arr in p.items()}
        if any(v is None for v in cwm.values()):continue
        delta=(cwm["BT"]-cwm["BI"]-cwm["ST"]+cwm["SI"])
        values.append((delta,cwm["BT"]-cwm["ST"],cwm["BI"]-cwm["SI"]))
    return values

def summary(vals):
    if not vals:return {"n":0,"mean":None,"min":None,"max":None}
    return {"n":len(vals),"mean":mean(vals),"min":min(vals),"max":max(vals)}

def weighted_panel(pairs,eive,metric,power,cap):
    rows=effect(pairs,eive,metric,power,cap)
    values=[v[0] for v in rows]
    return {"power":power,"cap":cap,"matched_pairs":len(rows),
            "mean_between_minus_station_track_relative_interior":mean(values) if values else None,
            "positive_pairs":sum(x>0 for x in values),
            "negative_pairs":sum(x<0 for x in values),
            "mean_between_minus_station_track_only":mean([v[1] for v in rows]) if rows else None,
            "mean_between_minus_station_interior_only":mean([v[2] for v in rows]) if rows else None,
            "matched_sign_flip_exploratory":exact_signflip(values),
            "NOT_physical_temperature_or_cover_percentage":True}

def single_taxon_impact(pairs,eive,metric="position"):
    base=effect(pairs,eive,metric,1.0,None)
    baseline=mean(x[0] for x in base)
    taxa=sorted(set().union(*(set().union(*item.values()) for item in pairs)))
    out=[]
    for taxon in taxa:
        rows=effect(pairs,eive,metric,1.0,None,excluded=taxon)
        if not rows:continue
        value=mean(x[0] for x in rows)
        out.append({"taxon":taxon,"remaining_pairs":len(rows),
                    "mean_after_one_taxon_removed":value,
                    "change_on_exclusion":value-baseline})
    same=[r for r in out if r["remaining_pairs"]==len(pairs)]
    return {"n_source_EIVE_taxa_in_source_matched_plots":len(taxa),
            "n_evaluable_one_taxon_omissions":len(out),
            "n_preserving_all_paired_geographic_sites":len(same),
            "mean_raw_weighted_effect_after_removing_one_taxon_full_pairs":summary(
                [r["mean_after_one_taxon_removed"] for r in same]),
            "number_full_pair_source_taxon_deletions_reversing_effect_sign":sum(
                (r["mean_after_one_taxon_removed"]>0)!=(baseline>0) for r in same),
            "top_five_source_taxa_whose_removal_most_increases_effect":[
                {"taxon":r["taxon"],"delta_effect_after_exclusion":r["change_on_exclusion"]}
                 for r in sorted(same,key=lambda a:-a["change_on_exclusion"])[:5]],
            "top_five_source_taxa_whose_removal_most_decreases_effect":[
                {"taxon":r["taxon"],"delta_effect_after_exclusion":r["change_on_exclusion"]}
                 for r in sorted(same,key=lambda a:a["change_on_exclusion"])[:5]],
            "source_species_NAMES_only_not_site_species_observation_rows":True}

def panel_source_outlier_trace(pairs):
    counts=Counter()
    for p in pairs:
        for zone,arr in p.items():
            counts[zone]+=sum(v>100 for v in arr.values())
    return dict(counts)

def analyze(original,eivepath):
    source,sourceqa=parse_xlsx_matrix(Path(original))
    names=species_names_from_publisher(original)
    eive,eiveqa=eive_table(eivepath)
    geo,geoqa=site_metadata(original)
    cover,coverqa=read_original_plant_values(original,source)
    results=[]
    for year,traitgate,matchgate in PANELS:
        matched,eligibility,matchqa=cohort(source,names,eive,geo,year,traitgate,matchgate)
        item={"year":year,"trait_gate":traitgate,"matching_gate":matchgate,
              "site_pair_eligibility":eligibility,"site_matching":matchqa}
        if not matched:
            item["estimable"]=False;results.append(item);continue
        system=TRAIT_GATES[traitgate][0]
        pairs=raw_matched_plot_values(cover,matched,names,eive,year,system)
        checks=group_comparison(cover,matched,names,eive,year,traitgate,
                               "original_positive_source_value","position")
        weighted=weighted_panel(pairs,eive,"position",1.,None)
        if not math.isclose(weighted["mean_between_minus_station_track_relative_interior"],
                            checks["mean_between_minus_station_track_relative_interior"],
                            rel_tol=0,abs_tol=1e-10):
            raise AssertionError("Original v334 source weighted EIVE baseline changed")
        first=group_comparison(cover,matched,names,eive,year,traitgate,"presence","position")
        presence=weighted_panel(pairs,eive,"position",0.,None)
        if not math.isclose(presence["mean_between_minus_station_track_relative_interior"],
                            first["mean_between_minus_station_track_relative_interior"],
                            rel_tol=0,abs_tol=1e-10):
            raise AssertionError("Original v334 source presence CWM baseline changed")
        half=group_comparison(cover,matched,names,eive,year,traitgate,"sqrt_original_source_value","position")
        sqrt=weighted_panel(pairs,eive,"position",.5,None)
        if not math.isclose(sqrt["mean_between_minus_station_track_relative_interior"],
                            half["mean_between_minus_station_track_relative_interior"],
                            rel_tol=0,abs_tol=1e-10):
            raise AssertionError("Original v334 sqrt source CWM baseline changed")
        exponents=[weighted_panel(pairs,eive,"position",x,None) for x in EXPONENTS]
        caps=[weighted_panel(pairs,eive,"position",1,x) for x in CAPS_AT_ONE]
        outliers=panel_source_outlier_trace(pairs)
        item.update({"estimable":True,
           "source_values_above100_in_matched_communities_no_site_ids":outliers,
           "weight_exponent_gradient":exponents,
           "positive_value_caps":caps,
           "single_source_taxon_omission_direct_weight":single_taxon_impact(pairs,eive)})
        results.append(item)
    return {"version":"v335","original_Dryad_SHA256":sourceqa["original_file_sha256"],
        "independent_EIVE_MD5":eiveqa["EIVE_original_md5"],
        "original_source_cover_QA":coverqa,"panels":results,
        "limitations":[
          "Source positive numeric cover values not verified as percentages, biomass or counts; no interpretation of their physical units",
          "A single original value above100 exists among all source species observations, and its matched-community category is audited",
          "Exponent and weight caps are exploratory technical sensitivities chosen AFTER learning of a sign reversal",
          "Source effect changes when dominance weights change, not a causal observation of thermal stress on abundance",
          "Selected geographic station-between pairs frozen from original v332 and not rematched for favorable effect",
          "All 31 sites in one railway corridor and EIVE calibrated from Europe: no claim of global railway effects",
          "Cannot infer temperature measurements, microclimate, train-carried propagules, species physiology or genetic adaptation"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze(a.kashmir_original,a.eive_original)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v335",
      "source_original_QA":out["original_source_cover_QA"],
      "panels":[{"year":x["year"],"trait_gate":x["trait_gate"],
         "match_gate":x["matching_gate"],"matched":x["site_matching"],
         **({k:x[k] for k in ("source_values_above100_in_matched_communities_no_site_ids",
              "weight_exponent_gradient","positive_value_caps",
              "single_source_taxon_omission_direct_weight")}
           if x.get("estimable") else {})} for x in out["panels"]]},
         ensure_ascii=False),flush=True)

if __name__=="__main__":main()

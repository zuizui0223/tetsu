#!/usr/bin/env python3
"""v337: source-abundance weighted vs occurrence EIVE CWM decomposition in Kashmir.

Original field matrix (true Dryad SHA) and independent EIVE temperature trait
(source MD5) checked by inherited source readers. Prior v332 matched
station–between-site lists are frozen BEFORE examining quantitative results.

Null: permute nonnegative original *positive species abundance values* among
the already observed species within EACH field plot independently. This
preserves all source species identity/plot richness/abundance spectrum/year
and geographic 1:1 pairs, but breaks species temperature–abundance pairing.
No source species-by-site and geographical rows are exported.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_geographic_station_between_EIVE_v332 import (
    site_metadata, TRAIT_GATES
)
from analyze_Kashmir_cover_weighted_matched_temperature_v334 import (
    read_original_plant_values, cohort, group_comparison
)

PANELS=(
  (2014,"relaxed_1system_3species_30pct","G1_geographic_altitude_primary"),
  (2014,"strict_3systems_5species_50pct","G1_geographic_altitude_primary"),
  (2014,"relaxed_1system_3species_30pct","G3_geo_altitude_richness_coverage"),
  (2017,"sparse_3systems_1species","G1_geographic_altitude_primary"),
)
SOURCE_DOI="10.5061/dryad.0p2ngf20h"
PAPER_DOI="10.1111/1365-2664.13961"
EIVE_DOI="10.5281/zenodo.7534792"
REFERENCE="Interior Plot 2"
METRIC="position"
PERM_DRAWS=999
SEED=337

def plot_traits(abundance,accepted_names,eive,minsystems,metric=METRIC):
    if abundance is None:raise ValueError("Original matched plot unavailable")
    out={}
    for index,amount in abundance.items():
        name=accepted_names[index]
        if name is None or name not in eive or eive[name]["systems"]<minsystems:
            continue
        if not math.isfinite(amount) or amount<=0:raise ValueError("Invalid source abundance")
        if name in out:raise ValueError("Duplicate accepted taxon")
        out[name]=(float(eive[name][metric]),float(amount))
    if not out:raise ValueError("Matched plot has no temperature-known taxa")
    return out

def cwm_from_pairs(trait_amount,exponent=1.0,shuffle=False,rng=None):
    if not trait_amount:raise ValueError("Empty true source plot")
    items=sorted(trait_amount)
    temps=[trait_amount[s][0] for s in items]
    values=[trait_amount[s][1] for s in items]
    if shuffle:
        if rng is None:raise ValueError("RNG required for within-source-plot label shuffle")
        rng.shuffle(values)
    if not all(x>0 for x in values):
        raise ValueError("No source positive species abundance")
    base=mean(temps)
    weights=[x**exponent for x in values]
    weighted=sum(t*w for t,w in zip(temps,weights))/sum(weights)
    mu_a=mean(weights)
    cov=mean([(t-base)*(w-mu_a) for t,w in zip(temps,weights)])
    if not math.isclose(weighted-base,cov/mu_a,rel_tol=0,abs_tol=1e-10):
        raise AssertionError("Community mean abundance–EIVE covariance identity failed")
    pp=[w/sum(weights) for w in weights]
    return {"presence_CWM":base,"abundance_CWM":weighted,
        "abundance_alignment":weighted-base,
        "abundance_trait_covariance_over_mean":cov/mu_a,
        "effective_number_matched_species":1/sum(p*p for p in pp),
        "largest_matched_species_share":max(pp),
        "source_taxa_n":len(items),
        "source_total_positive_abundance":sum(values)}

def source_site_pairs(plot_records,paired):
    keys={}
    for a in paired:
        for cls,site in (("B",a["between_source_site"]),("S",a["station_source_site"])):
            for typ,z in (("T","Track"),("I",REFERENCE)):
                k=(site,a["source_year"],z)
                keys[cls+typ]=plot_records.get(k)
    return keys

def matched_source_plots(raw,source_matched,names,eive,year,minsystems):
    output=[]
    for m in source_matched:
        z={}
        for key,site,zone in (
          ("BT",m["between_source_site"],"Track"),
          ("BI",m["between_source_site"],REFERENCE),
          ("ST",m["station_source_site"],"Track"),
          ("SI",m["station_source_site"],REFERENCE)):
            z[key]=plot_traits(raw.get((site,year,zone)),names,eive,minsystems)
        output.append(z)
    return output

def signed_linear(component):
    return component["BT"]-component["BI"]-component["ST"]+component["SI"]

def calc_paired(pairs,exponent,rng=None,permute=False):
    output=[]
    for p in pairs:
        stats={key:cwm_from_pairs(v,exponent,shuffle=permute,rng=rng)
               for key,v in p.items()}
        alignment=signed_linear({k:v["abundance_alignment"] for k,v in stats.items()})
        abundance=signed_linear({k:v["abundance_CWM"] for k,v in stats.items()})
        incidence=signed_linear({k:v["presence_CWM"] for k,v in stats.items()})
        if not math.isclose(abundance-incidence,alignment,rel_tol=0,abs_tol=1e-10):
            raise AssertionError("Between-vs-station four-community decomposition not additive")
        output.append({"four":stats,"relative_presence":incidence,
                       "relative_abundance":abundance,
                       "relative_abundance_alignment":alignment})
    return output

def means_of_pairs(rows):
    if not rows:return None
    keys=("relative_presence","relative_abundance","relative_abundance_alignment")
    obs={k:mean(x[k] for x in rows) for k in keys}
    obs["plot_alignments"]={s:mean(x["four"][s]["abundance_alignment"] for x in rows)
         for s in ("BT","BI","ST","SI")}
    obs["plot_effective_species"]={s:mean(x["four"][s]["effective_number_matched_species"] for x in rows)
         for s in ("BT","BI","ST","SI")}
    obs["plot_max_species_abundance_share"]={s:mean(x["four"][s]["largest_matched_species_share"] for x in rows)
         for s in ("BT","BI","ST","SI")}
    obs["n_matched_station_between_pairs"]=len(rows)
    obs["positive_relative_abundance_alignment_pairs"]=sum(
         x["relative_abundance_alignment"]>0 for x in rows)
    if not math.isclose(obs["relative_abundance"]-obs["relative_presence"],
                        obs["relative_abundance_alignment"],rel_tol=0,abs_tol=1e-10):
        raise AssertionError("Mean four-pair additive identity violated")
    return obs

def conditional_null(pairs,exponent=1.0,draws=PERM_DRAWS,seed=SEED):
    if not pairs:return None
    observed=means_of_pairs(calc_paired(pairs,exponent))
    rng=random.Random(seed)
    keys=("relative_abundance_alignment", "relative_abundance")
    samples={k:[] for k in keys}
    plot_fields=["BT","BI","ST","SI"]
    by_plot={f:[] for f in plot_fields}
    for i in range(draws):
        res=means_of_pairs(calc_paired(pairs,exponent,rng,permute=True))
        for k in keys:samples[k].append(res[k])
        for f in plot_fields:by_plot[f].append(res["plot_alignments"][f])
    result={"weight_exponent":exponent,"n_draws":draws,"seed":seed,
            "observed":observed,"technical_null":{}}
    def summarize(obs,v):
        v=sorted(v)
        n=len(v)
        greater=sum(abs(x)>=abs(obs)-1e-12 for x in v)
        return {"observed":obs,"null_mean":mean(v),
            "null_95pct":[v[int(.025*n)],v[int(.975*n)]],
            "two_sided_fraction_plus_one":(greater+1)/(n+1)}
    for k,v in samples.items():
        result["technical_null"][k]=summarize(observed[k],v)
    result["plot_abundance_alignment_null"]={
         f:summarize(observed["plot_alignments"][f],by_plot[f])
         for f in plot_fields}
    # Leave an entire matched between-station/station spatial pair out.
    if len(pairs)>=3:
        jack=[means_of_pairs(calc_paired(pairs[:i]+pairs[i+1:],exponent))
              for i in range(len(pairs))]
        result["leave_one_true_spatial_pair_out"]={
            k:{"min":min(j[k] for j in jack),"max":max(j[k] for j in jack),
               "all_positive":all(j[k]>0 for j in jack),
               "all_negative":all(j[k]<0 for j in jack)}
            for k in keys}
    return result

def analyze(original,eive_xlsx):
    source,sourceqa=parse_xlsx_matrix(Path(original))
    names=species_names_from_publisher(original)
    eive,eiveqa=eive_table(eive_xlsx)
    meta,metaqa=site_metadata(original)
    abundance,abundqa=read_original_plant_values(original,source)
    out={"version":"v337",
        "source_original_kashmir_sha256":sourceqa["original_file_sha256"],
        "source_original_kashmir_paper":PAPER_DOI,
        "independent_temperature_indicator_EIVE":EIVE_DOI,
        "EIVE_source_original_md5":eiveqa["EIVE_original_md5"],
        "source_observed_abundance_numeric_QA":abundqa,
        "original_paper_measurement":"At each site, twelve 1x1m quadrats for species abundance; units of each Dryad cell not further defined in manuscript sampling method",
        "field_site_source_geoQA":metaqa,
        "panels":[]}
    for index,(year,traitgate,matchgate) in enumerate(PANELS):
        matched,srcqa,matchqa=cohort(source,names,eive,meta,year,traitgate,matchgate)
        item={"year":year,"source_trait_gate":traitgate,
              "source_geographical_match_gate":matchgate,
              "source_mapped_site_eligibility":srcqa,
              "source_matching_QA":matchqa}
        if not matched:
            item["status"]="no genuine same-source geographic pairs"
            out["panels"].append(item);continue
        pairs=matched_source_plots(
            abundance,matched,names,eive,year,TRAIT_GATES[traitgate][0])
        for alpha in (1.0,0.5):
            result=conditional_null(pairs,alpha,PERM_DRAWS,SEED+100*index+(0 if alpha==1 else 17))
            if alpha==1:
                previous=group_comparison(abundance,matched,names,eive,year,traitgate,
                         "original_positive_source_value",METRIC)
            else:
                previous=group_comparison(abundance,matched,names,eive,year,traitgate,
                         "sqrt_original_source_value",METRIC)
            original_inc=group_comparison(abundance,matched,names,eive,year,traitgate,"presence",METRIC)
            obs=result["observed"]
            if not math.isclose(obs["relative_presence"],
                                original_inc["mean_between_minus_station_track_relative_interior"],
                                abs_tol=1e-10,rel_tol=0):
                raise AssertionError("v334 source presence CWM baseline differs in v337")
            if not math.isclose(obs["relative_abundance"],
                                previous["mean_between_minus_station_track_relative_interior"],
                                abs_tol=1e-10,rel_tol=0):
                raise AssertionError("v334 source abundance CWM baseline differs in v337")
            item.setdefault("exponent_sensitivity",[]).append(result)
        out["panels"].append(item)
    out["interpretation_gates"]=[
       "The genuine publisher paper explicitly calls source quantity species abundance, but does not specify whether Excel values are raw quadrat counts, sums, means or transformations",
       "Species-presence mean and abundance-weighted CWM difference is exactly the covariance of numeric source abundance with species EIVE, divided by mean source abundance",
       "Within-plot species abundance shuffle conditions on source taxa, observed original species richness, original abundances, original 1:1 geographic matches and year; it does not randomize railway or station treatment",
       "Within-plot abundance shuffling does not preserve intrinsic species-specific abundance preferences across plots and is a liberal null conditional on each plot; cannot alone prove temperature filter causes dominance",
       "Source EIVE indicator is European-trained taxon temperature preference, not Himalayan physiological tolerance, measured air temperature or actual ballast microclimate",
       "Only Kashmir one observation programme, 31 spatially dependent original field sites across 2014 and 2017",
       "No actual individual-species seed movement, railway traffic frequencies, train winds, soil temperature, successional mechanisms or seed fitness directly measured"]
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze(a.kashmir_original,a.eive_original)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    compact=[{"year":row["year"],"source_gate":row["source_trait_gate"],
              "match":row["source_geographical_match_gate"],
              "n":row["source_matching_QA"]["matched_count"],
              **({"exponent_sensitivity":row["exponent_sensitivity"]}
                 if "exponent_sensitivity" in row else {})} for row in out["panels"]]
    print(json.dumps({"version":"v337","results":compact},
                     ensure_ascii=False),flush=True)
if __name__=="__main__":main()

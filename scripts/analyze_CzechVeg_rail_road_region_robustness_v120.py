#!/usr/bin/env python3
"""v12.0 genuine CzechVeg–Pladias original-source phyto-region diagnostic.

Uses only the v11.9 verified source JSON aggregate, with the exact
public CzechVeg v1.1 117,739 plot/header + 2,510,340 plant rows and
published Pladias diaspore/lifeform original-source validated.
No raw source redistribution, no invented within-site railroad zones.
For each source phyto-region with >=5 text-labelled railway and road
candidate plots AND evaluable 3+ wind & nonwind taxa per lifeform,
compare wind-vs-other rail-minus-road prevalence under equal-species
and species source-plot opportunity weighting.

Report region-weighted and source plot-opportunity-weighted means,
region bootstrap and leave-one-region-out, plus stricter nrail,nroad
thresholds, and same-region strict-vs-broad site labeling changes.
These are observational sensitivity checks, NOT valid causal
RAIL-vs-ROAD or genuine Slovenia BS−BT within-railway replication.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import Counter
from pathlib import Path
from statistics import mean,median

SOURCE_V119="v11.9"
SOURCE_ZENODO_HEADER_MD5="0bc920462d1f9eac473e7b2a231e6782"
SOURCE_ZENODO_SPECIES_MD5="b296e5e19f03fb3bad0a8347c317f8a3"
BOOT_N=8000
RNG_SEED=20261114


def validated_source(report):
    if report.get("complete") is not True or report.get("version")!=SOURCE_V119:
        raise ValueError("Authenticated original Czech source v11.9 result not present")
    a=report.get("original_source_Czechveg_header_provenance") or {}
    b=report.get("original_source_Czechveg_species_provenance") or {}
    if a.get("checksum")!=SOURCE_ZENODO_HEADER_MD5 or b.get("checksum")!=SOURCE_ZENODO_SPECIES_MD5:
        raise ValueError("Original source Zenodo Czech vegetation record MD5 mismatch")
    tr=report.get("original_source_Pladias_trait_quality") or {}
    if (tr.get("source_Pladias_original_diaspore_rows")!=3147 or
        tr.get("source_Pladias_original_lifeform_rows")!=4599 or
        not tr.get("dynamic_publisher_XLSX_semantic_original_source_validation",{}).
            get("dispersal",{}).get("all_five_v118_actual_taxon_and_strategy_schema_anchors_equal",False)):
        raise ValueError("Original Pladias v11.8 verified strategy source not available")
    return report


def pct(xs,p):
    if not xs:return None
    s=sorted(xs);idx=(len(s)-1)*p
    lo=int(idx);hi=min(lo+1,len(s)-1)
    return s[lo]+(s[hi]-s[lo])*(idx-lo)


def summary(rows,draws=BOOT_N,seed=RNG_SEED):
    if not rows:return {"identified":False,"n_original_source_eligible_phyto_regions":0}
    ids=[r["source_phytogeographic_region"] for r in rows]
    if len(set(ids))!=len(ids):
        raise ValueError("Original CzechVeg region code repeated in same sensitivity cohort")
    a=[float(r["source_equal"]) for r in rows]
    b=[float(r["source_occupancy"]) for r in rows]
    if any(not math.isfinite(v) for v in a+b):
        raise ValueError("Missing nonfinite original biological source contrast")
    w=[int(r["nrail"])+int(r["nroad"]) for r in rows]
    if min(w)<=0:raise ValueError("At least one original region has no sampled source sites")
    mean_eq=mean(a);mean_oc=mean(b)
    weighted_eq=sum(z*wt for z,wt in zip(a,w))/sum(w)
    weighted_oc=sum(z*wt for z,wt in zip(b,w))/sum(w)
    rng=random.Random(seed)
    draws_eq=[];draws_oc=[];gap=[]
    for _ in range(draws):
        ind=[rng.randrange(len(rows)) for _ in rows]
        e=mean(a[i] for i in ind)
        o=mean(b[i] for i in ind)
        draws_eq.append(e);draws_oc.append(o);gap.append(e-o)
    loo=[]
    if len(rows)>1:
        for i,r in enumerate(rows):
            q=mean(a[j] for j in range(len(a)) if j!=i)
            u=mean(b[j] for j in range(len(b)) if j!=i)
            loo.append({
                "omitted_original_source_phyto_region":r["source_phytogeographic_region"],
                "mean_species_equal":q,"mean_occ_weighted":u,
                "mean_signs_opposite":q*u<0,
            })
    return {
        "identified":True,
        "n_original_source_eligible_phyto_regions":len(rows),
        "n_original_total_original_source_candidate_rail_sites":sum(int(r["nrail"]) for r in rows),
        "n_original_total_original_source_candidate_road_sites":sum(int(r["nroad"]) for r in rows),
        "n_original_regions_positive_equal":sum(v>0 for v in a),
        "n_original_regions_positive_occupancy_weight":sum(v>0 for v in b),
        "n_original_regions_positive_equal_negative_occ":sum(x>0 and y<0 for x,y in zip(a,b)),
        "n_original_regions_negative_equal_positive_occ":sum(x<0 and y>0 for x,y in zip(a,b)),
        "n_original_regions_any_opposite_sign":sum(x*y<0 for x,y in zip(a,b)),
        "original_region_equal_mean_species_equal":mean_eq,
        "original_region_equal_mean_occupancy_weighted":mean_oc,
        "original_mean_gap_equal_minus_occupancy":mean_eq-mean_oc,
        "original_original_plots_weight_region_mean_species_equal":weighted_eq,
        "original_original_plots_weight_region_mean_occupancy":weighted_oc,
        "original_region_equal_means_sign_opposite":mean_eq*mean_oc<0,
        "region_resampling_95_species_equal":[pct(draws_eq,.025),pct(draws_eq,.975)],
        "region_resampling_95_occ_weighted":[pct(draws_oc,.025),pct(draws_oc,.975)],
        "region_resampling_95_gap":[pct(gap,.025),pct(gap,.975)],
        "n_region_resample_draws":draws,
        "n_leave_one_source_region_out":len(loo),
        "n_leave_one_source_region_out_retain_mean_opposite_signs":
            sum(z["mean_signs_opposite"] for z in loo),
        "leave_one_region_source_values":loo,
        "conditional_on_nonrandom_source_region_codes_no_random_country_CI":True,
    }


def conditional_regions(panel,min_rail,min_road):
    cand=panel["source_region_comparisons_not_pairwise_causal_rail_vs_road"]
    return [
        x for x in cand if x.get("identified") is True
        and int(x["nrail"])>=min_rail
        and int(x["nroad"])>=min_road
        and x.get("source_equal") is not None
        and x.get("source_occupancy") is not None
    ]


def run(source):
    a=validated_source(source)
    strictr=a["region_stratified_high_precision_source_comparison"]
    broader=a["region_stratified_broader_source_comparison"]
    strict=conditional_regions(strictr,5,5)
    broad=conditional_regions(broader,5,5)
    if len(strict)!=strictr["n_original_regions_with_guild_comparison_identified"]:
        raise ValueError("Strict region source count does not reconcile")
    if len(broad)!=broader["n_original_regions_with_guild_comparison_identified"]:
        raise ValueError("Broad region source count does not reconcile")
    all_overlaps=[]
    m1={z["source_phytogeographic_region"]:z for z in strict}
    m2={z["source_phytogeographic_region"]:z for z in broad}
    for code in sorted(set(m1)&set(m2)):
        x=m1[code];y=m2[code]
        all_overlaps.append({
            "source_phyto_region":code,"strict_nrail":x["nrail"],
            "broad_nrail":y["nrail"],
            "strict_nroad":x["nroad"],"broad_nroad":y["nroad"],
            "strict_equal":x["source_equal"],"broad_equal":y["source_equal"],
            "strict_occ":x["source_occupancy"],"broad_occ":y["source_occupancy"],
            "species_equal_sign_agrees":x["source_equal"]*y["source_equal"]>0,
            "occ_weight_sign_agrees":x["source_occupancy"]*y["source_occupancy"]>0,
        })
    return {
        "complete":True,"version":"v12.0",
        "original_source_CzechVeg_DOI":"10.5281/zenodo.19231152",
        "original_source_Pl adias_DOI".replace(" ",""):"10.23855/preslia.2021.001",
        "source_derived_v119_verified":True,
        "n_original_veg_plot_rows_verified":117739,
        "n_original_species_plot_rows_scanned":source[
            "original_source_species_plot_and_trait_join_quality"][
            "n_original_source_species_csv_records_scanned"],
        "original_high_precision_HABITAT_source_railroad_vs_road":{
            "n_original_source_categorical_rail":a["high_precision_original_HABITAT_text_label_RAIL_vs_ROAD"][
                "n_original_rail_candidate_plots"],
            "n_original_source_categorical_road":a["high_precision_original_HABITAT_text_label_RAIL_vs_ROAD"][
                "n_original_road_control_plots"],
            "all_eligible_phyto_regions":summary(strict),
            "high_source_replication_min10_rail10_road":summary(
                conditional_regions(strictr,10,10)),
        },
        "original_broad_HABITAT_and_LOCALITY_source_railroad_vs_road":{
            "n_original_source_categorical_rail":a["broader_source_site_HABITAT_or_LOCALITY_label_RAIL_vs_ROAD"][
                "n_original_rail_candidate_plots"],
            "n_original_source_categorical_road":a["broader_source_site_HABITAT_or_LOCALITY_label_RAIL_vs_ROAD"][
                "n_original_road_control_plots"],
            "all_eligible_phyto_regions":summary(broad),
            "high_source_replication_min10_rail10_road":summary(
                conditional_regions(broader,10,10)),
            "high_source_replication_min20_rail20_road":summary(
                conditional_regions(broader,20,20)),
        },
        "overlap_strict_broad_exactly_same_phyto_source_regions":{
            "n_overlap_regions":len(all_overlaps),
            "n_sign_agree_equal":sum(z["species_equal_sign_agrees"] for z in all_overlaps),
            "n_sign_agree_occ":sum(z["occ_weight_sign_agrees"] for z in all_overlaps),
            "regions":all_overlaps,
        },
        "scientific_interpretation_limitations":[
            "Czech original study comparison is railway-text vs road-text across separate source vegetation plots, never equal-area same-site BS/BT source contrast.",
            "Year/area/habitat and effort not harmonized within phyto regions, and strict vs broad labels may include off-track vegetation.",
            "Species dispersal archetype and source Raunkiaer lifeform proxies derive from publicly exported Pladias and are not observed actual propagule travel or exact Slovenian annual lifespan coding.",
            "Observed candidate site frequencies are neither individual counts nor biomass or train wind/mowing measurements.",
            "Resampling 28 source regional codes does not estimate independent random countries; ecological causal inference is NOT established.",
            "No data from original author's raw CzechVeg or Pladias files redistributed.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--verified119",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=json.loads(Path(a.verified119).read_text(encoding="utf8"))
    result=run(data)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    sections=result["original_broad_HABITAT_and_LOCALITY_source_railroad_vs_road"]
    strict=result["original_high_precision_HABITAT_source_railroad_vs_road"]
    print(json.dumps({
        "complete":True,"version":"v12.0",
        "source_veg_n":result["n_original_veg_plot_rows_verified"],
        "strict":{k:v for k,v in strict.items()},
        "broad":{k:v for k,v in sections.items()},
        "same_regions_strict_broad":result["overlap_strict_broad_exactly_same_phyto_source_regions"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__": main()

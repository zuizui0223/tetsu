#!/usr/bin/env python3
"""v343: within ONE original Czech vegetation publication railway substrata.

After full v341 pooled Czech source analysis identified apparent opposite rail
yard vs embankment thermal differences, v342 demonstrated only one original
REFERENCE has both habitat classes. This source-level descriptive check uses
the same ORIGINAL G2 outcome-blind geographic plot-control matches; only the
study possessing both environments is selected BEFORE examining plant traits.

True CzechVeg source headers/species and independent EIVE verified by MD5.
No source raw plot locations, identities or individual source species rows
are redistributed. No railway-causal habitat contrast can be estimated from
nonrandom plots in just one original vegetation study.
"""
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

import numpy as np
from scipy.optimize import linear_sum_assignment

from audit_CzechVeg_actual_railways_nearby_controls_v340 import distance_km
from analyze_CzechVeg_rail_nearby_EIVE_source_v341 import (
    load_headers,match_source_outcome_blind,load_selected_species,
    compute_community)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table

GATE="G2_same_source_extended"
MAIN=(3,5,.5)
RELAXED=(1,3,.3)
HABITATS=("railway_track_yard_or_structure","railway_embankment_ditch_or_cut")

def original_study_with_both(pairs):
    by=defaultdict(Counter)
    for p in pairs:
        by[p["rail"]["reference"]][p["rail_habitat"]]+=1
    dual=[study for study,counts in by.items() if all(counts[h]>0 for h in HABITATS)]
    if len(dual)!=1:raise ValueError("v342 unique both-habitat study support changed")
    ref=dual[0]
    sub=[p for p in pairs if p["rail"]["reference"]==ref]
    if len(sub)!=27 or Counter(p["rail_habitat"] for p in sub)!=Counter({
         HABITATS[0]:18,HABITATS[1]:9}):
        raise ValueError("v342 source strata drift")
    return sub,hashlib.sha256(ref.encode()).hexdigest()[:10]

def observed_effects(pairs,plants,eive,metric,systems,min_species,min_cov):
    results=[]
    for p in pairs:
        a=compute_community(plants.get(p["rail"]["id"]),eive,metric,systems)
        b=compute_community(plants.get(p["reference"]["id"]),eive,metric,systems)
        if any(z["trait_coverage"] is None for z in (a,b)):continue
        if min(a["trait_n"],b["trait_n"])<min_species or min(
              a["trait_coverage"],b["trait_coverage"])<min_cov:
            continue
        effect=a["equal_species_cwm"]-b["equal_species_cwm"]
        cover=None
        if a["cover_weighted_cwm"] is not None and b["cover_weighted_cwm"] is not None:
            if min(a["trait_source_cover_nonmissing_n"],b["trait_source_cover_nonmissing_n"])>=min_species and min(
               a["source_cover_known_trait_fraction"],b["source_cover_known_trait_fraction"])>=.8:
                cover=a["cover_weighted_cwm"]-b["cover_weighted_cwm"]
        results.append({"rail":p["rail"],"habitat":p["rail_habitat"],
                        "effect":effect,"cover_effect":cover,
                        "dist_to_reference_km":p["distance_km"],
                        "eive_frac_rail":a["trait_coverage"],
                        "eive_frac_reference":b["trait_coverage"],
                        "source_rail_species":a["source_species_n"],
                        "source_ref_species":b["source_species_n"]})
    return results

def type_diffs(items):
    by=defaultdict(list)
    for v in items:by[v["habitat"]].append(v)
    out={}
    for category in HABITATS:
        z=by[category]
        values=[x["effect"] for x in z]
        cv=[x["cover_effect"] for x in z if x["cover_effect"] is not None]
        out[category]={
          "n_site_vs_nearby_same_study_pairs":len(z),
          "EIVE_centre_or_width_rail_minus_nearby_mean":mean(values) if values else None,
          "positive_site_pairs":sum(x>0 for x in values),
          "negative_site_pairs":sum(x<0 for x in values),
          "source_cover_weighted_EIVE_difference_mean":mean(cv) if cv else None,
          "n_usable_source_cover_species_pairs":len(cv),
          "mean_source_species_richness_rail":mean([x["source_rail_species"] for x in z]) if z else None,
          "mean_source_species_richness_reference":mean([x["source_ref_species"] for x in z]) if z else None}
    lhs=out[HABITATS[0]]["EIVE_centre_or_width_rail_minus_nearby_mean"]
    rhs=out[HABITATS[1]]["EIVE_centre_or_width_rail_minus_nearby_mean"]
    return out,lhs-rhs if lhs is not None and rhs is not None else None

def nearby_rail_type_pairs(items,max_dist_km=2,max_alt_m=50):
    yards=sorted([x for x in items if x["habitat"]==HABITATS[0]],
        key=lambda x:x["rail"]["id"])
    verges=sorted([x for x in items if x["habitat"]==HABITATS[1]],
        key=lambda x:x["rail"]["id"])
    if not yards or not verges:return {"n_source_rail_types_nearby":0}
    n=len(yards);m=len(verges)
    mat=np.full((n,m+n),1e8)
    mat[:,m:]=1e6
    physical={}
    for i,a in enumerate(yards):
        for j,b in enumerate(verges):
            x=a["rail"];y=b["rail"]
            if any(z[k] is None for z in (x,y) for k in ("lat","lon","altitude")):
                continue
            dist=distance_km(x,y)
            elev=abs(x["altitude"]-y["altitude"])
            if dist<=max_dist_km and elev<=max_alt_m:
                mat[i,j]=dist+0.02*elev+1e-8*(i*m+j)
                physical[i,j]=(dist,elev)
    row,col=linear_sum_assignment(mat)
    vals=[]
    for i,j in zip(row,col):
        if j>=m or mat[i,j]>=1e6:continue
        a,b=yards[i],verges[j]
        dist,alt=physical[i,j]
        vals.append({"rail_yard_vs_embankment_difference_of_rail_nearby_effect":
                       a["effect"]-b["effect"],"distance":dist,"alt":alt})
    return {"n_source_rail_types_nearby":len(vals),
        "source_rail_geodesic_max_km":max_dist_km,
        "source_altitude_max_m":max_alt_m,
        "mean_rail_yard_minus_embankment_thermal_effect_for_nearby_source_sites":
            mean([v["rail_yard_vs_embankment_difference_of_rail_nearby_effect"]
                  for v in vals]) if vals else None,
        "positive_spatially_matched_difference_count":sum(
           v["rail_yard_vs_embankment_difference_of_rail_nearby_effect"]>0 for v in vals),
        "median_source_track_to_embankment_geodesic_km":median(
             v["distance"] for v in vals) if vals else None,
        "median_source_track_to_embankment_altitude_m":median(
             v["alt"] for v in vals) if vals else None,
        "spatial_matches_not_independent_railway_studies_or_random_experiment":True}

def analyze(header,species,eive):
    rows,headqa=load_headers(header)
    all_pairs,matchqa=match_source_outcome_blind(rows,GATE)
    if len(all_pairs)!=39:raise ValueError("Original v341 39 matched pairs drift")
    source_pairs,study_hash=original_study_with_both(all_pairs)
    plot_ids={r[k]["id"] for r in source_pairs for k in ("rail","reference")}
    flora,florqa=load_selected_species(species,plot_ids)
    traits,traitqa=eive_table(eive)
    result={"version":"v343","original_source_study_hash_not_id":study_hash,
            "original_G2_matched_pairs":len(all_pairs),
            "original_study_matched_pairs_both_habitat_types":len(source_pairs),
            "number_original_independent_studies_in_this_subanalysis":1,
            "study_initial_rail_habitat_counts":dict(Counter(p["rail_habitat"] for p in source_pairs)),
            "source_species_rows_qa":florqa,
            "independent_EIVE_md5":traitqa["EIVE_original_md5"],
            "panels":[]}
    for metric in ("position","width"):
        for systems,n,coverage in (MAIN,RELAXED):
            selected=observed_effects(source_pairs,flora,traits,metric,systems,n,coverage)
            types,delta=type_diffs(selected)
            result["panels"].append({"metric":metric,"minimum_EIVE_systems":systems,
                "min_numeric_species_per_zone":n,"minimum_EIVE_trait_coverage":coverage,
                "n_source_plot_pair_EIVE_eligible":len(selected),
                "within_one_original_publication_source_rail_habitat_types":types,
                "unmatched_within_publication_yard_minus_embankment_difference":delta,
                "geographically_near_rail_type_pair_sensitivity":
                    nearby_rail_type_pairs(selected)})
    result["limits"]=[
      "One original Czech vegetation source only, and selection of this study occurred after opposite pooled habitat effects were observed",
      "Control plots are nearby vegetation with no rail HABITAT/LOCALITY term, NOT proven unaffected natural controls",
      "Original rail structure includes different wall, yard and track substrates, not universally active ballast",
      "Within-source rail-yard and embankment plots were NOT randomly located, and unrecorded soil/shading/management remain confounded",
      "Geographical distance <=2 km and altitude <=50m near-type pairing does not ensure same vegetation community source pool or causal microclimate",
      "Original EIVE traits are plant temperature preference indices, not physiological heat tolerance or temperature degrees C",
      "No study-independent railway climatic range-limit release, train dispersal, plant reproduction or evolution measured",
      "No original plot GPS, species observations or original plot IDs republished"]
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-headers",required=True)
    p.add_argument("--source-species",required=True)
    p.add_argument("--eive-xlsx",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=analyze(a.source_headers,a.source_species,a.eive_xlsx)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v343","study_source_hash":obj["original_source_study_hash_not_id"],
       "original_study_matched_pairs_both_habitat_types":obj["original_study_matched_pairs_both_habitat_types"],
       "panels":obj["panels"]},ensure_ascii=False),flush=True)
if __name__=="__main__":main()

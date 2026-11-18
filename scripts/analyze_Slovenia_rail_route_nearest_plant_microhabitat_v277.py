#!/usr/bin/env python3
"""v27.7 — genuine Slovenian RAIL plants: same-location vs nearest ORIGINAL
railway location on the SAME national line, using author's TRUE route chainage.

Reanalysis of Šabić et al. (2026), DOI 10.1007/s44473-026-00155-z.
Original 401279-byte publisher XLSX, SHA256 immutable verified upstream.
Genuine 89 matched BS shoulder/BT toe plant sites, 173 original left/right
side pairs, 15 actual national line sections; union left/right *within*
location before any comparison, NEVER 173 independent sites.

For each source site i on rail section L:
  S0_i = Jaccard similarity(BS_i, BT_i)           # same original location
  S1_i = Jaccard similarity(BS_i, BT_j)           # nearest OTHER location
         j chosen only by original midpoint rail-route kilometre, NEVER
         by botanical outcome.
  delta_i = S0_i - S1_i (positive: local floristic association).
Distances are absolute author original source national-line chainages
(start_km+end_km)/2, NOT a global network shortest path or train route
connection, not a train service. Cross-check geodesic site GPS to diagnose
line-km mismatches, not to turn straight-line distance into service flow.

Pre-fixed nearby source sensitivity radii 5/10/20/50/ALL km. If no
real other location within a radius, exclude (not zero-impute). Additional
same-year, same-SINGLE/DOUBLE-layout, both matched candidate panels.
Each source site once => equal site within original national line => equal
original national line. Require >=2 sites for line to estimate and >=4
independent original sections AND >=12 actual sites for a population-level
interpretation of any radius. LOO line-exclusion sensitivity reported.
Zero-distance same-source paired field habitat dissimilarity is not
causal gene, seed origin, local reproductive isolation, railway-mediated
dispersal, environmental filter, mowing/wind adaptation, or matched
non-rail control. Both source habitats lie in real railway sections.

A second descriptive plot shows observed cross-zone species similarity in
distance bins among genuine different-location directed BS_i vs BT_j
comparisons within original line. Dyads are NOT independent field samples;
per-line means and bin support reported without a fabricated p value.

Original plant-by-site lists, GPS, km, exact site names and individual
pairwise records are not exported; aggregated section-class/source count
and biological dissimilarity metrics only.
"""
from __future__ import annotations
import argparse, json, math
from collections import defaultdict,Counter
from datetime import datetime,timedelta
from pathlib import Path
from statistics import mean,median
from itertools import permutations

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx, validate_original_keys
from analyze_rail_local_site_coupling_v257 import slovenia_location_sets
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from analyze_Slovenia_gov_official_main_regional_traits_v225 import official_class

RADII=(5.,10.,20.,50.,math.inf)
MODES=("all","same_year","same_layout","same_year_and_layout")
DISTANCE_BINS=((0.,5.),(5.,20.),(20.,50.),(50.,math.inf))
MIN_LINES_REPORT=4
MIN_SITES_REPORT=12
EXPECTED={"source_actual_sites":89,"source_sections":15,
          "same_section_pairs_undirected":340,
          "source_eligible_multisite_sections":14,
          "source_eligible_multisite_locations":88}


def jaccard_similarity(a:set,b:set)->float:
    if not a or not b:raise ValueError("Both authentic original plant plots must be nonempty")
    return len(a&b)/len(a|b)


def to_original_sites(original:dict)->tuple[dict,dict]:
    source,locqa=slovenia_location_sets(original)
    loc,zones,taxa,traits,names,audit=validate_original_keys(original)
    result=defaultdict(list)
    source_bad=Counter()
    for lid,v in source.items():
        g=loc[lid]
        section=str(v["section"])
        typ=str(g["track_type"])
        if typ not in ("SINGLE","DOUBLE") or official_class(section) is None:
            raise ValueError("Source real track layout or original national line missing")
        try:
            left,right=float(g["start_km"]),float(g["end_km"])
            lat0,lon0,lat1,lon1=(float(g[k]) for k in
                                   ("start_lat","start_lon","end_lat","end_lon"))
            stamp=float(g["date"])
        except (KeyError,TypeError,ValueError) as err:
            raise ValueError("Missing exact original publisher route km/GPS/survey date") from err
        if not all(map(math.isfinite,(left,right,lat0,lon0,lat1,lon1,stamp))):
            raise ValueError("Original route km/GPS not finite")
        if min(left,right)<0 or not (44<=lat0<=48 and 44<=lat1<=48 and
                                     12<=lon0<=19 and 12<=lon1<=19):
            raise ValueError("Original route source km/GPS malformed")
        year=(datetime(1899,12,30)+timedelta(days=stamp)).year
        if year not in (2022,2023):raise ValueError("Original field collection year unsupported")
        if haversine_km(lat0,lon0,lat1,lon1)>5:
            raise ValueError("Author source same 100m transect endpoints too far apart")
        result[section].append({
            "source_id":str(lid),"section":section,"class":official_class(section),
            "year":year,"layout":typ,
            "routekm":(left+right)/2,
            "lat":(lat0+lat1)/2,"lon":(lon0+lon1)/2,
            "BS":v["BS"],"BT":v["BT"]
        })
    for sec in result:result[sec].sort(key=lambda x:(x["routekm"],x["source_id"]))
    sum_sites=sum(map(len,result.values()))
    sum_pairs=sum(len(x)*(len(x)-1)//2 for x in result.values())
    cohort={
        "source_actual_sites":sum_sites,
        "source_sections":len(result),
        "same_section_pairs_undirected":sum_pairs,
        "source_eligible_multisite_sections":sum(len(x)>=2 for x in result.values()),
        "source_eligible_multisite_locations":sum(len(x) for x in result.values() if len(x)>=2),
    }
    if cohort!=EXPECTED:raise ValueError(f"Strict source 89 sites 15 lines and 340 true within-line original location pairs drift: {cohort}")
    return dict(result),{
        "source_cohort":cohort,
        "source_publisher_QA":locqa,
        "author_original_2022_2023_source_plots":True,
        "exact_original_route_chainage_source_in_km":True,
        "geographic_GPS_not_train_service_path":True,
        "source_by_original_G_R":{"G":sum(len(x) for x in result.values() if x[0]["class"]=="main_G"),
                                 "R":sum(len(x) for x in result.values() if x[0]["class"]=="regional_R")},
    }


def source_neighbors(items:list[dict], radius:float, mode:str)->list[dict]:
    if mode not in MODES:raise ValueError("Unknown source matching constraint")
    if radius<=0:raise ValueError("Source radius must be positive")
    pairs=[]
    for this in items:
        eligible=[]
        for other in items:
            if other["source_id"]==this["source_id"]:continue
            distance=abs(this["routekm"]-other["routekm"])
            if distance>radius:continue
            if mode in ("same_year","same_year_and_layout") and this["year"]!=other["year"]:
                continue
            if mode in ("same_layout","same_year_and_layout") and this["layout"]!=other["layout"]:
                continue
            geo=haversine_km(this["lat"],this["lon"],other["lat"],other["lon"])
            eligible.append((distance,geo,other["source_id"],other))
        if not eligible:continue
        # Species-independent nearest-neighbor matching, no tuning by outcome.
        d,geo,idv,other=min(eligible,key=lambda x:x[:3])
        same=jaccard_similarity(this["BS"],this["BT"])
        neighbor=jaccard_similarity(this["BS"],other["BT"])
        pairs.append({
            "source_id":this["source_id"],
            "class":this["class"],"year":this["year"],
            "same_similarity":same,
            "neighbor_similarity":neighbor,
            "same_minus_nearest_other":same-neighbor,
            "route_gap_km":d,"geo_gap_km":geo,
            "route_geodesic_gap_km":d-geo
        })
    return pairs


def bin_for_distance(gap:float)->str:
    if not math.isfinite(gap) or gap<0:raise ValueError("Original source separation invalid")
    for lo,hi in DISTANCE_BINS:
        if lo<=gap<hi:return f"{lo:g}-{hi:g}" if math.isfinite(hi) else f"{lo:g}+"
    raise ValueError("Original source distance outside allowed bins")


def section_distance_profile(items:list[dict])->dict:
    bybin=defaultdict(list)
    geodiff=[]
    for i,j in permutations(range(len(items)),2):
        a,b=items[i],items[j]
        gap=abs(a["routekm"]-b["routekm"])
        geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
        # Different locations can have source coordinates nearly identical;
        # this is reported as source spatial precision limitation.
        geodiff.append(gap-geo)
        bybin[bin_for_distance(gap)].append(jaccard_similarity(a["BS"],b["BT"]))
    return {
        "n_different_original_location_directed_source_pairs":sum(len(v) for v in bybin.values()),
        "n_same_section_original_sites":len(items),
        "geodesic_distance_exceeds_original_route_km_count":sum(v < -0.5 for v in geodiff),
        "route_minus_geodesic_distance_km_median":median(geodiff) if geodiff else None,
        "bins":{k:{"n_dependent_directed_pairs_NOT_independent":len(v),
                   "mean_species_Jaccard_similarity":mean(v)}
                for k,v in sorted(bybin.items())}
    }


def summarise(items_by_section:dict,mode:str,radius:float)->dict:
    per_section={}
    qualified_pairs=0
    for section,items in sorted(items_by_section.items()):
        if len(items)<2:continue
        records=source_neighbors(items,radius,mode)
        if len(records)<2:continue # source >=2 real field sites per original section
        per_section[section]={
            "source_actual_sites_in_section":len(items),
            "source_nearest_other_eligible_site_count":len(records),
            "source_original_line_class":items[0]["class"],
            "mean_true_same_site_similarity":mean(x["same_similarity"] for x in records),
            "mean_other_nearest_site_similarity":mean(x["neighbor_similarity"] for x in records),
            "mean_same_minus_nearest":mean(x["same_minus_nearest_other"] for x in records),
            "median_nearest_route_separation_km":median(x["route_gap_km"] for x in records),
            "n_original_chains_with_geodesic_exceeding_source_km_by_halfkm":
                sum(x["route_geodesic_gap_km"]<-.5 for x in records),
        }
        qualified_pairs+=len(records)
    n=len(per_section)
    if not n:
        return {"identified":False,"source_sections":0,"source_sites":0}
    delta=[r["mean_same_minus_nearest"] for r in per_section.values()]
    loorange=[mean(delta[:i]+delta[i+1:]) for i in range(n)] if n>1 else []
    return {
        "identified": n>=MIN_LINES_REPORT and qualified_pairs>=MIN_SITES_REPORT,
        "source_real_rail_sections":n,
        "source_eligible_original_spatial_sites":qualified_pairs,
        "source_original_G_main_sections":sum(r["source_original_line_class"]=="main_G" for r in per_section.values()),
        "source_original_R_regional_sections":sum(r["source_original_line_class"]=="regional_R" for r in per_section.values()),
        "original_section_equal_same_site_species_similarity":mean(x["mean_true_same_site_similarity"] for x in per_section.values()),
        "original_section_equal_nearest_other_site_species_similarity":mean(x["mean_other_nearest_site_similarity"] for x in per_section.values()),
        "original_section_equal_same_site_minus_nearest_species_similarity":mean(delta),
        "source_original_line_sections_with_positive_same_nearest_gap":sum(v>0 for v in delta),
        "source_original_line_sections_with_negative_same_nearest_gap":sum(v<0 for v in delta),
        "leave_one_original_rail_section_out_min_max":[min(loorange),max(loorange)] if loorange else None,
        "median_of_section_median_nearest_route_km":median(x["median_nearest_route_separation_km"] for x in per_section.values()),
        "source_route_km_minus_geographic_distance_large_mismatch_record_count":
            sum(x["n_original_chains_with_geodesic_exceeding_source_km_by_halfkm"] for x in per_section.values()),
        "sampling_design":"original site -> national line equal; true same-site BS/BT vs nearest other original BT on SAME official line; no source species list imputed",
        "source_section_results_aggregate_only":per_section,
        "no_randomized_field_class_effect":True
    }


def distance_bins(items_by_section:dict)->dict:
    # Each source national section gets one vote per distance bin, no
    # duplicated site-to-site dyad counts as independent field experiments.
    by_bin=defaultdict(list)
    all_sections={}
    for section,items in sorted(items_by_section.items()):
        if len(items)<2:continue
        result=section_distance_profile(items)
        all_sections[section]={
            "source_sites":len(items),
            "source_within_rail_line_directed_pair_count":
                result["n_different_original_location_directed_source_pairs"],
            "source_route_geodesic_mismatch_count":
                result["geodesic_distance_exceeds_original_route_km_count"]
        }
        for bin_id,q in result["bins"].items():
            by_bin[bin_id].append(q["mean_species_Jaccard_similarity"])
    bins={}
    for bin_id,values in sorted(by_bin.items()):
        bins[bin_id]={
            "source_original_sections_with_different_site_opportunities":len(values),
            "original_section_equal_mean_cross_habitat_species_similarity":
                mean(values),
            "not_independent_source_site_pairs":True,
        }
    return {
        "source_multisite_rail_sections":len(all_sections),
        "source_different_location_directed_pairs_sum_not_independent":
            sum(x["source_within_rail_line_directed_pair_count"] for x in all_sections.values()),
        "source_pairs_route_distance_less_than_geodesic_halfkm_mismatch":
            sum(x["source_route_geodesic_mismatch_count"] for x in all_sections.values()),
        "cross_site_distance_bin_equal_section":bins
    }


def run(xlsx:Path,audit:Path)->dict:
    original=extract_original_xlsx(xlsx,audit)
    sections,source=to_original_sites(original)
    values={}
    for mode in MODES:
        values[mode]={}
        for radius in RADII:
            label=("ALL" if math.isinf(radius) else f"{radius:g}km")
            values[mode][label]=summarise(sections,mode,radius)
    bins=distance_bins(sections)
    principal=values["all"]["ALL"]
    if not principal["identified"] or principal["source_real_rail_sections"]!=14 or principal["source_eligible_original_spatial_sites"]!=88:
        raise ValueError("Original Slovenian 14section 88site nearest other SOURCE rail line cohort drift")
    if abs(principal["original_section_equal_same_site_species_similarity"]-.3632363548846349)>1e-10:
        raise ValueError("Source previously proven genuine original BS/BT same-site Jaccard similarity changed")
    return {
        "version":"v27.7",
        "source_original_rail_vegetation_DOI":"10.1007/s44473-026-00155-z",
        "source_geographic_and_official_rail_line_data":source,
        "source_true_rail_route_kilometrage_nearest_different_site_within_line":values,
        "source_original_section_equal_cross_zone_species_similarity_by_route_km":bins,
        "true_railway_main_branch_or_traffic_NOT_causal_inference":True,
        "source_plant_taxon_GPS_or_original_location_ids_NOT_exported":True,
        "both_original_BS_BT_source_zones_are_railway_habitats":True,
        "zero_distance_site_vs_nonzero_nearest_does_not_separate_environment_from_seed_supply":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",type=Path,required=True)
    p.add_argument("--source-audit",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    out=run(args.original_xlsx,args.source_audit)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for mode in MODES:
        for radius in RADII:
            panel=out["source_true_rail_route_kilometrage_nearest_different_site_within_line"][mode]["ALL" if math.isinf(radius) else f"{radius:g}km"]
            print(json.dumps({
                "version":"v27.7","source_mode":mode,
                "within_source_route_radius":"ALL" if math.isinf(radius) else f"{radius:g}km",
                "identified":panel["identified"],"sections":panel.get("source_real_rail_sections"),
                "source_actual_sites":panel.get("source_eligible_original_spatial_sites"),
                "source_same_site_similarity":panel.get("original_section_equal_same_site_species_similarity"),
                "source_nearest_other_site_similarity":panel.get("original_section_equal_nearest_other_site_species_similarity"),
                "source_same_minus_nearest":panel.get("original_section_equal_same_site_minus_nearest_species_similarity"),
                "n_source_sections_same_gt_nearest":panel.get("source_original_line_sections_with_positive_same_nearest_gap"),
                "source_route_km_median":panel.get("median_of_section_median_nearest_route_km")
            },ensure_ascii=False))
    print("V277_DISTANCE_BINS="+json.dumps(out["source_original_section_equal_cross_zone_species_similarity_by_route_km"],ensure_ascii=False))


if __name__=="__main__":
    main()

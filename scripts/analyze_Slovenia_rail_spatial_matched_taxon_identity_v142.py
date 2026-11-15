#!/usr/bin/env python3
"""v14.2: spatial-distance matched Slovenia real-source same-section plant
taxonomic identity excess beyond observed functional composition.

Original 2026 Slovenian railway field source, author Šabić et al.
DOI 10.1007/s44473-026-00155-z.
89 actual original matched BS/BT source locations; 15 route sections.
Use original publisher XLSX START/END GPS (lat,lon), source recorded
start/end railway chainage kilometers for identity and spatial support.

Predeclare geodesic distance bins (km):
 [0,5),[5,10),[10,20),[20,40),[40,80),[80,160),[160,500).
Predeclare 4 guild functional-overlap bins:
 [0,.7),[.7,.8),[.8,.9),[.9,1.000001).
No flexible bin retuning after seeing the source result.
A bin qualifies only if >=3 SAME-section and >=3 DIFFERENT-section
pairs exist. Compare within/between pairs WITHIN EACH identical
bin, weighted by min(within,between) source pair-opportunity count.
Other pairs are nonoverlap, excluded explicitly, not zero outcomes.

The fixed 4-guild checkerboard swaps from v14.1 conserve within-source
site functional profile and species original field frequency, optionally
within original survey-year × SINGLE/DOUBLE-track strata. Apply the
same FIXED geographic and functional pair windows to every null draw.
Two random-walk null chains are correlated. They are not randomized
train experiments. Likewise source same-section binary only = author
railway segment membership, not physical rail graph/path/traffic.
This is a strictly descriptive one-country source spatial identity
clustering and common-support diagnostic, NEVER proof of carriage.

No road plots or rail vs non-rail comparisons exist in this source.
Geographical matching adjusts physical distance BINS but not detailed
geology, climate, local habitat, plant pools or freight operations.
Strict 2022/23 field vegetation selection/unknown management applies.
No new observations, no within-species evolved railway stress adaptation.
"""
from __future__ import annotations
import argparse,copy,json,math,random
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path
from statistics import median,mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_DOI,SOURCE_SHA256,extract_original_xlsx,validate_original_keys
)
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_sites_and_guilds,functional_overlap,jaccard,
    make_matrices,swap_proposals,matrix_taxon_lists
)

DISTANCE_CUTS_KM=(0,5,10,20,40,80,160,500)
FUNCTIONAL_OVERLAP_CUTS=(0,.7,.8,.9,1.000001)
MIN_OPPORTUNITIES_EACH_ARM=3
DRAWS_PER_CHAIN=160
CHAINS=2
BURNIN=24000
THIN_PROPOSALS=5000
SEED=20261115


def original_gps_coordinates(original,sites):
    locs,_,_,_,_,_=validate_original_keys(original)
    byname={str(k):v for k,v in locs.items()}
    attached=[]
    transect_lengths=[]
    for one in sites:
        g=byname[one["site"]]
        fields=("start_lat","start_lon","end_lat","end_lon")
        try:c=[float(g[k]) for k in fields]
        except (TypeError,ValueError,KeyError) as e:
            raise ValueError("Original Slovenian XLSX site endpoint GPS missing") from e
        a,b,cend,d=c
        if not all(math.isfinite(x) for x in (a,b,cend,d)):
            raise ValueError("Original publisher site GPS values non-finite")
        if not(44<=a<=48 and 44<=cend<=48 and 12<=b<=19 and 12<=d<=19):
            raise ValueError("Original publisher WGS84 coordinate outside broad Slovenia physical bounds")
        p=((a+cend)/2.,(b+d)/2.)
        sl=haversine_km(a,b,cend,d)
        transect_lengths.append(sl)
        if sl>5:
            raise ValueError("Original author source 100m transect endpoints more than 5 km apart; unsafe GIS")
        km=None
        try:
            start,end=float(g["start_km"]),float(g["end_km"])
            if math.isfinite(start) and math.isfinite(end) and start>=0 and end>=0:
                km=(start+end)/2
        except (KeyError,TypeError,ValueError):pass
        attached.append({**one,
          "start_source_lat_lon":(a,b),
          "end_source_lat_lon":(cend,d),
          "lat":p[0],"lon":p[1],"mid_source_route_chainage_km":km})
    if len(attached)!=89:raise ValueError("GPS source cohort must preserve 89 real locations")
    return attached,{
        "verified_original_geographical_WGS84_site_GPS_endpoint_records":len(attached),
        "n_valid_source_mid_kilometrage":sum(x["mid_source_route_chainage_km"] is not None for x in attached),
        "original_author_100m_transect_endpoint_distance_m_median":1000*median(transect_lengths),
        "original_author_100m_transect_endpoint_distance_m_max":1000*max(transect_lengths),
        "source_latitude_minmax":[min(s["lat"] for s in attached),max(s["lat"] for s in attached)],
        "source_longitude_minmax":[min(s["lon"] for s in attached),max(s["lon"] for s in attached)],
        "no_unreported_road_distance_or_freight_flow_converted_from_GPS":True,
    }


def haversine_km(lat1,lon1,lat2,lon2):
    p1,p2=map(math.radians,(lat1,lat2))
    dlat=p2-p1
    dlon=math.radians(lon2-lon1)
    a=math.sin(dlat/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dlon/2)**2
    return 6371.0088*2*math.atan2(math.sqrt(a),math.sqrt(max(0,1-a)))


def into_bin(v,cuts):
    if not math.isfinite(v):return None
    for i,(lo,hi) in enumerate(zip(cuts[:-1],cuts[1:])):
        if lo<=v<hi:return i
    return None


def source_pairs(sites,guild):
    out=[]
    for i,j in combinations(range(len(sites)),2):
        x,y=sites[i],sites[j]
        geo=haversine_km(x["lat"],x["lon"],y["lat"],y["lon"])
        f=functional_overlap(x["known"],y["known"],guild)
        gb=into_bin(geo,DISTANCE_CUTS_KM)
        fb=into_bin(f,FUNCTIONAL_OVERLAP_CUTS)
        if gb is None or fb is None:
            raise ValueError("Original source 89 sites contain unexpected geodesic or functional similarity")
        same=x["section"]==y["section"]
        chainage=(abs(x["mid_source_route_chainage_km"]-
                     y["mid_source_route_chainage_km"])
                  if same and x["mid_source_route_chainage_km"] is not None
                  and y["mid_source_route_chainage_km"] is not None else None)
        out.append({
            "i":i,"j":j,"same":same,"geo_km":geo,"fun":f,
            "geo_bin":gb,"functional_bin":fb,
            "same_year":x["year"]==y["year"],
            "same_track_type":x["track_type"]==y["track_type"],
            "route_km_separation_same_source_section":chainage,
            "observed_known_species_shared_jaccard":jaccard(x["known"],y["known"]),
        })
    if len(out)!=3916 or sum(r["same"] for r in out)!=340:
        raise ValueError("Original 89-site / 15 section fixed 3916 pairs, 340 same section pairs source drift")
    return out


def spatial_common_support(pairs,sets=None,mode="geographic_plus_functional",min_each=MIN_OPPORTUNITIES_EACH_ARM):
    if mode not in ("geographic_only","geographic_plus_functional","geographic_functional_same_year_track"):
        raise ValueError("Source predeclared support mode unknown")
    grouped=defaultdict(lambda:[[],[]])
    for r in pairs:
        k=(r["geo_bin"],)
        if mode!="geographic_only":k+=(r["functional_bin"],)
        if mode=="geographic_functional_same_year_track":
            k+=(r["same_year"],r["same_track_type"])
        sim=(r["observed_known_species_shared_jaccard"] if sets is None
             else jaccard(sets[r["i"]],sets[r["j"]]))
        grouped[k][0 if r["same"] else 1].append(sim)
    common={k:v for k,v in grouped.items()
            if min(len(v[0]),len(v[1]))>=min_each}
    if not common:
        return {
            "identified":False,
            "reason":"NO original same source section versus different-section comparisons in identical geographical/functional strata with >=3 each",
            "source_support":summarize_support(pairs,common,mode),
            "weighted_source_identity_gap":None,
        }
    weights={k:min(len(v[0]),len(v[1])) for k,v in common.items()}
    w_total=sum(weights.values())
    val=sum(weights[k]*(mean(a)-mean(b)) for k,(a,b) in common.items())/w_total
    return {
        "identified":True,
        "weighted_source_identity_gap":val,
        "source_support":summarize_support(pairs,common,mode),
        "component_source_bins":[
          {"match_stratum_key":list(k),
           "within_pairs":len(v[0]),"other_pairs":len(v[1]),
           "source_pair_weight":weights[k],
           "source_same_section_identity_mean":mean(v[0]),
           "source_different_section_identity_mean":mean(v[1]),
           "source_identity_gap":mean(v[0])-mean(v[1])}
          for k,v in sorted(common.items(),key=lambda x:str(x[0]))
        ],
    }


def summarize_support(pairs,common,mode):
    within=sum(r["same"] for r in pairs)
    between=len(pairs)-within
    win=sum(len(v[0]) for v in common.values())
    b=sum(len(v[1]) for v in common.values())
    return {
       "matching_design":mode,
       "matching_original_geography_distance_breakpoints_km":list(DISTANCE_CUTS_KM),
       "matching_original_guild_profile_similarity_breakpoints":list(FUNCTIONAL_OVERLAP_CUTS),
       "predeclared_min_site_pairs_per_comparison_arm_per_stratum":MIN_OPPORTUNITIES_EACH_ARM,
       "n_original_eligible_geography_functional_strata":len(common),
       "n_total_same_section_pairs":within,
       "n_total_different_section_pairs":between,
       "n_same_section_pairs_on_identical_source_support":win,
       "n_different_section_pairs_on_identical_source_support":b,
       "proportion_same_section_pairs_represented":win/within if within else None,
       "source_pair_contrast_aggregate_weight_mass":
           sum(min(len(v[0]),len(v[1])) for v in common.values()),
       "original_site_pairs_NOT_independent_replicates":True,
    }


def source_geo_confounds(pairs):
    for_same=[r["geo_km"] for r in pairs if r["same"]]
    for_other=[r["geo_km"] for r in pairs if not r["same"]]
    return {
      "source_same_section_pair_geodesic_distance_km_median":median(for_same),
      "source_different_section_pair_geodesic_distance_km_median":median(for_other),
      "source_same_section_geodesic_distance_km_minmax":[min(for_same),max(for_same)],
      "source_different_section_geodesic_distance_km_minmax":[min(for_other),max(for_other)],
      "source_same_section_guild_profile_similarity_median":
         median(r["fun"] for r in pairs if r["same"]),
      "source_different_section_guild_profile_similarity_median":
         median(r["fun"] for r in pairs if not r["same"]),
      "n_original_source_same_section_pair_with_rail_kilometrage":sum(
         r["same"] and r["route_km_separation_same_source_section"] is not None
         for r in pairs),
    }


def spatial_null(sites,guild,pairs,mode,draws_per_chain=DRAWS_PER_CHAIN,
                 nchains=CHAINS,burnin=BURNIN,spacing=THIN_PROPOSALS,seed=SEED):
    original=spatial_common_support(pairs,mode=mode)
    if not original["identified"]:
        return {"identified":False,"original_support":original["source_support"],
                "reason":original["reason"]}
    for stratify in (False,True):
        # Report both original source guild occupancy fixed-margin
        # controls and stricter year×track conditional exchangeability.
        vals=[]
        acceptance=[]
        for chain in range(nchains):
            mats=make_matrices(sites,guild)
            rng=random.Random(seed+9901*chain+(7 if stratify else 0))
            acc=swap_proposals(mats,sites,rng,burnin,stratify)
            for draw in range(draws_per_chain):
                acc+=swap_proposals(mats,sites,rng,spacing,stratify)
                ss=matrix_taxon_lists(mats,len(sites))
                res=spatial_common_support(pairs,ss,mode=mode)
                if not res["identified"] or res["source_support"]!=original["source_support"]:
                    raise ValueError("Original geographic source pair support must never drift under null swaps")
                vals.append(res["weighted_source_identity_gap"])
            total=burnin+spacing*draws_per_chain
            acceptance.append({"chain":chain+1,
                               "proposals":total,
                               "accepted_source_checkerboard_swaps":acc,
                               "fraction":acc/total})
        if not vals:raise ValueError("Expected real original source null draws missing")
        sorted_vals=sorted(vals)
        def q(p):
            u=p*(len(sorted_vals)-1);lo=math.floor(u);hi=math.ceil(u)
            return sorted_vals[lo]+(sorted_vals[hi]-sorted_vals[lo])*(u-lo)
        yield_mode=("source_guild_species_degree_fixed_ALL_STUDY"
                    if not stratify else "source_guild_species_degree_fixed_WITHIN_YEAR_TRACK_TYPE")
        yield yield_mode,{
          "observed_spatially_matched_taxonomic_identity_section_gap":
             original["weighted_source_identity_gap"],
          "null_expected_spatially_matched_gap":mean(vals),
          "excess_above_guild_frequency_null":
             original["weighted_source_identity_gap"]-mean(vals),
          "null_MC_2_5":q(.025),"null_MC_97_5":q(.975),
          "MC_upper_tail_fraction_plus_one_not_field_causal_p":
             (1+sum(v>=original["weighted_source_identity_gap"] for v in vals))/(1+len(vals)),
          "n_correlated_chain_draws":len(vals),
          "chain_mixing_diagnostics":acceptance,
          "year_track_type_frequencies_preserved":stratify,
          "not_a_physical_rail_graph_distance_or_freight_seed_transport_trial":True,
        }


def run(original):
    sites,guild,prior=source_sites_and_guilds(original)
    coords,coordproof=original_gps_coordinates(original,sites)
    pairs=source_pairs(coords,guild)
    support={m:spatial_common_support(pairs,mode=m)
             for m in ("geographic_only","geographic_plus_functional",
                       "geographic_functional_same_year_track")}
    primary=support["geographic_plus_functional"]
    if not primary["identified"]:
        null={"identified":False,
              "reason":"geographic_plus_functional matching not possible with predeclared bins"}
    else:
        null=dict(spatial_null(coords,guild,pairs,"geographic_plus_functional"))
    loo=[]
    for sec in sorted({s["section"] for s in coords}):
        include=[s for s in coords if s["section"]!=sec]
        # source pairs must be freshly recomputed with different n
        # rather than using unchanged 89 index ids.
        mini=[]
        for i,j in combinations(range(len(include)),2):
            a,b=include[i],include[j]
            geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
            fun=functional_overlap(a["known"],b["known"],guild)
            mini.append({"i":i,"j":j,"same":a["section"]==b["section"],
                         "geo_bin":into_bin(geo,DISTANCE_CUTS_KM),
                         "functional_bin":into_bin(fun,FUNCTIONAL_OVERLAP_CUTS),
                         "same_year":a["year"]==b["year"],
                         "same_track_type":a["track_type"]==b["track_type"],
                         "observed_known_species_shared_jaccard":
                            jaccard(a["known"],b["known"])})
        res=spatial_common_support(mini,mode="geographic_plus_functional")
        loo.append({"excluded_source_section":sec,
                    "original_site_count_remaining":len(include),
                    "estimated":res["identified"],
                    "source_section_identity_gap_after_exclusion":
                        res["weighted_source_identity_gap"]})
    return {
      "complete":True,"version":"v14.2",
      "original_Slovenia_field_source_DOI":SOURCE_DOI,
      "source_2026_plant_supplement_SHA256":SOURCE_SHA256,
      "original_v141_source_field_cohort_reproduced":{
          "source_sites":len(coords),"original_sections":len({s["section"] for s in coords}),
          "actual_dual_known_guild_taxa":prior[
              "original_strategy_lifespan_known_taxa_present_in_matched_sites"]},
      "author_original_100m_GPS_transect_coordinate_audit":coordproof,
      "source_geo_site_section_spatial_confounding":source_geo_confounds(pairs),
      "predeclared_spatial_and_functional_matching_results":support,
      "source_null_spatially_matched_identity":null,
      "leave_one_source_rail_section_out_spatially_matched":{
          "n_source_sections":len(loo),
          "n_estimable":sum(r["estimated"] for r in loo),
          "n_positive":sum(r["estimated"] and
                            r["source_section_identity_gap_after_exclusion"]>0
                            for r in loo),
          "individual_source_leaveout_results":loo,
      },
      "scientific_limits":[
        "This uses true published ORIGINAL start/end GPS coordinates to match geodesic distance, but not actual railway-network route distance, continuity or train activity.",
        "Same ORIGINAL route section means shared author source line label. It is not freight connectivity, not shipping OD, not an observation of seeds or dispersal.",
        "Geo bins must overlap on BOTH same-section and different-section source pairs; unsupported same-section pairs are excluded as nonidentifiable rather than modeled outside source overlap.",
        "Four coarse source lifespan×diaspore strategy archetypes are not a measured functional beta hypervolume, and the author source study has no non-RAIL field survey control.",
        "Checkerboard conditional null preserves guild species counts per site and species frequency, with year/track type as sensitivity. It does not preserve regional species pools, environment or landscape features.",
        "All 3916 pairs share original 89 sites nested in only 15 source railway line segments and one country. MC swaps are not independent scientific studies or exact causal p-values.",
        "No actual railway network distance, source freight volumes, train schedules or road parallel corridor distances computed here; no train seed transport causal claim.",
        "Original observational site selection preferred visible vegetation and survey years 2022/23 after herbicide changes; cross-sectional data cannot infer temporal colonization, plasticity or adaptation.",
      ],
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    raw=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    out=run(raw)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    summ={
      "complete":True,"version":"v14.2",
      "GPS":out["author_original_100m_GPS_transect_coordinate_audit"],
      "geographic_source_confounding":out["source_geo_site_section_spatial_confounding"],
      "matched":{m:{"identified":v["identified"],"result":v.get("weighted_source_identity_gap"),
                     "support":v["source_support"]}
                 for m,v in out["predeclared_spatial_and_functional_matching_results"].items()},
      "conditional_null":out["source_null_spatially_matched_identity"],
      "section_loo":{k:v for k,v in out["leave_one_source_rail_section_out_spatially_matched"].items()
                     if k!="individual_source_leaveout_results"},
      "do_not_conflate_section_label_with_train_seed_dispersal":True,
    }
    print(json.dumps(summ,ensure_ascii=False),flush=True)


if __name__=="__main__":main()

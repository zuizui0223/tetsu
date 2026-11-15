#!/usr/bin/env python3
"""v15.4 independent REAL Czech railway botanical survey:
geographic conditioning of same-stretch species sharing using exact station
name identity and externally derived 2026 SR70 railway-station CENTROIDS.

Authors Kutlvašr et al. (2024), Preslia DOI 10.23855/preslia.2024.247,
39 field stations, 2020-22, S1 current station/taxon incidence and corrected
S2 metadata source-pinned in v14.3. Pladias source strategy+lifeform with
non-imputation. Full names matched conservatively v15.3 (35 of 39).
Czech coordinate source is third-party-derived official 2026 SR70 point, NOT
the author field GPS, railway-freight route, or geolocation of each plant.
Exactly 4 unlocated source stations remain unknown, not fabricated.

Use 35 source sites as independent botanical units and 5 original route
stretches as source clusters, never 595 independent field studies.
Match original same-stretch vs different-stretch source pairs by
FROZEN original Slovenian spatial bins and source guild-composition bins,
optionally independent source station size. Do not pool country's effect.

Null swaps preserve each source station's source functional-guild plant
richness and each original taxon's occupancy, with descriptive stricter
within-original-stretch swaps (which preserves geographic source pools).
The latter can overcondition the mechanism and is a sensitivity only.
No causal seed carriage, mowing, aerodynamic adaptation or global beta.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median
from analyze_Czech_2024_39station_functional_identity_v143 import (
    read_original_czech_publisher,current_sites,site_pairs,
    function_sim,jaccard,GROUPS
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import parse_traits
from audit_Czech_2024_S2_to_SR70_original_station_points_v153 import (
    original_s2_names,sr70_records,match_czech_stations
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    haversine_km,into_bin,DISTANCE_CUTS_KM,FUNCTIONAL_OVERLAP_CUTS
)

AREA_LOG_RATIO_CUTS=(0.,.5,1.,2.,20.)
MIN_SAME_AND_DIFF_PAIRS_EACH_STRATUM=3
NCHAINS=2
DRAWS_PER_CHAIN=130
BURNIN=16000
THIN=3000
SEED=20261115
EXPECTED_ORIGINAL_STATION_MATCHES=35
EXPECTED_PUBLISHER_ALL_STATIONS=39
EXPECTED_PUBLISHER_VALID_TAXA=544
EXPECTED_PUBLISHER_VALID_SITE_TAXON_INCIDENCE=5081


def original_35_station_sites(s1,s2,dispersal,lifeform,sr70):
    floras,meta,src=read_original_czech_publisher(s1,s2)
    d,l,traits_prov=parse_traits(dispersal,lifeform)
    original_sites,guild,cov=current_sites(floras,meta,d,l)
    if len(original_sites)!=39 or len(guild)!=EXPECTED_PUBLISHER_VALID_TAXA:
        raise ValueError("Publisher original 39 stations / 544 source matched Pladias taxa drifted")
    if cov["n_original_incidence_with_both_exact_Pladias_traits"]!=EXPECTED_PUBLISHER_VALID_SITE_TAXON_INCIDENCE:
        raise ValueError("Publisher original exact 5081 Pladias matched incidence drifted")
    sr,spatialprov=sr70_records(sr70)
    exact=match_czech_stations(original_s2_names(s2),sr)
    lookup={s["original_locality_code"]:s for s in exact}
    if len(lookup)!=39:
        raise ValueError("Original S2 site code coverage unexpectedly changed")
    joined=[]
    excluded=[]
    for s in original_sites:
        original=lookup[s["site"]]
        loc=original["matched_czech_2026_source_station_centroid_WGS84"]
        if loc is None:
            excluded.append({
                "original_station_code":s["site"],
                "original_station_name":original["original_source_s2_col1_label"],
                "reason":original["match_status"],
            })
            continue
        if not (s["section"]==original["source_section"]):
            raise ValueError("Original S2 section label inconsistent with publisher station flora")
        joined.append({**s,"lat":loc["lat"],"lon":loc["lon"],
                       "station_name":original["original_source_s2_col1_label"],
                       "matched_official_derived_sr70_codes":loc["source_sr70_registry_codes"]})
    if len(joined)!=EXPECTED_ORIGINAL_STATION_MATCHES or len(excluded)!=4:
        raise ValueError("Unexpected exact original-name station geographic source availability")
    if len({x["site"] for x in joined})!=len(joined):
        raise ValueError("Duplicate source station after matching")
    return joined,guild,{"original_botanical_publisher_source":src,
          "source_Pladias_trait_coverage":cov,
          "Pladias_source_read":traits_prov,
          "2026_SR70_source_provenance":spatialprov,
          "n_original_39_station_centroids_matched_by_exact_full_name":len(joined),
          "unmatched_4_authentic_original_station_codes_not_imputed":excluded,
          "not_source_vegetation_plot_coordinates":True,
          "field_2020_2022_vs_SR70_2026_time_mismatch":True}


def real_site_pairs(sites,guild):
    out=[]
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
        fun=function_sim(a["known"],b["known"],guild)
        if fun is None or not(geo>0 and math.isfinite(geo)):
            raise ValueError("Czech source sites not geographic distinct or missing original guild")
        gi=into_bin(geo,DISTANCE_CUTS_KM)
        fi=into_bin(fun,FUNCTIONAL_OVERLAP_CUTS)
        area_diff=(abs(math.log(a["area_ha"])-math.log(b["area_ha"]))
                   if a["area_ha"] is not None and b["area_ha"] is not None
                   and a["area_ha"]>0 and b["area_ha"]>0 else None)
        ai=into_bin(area_diff,AREA_LOG_RATIO_CUTS) if area_diff is not None else None
        out.append({"i":i,"j":j,
                    "same_stretch":a["section"]==b["section"],
                    "geodesic_station_centroid_km":geo,
                    "functional_coarse_overlap":fun,
                    "geo_bin":gi,"function_bin":fi,"area_log_diff_bin":ai,
                    "taxonomic_Jaccard_original_544":jaccard(a["known"],b["known"]),
                    "taxonomic_Jaccard_ALL_author_source_761":jaccard(
                       a["all_original_taxa"],b["all_original_taxa"])})
    if len(out)!=595:
        raise ValueError("35 original stations must yield exactly 595 dependent pairs")
    return out


def source_gap(pairs,known=True):
    k="taxonomic_Jaccard_original_544" if known else "taxonomic_Jaccard_ALL_author_source_761"
    same=[p[k] for p in pairs if p["same_stretch"]]
    other=[p[k] for p in pairs if not p["same_stretch"]]
    if not same or not other:raise ValueError("Original Czech study sections have no comparison")
    return {"n_same_source_section_pairs":len(same),
            "n_different_source_section_pairs":len(other),
            "same_stretch_jaccard":mean(same),
            "different_stretch_jaccard":mean(other),
            "same_minus_different_jaccard":mean(same)-mean(other)}


def conditional_source_gap(pairs,sets=None,mode="geography_and_function",
                           min_each=MIN_SAME_AND_DIFF_PAIRS_EACH_STRATUM):
    if mode not in ("geography_only","geography_and_function","geography_function_and_area"):
        raise ValueError("Undeclared Czech source original station centroid matching method")
    fields=("geo_bin",) if mode=="geography_only" else ("geo_bin","function_bin")
    if mode=="geography_function_and_area":fields+=("area_log_diff_bin",)
    buckets=defaultdict(lambda:[[],[]])
    unknown_area=0
    for p in pairs:
        key=tuple(p.get(k) for k in fields)
        if None in key:
            unknown_area+=1
            continue
        j=(p["taxonomic_Jaccard_original_544"] if sets is None
           else jaccard(sets[p["i"]],sets[p["j"]]))
        if j is None:raise ValueError("Source taxonomic null lost all species from a station")
        buckets[key][0 if p["same_stretch"] else 1].append(j)
    common={key:v for key,v in buckets.items()
            if min(len(v[0]),len(v[1]))>=min_each}
    within=sum(len(v[0]) for v in common.values())
    between=sum(len(v[1]) for v in common.values())
    weight=sum(min(len(v[0]),len(v[1])) for v in common.values())
    source={"supported_strata_count":len(common),
            "within_source_stretch_pairs_on_common_support":within,
            "different_source_stretch_pairs_on_common_support":between,
            "matched_source_pair_weight":weight,
            "excluded_missing_area_or_outside_original_declared_bins":unknown_area,
            "no_extrapolation_outside_common_support":True}
    if not common:
        return {"identified":False,
                "reason":"Original source 35 station centroids do not yield >=3 within and between section pairs per matching stratum",
                "source_support":source,
                "weighted_same_minus_different_jaccard":None}
    val=sum(min(len(a),len(b))*(mean(a)-mean(b))
            for a,b in common.values())/weight
    return {"identified":True,"source_support":source,
            "weighted_same_minus_different_jaccard":val}


def make_source_matrices(sites,guild):
    return {g:[x["known"] & {tax for tax,c in guild.items() if c==g}
               for x in sites] for g in GROUPS}


def swap(m,sites,rng,proposal_count,within_section):
    strata=defaultdict(list)
    for i,s in enumerate(sites):
        strata[s["section"] if within_section else "full"].append(i)
    allowed=[v for v in strata.values() if len(v)>=2]
    if not allowed:raise ValueError("No source original rail section contains two botanical stations")
    accepted=0
    for _ in range(proposal_count):
        idx=rng.choice(allowed)
        i,j=rng.sample(idx,2)
        rows=m[rng.choice(GROUPS)]
        left=rows[i]-rows[j];right=rows[j]-rows[i]
        if not left or not right:continue
        a=rng.choice(tuple(left));b=rng.choice(tuple(right))
        rows[i].remove(a);rows[i].add(b)
        rows[j].remove(b);rows[j].add(a)
        accepted+=1
    return accepted


def reconstructed_sets(m,n):
    result=[set() for _ in range(n)]
    for rows in m.values():
        for i,s in enumerate(rows):result[i].update(s)
    return result


def fixed_source_margin_null(sites,guild,pairs,mode="geography_and_function",
                             within_stretch=False,nchains=NCHAINS,draws=DRAWS_PER_CHAIN,
                             burnin=BURNIN,thin=THIN,seed=SEED):
    primary=conditional_source_gap(pairs,mode=mode)
    if not primary["identified"]:
        return {"identified":False,"reason":primary["reason"],
                "source_support":primary["source_support"]}
    original=primary["weighted_same_minus_different_jaccard"]
    rng_results=[];diagnostics=[]
    for chain in range(nchains):
        m=make_source_matrices(sites,guild)
        expected_count={g:[len(row) for row in m[g]] for g in GROUPS}
        source_prevalence=Counter(x for rows in m.values() for site in rows for x in site)
        within_group_prevalence=Counter(
            (sites[i]["section"],x) for rows in m.values()
            for i,site in enumerate(rows) for x in site)
        rng=random.Random(seed+1009*chain+(31 if within_stretch else 0))
        acc=swap(m,sites,rng,burnin,within_stretch)
        for _ in range(draws):
            acc+=swap(m,sites,rng,thin,within_stretch)
            if any([len(row) for row in m[g]]!=expected_count[g] for g in GROUPS):
                raise ValueError("Czech source each station×guild original richness not conserved")
            if Counter(x for rows in m.values() for site in rows for x in site)!=source_prevalence:
                raise ValueError("Czech source species original site frequency not conserved")
            if within_stretch and Counter(
                    (sites[i]["section"],x) for rows in m.values()
                    for i,site in enumerate(rows) for x in site)!=within_group_prevalence:
                raise ValueError("Original study stretch specific species prevalence not conserved")
            drawn=conditional_source_gap(pairs,reconstructed_sets(m,len(sites)),mode=mode)
            if not drawn["identified"] or drawn["source_support"]!=primary["source_support"]:
                raise ValueError("Original matched geographic pair window shifted under null")
            rng_results.append(drawn["weighted_same_minus_different_jaccard"])
        diagnostics.append({"source_chain":chain+1,
                            "accepted_checkerboard_swaps":acc,
                            "proposed":burnin+draws*thin})
    z=sorted(rng_results)
    return {"identified":True,
            "source_matched_35_station_observed_jaccard_gap":original,
            "null_mean":mean(rng_results),
            "observed_minus_null_mean":original-mean(rng_results),
            "mcmc_percentile_2_5":z[round(.025*(len(z)-1))],
            "mcmc_percentile_97_5":z[round(.975*(len(z)-1))],
            "n_correlated_source_mcmc_draws":len(z),
            "fixed_global_species_prevalence_and_station_guild_richness":True,
            "within_original_stretch_species_prevalence_additionally_fixed":within_stretch,
            "within_stretch_fixed_sensitivity_can_overcondition_true_effect":True,
            "no_per_pair_experimental_p":True,
            "source_chain_diagnostics":diagnostics}


def section_leaveout(sites,guild):
    out=[]
    for group in sorted({x["section"] for x in sites}):
        original=[s for s in sites if s["section"]!=group]
        pairs=real_site_pairs_without_n_gate(original,guild)
        stat=conditional_source_gap(pairs,mode="geography_and_function")
        out.append({"omitted_original_railway_stretch":group,
                    "remaining_source_stations":len(original),
                    "identified":stat["identified"],
                    "geographic_guild_conditioned_jaccard_gap":
                         stat["weighted_same_minus_different_jaccard"],
                    "source_support":stat["source_support"]})
    return out

def real_site_pairs_without_n_gate(sites,guild):
    # Rebuild the original same-source method without pretending a
    # 30-site/4-stretch sensitivity should also have 595 pairs.
    pairs=[]
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        geo=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
        f=function_sim(a["known"],b["known"],guild)
        area=(abs(math.log(a["area_ha"])-math.log(b["area_ha"]))
             if all(x is not None and x>0 for x in (a["area_ha"],b["area_ha"])) else None)
        pairs.append({"i":i,"j":j,
                      "same_stretch":a["section"]==b["section"],
                      "geo_bin":into_bin(geo,DISTANCE_CUTS_KM),
                      "function_bin":into_bin(f,FUNCTIONAL_OVERLAP_CUTS),
                      "area_log_diff_bin":into_bin(area,AREA_LOG_RATIO_CUTS) if area is not None else None,
                      "taxonomic_Jaccard_original_544":jaccard(a["known"],b["known"])})
    return pairs


def run(s1,s2,disp,life,sr70):
    sites,guild,source=original_35_station_sites(s1,s2,disp,life,sr70)
    pairs=real_site_pairs(sites,guild)
    geo_same=[x["geodesic_station_centroid_km"] for x in pairs if x["same_stretch"]]
    geo_other=[x["geodesic_station_centroid_km"] for x in pairs if not x["same_stretch"]]
    geocons={m:conditional_source_gap(pairs,mode=m) for m in (
        "geography_only","geography_and_function","geography_function_and_area")}
    primary=geocons["geography_and_function"]
    null={}
    if primary["identified"]:
        for restricted in (False,True):
            key="within_source_stretch_strict_frequency" if restricted else "whole_35site_source_global_frequency"
            null[key]=fixed_source_margin_null(sites,guild,pairs,within_stretch=restricted)
    else:
        null={"identified":False,"reason":"Primary source geographic+guild matching no common support"}
    return {
        "version":"v15.4_REAL_Czech_original_station_species_sharing_geographic_source",
        "original_study_country":"Czech Republic",
        "original_source_DOr_study_ID":"10.23855/preslia.2024.247",
        "source_2020_2022_plant_publisher_and_2026_station_centroid":source,
        "n_strict_exact_station_identity_matched":len(sites),
        "n_original_source_stretches":len({s["section"] for s in sites}),
        "n_dependent_real_station_pairs":len(pairs),
        "n_same_source_stretch_pairs":len(geo_same),
        "n_other_stretch_pairs":len(geo_other),
        "same_stretch_station_centroid_geodesic_km_median":median(geo_same),
        "other_stretch_station_centroid_geodesic_km_median":median(geo_other),
        "raw_544_both_trait_known_source_taxonomic_jaccard":source_gap(pairs),
        "raw_ALL_761_source_taxonomic_jaccard":source_gap(pairs,known=False),
        "fixed_original_slovenia_geographic_bins_and_Czech_original_trait_area_conditioning":geocons,
        "fixed_source_species_occupancy_checkerboard_null":null,
        "leave_one_original_Czech_stretch_out_not_country_replication":section_leaveout(sites,guild),
        "source_research_bounds":[
            "SR70 2026 official-derived station POINTS not the true 2020-22 original field vegetation plot coordinates. 35/39 exact name matching; 4 not imputed.",
            "Original 39 Czech stations are unequal multi-habitat areas unlike Slovenian 100m vegetation transects; may have source survey effort/area confounding.",
            "Geodesic distance cannot substitute for actual railway line physical graph distance, direction or freight service history; NO actual rail-to-road path yet in this version.",
            "Five original railway stretches are nested in ONE country and species pairs are dependent, no 595 independent ecological experiments.",
            "Geographic and guild bins are fixed from preceding Slovenian source method, not selected to maximize new Czech observed effect.",
            "Strict within-stretch species-frequency fixed swaps may remove a real mechanism along with common plant-pool structure; sensitivity not a causal null.",
            "No non-rail vegetation field samples, real seed-trap germination/colonization, measured train wind, grass cutting treatment or genetics.",
            "The expected original SR70 2026 station source still requires author original GPS for true same-grain crosscountry metaanalysis."
        ],
        "NO_pooled_Czech_Slovenia_beta_or_train_causal_claim":True
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("s1","s2","pladias-diaspore","pladias-lifeform","sr70","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.pladias_diaspore,a.pladias_lifeform,a.sr70)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    compact={k:v for k,v in result.items() if k in (
        "version","n_strict_exact_station_identity_matched","n_original_source_stretches",
        "n_dependent_real_station_pairs","n_same_source_stretch_pairs","n_other_stretch_pairs",
        "same_stretch_station_centroid_geodesic_km_median",
        "other_stretch_station_centroid_geodesic_km_median",
        "raw_544_both_trait_known_source_taxonomic_jaccard",
        "raw_ALL_761_source_taxonomic_jaccard",
        "fixed_original_slovenia_geographic_bins_and_Czech_original_trait_area_conditioning",
        "fixed_source_species_occupancy_checkerboard_null",
        "leave_one_original_Czech_stretch_out_not_country_replication")}
    print(json.dumps(compact,ensure_ascii=False))

if __name__=="__main__":main()

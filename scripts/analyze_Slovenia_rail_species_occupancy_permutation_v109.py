#!/usr/bin/env python3
"""v10.9: Is the Slovenian railway wind-diaspore deficit robust to species identity?

REAL source: Šabić, Grašič & Rovanšek 2026, DOI 10.1007/s44473-026-00155-z.
Reanalyze original 312 species-code traits and 89 actual rail locations:
  - exactly matched BS/BT zone-side pairs pooled WITHIN location,
  - species equal weighted BS-vs-BT occupancy preferences,
  - strict annual vs other-known-lifespan stratified wind-label permutation,
  - uncommon-vs-common species sensitivity,
  - leave-one-section-out and leave-one-taxon-out diagnostics.
Permutation labels are a diagnostic under an exchangeability assumption,
NOT a phylogenetic null, train-carriage test or causal wind/mowing effect.
Original publication field observations (NOT new tetsu field collection).
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_DOI, SOURCE_SHA256, extract_original_xlsx,
    validate_original_keys, pair_source
)

WIND=frozenset(("Epilobium", "Phragmites"))
STRATA=("strict_annual", "other_known_lifespan")
PERMUTATIONS=9999
SEED=20261114
THRESHOLDS=(1,3,5,10)


def source_guild(trait):
    d=trait.get("Diaspore.strategy")
    life=trait.get("Lifespan.categories")
    if d is None or d=="" or life is None or life=="":
        return None
    return {
        "wind":d in WIND,
        "stratum":"strict_annual" if life=="annual" else "other_known_lifespan",
    }


def original_location_matched_sets(original):
    geos,zones,species,traits,codebook,integrity=validate_original_keys(original)
    paired=pair_source(geos,zones,species,traits)
    byloc={}
    for row in paired:
        lid=row["location_id"]
        side=row["track_side"]
        if lid not in byloc:
            byloc[lid]={"section_number":str(row["section_number"]),
                         "BS":set(),"BT":set(),"n_matched_sides":0}
        obj=byloc[lid]
        if obj["section_number"]!=str(row["section_number"]):
            raise ValueError("Original site cannot belong to two railway segments")
        obj["BS"].update(species[(lid,side+"BS")])
        obj["BT"].update(species[(lid,side+"BT")])
        obj["n_matched_sides"]+=1
    if len(byloc)!=89 or len(paired)!=173 or len(set(x["section_number"] for x in byloc.values()))!=15:
        raise ValueError("Original matched 89-location 15-section source drift")
    if any(x["n_matched_sides"] not in (1,2) for x in byloc.values()):
        raise ValueError("Source left/right sides improperly grouped")
    return byloc,traits,integrity


def taxa_preference(source_sites,traits,exclude_section=None):
    accum=defaultdict(lambda: [0,0,0])
    nloc=0
    for site in source_sites.values():
        if site["section_number"]==exclude_section:
            continue
        nloc+=1
        bs=site["BS"]; bt=site["BT"]
        for code in sorted(bs|bt):
            v=accum[code]
            if code in bs and code in bt:v[2]+=1
            elif code in bs:v[0]+=1
            else:v[1]+=1
    out=[]
    for code,v in sorted(accum.items()):
        meta=source_guild(traits[code])
        if meta is None:
            continue
        bs_only,bt_only,both=v
        union=sum(v)
        if union==0:continue
        out.append({
            "taxon_code":code,"wind":meta["wind"],"stratum":meta["stratum"],
            "BS_only":bs_only,"BT_only":bt_only,"both":both,
            "n_original_locations_present_either_zone":union,
            "occupancy_preference_BS_minus_BT":(bs_only-bt_only)/union,
        })
    return out,nloc


def contrast(rows,min_sites,excluded_taxon=None):
    by={s:{"wind":[],"other":[]} for s in STRATA}
    for r in rows:
        if r["taxon_code"]==excluded_taxon or r["n_original_locations_present_either_zone"]<min_sites:
            continue
        by[r["stratum"]]["wind" if r["wind"] else "other"].append(
            r["occupancy_preference_BS_minus_BT"])
    counts={}
    for s in STRATA:
        counts[s]={"wind":len(by[s]["wind"]),"other":len(by[s]["other"])}
    if any(min(counts[s].values())<3 for s in STRATA):
        return {"identified":False,"min_original_locations_per_taxon":min_sites,
                "n_qualifying_species_by_guild":counts}
    n=sum(sum(counts[s].values()) for s in STRATA)
    details={}
    value=0.0
    for s in STRATA:
        diff=mean(by[s]["wind"])-mean(by[s]["other"])
        wt=sum(counts[s].values())/n
        details[s]={"wind_minus_other_species_mean_BS_preference":diff,
                    "n_wind_taxa":counts[s]["wind"],
                    "n_other_taxa":counts[s]["other"],
                    "stratum_taxon_weight":wt,
                    "wind_preference_mean":mean(by[s]["wind"]),
                    "other_preference_mean":mean(by[s]["other"])}
        value+=wt*diff
    return {
        "identified":True,"min_original_locations_per_taxon":min_sites,
        "n_qualifying_species":n,
        "n_qualifying_species_by_guild":counts,
        "equal_species_within_lifespan_wind_minus_other_preference":value,
        "by_lifespan_stratum":details,
    }


def taxon_label_permutation(rows,min_sites,draws=PERMUTATIONS,seed=SEED):
    c=contrast(rows,min_sites)
    if not c["identified"]:
        return c
    groups={}
    for s in STRATA:
        a=[r for r in rows if r["stratum"]==s and
           r["n_original_locations_present_either_zone"]>=min_sites]
        groups[s]={"values":[r["occupancy_preference_BS_minus_BT"] for r in a],
                   "n_wind":sum(r["wind"] for r in a),
                   "n":len(a)}
    obs=c["equal_species_within_lifespan_wind_minus_other_preference"]
    rng=random.Random(seed)
    n_all=sum(groups[s]["n"] for s in STRATA)
    extreme=0
    less=0
    for _ in range(draws):
        stat=0
        for s in STRATA:
            d=groups[s]
            n=d["n"]; nw=d["n_wind"]
            picked=rng.sample(range(n),nw)
            selected=sum(d["values"][i] for i in picked)
            total=sum(d["values"])
            diff=selected/nw-(total-selected)/(n-nw)
            stat+=(n/n_all)*diff
        if abs(stat)+1e-12>=abs(obs):extreme+=1
        if stat<=obs+1e-12:less+=1
    return {**c,
        "source_lifespan_stratified_species_label_permutations":draws,
        "null_unit":"source_plant_species_trait_label_preserving_strict_annual_strata",
        "empirical_two_sided_permutation_p":(extreme+1)/(draws+1),
        "empirical_one_sided_wind_more_BT_p":(less+1)/(draws+1),
        "exchangeability_and_phylogenetic_limitation":
            "Species trait labels are fixed biologically and may cluster by genus/family; "
            "label permutation is a species-level diagnostic, not a causal or phylogenetic null.",
    }


def run(original,verified106):
    if not (verified106.get("complete") is True and
            verified106.get("version")=="v10.6" and
            verified106.get("source_original_xlsx_sha256")==SOURCE_SHA256 and
            verified106.get("original_same_site_means_15_section_bootstrap",{}).get(
                "within_lifespan_wind_component",{}).get(
                "equal_field_location_mean_BS_minus_BT")==-0.040424):
        raise ValueError("Required genuine original Slovenia source v10.6 receipt missing")
    sites,traits,integrity=original_location_matched_sets(original)
    preference,nsite=taxa_preference(sites,traits)
    if nsite!=89:
        raise ValueError("Original site count drift")
    results={str(k):taxon_label_permutation(preference,k) for k in THRESHOLDS}
    if not results["1"]["identified"]:
        raise ValueError("Original wind guild failed basic species-replicate support")
    preferred=3 if results["3"]["identified"] else 1
    original_contrast=contrast(preference,preferred)
    sections=sorted({s["section_number"] for s in sites.values()})
    loo=[]
    for sect in sections:
        reduced,_=taxa_preference(sites,traits,exclude_section=sect)
        c=contrast(reduced,preferred)
        loo.append({"omitted_original_railway_section":sect,
                    "identified":c["identified"],
                    "n_qualifying_species":c.get("n_qualifying_species"),
                    "species_equal_within_lifespan_wind_minus_other_BS_preference":
                        c.get("equal_species_within_lifespan_wind_minus_other_preference")})
    original_negative=original_contrast["equal_species_within_lifespan_wind_minus_other_preference"]<0
    influence=[]
    for t in preference:
        if t["n_original_locations_present_either_zone"]<preferred:continue
        leave=contrast(preference,preferred,excluded_taxon=t["taxon_code"])
        if not leave["identified"]:continue
        influence.append({
            "source_EPPO_taxon_code":t["taxon_code"],
            "guild":"wind" if t["wind"] else "other_source_strategy",
            "lifespan_stratum":t["stratum"],
            "source_locations_present":t["n_original_locations_present_either_zone"],
            "BS_only":t["BS_only"],"BT_only":t["BT_only"],"both":t["both"],
            "species_BS_BT_preference":round(t["occupancy_preference_BS_minus_BT"],6),
            "delta_contrast_from_omitting_this_taxon":round(
                leave["equal_species_within_lifespan_wind_minus_other_preference"]-
                original_contrast["equal_species_within_lifespan_wind_minus_other_preference"],6),
            "sign_after_taxon_leaveout":
                "negative" if leave["equal_species_within_lifespan_wind_minus_other_preference"]<0
                else "positive_or_zero",
        })
    by_influence=sorted(influence,key=lambda x:abs(x["delta_contrast_from_omitting_this_taxon"]),reverse=True)
    windtaxa=[t for t in preference if t["wind"]]
    all_other=[t for t in preference if not t["wind"]]
    return {
        "complete":True,"version":"v10.9",
        "original_article_DOI":SOURCE_DOI,
        "original_publisher_field_XLSX_sha256":SOURCE_SHA256,
        "n_original_railway_field_locations":len(sites),
        "n_original_matched_BS_BT_sides":sum(x["n_matched_sides"] for x in sites.values()),
        "n_original_source_railway_sections":len(sections),
        "n_original_taxon_trait_rows":len(traits),
        "n_source_wind_taxa_with_any_original_BS_BT_occupancy":len(windtaxa),
        "n_source_other_strategy_taxa_with_any_original_BS_BT_occupancy":len(all_other),
        "n_source_taxa_with_both_trait_codes_present":len(preference),
        "species_equal_diagnostics_by_min_site_occupancy":results,
        "primary_predeclared_min_3_or_fallback_1_threshold":preferred,
        "primary_taxon_label_permutation":results[str(preferred)],
        "leave_one_railway_section_out":loo,
        "n_leave_one_section_out_identified":sum(x["identified"] for x in loo),
        "n_leave_one_section_out_negative":sum(
            x["identified"] and x["species_equal_within_lifespan_wind_minus_other_BS_preference"]<0
            for x in loo),
        "n_original_taxa_with_leave_one_out_estimable":len(influence),
        "n_taxa_whose_removal_changes_sign":sum(
            (x["sign_after_taxon_leaveout"]=="negative")!=original_negative
            for x in influence),
        "original_primary_species_equal_contrast_sign":
            "negative" if original_negative else "positive_or_zero",
        "top_20_influential_species_coded_without_full_raw_source_records":by_influence[:20],
        "source_integrity":{
            "original_source_taxa_incidence_deduplicated":integrity[
                "n_unique_original_EPPO_species_times_zone_source_incidence"],
            "site_track_sides_are_pooled_before_species_scoring":True,
            "no_unobserved_original_taxon_zone_presences_added":True,
        },
        "limitation":(
            "This is an observational species equal-weighted robustness diagnostic. "
            "Source category is published Diaspore.strategy, not actual wind or train dispersal. "
            "Study sites nonrandom one-country sample; paired zone sizes differ. "
            "Taxon categories may be phylogenetically clustered and occupation frequencies "
            "confounded by taxonomy/ecology. No train freight seeds, gusts, mowing, "
            "fitness, heritable selection, real species turnover through time "
            "or global rail-vs-road causal effects measured."
        ),
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","verified106","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    source=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    v106=json.loads(Path(a.verified106).read_text(encoding="utf-8"))
    result=run(source,v106)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=="__main__":main()

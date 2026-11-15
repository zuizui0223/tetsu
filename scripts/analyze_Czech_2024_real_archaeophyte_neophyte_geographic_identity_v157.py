#!/usr/bin/env python3
"""v15.7: original publisher Czech author source archaeo-/neo-phyte explicitly
tagged RAILWAY flora identity across 35 strictly located stations.

Exactly original Kutlvašr et al 2024 Preslia appendix S1 station-by-species
REC incidences and author original 'cat.' source status strings. Do not assume
untagged species are confirmed natives. Original historical a/b/c/d source
records are not equivalent observed absence, ignored.
2026 third-party SR70 official-derived exact station centroid names, not
original 2020–22 author vegetation plot GPS (v15.3).
Source geographic/guild/area bins fixed in Czech v15.4, no retuning.

For each explicitly recorded archaeophyte, neophyte, or union:
- 35-site source within-vs-between original stretch Jaccard,
- same fixed geographic+coarse guild similarity pair matching,
- station-area matching sensitivity,
- site class richness and class-species global occupancy-preserving
  checkerboard shuffle + original-stretch prevalence-preserving stricter
  sensitivity (latter overconditions possible train transport).
No true rail cargo/seed arrival, wind exposure or heritable adaptation.
"""
from __future__ import annotations
import argparse,json,random
from pathlib import Path
from collections import Counter,defaultdict
from statistics import mean
from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    original_35_station_sites,real_site_pairs,
    conditional_source_gap,MIN_SAME_AND_DIFF_PAIRS_EACH_STRATUM
)
from analyze_Czech_2024_39station_functional_identity_v143 import jaccard
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import normalize_taxon
from analyze_Czech_2024_rail_station_area_aliens_v108 import EXPECTED_S1_SHA
from itertools import combinations
import hashlib,zipfile

NCHAINS=2
DRAWS_PER_CHAIN=120
BURNIN=15000
THIN=2800
SEED=20261115

def read_source_status(s1,expected=1053):
    raw=Path(s1).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=EXPECTED_S1_SHA or not zipfile.is_zipfile(s1):
        raise ValueError("Original Czech Preslia S1 publisher source hash mismatch")
    from openpyxl import load_workbook
    wb=load_workbook(s1,read_only=True,data_only=True)
    try:rows=list(next(iter(wb.worksheets)).values)
    finally:wb.close()
    d={}
    for row in rows[7:]:
        if len(row)<2 or row[1] is None or not str(row[1]).strip():continue
        name=normalize_taxon(row[1])
        if name in d:raise ValueError("Duplicated original plant normalized taxonomic name")
        cat=str(row[0] or "").casefold()
        tags={x.strip() for x in cat.split(",") if x.strip()}
        # Author categories 'nat, arch' = naturalized archaeophyte, NOT native;
        # 'inv, neo' = invasive neophyte. Untagged may carry red-list codes.
        alien_tags=tags & {"arch","neo"}
        if len(alien_tags)>1:raise ValueError("Original source one plant cannot have arch and neo simultaneously")
        d[name]=next(iter(alien_tags)) if alien_tags else None
    if len(d)!=expected:raise ValueError("Original Preslia 1053 distinct taxonomy source rows drifted")
    return d,{"original_S1_publisher_SHA256":EXPECTED_S1_SHA,
              "original_all_recent_and_historical_source_taxonomy_labels":len(d),
              "original_taxa_source_tag_arch":sum(v=="arch" for v in d.values()),
              "original_taxa_source_tag_neo":sum(v=="neo" for v in d.values()),
              "original_source_unlabelled_does_NOT_mean_proven_native":True}

def group_sets(sites,status):
    keys={"archaeophyte":{"arch"},"neophyte":{"neo"},"explicit_aliens_union":{"arch","neo"}}
    out={}
    for k,tags in keys.items():
        codes={name for name,code in status.items() if code in tags}
        sets=[s["all_original_taxa"] & codes for s in sites]
        out[k]={"source_site_sets":sets,
                "n_source_distinct_observed_class_taxa":len(set().union(*sets)),
                "n_source_site_taxon_incidence":sum(len(s) for s in sets),
                "n_actual_original_site_stations_with_class":sum(bool(s) for s in sets)}
    return out


def source_class_pairs(pair_rows,sets):
    # Shared source pair matching strata remain as in v15.4. If neither
    # site bears a member of the class, group-specific Jaccard undefined;
    # drop only double-empty. One-empty pair has Jaccard=0 legitimately.
    return [p for p in pair_rows if sets[p["i"]] or sets[p["j"]]]


def source_subgroup_gap(pairs,sets):
    same=[jaccard(sets[p["i"]],sets[p["j"]]) for p in pairs if p["same_stretch"]]
    other=[jaccard(sets[p["i"]],sets[p["j"]]) for p in pairs if not p["same_stretch"]]
    if not same or not other:return {"identified":False,"reason":"No original station pair support"}
    if any(x is None for x in same+other):raise ValueError("Zero/zero original plant group not filtered")
    return {"identified":True,"n_same_stretch_pairs":len(same),
            "n_other_stretch_pairs":len(other),
            "mean_jaccard_same":mean(same),
            "mean_jaccard_other":mean(other),
            "same_minus_other":mean(same)-mean(other)}


def source_single_group_swap(rows,sites,rng,proposals,within_stretch):
    strata=defaultdict(list)
    for i,s in enumerate(sites):
        strata[s["section"] if within_stretch else "all"].append(i)
    eligible=[v for v in strata.values() if len(v)>=2]
    nacc=0
    for _ in range(proposals):
        ind=rng.choice(eligible)
        i,j=rng.sample(ind,2)
        left=rows[i]-rows[j]
        right=rows[j]-rows[i]
        if not left or not right:continue
        a=rng.choice(tuple(left));b=rng.choice(tuple(right))
        rows[i].remove(a);rows[i].add(b)
        rows[j].remove(b);rows[j].add(a)
        nacc+=1
    return nacc


def conditional_group_null(sites,pairs,sets,mode="geography_and_function",
                           within_stretch=False,nchains=NCHAINS,
                           draws=DRAWS_PER_CHAIN,burnin=BURNIN,thin=THIN,seed=SEED):
    original=conditional_source_gap(pairs,sets=sets,mode=mode)
    if not original["identified"]:
        return {"identified":False,"source_support":original["source_support"],
                "reason":original["reason"]}
    obs=original["weighted_same_minus_different_jaccard"]
    null=[];diag=[]
    for chain in range(nchains):
        x=[set(s) for s in sets]
        initial_rich=[len(s) for s in x]
        frequency=Counter(y for s in x for y in s)
        stretch_frequency=Counter((sites[i]["section"],y) for i,s in enumerate(x) for y in s)
        rng=random.Random(seed+1009*chain+(77 if within_stretch else 0))
        accepted=source_single_group_swap(x,sites,rng,burnin,within_stretch)
        for _ in range(draws):
            accepted+=source_single_group_swap(x,sites,rng,thin,within_stretch)
            if [len(s) for s in x]!=initial_rich:
                raise ValueError("Source station × explicit alien subgroup richness changed")
            if Counter(y for s in x for y in s)!=frequency:
                raise ValueError("Original exact source taxon total frequency not held")
            if within_stretch and Counter(
                (sites[i]["section"],y) for i,s in enumerate(x) for y in s)!=stretch_frequency:
                raise ValueError("Original source within-stretch taxon frequency not preserved")
            stat=conditional_source_gap(pairs,sets=x,mode=mode)
            if not stat["identified"] or stat["source_support"]!=original["source_support"]:
                raise ValueError("Matched source opportunity changed under fixed richness swaps")
            null.append(stat["weighted_same_minus_different_jaccard"])
        diag.append({"chain":chain+1,"accepted":accepted,
                    "proposed":burnin+draws*thin})
    z=sorted(null)
    return {"identified":True,"original_jaccard_gap":obs,
            "null_expected_gap":mean(null),"observed_minus_null":obs-mean(null),
            "chain_percentile_2_5":z[round(.025*(len(z)-1))],
            "chain_percentile_97_5":z[round(.975*(len(z)-1))],
            "n_correlated_MCMC_draws":len(z),
            "station_subgroup_richness_and_taxon_total_frequency_preserved":True,
            "within_original_source_stretch_taxon_frequency_also_fixed":within_stretch,
            "strict_null_may_overcondition_real_corridor_effect":True,
            "not_independent_field_pair_p_values":True,
            "source_chain_diagnostics":diag}

def run(s1,s2,disp,life,sr70):
    sites,guild,source=original_35_station_sites(s1,s2,disp,life,sr70)
    allpairs=real_site_pairs(sites,guild)
    labels,labelprov=read_source_status(s1)
    all_groups=group_sets(sites,labels)
    group_results={}
    for label,data in all_groups.items():
        sets=data["source_site_sets"]
        pairs=source_class_pairs(allpairs,sets)
        valid={
            "source_class_original_taxa_in_35_locations":data["n_source_distinct_observed_class_taxa"],
            "original_subgroup_site_taxon_incidence":data["n_source_site_taxon_incidence"],
            "n_sites_with_at_least_one_subgroup_taxon":data["n_actual_original_site_stations_with_class"],
            "n_actual_pairwise_opportunities_not_double_zero":len(pairs),
            "raw":source_subgroup_gap(pairs,sets),
            "conditioning":{mode:conditional_source_gap(pairs,sets=sets,mode=mode)
                for mode in ("geography_only","geography_and_function",
                             "geography_function_and_area")},
        }
        matched=valid["conditioning"]["geography_and_function"]
        if matched["identified"] and valid["n_sites_with_at_least_one_subgroup_taxon"]>=20:
            valid["source_fixed_margins_null_global"]=conditional_group_null(
                sites,pairs,sets,within_stretch=False)
            valid["source_fixed_margins_null_by_original_stretch"]=conditional_group_null(
                sites,pairs,sets,within_stretch=True)
        else:
            valid["source_fixed_margins_null_global"]={"identified":False}
            valid["source_fixed_margins_null_by_original_stretch"]={"identified":False}
        group_results[label]=valid
    return {
      "version":"v15.7_REAL_Czech35_original_arch_vs_neo_station_identity_geographic",
      "original_study_source_DOI":"10.23855/preslia.2024.247",
      "original_source_35station_Pladias_SR70_and_preslia_provenance":source,
      "original_publisher_source_explicit_taxonomic_status":labelprov,
      "actual_source_original_station_class_specific_taxon_identity":group_results,
      "strict_limitations":[
        "Source S1 explicit arch/neo labels only; unlabelled original taxa may be native or other status and never imputed as confirmed native.",
        "Station point source is official-derived 2026 SR70 and NOT author 2020-22 real field plot coordinates.",
        "Original author 35 stations / five stretches / one country; 595 potential dependent station pairs not independent field studies.",
        "Explicit archaeophyte/neophyte groups are time-of-arrival categories, not a direct seed carriage observation or biological train vector classification.",
        "Class source Jaccard is taxonomic sharing, not abundance, plant recruitment, annual species rate or proven gene flow.",
        "Species frequency fixed-margin null may preserve actual rail effect if it is organized by railway stretches; all stochastic draws correlated.",
        "Original source station area/microhabitat management differs, and no actual train traffic, freight OD, mowing schedule, wind or fitness measured.",
      ],
      "never_pool_with_Slovenia_or_claim_movement_causality":True
    }

def main():
    p=argparse.ArgumentParser()
    for n in ("s1","s2","pladias-diaspore","pladias-lifeform","sr70","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.pladias_diaspore,a.pladias_lifeform,a.sr70)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"version":result["version"],
       "source":result["original_publisher_source_explicit_taxonomic_status"],
       "groups":result["actual_source_original_station_class_specific_taxon_identity"]},ensure_ascii=False))

if __name__=="__main__":main()

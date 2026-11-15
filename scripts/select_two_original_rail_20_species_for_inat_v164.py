#!/usr/bin/env python3
"""v16.4 preselect 20 real plant species frequent in two ORIGINAL author railway
floras, before looking at modern iNaturalist photo location results.

These are candidate rail-RECORDED species, NOT demonstrated railway-preferring
species (no matched non-rail species-site inventory in these original studies).
Slovenia: Šabić et al 2026 author 89 real 100m railway vegetation sites,
authentic original source supplement SHA pinned in v15.0 and v16.2.
Czech: Kutlvašr et al 2024 author 39 railway whole stations with S1 REC
actual source observations (not historical), publisher S1/S2 SHA checked.
Use site-occupancy rank, not taxon abundance or photo frequency.
Alternate source study ranks: 10/10 initial country candidates,
skip duplicate scientific binomial names, fill until exactly 20 unique.
Exclude infraspecific, aggregate and genus taxa from this *first exact
species* feasibility, explicitly record this preselection limitation.

The 30-m study is NOT a global exhaustive count at this stage: a complete
photo observation census is performed in fixed original SI/CZ bounding boxes
using real 2024/2023 OSM railway polylines. Global iNaturalist totals only
serve as opportunity denominators, NEVER global 30-m railroad counts.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from collections import Counter
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
from analyze_Czech_2024_39station_functional_identity_v143 import read_original_czech_publisher

DOIS={"SI":"10.1007/s44473-026-00155-z","CZ":"10.23855/preslia.2024.247"}
SITE_COUNTS={"SI":89,"CZ":39}
N_CANDIDATE=20
SPECIES_RE=re.compile(r"^[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+$")

def rank_observed_source_floras(site_sets):
    source=Counter()
    for site in site_sets:
        for raw in site:
            name=str(raw or "").strip()
            if SPECIES_RE.fullmatch(name):
                source[name]+=1
    return sorted(source.items(),key=lambda x:(-x[1],x[0]))

def source_real_floras(si_original,si_audit,cz_s1,cz_s2):
    orig=extract_original_xlsx(si_original,si_audit)
    si_sites,_,si_meta=source_sites_and_guilds(orig)
    cz_floras,cz_meta,cz_source=read_original_czech_publisher(cz_s1,cz_s2)
    if (len(si_sites)!=SITE_COUNTS["SI"] or len(cz_floras)!=SITE_COUNTS["CZ"] or
        si_meta["original_sampled_locations"]!=SITE_COUNTS["SI"] or
        cz_source["original_source_Czech_railway_Preslia_DOI"]!=DOIS["CZ"]):
        raise ValueError("Original source authors' site count or DOI changed")
    if si_meta["original_authors_field_study_DOI"]!=DOIS["SI"]:
        raise ValueError("Slovenia source original DOI mismatch")
    return {
      "SI":[set(site["all_original_taxa"]) for site in si_sites],
      "CZ":[set(x) for x in cz_floras.values()],
    },{
      "SI":{"source_DOI":DOIS["SI"],"original_source_sites":len(si_sites),
            "real_author_original_supplement_SHA":si_meta["original_source_supplement_SHA256"],
            "original_field_site_grain":"89 original 100m rail locations"},
      "CZ":{"source_DOI":DOIS["CZ"],"original_source_sites":len(cz_floras),
            "real_author_original_publisher_S1_SHA":cz_source["source_S1_SHA256"],
            "real_author_original_publisher_S2_SHA":cz_source["source_S2_corrected_erratum_SHA256"],
            "original_field_site_grain":"39 original unequal whole railway stations"},
    }

def top_candidates(floras,n=N_CANDIDATE):
    if n<2 or n%2:
        raise ValueError("Need an even number of source-selected species")
    ranked={country:rank_observed_source_floras(floras[country])
            for country in ("SI","CZ")}
    if any(len(ranked[country])<n for country in ("SI","CZ")):
        raise ValueError("Not enough original exact binomial species for no-guessed source 20")
    counts={country:dict(ranked[country]) for country in ("SI","CZ")}
    chosen=[];seen=set();idx={"SI":0,"CZ":0}
    # First 10 from each author study, distinct species if overlapping names.
    # If cross-country duplicates reduce 20, advance the same fixed
    # alternating sources until target n is reached, without seeing photo data.
    while len(chosen)<n:
        country="SI" if len(chosen)%2==0 else "CZ"
        while idx[country]<len(ranked[country]) and ranked[country][idx[country]][0] in seen:
            idx[country]+=1
        if idx[country]>=len(ranked[country]):
            raise ValueError("Original source taxonomic names overlap too heavily to fill top20")
        name,site_n=ranked[country][idx[country]]
        idx[country]+=1
        seen.add(name)
        chosen.append({
          "original_source_rank":len(chosen)+1,
          "original_source_selected_from":country,
          "scientific_name_original_exact_binomial":name,
          "original_rail_site_occupancy_n_in_source_country":site_n,
          "original_rail_site_occupancy_fraction":site_n/len(floras[country]),
          "other_original_country_same_source_site_count":counts["CZ" if country=="SI" else "SI"].get(name,0),
          "selection_not_based_on_iNaturalist_photos":True,
          "NOT_railway_affinity_vs_nonnear_railway":True,
        })
    if len(chosen)!=n or len({x["scientific_name_original_exact_binomial"] for x in chosen})!=n:
        raise ValueError("Original 20 source candidate distinct scientific species requirement failed")
    return chosen,{c:{"original_distinct_binomial_species":len(ranked[c]),
                     "source_ranked_binomials":len(ranked[c]),
                     "total_original_source_rail_sites":len(floras[c])} for c in ("SI","CZ")}

def run(si_original,si_audit,cz_s1,cz_s2):
    flora,source=source_real_floras(si_original,si_audit,cz_s1,cz_s2)
    chosen,coverage=top_candidates(flora)
    return {
      "complete":True,"version":"v16.4_AUTHOR_ORIGINAL_REAL_RAIL_FLORA_20_SPECIES_FIRST_PHOTO_FEASIBILITY",
      "selection_before_iNaturalist_2026_11_15_photo_results":True,
      "source_original_study_provenance":source,
      "source_original_research_ranking_and_taxon_eligibility":coverage,
      "selection":"alternate SI/CZ original observed station-occupancy ranks, exact unranked-photo binomials only; skip duplicates; 20 unique species",
      "n_distinct_original_species":len(chosen),
      "original_selected_20_species":chosen,
      "not_an_estimate_of_railway_affinity_vs_matched_nonnear_railway_sites":True,
      "not_original_phylogenetic_adaptation_or_train_carriage_evidence":True
    }

def main():
    p=argparse.ArgumentParser()
    for k in ("si-xlsx","si-audit","cz-s1","cz-s2","output"):
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    r=run(a.si_xlsx,a.si_audit,a.cz_s1,a.cz_s2)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":r["version"],
        "n":r["n_distinct_original_species"],
        "original_species":[
          {"rank":x["original_source_rank"],"country":x["original_source_selected_from"],
           "taxon":x["scientific_name_original_exact_binomial"],
           "original_rail_site_occupancy_n":x["original_rail_site_occupancy_n_in_source_country"]}
          for x in r["original_selected_20_species"]]
    },ensure_ascii=False))

if __name__=="__main__":main()

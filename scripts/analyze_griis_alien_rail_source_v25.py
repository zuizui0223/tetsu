#!/usr/bin/env python3
"""v2.5: country-specific introduced-plant species keys vs mapped railway graph.

TRIPLE REAL-SOURCE JOIN:
- GBIF full geographic-cell speciesKey indexed record counts (v2.0)
- Overture mapped railway physical connector graph (v2.1)
- GBIF indexed GRIIS country-introduced plant backbone keys (v2.4)

Unknown introduction status is NEVER called native. Taiwan unassessed.
This is a noncausal observation-bias-sensitive MACRO feasibility/association
pilot, NOT evidence that trains transport seeds or homogenize vegetation.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from datetime import datetime,timezone
from itertools import combinations
from math import isfinite, log
from pathlib import Path
from statistics import median

from analyze_macro_rail_assemblages_v22 import (
    COUNTRIES,PERIODS,source_species,incidence_probabilities,
    expected_rarefied_overlap,jaccard,geographic_distance_km,matched_contrasts
)

ALIEN_RAREFACTION_RECORDS=20
MIN_ALIEN_RECORDS=20
MIN_ALIEN_SPECIES=5
MIN_CHECKLIST_MATCHED_PLANT_SPECIES=100


def valid_griis(aliases):
    if not aliases.get("complete") or aliases.get("errors"):
        raise ValueError("Incomplete GRIIS source-country list")
    expected=set(COUNTRIES.values())-{"TW"}
    if set(aliases.get("countries",{}))!=expected:
        raise ValueError("Griis country frame must include nine exact countries")
    out={}
    for country,row in aliases["countries"].items():
        if row.get("country")!=country or not row.get("full_source_pages_verified"):
            raise ValueError("Country source provenance invalid")
        raw=row.get("gbif_backbone_plant_species_key_set")
        if not isinstance(raw,list) or any(not isinstance(x,str) or not x for x in raw):
            raise ValueError("Backbone source key set invalid")
        if len(set(raw))!=len(raw):
            raise ValueError("Duplicate GRIIS accepted species key")
        if len(raw)!=row.get("n_directly_matched_plant_backbone_species_keys"):
            raise ValueError("Backbone source key count mismatch")
        out[country]=set(raw) if len(raw)>=MIN_CHECKLIST_MATCHED_PLANT_SPECIES else None
    return out


def alien_source_species(source,alien_key_set):
    counts=source_species(source)
    if counts is None or alien_key_set is None:
        return None
    intro={key:n for key,n in counts.items() if key in alien_key_set}
    if len(intro)<MIN_ALIEN_SPECIES or sum(intro.values())<MIN_ALIEN_RECORDS:
        return None
    return intro


def cell_observations(city,aliases,period):
    result={}
    all_qualified_union=set()
    for name,cell in city["cells"].items():
        source=cell["by_period"][period]
        all_plants=source_species(source)
        if all_plants:
            all_qualified_union.update(set(all_plants)&aliases)
        intro=alien_source_species(source,aliases)
        if intro:
            result[name]={
                "counts":intro,
                "probs":incidence_probabilities(intro,k=ALIEN_RAREFACTION_RECORDS),
                "bbox":cell["bbox"],
                "n_aliens":len(intro),
                "n_alien_keyed_records":sum(intro.values()),
                "n_total_source_records":source["indexed_occurrence_count"]
            }
    return result,all_qualified_union


def city_pair_analysis(plant,rail,aliens,period):
    cell_data,source_union=cell_observations(plant,aliens,period)
    rail_map=rail["cell_component_cover"]
    if set(rail_map)!=set(plant["cells"]):
        raise ValueError("Rail cells and GBIF city cells differ")
    pair_lookup={}
    for row in rail["cell_pair_connection_descriptors"]:
        key=tuple(sorted((row["a"],row["b"])))
        if key in pair_lookup or key[0]==key[1]:
            raise ValueError("Duplicate graph rail pair")
        pair_lookup[key]=row["within_clipped_graph_connected"]
    if len(pair_lookup)!=36:
        raise ValueError("Missing mapped rail-cell source pairs")
    rows=[]
    for a,b in combinations(sorted(cell_data),2):
        x,y=cell_data[a],cell_data[b]
        rail_connected=pair_lookup.get((a,b))
        if not isinstance(rail_connected,bool):
            raise ValueError("Unverifiable rail connection status")
        obs=jaccard(x["counts"],y["counts"])
        rare=expected_rarefied_overlap(x["probs"],y["probs"])
        if rare is None or not isfinite(rare):
            continue
        rows.append({
            "a":a,"b":b,
            "mapped_rail_connected_within_window":rail_connected,
            "both_cells_have_source_rail":bool(rail_map[a]["n_source_components"] and rail_map[b]["n_source_components"]),
            "geographic_distance_km":geographic_distance_km(x["bbox"],y["bbox"]),
            "log_occurrence_effort_ratio":abs(log(x["n_total_source_records"]/y["n_total_source_records"])),
            "log_geomean_records":(log(x["n_total_source_records"])+log(y["n_total_source_records"]))/2,
            "observed_introduced_species_key_jaccard":obs,
            "ratio_expected_rarefied_100_record_intersection_to_union":rare
        })
    result={
        "n_qualified_all_plant_cell_sources_with_at_least_one_griis_alien":sum(bool(
            source_species(cell["by_period"][period]) and set(
                source_species(cell["by_period"][period]) or {})&aliens
            ) for cell in plant["cells"].values()),
        "n_source_eligible_introduced_plant_cells":len(cell_data),
        "n_griis_introduced_species_keys_observed_in_qualified_total_plant_cells":len(source_union),
        "n_independent_cell_pair_source_rows":len(rows),
        "n_connected_cell_pair_sources":sum(r["mapped_rail_connected_within_window"] for r in rows),
        "all_pairs_comparison":matched_contrasts(rows,"all"),
        "strict_both_rail_comparison":matched_contrasts(rows,"both")
    }
    return result


def summarize_regions(result,period,group):
    values=[r["periods"][period][group]["median_matched_similarity_delta"]
        for r in result.values()
        if r.get("country_list_adequate") and
           r["periods"][period][group]["source_comparison_adequate"]]
    return {
        "n_independent_geographic_windows_with_matched_rail_contrasts":len(values),
        "n_positive_windows":sum(v>0 for v in values),
        "n_negative_windows":sum(v<0 for v in values),
        "median_city_level_descriptive_delta":round(median(values),6) if values else None,
        "five_city_minimum_screen":len(values)>=5,
        "no_causal_or_statistical_significance_test":True
    }


def analyze(plant,rail,griis):
    if not plant.get("source_complete") or not rail.get("complete"):
        raise ValueError("Incomplete true-source geography or railway graph")
    if set(plant.get("regions",{}))!=set(COUNTRIES) or set(rail.get("regions",{}))!=set(COUNTRIES):
        raise ValueError("Source region frame mismatch")
    aliases=valid_griis(griis)
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "study":"Exploratory 9-country GRIIS-introduced only speciesKey × observed records × mapped Overture rail graph",
        "unassessed_country":{"TW":"No verified comparable GRIIS country checklist, source status unknown"},
        "country_lists":{},"regions":{},"independent_city_direction":{},
        "not_train_seed_transport_effect":True,
        "no_country_nonlisted_species_treated_as_native":True,
    }
    for name,country in COUNTRIES.items():
        if plant["regions"][name].get("country")!=country or rail["regions"][name].get("country")!=country:
            raise ValueError("Country mismatch in actual sources")
        if country=="TW":
            out["regions"][name]={"country":country,"alien_status":"unassessed","periods":{}}
            continue
        alias=aliases[country]
        row=griis["countries"][country]
        out["country_lists"][country]={
            "source_dataset":row["source_dataset_reference"],
            "n_griis_plant_backbone_species_keys":row["n_directly_matched_plant_backbone_species_keys"],
            "n_name_usages_unmatched":row["n_unmatched_plant_name_usages"],
            "usable_checklist":alias is not None
        }
        result={"country":country,"country_list_adequate":alias is not None,"periods":{}}
        if alias is not None:
            for period in PERIODS:
                ans=city_pair_analysis(plant["regions"][name],rail["regions"][name],alias,period)
                result["periods"][period]=ans
                print(json.dumps({
                    "city":name,"period":period,
                    "observed_alien_taxon_union":ans["n_griis_introduced_species_keys_observed_in_qualified_total_plant_cells"],
                    "alien_cells_eligible":ans["n_source_eligible_introduced_plant_cells"],
                    "source_pairs":ans["n_independent_cell_pair_source_rows"],
                    "rail_connected_pairs":ans["n_connected_cell_pair_sources"],
                    "matched_pairs":ans["all_pairs_comparison"]["n_distance_and_record_effort_matched_pairs"],
                    "descriptive_rarefied20_delta":ans["all_pairs_comparison"]["median_matched_similarity_delta"]
                },ensure_ascii=False),flush=True)
        else:
            result["alien_status"]="source checklisted but too small (<100 Plantae accepted keys)"
        out["regions"][name]=result
    for period in PERIODS:
        out["independent_city_direction"][period]={
            group:summarize_regions(out["regions"],period,group)
            for group in ("all_pairs_comparison","strict_both_rail_comparison")
        }
    out["scientific_limits"]=[
        "GRIIS country checklist membership is evidence for recorded introduction status in that country, not establishment at a mapped railway.",
        "The source checklists are updated in different years and can be incomplete. Non-membership cannot classify native species.",
        "Species keys are only observed indexed GBIF records, not full floristic or taxon-by-survey plot absence data.",
        "Twenty-record rarefaction makes record count more comparable but does not standardize observer/data-provider geography, survey routes or effort.",
        "The railway graph is mapped physical connectivity, not verified historical active freight flow, and is cut at each city window.",
        "Current study lacks comparable road graph, climate, port activity, urban land use and validated station/trackside plant observations.",
        "Country samples are 1 city each and do not establish biogeographic causal generalization even when multiple city comparisons share a sign.",
        "No p-values are computed because cell pairs sharing sites are not independent; estimation is descriptive and limited by positivity."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--gbif-source",required=True)
    p.add_argument("--rail-source",required=True)
    p.add_argument("--griis-source",required=True)
    p.add_argument("--output",default="results/audit/alien_rail_source_cross9country_v25.json")
    a=p.parse_args()
    data=[json.loads(Path(x).read_text(encoding="utf-8")) for x in (
        a.gbif_source,a.rail_source,a.griis_source)]
    result=analyze(*data)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"n_countries_with_source_list":len(result["country_lists"]),
        "direction":result["independent_city_direction"],
        "noncausal":result["not_train_seed_transport_effect"]},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

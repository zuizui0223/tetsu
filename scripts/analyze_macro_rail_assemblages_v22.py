#!/usr/bin/env python3
"""v2.2: real-source GBIF observed-species vs mapped rail connectivity diagnostic.

INPUTS: immutable GitHub Actions aggregate JSONs, 10 country city windows.
Ecological scope: 2010–17 / 2018–25 observed speciesKey composition only.
NO country-specific alien-status, active freight services, full road topology,
or observer/source-standardized biological surveys, hence no causal inference.

Replicates are COUNTRIES/WINDOWS, never within-city site pairs.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from math import asin, cos, expm1, isfinite, lgamma, log, radians, sin, sqrt
from pathlib import Path
from statistics import median

PERIODS = ("2010_2017","2018_2025")
FIXED_RAREFACTION_RECORDS = 100
MIN_RECORDS = 100
MIN_DATASETS = 3
MIN_PAIRED_CATEGORIES = 3
COUNTRIES = {
    "jp_sendai":"JP","de_leipzig":"DE","kr_daejeon":"KR","tw_hsinchu":"TW",
    "us_madison":"US","ca_ottawa":"CA","gb_york":"GB","fr_dijon":"FR",
    "au_canberra":"AU","nz_hamilton":"NZ",
}


def source_species(s):
    """Count only positively recorded canonical species-key categories."""
    if not isinstance(s,dict) or not s.get("full_facet_pagination_completed"):
        return None
    if not s.get("source_eligible_for_recording_effort_screen"):
        return None
    if s.get("n_source_datasets_at_least",0)<MIN_DATASETS:
        return None
    if s.get("indexed_occurrence_count",0)<MIN_RECORDS:
        return None
    raw=s.get("species_key_occurrence_counts")
    if not isinstance(raw,dict):
        return None
    counts={}
    for k,v in raw.items():
        if not isinstance(k,str) or not k or isinstance(v,bool) or not isinstance(v,int) or v<0:
            raise ValueError("Malformed source GBIF species facet value")
        if v>0:
            counts[k]=v
    if len(counts)!=s.get("n_species_key_categories"):
        raise ValueError("Taxon source category count mismatch")
    if sum(counts.values())<FIXED_RAREFACTION_RECORDS:
        return None
    if sum(counts.values())>s["indexed_occurrence_count"]:
        raise ValueError("Keyed source count exceeds source records")
    return counts


def incidence_probabilities(counts, k=FIXED_RAREFACTION_RECORDS):
    """P(at least one species record in simple random k-without-replacement draw).

    This is a model of rarefied INDEXED RECORDS, not individual organisms or
    independent samples. Does not correct observer, dataset or spatial biases.
    """
    n=sum(counts.values())
    if not isinstance(k,int) or k<1 or n<k:
        raise ValueError("Cannot rarefy inadequate source material")
    result={}
    for taxon,c in counts.items():
        if c<=0 or c>n:
            raise ValueError("Bad species-record source count")
        if n-c<k:
            p=1.
        else:
            log_zero=(lgamma(n-c+1)-lgamma(n-c-k+1)
                      -lgamma(n+1)+lgamma(n-k+1))
            p=-expm1(min(0.,log_zero))
        result[taxon]=max(0.,min(1.,p))
    return result


def expected_rarefied_overlap(p,q):
    """Ratio of expected intersection and union, NOT E[Jaccard]."""
    keys=set(p)|set(q)
    intersection=sum(p.get(k,0.)*q.get(k,0.) for k in keys)
    union=sum(p.get(k,0.)+q.get(k,0.)-p.get(k,0.)*q.get(k,0.) for k in keys)
    return intersection/union if union>0 else None


def jaccard(a,b):
    x=set(a);y=set(b)
    union=x|y
    return len(x&y)/len(union) if union else None


def geographic_distance_km(box1,box2):
    if len(box1)!=4 or len(box2)!=4:
        raise ValueError("Malformed study cells")
    x1=(box1[0]+box1[2])/2; y1=(box1[1]+box1[3])/2
    x2=(box2[0]+box2[2])/2; y2=(box2[1]+box2[3])/2
    lat1,lat2=radians(y1),radians(y2)
    dlat=lat2-lat1
    dlon=radians(x2-x1)
    a=sin(dlat/2)**2+cos(lat1)*cos(lat2)*sin(dlon/2)**2
    return 6371.0088*2*asin(min(1,sqrt(max(0,a))))


def calculate_pairs(bio,graph,period):
    """Pre-outcome masks: recording opportunity and graph data availability."""
    if len(bio["cells"])!=9 or len(graph["cell_pair_connection_descriptors"])!=36:
        raise ValueError("Expected fixed 9 cells / 36 source rail-cell pairs")
    cell_species={}
    for cid,cell in bio["cells"].items():
        if cid not in graph["cell_component_cover"]:
            raise ValueError("Rail and GBIF cell IDs disagree")
        if period not in cell["by_period"]:
            raise ValueError("Missing specified period")
        counts=source_species(cell["by_period"][period])
        if counts is not None:
            cell_species[cid]={
                "counts":counts,
                "prob":incidence_probabilities(counts),
                "bbox":cell["bbox"],
                "source_n":cell["by_period"][period]["indexed_occurrence_count"],
                "keyed_source_n":sum(counts.values()),
                "rail_present":graph["cell_component_cover"][cid]["n_source_components"]>0,
            }
    rail_pairs={}
    for r in graph["cell_pair_connection_descriptors"]:
        key=tuple(sorted([r["a"],r["b"]]))
        if key in rail_pairs or key[0]==key[1]:
            raise ValueError("Duplicate/invalid mapped rail pair")
        rail_pairs[key]=r["within_clipped_graph_connected"]
    if len(rail_pairs)!=36:
        raise ValueError("Missing graph pair")
    rows=[]
    for a,b in combinations(sorted(cell_species),2):
        u,v=cell_species[a],cell_species[b]
        conn=rail_pairs.get((a,b))
        if not isinstance(conn,bool):
            raise ValueError("No physical rail-pair source status")
        raw=jaccard(u["counts"],v["counts"])
        r100=expected_rarefied_overlap(u["prob"],v["prob"])
        if raw is None or r100 is None or not isfinite(r100):
            continue
        rows.append({
            "a":a,"b":b,
            "mapped_rail_connected_within_window":conn,
            "both_cells_have_source_rail":u["rail_present"] and v["rail_present"],
            "geographic_distance_km":geographic_distance_km(u["bbox"],v["bbox"]),
            "log_occurrence_effort_ratio":abs(log(u["source_n"]/v["source_n"])),
            "log_geomean_records":(log(u["source_n"])+log(v["source_n"]))/2,
            "observed_species_key_jaccard":raw,
            "ratio_expected_rarefied_100_record_intersection_to_union":r100,
        })
    return rows,len(cell_species)


def matched_contrasts(rows,subset="all",metric="ratio_expected_rarefied_100_record_intersection_to_union"):
    """Source-matched descriptive pairs; NO hypothesis p-values.

    Matched on distance and record-count geometry, never on species overlap.
    Pairs sharing cells remain dependent; must aggregate at city block level.
    """
    source=rows if subset=="all" else [r for r in rows if r["both_cells_have_source_rail"]]
    exposed=[r for r in source if r["mapped_rail_connected_within_window"]]
    control=[r for r in source if not r["mapped_rail_connected_within_window"]]
    sorted_exposed=sorted(exposed,key=lambda r:(r["geographic_distance_km"],r["a"],r["b"]))
    used=set();matches=[]
    for positive in sorted_exposed:
        scored=[]
        for j,candidate in enumerate(control):
            if j in used: continue
            dd=abs(positive["geographic_distance_km"]-candidate["geographic_distance_km"])
            max_dd=max(5.,.33*max(positive["geographic_distance_km"],candidate["geographic_distance_km"]))
            if dd>max_dd: continue
            dg=abs(positive["log_geomean_records"]-candidate["log_geomean_records"])
            if dg>1.:continue
            de=abs(positive["log_occurrence_effort_ratio"]-candidate["log_occurrence_effort_ratio"])
            if de>1.:continue
            score=dd/max_dd+dg+de
            scored.append((score,j))
        if scored:
            _,j=min(scored)
            used.add(j)
            negative=control[j]
            matches.append(positive[metric]-negative[metric])
    adequate=(len(matches)>=MIN_PAIRED_CATEGORIES
              and len(exposed)>=MIN_PAIRED_CATEGORIES
              and len(control)>=MIN_PAIRED_CATEGORIES)
    return {
        "n_pairs_with_eligible_species_source":len(source),
        "n_rail_connected_pairs":len(exposed),
        "n_unconnected_pairs":len(control),
        "n_distance_and_record_effort_matched_pairs":len(matches),
        "source_comparison_adequate":adequate,
        "median_matched_similarity_delta":round(median(matches),6) if adequate else None,
        "median_unadjusted_source_delta":round(median(x[metric] for x in exposed)-median(x[metric] for x in control),6)
             if exposed and control else None,
        "not_a_railway_biological_effect":True,
    }


def analyze(plant,rail):
    if not plant.get("source_complete") or not rail.get("complete"):
        raise ValueError("Incomplete source; cannot report cross-source comparisons")
    if set(plant.get("regions",{}))!=set(COUNTRIES) or set(rail.get("regions",{}))!=set(COUNTRIES):
        raise ValueError("Country frame missing or changed")
    out={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "study_kind":"Exploratory raw-plant-record speciesKey compositional diagnostic vs mapped rail physical connector network",
        "source_status":"geographically joined 10 fixed cities x 3x3 source cells",
        "national_status":"ten geographic source windows, NOT ten independently sampled national populations",
        "regions":{}, "country_equal_weight_direction":{},
        "source_significance_p_value":None,
        "causal_train_or_alien_plant_effect_unidentified":True,
    }
    for name,country in COUNTRIES.items():
        b=plant["regions"][name]
        g=rail["regions"][name]
        if b.get("country")!=country or g.get("country")!=country:
            raise ValueError("Country or geographic-cell source mismatch")
        summary={"country":country,"source_rail_components":g["n_distinct_components_in_source_window"],
                 "periods":{}}
        for period in PERIODS:
            pairs,n=calculate_pairs(b,g,period)
            result={
                "n_eligible_independent_cell_sources":n,
                "n_pairwise_observations_that_share_cells":len(pairs),
                "unadjusted":matched_contrasts(pairs,"all"),
                "both_cells_mapped_rail":matched_contrasts(pairs,"both")
            }
            summary["periods"][period]=result
            print(json.dumps({"region":name,"period":period,"cells":n,"pairs":len(pairs),
                "connected":result["unadjusted"]["n_rail_connected_pairs"],
                "unconnected":result["unadjusted"]["n_unconnected_pairs"],
                "matched":result["unadjusted"]["n_distance_and_record_effort_matched_pairs"],
                "delta_rarefied":result["unadjusted"]["median_matched_similarity_delta"],
                "both_rail_delta":result["both_cells_mapped_rail"]["median_matched_similarity_delta"]}),flush=True)
        out["regions"][name]=summary
    for period in PERIODS:
        valid=[(n,r["periods"][period]["unadjusted"]["median_matched_similarity_delta"]) for n,r in out["regions"].items()
            if r["periods"][period]["unadjusted"]["source_comparison_adequate"]]
        strict=[(n,r["periods"][period]["both_cells_mapped_rail"]["median_matched_similarity_delta"]) for n,r in out["regions"].items()
            if r["periods"][period]["both_cells_mapped_rail"]["source_comparison_adequate"]]
        def summarize(group):
            vals=[v for _,v in group]
            return {"n_informative_independent_city_blocks":len(vals),
                    "positive_blocks":sum(v>0 for v in vals),
                    "negative_blocks":sum(v<0 for v in vals),
                    "median_across_city_blocks":round(median(vals),6) if vals else None,
                    "minimum_blocks_for_interpretation":5,
                    "directional_result_supported":len(vals)>=5,
                    "not_survey_based_railway_causal_inference":True}
        out["country_equal_weight_direction"][period]={"all_cells":summarize(valid),
                    "both_cells_with_rail":summarize(strict)}
    out["interpretation_limits"]=[
        "The response is observed indexed GBIF species-key composition, not independently verified plant community occupancy or non-native species.",
        "100-record rarefaction addresses unequal RECORD count only, not observer effort, spatial reporting heterogeneity or dataset composition.",
        "2026 Overture connectors include mapped and potentially disused/nonoperating rail; do not claim freight flows.",
        "Physical connectivity is truncated at chosen city window, and lack of within-window path is not proof of geographical absence of rail connection.",
        "Country-specific alien status, matched road network connectivity, climate and ports not yet joined. Any positive direction is noncausal exploratory.",
        "City pair comparisons are dependent; no per-pair p-values, or fake precision from hundreds of non-independent pairs.",
        "Minimum 5 blocks is a design screen, not a significance threshold.",
    ]
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gbif-source",required=True)
    ap.add_argument("--rail-source",required=True)
    ap.add_argument("--output",default="results/audit/rail_gbif_observed_composition_joint_v22.json")
    args=ap.parse_args()
    p=json.loads(Path(args.gbif_source).read_text(encoding="utf-8"))
    q=json.loads(Path(args.rail_source).read_text(encoding="utf-8"))
    result=analyze(p,q)
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"n_regions":len(result["regions"]),"country_block_direction":result["country_equal_weight_direction"],
                      "no_causal_effect":result["causal_train_or_alien_plant_effect_unidentified"]},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

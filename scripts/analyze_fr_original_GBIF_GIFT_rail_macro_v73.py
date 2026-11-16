#!/usr/bin/env python3
"""v7.3: REAL FR Burgundy original GBIF event vs old indexed facet functional
RAIL/ROAD macro source diagnostic, on identical physical 50km cell pairs.

Source contains country-GRIIS strict COL-accepted species, independently
matched by legacy accepted GBIF species keys to GIFT 1.2.1 growth form.
Real 81-file 2018-25 original GBIF A5 source is NONEXHAUSTIVE. This is
not observed railway-carried seed transfer, stress, adaptation or evolution.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from itertools import combinations
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from analyze_macro_50km_rail_road_aliens_v34 import _routes, partial_rank_r
from analyze_macro_rail_assemblages_v22 import incidence_probabilities, expected_rarefied_overlap

BLOCK="fr_burgundy"
PERIOD="2018_2025"
GROUPS=("herbaceous","woody")
SOURCES=("original_event_all","original_event_precise5km","legacy_indexed_facet")
ANCHORS=("nearest_in_cell","nearest_largest_component_in_cell")
K=20
MIN_GROUP_TAXA=3
MIN_GROUP_RECORDS=20
MIN_CELLS=8
MIN_PAIRS=20

def growth_bridge(crosswalk,gift):
    if gift.get("complete") is not True or gift.get("source",{}).get("GIFT_meta_ID")!="1.2.1":
        raise ValueError("Original GIFT 1.2.1 source not verified")
    fr=crosswalk.get("countries",{}).get("FR")
    if not fr or not fr.get("old_GRIIS_to_COL_strict_accepted_map"):
        raise ValueError("Original six-country official FR GBIF-COL strict map missing")
    legacy=fr["old_GRIIS_to_COL_strict_accepted_map"]
    trait=gift.get("growth_form_by_source_gbif_species_key",{})
    groups=defaultdict(set)
    for old,col in legacy.items():
        if old in trait:
            if trait[old] not in GROUPS:raise ValueError("Unexpected GIFT life form")
            groups[col].add(trait[old])
    return {col:next(iter(g)) for col,g in groups.items() if len(g)==1}, \
           {col for col,g in groups.items() if len(g)>1},legacy

def grouped(counts,traits):
    out={g:{} for g in GROUPS}
    for key,count in counts.items():
        if type(count) is not int or count<=0:
            raise ValueError("Original source species counts invalid")
        g=traits.get(key)
        if g in out:
            out[g][key]=out[g].get(key,0)+count
    return out

def prepared_cells(event,crosswalk,gift,flora):
    if (event.get("source_type","").startswith("REAL interpreted original GBIF") is not True
        or event.get("frozen_geographical_block")!=BLOCK
        or event.get("source_period_analyzed")!=PERIOD
        or event.get("n_original_equal_area_50km_cells")!=16
        or event.get("source_independently_verified_same_all_Plantae_bbox_row_count") is not True
        or event.get("source_archive_NOT_complete_regional_original_events_until_A5_partition_gap_resolved") is not True
        or event.get("NO_RAIL_effect_or_train_dispersal_result") is not True):
        raise ValueError("Original 81-source experiment provenance/inference gates invalid")
    e=event["by_period"][PERIOD]["cells"]
    f=flora["regions"][BLOCK]["cells"]
    if set(e)!=set(f) or len(e)!=16:
        raise ValueError("Original and indexed FR 50km cell identity mismatch")
    traits,ambiguous,legacy=growth_bridge(crosswalk,gift)
    out={}
    for cid,ev in e.items():
        fp=f[cid]["periods"][PERIOD]
        if fp.get("full_facet_source_pagination_completed") is not True:
            raise ValueError("Original GBIF indexed species facet source incomplete")
        indexed=defaultdict(int)
        for old,n in fp["observed_alien_specieskey_counts"].items():
            col=legacy.get(old)
            if col:indexed[col]+=n
        out[cid]={
            "original_event_all":grouped(ev["COL_strict_country_GRIIS_speciesKey_event_counts"],traits),
            "original_event_precise5km":grouped(ev["COL_strict_country_GRIIS_speciesKey_known_precision_5km_counts"],traits),
            "legacy_indexed_facet":grouped(dict(indexed),traits)}
    return out,traits,ambiguous

def source_measure(counts):
    if len(counts)<MIN_GROUP_TAXA or sum(counts.values())<MIN_GROUP_RECORDS:
        return None
    return incidence_probabilities(counts,k=K)

def paired_rows(cells,flora,rail,road,anchor):
    rails=_routes(rail,BLOCK,anchor,"pair_routes")
    roads=_routes(road,BLOCK,anchor,"pair_paths")
    if set(rails)!=set(roads):raise ValueError("Physical RAIL and ROAD cell-pair identities differ")
    probs={cid:{s:{g:source_measure(cells[cid][s][g]) for g in GROUPS}
                for s in SOURCES} for cid in cells}
    pairs=[]
    for a,b in combinations(sorted(cells),2):
        r=rails[(a,b)]
        d=roads[(a,b)]
        rp=r.get("cell_centroid_snap_augmented_rail_path_km")
        dp=d.get("centroid_snap_augmented_road_path_km")
        geo=r.get("source_cell_centroid_geodesic_km")
        rr=r.get("rail_to_geodesic_ratio")
        dr=d.get("road_to_geodesic_ratio")
        vals=(rp,dp,geo,rr,dr)
        if not all(type(z) in (int,float) and math.isfinite(z) and z>0 for z in vals):continue
        if rr<.95 or dr<.95:continue
        fa=flora["regions"][BLOCK]["cells"][a]["periods"][PERIOD]
        fb=flora["regions"][BLOCK]["cells"][b]["periods"][PERIOD]
        na=fa["indexed_all_plantae_source_records"]
        nb=fb["indexed_all_plantae_source_records"]
        if not (na>0 and nb>0):raise ValueError("GBIF source record opportunity missing")
        rows={}
        for s in SOURCES:
            rows[s]={}
            for g in GROUPS:
                pa,pb=probs[a][s][g],probs[b][s][g]
                if pa is None or pb is None:
                    rows[s][g]=None
                    continue
                expected=expected_rarefied_overlap(pa,pb)
                aa,bb=cells[a][s][g],cells[b][s][g]
                raw=len(set(aa)&set(bb))/len(set(aa)|set(bb))
                rows[s][g]={"rarefied20":expected,"raw_jaccard":raw}
        if any(rows[s][g] is None or rows[s][g]["rarefied20"] is None
               for s in SOURCES for g in GROUPS):
            continue
        pairs.append({"a":a,"b":b,"relative_rail_path":math.log(rp/dp),
                      "covariates":[math.log(geo),math.log(dp),abs(math.log(na/nb)),
                                    (math.log(na)+math.log(nb))/2],
                      "outcomes":rows})
    return pairs

def analyze(event,crosswalk,gift,flora,rail,road):
    cells,traits,ambiguous=prepared_cells(event,crosswalk,gift,flora)
    result={"complete":True,"utc":datetime.now(timezone.utc).isoformat(),
            "study":"Original 81-file GBIF event vs legacy GBIF facet GRIIS-COL-GIFT growth form × mapped rail/road FR macro 50km",
            "one_region_only_FR_no_independent_country_replication":True,
            "period":PERIOD,"n_original_geographical_cells":len(cells),"rarefaction_indexed_event_k":K,
            "strict_direct_COL_growth_form_concepts":len(traits),
            "n_conflicting_source_COL_growth_form_concepts_excluded":len(ambiguous),
            "source_group_counts":{},"anchors":{}}
    for s in SOURCES:
        result["source_group_counts"][s]={}
        for g in GROUPS:
            count=sum(sum(cells[cid][s][g].values()) for cid in cells)
            species={k for cid in cells for k in cells[cid][s][g]}
            result["source_group_counts"][s][g]={"n_record_or_event_species_counts":count,
                                                  "n_distinct_COL_species_concepts":len(species)}
    anchor_raw={}
    for anchor in ANCHORS:
        pairs=paired_rows(cells,flora,rail,road,anchor)
        n_unique=len({c for row in pairs for c in (row["a"],row["b"])})
        stats={}
        for s in SOURCES:
            stats[s]={}
            for g in GROUPS:
                stats[s][g]={}
                for metric in ("rarefied20","raw_jaccard"):
                    x=[p["relative_rail_path"] for p in pairs]
                    y=[p["outcomes"][s][g][metric] for p in pairs]
                    z=[p["covariates"] for p in pairs]
                    rho=partial_rank_r(x,y,z) if n_unique>=MIN_CELLS and len(pairs)>=MIN_PAIRS else None
                    stats[s][g][metric]={"partial_rank_rho":rho}
        result["anchors"][anchor]={"n_same_source_cells_in_qualified_pairs":n_unique,
                                   "n_same_physical_cell_pairs_NOT_independent":len(pairs),
                                   "source_response_comparison":stats}
        anchor_raw[anchor]=[(p["a"],p["b"],p["relative_rail_path"]) for p in pairs]
    result["two_anchor_methods_independent_network_variation"]=anchor_raw[ANCHORS[0]]!=anchor_raw[ANCHORS[1]]
    result["limitations"]=[
        "ONE selected French 200km geographic macro block; 120 pairs SHARE only 16 cells, no pairwise p value or six-country generalization.",
        "2026 physical RAIL/ROAD networks are not real train journeys, dated active freight links, seed flux or historical rail openings.",
        "August 2026 A5 GBIF original interpreted Parquet 81-source candidate union is NOT a verified complete regional GBIF archive.",
        "Original event dedup does not guarantee observer/dataset independence, standardized plant absences or taxonomically complete flora.",
        "The full vs precise-5km source groups overlap; they are not independent external validation sets.",
        "Legacy GBIF indexed species facets differ in source vintage/taxonomic concept from experimental August 2026 original events.",
        "GIFT growth form does NOT measure grass cutting tolerance, slipstream, wind exposure, plasticity or genetic adaptation.",
        "No p-values are reported because original 16 cells generate dependent 120 pair comparisons."
    ]
    return result

def main():
    ap=argparse.ArgumentParser()
    for a in ("event","crosswalk","gift","flora","rail","road","output"):
        ap.add_argument("--"+a,required=True)
    v=ap.parse_args()
    load=lambda path:json.loads(Path(path).read_text(encoding="utf-8"))
    result=analyze(*[load(getattr(v,x)) for x in ("event","crosswalk","gift","flora","rail","road")])
    dest=Path(v.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,
        "source_groups":result["source_group_counts"],
        "anchors":result["anchors"],
        "no_seed_or_mowing_or_wind_causality":True},ensure_ascii=False),flush=True)
if __name__=="__main__":main()

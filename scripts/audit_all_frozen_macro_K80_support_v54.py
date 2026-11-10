#!/usr/bin/env python3
"""v5.4 feasibility audit across ALL TEN pre-frozen 50km source blocks.

No railway outcome or physical path is used here. Reuse one immutable
country-masked 8-country GRIIS-plant GBIF speciesKey source run v3.6:
10 x 16 = 160 mapped 2500km² cells, two fixed periods, one source
country each.

At recording depths 20,40,80,160,320 indexed introduced plant GBIF
RECORDS, count source-eligible geographic CELLS in BOTH periods.
No data-driven selection of country blocks or rail outcome signs.
This quantifies how many geographically independent macro regions
can feasibly test K80 *before* new biological route modeling.
"""
from __future__ import annotations
import argparse,json
from datetime import datetime,timezone
from pathlib import Path

from audit_gbif_equalarea_macroblocks_v29 import PILOT,BLOCKS
from analyze_6country_50km_network_plants_v44 import PERIODS

DEPTHS=(20,40,80,160,320)
MIN_INDEPENDENT_SOURCE_CELLS=8
MIN_MAX_PAIR_OPPORTUNITIES=20


def eligible_count_at_depth(row,k):
    if row.get("source_pre_outcome_both_period_recording_eligible") is not True:
        return False
    for p in PERIODS:
        d=row.get("periods",{}).get(p)
        if not isinstance(d,dict) or not d.get("full_facet_source_pagination_completed"):
            return False
        if d.get("alien_observed_minimum_5_species_and_20_records") is not True:
            return False
        counts=d.get("observed_alien_specieskey_counts")
        if not isinstance(counts,dict) or len(counts)<5 or any(
            not isinstance(v,int) or isinstance(v,bool) or v<=0 for v in counts.values()
        ):
            raise ValueError("Malformed original country-GRIIS indexed speciesKey source")
        n=sum(counts.values())
        if n!=d.get("n_indexed_GRIIS_introduced_plant_records"):
            raise ValueError("Original indexed GRIIS species record count inconsistent")
        if n<k:
            return False
    return True


def analyze(flora):
    if not flora.get("complete") or flora.get("source_errors"):
        raise ValueError("Original 10 source country regions incomplete")
    if set(flora.get("regions",{}))!=set(PILOT):
        raise ValueError("Pre-frozen 10 country/region source coverage changed")
    output={
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_source":"griis-50km-10block-8country-introduced-species-source-v36, run #42761140912",
        "pre_outcome_frozen_region_frame":list(PILOT),
        "depths_indexed_source_RECORDS":list(DEPTHS),
        "n_source_country_blocks":len(PILOT),
        "n_distinct_source_countries":len({BLOCKS[b][0] for b in PILOT}),
        "total_source_geographic_cells":len(PILOT)*16,
        "regions":{},"source_summary":{}
    }
    for block in PILOT:
        region=flora["regions"][block]
        if region.get("country")!=BLOCKS[block][0]:
            raise ValueError("Country GRIIS source identity changed")
        cells=region.get("cells")
        if not isinstance(cells,dict) or len(cells)!=16:
            raise ValueError("Original 16-cell geographic block incomplete")
        eligible={}
        for k in DEPTHS:
            names=[cid for cid,v in cells.items() if eligible_count_at_depth(v,k)]
            n=len(names)
            eligible[str(k)]={
                "n_original_source_eligible_both_periods":n,
                "max_possible_source_pair_comparisons_before_rail_or_road_route_QA":n*(n-1)//2,
                "satisfies_min8_cells_and_min20_upper_pair_opportunity":
                    n>=MIN_INDEPENDENT_SOURCE_CELLS and n*(n-1)//2>=MIN_MAX_PAIR_OPPORTUNITIES,
            }
        counts=[eligible[str(k)]["n_original_source_eligible_both_periods"] for k in DEPTHS]
        if any(b>a for a,b in zip(counts,counts[1:])):
            raise ValueError("Deeper source recording opportunities cannot have more eligible geographic cells")
        output["regions"][block]={
            "country":BLOCKS[block][0],
            "n_original_50km_equal_area_cells":16,
            "by_depth":eligible,
            "NO_train_rail_biological_effect_computed":True,
        }
        print(json.dumps({"block":block,"country":BLOCKS[block][0],
            **{"K"+str(k):eligible[str(k)]["n_original_source_eligible_both_periods"] for k in DEPTHS},
            "K80_min8_source_support":eligible["80"]["satisfies_min8_cells_and_min20_upper_pair_opportunity"]
        },ensure_ascii=False),flush=True)
    for k in DEPTHS:
        counts=[r["by_depth"][str(k)] for r in output["regions"].values()]
        supported=[n for n,x in output["regions"].items() if x["by_depth"][str(k)]["satisfies_min8_cells_and_min20_upper_pair_opportunity"]]
        output["source_summary"][str(k)]={
            "n_geographic_cells_quality_eligible":sum(x["n_original_source_eligible_both_periods"] for x in counts),
            "n_geographical_country_blocks_passing_pre_rail_source_minimum":len(supported),
            "n_distinct_countries_among_eligible_blocks":len({BLOCKS[b][0] for b in supported}),
            "supported_predefined_country_blocks":supported,
            "cannot_infer_pairwise_rail_or_road_network_positivity_from_flora_coverage":True
        }
    output["complete"]=True
    output["limits"]=[
        "These are GBIF positive indexed source RECORD opportunity counts, not true species observed absence, standardized survey, population abundance or effect.",
        "K80 country sample size can be a maximum source feasibility: actual physical RAIL and ROAD path overlap may remove additional cells/pairs.",
        "Original 24 geographical blocks were frozen before outcome; only the 10 existing source-harvested regions are measured in this coverage audit.",
        "Countries with limited GRIIS checklist completeness (JP) and sparse eligible source records (NZ) are not concluded ecologically low in introduced taxa.",
        "Source data provenance is v3.6 historical index; full DOI-event raw occurrence and current recorded-by/dataset effort are not deduplicated.",
        "No physical railway train service, freight or propagule transfer is used or estimated in this source-only feasibility audit.",
    ]
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",required=True)
    p.add_argument("--output",default="results/audit/frozen_10macroblock_K80_source_feasibility_v54.json")
    a=p.parse_args()
    flora=json.loads(Path(a.flora).read_text(encoding="utf-8"))
    result=analyze(flora)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"source_summary":result["source_summary"]},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

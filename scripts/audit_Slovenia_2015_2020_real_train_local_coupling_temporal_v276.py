#!/usr/bin/env python3
"""v27.6 2015 vs 2020 authentic historical TRAIN counts against IDENTICAL
2022/23 Slovenian author-original railway microhabitat local species coupling.

Source author vegetation 89 real sites, 15 true line sections (173 BS/BT
matched sides) and official 2023 station chainage PDF are SHA authenticated.
Only 27 sites on 9 original true undirected 2020 segments and 5 actual
national lines have traffic source exposure. Original Eurostat 2015 and
2020 directional network code + original text label must match exactly.
2015 lacks exact source segments => abort; do NOT impute 2015=2020.

The ecological response is the same exact within-line frequency-preserving
checkerboard-null local species coupling from v27.5, computed over all
89 true Slovenian sites, not re-derived separately by historical year.
For each historical predictor year, physically same 9 interval × 5-line
source plant cohort, evaluate:
  - raw source segment-level Spearman of log goods/passenger vs coupling;
  - true source within-national-line de-meaned freight slope,
    unadjusted or adjusted for within-line passenger log-count;
  - all 48 within-national-line exact goods-label reassignment
    diagnostic, not experimental p.
Changed source year is technical temporal sensitivity, NEVER second
independent vegetation study, and 2015/2020 both precede 2022/23
vegetation. No traffic origin-destination, seed carriage, air gust, mowing
date, genetic/phenotypic response or actual 2022/23 train movements.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_rail_local_site_coupling_v257 import slovenia_location_sets
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import (
    query_series,unpack_jsonstat,spearman_ties_pairs
)
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
    eurostat_2020_segments
)
from validate_SI_historical_2015_2020_train_traffic_vs_flora_v235 import (
    matched_traffic_for_2020_site_segments,EUROSTAT_SHA
)
from analyze_Slovenia_historical_train_frequency_local_species_coupling_v275 import (
    site_specific_null,mapped_segment_data,join_and_group,traffic_diagnostic
)

DRAWS=99
FREEZE_JACCARD_SIMILARITY=.3632363548846349


def historical_site_rematch(rows2020:list[dict],
                            original2020:list[dict],
                            original2015:list[dict])->tuple[list[dict],dict]:
    ids=[{"eurostat_segment_code":x["source_segment"]} for x in rows2020]
    sources,qa=matched_traffic_for_2020_site_segments(ids,original2020,original2015)
    if len(sources)!=9:
        raise ValueError("Cannot match all nine true 2015 and 2020 physical source segments")
    original20={x["source_segment"] for x in rows2020}
    if set(sources)!=original20:
        raise ValueError("Source 2015 exactly matched physical segment set drifted")
    updated=[{
        **s,"freight":sources[s["source_segment"]]["goods"],
        "passenger":sources[s["source_segment"]]["passenger"]
    } for s in rows2020]
    if len(updated)!=27:
        raise ValueError("Source mapped historical site count must remain 27")
    return updated,{
        "actual_2015_true_2direction_source_segments":len(sources),
        "actual_2015_original_vegetation_sites_with_exact_same_physical_segment":len(updated),
        "2015_and_2020_source_codes_plus_original_labels_strict":True,
        "no_other_line_unmatched_2015_freight_imputed":True
    }


def paired_historical(segments2015:list[dict],segments2020:list[dict])->dict:
    if len(segments2015)!=9 or len(segments2020)!=9:
        raise ValueError("Both historical years must share the exact nine physical source intervals")
    key=lambda x:(x["original_line"],x["source_plant_site_count"],
                  round(x["mean_local_excess"],12),round(x["mean_observed_Jaccard_beta"],12))
    q15=sorted(segments2015,key=key)
    q20=sorted(segments2020,key=key)
    # When multiple identical plant score/location counts coincide in same
    # line, sorted key may be ambiguous. Require unique source score keys.
    if len({key(x) for x in q20})!=9 or len({key(x) for x in q15})!=9:
        raise ValueError("Cannot uniquely pair nine original physical segments by source line and fixed plant score")
    if [key(x) for x in q15]!=[key(x) for x in q20]:
        raise ValueError("2015 and 2020 plant outcome drift: only historical traffic predictors may change")
    goods=[(x["goods_2020"],y["goods_2020"]) for x,y in zip(q15,q20)]
    pax=[(x["passenger_2020"],y["passenger_2020"]) for x,y in zip(q15,q20)]
    out={
        "original_identical_physical_segment_plant_response_in_two_years":True,
        "original_source_physical_segments":9,
        "2015_2020_same_physical_segments_goods_spearman":spearman_ties_pairs(goods),
        "2015_2020_same_physical_segments_passenger_spearman":spearman_ties_pairs(pax),
        "actual_freight_2015_min_max":[min(a for a,b in goods),max(a for a,b in goods)],
        "actual_freight_2020_min_max":[min(b for a,b in goods),max(b for a,b in goods)],
        "no_2022_2023_field_year_train_intensity_measured":True,
    }
    return out


def run(original_xlsx:Path,source_audit:Path,station_pdf:Path,draws:int=DRAWS)->dict:
    author=extract_original_xlsx(original_xlsx,source_audit)
    original_site,qa=slovenia_location_sets(author)
    local_sites,localqa=site_specific_null(original_site,draws=draws)
    field20,siteqa=mapped_segment_data(author,station_pdf)
    segment20,join20=join_and_group(field20,local_sites)
    if len(segment20)!=9 or join20["original_mapped_sites_with_valid_frequency_fixed_null"]!=27:
        raise ValueError("Original 2020 27 plant sites/9 traffic segments plus local null unsupported")
    raw2020,qa20=eurostat_2020_segments()
    if qa20["source"]["sha256_original_json_response"]!=EUROSTAT_SHA["2020"]:
        raise ValueError("Frozen true 2020 Eurostat original directional JSON source hash changed")
    raw2015,qa15=query_series("rail_tf_ns15_si")
    if qa15["sha256_original_json_response"]!=EUROSTAT_SHA["2015"]:
        raise ValueError("Frozen true 2015 Eurostat original directional JSON source hash changed")
    direction15,dirqa=unpack_jsonstat(raw2015,"2015")
    site15,matchqa=historical_site_rematch(field20,raw2020,direction15)
    segment15,join15=join_and_group(site15,local_sites)
    if len(segment15)!=9 or join15["original_mapped_sites_with_valid_frequency_fixed_null"]!=27:
        raise ValueError("Strict same 27 original plant sites/9 physical traffic intervals across 2015/2020 failed")
    timeline=paired_historical(segment15,segment20)
    years={
        "2015":traffic_diagnostic(segment15),
        "2020":traffic_diagnostic(segment20)
    }
    return {
        "version":"v27.6",
        "source_genuine_Slovenian_plant_DOI":"10.1007/s44473-026-00155-z",
        "original_author_rail_site_frequency_local_beta_source_QA":{k:v for k,v in localqa.items() if k!="original_line_source_support"},
        "actual_2020_official_2023_station_pdf_and_Eurostat_source_QA":siteqa,
        "historical_directional_source_2015_2020_SHA256":EUROSTAT_SHA,
        "same_2015_2020_true_physical_segment_crosswalk_QA":matchqa,
        "same_original_botany_2015_2020_source_comparison":timeline,
        "exact_same_9segment_5line_source_local_coupling_by_historical_train_year":years,
        "source_original_plant_2022_2023_not_direct_2015_or_2020_traffic_response":True,
        "source_historical_2015_2020_comparison_not_two_independent_ecological_replications":True,
        "not_a_causal_train_frequency_local_seed_coupling_or_gust_effect":True,
        "no_original_field_site_species_GPS_or_Eurostat_direction_table_reexport":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",type=Path,required=True)
    p.add_argument("--source-audit",type=Path,required=True)
    p.add_argument("--official-station-pdf",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    args=p.parse_args()
    result=run(args.original_xlsx,args.source_audit,args.official_station_pdf,args.draws)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
        "version":"v27.6",
        "actual_historical_traffic_site_correspondence":result["same_2015_2020_true_physical_segment_crosswalk_QA"],
        "genuine_2015_2020_nine_segment_rank_QA":result["same_original_botany_2015_2020_source_comparison"],
        "actual_source_traffic_year_diagnostic":{
            year:{
                "historical_goods_Spearman":
                    v["source_2020_historical_goods_trains_year_9segment_spearman_local_excess"],
                "historical_passenger_Spearman":
                    v["source_2020_historical_passenger_trains_year_9segment_spearman_local_excess"],
                "within_line_freight_source_reassignment":
                    v["within_line_goods_slope_with_source_exact_48_label_redistributions"],
                "within_line_passenger_adjusted_freight_source_reassignment":
                    v["within_line_goods_adjusted_for_passenger_slope_with_exact_48"]
            } for year,v in result["exact_same_9segment_5line_source_local_coupling_by_historical_train_year"].items()
        }
    },ensure_ascii=False))
if __name__=="__main__":main()

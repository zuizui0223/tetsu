#!/usr/bin/env python3
"""v12.2: original CzechVeg RAIL-text vs ROAD-text data provenance
confounding by ORIGINAL source study/reference title.

Critical source-unit question: did rail- and road-labelled vegetation
quadrats appear WITHIN the same original published field survey, rather
than being assembled from incompatible studies? Source original
CzechVeg-Open v1.1 (2026) 117739 plot metadata, authentic Zenodo MD5,
Chytrý et al. DOI 10.5281/zenodo.19231152. Original publication
REFERENCE_TITLE metadata is not actual field management or vegetation
treatment assignment. No source species rows or physical seed movements
used. All matched and nonmatched source exclusions disclosed.
"""
from __future__ import annotations
import argparse,csv,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    plot_labels,verify_file,HEADERS_NAME,
)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    year_from_source,area_key,
)
from analyze_CzechVeg_rail_road_region_robustness_v120 import validated_source

BAD_TITLE=frozenset(("", "NA", "N/A","NULL","UNKNOWN","UNKNOWN SOURCE","NONE","NAN"))
THRESHOLDS=(1,2,3,5)


def source_title(t):
    if t is None:return None
    s=" ".join(str(t).strip().split())
    if s.upper() in BAD_TITLE:return None
    if len(s)<12:return None
    return s


def provenance_blocks(headers_file,plots):
    by={mode:{
        "study":defaultdict(lambda:{"rail":set(),"road":set()}),
        "study_region":defaultdict(lambda:{"rail":set(),"road":set()}),
        "study_region_decade_area":defaultdict(lambda:{"rail":set(),"road":set()}),
        "no_study_title":Counter(),
        "surveyed_source_count":Counter(),
    } for mode in ("strict","broad")}
    columns=[]
    with open(headers_file,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        columns=reader.fieldnames or []
        required={"RELEVE_NR","REFERENCE_TITLE","YEAR","SKALICKY1988"}
        if required-set(columns):
            raise ValueError("Original Czech vegetation public-study bibliographic keys missing")
        n=0
        for original in reader:
            n+=1
            lid=str(original.get("RELEVE_NR") or "").strip()
            if lid not in plots:continue
            plot=plots[lid]
            study=source_title(original.get("REFERENCE_TITLE"))
            region=str(original.get("SKALICKY1988") or "").strip()
            year=year_from_source(original.get("YEAR"))
            for mode in ("strict","broad"):
                cat=plot[mode]
                if cat not in ("rail","road"):continue
                a=by[mode]
                a["surveyed_source_count"][cat]+=1
                if study is None:
                    a["no_study_title"][cat]+=1
                    continue
                a["study"][study][cat].add(lid)
                if region not in ("","NA"):
                    a["study_region"][(study,region)][cat].add(lid)
                    if year is not None:
                        a["study_region_decade_area"][
                            (study,region,year//10*10,area_key(plot["area"],"class4"))
                        ][cat].add(lid)
    if n!=117739:
        raise ValueError("Original CzechVeg 117739 source plot headers changed")
    return by,{"n_original_full_header_rows":n,
               "source_original_reference_title_column_confirmed":True,
               "original_source_year_and_area_metadata_only_no_outcome":True}


def summarize(blocks,mode):
    totals=Counter()
    summaries={}
    for key in ("study","study_region","study_region_decade_area"):
        groups=blocks[key]
        possible=len(groups)
        subset={}
        for minimum in THRESHOLDS:
            ok={k:v for k,v in groups.items()
                if len(v["rail"])>=minimum and len(v["road"])>=minimum}
            ids={site for r in ok.values() for site in r["rail"]+r["road"]} if False else None
            summaries[key+"__min"+str(minimum)]={
                "n_source_study_strata_with_both_original_rail_and_road":len(ok),
                "source_rail_candidate_plot_memberships_in_shared_source_study":sum(
                    len(v["rail"]) for v in ok.values()),
                "source_road_candidate_plot_memberships_in_shared_source_study":sum(
                    len(v["road"]) for v in ok.values()),
                "original_rail_candidate_plot_distinct_ids_in_shared_source_study":len(
                    {x for v in ok.values() for x in v["rail"]}),
                "original_road_candidate_plot_distinct_ids_in_shared_source_study":len(
                    {x for v in ok.values() for x in v["road"]}),
                "source_same_original_author_study_enough_for_both_source_categories":len(ok)>0,
            }
        summaries[key+"__candidate_study_strata_total"]=possible
    both=[(study,v) for study,v in blocks["study"].items()
          if len(v["rail"])>=2 and len(v["road"])>=2]
    both.sort(key=lambda z:(-min(len(z[1]["rail"]),len(z[1]["road"])),-sum(len(v) for v in z[1].values())))
    return {
        "source_original_site_label_type":mode,
        "n_source_rail_candidate_sites_with_original_valid_plot_area":
            blocks["surveyed_source_count"]["rail"],
        "n_source_road_candidate_sites_with_original_valid_plot_area":
            blocks["surveyed_source_count"]["road"],
        "n_source_rail_site_plots_without_original_study_reference_title":
            blocks["no_study_title"]["rail"],
        "n_source_road_site_plots_without_original_study_reference_title":
            blocks["no_study_title"]["road"],
        "matching_tier_counts":summaries,
        "original_same_study_rail_vs_road_top10_metadata_only":[
            {"original_source_reference_title_excerpt":study[:140],
             "n_source_rail_sites":len(v["rail"]),
             "n_source_road_sites":len(v["road"])}
             for study,v in both[:10]],
        "original_study_title_is_not_a_randomization_or_site_distance_control":True,
    }


def run(headers,meta,verified119):
    prev=validated_source(verified119)
    true=verify_file(meta,HEADERS_NAME,headers)
    if true["checksum"]!=prev["original_source_Czechveg_header_provenance"]["checksum"]:
        raise ValueError("Real original Czech vegetation plot header source no longer matches v11.9")
    plots,locqc=plot_labels(headers)
    a,details=provenance_blocks(headers,plots)
    return {
        "complete":True,"version":"v12.2",
        "original_Zenodo_DOI":"10.5281/zenodo.19231152",
        "original_verified_source_headers":true,
        "original_site_area_and_label_audit":locqc,
        "original_study_provenance_audit":details,
        "strict_HABITAT_rail_vs_road":summarize(a["strict"],"strict"),
        "broad_HABITAT_LOCALITY_rail_vs_road":summarize(a["broad"],"broad"),
        "no_biological_estimation":True,
        "source_study_validation_limits":[
            "Source REFERENCE_TITLE can be missing or generic even when survey was performed properly.",
            "Being in same paper does not mean rail and road plots match vegetation/habitat/year/area.",
            "Different source papers with no common rail and road plots cannot identify a railway-caused wind-dispersal contrast.",
            "Original published source 117739 plot metadata scanned, no raw source files redistributed.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for opt in ("zenodo-record","source-headers","verified119","output"):
        p.add_argument("--"+opt,required=True)
    a=p.parse_args()
    res=run(a.source_headers,
            json.loads(Path(a.zenodo_record).read_text(encoding="utf-8")),
            json.loads(Path(a.verified119).read_text(encoding="utf-8")))
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"version":"v12.2",
        "original_117739_plot_study":res["original_study_provenance_audit"],
        "strict":res["strict_HABITAT_rail_vs_road"],
        "broad":res["broad_HABITAT_LOCALITY_rail_vs_road"]},
        ensure_ascii=False),flush=True)
if __name__=="__main__":main()

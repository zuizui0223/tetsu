#!/usr/bin/env python3
"""v12.4 genuine original CzechVeg vegetation SYNTAXON rank/source audit.

Original external field source CzechVeg-Open v1.1 2026, Zenodo
DOI 10.5281/zenodo.19231152, 117739 original phytosociology plot
headers verified against original publisher MD5 + v11.9 prior receipt.

Question: Are publication-matched RAIL and ROAD source text candidate
plots also comparable by ORIGINAL SOURCE community syntaxon, without
inventing or converting unknown Czech vegetation classification ranks?

Do NOT interpret raw syntaxon ID as a common class, phytosociological
rank, habitat control, managed railway or plant fitness. Report exact
raw code & name availability and exact source label overlap only.
"""
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path
from collections import Counter,defaultdict
from audit_CzechVeg_shared_railroad_reference_studies_v122 import source_title
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file,HEADERS_NAME,plot_labels
)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    year_from_source,area_key
)
from analyze_CzechVeg_rail_road_region_robustness_v120 import validated_source

SYNTAX_COLUMNS=("SYNTAXON","SYNTAXON_NAME","SYNTAX_OLD","SYNTAXON_OLD_NAME")
NULLS={"","NA","N/A","NONE","UNKNOWN","NULL","NAN","-","?"}
MATCH_LEVELS=(
    ("same_paper_phyto_decade_area_SYNTAXON",True,True,True,"SYNTAXON"),
    ("same_paper_phyto_decade_area_SYNTAXON_NAME",True,True,True,"SYNTAXON_NAME"),
    ("same_paper_phyto_decade_area_SYNTAX_OLD",True,True,True,"SYNTAX_OLD"),
    ("same_paper_phyto_decade_area_SYNTAXON_OLD_NAME",True,True,True,"SYNTAXON_OLD_NAME"),
    ("same_paper_phyto_decade_SYNTAXON",True,True,False,"SYNTAXON"),
    ("same_paper_phyto_SYNTAXON",True,False,False,"SYNTAXON"),
    ("same_paper_phyto_decade_area_no_syntax",True,True,True,None),
)
MIN_SIDES=(1,2,3)


def clean_syntaxon(value):
    if value is None:return None
    s=" ".join(str(value).strip().split())
    if s.upper() in NULLS or len(s)>400:return None
    return s


def source_meta_with_taxonomic_vegetation(header, eligible):
    active={}
    bycolumn={col:Counter() for col in SYNTAX_COLUMNS}
    ex={col:[] for col in SYNTAX_COLUMNS}
    quality=Counter()
    with open(header,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        fields=reader.fieldnames or []
        if set(SYNTAX_COLUMNS)-set(fields):
            raise ValueError("Original published Czech vegetation syntaxon columns unavailable")
        if "REFERENCE_TITLE" not in fields:
            raise ValueError("No original paper-studies source metadata")
        for row in reader:
            quality["n_full_original_CzechVeg_source_header_rows"]+=1
            rid=str(row.get("RELEVE_NR") or "").strip()
            if rid not in eligible:continue
            src=eligible[rid]
            if src["strict"] not in ("rail","road") and src["broad"] not in ("rail","road"):
                continue
            study=source_title(row.get("REFERENCE_TITLE"))
            region=str(row.get("SKALICKY1988") or "").strip()
            year=year_from_source(row.get("YEAR"))
            if study is None:
                quality["rail_or_road_source_plot_reference_study_missing"]+=1
                continue
            if not region or region=="NA":
                quality["rail_or_road_source_plot_phyto_region_missing"]+=1
                continue
            if year is None:
                quality["rail_or_road_source_plot_valid_survey_year_missing"]+=1
                continue
            syntax={k:clean_syntaxon(row.get(k)) for k in SYNTAX_COLUMNS}
            active[rid]={
                "study":study,"source_region":region,"survey_year":year,
                "area":src["area"],"strict":src["strict"],"broad":src["broad"],
                **syntax
            }
            quality["n_source_candidate_plot_with_paper_region_year_and_area"]+=1
            for k,s in syntax.items():
                if s is not None:
                    bycolumn[k][s]+=1
                    quality["n_source_candidate_plot_with_"+k]+=1
                    if len(ex[k])<8:ex[k].append(s[:120])
                else:quality["n_source_candidate_plot_missing_"+k]+=1
    if quality["n_full_original_CzechVeg_source_header_rows"]!=117739:
        raise ValueError("Original source 117739 plot rows changed")
    return active,{
        "source_header_rows":quality["n_full_original_CzechVeg_source_header_rows"],
        "retained_source_plot_with_original_study_phyto_year_area":len(active),
        "source_plot_eligibility_and_missing_syntaxon_counts":dict(quality),
        "source_vegetation_column_distinct_label_counts":{
            k:len(bycolumn[k]) for k in SYNTAX_COLUMNS},
        "source_vegetation_column_original_raw_top12":{
            k:[{"original_raw_code_or_label":v[:120],"n_source_candidate_plots":n}
               for v,n in count.most_common(12)] for k,count in bycolumn.items()},
        "original_source_vegetation_column_examples_first8":ex,
        "NO_unverified_conversion_to_common_Czech_vegetation_class":True,
    }


def block_key(row,level):
    title,with_region,with_year,with_area,syn=level
    return (row["study"],
            row["source_region"] if with_region else None,
            row["survey_year"]//10*10 if with_year else None,
            area_key(row["area"],"class4") if with_area else None,
            row[syn] if syn else None)


def original_source_cross_arm_support(active,mode,level,minimum):
    bins=defaultdict(lambda:{"rail":set(),"road":set()})
    counts=Counter()
    syn=level[4]
    for rid,rec in active.items():
        label=rec[mode]
        if label not in ("rail","road"):continue
        counts["candidate_source_"+label]+=1
        if syn is not None and rec[syn] is None:
            counts["candidate_"+label+"_missing_source_syntaxon"]+=1
            continue
        k=block_key(rec,level)
        bins[k][label].add(rid)
    keep={k:v for k,v in bins.items() if
          len(v["rail"])>=minimum and len(v["road"])>=minimum}
    studies={k[0] for k in keep}
    regions={k[1] for k in keep if k[1] is not None}
    return {
        "n_original_source_both_arm_blocks":len(keep),
        "n_distinct_original_publication_studies_retained":len(studies),
        "n_original_phyto_regions_retained":len(regions),
        "source_rail_sites_in_common_vegetation_support":
            len({p for v in keep.values() for p in v["rail"]}),
        "source_road_sites_in_common_vegetation_support":
            len({p for v in keep.values() for p in v["road"]}),
        "balanced_plot_pair_mass_minimum":
            sum(min(len(v["rail"]),len(v["road"])) for v in keep.values()),
        "candidate_original_rail_syntaxon_missing":
            counts["candidate_rail_missing_source_syntaxon"],
        "candidate_original_road_syntaxon_missing":
            counts["candidate_road_missing_source_syntaxon"],
        "all_original_plot_sampling_study_area_year_source_unchanged":True
    }


def run(meta,header,verified119):
    report=validated_source(verified119)
    hdr=verify_file(meta,HEADERS_NAME,header)
    if hdr["checksum"]!=report["original_source_Czechveg_header_provenance"]["checksum"]:
        raise ValueError("Original CzechVeg source header MD5 differs from prior verified source")
    plots,plot_quality=plot_labels(header)
    actual,quality=source_meta_with_taxonomic_vegetation(header,plots)
    output={}
    for mode in ("strict","broad"):
        output[mode]={}
        for level in MATCH_LEVELS:
            key=level[0]
            output[mode][key]={}
            for min_each in MIN_SIDES:
                output[mode][key][str(min_each)]=original_source_cross_arm_support(
                    actual,mode,level,min_each)
    return {
        "complete":True,"version":"v12.4",
        "original_117739_CzechVeg_survey_DOI":"10.5281/zenodo.19231152",
        "original_original_author_pub_type":"CzechVeg-Open v1.1 2026",
        "original_source_checksum":hdr,
        "original_source_syntaxon_label_audit":quality,
        "source_study_plot_crosstype_common_vegetation_support":output,
        "not_plant_trait_biological_model":True,
        "source_vegetation_rank_not_invented":True,
        "limitations":[
            "CzechVeg source SYNTAXON codes and names have not been decoded into phyto sociological class, alliance or community rank; they are ONLY exact original text strata.",
            "Matching exact source vegetation codes may over-restrict original sampling and differ in classification level across original botanical papers.",
            "Original source RAIL and ROAD site text is an unvalidated site-context proxy for railway right-of-way.",
            "Same original publication does NOT guarantee identical original field survey design, randomized exposure or actual train seed transport.",
            "No source species physiology, mowing treatment, train slipstream, seed-trap, gene flow, flower phenology measured.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("zenodo-record","source-headers","verified119","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    o=run(json.loads(Path(a.zenodo_record).read_text(encoding="utf8")),
          a.source_headers,
          json.loads(Path(a.verified119).read_text(encoding="utf8")))
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(o,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({
        "complete":True,"version":"v12.4",
        "source":o["original_source_syntaxon_label_audit"],
        "eligible_blocks":{
          mode:{
            level:{
              min_each:{
                k:v for k,v in val.items() if
                k in ("n_original_source_both_arm_blocks",
                      "n_distinct_original_publication_studies_retained",
                      "n_original_phyto_regions_retained",
                      "source_rail_sites_in_common_vegetation_support",
                      "source_road_sites_in_common_vegetation_support",
                      "balanced_plot_pair_mass_minimum",
                      "candidate_original_rail_syntaxon_missing",
                      "candidate_original_road_syntaxon_missing")
              } for min_each,val in limits.items()}
              for level,limits in group.items()}
              for mode,group in o["source_study_plot_crosstype_common_vegetation_support"].items()}
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

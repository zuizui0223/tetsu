#!/usr/bin/env python3
"""v19.4: source-original abiotic vs vegetation-class regional reference audit.

This is an opportunity test ONLY. The CzechVeg SYNTAXON TH/TF strata
are derived from botanical composition; choosing them as reference
could mechanically make reference communities similar. To prevent
selection on the plant-response, count alternative reference plots
selected by field-measured ALTITUDE and INCLINATIO, with optional
bare rock and non-standardized SOIL/GEOLOGY metadata audits.
COV_ROCK is a substrate cover measurement but may mediate vegetation.
No source plant abundance or community effects estimated here.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    get_meta,matched_rail_source,ANALOG_CODES
)

EXPECTED_ORIGINAL_HEADERS=117739
NUM_FIELDS=("ALTITUDE","INCLINATIO","COV_ROCK","PRECISION")
STRING_FIELDS=("GEOLOGY","SOIL")
MIN_REF=5
SOURCE_STUDY_QUASI_INDEPENDENCE=True

def number(raw,min_=None,max_=None):
    s=str(raw or "").strip()
    if s in {"","-1","NA","N/A","NaN","nan","NULL"}:return None
    try:n=float(s.replace(",",".")) 
    except ValueError:return None
    if not math.isfinite(n):return None
    if min_ is not None and n<min_:return None
    if max_ is not None and n>max_:return None
    return n

def abiotic_rows(path,source):
    measured={};audit=Counter()
    with open(path,encoding="utf-8-sig",newline="") as f:
        rows=csv.DictReader(f)
        required=set(NUM_FIELDS)|set(STRING_FIELDS)|{"RELEVE_NR"}
        if required-set(rows.fieldnames or ()):
            raise ValueError("Original CzechVeg-Open missing published environmental header columns")
        for r in rows:
            audit["total_headers"]+=1
            pid=str(r.get("RELEVE_NR") or "").strip()
            if pid not in source:continue
            vals={
                "elev_m":number(r.get("ALTITUDE"),0,2500),
                "slope_deg":number(r.get("INCLINATIO"),0,90),
                "rock_cover_pct":number(r.get("COV_ROCK"),0,100),
                "position_uncertainty_m":number(r.get("PRECISION"),0,1_000_000),
                "geology_recorded":str(r.get("GEOLOGY") or "").strip() not in ("","NA"),
                "soil_recorded":str(r.get("SOIL") or "").strip() not in ("","NA"),
            }
            measured[pid]=vals
            audit["eligible_original_source_headers"]+=1
            for fld,v in vals.items():
                if v is not None and v is not False:audit["recorded_"+fld]+=1
    if audit["total_headers"]!=EXPECTED_ORIGINAL_HEADERS or len(measured)!=len(source):
        raise ValueError("Publisher original header count or source ID join mismatch")
    return measured,dict(audit)

def source_reference_options(meta,blocks):
    candidate=defaultdict(list)
    for pid,m in meta.items():
        if m["label"] is None and m["strict"] is None:
            candidate[(m["region"],m["decade"],m["area"])].append(pid)
    out={}
    for key,v in blocks.items():
        study,reg,dec,area=key
        options=[s for s in candidate[(reg,dec,area)]
                 if meta[s]["study"]!=study]
        for rid in v["rail"]:
            out[rid]=options
    return out

def matching(source,ctrl):
    """Predeclared abiotic matching, without access to vegetation/taxon outcomes."""
    if source["elev_m"] is None:return []
    out=[]
    diff=abs(source["elev_m"]-ctrl["elev_m"]) if ctrl["elev_m"] is not None else float("inf")
    if diff<=200:out.append("elev200")
    if diff<=100:
        out.append("elev100")
        if source["slope_deg"] is not None and ctrl["slope_deg"] is not None:
            if abs(source["slope_deg"]-ctrl["slope_deg"])<=10:
                out.append("elev100_slope10")
                if (source["rock_cover_pct"] is not None
                    and ctrl["rock_cover_pct"] is not None
                    and abs(source["rock_cover_pct"]-ctrl["rock_cover_pct"])<=20):
                    out.append("elev100_slope10_rock20")
        if ctrl["rock_cover_pct"] is not None and ctrl["rock_cover_pct"]>=10:
            out.append("elev100_reference_rock_ge10")
    return out

MATCHES=("elev200","elev100","elev100_slope10",
         "elev100_slope10_rock20","elev100_reference_rock_ge10")

def audit(meta,blocks,measured):
    ropts=source_reference_options(meta,blocks)
    stats={typ:Counter() for typ in MATCHES}
    by_region={typ:defaultdict(Counter) for typ in MATCHES}
    metaq=Counter()
    for rid,opts in ropts.items():
        orig=measured[rid]
        metaq["original_rail_text_plots"]+=1
        for col in ("elev_m","slope_deg","rock_cover_pct","position_uncertainty_m","geology_recorded","soil_recorded"):
            if orig[col] is not None and orig[col] is not False:metaq["rail_has_"+col]+=1
        qualified=defaultdict(set)
        for r in opts:
            for typ in matching(orig,measured[r]):
                qualified[typ].add(r)
        for typ in MATCHES:
            counts=stats[typ];reg=meta[rid]["region"]
            values=qualified[typ]
            counts["rail_plots_scanned"]+=1
            counts["total_source_paired_control_opportunities"]+=len(values)
            if values:
                counts["rail_plots_with_at_least_one_reference"]+=1
            if len(values)>=MIN_REF:
                counts["rail_plots_with_min5_references"]+=1
                by_region[typ][reg]["rail_plots_min5"]+=1
                counts["source_five_reference_sites_from_distinct_papers"]+=len(
                    {meta[c]["study"] for c in values})>=2
                counts["source_analog_class_reference_opportunities"]+=sum(
                    meta[c]["source_class"] in ANALOG_CODES for c in values)
            by_region[typ][reg]["all_rail_plots"]+=1
    summary={}
    for typ in MATCHES:
        counts=stats[typ]
        n_regions=sum(x["rail_plots_min5"]>0 for x in by_region[typ].values())
        summary[typ]={
            **dict(counts),"original_phyto_regions_supported_min5":n_regions,
            "min_reference_opportunity":MIN_REF,
            "rail_sites_are_original_text_not_verified_right_of_way":True,
            "source_environmental_match_independent_of_taxon_or_syntaxon":True,
            "vegetation_source_THTF_label_used_for_filter":False,
            "source_region_rail_plots_min5_no_plot_ids":{
                str(k):v["rail_plots_min5"] for k,v in sorted(by_region[typ].items()) if v["rail_plots_min5"]>0
            }
        }
    return {"meta_rail_numeric_coverage":dict(metaq),
            "source_reference_type_support":summary,
            "all_original_source_rail_plots_are_from_same_original_publication_as_road_arm":True,
            "reference_plots_are_from_different_original_publication_than_rail":True,
            "nonrail_reference_status_not_confirmed_by_absent_transport_text":True,
            "geology_and_soil_text_are_unstandardized":True,
            "rock_cover_percent_may_be_ecological_response_or_management_mediator":True,
            "species_occurrence_response_NOT_used":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-headers",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    meta,proof=get_meta(a.original_headers)
    blocks,support=matched_rail_source(meta)
    measured,qa=abiotic_rows(a.original_headers,meta)
    result={"version":"v19.4", "original_source_DOI":"10.5281/zenodo.19231152",
            "original_source_headers_md5":proof["publisher_original_headers_md5"],
            "original_header_count":qa["total_headers"],
            "source_eligible_header_count":len(meta),
            "source_matched_railroad_frame":support,
            "numeric_field_quality":qa,
            "abiotic_ref_opportunity":audit(meta,blocks,measured),
            "no_taxon_or_phenotype_effect_estimated":True}
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({
        "version":result["version"],
        "real_original_headers":qa["total_headers"],
        "abiotic_rail_data_coverage":result["abiotic_ref_opportunity"]["meta_rail_numeric_coverage"],
        "support":{k:{kk:vv for kk,vv in v.items() if kk!="source_region_rail_plots_min5_no_plot_ids"}
             for k,v in result["abiotic_ref_opportunity"]["source_reference_type_support"].items()}
    },ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()

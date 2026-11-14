#!/usr/bin/env python3
"""v12.6 audit: original CzechVeg hierarchical SYNTAXON-to-CLASS exact
source-code prefix, preserving source country botanical classification.

The original published Vegetation of Czech Republic syntaxon codes
have first character formation group, second character vegetation
class, third character alliance, subsequent numeric association.
Published Pladias official vegetation catalogue displays XD=Galio-
Urticetea, XDA=Senecionion fluviatilis, XDA01=Cuscuto... etc.
Source publisher CzechVeg v1.1 original 117739 plot metadata
MD5 verified. Interpret 2-character canonical prefix as class
ONLY if that prefix is independently observed as an EXACT 2-letter
standalone code somewhere in this SAME original 117739-row source
header table. Missing/invalid codes remain unknown. This avoids
inventing genus/vegetation codes and tests sample eligibility before
re-fitting plant guild effects.

This does NOT prove source plots are active railway right-of-way nor
direct management, train-gust, viable freight seed transport or genetic
adaptation. Conditioning on vegetation class might condition on a
railway-disturbance consequence; no total causal effect identified.
"""
from __future__ import annotations
import argparse,csv,json,re
from collections import Counter
from pathlib import Path
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file,HEADERS_NAME,plot_labels,
)
from analyze_CzechVeg_rail_road_region_robustness_v120 import validated_source
from audit_CzechVeg_syntaxon_same_study_source_v124 import (
    source_meta_with_taxonomic_vegetation,original_source_cross_arm_support,
)

SOURCE_CODE_RE=re.compile(r"^[A-Z]{2}(?:[A-Z](?:\d{2}[a-z]?)?)?$")
CLASS_CODE_RE=re.compile(r"^[A-Z]{2}$")
SOURCE_MISSING={"","NA","N/A","UNKNOWN","NULL","NONE","NAN","-","?"}


def source_canonical_syntaxon_class(code,anchored_roots):
    if code is None:return None
    v=str(code).strip()
    if v.upper() in SOURCE_MISSING:return None
    if not SOURCE_CODE_RE.fullmatch(v):return None
    group=v[:2]
    return group if group in anchored_roots else None


def canonical_roots_from_original_full_publisher_csv(original_header):
    observed=Counter()
    counts=Counter()
    with open(original_header,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if "SYNTAXON" not in (reader.fieldnames or []):
            raise ValueError("Original CzechVeg SYNTAXON source code missing")
        for row in reader:
            counts["n_source_vegetation_plot_header_rows"]+=1
            code=str(row.get("SYNTAXON") or "").strip()
            if CLASS_CODE_RE.fullmatch(code):
                observed[code]+=1
            if code and code.upper() not in SOURCE_MISSING:
                counts["n_source_vegetation_syntaxon_code_nonmissing"]+=1
    if counts["n_source_vegetation_plot_header_rows"]!=117739:
        raise ValueError("Original 117739 source plot metadata count invalid")
    return set(observed),{
        "n_independent_source_class_code_root_labels":len(observed),
        "n_source_vegetation_plot_headers":counts["n_source_vegetation_plot_header_rows"],
        "n_source_syntaxon_nonmissing_all_original":counts[
            "n_source_vegetation_syntaxon_code_nonmissing"],
        "top_20_original_source_published_class_code_roots":[
            {"source_published_vegetation_class_code":k,
             "n_original_plot_headers_with_direct_class_code":v}
            for k,v in observed.most_common(20)],
    }


def convert_source_syntaxon_meta(meta,roots):
    output={}
    counts=Counter()
    for plot,m in meta.items():
        source=m["SYNTAXON"]
        class_code=source_canonical_syntaxon_class(source,roots)
        if source is not None and class_code is None:
            counts["candidate_with_known_source_syntaxon_code_but_unmapped_class"]+=1
        if class_code:
            counts["candidate_with_verified_original_source_class_code"]+=1
        else:
            counts["candidate_unknown_class_excluded_from_class_match"]+=1
        output[plot]={**m,"SOURCE_CLASS_CODE":class_code}
    return output,dict(counts)


def matched_class_summary(meta,mode,condition,minimum):
    from collections import defaultdict
    d=defaultdict(lambda:{"rail":set(),"road":set()})
    for rid,p in meta.items():
        lab=p[mode]
        if lab not in ("rail","road"):continue
        code=p["SOURCE_CLASS_CODE"]
        if code is None:continue
        if condition=="same_study_region_decade_area_class":
            key=(p["study"],p["source_region"],
                 p["survey_year"]//10*10,
                 ("1-4" if p["area"]<=4 else "5-16" if p["area"]<=16
                  else "17-64" if p["area"]<=64 else "65-256" if p["area"]<=256
                  else "257-1000"),code)
        elif condition=="same_study_region_decade_class":
            key=(p["study"],p["source_region"],p["survey_year"]//10*10,code)
        elif condition=="same_study_region_class":
            key=(p["study"],p["source_region"],code)
        else:raise ValueError("Source class matching condition not explicitly declared")
        d[key][lab].add(rid)
    ok={k:v for k,v in d.items() if min(len(v["rail"]),len(v["road"]))>=minimum}
    return {
        "n_original_same_source_study_vegetation_class_comparison_blocks":len(ok),
        "n_original_studies_with_actual_rail_road_class_overlap":len({k[0] for k in ok}),
        "n_original_source_veg_classes_with_valid_both_arm_blocks":len({k[-1] for k in ok}),
        "n_original_rail_plots_on_class_comparable_support":
            len({site for v in ok.values() for site in v["rail"]}),
        "n_original_road_plots_on_class_comparable_support":
            len({site for v in ok.values() for site in v["road"]}),
        "n_original_phytoregions_with_comparable_source_vegetation_classes":
            len({k[1] for k in ok}),
        "effective_original_block_pair_minimum_mass":
            sum(min(len(v["rail"]),len(v["road"])) for v in ok.values()),
        "matched_source_class_codes":sorted({k[-1] for k in ok}),
    }


def run(meta,header,receipt119,receipt124):
    if not(receipt124.get("complete") is True and receipt124.get("version")=="v12.4"):
        raise ValueError("Require real original CzechVeg v12.4 vegetation source receipt")
    prev=validated_source(receipt119)
    source=verify_file(meta,HEADERS_NAME,header)
    if source["checksum"]!=prev["original_source_Czechveg_header_provenance"]["checksum"]:
        raise ValueError("Original CzechVeg original header MD5 not authenticated")
    roots,rootqc=canonical_roots_from_original_full_publisher_csv(header)
    plots,plotqc=plot_labels(header)
    core,originalqc=source_meta_with_taxonomic_vegetation(header,plots)
    if (len(core)!=receipt124["original_source_syntaxon_label_audit"][
       "retained_source_plot_with_original_study_phyto_year_area"]):
        raise ValueError("Real original publication/study-year candidate source drift")
    enriched,audit=convert_source_syntaxon_meta(core,roots)
    conditions=("same_study_region_decade_area_class",
                "same_study_region_decade_class",
                "same_study_region_class")
    results={
        mode:{condition:{
            str(minimum):matched_class_summary(enriched,mode,condition,minimum)
            for minimum in (1,2,3)
        } for condition in conditions}
        for mode in ("strict","broad")
    }
    return {
        "complete":True,"version":"v12.6",
        "original_source_CzechVeg_DOI":"10.5281/zenodo.19231152",
        "source_vegetation_class_hierarchy_evidence":"https://pladias.cz/en/vegetation/",
        "original_source_file":source,
        "original_published_standalone_two_letter_class_root_audit":rootqc,
        "original_source_study_region_year_metadata_audit":originalqc,
        "source_railroad_candidate_class_availability":audit,
        "original_same_paper_region_year_and_VERIFIED_class_rail_road":results,
        "not_a_genuine_train_microhabitat_experiment":True,
        "classification_and_causal_limits":[
            "First two letters of standard Czech vegetation syntaxon classification identify vegetation class, but original code must match published source pattern and its root class must be found independently as an original standalone class code.",
            "Some original SYNTAXON records use nonstandard codes or missing labels; these are unclassified, not assigned class by intuition.",
            "Vegetation class itself could be a consequence of railway environmental filtering: conditioning on it is not a controlled estimate of the total effect of railway land use.",
            "Source railway/road habitat and locality TEXT clues are unvalidated operational railway boundaries, not measured train wind/seed flux/management.",
            "Original source botanical studies were not randomized into experimental railway and roadside treatments.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("zenodo-record","source-headers","verified119","verified124","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    out=run(
        json.loads(Path(a.zenodo_record).read_text(encoding="utf-8")),
        a.source_headers,
        json.loads(Path(a.verified119).read_text(encoding="utf-8")),
        json.loads(Path(a.verified124).read_text(encoding="utf-8")))
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v12.6",
        "class_root":out["original_published_standalone_two_letter_class_root_audit"],
        "availability":out["source_railroad_candidate_class_availability"],
        "comparison":{
            mode:{typ:{
                minimum:{k:v for k,v in stat.items() if k!="matched_source_class_codes"}
                for minimum,stat in groups.items()}
                for typ,groups in variants.items()}
                for mode,variants in out[
                    "original_same_paper_region_year_and_VERIFIED_class_rail_road"].items()}
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

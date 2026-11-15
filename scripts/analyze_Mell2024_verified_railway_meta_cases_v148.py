#!/usr/bin/env python3
"""v14.8 exact original Springer Mell 2024 railway meta case and Penone bias audit.
Do not pool railway habitat richness/cover contrasts as train seed transport.
Publisher originals sha256-pinned; excerpts only, not workbooks, retained.
"""
import argparse, hashlib, json, math, zipfile
from pathlib import Path
from audit_Mell2024_supplementary_railway_meta_v147 import (
    sheet_files,get_shared_strings,worksheet_rows)
SHA5="85042f37984a6f8f316f1f578e7f3093334ded057ebdc0c72437ac8b34579c90"
SHA8="0fd8ff9922c1138d5abd5289e743058b10be3b52c73d4d338cc9b9d2505def17"
DOI="10.1186/s13750-024-00328-3"

def source_tables(path,sha,required):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:
        raise ValueError("Original Springer publisher supplement SHA changed")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError("Publisher ZIP CRC failure")
        paths=dict(sheet_files(z))
        if not set(required).issubset(paths):
            raise ValueError("Original publisher worksheet schema changed")
        shared=get_shared_strings(z)
        return {name:list(worksheet_rows(z,paths[name],shared)) for name in required}

def as_records(rows,required):
    if not rows or rows[0][0]!=1:raise ValueError("Missing source header")
    h={v.strip():k for k,v in rows[0][1].items()}
    if not set(required).issubset(h):
        raise ValueError("Original source header keys missing "+str(set(required)-h.keys()))
    return [{name:raw.get(i,"") for name,i in h.items()} for _,raw in rows[1:]]

def rail_cases(cases):
    if len(cases)!=205 or len({c["Case.ID"] for c in cases})!=205:
        raise ValueError("Original 205 meta analysis cases missing or duplicated")
    rail=[c for c in cases if c["LTI"].strip().lower()=="railway"]
    expected={("EM_8","Agence_Gaiadomo_2012_1"),("11955","Everingham_2019_3"),
              ("11955","Everingham_2019_6")}
    if {(c["Study.ID"],c["Case.ID"]) for c in rail}!=expected:
        raise ValueError("Original 3 railway cases and 2 studies changed")
    summaries=[]
    for c in sorted(rail,key=lambda x:x["Case.ID"]):
        n1,n0=int(float(c["N_treatment"])),int(float(c["N_control"]))
        m1,m0=float(c["Mean_treatment"]),float(c["Mean_control"])
        sd1,sd0=float(c["Var_treatment"]),float(c["Var_control"])
        if (c["Measure_var_treatment"],c["Measure_var_control"])!=("SD","SD"):
            raise ValueError("Source variability does not represent standard deviations")
        if min(n1,n0)<2 or not all(math.isfinite(x) and x>=0 for x in (sd1,sd0)):
            raise ValueError("Unsafe source counts or SD")
        pool=((n1-1)*sd1**2+(n0-1)*sd0**2)/(n1+n0-2)
        g=(1-3/(4*(n1+n0-2)-1))*(m1-m0)/math.sqrt(pool) if pool>0 else None
        summaries.append({
            "source_study_id":c["Study.ID"],"source_case_id":c["Case.ID"],
            "author_year":c["Author_date"],"country":c["Country"],
            "original_continent":c["Continent"],"bias_rating":c["Bias"],
            "response":c["Response_variable"],"treatment":c["Treatment"],"control":c["Control"],
            "n_treatment":n1,"mean_treatment":m1,"SD_treatment":sd1,
            "n_control":n0,"mean_control":m0,"SD_control":sd0,
            "raw_mean_difference_rail_minus_comparator":m1-m0,
            "descriptive_study_case_Hedges_g_NOT_POOLED":g,
            "original_geographic_metadata_conflict":c["Country"]=="France" and c["Continent"]=="North_America"})
    return {"all_original_meta_cases":205,"railway_cases":3,"unique_original_rail_studies":2,
            "rail_percent":300/205,
            "cases":summaries,"NO_independent_3_study_replication":True,
            "NO_rail_seed_or_management_specific_causal_effect":True}

def penone_audit(accepted,critical):
    doi="10.1016/j.biocon.2012.01.041"
    a=[r for r in accepted if r["DOI"]==doi]
    c=[r for r in critical if r["DOI"]==doi]
    if len(a)!=1 or len(c)!=1:raise ValueError("Original Penone review cohort not uniquely matched")
    r=c[0]
    if r["Risk_of_bias"]!="High" or "paved road" not in r["Other_bias"].lower():
        raise ValueError("Original reviewer Penone bias appraisal changed")
    return {"DOI":doi,"accepted_at_full_text":True,
            "Mell_reviewers_bias_rating":r["Risk_of_bias"],
            "Mell_reviewers_confound":r["Other_bias"],
            "Mell_reviewers_reason":r["Reasons_for_rejection"],
            "excluded_after_critical_appraisal":True,
            "not_a_new_experimental_finding":True}

def run(file5,file8):
    f8=source_tables(file8,SHA8,["File"])
    req={"Study.ID","Case.ID","Author_date","LTI","Country","Continent","Bias",
         "Response_variable","Treatment","Control","N_treatment","N_control",
         "Mean_treatment","Mean_control","Var_treatment","Var_control",
         "Measure_var_treatment","Measure_var_control"}
    cases=as_records(f8["File"],req)
    f5=source_tables(file5,SHA5,["Sheet B - Accepted ful-text",
                                 "Sheet C - Critical appraisal",
                                 "Sheet D - Systematic review table"])
    accepted=as_records(f5["Sheet B - Accepted ful-text"],{"DOI","ID_doc","Publication_title"})
    critical=as_records(f5["Sheet C - Critical appraisal"],
        {"DOI","ID_doc","Risk_of_bias","Other_bias","Reasons_for_rejection"})
    reviewed=as_records(f5["Sheet D - Systematic review table"],
        {"ID_doc","DOI","LTI","Included_in_metaanalysis"})
    if len(reviewed)!=316:raise ValueError("Original review 316 included study rows changed")
    return {"version":"v14.8_verified_original_springer_meta","review_DOI":DOI,
            "appendix5_SHA256":SHA5,"appendix8_SHA256":SHA8,
            "source_316_systematic_review_included_studies":316,
            "rail_meta_cases":rail_cases(cases),
            "Penone_2012_review_risk":penone_audit(accepted,critical),
            "no_train_seed_dispersal_effect":True,
            "no_pooled_railway_meta_effect":True,
            "source_original_France_NorthAmerica_metadata_error_not_silently_fixed":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--study-metadata-5",required=True)
    p.add_argument("--effect-cases-8",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    d=run(a.study_metadata_5,a.effect_cases_8)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(d,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps(d,ensure_ascii=False))

if __name__=="__main__":main()

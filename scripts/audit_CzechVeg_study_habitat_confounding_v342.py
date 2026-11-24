#!/usr/bin/env python3
"""v342: CzechVeg original rail habitat x source-study cross-tabs and canopy audit.

After v341 observed opposite temperature CWM tendencies for rail yards and
embankments, this diagnostic tests whether category is confounded with
original bibliographic study/vegetation physiognomy. It DOES NOT refit EIVE
or output new effects. Original v1.1 headers MD5 source verified, no original
species-by-plot/GPS/site ID or LOCALITY text rehosted.
"""
import argparse,csv,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_CzechVeg_railway_headers_source_v339 import (
    EXPECTED_MD5,ascii_cz,good_number)
from analyze_CzechVeg_rail_nearby_EIVE_source_v341 import (
    load_headers,match_source_outcome_blind)

EXPECTED_MATCHED={"G1_same_source_high_precision":2,
                  "G2_same_source_extended":39,
                  "G4_cross_source_extended":49}
GATES=tuple(EXPECTED_MATCHED)

def small_hash(value):
    return hashlib.sha256(str(value).encode()).hexdigest()[:10]

def cover_fraction(x):
    return x if x is not None and 0<=x<=100 else None

def land_cover_group(habitat):
    h=ascii_cz(habitat)
    if not h or h in ("na","nan","null"):return "unknown_or_unlabelled"
    # Locality/railway labels already filtered by source railway mask.
    if re.search(r"\bles\w*|\bbukov\w*|\bdoubrav\w*|\bborov\w*|\bhabrov\w*",h):
        return "forest_like"
    if re.search(r"\bmokr\w*|\braselin\w*|\brybnik\w*|\brakosin\w*|\bvodn\w*",h):
        return "wetland_like"
    if re.search(r"\blouk\w*|\btravnik\w*|\btravnat\w*|\bpastvin\w*|\bstep\w*",h):
        return "grassland_like"
    if re.search(r"\bruder\w*|\bsterkov\w*|\bskladk\w*|\bcest\w*|\bsilnic\w*|\bparkovist\w*",h):
        return "ruderal_or_transport_like"
    if re.search(r"\bpol\w*|\born\w*|\bpole\b|\bplodin\w*",h):
        return "arable_like"
    return "other_habitat_not_confirmed"

def plot_original_habitat_subset(source,ids):
    found={}
    with Path(source).open("r",encoding="utf-8-sig",newline="") as fh:
        for row in csv.DictReader(fh):
            ident=str(row.get("RELEVE_NR") or "").strip()
            if ident in ids:
                if ident in found:raise ValueError("Plot id repeated in source header")
                found[ident]=str(row.get("HABITAT") or "")
    if set(found)!=set(ids):
        raise ValueError("Original selected reference ID missing source metadata")
    return found

def group_source_pairs(pairs):
    studies=defaultdict(Counter)
    for p in pairs:
        k=small_hash(p["rail"]["reference"])
        studies[k][p["rail_habitat"]]+=1
    study_both=sum(all(s[k]>0 for k in (
        "railway_track_yard_or_structure","railway_embankment_ditch_or_cut"))
        for s in studies.values())
    return {"n_rail_studies_in_matched_pairs":len(studies),
        "n_original_rail_studies_with_both_track_and_embankment":study_both,
        "study_habitat_count_cross_tab_anonymous":{
         key:{"track_yard_or_structure":v["railway_track_yard_or_structure"],
              "embankment_ditch_cut":v["railway_embankment_ditch_or_cut"]}
         for key,v in sorted(studies.items())}}

def numerically_valid_veg_open_pair(p):
    vals=[cover_fraction(q[k]) for q in (p["rail"],p["reference"])
          for k in ("cov_trees","cov_shrubs")]
    if any(x is None for x in vals):
        return False,"missing_original_tree_or_shrub_coverage"
    return (max(vals[0],vals[2])<=10 and max(vals[1],vals[3])<=20),(
        "physiognomy_open_both" if max(vals[0],vals[2])<=10 and max(vals[1],vals[3])<=20
        else "one_or_both_original_plots_too_woody")

def summarize_pairs(pairs,headers_source):
    ids={p["reference"]["id"] for p in pairs}
    habitats=plot_original_habitat_subset(headers_source,ids) if ids else {}
    covercounts=Counter()
    refs=Counter()
    per_source=defaultdict(Counter)
    alt=[]
    tree=[]
    shrub=[]
    open_study=set()
    for p in pairs:
        ref=land_cover_group(habitats[p["reference"]["id"]])
        refs[ref]+=1
        studykey=small_hash(p["rail"]["reference"])
        per_source[studykey][ref]+=1
        rail=p["rail"];near=p["reference"]
        if rail["altitude"] is not None and near["altitude"] is not None:
            alt.append(abs(rail["altitude"]-near["altitude"]))
        tt=[cover_fraction(q["cov_trees"]) for q in (rail,near)]
        ss=[cover_fraction(q["cov_shrubs"]) for q in (rail,near)]
        if None not in tt:tree.append(tt[0]-tt[1])
        if None not in ss:shrub.append(ss[0]-ss[1])
        good,why=numerically_valid_veg_open_pair(p)
        covercounts[why]+=1
        if good:
            covercounts["open_pairs"]+=1
            open_study.add(studykey)
    return {"source_reference_habitat_types":dict(refs),
        "source_reference_habitat_types_by_study_anonymous":{
            k:dict(v) for k,v in sorted(per_source.items())},
        "original_pair_woody_vegetation_coverage_QA":dict(covercounts),
        "n_source_references_with_open_both_plots":len(open_study),
        "source_altitude_difference_abs_median":median(alt) if alt else None,
        "source_altitude_difference_abs_max":max(alt) if alt else None,
        "original_source_tree_cover_rail_minus_nearby_mean_pctpoints":mean(tree) if tree else None,
        "original_source_shrub_cover_rail_minus_nearby_mean_pctpoints":mean(shrub) if shrub else None,
        "source_pairs_with_numeric_altitude_both":len(alt),
        "source_pairs_with_numeric_tree_cover_both":len(tree),
        "source_pairs_with_numeric_shrub_cover_both":len(shrub)}

def audit(source):
    rows,sourceqa=load_headers(source)
    results={}
    for gate,expected in EXPECTED_MATCHED.items():
        pairs,meta=match_source_outcome_blind(rows,gate)
        if len(pairs)!=expected:
            raise ValueError(f"Original v341 outcome-blind matched cohort changed: {gate}")
        results[gate]={
            "original_matched_source_qa":meta,
            **group_source_pairs(pairs),
            "local_reference_habitat_and_tree_cover":summarize_pairs(pairs,source)}
    return {"version":"v342","source_original_full_CzechVeg_MD5":EXPECTED_MD5,
            "source_original_plot_population":sourceqa["original_header_rows"],
            "source_2014_2026_rail_author_studies_independence_gate":results,
            "scientific_gate":[
                "Opposite v341 railway yard and embankment EIVE signs are posthoc and can be confounded by source study unless both habitats occur within same reference study",
                "At most 4 source studies pass strict EIVE numeric trait gate in v341 G2 and G1 has only one eligible plot",
                "Original site vegetation tree and shrub cover can be missing (-1), and missing does NOT mean treeless",
                "Reference HABITAT without railway text is not necessarily actual unaffected or equivalent nonrail plant community",
                "Open-plot conditioning is a potential mediator control, not a randomized habitat exposure",
                "No original plot GPS, botanical source records or raw locality descriptions are republished",
                "No test of railway impact on climate limits, real ballast heat or train-mediated propagule dispersal"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-headers",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    obj=audit(a.source_headers)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(obj,ensure_ascii=False),flush=True)
if __name__=="__main__":main()

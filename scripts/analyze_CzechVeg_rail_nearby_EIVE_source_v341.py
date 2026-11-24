#!/usr/bin/env python3
"""v341 CzechVeg genuine source vegetation species join to independent EIVE.

Source DOI 10.5281/zenodo.19231152: exact v1.1 full 117739 plot header
and ~383MB original vascular-plant species rows. Source-verified MD5s.
HABITAT-derived explicit railway-substrate candidates only, matched to
nearby plots with NO rail text in HABITAT or LOCALITY. These controls are
UNLABELLED nearby vegetation, *not* proved unaffected nonrail habitats.

All same-source/cross-source nearest one-to-one geographic/area/year
matching is outcome-blind and fixed from prior v340 matching feasibility.
Source plot row identities and GPS never exported. Railway EIVE thermal
preference is not measured physical microclimate or seed carriage.
"""
import argparse,csv,hashlib,itertools,json,math,re,unicodedata
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

import numpy as np
from scipy.optimize import linear_sum_assignment

from audit_CzechVeg_railway_headers_source_v339 import (
    EXPECTED_MD5,EXPECTED_PLOTS,ascii_cz,classify,good_number)
from audit_CzechVeg_actual_railways_nearby_controls_v340 import (
    CERTIFIED,MATCH_GATES,rail_certainty,src_metadata,in_bounds,distance_km)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table,normal_name

SOURCE_SPECIES_MD5="b296e5e19f03fb3bad0a8347c317f8a3"
MIN_EXPECTED_SOURCE_SPECIES_ROWS=1000000
SOURCE_DOI="10.5281/zenodo.19231152"
EIVE_DOI="10.5281/zenodo.7534792"
MAIN_GATES=("G1_same_source_high_precision","G2_same_source_extended",
            "G4_cross_source_extended")
PANELS=((3,5,0.5),(1,3,0.3))
MEASURES=("position","width")
HABITATS=("railway_track_yard_or_structure","railway_embankment_ditch_or_cut")

def safe_md5(path,expected):
    obj=hashlib.md5()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):
            obj.update(chunk)
    if obj.hexdigest()!=expected:
        raise ValueError(f"Not original checksum-protected study file: {Path(path).name}")
    return obj.hexdigest()

def load_headers(path):
    safe_md5(path,EXPECTED_MD5)
    rows=[]
    cnt=Counter()
    with Path(path).open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        needed={"RELEVE_NR","HABITAT","LOCALITY","REFERENCE","YEAR","SURF_AREA",
                "PRECISION","DEG_LAT","DEG_LON"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError("Publisher original header schema mismatch")
        seen=set()
        for data in reader:
            cnt["original_header_rows"]+=1
            ident=str(data["RELEVE_NR"]).strip()
            if not ident or ident in seen:
                raise ValueError("Source original vegetation plot ID missing/duplicated")
            seen.add(ident)
            typ,_=classify(data["HABITAT"],data["LOCALITY"])
            cat=(rail_certainty(data["HABITAT"]) if typ=="rail_habitat_text_candidate"
                 else "not_certified_rail_habitat")
            flag="rail" if cat in CERTIFIED else (
                 "potential_reference" if typ=="no_rail_text_evidence" else "excluded")
            if flag=="excluded":continue
            meta=src_metadata(data)
            meta["id"]=ident;meta["habitat_type"]=cat
            meta["role"]=flag
            if flag=="rail":
                cnt["certified_rail_candidates"]+=1
            else:
                cnt["unlabelled_potential_references"]+=1
            rows.append(meta)
    if cnt["original_header_rows"]!=EXPECTED_PLOTS:
        raise ValueError("CzechVeg original population 117739 plot drift")
    if cnt["certified_rail_candidates"]!=84:
        raise ValueError("v340 certified rail candidate population changed")
    return rows,dict(cnt)

def index_nonrail(rows):
    grid=defaultdict(list)
    for x in rows:
        if x["role"]!="potential_reference" or x["lat"] is None:continue
        key=(int(x["lat"]/.02),int(x["lon"]/.02))
        grid[key].append(x)
    return grid

def possible_controls(rail,grid,gate):
    if any(rail[k] is None for k in ("lat","lon","year","area","precision")):
        return []
    if rail["precision"]>gate["max_precision_m"]:return []
    i=int(rail["lat"]/.02);j=int(rail["lon"]/.02)
    out=[]
    for ii in range(i-2,i+3):
        for jj in range(j-3,j+4):
            for candidate in grid.get((ii,jj),[]):
                valid,dist=in_bounds(rail,candidate,gate)
                if not valid:continue
                if candidate["id"]==rail["id"]:
                    raise AssertionError("Identical source plot used as its own control")
                cost=dist+0.1*abs(rail["year"]-candidate["year"])+0.25*abs(
                    math.log(rail["area"]/candidate["area"]))
                out.append((candidate,dist,cost))
    return out

def match_source_outcome_blind(rows,gate_name):
    gate=MATCH_GATES[gate_name]
    grid=index_nonrail(rows)
    rail=sorted((r for r in rows if r["role"]=="rail"),key=lambda x:x["id"])
    # For cross-reference matching, explicitly do not count the references
    # as same-study causal controls.
    groups=defaultdict(list)
    for r in rail:
        group=r["reference"] if gate["same_reference"] else "ALL_STUDIES"
        groups[group].append(r)
    all_pairs=[]
    source_counts=Counter()
    feasible_edges=0
    for ref,rail_group in sorted(groups.items()):
        # Missing source REF is not a valid exact-study comparison.
        if gate["same_reference"] and not ref:continue
        candidates_by_rail=[]
        control_pool={}
        for r in rail_group:
            current=possible_controls(r,grid,gate)
            candidates_by_rail.append(current)
            for c,d,cost in current:
                control_pool[c["id"]]=c
                feasible_edges+=1
        controlids=sorted(control_pool)
        n,m=len(rail_group),len(controlids)
        if m==0:continue
        lookup={x:i for i,x in enumerate(controlids)}
        values=np.full((n,m+n),1e8,dtype=float)
        values[:,m:]=1e6  # 1 dummy unmatched slot per rail plot
        pair_info={}
        for i,poss in enumerate(candidates_by_rail):
            for candidate,dist,cost in poss:
                j=lookup[candidate["id"]]
                values[i,j]=cost+1e-8*(i*m+j)
                pair_info[i,j]=(dist,cost)
        ix,jx=linear_sum_assignment(values)
        for i,j in zip(ix,jx):
            if j>=m or values[i,j]>=1e6:continue
            a=rail_group[i];b=control_pool[controlids[j]]
            dist,cost=pair_info[i,j]
            all_pairs.append({"rail":a,"reference":b,
                "distance_km":dist,"year_diff":abs(a["year"]-b["year"]),
                "area_ratio":max(a["area"],b["area"])/min(a["area"],b["area"]),
                "same_study":a["reference"]==b["reference"] and bool(a["reference"]),
                "rail_habitat":a["habitat_type"]})
            source_counts[a["reference"]]+=1
    if len({x["rail"]["id"] for x in all_pairs})!=len(all_pairs):
        raise AssertionError("Rail sample reuse in 1:1 matching")
    if len({x["reference"]["id"] for x in all_pairs})!=len(all_pairs):
        raise AssertionError("Nonrail candidate control reuse in 1:1 matching")
    result={"matching_source_gate":gate_name,
            "n_original_certified_railway_substrate_candidates":len(rail),
            "n_1to1_matched_source_plot_pairs":len(all_pairs),
            "n_explicit_railway_study_REFERENCES_in_matched_sample":len(source_counts),
            "source_matched_habitat_counts":dict(Counter(x["rail_habitat"] for x in all_pairs)),
            "same_study_matched_pairs":sum(x["same_study"] for x in all_pairs),
            "all_potential_feasible_geographic_edges":feasible_edges,
            "n_distinct_original_source_control_studies":len(set(x["reference"]["reference"] for x in all_pairs)),
            "median_match_distance_km":median(x["distance_km"] for x in all_pairs) if all_pairs else None,
            "median_original_source_year_difference":median(x["year_diff"] for x in all_pairs) if all_pairs else None,
            "median_original_plot_area_ratio":median(x["area_ratio"] for x in all_pairs) if all_pairs else None}
    return all_pairs,result

def source_exact_species(raw):
    if not raw:return None
    value=normal_name(raw)
    if re.fullmatch(r"[a-z][a-z-]+ [a-z][a-z-]+",value):
        return value
    return None

def load_selected_species(path,plot_ids):
    safe_md5(path,SOURCE_SPECIES_MD5)
    plants=defaultdict(dict)
    checks=Counter()
    with Path(path).open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        needed={"RELEVE_NR","KAPLAN_PLUS_NAME_SIMPLE","PLANT_GROUP","COVER_PERC"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError(f"Original Czech species table lacks columns: {needed-set(reader.fieldnames or [])}")
        for r in reader:
            checks["original_source_taxon_rows"]+=1
            plot=str(r["RELEVE_NR"]).strip()
            if plot not in plot_ids:continue
            checks["target_source_species_rows"]+=1
            if str(r["PLANT_GROUP"] or "").strip().casefold()!="vascular":
                checks["nonvascular_records_excluded"]+=1
                continue
            label=normal_name(r["KAPLAN_PLUS_NAME_SIMPLE"])
            if not label:raise ValueError("Original vascular source row lacks taxonomic name")
            taxon=source_exact_species(label)
            key=taxon if taxon else "unresolved_source_taxon:"+label
            cov=good_number(r["COVER_PERC"])
            if cov is not None and (cov<0 or cov>100):
                checks["source_cover_out_of_range"]+=1
                cov=None
            if cov is not None and cov==0:
                cov=None
            prior=plants[plot].get(key)
            if prior is not None:
                checks["duplicate_original_accepted_taxon_layers"]+=1
                if cov is not None:
                    plants[plot][key]=max(prior,cov) if prior is not None else cov
            else:
                plants[plot][key]=cov
            checks["source_selected_vascular_species_records"]+=1
    if checks["original_source_taxon_rows"]<MIN_EXPECTED_SOURCE_SPECIES_ROWS:
        raise ValueError("Source species table row count unexpectedly low")
    if len(plants)<len(plot_ids)*.8:
        raise ValueError("Most matched original vegetation plots absent from source species table")
    return plants,dict(checks)

def compute_community(taxa,eive,metric,systems):
    if not taxa:return {"source_species_n":0,"trait_n":0,"trait_coverage":None,
                        "equal_species_cwm":None,"cover_weighted_cwm":None,
                        "trait_source_cover_nonmissing_n":0}
    eligible={k:(eive[k][metric],cover) for k,cover in taxa.items()
              if k in eive and eive[k]["systems"]>=systems}
    valid_cover=[(v[0],v[1]) for v in eligible.values() if v[1] is not None]
    cov_cwm=(sum(t*c for t,c in valid_cover)/sum(c for t,c in valid_cover)
             if valid_cover else None)
    return {"source_species_n":len(taxa),"trait_n":len(eligible),
            "trait_coverage":len(eligible)/len(taxa),
            "equal_species_cwm":mean(x[0] for x in eligible.values()) if eligible else None,
            "cover_weighted_cwm":cov_cwm,
            "trait_source_cover_nonmissing_n":len(valid_cover),
            "source_cover_known_trait_fraction":len(valid_cover)/len(eligible) if eligible else None}

def panel_stats(pairs,plants,eive,systems,minn,coverage,metric):
    study=defaultdict(list)
    by_habitat=defaultdict(list)
    qa=Counter()
    coverage_measure=[]
    for p in pairs:
        a=compute_community(plants.get(p["rail"]["id"]),eive,metric,systems)
        b=compute_community(plants.get(p["reference"]["id"]),eive,metric,systems)
        qa["source_geo_matched_pairs"]+=1
        if a["trait_coverage"] is None or b["trait_coverage"] is None:
            qa["missing_matched_plots_or_empty_vascular_flora"]+=1;continue
        coverage_measure.append((a["trait_coverage"],b["trait_coverage"]))
        if min(a["trait_n"],b["trait_n"])<minn or min(
                    a["trait_coverage"],b["trait_coverage"])<coverage:
            qa["pairs_below_EIVE_numeric_trait_gate"]+=1;continue
        qa["eligible_source_vegetation_pairs"]+=1
        difference=a["equal_species_cwm"]-b["equal_species_cwm"]
        reference=p["rail"]["reference"]
        study[reference].append(difference)
        by_habitat[p["rail_habitat"]].append(difference)
        if a["cover_weighted_cwm"] is not None and b["cover_weighted_cwm"] is not None:
            if min(a["trait_source_cover_nonmissing_n"],b["trait_source_cover_nonmissing_n"])>=minn and min(
                a["source_cover_known_trait_fraction"],b["source_cover_known_trait_fraction"])>=.8:
                by_habitat[p["rail_habitat"]+"_cover_eligible"].append(
                    a["cover_weighted_cwm"]-b["cover_weighted_cwm"])
                qa["eligible_80pct_source_cover_pairs"]+=1
    vals=[v for arr in study.values() for v in arr]
    studymeans=[mean(v) for v in study.values()]
    byhab={k:{"n_source_pairs":len(v),"site_equal_mean":mean(v),
              "positive_pairs":sum(x>0 for x in v)}
           for k,v in sorted(by_habitat.items()) if v}
    directions={}
    for k,v in study.items():
        token=hashlib.sha1(str(k).encode()).hexdigest()[:10]
        directions[token]={"n_source_pairs":len(v),
                           "study_mean":mean(v),"positive_source_pairs":sum(x>0 for x in v)}
    if 1<=len(studymeans)<=12:
        obs=abs(mean(studymeans))
        samples=[abs(sum(s*v for s,v in zip(studymeans,sign))/len(studymeans))
                 for sign in itertools.product((-1,1),repeat=len(studymeans))]
        pvalue=sum(x>=obs-1e-12 for x in samples)/len(samples)
    else:pvalue=None
    return {"metric":metric,"minimum_EIVE_contributing_systems":systems,
            "min_numeric_species_per_habitat":minn,
            "min_trait_coverage_each_plot":coverage,
            "source_qa":dict(qa),
            "all_source_plot_rail_reference_numeric_trait_coverage_means":{
                "rail":mean(x[0] for x in coverage_measure) if coverage_measure else None,
                "reference":mean(x[1] for x in coverage_measure) if coverage_measure else None},
            "n_original_studies_with_eligible_pairs":len(study),
            "n_original_field_plot_pairs_with_EIVE":len(vals),
            "site_equal_mean_rail_minus_nearby_unlabelled":mean(vals) if vals else None,
            "study_equal_mean_rail_minus_nearby_unlabelled":mean(studymeans) if studymeans else None,
            "study_mean_positive_count":sum(x>0 for x in studymeans),
            "study_mean_negative_count":sum(x<0 for x in studymeans),
            "study_sign_flip_exploratory_two_sided_fraction":pvalue,
            "source_study_aggregates_anonymous":directions,
            "rail_habitat_label_stratification":byhab}

def analyze(headers,species,eivefile):
    rows,hq=load_headers(headers)
    matches={}
    all_ids=set()
    for gate in MAIN_GATES:
        paired,qa=match_source_outcome_blind(rows,gate)
        matches[gate]=(paired,qa)
        for p in paired:
            all_ids.add(p["rail"]["id"]);all_ids.add(p["reference"]["id"])
    if not all_ids:
        raise ValueError("Original CzechVeg source has no feasible 1:1 railway-reference plots")
    plant_data,plantqa=load_selected_species(species,all_ids)
    eive,eiveqa=eive_table(eivefile)
    panels=[]
    for gate in MAIN_GATES:
        pairs,matchqa=matches[gate]
        obj={"gate":gate,"pretrait_matching_qa":matchqa,"panels":[]}
        for systems,minn,cover in PANELS:
            for metric in MEASURES:
                obj["panels"].append(panel_stats(pairs,plant_data,eive,systems,minn,cover,metric))
        panels.append(obj)
    return {"version":"v341","source":"genuine CzechVeg-Open entire 117739 vegetation plot study",
        "original_2026_03_26_headers_md5":EXPECTED_MD5,
        "original_2026_03_26_species_md5":SOURCE_SPECIES_MD5,
        "independent_EIVE_1_0_source_md5":eiveqa["EIVE_original_md5"],
        "original_source_header_population":hq,
        "original_selected_species_source_qa":plantqa,
        "selected_unique_original_source_plots_for_species_join":len(all_ids),
        "results":panels,
        "cautions":[
           "Czech railway samples are selected from HABITAT descriptions only, not direct field verified active ballast tracks",
           "UNLABELLED nearby controls are not experimentally demonstrated nonrail habitats or unaffected community counterfactuals",
           "The matched observations often have up to 2km physical separation, 7 years difference and coordinate uncertainty up to 1km",
           "Different original vegetation studies and habitat age/manageability may explain selection, not train operating ecology",
           "Same-study reference equal weighting does not make studies independent climate regions",
           "EIVE matched Europe taxa temperature indicator, not measured ballast or regional air temperature",
           "CzechVeg original taxon cover percent converted from multiple ordinal scales; treat as estimated relative cover, not directly observed biomass",
           "No original source location, plot ID, species-by-plot source rows or actual GPS provided in derived outputs",
           "Do not call this a climate-offset railway effect or a meta-analysis of genuinely paired studies"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-headers",required=True)
    p.add_argument("--source-species",required=True)
    p.add_argument("--eive-xlsx",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=analyze(a.source_headers,a.source_species,a.eive_xlsx)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"version":"v341",
        "selected_unique_original_source_plots_for_species_join":data["selected_unique_original_source_plots_for_species_join"],
        "original_selected_species_source_qa":data["original_selected_species_source_qa"],
        "results":[{"gate":x["gate"],"pretrait_matching_qa":x["pretrait_matching_qa"],
                    "panels":[{k:p[k] for k in (
                        "metric","minimum_EIVE_contributing_systems",
                        "min_numeric_species_per_habitat",
                        "n_original_studies_with_eligible_pairs",
                        "n_original_field_plot_pairs_with_EIVE",
                        "site_equal_mean_rail_minus_nearby_unlabelled",
                        "study_equal_mean_rail_minus_nearby_unlabelled",
                        "study_mean_positive_count",
                        "study_mean_negative_count",
                        "study_sign_flip_exploratory_two_sided_fraction",
                        "source_qa","rail_habitat_label_stratification")}
                        for p in x["panels"]]} for x in data["results"]]},
        ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()

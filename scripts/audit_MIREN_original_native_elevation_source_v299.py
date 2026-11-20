#!/usr/bin/env python3
"""Source-only gate for ORIGINAL MIREN Elevation metadata (metres).
Read original 2007-2023 Zenodo CSV, validate original Plot1/Plot3 altitude,
and audit true first-eligible-year MIREN field cohort, before biological test.
Do not release original site-level altitudes, GPS, site/species identities.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
 read_source_miren,traits_from_source,MIREN_MD5)
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 import PANEL_DEFINITIONS

MAX_WITHIN_PLOT_ELEVATION_SPREAD_M=5.
MAX_PAIRED_ROADSIDE_INLAND_ALTITUDE_DIFFERENCE_M=100.
MIN_NUMERIC_ELEVATION_FRACTION=.90
MIN_ORIGINAL_ROAD_ELEVATION_SPAN_M=150.
MIN_ORIGINAL_ROAD_SITES=4
MIN_ORIGINAL_REGION_ROADS=2
MIN_ORIGINAL_SITE_RETENTION=.80

def read_original_elevation_csv(src:Path,original_records:dict):
    keys=set(original_records)
    raw=defaultdict(list)
    pids=defaultdict(set)
    errors=Counter()
    with src.open(encoding="utf-8-sig",newline="") as f:
        it=csv.DictReader(f)
        fields={"Region","Road","Site","Year","Plot","Plot_id","Elevation"}
        if not fields.issubset(it.fieldnames or []):
            raise ValueError("Genuine MIREN original author Elevation field missing")
        for row in it:
            year=row["Year"].strip()
            plot=row["Plot"].strip()
            m=re.fullmatch(r"([123])([A-Za-z]?)",plot)
            if not m or not year.isdigit():continue
            key=(row["Region"].strip(),row["Road"].strip(),
                 row["Site"].strip(),int(year),int(m.group(1)))
            if key not in keys:continue
            if not row["Plot_id"].strip().endswith("."+plot):
                errors["original_plot_id_suffix_invalid"]+=1
                continue
            pids[key].add(row["Plot_id"].strip())
            try:
                e=float(row["Elevation"].strip())
                if not (math.isfinite(e) and -500<=e<=9000):
                    raise ValueError("Original terrestrial altitude invalid")
            except (ValueError,TypeError):
                e=None
                errors["original_plant_row_invalid_elevation"]+=1
            raw[key].append(e)
    meta={}
    for key,vals in raw.items():
        good=[x for x in vals if x is not None]
        frac=len(good)/len(vals)
        spread=max(good)-min(good) if good else None
        okay=(len(pids[key])==1 and
              frac>=MIN_NUMERIC_ELEVATION_FRACTION and
              spread is not None and spread<=MAX_WITHIN_PLOT_ELEVATION_SPREAD_M)
        meta[key]=float(median(good)) if okay else None
    qa={
     "n_source_original_plots_with_elevation_metadata":len(raw),
     "n_source_plant_rows_with_author_plot_elevation":sum(map(len,raw.values())),
     "n_source_plots_with_qualified_original_Elevation":sum(x is not None for x in meta.values()),
     "n_source_plots_without_usable_author_elevation":sum(x is None for x in meta.values()),
     "original_invalid_source_row_QA":dict(errors),
     "original_site_specific_elevation_GPS_and_species_NEVER_exported":True}
    return meta,qa

def original_paired_site_elevation(row:dict,region:str,road:str,meta:dict):
    site=row["site"];year=int(row["year"])
    a=meta.get((region,road,site,year,1))
    b=meta.get((region,road,site,year,3))
    if a is None or b is None:
        return None
    if abs(a-b)>MAX_PAIRED_ROADSIDE_INLAND_ALTITUDE_DIFFERENCE_M:
        return None
    return (a+b)/2

def evaluate_original_region_support(roads:list,region:str,meta:dict):
    n_original=sum(len(sites) for _,sites in roads)
    byroad=[]
    for road,sites in roads:
        matched=[(row,original_paired_site_elevation(row,region,road,meta))
                 for row in sites]
        ok=[(row,e) for row,e in matched if e is not None]
        elev=[e for row,e in ok]
        span=(max(elev)-min(elev)) if len(elev)>=2 else 0.
        eligible=len(ok)>=MIN_ORIGINAL_ROAD_SITES and span>=MIN_ORIGINAL_ROAD_ELEVATION_SPAN_M
        byroad.append({
          "source_n_original_paired_sites":len(sites),
          "source_n_valid_paired_original_Elevation_sites":len(ok),
          "original_source_road_elevation_range_m":float(span),
          "source_n_survey_years_among_original_sites":len({int(row["year"]) for row in sites}),
          "source_author_road_elevation_slope_identifiable":bool(eligible)})
    supported=[x for x in byroad if x["source_author_road_elevation_slope_identifiable"]]
    n_sites=sum(x["source_n_valid_paired_original_Elevation_sites"] for x in supported)
    qualified=(len(supported)>=MIN_ORIGINAL_REGION_ROADS and n_original>0 and
               n_sites/n_original>=MIN_ORIGINAL_SITE_RETENTION)
    return {
       "source_original_roads":len(roads),"source_original_paired_sites":n_original,
       "source_road_n_with_ge4_sites_and_ge150m_elevation":len(supported),
       "source_n_paired_sites_retained_with_valid_elevation":n_sites,
       "source_site_retention_fraction":n_sites/n_original if n_original else 0,
       "source_road_n_with_mixed_original_first_eligible_survey_years":
         sum(x["source_n_survey_years_among_original_sites"]>1 for x in byroad),
       "source_environmental_gradient_support_PASS":bool(qualified),
       "by_road_anonymous_original_metadata_support":byroad}

def run(src:Path,height:Path,seed:Path,output:Path):
    if hashlib.md5(src.read_bytes()).hexdigest()!=MIREN_MD5:
        raise ValueError("Original MIREN author CSV MD5 drift")
    for p,key in ((height,"height"),(seed,"seed")):
        if hashlib.sha256(p.read_bytes()).hexdigest()!=LEDASHA[key]:
            raise ValueError("Original Oldenburg LE DA SHA drift")
    original,author_qa=read_source_miren(src)
    h,_=read_leda(height,"height_m");s,_=read_leda(seed,"seed_mass_mg")
    traits,trait_qa=traits_from_source(original,h,s)
    if len(traits)!=1096 or author_qa["original_source_record_QA"][
        "original_source_species_observation_rows"]!=167804:
        raise ValueError("Original source MIREN flora/source taxonomy support drift")
    meta,meta_qa=read_original_elevation_csv(src,original)
    panels={}
    for panel,(minimum,fraction,year,nregion,nroad,nsite) in PANEL_DEFINITIONS.items():
        regions,source_qa=valid_region_panels(original,traits,minimum,fraction,year)
        actual=(len(regions),sum(len(a) for a in regions.values()),
                sum(len(sites) for arr in regions.values() for _,sites in arr))
        if actual!=(nregion,nroad,nsite):
            raise ValueError("Original true first-year road/inland source cohort drift")
        r={reg:evaluate_original_region_support(roads,reg,meta)
           for reg,roads in sorted(regions.items())}
        panels[panel]={
          "n_original_regions":nregion,"n_original_road_samples":nroad,
          "n_original_road_inland_same_site_pairs":nsite,
          "source_original_region_code_elevation_support":r,
          "n_original_region_codes_with_valid_environmental_gradients":
            sum(q["source_environmental_gradient_support_PASS"] for q in r.values())}
    out={"version":"v29.9_SOURCE_ONLY_ELEVATION_QA",
         "source_original_dataset":"10.5281/zenodo.15749536",
         "source_original_pub":"10.1111/geb.70284",
         "source_species_rows":167804,
         "source_original_plant_Elevation_data":meta_qa,
         "source_true_MIREN_first_survey_field_support":panels,
         "source_analysed_field_type":"ROAD_AND_SAME_SITE_INLAND_NOT_RAIL",
         "new_community_stature_environmental_effect_computed":False,
         "original_per_site_elevation_and_GPS_NOT_exported":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"source_original_elevation_metadata":meta_qa,
      "original_source_elevation_supported_regions":{
        k:{"eligible":v["n_original_region_codes_with_valid_environmental_gradients"],
           "total":v["n_original_regions"],
           "region_codes_and_counts":{reg:{
             "original_source_elevation_valid_roads":q["source_road_n_with_ge4_sites_and_ge150m_elevation"],
             "original_sites_retained":q["source_n_paired_sites_retained_with_valid_elevation"],
             "source_survey_year_mixed_roads":q["source_road_n_with_mixed_original_first_eligible_survey_years"],
             "source_elevation_support":q["source_environmental_gradient_support_PASS"]}
             for reg,q in v["source_original_region_code_elevation_support"].items()}}
             for k,v in panels.items()}},indent=2))
    return out
def main():
    p=argparse.ArgumentParser()
    for name in ("miren","height","seed","output"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    run(a.miren,a.height,a.seed,a.output)
if __name__=="__main__":main()

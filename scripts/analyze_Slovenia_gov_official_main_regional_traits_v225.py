#!/usr/bin/env python3
"""v22.5 source-verified SI railway main(G) vs regional(R), robustness of v22.4.

Official 2013 decree G/R line designation copied from Slovenian government
2024 network national emergency plan (July 2024 update), p. 7:
https://www.gov.si/assets/organi-v-sestavi/URSZR/Datoteke/Drzavni-nacrti/dn_zelezniska_azuriran_2024.pdf
Original vegetation author 2022/2023 site field 'track_section_number' is
national line ID as found in original source. No OSM dependence.
These are administrative classes, NOT measured passenger/freight volumes
or wind, mowing pressure, station departure frequency or true adaptation.

Primary comparison is source 89 original sites, annual BS-BT, and source
82 original sites meeting joint species-median LEDA >=5 per BS and BT;
line-level effects G(6 lines) vs R(9 lines) are descriptive and
confounded by climate, geographic region, substrate, field site selection.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timedelta
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx,validate_original_keys,pair_source,locations_independent
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds,jaccard
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import original_gps_coordinates,haversine_km
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    HEIGHT_SHA,SEED_SHA,audit_si_source,all_original_matched_taxa,ztransform,
    create_sites,site_effects,section_estimates,DIMENSIONS
)

OFFICIAL_SOURCE_URL="https://www.gov.si/assets/organi-v-sestavi/URSZR/Datoteke/Drzavni-nacrti/dn_zelezniska_azuriran_2024.pdf"
# Official 2024 plan, page 7/36, table 1: G=glavne proge, R=regionalne proge.
MAIN=frozenset(("10","20","30","40","41","44","50","60","62","64"))
REGIONAL=frozenset(("21","31","32","33","34","35","42","43","61","70","71","72","73","80","81","82"))
KNOWN=("10","20","21","31","32","34","41","42","44","50","64","70","72","80","81")
CUTS=(0,5,10,20,40,80,160,500)

def official_class(section):
    section=str(section)
    if section in MAIN:return "main_G"
    if section in REGIONAL:return "regional_R"
    return None

def aggregate(values):
    return {"n":len(values),"mean":mean(values) if values else None,
            "median":median(values) if values else None,
            "min":min(values) if values else None,
            "max":max(values) if values else None}

def mean_by_section(records,metric):
    by=defaultdict(list)
    for r in records:
        val=r.get(metric)
        if val is not None and math.isfinite(val):
            by[r["section"]].append(float(val))
    return {s:mean(v) for s,v in by.items() if v}

def categorical_contrast(source_section_means):
    g=[v for k,v in source_section_means.items() if official_class(k)=="main_G"]
    r=[v for k,v in source_section_means.items() if official_class(k)=="regional_R"]
    out={"G_main_section_n":len(g),"R_regional_section_n":len(r),
         "G_main_section_equal_mean":mean(g) if g else None,
         "R_regional_section_equal_mean":mean(r) if r else None,
         "main_minus_regional":mean(g)-mean(r) if g and r else None}
    if len(g)>=3 and len(r)>=3:
        leave=[]
        for s in sorted(source_section_means):
            keep={k:v for k,v in source_section_means.items() if k!=s}
            gg=[v for k,v in keep.items() if official_class(k)=="main_G"]
            rr=[v for k,v in keep.items() if official_class(k)=="regional_R"]
            leave.append({"leftout_original_line_code":s,"G_minus_R":mean(gg)-mean(rr)})
        out["leave_one_section_out_min_max"]=[min(x["G_minus_R"] for x in leave),
                                             max(x["G_minus_R"] for x in leave)]
        out["leave_one_section_out_sign_positive_n"]=sum(x["G_minus_R"]>0 for x in leave)
        out["leave_one_section_out_n"]=len(leave)
    return out

def official_summary(original,raw_height,raw_seed):
    loc,zone,sp,traits,cb,quality=validate_original_keys(original)
    pairs=pair_source(loc,zone,sp,traits)
    full=locations_independent(pairs)
    observed={str(x["section_number"]) for x in full}
    if observed!=set(KNOWN) or len(full)!=89 or len(observed)!=15:
        raise ValueError("National original 15 line codes do not match published main/regional registry")
    for row in full:
        row["section"]=str(row["section_number"])
        row["class"]=official_class(row["section"])
        if row["class"] is None:raise ValueError("Unknown government class for original source railway")
    attrs={
       "strict_annual_BS_minus_BT":"BS_minus_BT_annual_strict",
       "strict_therophyte_BS_minus_BT":"BS_minus_BT_therophyte_strict",
       "original_cover_BS_minus_BT":"BS_minus_BT_cover",
       "original_richness_BS_minus_BT":"BS_minus_BT_richness"
    }
    source_full={}
    for label,key in attrs.items():
        sub=[{**x,"effect":x.get(key)} for x in full]
        per=mean_by_section(sub,"effect")
        source_full[label]={
           "all_locations_89_site_equal":aggregate([x["effect"] for x in sub if x["effect"] is not None]),
           "main_vs_regional_official_line_equal":categorical_contrast(per)}
    counts={}
    for typ in ("main_G","regional_R"):
        subset=[x for x in full if x["class"]==typ]
        counts[typ]={"n_source_original_field_sites":len(subset),
                     "n_source_original_lines":len({x["section"] for x in subset}),
                     "line_code_site_counts":dict(Counter(x["section"] for x in subset)),
                     "track_type_in_source":dict(Counter(str(x["track_type"]) for x in subset))}
    qc=audit_si_source(original,raw_height,raw_seed)
    names,codebook=original_slovenian_EPPO_codebook_binomial_mapping(original)
    taxa=all_original_matched_taxa(zone,sp)
    exact={t:binomial(names.get(t)) for t in taxa}
    values,transform=ztransform(exact,raw_height,raw_seed)
    if len(values)!=224:raise ValueError("v22.4 known joint taxa gate")
    eligible,source_pairs=create_sites(original,values)
    scored=site_effects(eligible,values)
    section,overall=section_estimates(scored)
    original_locations={str(x["rail_location_id"]):x for x in full}
    for site,s in scored.items():
        s["class"]=official_class(s["section"])
        src=original_locations[site]
        s["annual_delta"]=src["BS_minus_BT_annual_strict"]
        s["therophyte_delta"]=src["BS_minus_BT_therophyte_strict"]
        s["original_cover_high_value_flag"]=bool(src["has_any_original_cover_above100"])
        stamp=loc[src["rail_location_id"]]["date"]
        s["year"]=(datetime(1899,12,30)+timedelta(days=float(stamp))).year
        s["track_type"]=str(src["track_type"])
    contrast={}
    for m in DIMENSIONS:
        per={s:r[m] for s,r in section.items()}
        contrast[m]=categorical_contrast(per)
        contrast[m]["section_equal_all_original"]=overall[m]
        contrast[m]["site_equal_all_original"]=mean(x["delta"][m] for x in scored.values())
        contrast[m]["n_negative_sections"]=sum(v<0 for v in per.values())
        contrast[m]["source_lines_with_at_least_3_eligible_sites"]={}
        th={s:v[m] for s,v in section.items() if sum(r["section"]==s for r in scored.values())>=3}
        contrast[m]["source_lines_with_at_least_3_eligible_sites"]={
           "source_sections":len(th),"G_R":categorical_contrast(th)}
    # Strict matching >=10 known species on each source side.
    elig10={k:[r for r in v if r["nBS"]>=10 and r["nBT"]>=10]
            for k,v in eligible.items()}
    elig10={k:v for k,v in elig10.items() if v}
    e10=site_effects(elig10,values)
    sec10,mean10=section_estimates(e10)
    checks={"paired_at_least_10_each_zone":{
             "real_same_side_pairs":sum(map(len,elig10.values())),
             "n_sites":len(e10),"n_sections":len(sec10),
             "joint_trace":mean10.get("joint_trace"),"height_var":mean10.get("height_var"),
             "seed_var":mean10.get("seed_var")},
       "exclude_any_unusual_original_cover_gt100":{
             "source_sites_excluded":sum(x["original_cover_high_value_flag"] for x in scored.values())}}
    no_bad={site:x for site,x in scored.items() if not x["original_cover_high_value_flag"]}
    sec_ok,mean_ok=section_estimates(no_bad)
    checks["exclude_any_unusual_original_cover_gt100"].update({
        "n_sites":len(no_bad),"n_sections":len(sec_ok),"joint_trace":mean_ok["joint_trace"],
        "height_var":mean_ok["height_var"],"seed_var":mean_ok["seed_var"]})
    checks["by_original_year"]={}
    for year in (2022,2023):
        sub={site:x for site,x in scored.items() if x["year"]==year}
        secs,avg=section_estimates(sub)
        checks["by_original_year"][str(year)]={"n_sites":len(sub),
              "n_source_lines":len(secs),"section_equal":avg}
    checks["by_original_track_type"]={}
    for t in sorted({x["track_type"] for x in scored.values()}):
        sub={site:x for site,x in scored.items() if x["track_type"]==t}
        secs,avg=section_estimates(sub)
        checks["by_original_track_type"][t]={"n_sites":len(sub),
              "n_source_lines":len(secs),"section_equal":avg}
    checks["leave_one_original_section_out"]={}
    for metric in ("joint_trace","height_var","seed_var"):
        values_per={s:v[metric] for s,v in section.items()}
        remaining=[mean(v for ss,v in values_per.items() if ss!=s) for s in sorted(values_per)]
        checks["leave_one_original_section_out"][metric]={
          "min":min(remaining),"max":max(remaining),
          "positive_n":sum(v>0 for v in remaining),"n":len(remaining)}
    # Original fixed geo: all pairs scored as descriptive only, do not
    # report naive pairwise p-values. Exclude same original source lines.
    source_sites,src_guild,meta=source_sites_and_guilds(original)
    attached,gpsqa=original_gps_coordinates(original,source_sites)
    pair=defaultdict(list)
    for a,b in combinations(attached,2):
        if a["section"]==b["section"]:continue
        ca=official_class(a["section"]);cb=official_class(b["section"])
        group="G_G" if ca==cb=="main_G" else ("R_R" if ca==cb=="regional_R" else "G_R")
        dist=haversine_km(a["lat"],a["lon"],b["lat"],b["lon"])
        band=next((f"{lo}-{hi}" for lo,hi in zip(CUTS[:-1],CUTS[1:]) if lo<=dist<hi),None)
        if band is None:continue
        sim=jaccard(a["known"],b["known"])
        pair[(group,band)].append((sim,dist))
    geography={}
    for (g,band),v in sorted(pair.items()):
        geography.setdefault(g,{})[band]={"dependent_site_pairs":len(v),
           "mean_known_species_Jaccard":mean(x for x,d in v),
           "mean_greatcircle_km":mean(d for x,d in v)}
    # Unmatched line pairs can be geographically very differently distributed.
    return {
      "status":"REAL_SI_ORIGINAL_GOVERNMENT_MAIN_VS_REGIONAL_DESCRIPTIVE",
      "source_2026_author_rail_vegetation_xlsx_sha":meta["original_source_supplement_SHA256"],
      "source_main_regional_official_government_2024_table_p7":OFFICIAL_SOURCE_URL,
      "national_main_G_line_codes_in_source":sorted(x for x in observed if official_class(x)=="main_G"),
      "national_regional_R_line_codes_in_source":sorted(x for x in observed if official_class(x)=="regional_R"),
      "source_full_89_sites_original_class_counts":counts,
      "original_89site_source_traits_by_main_vs_regional":source_full,
      "source_LEDA_82sites_15lines_support":{
         "qualified_original_site_n":len(scored),
         "qualified_original_same_side_pairs":sum(map(len,eligible.values())),
         "main_site_n":sum(x["class"]=="main_G" for x in scored.values()),
         "regional_site_n":sum(x["class"]=="regional_R" for x in scored.values()),
         "main_section_n":sum(official_class(k)=="main_G" for k in section),
         "regional_section_n":sum(official_class(k)=="regional_R" for k in section)},
      "original_82site_BS_minus_BT_functional_spread_main_vs_regional":contrast,
      "source_year_tracktype_and_richness_robustness":checks,
      "between_distinct_line_taxonomic_similarity_by_G_R_and_distance":geography,
      "no_measured_train_frequency":True,
      "limitations":[
        "G/R main versus regional are national legal categories, not actual 2022/2023 trains/day, freight train OD or train gusts",
        "source original railway route identifier must refer to government national line code; confirmed code universe but full location chainage independent GIS not rechecked",
        "source field sampling not random; G/R share geography, background species pool, urbanization, management, track geometry and climate",
        "only one national original vegetation study and 15 original route lines; pairwise site comparisons heavily dependent",
        "LEDA species medians are not locally measured phenotypes, cannot prove plant acclimation/evolution",
        "different zones' sampled area and detectability remain unknown",
        "no causally interpretable G/R difference and no train frequency relationship"
      ]}

def main():
    p=argparse.ArgumentParser()
    for arg in ("source_xlsx","source_audit","height","seed","output"):
        p.add_argument("--"+arg.replace("_","-"),required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.source_audit)
    for file,digest in ((a.height,HEIGHT_SHA),(a.seed,SEED_SHA)):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=digest:
            raise ValueError("Source original LEDA SHA mismatch")
    height,hq=read_leda(a.height,"height_m")
    seed,sq=read_leda(a.seed,"seed_mass_mg")
    ans=official_summary(original,height,seed)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ans,ensure_ascii=False,indent=2),encoding="utf8")
    slim={x:ans[x] for x in ("status","national_main_G_line_codes_in_source",
       "national_regional_R_line_codes_in_source","source_full_89_sites_original_class_counts",
       "original_89site_source_traits_by_main_vs_regional",
       "source_LEDA_82sites_15lines_support",
       "original_82site_BS_minus_BT_functional_spread_main_vs_regional",
       "source_year_tracktype_and_richness_robustness")}
    print("V225_OFFICIAL_G_R_SUMMARY="+json.dumps(slim,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()

#!/usr/bin/env python3
"""v323 real MIREN mountain ROAD Plot1 vs Plot3 EIVE thermal profile.

MIREN roads are independent observational habitat benchmarks, NOT railways.
Only source-verified 100m2 same-site/year records, explicit European
region-code cohort AUT/CHE/ITA/NOR, original first complete visit >=2010.
Species lists have nonrandom survey protocols: do not infer climate
causation, rail-specific changes or native-flora completeness.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table, normal_name
from analyze_Slovenia_thermal_generalist_section_v319 import section_test, section_drop
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial

MIREN_MD5="b640195cc44f252ccb6ddb28eb6fff13"
EXPECTED_ROWS=167804
EUROPE_SOURCE_REGIONS=("AUT","CHE","ITA","NOR")
EXPECTED={"Region","Road","Site","Plot","Plot_id","Year",
          "Accepted.Name.MIREN","Protocol"}
MIN_ROADS_PER_REGION=2
MIN_SITES_PER_ROAD=4

def read_miren(path):
    raw=Path(path).read_bytes()
    if hashlib.md5(raw).hexdigest()!=MIREN_MD5:
        raise ValueError("Original MIREN Zenodo source hash mismatch")
    rows=0
    plot_taxa=defaultdict(set)
    protocols=Counter()
    qc=Counter()
    with Path(path).open("r",newline="",encoding="utf-8-sig") as fh:
        reader=csv.DictReader(fh)
        if not EXPECTED.issubset(reader.fieldnames or []):
            raise ValueError("Original MIREN columns missing")
        for row in reader:
            rows+=1
            region=(row["Region"] or "").strip()
            if region not in EUROPE_SOURCE_REGIONS:continue
            qc["records_in_fixed_europe_regions"]+=1
            plot=(row["Plot"] or "").strip()
            if plot not in ("1","3"):
                qc["records_not_strict_plot1_plot3"]+=1
                continue
            p=(row["Plot_id"] or "").strip()
            if not p.endswith("."+plot):
                qc["plot_id_mismatch"]+=1;continue
            year=(row["Year"] or "").strip()
            road=(row["Road"] or "").strip()
            site=(row["Site"] or "").strip()
            if not year.isdigit() or not 2007<=int(year)<=2023 or not road or not site:
                qc["invalid_site_year"]+=1;continue
            prot=str(row.get("Protocol") or "").strip()
            protocols[(region,prot)]+=1
            species=binomial((row["Accepted.Name.MIREN"] or "").strip())
            if species is None:
                qc["not_unambiguous_binomial"]+=1;continue
            name=normal_name(species)
            key=(region,road,site,int(year),plot)
            if name in plot_taxa[key]:
                qc["duplicate_plot_year_species"]=qc.get("duplicate_plot_year_species",0)+1
            plot_taxa[key].add(name)
    if rows!=EXPECTED_ROWS:
        raise ValueError("Original source record count changed")
    return plot_taxa,{
        "source_md5":MIREN_MD5,"full_original_rows":rows,
        "fixed_regions":list(EUROPE_SOURCE_REGIONS),
        "filtered_qc":dict(qc),
        "source_protocol_record_counts":[{"region":r,"protocol":p,"records":n}
              for (r,p),n in sorted(protocols.items())],
        "observed_strict_plot1_and3_keys":len(plot_taxa)}

def select_first_complete(plot_taxa,min_year):
    siteyears=defaultdict(dict)
    for (region,road,site,year,plot),taxa in plot_taxa.items():
        if year<min_year:continue
        siteyears[(region,road,site,year)][plot]=taxa
    bysite=defaultdict(list)
    for (region,road,site,year),plots in siteyears.items():
        if set(plots)!={"1","3"}:continue
        if not plots["1"] or not plots["3"]:continue
        bysite[(region,road,site)].append((year,plots))
    pairs=[]
    for (region,road,site),options in bysite.items():
        year,plots=min(options,key=lambda x:x[0])
        pairs.append({"region":region,"road":road,"site":site,
                      "year":year,"roadside":plots["1"],"inland":plots["3"]})
    return sorted(pairs,key=lambda x:(x["region"],x["road"],x["site"]))

def traitmean(values,eive,metric,systems):
    known=[eive[s][metric] for s in values
           if s in eive and eive[s]["systems"]>=systems]
    return (sum(known)/len(known) if known else None,
            len(known)/len(values) if values else 0,len(known))

def analyze_panel(pairs,eive,metric,systems,minimum_coverage,min_species):
    groups=defaultdict(list)
    cohort=Counter()
    for pair in pairs:
        x,cx,nx=traitmean(pair["roadside"],eive,metric,systems)
        y,cy,ny=traitmean(pair["inland"],eive,metric,systems)
        region,road=pair["region"],pair["road"]
        cohort["original_complete_site_pairs"]+=1
        if nx<min_species or ny<min_species or min(cx,cy)<minimum_coverage:
            continue
        cohort["eligible_sites_before_road_region_gate"]+=1
        groups[(region,road)].append(x-y)
    road_means={k:sum(v)/len(v) for k,v in groups.items()
                if len(v)>=MIN_SITES_PER_ROAD}
    byregion=defaultdict(list)
    regional_sites=defaultdict(int)
    for (region,road),v in sorted(road_means.items()):
        byregion[region].append(v)
        regional_sites[region]+=len(groups[(region,road)])
    eligible={k:v for k,v in byregion.items()
              if len(v)>=MIN_ROADS_PER_REGION}
    road_level={region:{"roads":len(arr),"eligible_first_visit_sites":regional_sites[region],
           "road_equal_mean_roadside_minus_inland":sum(arr)/len(arr),
           "positive_roads":sum(v>0 for v in arr),
           "negative_roads":sum(v<0 for v in arr)}
           for region,arr in sorted(eligible.items())}
    region_means=[v["road_equal_mean_roadside_minus_inland"] for v in road_level.values()]
    return {"metric":metric,"minimum_EIVE_contributing_systems":systems,
            "earliest_allowed_survey_year":min(p["year"] for p in pairs) if pairs else None,
            "min_source_trait_coverage_each_plot":minimum_coverage,
            "min_numeric_taxa_each_plot":min_species,
            "qc":dict(cohort),
            "supported_region_codes":len(eligible),
            "supported_road_ids":sum(len(v) for v in eligible.values()),
            "supported_first_survey_sites":sum(regional_sites[k] for k in eligible),
            "regions":road_level,
            "region_equal_mean_roadside_minus_inland":(
                sum(region_means)/len(region_means) if region_means else None),
            "regions_roadside_higher":sum(v>0 for v in region_means),
            "regions_inland_higher":sum(v<0 for v in region_means),
            "section_direction_signflip_descriptive":section_test(eligible),
            "leave_one_region_out":section_drop(eligible)}

def analyze(csvpath,eivepath):
    plots,qa=read_miren(csvpath)
    eive,eivemeta=eive_table(eivepath)
    first10=select_first_complete(plots,2010)
    first17=select_first_complete(plots,2017)
    panels=[]
    for earliest,sub in ((2010,first10),(2017,first17)):
        for systems in (1,3):
            for coverage,minimum_n in ((.25,3),(.50,5)):
                for metric in ("position","width"):
                    p=analyze_panel(sub,eive,metric,systems,coverage,minimum_n)
                    p["min_source_year"]=earliest
                    panels.append(p)
    return {"version":"v323","source":"MIREN GLOBAL MOUNTAIN ROAD survey v3",
      "source_doi":"10.5281/zenodo.15749536",
      "source_data_paper":"10.1111/geb.70284",
      "eive_source":"10.5281/zenodo.7534792",
      "original_miren_audit":qa,"original_eive_audit":eivemeta,
      "first_complete_MIREN_site_pairs_after_2010":len(first10),
      "first_complete_MIREN_site_pairs_after_2017":len(first17),
      "fixed_europe_source_region_codes":list(EUROPE_SOURCE_REGIONS),
      "panels":panels,
      "limits":[
       "MIREN records are ROAD vegetation, NEVER independent railway observations",
       "Four European region codes were fixed before observing thermal outcomes, not all 25 MIREN global regions",
       "First observed both plots in each original site/year used to avoid repeated-site pseudo-replication",
       "Only fully source-recorded Plot 1/3 exist; positive-only lists with missing plot-year are not absence",
       "2007-08 MIREN source protocol often targets aliens, so primary source years >=2010; Protocol labels are audited but not interpreted as complete flora",
       "Roadside/interior both 100 m2 but plant pool, altitude and management differ from Slovenian ballast shoulder/toe",
       "Road-region nonrandom assignment, sparse independent regions and EIVE Europe-biased traits prevent causal/continental generality claims",
       "This does not test rail vs natural climate offset, measured track heat, temperature tolerance or species dispersal"]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--miren-csv",required=True)
    ap.add_argument("--eive-xlsx",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    d=analyze(a.miren_csv,a.eive_xlsx)
    p=Path(a.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"version":d["version"],
        "first_complete_MIREN_site_pairs_after_2010":d["first_complete_MIREN_site_pairs_after_2010"],
        "first_complete_MIREN_site_pairs_after_2017":d["first_complete_MIREN_site_pairs_after_2017"],
        "source_protocols":d["original_miren_audit"]["source_protocol_record_counts"],
        "panels":[{k:p[k] for k in (
            "metric","minimum_EIVE_contributing_systems","min_source_year",
            "min_source_trait_coverage_each_plot","min_numeric_taxa_each_plot",
            "supported_region_codes","supported_road_ids","supported_first_survey_sites",
            "regions","region_equal_mean_roadside_minus_inland",
            "section_direction_signflip_descriptive")} for p in d["panels"]]},
         ensure_ascii=False),flush=True)

if __name__=="__main__":main()

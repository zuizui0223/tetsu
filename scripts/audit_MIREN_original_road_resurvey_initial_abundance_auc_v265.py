#!/usr/bin/env python3
"""v26.5 — ORIGINAL MIREN mountain ROAD (not railway) repeated vegetation:
does initial within-plot abundance class predict later species re-recording,
and is that predictive power greater at the road edge than inland?

Original Zenodo 15749536 v3 MIREN 2007–23 CSV, exact md5 and 167804 rows.
A source physical site is (Region,Road,Site); only *strict original*
Plot1 roadside vs Plot3 52–102m interior, both observed in >=2 years,
first and last available complete paired survey. Replacement Plot1A/3A
is NEVER treated as the same permanent original plot.
Predefined PRIMARY: initial year >=2012, when full-flora surveys were
generally used. 2007–08 alien-focused cohorts are sensitivity ONLY.

AUC within each initial original plot: probability that the initial
ordinal Abundance 1–3 of a taxon re-recorded later exceeds that of a
taxon absent from later record, with ties worth 0.5. Each original site
is paired roadside–interior, then equal road means, then equal Region
code means. Sample-support gate >=4 paired sites per road and >=2
roads per region. SECONDARY sensitivity: Cover ordinal 1–8 with
matched Abundance+Cover same-site supports, not a plot-cover percentage.
Data cannot distinguish survival, mortality, field detection or
protocol differences. This is ROAD, not extra RAIL field replication.

Every first/last site necessarily changes source Protocol versions;
therefore the time-transitions are observational and not an experiment.
Do not include positive-only GBIF or infer a true taxon absence.
No raw species/location/plot-year records redistributed in output.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math
from collections import defaultdict, Counter
from pathlib import Path
from statistics import mean, median

SOURCE_MD5="b640195cc44f252ccb6ddb28eb6fff13"
SOURCE_ROWS=167804
SOURCE_DOI="10.5281/zenodo.15749536"
REQUIRED=("Plot_id","Region","Road","Site","Plot","Year","Accepted.Name.MIREN",
          "Abundance","Cover","Protocol")
PLOTS=("1","3")
METRICS=("Abundance","Cover")
PRIMARY_MIN_YEAR=2012
MIN_SITES_PER_ROAD=4
MIN_ROADS_PER_REGION=2


def numerical(x:str) -> float|None:
    if str(x).strip().upper() in ("","NA","N/A","NAN"):return None
    try: value=float(x)
    except ValueError:raise ValueError("Non-numeric ordinal source plant count or cover")
    if not math.isfinite(value) or value<=0:
        raise ValueError("Original positive source abundance/cover not positive finite")
    return value


def source_read(path:Path)->tuple[dict,dict]:
    raw=path.read_bytes()
    if hashlib.md5(raw).hexdigest()!=SOURCE_MD5:
        raise ValueError("Original official MIREN 2025 Zenodo v3 MD5 source changed")
    plant={}             # (region,road,site,year,plot) -> taxon -> (abundance,cover)
    protocol={}
    source_ids={}
    count=Counter()
    rows=0
    with io.StringIO(raw.decode("utf-8-sig"),newline="") as fh:
        reader=csv.DictReader(fh)
        if set(REQUIRED)-set(reader.fieldnames or ()):
            raise ValueError("Original author source field names changed")
        for row in reader:
            rows+=1
            plot=row["Plot"].strip()
            if plot not in PLOTS:
                count["original_excluded_plot_2_or_replacements"]+=1
                continue
            reg,road,site=(row[k].strip() for k in ("Region","Road","Site"))
            year=int(row["Year"])
            taxon=row["Accepted.Name.MIREN"].strip()
            pid=row["Plot_id"].strip()
            if not all((reg,road,site,taxon,pid)) or year<2007 or year>2023:
                raise ValueError("Missing or implausible original source plot/taxon/year")
            if not pid.endswith("."+plot):
                raise ValueError("Original plot ID conflicts with plot field")
            key=(reg,road,site,year,plot)
            if key in source_ids and source_ids[key]!=pid:
                raise ValueError("Multiple source original plot IDs at same site-year")
            source_ids[key]=pid
            ver=row["Protocol"].strip()
            if key in protocol and protocol[key]!=ver:
                raise ValueError("Contradictory protocol within source plot-year")
            protocol[key]=ver
            metrics=(numerical(row["Abundance"]),numerical(row["Cover"]))
            if key not in plant:plant[key]={}
            if taxon in plant[key]:
                if plant[key][taxon]!=metrics:
                    raise ValueError("Conflicting repeated species metric within source plot-year")
                count["source_same_taxon_plot_year_exact_duplicates_collapsed"]+=1
            else:
                plant[key][taxon]=metrics
            count["original_STRICT_plot_1_3_source_rows"]+=1
    if rows!=SOURCE_ROWS or not 1700<=len(plant)<=5500:
        raise ValueError("Genuine MIREN original author source records or plot support changed")
    all_siteyears=defaultdict(set)
    for reg,road,site,year,plot in plant:
        all_siteyears[(reg,road,site,year)].add(plot)
    complete=defaultdict(list)
    for (reg,road,site,year),plots in all_siteyears.items():
        if set(PLOTS)<=plots:
            complete[reg,road,site].append(year)
    repeated={site:sorted(ys) for site,ys in complete.items() if len(ys)>=2}
    if len(repeated)!=484:
        raise ValueError("Previously genuine MIREN 484 repeated original road site source drift")
    qc={
        "source_original_data_DOI":SOURCE_DOI,
        "source_original_CSV_MD5":SOURCE_MD5,
        "source_original_bytes":len(raw),
        "original_source_species_observation_rows":rows,
        **dict(count),
        "original_strict_1_3_plot_years":len(plant),
        "original_same_site_year_1_and_3_pair_count":sum(len(v) for v in complete.values()),
        "original_spatial_sites_with_repeated_strict_both_plots":len(repeated),
        "source_region_codes_repeated":len({site[0] for site in repeated}),
        "source_true_initial_year_2012_or_later_sites":sum(ys[0]>=PRIMARY_MIN_YEAR for ys in repeated.values()),
        "all_repeats_have_at_least_two_years":True,
        "replacement_original_A_suffix_plots_excluded":True,
    }
    return {"plants":plant,"protocol":protocol,"repeated":repeated},qc


def score_auc(initial:dict[str,tuple[float|None,float|None]],
              later:dict[str,tuple[float|None,float|None]],
              metric:str)->dict:
    ix=METRICS.index(metric)
    before={name:value[ix] for name,value in initial.items()}
    if len(before)<2:
        return {"identified":False,"reason":"fewer_than_two_source_initial_taxa"}
    if any(x is None for x in before.values()):
        return {"identified":False,"reason":"initial_source_abundance_or_cover_incomplete"}
    persist=[v for taxon,v in before.items() if taxon in later]
    missing=[v for taxon,v in before.items() if taxon not in later]
    if not persist or not missing:
        return {"identified":False,"reason":"no_both_source_recorded_fates"}
    n=len(persist)*len(missing)
    score=sum((1.0 if p>l else .5 if p==l else 0.0)
              for p in persist for l in missing)/n
    return {"identified":True,"auc":score,
            "initial_taxa":len(before),"later_recorded":len(later),
            "persistent_taxa":len(persist),"not_re_recorded_taxa":len(missing),
            "n_initial_ordinal_ties":sum(p==l for p in persist for l in missing),
            "n_between_fates_pairs":n}


def evaluate_pairs(data:dict)->tuple[list[dict],dict]:
    plant=data["plants"]
    proto=data["protocol"]
    source=data["repeated"]
    out=[]
    survey_protocol_pairs=Counter()
    for (reg,road,site),ys in sorted(source.items()):
        early,late=ys[0],ys[-1]
        record={"region":reg,"road":road,"site":site,"initial_year":early,
                "final_year":late,"years_elapsed":late-early,
                "initial_protocol_plot1":proto[(reg,road,site,early,"1")],
                "final_protocol_plot1":proto[(reg,road,site,late,"1")]}
        if proto[(reg,road,site,early,"1")]!=proto[(reg,road,site,early,"3")]:
            raise ValueError("Different protocol in originally paired roadside/interior year")
        if proto[(reg,road,site,late,"1")]!=proto[(reg,road,site,late,"3")]:
            raise ValueError("Different protocol in final originally paired roadside/interior year")
        survey_protocol_pairs[(record["initial_protocol_plot1"],record["final_protocol_plot1"])]+=1
        for plot in PLOTS:
            before=plant[(reg,road,site,early,plot)]
            after=plant[(reg,road,site,late,plot)]
            for metric in METRICS:
                record[(plot,metric)]=score_auc(before,after,metric)
        out.append(record)
    if len(out)!=484:raise ValueError("Original repeated sites were duplicated")
    return out,{
        "source_protocol_transition_count":{
            f"{initial} -> {last}":int(n)
            for (initial,last),n in sorted(survey_protocol_pairs.items())
        },
        "n_original_sites_with_identical_source_protocol_at_both_observations":
            sum(row["initial_protocol_plot1"]==row["final_protocol_plot1"] for row in out),
        "all_484_complete_sites_first_and_last_used_only_once":True
    }


def aggregate(rows:list[dict],metric:str,period:str)->dict:
    if period=="post2012":valid_period=lambda row:row["initial_year"]>=PRIMARY_MIN_YEAR
    elif period=="early2007_08":valid_period=lambda row:row["initial_year"]<PRIMARY_MIN_YEAR
    elif period=="all_first_last":valid_period=lambda row:True
    else:raise ValueError("Source survey period not predeclared")
    paired=[]
    for row in rows:
        if not valid_period(row):continue
        edge=row[("1",metric)];inside=row[("3",metric)]
        if not edge["identified"] or not inside["identified"]:continue
        paired.append({"region":row["region"],"road":row["road"],
                       "site":row["site"],"year_initial":row["initial_year"],
                       "auc_roadside":edge["auc"],"auc_interior":inside["auc"],
                       "delta":edge["auc"]-inside["auc"]})
    byroad=defaultdict(list)
    for row in paired:
        byroad[row["region"],row["road"]].append(row)
    region_roads=defaultdict(list)
    for (region,road),rs in byroad.items():
        if len(rs)<MIN_SITES_PER_ROAD:continue
        region_roads[region].append({
            "n_source_sites":len(rs),
            "auc_roadside":mean(x["auc_roadside"] for x in rs),
            "auc_interior":mean(x["auc_interior"] for x in rs),
            "roadside_minus_interior":mean(x["delta"] for x in rs)
        })
    regions={}
    for region,roads in sorted(region_roads.items()):
        if len(roads)<MIN_ROADS_PER_REGION:continue
        regions[region]={
            "source_real_roads":len(roads),
            "source_real_spatial_sites":sum(x["n_source_sites"] for x in roads),
            "source_road_equal_auc_roadside":mean(x["auc_roadside"] for x in roads),
            "source_road_equal_auc_interior":mean(x["auc_interior"] for x in roads),
            "source_road_equal_delta":mean(x["roadside_minus_interior"] for x in roads)
        }
    vals=[v["source_road_equal_delta"] for v in regions.values()]
    supported=set(regions)
    return {
        "period":period,"metric":metric,
        "source_original_sites_with_two_plot_AUC_identified_all_road_sizes":len(paired),
        "source_regions_in_supported_equal_region_summary":len(regions),
        "source_roads_in_supported_equal_region_summary":sum(v["source_real_roads"] for v in regions.values()),
        "source_sites_in_supported_equal_region_summary":sum(v["source_real_spatial_sites"] for v in regions.values()),
        "mean_original_region_equal_roadside_minus_interior_AUC":mean(vals) if vals else None,
        "regions_roadside_AUC_greater":sum(x>0 for x in vals),
        "regions_roadside_AUC_lower":sum(x<0 for x in vals),
        "regions_exactly_tied":sum(x==0 for x in vals),
        "delete_one_source_region_out_range":
            [min((sum(vals)-v)/(len(vals)-1) for v in vals),
             max((sum(vals)-v)/(len(vals)-1) for v in vals)] if len(vals)>1 else None,
        "original_region_code_aggregate_only_not_country_replicates":regions,
        "regional_mean_is_not_unbiased_world_population_inference":True,
        "experimental_traffic_mowing_or_detection_causality_not_identified":True
    }


def matched_metric_comparison(rows:list[dict],period="post2012")->dict:
    # Same 2014-style risk set and same ORIGINAL SOURCE spatial sites for the
    # initial ordinal abundance and ordinal cover metric: no changing cohort.
    base=[]
    for row in rows:
        if row["initial_year"]<PRIMARY_MIN_YEAR:continue
        if all(row[(plot,metric)]["identified"] for plot in PLOTS for metric in METRICS):
            base.append(row)
    return {
        "same_source_sites_for_Abundance_and_Cover":len(base),
        "ordinal_abundance_on_identical_source_sites":aggregate(base,"Abundance",period),
        "ordinal_cover_on_identical_source_sites":aggregate(base,"Cover",period),
        "not_two_independent_field_replications":True
    }


def run(path:Path)->dict:
    source,qc=source_read(path)
    rows,periodqa=evaluate_pairs(source)
    all_period=aggregate(rows,"Abundance","all_first_last")
    main=aggregate(rows,"Abundance","post2012")
    old=aggregate(rows,"Abundance","early2007_08")
    match=matched_metric_comparison(rows)
    if (main["source_sites_in_supported_equal_region_summary"]<100 or
        main["source_regions_in_supported_equal_region_summary"]<4):
        raise ValueError("Modern ROAD source primary support unexpectedly poor")
    return {
        "version":"v26.5","status":"SOURCE_VERIFIED_MIREN_2007_2023_ROADS_NOT_RAIL",
        "author_original_source_QA":qc,"protocol_QA":periodqa,
        "primary_modern_abundance_ordinal_road_vs_interior":main,
        "sensitivity_all_years_alien_focus_early_included":all_period,
        "sensitivity_early2007_08_only_not_full_flora":old,
        "same_original_site_cover_and_abundance_measurement_sensitivity":match,
        "interpretation":"within-initial-plot ordinal 1-3 Abundance predicts re-recording versus not re-recording, not a measure of field individual survival",
        "independent_road_database_not_railway_field_replication":True,
        "MIREN_regions_not_random_world_countries":True,
        "v2007_protocol_alien_focused_not_modern_full_species_frame":True,
        "none_of_484_repeated_site_pairs_has_identical_Protocol_version":periodqa[
            "n_original_sites_with_identical_source_protocol_at_both_observations"]==0,
        "cannot_compare_absolute_Kashmir_numeric_AUC_directly_to_MIREN_ordinal_AUC":True,
        "original_source_plot_species_lists_and_source_GPS_not_redistributed":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    report=run(args.original_csv)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
        "version":report["version"],"status":report["status"],
        "original_source_QA":report["author_original_source_QA"],
        "modern":{key:value for key,value in report[
            "primary_modern_abundance_ordinal_road_vs_interior"].items()
            if key!="original_region_code_aggregate_only_not_country_replicates"},
        "paired_measurement_abundance":{
            key:value for key,value in report[
                "same_original_site_cover_and_abundance_measurement_sensitivity"][
                    "ordinal_abundance_on_identical_source_sites"].items()
            if key!="original_region_code_aggregate_only_not_country_replicates"},
        "paired_measurement_cover":{
            key:value for key,value in report[
                "same_original_site_cover_and_abundance_measurement_sensitivity"][
                    "ordinal_cover_on_identical_source_sites"].items()
            if key!="original_region_code_aggregate_only_not_country_replicates"}
    },ensure_ascii=False))


if __name__=="__main__":main()

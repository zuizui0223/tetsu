#!/usr/bin/env python3
"""v27.2 original 2012+ MIREN mountain ROAD flora native-vs-alien
same-site assemblage coupling under source richness+frequency fixed null.

The genuine MIREN 2025 Zenodo CSV (MD5 b640195cc44f252ccb6ddb28eb6fff13)
has native/alien/unknown STATUS per region x species. Two conflicting status
assignments in prior source v9.9 were safely classified as UNKNOWN, never
silently treated as native or alien. This source re-read independently
verifies that each actual original 2012+ Plot1/Plot3 observed taxon has a
classification in its native REGION source pool.

For 669 original physical sites/48 roads/16 source Region codes that
previously qualified in v26.9/27.0, stratify into:
  NATIVE (both source plots contain >=1 source Native taxon)
  ALIEN (both source plots contain >=1 source Alien taxon)
  BOTH_STATUS (each source plot contains >=1 Native AND >=1 Alien),
               paired-site sensitivity.
For each category and original road, retain sites with positive recorded
taxa in *both zones*, min4 eligible sites per road, min2 roads per region.
Species are never imputed to fill source missingness. The corridor-side
species sets are fixed; reference-zone site x species checkerboard trades
preserve exact reference site status-specific richness and same-ROAD
species occurrence frequencies. The true original species and environment
remain separate from any railway study.

Descriptive source mode also compares site-equal fraction of ALIEN among
original corridor-only versus original interior-only taxa, excluding source
Unknown status. Region/road-balanced, not per-taxonomic-record pseudoreps.

All resulting p-like ranks are technical correlated MCMC diagnostics,
NOT experimental causal p values. Native/alien support cohorts may differ,
so their beta differences cannot be directly interpreted as a
status-specific effect unless same-site matched support is reported.
No actual vehicle dispersal, train traffic, mowing or evolution measured.
Never redistribute source locality/species rows.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 import source_read
from analyze_original_RAIL_vs_global_ROAD_taxonomic_microhabitat_beta_v269 import road_first_surveys
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    modern_roads,sampler,finish_null,METRICS
)

SOURCE_MD5="b640195cc44f252ccb6ddb28eb6fff13"
CATEGORIES=("Native","Alien")
DRAWS=79
BURNIN=2500
SPACING=300
SEED=272118
MIN_SITES=4
MIN_ROADS=2


def statuses(path:Path)->tuple[dict,dict]:
    if hashlib.md5(path.read_bytes()).hexdigest()!=SOURCE_MD5:
        raise ValueError("Author MIREN genuine original source checksum changed")
    value=defaultdict(set)
    nrows=0
    with path.open(encoding="utf-8-sig",newline="") as handle:
        csv_reader=csv.DictReader(handle)
        if not {"Region","Accepted.Name.MIREN","Status","Plot"}<=set(csv_reader.fieldnames or []):
            raise ValueError("Genuine MIREN source native status metadata missing")
        for row in csv_reader:
            nrows+=1
            code=row["Plot"].strip()
            if code not in ("1","2","3"):continue
            region=row["Region"].strip()
            taxon=row["Accepted.Name.MIREN"].strip()
            label=row["Status"].strip()
            if label.upper() in ("","NA","N/A"):
                label="Unknown"
            if label not in ("Native","Alien","Unknown"):
                raise ValueError("Original MIREN nativity has undocumented class")
            if not region or not taxon:
                continue
            value[region,taxon].add(label)
    if nrows!=167804:
        raise ValueError("MIREN source nativity scan population changed")
    disagree={key for key,v in value.items() if len(v)>1}
    resolved={key:(next(iter(v)) if len(v)==1 else "Unknown") for key,v in value.items()}
    qa={
        "source_rows":nrows,
        "classified_region_taxon_pairs":len(resolved),
        "conflicting_original_region_taxon_status_conservatively_UNKNOWN":len(disagree),
        "unknown_region_taxon_pairs":sum(x=="Unknown" for x in resolved.values()),
        "native_region_taxon_pairs":sum(x=="Native" for x in resolved.values()),
        "alien_region_taxon_pairs":sum(x=="Alien" for x in resolved.values()),
        "source_status_regional_not_global_plant_invasiveness":True
    }
    return resolved,qa


def filtered_groups(source_path:Path)->tuple[dict,dict,dict]:
    source,qc=source_read(source_path)
    first,firstqc=road_first_surveys(source)
    roadgroups,mainqa=modern_roads(first)
    status,statusqa=statuses(source_path)
    classified=defaultdict(lambda:defaultdict(list))
    site_record=[]
    site_unknown=set()
    for region,roads in sorted(roadgroups.items()):
        for road,sites in sorted(roads.items()):
            for site in sites:
                a=site["corridor"];b=site["reference"]
                if any((region,t) not in status for t in a|b):
                    raise ValueError("Original MIREN 2012+ plot source species has no author Region status")
                a_by={kind:{t for t in a if status[region,t]==kind} for kind in CATEGORIES}
                b_by={kind:{t for t in b if status[region,t]==kind} for kind in CATEGORIES}
                for kind in CATEGORIES:
                    if a_by[kind] and b_by[kind]:
                        classified[kind][region,road].append({
                            "site":site["site"],"corridor":a_by[kind],"reference":b_by[kind]})
                if all(a_by[k] and b_by[k] for k in CATEGORIES):
                    classified["both_Native"][region,road].append({
                        "site":site["site"],"corridor":a_by["Native"],"reference":b_by["Native"]})
                    classified["both_Alien"][region,road].append({
                        "site":site["site"],"corridor":a_by["Alien"],"reference":b_by["Alien"]})
                # Site-level alien fraction of original corridor-only / interior-only.
                unique_a=a-b;unique_b=b-a
                known_a=[status[region,t] for t in unique_a if status[region,t] in CATEGORIES]
                known_b=[status[region,t] for t in unique_b if status[region,t] in CATEGORIES]
                nunknown=sum(status[region,t]=="Unknown" for t in a|b)
                if nunknown:site_unknown.add((region,road,site["site"]))
                site_record.append({
                    "region":region,"road":road,
                    "corridor_only_n":len(unique_a),"reference_only_n":len(unique_b),
                    "known_corridor_only_n":len(known_a),"known_reference_only_n":len(known_b),
                    "alien_share_of_known_corridor_only":(sum(x=="Alien" for x in known_a)/len(known_a) if known_a else None),
                    "alien_share_of_known_reference_only":(sum(x=="Alien" for x in known_b)/len(known_b) if known_b else None)
                })
    qa={
        "source_original_MIREN_QA":qc,
        "source_original_first_year_QA":firstqc,
        "original_modern_source_road_region_sites":mainqa,
        "original_status_source_QA":statusqa,
        "source_original_sites_with_any_UNKNOWN_status_species":len(site_unknown),
        "status_specific_sites_positive_both_zones_per_status":{
            key:sum(map(len,byroad.values())) for key,byroad in classified.items()
        },
        "no_imputed_source_nonobserved_native_or_alien_taxa":True,
    }
    return classified,site_record,qa


def balance_by_roads(raw:dict)->dict:
    regional=defaultdict(dict)
    for (reg,road),values in raw.items():
        if len(values)>=MIN_SITES:
            regional[reg][road]=sorted(values,key=lambda x:str(x["site"]))
    return {reg:roads for reg,roads in regional.items() if len(roads)>=MIN_ROADS}


def analyze_status_panel(byroad:dict,*,draws:int=DRAWS,burnin:int=BURNIN,
                         spacing:int=SPACING,seed:int=SEED)->dict:
    region=balance_by_roads(byroad)
    if not region:
        return {"identified":False,"n_source_regions":0}
    results={}
    for i,(r,roads) in enumerate(sorted(region.items())):
        vals=[]
        for j,(road,sites) in enumerate(sorted(roads.items())):
            x=sampler([s["corridor"] for s in sites],
                      [s["reference"] for s in sites],
                      draws=draws,burnin=burnin,spacing=spacing,seed=seed+631*i+109*j)
            vals.append(x)
        observed={metric:mean(x["observed"][metric] for x in vals) for metric in METRICS}
        nul=[{metric:mean(x["null_draws"][k][metric] for x in vals)
              for metric in METRICS} for k in range(draws)]
        results[r]={
            "source_original_roads":len(vals),
            "source_original_site_pairs":sum(x["original_site_count"] for x in vals),
            "n_roads_without_accepted_swaps":sum(x["accepted_swaps_total"]==0 for x in vals),
            "observed_and_frequency_null":finish_null(observed,nul)
        }
    obs={m:mean(x["observed_and_frequency_null"][m]["observed"] for x in results.values())
         for m in METRICS}
    # Technical macro null sampling, allowing each original road MCMC draw;
    # note averaging per-Region expectations then comparing is descriptive,
    # not an exact synchronized independent null.
    overall={
        m:{
            "mean_source_region_observed":obs[m],
            "mean_region_equal_frequency_null":mean(x["observed_and_frequency_null"][m]["frequency_null_mean"] for x in results.values()),
            "mean_region_equal_observed_minus_null":mean(x["observed_and_frequency_null"][m]["observed_minus_frequency_null"] for x in results.values()),
            "regions_observed_below_null":sum(x["observed_and_frequency_null"][m]["observed_minus_frequency_null"]<0 for x in results.values()),
            "regions_observed_above_null":sum(x["observed_and_frequency_null"][m]["observed_minus_frequency_null"]>0 for x in results.values())
        } for m in METRICS
    }
    return {
        "identified":True,
        "n_source_regions":len(results),
        "n_original_roads":sum(x["source_original_roads"] for x in results.values()),
        "n_original_matched_source_sites":sum(x["source_original_site_pairs"] for x in results.values()),
        "source_original_roads_without_any_frequency_swap":sum(x["n_roads_without_accepted_swaps"] for x in results.values()),
        "equal_road_then_equal_region_results":overall,
        "source_region_codes_not_independent_study_replicates":True,
        "region_code_aggregate_support":results
    }


def alien_fraction_comparison(records:list[dict])->dict:
    within=defaultdict(list)
    for x in records:
        within[x["region"],x["road"]].append(x)
    roads=defaultdict(list)
    for (reg,road),sites in within.items():
        vals_a=[x["alien_share_of_known_corridor_only"] for x in sites if x["alien_share_of_known_corridor_only"] is not None]
        vals_b=[x["alien_share_of_known_reference_only"] for x in sites if x["alien_share_of_known_reference_only"] is not None]
        paired=[x for x in sites if x["alien_share_of_known_corridor_only"] is not None and x["alien_share_of_known_reference_only"] is not None]
        if len(paired)<MIN_SITES:continue
        roads[reg].append({
            "n_matched_source_sites":len(paired),
            "road_equal_site_alien_share_corridor_only":mean(x["alien_share_of_known_corridor_only"] for x in paired),
            "road_equal_site_alien_share_reference_only":mean(x["alien_share_of_known_reference_only"] for x in paired),
            "road_equal_site_alien_share_difference":mean(x["alien_share_of_known_corridor_only"]-x["alien_share_of_known_reference_only"] for x in paired)
        })
    eligible={reg:rs for reg,rs in roads.items() if len(rs)>=MIN_ROADS}
    regional={reg:{"n_roads":len(rs),
                   "n_original_sites":sum(x["n_matched_source_sites"] for x in rs),
                   "mean_site_road_alien_share_corridor_only":mean(x["road_equal_site_alien_share_corridor_only"] for x in rs),
                   "mean_site_road_alien_share_reference_only":mean(x["road_equal_site_alien_share_reference_only"] for x in rs),
                   "delta":mean(x["road_equal_site_alien_share_difference"] for x in rs)
                  } for reg,rs in sorted(eligible.items())}
    return {
        "source_region_codes":len(regional),
        "source_roads":sum(x["n_roads"] for x in regional.values()),
        "source_sites":sum(x["n_original_sites"] for x in regional.values()),
        "mean_region_equal_alien_fraction_among_corridor_only_taxa":mean(x["mean_site_road_alien_share_corridor_only"] for x in regional.values()) if regional else None,
        "mean_region_equal_alien_fraction_among_interior_only_taxa":mean(x["mean_site_road_alien_share_reference_only"] for x in regional.values()) if regional else None,
        "mean_region_equal_alien_fraction_diff":mean(x["delta"] for x in regional.values()) if regional else None,
        "regions_corridor_only_more_alien":sum(x["delta"]>0 for x in regional.values()),
        "source_native_alien_status_Unknown_excluded_not_imputed":True,
        "source_unique_taxon_averages_are_not_original_independent_field_sites":True
    }


def run(source_path:Path,draws:int=DRAWS,burnin:int=BURNIN,spacing:int=SPACING)->dict:
    bystatus,source_sites,qa=filtered_groups(source_path)
    outputs={}
    for i,status in enumerate(("Native","Alien","both_Native","both_Alien")):
        outputs[status]=analyze_status_panel(bystatus[status],draws=draws,burnin=burnin,
                                             spacing=spacing,seed=SEED+100001*i)
    return {
        "version":"v27.2",
        "source_study":"MIREN mountain ROADS, not railway vegetation",
        "original_author_status_and_biogeographical_QA":qa,
        "source_regional_status_group_fixed_margin_nulls":outputs,
        "original_corridor_only_vs_reference_only_alien_fraction":alien_fraction_comparison(source_sites),
        "status_categories_based_on_original_MIREN_Region_not_global_nativity":True,
        "source_Unknown_missing_and_conflicted_taxa_not_imputed_to_Alien_or_Native":True,
        "study_not_train_mowing_wind_or_local_fitness_experiment":True,
        "no_original_taxon_name_GPS_or_site_level_source_records_output":True
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--miren-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURNIN)
    p.add_argument("--spacing",type=int,default=SPACING)
    a=p.parse_args()
    out=run(a.miren_csv,a.draws,a.burnin,a.spacing)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for status,panel in out["source_regional_status_group_fixed_margin_nulls"].items():
        print(json.dumps({"version":"v27.2","status":status,"identified":panel["identified"],
          "regions":panel.get("n_source_regions"),
          "roads":panel.get("n_original_roads"),
          "sites":panel.get("n_original_matched_source_sites"),
          "raw_taxonomy_beta":panel.get("equal_road_then_equal_region_results",{}).get("beta_jac"),
          "raw_shared":panel.get("equal_road_then_equal_region_results",{}).get("shared_taxa")},ensure_ascii=False))
    print(json.dumps({"version":"v27.2","source_status":out["original_author_status_and_biogeographical_QA"]["original_status_source_QA"],
          "original_alien_nativity_exclusives":out["original_corridor_only_vs_reference_only_alien_fraction"]},ensure_ascii=False))


if __name__=="__main__":main()

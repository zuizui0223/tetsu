#!/usr/bin/env python3
"""v12.1: REAL CzechVeg RAIL-text vs ROAD-text plant functional comparisons
with geographically and temporally co-observed original plot-area strata.

Original external authors' CzechVeg-Open 2026 v1.1 full 117739 vegetation
plot headers, 2510340 original source species×plot records (Zenodo
10.5281/zenodo.19231152), original published Pladias Sádlo 2018 strategy
and original Raunkiaer life form. No tetsu original field collection.

Unlike v11.9's unbalanced pooled source class preference, include only
blocks with BOTH plot types actually recorded: original phytogeographic
region, sampling year epoch, and original recorded plot-area bin.
Source plot vegetation category, actual rail right-of-way, management,
functional selection and transport remain NOT matched/identified.

All six predeclared overlap classes run and ALL outputs retained:
A exact area / 10yr / region with >=2 rail and road source plots
B area class / 10yr / region, >=2 each
C area class / 20yr / region, >=2 each
D 10yr / region, >=3 each (area NOT controlled)
E area class / region, >=3 each (year NOT controlled)
F region only, >=5 each (year and area NOT controlled)
Run separately for STRICT HABITAT-only and BROAD HABITAT/LOCALITY rail
and road source categories, with overlapping area range 1..1000 m2.

For each retained matched block b, target propensity weight
omega_b=min(n_rail,b, n_road,b) (standardized to sum 1).
Each taxon source q_s = sum_b omega_b[
  observed rail source plot proportion_s,b -
  observed road source plot proportion_s,b
].
Only source taxon×lifeform×strategy exact matches retained.
Then compare strict wind-Epilobium/Phragmites vs remaining dispersal
strategy within therophyte vs other original source life-form labels,
with equal species versus observed matched-source plot-occupancy weight.
NO random treatment/exposure/field causal effect or cross-country same
within-railway BS/BT effect.
"""
from __future__ import annotations
import argparse,csv,json,hashlib,math,re
from collections import Counter,defaultdict
from pathlib import Path

from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    HEADERS_NAME,SPECIES_NAME,verify_file,parse_traits,plot_labels,
    plot_taxa,guild_panel,
)
from analyze_CzechVeg_rail_road_region_robustness_v120 import validated_source
from analyze_Slovenia_rail_species_occupancy_weighting_v110 import group_estimator


MATCHING_SPECS=(
    ("A_region_decade_exact_area_min2",10,"exact",2),
    ("B_region_decade_area4_min2",10,"class4",2),
    ("C_region_20year_area4_min2",20,"class4",2),
    ("D_region_decade_no_area_min3",10,"none",3),
    ("E_region_area4_no_year_min3",0,"class4",3),
    ("F_region_only_min5",0,"none",5),
)
NO_TAXON_FUZZY_OR_GENUS_MATCH=True


def year_from_source(v):
    s=str(v or "").strip()
    if not re.fullmatch(r"\d{4}(?:\.0)?",s):
        return None
    n=int(float(s))
    return n if 1850<=n<=2026 else None


def area_key(area,mode):
    if not 1<=area<=1000:
        raise ValueError("Plot surface outside original source positive area scope")
    if mode=="exact":return area
    if mode=="none":return None
    if mode=="class4":
        if area<=4:return "1_to_4"
        if area<=16:return "gt4_to_16"
        if area<=64:return "gt16_to_64"
        if area<=256:return "gt64_to_256"
        return "gt256_to_1000"
    raise ValueError("Unsupported plot-area balance definition")


def header_source_context(header_csv,plotmeta):
    """Only source named original fields, not surrogate date from upload."""
    stats=Counter()
    out={}
    with open(header_csv,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        columns=reader.fieldnames or []
        veg_names=[k for k in columns if re.search(r"syntax|veget|eunis|habitat|formation|class",k,re.I)]
        for row in reader:
            stats["original_csv_rows"]+=1
            rid=str(row.get("RELEVE_NR") or "").strip()
            if rid not in plotmeta:continue
            x=plotmeta[rid]
            yr=year_from_source(row.get("YEAR"))
            region=str(row.get("SKALICKY1988") or "").strip()
            if yr is None:stats["candidate_with_missing_or_unsupported_year"]+=1
            if not region or region=="NA":stats["candidate_with_missing_phyto_region"]+=1
            if yr is None or not region or region=="NA":continue
            out[rid]={**x,"source_year_int":yr,"source_phyto_region":region}
            stats["candidate_with_valid_year_and_region"]+=1
    if stats["original_csv_rows"]!=117739:
        raise ValueError("Source real CzechVeg full published plot header table changed")
    return out,{
        "original_full_header_rows":stats["original_csv_rows"],
        "original_area_valid_plot_rows":len(plotmeta),
        "candidate_with_missing_or_unsupported_year":stats["candidate_with_missing_or_unsupported_year"],
        "candidate_with_missing_phyto_region":stats["candidate_with_missing_phyto_region"],
        "eligible_original_source_plot_year_and_region_rows":len(out),
        "source_original_columns_with_possible_vegetation_class_labels":veg_names[:30],
        "vegetation_classification_not_harmonized_to_one_source_standard":True,
    }


def block_key(meta,spec):
    name,years,area_mode,min_per_arm=spec
    year=int(meta["source_year_int"])
    reg=meta["source_phyto_region"]
    return (reg,year//years*years if years else None,
            area_key(meta["area"],area_mode))


def matched_blocks(meta,source_taxa,label,spec):
    """No data-leakage outcome matching: source sample unit and original
    metadata ONLY before the plant traits/outcome is inspected."""
    name,_,_,min_per_arm=spec
    groups=defaultdict(lambda:{"rail":[],"road":[]})
    counts=Counter()
    for site,rec in sorted(meta.items()):
        grp=rec[label]
        if grp not in ("rail","road"):continue
        # Zero source vascular species rows may signal source survey completeness,
        # cannot assert a source plot genuinely surveyed zero vascular plants.
        if not source_taxa.get(site):
            counts["original_candidate_plot_without_any_vascular_records_excluded"]+=1
            continue
        k=block_key(rec,spec)
        groups[k][grp].append(site)
        counts["source_candidate_plots_with_vascular_data"]+=1
    valid={k:v for k,v in groups.items() if
           len(v["rail"])>=min_per_arm and len(v["road"])>=min_per_arm}
    diag={
        "source_comparability_rule":name,
        "minimum_sampled_source_rail_and_road_plots_per_original_block":min_per_arm,
        "original_source_candidate_plots_with_vascular_records":
            counts["source_candidate_plots_with_vascular_data"],
        "candidate_plots_with_no_source_vascular_records_excluded":
            counts["original_candidate_plot_without_any_vascular_records_excluded"],
        "n_original_source_region_year_area_candidate_blocks":len(groups),
        "n_original_source_blocks_retaining_both_arms":len(valid),
        "n_source_rail_plots_on_matched_common_support":
            sum(len(v["rail"]) for v in valid.values()),
        "n_source_road_plots_on_matched_common_support":
            sum(len(v["road"]) for v in valid.values()),
        "n_source_rail_plot_loss_from_unmatched_strata":
            sum(len(v["rail"]) for v in groups.values())-
            sum(len(v["rail"]) for v in valid.values()),
        "n_source_road_plot_loss_from_unmatched_strata":
            sum(len(v["road"]) for v in groups.values())-
            sum(len(v["road"]) for v in valid.values()),
        "n_distinct_original_phyto_regions_on_matched_common_support":
            len({k[0] for k in valid}),
        "n_original_block_weights_effective_n_min_arm":
            sum(min(len(v["rail"]),len(v["road"])) for v in valid.values()),
        "top_original_geo_year_area_blocks_counts_only":[
            {"source_region":str(k[0]),"source_year_bin":k[1],
             "source_area_bin":str(k[2]),"n_rail":len(v["rail"]),
             "n_road":len(v["road"])}
            for k,v in sorted(valid.items(),
                              key=lambda x:(-min(len(x[1]["rail"]),len(x[1]["road"])),
                                            str(x[0])))[:15]],
    }
    return valid,diag


def block_standardized_guild(blocks,plot_taxon):
    """The source biological score has common-support standardization:
    common block mass and equal observed source plot sampling fraction
    per site type within each region-year-area block."""
    if not blocks:
        return {"identified":False,"cause":"No observed original same-region/year/area support"}
    weights={k:min(len(v["rail"]),len(v["road"])) for k,v in blocks.items()}
    total=sum(weights.values())
    if not total:
        raise ValueError("Zero block overlap weight on valid source data")
    q=defaultdict(lambda:{"score":0.,"occ":0,"wind":None,"life":None})
    guild_coverage=Counter()
    total_original_incidence=0
    for block in sorted(blocks,key=str):
        siteclass=blocks[block]
        a=len(siteclass["rail"]);b=len(siteclass["road"])
        wt=weights[block]/total
        counts=defaultdict(lambda:[0,0,None,None])
        for group,idx in (("rail",0),("road",1)):
            for site in siteclass[group]:
                for species,disp,life in plot_taxon.get(site,[]):
                    total_original_incidence+=1
                    if disp is None or life is None:continue
                    record=counts[species]
                    if (record[2] is not None and record[2]!=disp or
                        record[3] is not None and record[3]!=life):
                        raise ValueError("Source Pladias trait changed within matched source species")
                    record[idx]+=1
                    record[2]=disp;record[3]=life
                    guild_coverage["plot_species_original_joint_trait_incidence"]+=1
        for species,(nrail,nroad,d,l) in counts.items():
            cell=q[species]
            wind=d in ("Epilobium","Phragmites")
            lf="strict_annual" if l=="therophyte" else "other_known_lifespan"
            if cell["wind"] not in (None,wind) or cell["life"] not in (None,lf):
                raise ValueError("One source name cannot have multiple source strategies")
            cell["wind"]=wind
            cell["life"]=lf
            cell["score"]+=wt*(nrail/a-nroad/b)
            cell["occ"]+=nrail+nroad
    rows=[{
        "taxon_code":species,"wind":v["wind"],"stratum":v["life"],
        "n_original_locations_present_either_zone":v["occ"],
        "occupancy_preference_BS_minus_BT":v["score"],
    } for species,v in sorted(q.items()) if v["occ"]>0]
    equal=group_estimator(rows,1,"species")
    freq=group_estimator(rows,1,"opportunities")
    if equal["identified"]!=freq["identified"]:
        raise ValueError("Source same taxon universe not recognized with two weights")
    base={
        "identified":equal["identified"],
        "n_original_joint_trait_species_source_common_support":len(rows),
        "n_original_joint_trait_species_with_source_strict_wind_code":
            sum(r["wind"] for r in rows),
        "n_original_matched_plot_by_taxon_trait_recognized_incidence":
            guild_coverage["plot_species_original_joint_trait_incidence"],
        "all_original_incidence_matching_blocks_including_unknown_traits":
            total_original_incidence,
        "fraction_joint_trait_incidence_coverage":(
            guild_coverage["plot_species_original_joint_trait_incidence"]/
            total_original_incidence if total_original_incidence else None),
        "n_source_standardization_blocks":len(blocks),
        "effective_overlapped_source_plot_pairs_from_block_minimums":total,
        "source_species_prevalence_per_block_standardized_before_guild_mean":True,
        "RAIL_vs_ROAD_not_within_railway_ballast_shoulder_vs_toe":True,
    }
    if base["identified"]:
        e=equal["wind_minus_other_preference_lifespan_adjusted"]
        w=freq["wind_minus_other_preference_lifespan_adjusted"]
        base.update({
            "original_source_species_equal_strict_wind_vs_other_prevalence_delta":e,
            "original_source_occ_weight_strict_wind_vs_other_prevalence_delta":w,
            "original_source_weighting_gap_equal_minus_occ":e-w,
            "opposite_source_estimand_signs":e*w<0,
            "original_source_known_traits_by_lifespan_and_diaspore":
                equal["by_stratum"],
        })
    return base


def diagnostic(original,plot_taxon,source_label,spec):
    matched,summary=matched_blocks(original,plot_taxon,source_label,spec)
    matched_sites={p for v in matched.values() for p in v["rail"]+v["road"]}
    pooled=guild_panel(original,plot_taxon,source_label,matched_sites)
    standardized=block_standardized_guild(matched,plot_taxon)
    return {
        "source_label_definition":source_label,
        "source_common_support_matching_spec":spec[0],
        "original_source_sample_restriction":summary,
        "original_pooled_prevalence_on_EXACT_matched_plots_BEFORE_block_standardizing":{
            "identified":pooled["identified"],
            "n_source_taxa":pooled.get("n_source_unique_joint_trait_taxa_in_original_rail_or_road_vegetation"),
            "species_equal":pooled.get("strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"),
            "source_occ_weight":pooled.get("strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap"),
        },
        "within_original_source_region_time_area_block_standardized_trait_comparison":
            standardized,
    }


def validate_original_joins(zenodo,source_head,source_sp,disp_file,life_file,
                            original_source119):
    verified=validated_source(original_source119)
    a=verify_file(zenodo,HEADERS_NAME,source_head)
    b=verify_file(zenodo,SPECIES_NAME,source_sp)
    previous=verified["original_source_Czechveg_header_provenance"]
    earlier=verified["original_source_Czechveg_species_provenance"]
    if a["checksum"]!=previous["checksum"] or b["checksum"]!=earlier["checksum"]:
        raise ValueError("Original full Czech source differs from verified v11.9")
    disp,life,traits=parse_traits(disp_file,life_file)
    baseline=verified["original_source_Pladias_trait_quality"]
    keys=("source_Pladias_original_diaspore_rows","source_Pladias_original_lifeform_rows",
          "n_unique_strategy_taxon_names","n_unique_strict_unambiguous_strategy_names",
          "n_source_strategy_conflicting_taxa_excluded",
          "n_source_lifeform_unambiguous_names","n_source_lifeform_ambiguous_taxa_excluded",
          "source_dispersal_code_counts")
    changed=[k for k in keys if traits.get(k)!=baseline.get(k)]
    if changed:
        raise ValueError("Original Pladias strategy/ life form biological inventory changed "+
                         ",".join(changed))
    return disp,life,{"original_vegetation_header_source":a,
                      "original_vegetation_species_source":b,
                      "Pladias_all_original_source_trait_table_aggregates_unchanged":True,
                      "source_trait_file_binary_SHA_fresh": {
                          "strategy":traits["source_strategy_SHA256"],
                          "lifeform":traits["source_lifeform_SHA256"]},
                      "original_external_plant_trait_rights_respected":True}


def run(meta,source_headers,source_species,source_disp,source_life,source119):
    disp,life,proof=validate_original_joins(meta,source_headers,source_species,
                                             source_disp,source_life,source119)
    plots,original_plot_quality=plot_labels(source_headers)
    source_context,source_context_quality=header_source_context(source_headers,plots)
    taxa,source_taxa_quality=plot_taxa(source_species,source_context,disp,life)
    if source_taxa_quality["n_original_source_species_csv_records_scanned"]!=2510340:
        raise ValueError("Real CzechVeg original species row count changed")
    outcomes={}
    for kind in ("strict","broad"):
        outcomes[kind]={}
        for spec in MATCHING_SPECS:
            outcomes[kind][spec[0]]=diagnostic(source_context,taxa,kind,spec)
    return {
        "complete":True,"version":"v12.1",
        "published_source_authors":"Chytrý et al. (2026) CzechVeg-Open v1.1 and Pladias Sádlo et al. 2018",
        "original_CzechVeg_DOI":"10.5281/zenodo.19231152",
        "original_Pladias_DOI":"10.23855/preslia.2021.001",
        "published_source_checks":proof,
        "original_full_plot_audit":original_plot_quality,
        "original_temporal_and_area_candidate_audit":source_context_quality,
        "original_joint_trait_source_sample_matching":source_taxa_quality,
        "six_declared_original_source_common_support_specs":list(x[0] for x in MATCHING_SPECS),
        "results_by_strict_vs_broad_source_railway_label_and_support":outcomes,
        "interpretation_limits":[
            "Each candidate RAIL plot and ROAD plot was historically selected nonrandomly by different original botanical studies.",
            "Matching original phytogeography, decade/20yr and plot-area bins cannot guarantee species-level habitat/soil/vegetation formation/mowing/traffic equal.",
            "Source vegetation formation not fully harmonized: no causal treatment allocation or direct train freight wind measurement.",
            "Source phytic region/year/area supports are not actual same-side railway BS-BT pairs; do NOT compare coefficients directly with Slovenia.",
            "Original plant species trait-known subset incomplete, and therophyte life form only a proxy for strict annual lifespan.",
            "Stratification diagnostic posthoc; p-values not calculated. No globally representative certainty.",
            "Source species names and full source raw records not redistributed. Any substantial Pladias paper should consult Governing Board.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("zenodo-record","header-csv","species-csv","dispersal-xlsx",
              "lifeform-xlsx","verified119","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    source=json.loads(Path(a.zenodo_record).read_text(encoding="utf-8"))
    last=json.loads(Path(a.verified119).read_text(encoding="utf-8"))
    out=run(source,a.header_csv,a.species_csv,a.dispersal_xlsx,
            a.lifeform_xlsx,last)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf8")
    print(json.dumps({
        "complete":True,"version":"v12.1",
        "n_original_source_veg_rows":out[
            "original_joint_trait_source_sample_matching"]["n_original_source_species_csv_records_scanned"],
        "n_original_source_plots":out["original_full_plot_audit"]["counts"]["original_header_rows"],
        "header_context":out["original_temporal_and_area_candidate_audit"],
        "strategy_counts":out["published_source_checks"],
        "results":{
            label:{key:{
                "balance":d["original_source_sample_restriction"],
                "raw_on_same_plots":d["original_pooled_prevalence_on_EXACT_matched_plots_BEFORE_block_standardizing"],
                "standardized":{k:v for k,v in d[
                   "within_original_source_region_time_area_block_standardized_trait_comparison"].items()
                  if k!="original_source_known_traits_by_lifespan_and_diaspore"}
            }for key,d in blocks.items()}
            for label,blocks in out["results_by_strict_vs_broad_source_railway_label_and_support"].items()
        }
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()

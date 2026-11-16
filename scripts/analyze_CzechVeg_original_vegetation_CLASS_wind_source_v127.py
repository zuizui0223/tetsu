#!/usr/bin/env python3
"""v12.7 source-authenticated CzechVeg vegetation CLASS-matched wind trait.

Original external authors' CzechVeg-Open v1.1 2026 Zenodo source
117739 vegetation headers and 2510340 source species×plot records,
original Pladias Sádlo (2018) seed dispersal strategy / life form
on exact plant species names (no genus guesses or trait imputation).
Original source 2-letter Czech syntaxon class prefix validated against
standalone original source class codes AND official Pladias vegetation
class hierarchy (e.g. XD, XDA, XDA01).

Compare RAIL-text vs ROAD-text vegetation plots within same ORIGINAL
publication, source phytoregion, decade, area class, and common
original vegetation CLASS. This is a *conditional descriptive study*
and NOT the total causal railway effect: vegetation class can itself be
affected by railway management and may be a mediator/collider.

Six scenarios: strict/broad source rail text classification times
common source class with matching conditions:
C1 source study+region+decade+area_class+vegetation_class
C2 source study+region+decade+vegetation_class
C3 source study+region+vegetation_class.
Requirements 1/2/3 source plots per rail and road group per block.
No original Slovenian matched BS vs BT equivalence, no true
freight-borne viable seed, train slipstream, mowing, adaptation
or flowering data.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path

from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    validate_original_joins,block_standardized_guild,area_key
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    plot_labels,plot_taxa,guild_panel,
)
from audit_CzechVeg_syntaxon_same_study_source_v124 import (
    source_meta_with_taxonomic_vegetation,
)
from audit_CzechVeg_original_vegetation_class_v126 import (
    canonical_roots_from_original_full_publisher_csv,
    convert_source_syntaxon_meta,
)
from analyze_CzechVeg_syntaxon_matched_rail_road_traits_v125 import (
    extract_estimate,source_study_cluster_resample,source_leave_one_publication_out
)

CLASS_CONDITIONS=(
    ("C1_same_paper_region_decade_area_and_vegetation_class",True,True),
    ("C2_same_paper_region_decade_vegetation_class",True,False),
    ("C3_same_paper_region_vegetation_class",False,False),
)
THRESHOLDS=(1,2,3)


def class_key(record,condition):
    _,decade,use_area=condition
    name=record["study"];region=record["source_region"]
    y=record["survey_year"]//10*10 if decade else None
    ac=area_key(record["area"],"class4") if use_area else None
    cc=record["SOURCE_CLASS_CODE"]
    if cc is None:raise ValueError("Original vegetation class unknown cannot be assigned")
    return (name,region,y,ac,cc)


def selected_source_class_blocks(source,plot_taxa,mode,condition,min_sites):
    out=defaultdict(lambda:{"rail":[],"road":[]})
    stats=Counter()
    for key,site in sorted(source.items()):
        lab=site[mode]
        if lab not in ("rail","road"):continue
        stats["original_source_candidate_"+lab]+=1
        if not plot_taxa.get(key):
            stats["source_"+lab+"_plot_no_recorded_vascular_taxa"]+=1
            continue
        if site["SOURCE_CLASS_CODE"] is None:
            stats["source_"+lab+"_plot_class_unknown"]+=1
            continue
        out[class_key(site,condition)][lab].append(key)
    common={k:v for k,v in out.items() if
            min(len(v["rail"]),len(v["road"]))>=min_sites}
    return common,{
        "original_source_class_matched_design":condition[0],
        "n_rail_plots_before_class_and_common_study_restrictions":
            stats["original_source_candidate_rail"],
        "n_road_plots_before_class_and_common_study_restrictions":
            stats["original_source_candidate_road"],
        "n_original_RAIL_source_plots_with_unknown_class":
            stats["source_rail_plot_class_unknown"],
        "n_original_ROAD_source_plots_with_unknown_class":
            stats["source_road_plot_class_unknown"],
        "n_source_class_matched_original_publication_study_blocks":len(common),
        "n_original_rail_matched_plots":len({x for v in common.values() for x in v["rail"]}),
        "n_original_road_matched_plots":len({x for v in common.values() for x in v["road"]}),
        "n_original_studies":len({k[0] for k in common}),
        "n_original_phyto_regions":len({k[1] for k in common}),
        "n_original_vegetation_class_codes":len({k[-1] for k in common}),
        "effective_original_source_rail_road_min_plot_pair_mass":sum(
            min(len(v["rail"]),len(v["road"])) for v in common.values()),
        "minimum_observed_original_rail_and_road_plots_each_block":min_sites
    }


def source_class_panel(source,plant_taxa,mode,condition,min_):
    blocks,frame=selected_source_class_blocks(source,plant_taxa,mode,condition,min_)
    selected={p for v in blocks.values() for arm in ("rail","road") for p in v[arm]}
    base=guild_panel(source,plant_taxa,mode,selected)
    estimate=block_standardized_guild(blocks,plant_taxa)
    out={
        "original_source_study_and_canonical_class_support":frame,
        "original_matched_source_plant_diaspore_response":extract_estimate(estimate),
        "original_plot_pooled_on_same_class_matched_source_sites":{
            "identified":base["identified"],
            "source_equal":base.get("strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"),
            "source_occ_weight":base.get("strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap"),
        },
    }
    if condition==CLASS_CONDITIONS[0] and min_==2:
        out["source_publication_cluster_resampling"] = source_study_cluster_resample(
            blocks,plant_taxa,seed=20261114,draws=1200)
        out["original_study_leave_one_out"] = source_leave_one_publication_out(
            blocks,plant_taxa)
    return out


def run(meta,orig_header,orig_species,disp_file,life_file,prev119,prev125,prev126):
    if not(prev125.get("complete") is True and prev125.get("version")=="v12.5"):
        raise ValueError("Require actual original v12.5 same syntaxon source")
    if not(prev126.get("complete") is True and prev126.get("version")=="v12.6"):
        raise ValueError("Require actual original v12.6 anchored source class audit")
    disp,life,proof=validate_original_joins(meta,orig_header,orig_species,
                                            disp_file,life_file,prev119)
    roots,rootqc=canonical_roots_from_original_full_publisher_csv(orig_header)
    if len(roots)!=prev126[
        "original_published_standalone_two_letter_class_root_audit"][
          "n_independent_source_class_code_root_labels"]:
        raise ValueError("Authentic Czech source vegetation class taxonomy root count drift")
    plots,plotqc=plot_labels(orig_header)
    original,exactqc=source_meta_with_taxonomic_vegetation(orig_header,plots)
    classmeta,classqc=convert_source_syntaxon_meta(original,roots)
    if classqc["candidate_with_verified_original_source_class_code"]!=prev126[
        "source_railroad_candidate_class_availability"][
            "candidate_with_verified_original_source_class_code"]:
        raise ValueError("Real v12.6 source class presence/validity drift")
    plant_taxa,original_species_qc=plot_taxa(orig_species,classmeta,disp,life)
    if original_species_qc["n_original_source_species_csv_records_scanned"]!=2510340:
        raise ValueError("Original Czech 2510340 species source records not fully scanned")
    tests={}
    checks=[]
    for mode in ("strict","broad"):
        tests[mode]={}
        for condition in CLASS_CONDITIONS:
            tests[mode][condition[0]]={}
            for minimum in THRESHOLDS:
                v=source_class_panel(classmeta,plant_taxa,mode,condition,minimum)
                tests[mode][condition[0]][str(minimum)]=v
                # Compare exact source v12.6 plot-class common-support counts,
                # except selected site/plant rows may be missing. Thus actual
                # biological samples must never exceed the header-only audit.
                sourceblock=prev126[
                    "original_same_paper_region_year_and_VERIFIED_class_rail_road"][
                    mode][{"C1":"same_study_region_decade_area_class",
                           "C2":"same_study_region_decade_class",
                           "C3":"same_study_region_class"}[condition[0][:2]]][str(minimum)]
                observed=v["original_source_study_and_canonical_class_support"]
                for x,y in (
                  ("n_original_rail_matched_plots",
                   "n_original_rail_plots_on_class_comparable_support"),
                  ("n_original_road_matched_plots",
                   "n_original_road_plots_on_class_comparable_support")
                ):
                    if observed[x]>sourceblock[y]:
                        raise ValueError("Original 117739 syntaxon class source plot coverage exceeds source header audit")
                checks.append({"group":mode,"condition":condition[0],
                               "minimum":minimum,"matches_v126_source_upper_bound":True})
    return {
        "complete":True,"version":"v12.7",
        "original_botanical_data_DOI":"10.5281/zenodo.19231152",
        "original_plant_strategy_source_DOI":"10.23855/preslia.2021.001",
        "verified_original_header_species_and_trait_sources":proof,
        "root_original_vegetation_class_provenance":rootqc,
        "original_source_class_match_coverage":classqc,
        "original_source_species_trait_match_coverage":original_species_qc,
        "provenance_checks_against_original_v126":checks,
        "original_RAIL_text_vs_ROAD_text_same_publication_vegetation_class_matched_plant_trait_results":tests,
        "notes":[
            "This is a CzechVeg original botanical study-site text comparison RAIL vs ROAD and is NOT the Slovenian matched BS−BT within railway design.",
            "Source vegetation class code is the canonical two-letter class prefix only when independently anchored in the source national standalone class code inventory, never inferred for unknown/malformed source codes.",
            "Vegetation class can itself respond to railway disturbance: conditioning on class can discard biologically mediated railway effects and is not a total causal railway effect.",
            "Even same published study may contain heterogeneous original plot collection methods/soil/climate/railway right-of-way site types.",
            "Pladias therophyte is life-form proxy not exactly annual lifespan; species strategy class is trait archetype not observation of actual wind-borne seeds.",
            "No original source train-moved viable seed, train slipstream gust, measured grass cutting, plant biomass/flowering/genetic adaptation data.",
            "No raw original author vegetation CSV nor Pladias source tables redistributed.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    for x in ("zenodo-record","header-csv","species-csv","dispersal-xlsx",
              "lifeform-xlsx","verified119","verified125","verified126","output"):
        p.add_argument("--"+x,required=True)
    a=p.parse_args()
    result=run(
        json.loads(Path(a.zenodo_record).read_text(encoding="utf-8")),
        a.header_csv,a.species_csv,a.dispersal_xlsx,a.lifeform_xlsx,
        json.loads(Path(a.verified119).read_text(encoding="utf-8")),
        json.loads(Path(a.verified125).read_text(encoding="utf-8")),
        json.loads(Path(a.verified126).read_text(encoding="utf-8")))
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    short={}
    for mode,groups in result["original_RAIL_text_vs_ROAD_text_same_publication_vegetation_class_matched_plant_trait_results"].items():
        short[mode]={}
        for name,d in groups.items():
            short[mode][name]={}
            for threshold,x in d.items():
                b=x["original_source_study_and_canonical_class_support"]
                e=x["original_matched_source_plant_diaspore_response"]
                short[mode][name][threshold]={
                    "blocks":b["n_source_class_matched_original_publication_study_blocks"],
                    "studies":b["n_original_studies"],
                    "rail":b["n_original_rail_matched_plots"],
                    "road":b["n_original_road_matched_plots"],
                    "classes":b["n_original_vegetation_class_codes"],
                    "estimate":e,
                    "bootstrap":x.get("source_publication_cluster_resampling"),
                    "n_study_leave_one_out":len(x.get("original_study_leave_one_out",[]))
                }
    print(json.dumps({"complete":True,"version":"v12.7","results":short},ensure_ascii=False),flush=True)


if __name__=="__main__":main()

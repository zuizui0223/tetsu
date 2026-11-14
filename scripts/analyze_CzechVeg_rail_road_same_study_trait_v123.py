#!/usr/bin/env python3
"""v12.3: original CzechVeg same-PUBLICATION railway-road plant-source
common-study x region x year x plot area conservative diaspore check.

NOT a Slovenia BS-v-BT source replication, NOT a randomized comparison.
Original 117739 source vegetation plots, 2510340 original species×plot
records from Chytrý et al. CzechVeg-Open 2026, original published source
DOI 10.5281/zenodo.19231152; Pladias original Sádlo et al. 2018.
Studies have heterogeneous botanical plot designs even within a paper.
Reference-title shared source is only a conservative design-control proxy.

Uses exact official source checksums and full Pladias source biological
inventory gate from v11.9 and source-style v12.1 six-cohort original
field data. Reports inability to estimate when <3 taxa per strategy and
therophyte group; never fills missing plant traits by genus inference.
No original source raw redistributed, no wind/carriage/mowing biology
causally inferred.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    validate_original_joins,header_source_context,area_key,
    block_standardized_guild,
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    plot_labels,plot_taxa,guild_panel,
)
from audit_CzechVeg_shared_railroad_reference_studies_v122 import source_title


SPECS=(
    ("S1_study_region_decade_area4_min1",10,"class4",1,True),
    ("S2_study_region_decade_area4_min2",10,"class4",2,True),
    ("S3_study_region_decade_no_area_min2",10,"none",2,True),
    ("S4_study_region_area4_no_year_min2",0,"class4",2,True),
    ("S5_study_region_only_min2",0,"none",2,True),
    ("S6_study_only_min3",0,"none",3,False),
)


def source_study_meta(original_header,sites):
    out={}
    source_stats=Counter()
    with open(original_header,encoding="utf-8-sig",newline="") as file:
        reader=csv.DictReader(file)
        if "REFERENCE_TITLE" not in (reader.fieldnames or []):
            raise ValueError("Original Czech source bibliographic study key absent")
        n=0
        for row in reader:
            n+=1
            plot=str(row.get("RELEVE_NR") or "").strip()
            if plot not in sites:continue
            s=source_title(row.get("REFERENCE_TITLE"))
            if s is None:
                source_stats["source_original_plot_no_meaningful_reference"]+=1
                continue
            out[plot]={**sites[plot],"original_reference_title":s}
            source_stats["source_original_valid_area_year_region_and_study"]+=1
    if n!=117739:raise ValueError("Original full CzechVeg published plot headers incomplete")
    return out,dict(source_stats)


def source_study_block_key(meta,spec):
    title,yearspan,area_mode,minimum,geographic=spec
    study=meta["original_reference_title"]
    region=meta["source_phyto_region"] if geographic else None
    dec=(meta["source_year_int"]//yearspan*yearspan) if yearspan else None
    area=area_key(meta["area"],area_mode)
    return (study,region,dec,area)


def common_original_study_blocks(meta,source_taxa,label,spec):
    kname,_,_,minimum,_=spec
    block=defaultdict(lambda:{"rail":[],"road":[]})
    candidate=Counter()
    for plot,m in sorted(meta.items()):
        group=m[label]
        if group not in ("rail","road"):continue
        if not source_taxa.get(plot):
            candidate["source_plot_without_original_vascular_records"]+=1
            continue
        candidate[group]+=1
        key=source_study_block_key(m,spec)
        block[key][group].append(plot)
    common={k:v for k,v in block.items()
            if len(v["rail"])>=minimum and len(v["road"])>=minimum}
    report={
        "predeclared_source_study_region_year_area_definition":kname,
        "source_rail_candidate_with_valid_paper_and_real_vascular_plants":candidate["rail"],
        "source_road_candidate_with_valid_paper_and_real_vascular_plants":candidate["road"],
        "candidate_plot_without_original_vascular_rows_excluded":
            candidate["source_plot_without_original_vascular_records"],
        "n_actual_original_published_study_based_source_blocks":len(block),
        "n_original_study_blocks_with_rail_and_road_common_support":len(common),
        "n_distinct_original_published_studies_in_common_support":len({x[0] for x in common}),
        "n_source_RAIL_plots_in_common_support":sum(len(x["rail"]) for x in common.values()),
        "n_source_ROAD_plots_in_common_support":sum(len(x["road"]) for x in common.values()),
        "n_original_rail_sites_removed_by_study_and_design_gates":
            candidate["rail"]-sum(len(x["rail"]) for x in common.values()),
        "n_original_road_sites_removed_by_study_and_design_gates":
            candidate["road"]-sum(len(x["road"]) for x in common.values()),
        "n_region_codes_within_shared_original_study_blocks":
            len({x[1] for x in common if x[1] is not None}),
        "effective_balanced_source_site_pair_count":
            sum(min(len(v["rail"]),len(v["road"])) for v in common.values()),
        "n_study_block_source_plot_count_examples_no_reference_titles":sorted([
            {"n_rail":len(v["rail"]),"n_road":len(v["road"]),
             "area_band":str(k[3]),"decade_start":k[2],
             "region":k[1]} for k,v in common.items()],
            key=lambda x:-min(x["n_rail"],x["n_road"]))[:12]
    }
    return common,report


def comparison(meta,source_taxa,label,spec):
    blocks,balance=common_original_study_blocks(meta,source_taxa,label,spec)
    matched={x for cells in blocks.values() for x in cells["rail"]+cells["road"]}
    raw=guild_panel(meta,source_taxa,label,matched)
    standard=block_standardized_guild(blocks,source_taxa)
    return {
        "source_common_support":balance,
        "pooled_on_same_actual_selected_study_plots":{
            "identified":raw["identified"],
            "n_source_unique_joint_trait_taxa":
                raw.get("n_source_unique_joint_trait_taxa_in_original_rail_or_road_vegetation"),
            "equal":raw.get("strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"),
            "occ_weighted":raw.get("strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap"),
        },
        "source_same_publication_geo_year_area_standardized":{
            k:v for k,v in standard.items()
            if k!="original_source_known_traits_by_lifespan_and_diaspore"},
        "no_railway_mechanical_causal_inference":True
    }


def run(meta,header_csv,species_csv,disp_source,life_source,prev119):
    disp,life,proof=validate_original_joins(meta,header_csv,species_csv,
                                             disp_source,life_source,prev119)
    plots,plotquality=plot_labels(header_csv)
    with_year,yearquality=header_source_context(header_csv,plots)
    m,studyquality=source_study_meta(header_csv,with_year)
    source_taxa,traitsquality=plot_taxa(species_csv,m,disp,life)
    if traitsquality["n_original_source_species_csv_records_scanned"]!=2510340:
        raise ValueError("Original 2510340 Czech source species-by-plot rows not recovered")
    panels={}
    for label in ("strict","broad"):
        panels[label]={}
        for spec in SPECS:
            panels[label][spec[0]]=comparison(m,source_taxa,label,spec)
    return {
        "complete":True,"version":"v12.3",
        "original_CzechVeg_DOI":"10.5281/zenodo.19231152",
        "original_Pladias_DOI":"10.23855/preslia.2021.001",
        "source_provenance":proof,
        "original_published_study_source_quality":studyquality,
        "original_117739_survey_header_quality":yearquality,
        "original_2510340_source_species_quality":traitsquality,
        "source_strict_and_broad_PUBLISHED_study_overlap_6_designs":panels,
        "scientific_caveats":[
            "Same REFERENCE_TITLE can still include fundamentally different rail vs road vegetation habitats and original sampling methods.",
            "Source studies/regions with too few wind, annual or other taxa are NOT a biological null effect.",
            "Original Czech railway and roadside plot positions are word matches not externally confirmed inside right-of-way.",
            "Plots not randomized; no effect of actual train seed transport, mowing, wind stress, phenology or genetic adaptation measured.",
            "The original Czech RAIL-text/ROAD-text estimand is NOT Slovenia's original within-field RAIL BS-versus-BT paired plot contrast.",
            "No original published full raw source or Pladias trait table redistributed in output.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for field in ("zenodo-record","header-csv","species-csv","dispersal-xlsx",
                  "lifeform-xlsx","verified119","output"):
        p.add_argument("--"+field,required=True)
    a=p.parse_args()
    res=run(json.loads(Path(a.zenodo_record).read_text(encoding="utf8")),
            a.header_csv,a.species_csv,a.dispersal_xlsx,a.lifeform_xlsx,
            json.loads(Path(a.verified119).read_text(encoding="utf8")))
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({
        "complete":True,"version":"v12.3",
        "published_source_study_quality":res["original_published_study_source_quality"],
        "original_species_rows":res["original_2510340_source_species_quality"][
            "n_original_source_species_csv_records_scanned"],
        "result":res["source_strict_and_broad_PUBLISHED_study_overlap_6_designs"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()

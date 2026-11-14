#!/usr/bin/env python3
"""v12.5 TRUE CzechVeg source vegetation-syntaxon-matched diaspore comparison.

Original external authors' CzechVeg-Open v1.1 2026, full 117739
releve metadata and 2510340 plant/releve presences from official
Zenodo 10.5281/zenodo.19231152; original Sádlo et al. (2018)
Pladias dispersal archetypes and Raunkiaer life forms. NOT tetsu
original field collection. This is RAIL-TEXT vs ROAD-TEXT across
distinct botanic plots, NOT within-RAIL Slovenian BS-vs-BT data.

Use original publication title, phytogeographic region, original year
decade, original plot area category, AND original source SYNTAXON
code equality without invented syntaxonomic class conversion.
Report strict vs broad site-text original two habitat definitions,
min 1 and 2 source original plots per arm per common-source block.
Retain original no-vegetation-code baseline and verify exact v12.3
numerical source estimates and sample counts.

This source-controlled analysis CONDITIONS on vegetation type, which
may itself be altered by railway management, so the resulting contrast
is a restricted composition comparison, NOT a total rail effect.
No general causal physical seed dispersal/adaptation claims.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    validate_original_joins,block_standardized_guild,
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    plot_labels,plot_taxa,guild_panel,
)
from audit_CzechVeg_syntaxon_same_study_source_v124 import (
    source_meta_with_taxonomic_vegetation,block_key,MATCH_LEVELS,
)
from analyze_CzechVeg_rail_road_region_robustness_v120 import validated_source


LEVELS=("same_paper_phyto_decade_area_no_syntax",
        "same_paper_phyto_decade_area_SYNTAXON",
        "same_paper_phyto_decade_area_SYNTAXON_NAME",
        "same_paper_phyto_decade_area_SYNTAX_OLD")
MINIMUMS=(1,2)
ROOT_RNG=20261114
STUDY_DRAWS=1200


def source_blocks(meta,taxa,mode,level,minimum):
    source_level=next(x for x in MATCH_LEVELS if x[0]==level)
    groups=defaultdict(lambda:{"rail":[],"road":[]})
    qc=Counter()
    for pid,m in sorted(meta.items()):
        arm=m[mode]
        if arm not in ("rail","road"):continue
        qc["source_"+arm+"_candidates"]+=1
        if not taxa.get(pid):
            qc["missing_original_vascular_taxa_"+arm]+=1
            continue
        if source_level[4] is not None and m[source_level[4]] is None:
            qc["missing_original_syntaxon_"+arm]+=1
            continue
        code=block_key(m,source_level)
        groups[code][arm].append(pid)
    valid={k:v for k,v in groups.items() if
           min(len(v["rail"]),len(v["road"]))>=minimum}
    return valid,{
        "original_vegetation_source_field":source_level[4],
        "original_rail_label":mode,
        "minimum_source_plots_per_original_arm":minimum,
        "n_source_blocks_with_original_rail_and_road":len(valid),
        "n_published_source_studies_in_matched_blocks":len({k[0] for k in valid}),
        "n_source_phytoregions_in_matched_blocks":len({k[1] for k in valid}),
        "n_source_rail_plots_matched":sum(len(v["rail"]) for v in valid.values()),
        "n_source_road_plots_matched":sum(len(v["road"]) for v in valid.values()),
        "balanced_original_study_block_mass":sum(
            min(len(v["rail"]),len(v["road"])) for v in valid.values()),
        "n_original_rail_candidates_before_source_gates":qc["source_rail_candidates"],
        "n_original_road_candidates_before_source_gates":qc["source_road_candidates"],
        "n_original_source_rail_plots_missing_syntaxon":qc["missing_original_syntaxon_rail"],
        "n_original_source_road_plots_missing_syntaxon":qc["missing_original_syntaxon_road"],
        "n_original_source_rail_plot_with_no_vascular_rows":
            qc["missing_original_vascular_taxa_rail"],
        "n_original_source_road_plot_with_no_vascular_rows":
            qc["missing_original_vascular_taxa_road"],
        "original_vegetation_source_ranks_not_harmonized":True,
    }


def extract_estimate(result):
    if not result.get("identified"):
        return {"identified":False,"n_source_taxa":result.get(
            "n_original_joint_trait_species_source_common_support"),
            "n_wind_taxa":result.get(
            "n_original_joint_trait_species_with_source_strict_wind_code")}
    return {
        "identified":True,
        "n_source_taxa":result["n_original_joint_trait_species_source_common_support"],
        "n_wind_taxa":result["n_original_joint_trait_species_with_source_strict_wind_code"],
        "species_equal":result["original_source_species_equal_strict_wind_vs_other_prevalence_delta"],
        "occurrence_weighted":result["original_source_occ_weight_strict_wind_vs_other_prevalence_delta"],
        "weighting_gap":result["original_source_weighting_gap_equal_minus_occ"],
        "source_opposite_signs":result["opposite_source_estimand_signs"],
        "fraction_matched_original_source_plant_incidence_with_joint_traits":
            result["fraction_joint_trait_incidence_coverage"],
    }


def numerical_percentile(values,p):
    arr=sorted(values)
    if not arr:return None
    rank=(len(arr)-1)*p
    lo=int(rank);hi=min(lo+1,len(arr)-1)
    return arr[lo]+(arr[hi]-arr[lo])*(rank-lo)


def source_study_cluster_resample(blocks,plot_taxa,seed=ROOT_RNG,draws=STUDY_DRAWS):
    """Cluster on original publication title. Critically tiny samples
    are described and never interpreted as country-population confidence.
    """
    if not blocks:
        return {"n_source_studies":0,"identified":False}
    studies=sorted({k[0] for k in blocks})
    if len(studies)<3:
        return {"n_source_studies":len(studies),"identified":False,
                "cause":"Fewer than three independent original publication study labels"}
    grouped=defaultdict(list)
    for key,v in blocks.items():
        grouped[key[0]].append((key,v))
    rng=random.Random(seed)
    eq=[];freq=[];gap=[]
    nonident=0
    signs=0
    for _ in range(draws):
        sample=[rng.choice(studies) for __ in studies]
        boot={}
        for i,s in enumerate(sample):
            for j,(key,v) in enumerate(grouped[s]):
                # A resampled original study must remain an intact source unit.
                boot[("source_study_resample",i,j,key)]=v
        result=block_standardized_guild(boot,plot_taxa)
        if result.get("identified"):
            a=result["original_source_species_equal_strict_wind_vs_other_prevalence_delta"]
            b=result["original_source_occ_weight_strict_wind_vs_other_prevalence_delta"]
            eq.append(a);freq.append(b);gap.append(a-b)
            signs+=int(a*b<0)
        else:nonident+=1
    return {
        "n_original_published_study_clusters":len(studies),
        "n_source_study_bootstrap_draws_requested":draws,
        "n_source_study_bootstrap_draws_identified":len(eq),
        "n_source_study_bootstrap_draws_nonidentifiable":nonident,
        "95pct_source_study_resampling_equal":(
            [numerical_percentile(eq,.025),numerical_percentile(eq,.975)]
            if eq else None),
        "95pct_source_study_resampling_occ":(
            [numerical_percentile(freq,.025),numerical_percentile(freq,.975)]
            if freq else None),
        "95pct_source_study_resampling_equal_minus_occ_gap":(
            [numerical_percentile(gap,.025),numerical_percentile(gap,.975)]
            if gap else None),
        "fraction_of_identifiable_source_bootstrap_with_opposite_signs":
            signs/len(eq) if eq else None,
        "source_region_and_study_not_randomized_country":True,
    }


def source_leave_one_publication_out(blocks,plot_taxa):
    papers=sorted({key[0] for key in blocks})
    out=[]
    for study in papers:
        reduced={k:v for k,v in blocks.items() if k[0]!=study}
        r=block_standardized_guild(reduced,plot_taxa)
        out.append({
            "omitted_study_title_excerpt":study[:110],
            "remaining_source_studies":len(papers)-1,
            **extract_estimate(r)})
    return out


def panel(meta,taxa,mode,level,minimum):
    blocks,count=source_blocks(meta,taxa,mode,level,minimum)
    matched={v for block in blocks.values() for arm in ("rail","road") for v in block[arm]}
    raw=guild_panel(meta,taxa,mode,matched)
    aligned=block_standardized_guild(blocks,taxa)
    estimate=extract_estimate(aligned)
    result={
        "source_original_sampling_frame":count,
        "source_original_study_area_year_and_syntaxon_standardized":estimate,
        "unstandardized_original_same_sites_composition":{
            "identified":raw["identified"],
            "n_original_trait_known_taxa":raw.get("n_source_unique_joint_trait_taxa_in_original_rail_or_road_vegetation"),
            "species_equal":raw.get("strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"),
            "occurrence_weighted":raw.get("strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap")
        }
    }
    if level=="same_paper_phyto_decade_area_SYNTAXON" and minimum==2:
        result["source_study_resampling_conditional"] = (
            source_study_cluster_resample(blocks,taxa))
        result["leave_one_original_study_out"] = (
            source_leave_one_publication_out(blocks,taxa))
    return result


def run(meta,header_csv,species_csv,disp_source,life_source,prior119,prior123,prior124):
    if prior123.get("complete") is not True or prior123.get("version")!="v12.3":
        raise ValueError("Actual source v12.3 same-study trait result required")
    if prior124.get("complete") is not True or prior124.get("version")!="v12.4":
        raise ValueError("Actual source v12.4 original syntaxon audit receipt required")
    disp,life,quality=validate_original_joins(meta,header_csv,species_csv,
                                                disp_source,life_source,prior119)
    plots,orig_quality=plot_labels(header_csv)
    meta_syntax,syntax_quality=source_meta_with_taxonomic_vegetation(header_csv,plots)
    if (syntax_quality["retained_source_plot_with_original_study_phyto_year_area"]!=
            prior124["original_source_syntaxon_label_audit"][
                "retained_source_plot_with_original_study_phyto_year_area"]):
        raise ValueError("Original Czech Veg source syntaxon 117739 plot cohort drift")
    taxa,source_taxa=plot_taxa(species_csv,meta_syntax,disp,life)
    if source_taxa["n_original_source_species_csv_records_scanned"]!=2510340:
        raise ValueError("Czech original full 2,510,340 plant rows missing")
    main={}
    for mode in ("strict","broad"):
        main[mode]={}
        for level in LEVELS:
            for minimum in MINIMUMS:
                name=level+"__min"+str(minimum)
                main[mode][name]=panel(meta_syntax,taxa,mode,level,minimum)
        base=main[mode]["same_paper_phyto_decade_area_no_syntax__min2"]
        actual=base["source_original_study_area_year_and_syntaxon_standardized"]
        baseline=prior123[
            "source_strict_and_broad_PUBLISHED_study_overlap_6_designs"][mode][
                "S2_study_region_decade_area4_min2"]
        prior=baseline["source_same_publication_geo_year_area_standardized"]
        sites=base["source_original_sampling_frame"]
        olds=baseline["source_common_support"]
        if (sites["n_source_rail_plots_matched"]!=olds["n_source_RAIL_plots_in_common_support"] or
            sites["n_source_road_plots_matched"]!=olds["n_source_ROAD_plots_in_common_support"]):
            raise ValueError("Original v12.3 matched source vegetation site counts do not reproduce")
        if actual["identified"]!=prior["identified"]:
            raise ValueError("Original v12.3 identified source design drift")
        if actual["identified"]:
            for newkey,oldkey in (
              ("species_equal","original_source_species_equal_strict_wind_vs_other_prevalence_delta"),
              ("occurrence_weighted","original_source_occ_weight_strict_wind_vs_other_prevalence_delta")):
                if not math.isclose(actual[newkey],prior[oldkey],rel_tol=0,abs_tol=1e-12):
                    raise ValueError("Original v12.3 source numerical wind trait mismatch: "+newkey)
    return {
        "complete":True,"version":"v12.5",
        "original_published_CzechVeg_DOI":"10.5281/zenodo.19231152",
        "original_external_Pladias_DOI":"10.23855/preslia.2021.001",
        "original_full_source_integrity":quality,
        "original_syntaxon_national_source_117739_plot_quality":syntax_quality,
        "original_2510340_source_species_trait_quality":source_taxa,
        "two_source_site_label_definitions_four_original_syntaxon_models":{
            "strict":main["strict"],"broad":main["broad"]},
        "matches_original_v12_3_source_biological_effects_exactly":True,
        "scientific_noncausal_boundary":[
            "No original strict HABITAT RAIL and ROAD plot shares a published study, region, decade, plot area class AND original SYNTAXON across source v12.4; such a biological comparison is structurally unidentifiable.",
            "Broader source habitat/locality labels may describe areas adjacent to but not within actual railway operational land, and source syntaxonomic codes span different ranks.",
            "Conditioning on vegetation type can block an effect of rail exposure itself, so this is NOT total effect of railway management; a survivor-subset compositional diagnostic only.",
            "Original source publication title does not identify independently replicated random study designs.",
            "No real train-induced seed movement, plant establishment fitness, mowing exposure or wind-speed measurement exists in these source files.",
            "Do not generalize selected 8-study same-syntaxon candidate comparisons to a worldwide railway dispersal strategy effect.",
            "Original published botany/Pladias database source provider citation and substantial reuse governance must be respected.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("zenodo-record","header-csv","species-csv","dispersal-xlsx",
                 "lifeform-xlsx","verified119","verified123","verified124","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    result=run(
        json.loads(Path(a.zenodo_record).read_text(encoding="utf8")),
        a.header_csv,a.species_csv,a.dispersal_xlsx,a.lifeform_xlsx,
        json.loads(Path(a.verified119).read_text(encoding="utf8")),
        json.loads(Path(a.verified123).read_text(encoding="utf8")),
        json.loads(Path(a.verified124).read_text(encoding="utf8")))
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    compact={}
    for mode,v in result["two_source_site_label_definitions_four_original_syntaxon_models"].items():
        compact[mode]={}
        for name,info in v.items():
            count=info["source_original_sampling_frame"]
            effect=info["source_original_study_area_year_and_syntaxon_standardized"]
            compact[mode][name]={
                "blocks":count["n_source_blocks_with_original_rail_and_road"],
                "studies":count["n_published_source_studies_in_matched_blocks"],
                "rail":count["n_source_rail_plots_matched"],
                "road":count["n_source_road_plots_matched"],
                "balanced":count["balanced_original_study_block_mass"],
                "effect":effect,
                "bootstrap":info.get("source_study_resampling_conditional"),
                "loo_count":len(info.get("leave_one_original_study_out",[])),
                "loo_sign_reversal":sum(
                    x.get("source_opposite_signs",False)
                    for x in info.get("leave_one_original_study_out",[]))
            }
    print(json.dumps({
        "complete":True,"version":"v12.5",
        "verified_117739_source_plots":True,
        "verified_2510340_source_rows":True,
        "by_mode":compact
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

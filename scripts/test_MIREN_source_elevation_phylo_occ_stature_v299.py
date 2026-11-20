#!/usr/bin/env python3
"""v29.9 MIREN real environmental-gradient mechanism candidate, author ORIGINAL
Elevation in metres, original field site/year and geography, trait-labelled
species on genuine V.PhyloMaker2 75%-clades x original occupancy quintile.

Outcomes:
  A. primary predeclared: within ORIGINAL ROAD, slope of community-weighted
     z(log canopy height) vs native site Elevation (per 1000 metres).
  B. secondary: within original road, slope of abs between-site height-CWM
     distance vs abs original elevation difference (per 1000 metres).

Both slopes are original-plot observations minus paired 2D trait vector
permutation expectation WITHIN real clade75 x species occurrence fifth.
MIREN ROADS and adjacent INLAND vegetation are a single field programme.
This is observational environmental association, not railway adaptation,
not full climate causality, and no field p-value from technical shuffles.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from statistics import mean
from collections import defaultdict
import numpy as np

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 import (
 PANEL_DEFINITIONS,PINNED_PRIOR_ROAD_RESULT
)
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
 read_source_miren,regional_comparison)
from validate_MIREN_original_roads_true_phylogeny_clade_negative_control_v297 import (
 trait_pool_original_matched_regions,construct_original_region_tree_panel,
 run_region_null,MODES as V297_MODES)
from analyze_MIREN_roads_inland_joint_phylo75_occupancy_null_v298 import (
 source_study_incidence,JOINT_MODE
)
from phylogeny_occupancy_joint_label_null_v298 import (
 original_global_taxon_occupancy_quintiles,
 joint_true_clade_and_occupancy_groups,original_source_joint_exchangeability)
from audit_VPhyloMaker2_RData_python_no_R_environment_v296 import (
 inspect_tree_rdata,TREES,COMMIT
)
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import (
 original_root_depths,original_true_clade_groups,permute_paired_traits,source_json_native)
from audit_MIREN_original_native_elevation_source_v299 import (
 read_original_elevation_csv,original_paired_site_elevation,
 MAX_PAIRED_ROADSIDE_INLAND_ALTITUDE_DIFFERENCE_M,
 MIN_ORIGINAL_ROAD_ELEVATION_SPAN_M)
 
ORIGINAL_PANEL_DESCRIPTORS={
 "primary_5species_25percent":"primary_original_5_species_25percent_trait_source",
 "strict_10species_50percent":"strict_original_10_species_50percent_trait_source",
}
SEEDS=(299201,299202)
DRAWS_EACH=199
MIN_WITHIN_ROAD_ALTITUDE_SPAN_M=150.
MIN_PAIRED_SITES_PER_ROAD=4
MIN_ROADS_PER_REGION=2
MIN_SOURCE_ORIGINAL_PHYLO_SITE_RETENTION=.80
MIN_ALTITUDE_YEAR_CONSISTENT_SUBSET_RETENTION=.60
MIN_REGIONS_FOR_DESCRIPTIVE_MEAN=3
ZONES=("roadside","interior","roadside_minus_interior")
METRICS=("height_signed_slope_per_km","height_distance_vs_elevation_gap_slope_per_km")

def _observed_plant_taxa(W:np.ndarray,names:list[str],row:int)->set[str]:
    return {names[int(x)] for x in np.flatnonzero(W[row])}

def construct_native_environment_cohort(original_roads:list,reg:str,source_phylo:dict,
                                        original_altitude:dict,min_species:int)->tuple[list,dict]:
    """Match original author site id/YEAR to row of immutable phylo-tip W.
    Prevent silent reordering or replacement of original source field units.
    """
    phylo_names=source_phylo["species"]
    source_name_set=set(phylo_names)
    original_altitude_road_list=[]
    source_road_index=0
    original_eligible_sites=sum(len(d["roadside"]) for d in source_phylo["roads"])
    all_original_rows=sum(len(rows) for road,rows in original_roads)
    failed_original_elevation=0
    n_paired_nonmatched=0
    n_source_mixed_year_roads=0
    for road_name,rows in original_roads:
        # A row is retained by original source phylogenetic study if both
        # habitats contain >= threshold existing ORIGINAL backbone tip taxa.
        # Source phylo pool includes only retained roads and may exclude some
        # single-road taxa; for site match the true current pool works only
        # if no new taxa occur on a different discarded road. Use source W
        # as the final source species universe; if true-site row drops, fail.
        elig=[row for row in rows if min(
            len(row[z]&source_name_set) for z in ("roadside","interior"))>=min_species]
        if len(elig)<MIN_PAIRED_SITES_PER_ROAD:
            continue
        if source_road_index>=len(source_phylo["roads"]):
            raise ValueError("Original source road count vs true phylogenetic retained cohort drift")
        W=source_phylo["roads"][source_road_index]
        source_road_index+=1
        if len(W["roadside"])!=len(elig) or len(W["interior"])!=len(elig):
            raise ValueError("Original source first-year paired sites vs phylogenetic W row order drift")
        if len({int(row["year"]) for row in elig})>1:
            n_source_mixed_year_roads+=1
        keep=[]
        for j,row in enumerate(elig):
            for z in ("roadside","interior"):
                actual=_observed_plant_taxa(W[z],phylo_names,j)
                expected=row[z]&source_name_set
                if actual!=expected:
                    raise ValueError("Original source site-by-species incidence misaligned during elevation join")
            alt=original_paired_site_elevation(row,reg,road_name,original_altitude)
            if alt is None:
                failed_original_elevation+=1
                continue
            keep.append((j,float(alt),int(row["year"])))
        if len(keep)<MIN_PAIRED_SITES_PER_ROAD:
            n_paired_nonmatched+=len(elig)
            continue
        site_alt=np.asarray([alt for j,alt,y in keep],dtype=float)
        if np.ptp(site_alt)<MIN_WITHIN_ROAD_ALTITUDE_SPAN_M:
            n_paired_nonmatched+=len(elig)
            continue
        selected_indices=np.asarray([j for j,alt,y in keep],dtype=int)
        years=[y for j,alt,y in keep]
        road_matrices={z:W[z][selected_indices] for z in ("roadside","interior")}
        original_altitude_road_list.append({
            "elev_m":site_alt,
            "source_original_years":years,
            "matrices":road_matrices,
            "n_phylo_paired_original_sites_in_road":len(elig),
        })
    if source_road_index!=len(source_phylo["roads"]):
        raise ValueError("An original source valid phylo road was lost in new native elevation pairing")
    kept=sum(len(r["elev_m"]) for r in original_altitude_road_list)
    frac=kept/original_eligible_sites if original_eligible_sites else 0.
    selected=(len(original_altitude_road_list)>=MIN_ROADS_PER_REGION and
              kept>=MIN_PAIRED_SITES_PER_ROAD*MIN_ROADS_PER_REGION and
              frac>=MIN_SOURCE_ORIGINAL_PHYLO_SITE_RETENTION)
    QA={
       "source_original_phylotip_first_year_paired_field_sites":original_eligible_sites,
       "original_field_sites_before_phylogeny_filter":all_original_rows,
       "original_source_roads_with_4_sites_150m_valid_elevation":
         len(original_altitude_road_list),
       "source_true_field_sites_on_usable_elevation_roads":kept,
       "fraction_original_phylotip_matched_sites_retained_for_elevation":frac,
       "source_original_phylo_roads_with_mixed_first_eligible_survey_year":
         n_source_mixed_year_roads,
       "source_plant_site_pair_metadata_missing_elevation_n":failed_original_elevation,
       "source_original_paired_sites_on_unsupported_elevation_roads_n":n_paired_nonmatched,
       "original_environmental_altitude_support":bool(selected),
       "predeclared_min_altitude_span_m":MIN_WITHIN_ROAD_ALTITUDE_SPAN_M,
       "predeclared_min_original_road_count":MIN_ROADS_PER_REGION,
       "predeclared_min_original_sites_per_road":MIN_PAIRED_SITES_PER_ROAD,
       "predeclared_min_phylo_original_site_retention":MIN_SOURCE_ORIGINAL_PHYLO_SITE_RETENTION,
       "original_site_Elevation_m_and_GPS_NEVER_exported":True,
    }
    return original_altitude_road_list,QA

def road_environment_slopes(W:np.ndarray,alt_m:np.ndarray,T:np.ndarray)->dict:
    if W.ndim!=2 or W.shape[1]!=len(T) or len(W)!=len(alt_m) or len(alt_m)<4:
        raise ValueError("Original site W and author-provided native site Elevation misaligned")
    if not np.all(np.isfinite(alt_m)) or np.ptp(alt_m)<MIN_WITHIN_ROAD_ALTITUDE_SPAN_M:
        raise ValueError("Insufficient original source environmental altitude gradient")
    H=W@T[:,0]
    if not np.isfinite(H).all():raise ValueError("Unusable original LEDA field CWM height")
    e=alt_m/1000.
    d=e-e.mean()
    signed=float(np.dot(d,H-H.mean())/np.dot(d,d))
    i,j=np.triu_indices(len(alt_m),1)
    elev_gap=np.abs(e[i]-e[j])
    trait_gap=np.abs(H[i]-H[j])
    x=elev_gap-elev_gap.mean()
    den=float(np.dot(x,x))
    if den<=1e-9:
        raise ValueError("Original between-site pairwise elevation gaps do not vary")
    gap_slope=float(np.dot(x,trait_gap-trait_gap.mean())/den)
    return {"height_signed_slope_per_km":signed,
            "height_distance_vs_elevation_gap_slope_per_km":gap_slope}

def region_environment_source_scores(roads:list,T:np.ndarray)->dict:
    by={z:{k:[] for k in METRICS} for z in ("roadside","interior")}
    if len(roads)<MIN_ROADS_PER_REGION:
        raise ValueError("Need >=2 independent original source roads for region environmental test")
    for road in roads:
        for z in ("roadside","interior"):
            out=road_environment_slopes(road["matrices"][z],road["elev_m"],T)
            for metric,val in out.items():by[z][metric].append(val)
    scores={z:{metric:mean(arr) for metric,arr in vals.items()}
            for z,vals in by.items()}
    scores["roadside_minus_interior"]={
        metric:scores["roadside"][metric]-scores["interior"][metric]
        for metric in METRICS}
    return scores

def label_null_environment_slopes(roads:list,T:np.ndarray,groups:list[np.ndarray],
                                  tree_index:int,region_index:int,
                                  source_scenario:str,draws:int=DRAWS_EACH)->dict:
    if draws<5:raise ValueError("Need enough technical permutations")
    observed=region_environment_source_scores(roads,T)
    outcomes={z:{m:[] for m in METRICS} for z in ZONES}
    for seed in SEEDS:
        rng=np.random.default_rng(seed+17011*tree_index+5009*region_index+
                                  sum(ord(x) for x in source_scenario))
        for rep in range(draws):
            shuffled=permute_paired_traits(T,groups,rng)
            candidate=region_environment_source_scores(roads,shuffled)
            for zone in ZONES:
                for met in METRICS:
                    outcomes[zone][met].append(candidate[zone][met])
    result={}
    for zone in ZONES:
        result[zone]={}
        for metric,null in outcomes[zone].items():
            actual=observed[zone][metric];ref=mean(null)
            result[zone][metric]={
              "source_original_observed_environmental_slope":actual,
              "source_true_phylo_and_occupancy_trait_null_mean_slope":ref,
              "original_observed_minus_conditional_null_environmental_slope":actual-ref,
              "source_technical_lower_rank_NOT_field_p":
                  (1+sum(v<=actual for v in null))/(1+len(null)),
              "source_technical_upper_rank_NOT_field_p":
                  (1+sum(v>=actual for v in null))/(1+len(null))}
    return result

def original_year_consistent_roads(roads:list)->list:
    return [r for r in roads if len(set(r["source_original_years"]))==1]

def source_environment_summary(by_region:dict)->dict:
    # original source regions, NOT independent botanical papers/countries
    if len(by_region)<MIN_REGIONS_FOR_DESCRIPTIVE_MEAN:
        return {"source_regions_n":len(by_region),"source_status":"INSUFFICIENT_SOURCE_REGIONAL_SUPPORT",
                "region_equal_estimates":None}
    result={}
    for z in ZONES:
        result[z]={}
        for metric in METRICS:
            x={reg:r[z][metric]["original_observed_minus_conditional_null_environmental_slope"]
               for reg,r in by_region.items()}
            obs={reg:r[z][metric]["source_original_observed_environmental_slope"]
                 for reg,r in by_region.items()}
            loo=[mean(v for reg,v in x.items() if reg!=removed) for removed in x]
            result[z][metric]={
                "original_source_region_equal_mean_observed_environmental_slope":mean(obs.values()),
                "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null":mean(x.values()),
                "original_region_codes_positive_excess":sum(v>0 for v in x.values()),
                "original_region_codes_negative_excess":sum(v<0 for v in x.values()),
                "source_original_region_loo_mean_range_not_field_CI":[min(loo),max(loo)],
                "source_original_region_code_aggregate_effects":x,
                "field_real_road_site_pairs_NOT_independent_publications":True}
    return {"source_regions_n":len(by_region),
            "source_status":"SOURCE_SUPPORTED_ORIGINAL_MIREN_REGION_EQUAL_DESCRIPTIVE",
            "region_equal_estimates":result}

def run(original_miren:Path,height:Path,seed:Path,tree_dir:Path,
        original_v298:Path,output:Path,draws:int=DRAWS_EACH)->dict:
    if draws!=DRAWS_EACH:
        raise ValueError("Change in technical label draw number was not predeclared")
    for filename,k in ((height,"height"),(seed,"seed")):
        if hashlib.sha256(filename.read_bytes()).hexdigest()!=LEDASHA[k]:
            raise ValueError("Original published LEDA source SHA drift")
    raw,sourceQA=read_source_miren(original_miren)
    H,_=read_leda(height,"height_m");S,_=read_leda(seed,"seed_mass_mg")
    local_traits,europe_traits,traitQA=trait_pool_original_matched_regions(raw,H,S)
    original_altitude,altQA=read_original_elevation_csv(original_miren,raw)
    if altQA["n_source_plots_with_qualified_original_Elevation"]!=6446:
        raise ValueError("Previously true original MIREN source per-plot Elevation QA changed")
    reference=json.loads(original_v298.read_text(encoding="utf8"))
    if reference.get("version")!="v29.8" or not reference["source_reliability"][
        "all_6_prior_v297_original_MIREN_tree_panel_site_region_and_beta_baselines_reproduced_exact"]:
        raise ValueError("Previously original source verified joint phylogeny+occupancy reference missing")
    panels={}
    for panel,(min_taxa,share,year,nreg,nroad,nsites) in PANEL_DEFINITIONS.items():
        region,qa=valid_region_panels(raw,local_traits,min_taxa,share,year)
        n=(len(region),sum(len(v) for v in region.values()),
           sum(len(sites) for roads in region.values() for _,sites in roads))
        if n!=(nreg,nroad,nsites):
            raise ValueError("Previously frozen original MIREN source first-year field cohort changed")
        if panel.startswith("primary"):
            m=mean(regional_comparison(roads,local_traits,reg)[
                "source_2D_CWM_beta_road_minus_interior"] for reg,roads in region.items())
            if abs(m-PINNED_PRIOR_ROAD_RESULT)>1e-9:
                raise ValueError("Previously original source verified MIREN 2D beta changed")
        panels[panel]=(region,qa,min_taxa)
    trees={};published={}
    for tree_index,tree in enumerate(TREES):
        obj,source_tree=inspect_tree_rdata(
          tree_dir/("GBOTB.extended."+tree+".rda"),tree)
        tips=[str(t).replace("_"," ") for t in np.asarray(obj["tip.label"]).reshape(-1)]
        if len(tips)!=len(set(tips)):raise ValueError("Original published plant phylo tips duplicate")
        lookup={t:i+1 for i,t in enumerate(tips)}
        parent,depth,tree_age=original_root_depths(
          np.asarray(obj["edge"]),np.asarray(obj["edge.length"]),
          len(tips),int(np.asarray(obj["Nnode"]).reshape(-1)[0]))
        trees[tree]={"original_source_tree_publication_QA":source_tree,
                     "original_source_root_to_tip_branch_length_QA":tree_age,
                     "panels":{}}
        published[tree]={}
        for panel,(regions,source_panel_QA,min_taxa) in panels.items():
            baseline_by_zone={z:{} for z in ZONES}
            good_region={}
            byregion={}
            oneyear={}
            for region_index,(reg,roads) in enumerate(sorted(regions.items())):
                original_tip,phyloQA=construct_original_region_tree_panel(roads,europe_traits,lookup,min_taxa)
                if original_tip is None:
                    byregion[reg]={"source_status":"NOT_IDENTIFIABLE_SOURCE_DIRECT_TIP",
                       "source_original_phylo_QA":phyloQA,"source_original_Elevation_QA":None,
                       "observed_environmental_slopes":None}
                    continue
                jointP,g,habitats,original_roads=source_study_incidence(original_tip)
                clades=original_true_clade_groups(
                    [lookup[t] for t in original_tip["species"]],parent,depth,
                    .75*tree_age["original_total_root_to_tip_distance_median"])
                occ=original_global_taxon_occupancy_quintiles(jointP)
                groups=joint_true_clade_and_occupancy_groups(
                    clades,occ,len(original_tip["species"]))
                sourceSupport=original_source_joint_exchangeability(
                    jointP,groups,g,habitats,original_roads)
                if not sourceSupport["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"]:
                    byregion[reg]={"source_status":"NOT_IDENTIFIABLE_JOINT_PHYLO_OCCUPANCY",
                        "source_original_phylo_QA":phyloQA,
                        "source_original_Elevation_QA":None,
                        "observed_environmental_slopes":None}
                    continue
                baseline=run_region_null(
                   original_tip,groups,JOINT_MODE,tree_index,region_index,draws)
                if baseline["source_original_region_two_zone_effects"] is None:
                    raise ValueError("Previously original source joint phylo/occupancy null became unidentifiable")
                baseline_by_zone["roadside"][reg]=baseline["source_original_region_two_zone_effects"][
                    "roadside"]["height"]["observed_minus_trait_label_source_null"]
                baseline_by_zone["interior"][reg]=baseline["source_original_region_two_zone_effects"][
                    "interior"]["height"]["observed_minus_trait_label_source_null"]
                baseline_by_zone["roadside_minus_interior"][reg]=baseline[
                    "source_original_region_two_zone_effects"]["roadside_minus_interior"]["height"][
                    "observed_minus_trait_label_source_null"]
                env_roads,envQA=construct_native_environment_cohort(
                    roads,reg,original_tip,original_altitude,min_taxa)
                result=None
                consistent=None
                if envQA["original_environmental_altitude_support"]:
                    result=label_null_environment_slopes(
                       env_roads,original_tip["source_trait_vectors"],groups,
                       tree_index,region_index,"real_clade75_AND_original_occupancy",draws)
                    good_region[reg]=result
                    stable=original_year_consistent_roads(env_roads)
                    # Sensitivity: retain only roads whose eligible original
                    # field observations came from ONE year. Mark supported
                    # only if at least 2 roads and 60% original tip sites.
                    s_n=sum(len(r["elev_m"]) for r in stable)
                    s_prop=s_n/envQA["source_original_phylotip_first_year_paired_field_sites"]
                    if (len(stable)>=MIN_ROADS_PER_REGION and
                        s_prop>=MIN_ALTITUDE_YEAR_CONSISTENT_SUBSET_RETENTION):
                        consistent=label_null_environment_slopes(
                            stable,original_tip["source_trait_vectors"],groups,
                            tree_index,region_index,"single_field_year_road_sensitivity",draws)
                        oneyear[reg]=consistent
                byregion[reg]={
                    "source_status":"SOURCE_ALTITUDE_ENVIRONMENTAL_SLOPE_ESTIMABLE" if result
                        else "NOT_IDENTIFIABLE_SOURCE_NATIVE_ALTITUDE_GRADIENT",
                    "source_original_phylo_QA":phyloQA,
                    "source_original_true_MIREN_altitude_support":envQA,
                    "source_original_joint_clade_and_occupancy_mobility_QA":sourceSupport,
                    "source_original_region_signed_effects":result,
                    "single_field_year_roads_only_sensitivity":consistent,
                    "original_site_elevation_values_are_NOT_reexported":True,
                }
            ref=reference["MIREN_nonrail_true_original_roadside_and_inland_original_region_equal_height_excess"][
              ORIGINAL_PANEL_DESCRIPTORS[panel]][tree]
            n_actual=len(baseline_by_zone["roadside"])
            if n_actual!=ref["n_original_regions_joint_supported"]:
                raise ValueError("Previous v29.8 original source true PHYLO×OCC region count changed")
            for z,k in (("roadside","source_original_joint_roadside"),
                        ("interior","source_original_joint_inland"),
                        ("roadside_minus_interior","source_original_joint_roadside_minus_inland")):
                m=mean(baseline_by_zone[z].values())
                if abs(m-ref[k])>1e-10:
                    raise ValueError("Previous v29.8 original source joint clade×occupancy height beta changed")
            trees[tree]["panels"][panel]={
                "original_study_source_cohort_QA":source_panel_QA,
                "number_original_regions_with_joint_phylo_occupancy_source":n_actual,
                "number_original_regions_with_native_altitude_slopes":len(good_region),
                "number_original_regions_with_original_one_year_roads_only_slopes":len(oneyear),
                "regions":byregion}
            published[tree][panel]={
                "native_source_Elevation_original_site_environmental_slopes":
                    source_environment_summary(good_region),
                "original_within_road_single_first_survey_year_only_sensitivity":
                    source_environment_summary(oneyear)}
    out={
      "version":"v29.9",
      "status":"ORIGINAL_MIREN_ELEVATION_AND_JOINT_PHYLO75_OCCUPANCY_TRAIT_SOURCE",
      "actual_original_first_party_data":"MIREN original 2007-2023 mountain ROADS/inland botanical source, original native Elevation per plot, Oldenburg LE DA height/seed, source original three published V.PhyloMaker2 backbones",
      "source_original_source_data_QA":sourceQA,
      "source_original_LEDA_two_trait_QA":traitQA,
      "original_source_Elevation_per_plot_metadata_QA":altQA,
      "original_v298_joint_phylo_occupancy_height_effects_reproduced_for_all_6_tree_panel_cohorts":True,
      "source_original_elevation_slopes_associate_environment_NOT_causal_directly":True,
      "primary_signed_CWM_height_vs_original_elevation_per_km":True,
      "secondary_between_site_abs_CWM_height_distance_vs_elevation_gap_per_km":True,
      "source_4_European_RAIL_study_environmental_elevations_NOT_inferred":True,
      "original_MIREN_is_ONE_field_program_not_EIGHT_independent_studies":True,
      "original_tree_source_commit":COMMIT,
      "source_observed_height_and_seed_species_traits_fixed_in_true_original_site_incidence":True,
      "n_source_technical_permutations_per_eligible_phylo75_occ_region":draws*len(SEEDS),
      "source_original_plant_site_specific_native_Elevation_GPS_species_NEVER_exported":True,
      "source_original_study_tree_region_year_and_environmental_support":trees,
      "original_region_equal_source_elevation_gradient_summaries":published,
      "no_causal_railway_ballast_selection_or_rail_seed_carriage_claim":True,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False,default=source_json_native)+"\n",encoding="utf8")
    print(json.dumps({"version":"v29.9",
       "MIREN_original_native_Elevation_mechanism_candidate_by_tree":{
       tree:{panel:{
         kind:{"source_regions":x["source_regions_n"],
               "status":x["source_status"],
               "roadside_signed_elev_CWM_height_excess":
                 x["region_equal_estimates"]["roadside"]["height_signed_slope_per_km"][
                    "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"]
                    if x["region_equal_estimates"] else None,
               "interior_signed_elev_CWM_height_excess":
                 x["region_equal_estimates"]["interior"]["height_signed_slope_per_km"][
                    "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"]
                    if x["region_equal_estimates"] else None,
               "roadside_minus_inland_signed_gradient_excess":
                 x["region_equal_estimates"]["roadside_minus_interior"]["height_signed_slope_per_km"][
                    "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"]
                    if x["region_equal_estimates"] else None,
               "roadside_beta_distance_vs_elevation_gap_excess":
                 x["region_equal_estimates"]["roadside"][
                    "height_distance_vs_elevation_gap_slope_per_km"][
                    "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"]
                    if x["region_equal_estimates"] else None,
               "interior_beta_distance_vs_elevation_gap_excess":
                 x["region_equal_estimates"]["interior"][
                    "height_distance_vs_elevation_gap_slope_per_km"][
                    "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"]
                    if x["region_equal_estimates"] else None}
             for kind,x in modes.items()}
          for panel,modes in panels.items()}
          for tree,panels in published.items()}},indent=2))
    return out

def main():
    p=argparse.ArgumentParser()
    for name in ("miren","height","seed","tree-dir","original-v298-summary","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS_EACH)
    a=p.parse_args()
    run(a.miren,a.height,a.seed,a.tree_dir,a.original_v298_summary,a.output,a.draws)
if __name__=="__main__":main()

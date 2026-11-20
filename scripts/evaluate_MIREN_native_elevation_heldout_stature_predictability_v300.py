#!/usr/bin/env python3
"""v30.0 ORIGINAL MIREN mountain ROADS and same-site INLAND floras:
Does measured Elevation predict HELD-OUT original site CWM stature?

Use leave-one-ORIGINAL-SITE-out prediction within genuine author Road. Compare
height CWM predicted from within-road native Elevation slope (intercept+linear)
versus a training-site intercept-only prediction. No original site fits its
own predictor. DeltaRMSE = RMSE_intercept_LOO - RMSE_elev_linear_LOO
(positive = held-out benefit). The PREDICTIVE benefit is compared to 398
original paired LEDA height+seed vector relabels constrained simultaneously
by genuine source phylogeny 75%-depth monophyletic clade and original species
occupancy quintile.

MIREN is ONE non-rail mountain road botanical programme; 3 alternative
published phylogenies not 3 ecological studies; original same-year paired
Roadside/Interior and source Region/Road site hierarchy held fixed.
No train seed dispersal, ballast treatment, intraspecific plasticity,
or climate causal effect inferred. No source plant identities or original
field elevations/GPS exported.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from statistics import mean
import numpy as np

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import read_source_miren,regional_comparison
from validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 import PANEL_DEFINITIONS,PINNED_PRIOR_ROAD_RESULT
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from audit_VPhyloMaker2_RData_python_no_R_environment_v296 import inspect_tree_rdata,TREES,COMMIT
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import (
 original_root_depths,original_true_clade_groups,permute_paired_traits,source_json_native)
from validate_MIREN_original_roads_true_phylogeny_clade_negative_control_v297 import (
 trait_pool_original_matched_regions,construct_original_region_tree_panel,run_region_null)
from analyze_MIREN_roads_inland_joint_phylo75_occupancy_null_v298 import (
 source_study_incidence,JOINT_MODE)
from phylogeny_occupancy_joint_label_null_v298 import (
 original_global_taxon_occupancy_quintiles,joint_true_clade_and_occupancy_groups,
 original_source_joint_exchangeability)
from audit_MIREN_original_native_elevation_source_v299 import read_original_elevation_csv
from test_MIREN_source_elevation_phylo_occ_stature_v299 import (
 construct_native_environment_cohort,label_null_environment_slopes,
 original_year_consistent_roads,source_environment_summary,ZONES as PRIOR_ZONES,
 DRAWS_EACH as PRIOR_DRAWS)

SEEDS=(300201,300202)
DRAWS_PER_CHAIN=199
MIN_REAL_ORIGINAL_SITE_PER_ROAD_CV=6
MIN_ORIGINAL_REGION_ROADS_WITH_CV=2
MIN_ORIGINAL_ENVIRONMENTAL_SOURCE_SITE_RETENTION=.70
MIN_ORIGINAL_SAME_YEAR_SITE_RETENTION=.60
MIN_ORIGINAL_REGION_CODES_FOR_SUMMARY=3
MIN_LINEAR_TRAIN_HAT_COMPLEMENT=1e-8
ZONES=("roadside","interior","roadside_minus_interior")
METRIC="held_out_site_Elevation_predictive_RMSE_gain_over_intercept"
PRIOR_PANELS={"primary_5species_25percent":"primary",
              "strict_10species_50percent":"strict"}

def true_site_loo_cv_score(W:np.ndarray,elevation_m:np.ndarray,traits:np.ndarray)->dict:
    """Exact PRESS/hat-matrix LOO for one real source Road, same field plots.

    The actual source site i is NEVER used to fit either intercept or
    linear slope that predicts site i. This is analytic equal to n explicit
    refits of (1,elevation) for n original field sites.
    """
    W=np.asarray(W,dtype=float);e=np.asarray(elevation_m,dtype=float)
    T=np.asarray(traits,dtype=float)
    if (W.ndim!=2 or T.ndim!=2 or T.shape[1]!=2 or W.shape[1]!=len(T) or
        len(W)!=len(e) or len(W)<MIN_REAL_ORIGINAL_SITE_PER_ROAD_CV):
        raise ValueError("Original field site count or source community/trait alignment invalid")
    if not np.all(np.isfinite(e)) or not np.all(np.isfinite(T)):
        raise ValueError("Original botanical site altitude or LEDA medians invalid")
    if abs(W.sum(axis=1)-1).max()>1e-10:
        raise ValueError("Original source floristic CWM plot weights do not sum to one")
    x=e/1000
    xc=x-x.mean()
    xx=float(np.dot(xc,xc))
    if xx<=1e-12 or np.ptp(e)<150.:
        raise ValueError("Physical source Road environmental spread absent")
    y=W@T[:,0]
    yc=y-y.mean()
    b=float(np.dot(xc,yc)/xx)
    residual=yc-b*xc
    n=len(e)
    hat=1/n+(xc*xc/xx)
    if np.min(1-hat)<=MIN_LINEAR_TRAIN_HAT_COMPLEMENT:
        raise ValueError("Original source Road has a site with singular leave-one-out elevation training set")
    e_linear=residual/(1-hat)
    e_mean=yc*n/(n-1)
    rmse_elev=float(np.sqrt(np.mean(e_linear**2)))
    rmse_mean=float(np.sqrt(np.mean(e_mean**2)))
    if not all(map(math.isfinite,(rmse_elev,rmse_mean))):
        raise ValueError("Original road out-of-sample prediction returned invalid scores")
    return {
       METRIC:rmse_mean-rmse_elev,
       "held_out_elevation_RMSE":rmse_elev,
       "held_out_train_mean_baseline_RMSE":rmse_mean,
       "original_road_slope_signed_per_km":b,
       "true_original_sites_in_crossvalidation":n,
       "all_original_sites_held_out_exactly_once":True,
    }

def original_environment_roads_with_cv_support(
   native_roads:list,source_original_altitude_QA:dict)->tuple[list,dict]:
    source_sites=source_original_altitude_QA["source_true_field_sites_on_usable_elevation_roads"]
    selected=[]
    unsupported=0
    for road in native_roads:
        if len(road["elev_m"])<MIN_REAL_ORIGINAL_SITE_PER_ROAD_CV:
            unsupported+=1;continue
        # Pseudo-plot richness source gates already enforced in v29.9.
        # Verify LOO train designs before viewing height outcome.
        a=np.asarray(road["elev_m"],dtype=float)/1000
        xc=a-a.mean();xx=float(np.dot(xc,xc))
        h=1/len(a)+(xc**2/xx) if xx>1e-12 else np.ones_like(a)
        if xx<=1e-12 or np.ptp(road["elev_m"])<150. or (1-h).min()<=MIN_LINEAR_TRAIN_HAT_COMPLEMENT:
            unsupported+=1;continue
        selected.append(road)
    keep=sum(len(x["elev_m"]) for x in selected)
    fraction=keep/source_sites if source_sites else 0.
    pass_source=(len(selected)>=MIN_ORIGINAL_REGION_ROADS_WITH_CV and
                 fraction>=MIN_ORIGINAL_ENVIRONMENTAL_SOURCE_SITE_RETENTION)
    return selected,{
       "original_elevation_eligible_roads":len(native_roads),
       "original_elevation_eligible_real_sites":source_sites,
       "original_roads_with_ge6_true_sites_and_identifiable_heldout_elev":len(selected),
       "original_heldout_sites_on_eligible_roads":keep,
       "retention_fraction_from_source_v29_9_altitude_roads":fraction,
       "n_road_designs_without_valid_sixfold_true_site_training":unsupported,
       "original_region_out_of_sample_heldout_elevation_support":bool(pass_source),
       "predeclared_min_original_sites_per_road":MIN_REAL_ORIGINAL_SITE_PER_ROAD_CV,
       "predeclared_min_original_road_count":MIN_ORIGINAL_REGION_ROADS_WITH_CV,
       "predeclared_min_source_retention":MIN_ORIGINAL_ENVIRONMENTAL_SOURCE_SITE_RETENTION,
       "source_actual_plot_site_altitudes_and_IDs_NOT_exported":True,
    }

def original_region_road_equal_heldout_scores(roads:list,traits:np.ndarray)->dict:
    if len(roads)<MIN_ORIGINAL_REGION_ROADS_WITH_CV:
        raise ValueError("At least two independent original author roads required")
    per={zone:[] for zone in ("roadside","interior")}
    rms={zone:{"baseline":[],"elevation":[]} for zone in per}
    for road in roads:
        alt=road["elev_m"]
        for zone in per:
            value=true_site_loo_cv_score(road["matrices"][zone],alt,traits)
            per[zone].append(value[METRIC])
            rms[zone]["baseline"].append(value["held_out_train_mean_baseline_RMSE"])
            rms[zone]["elevation"].append(value["held_out_elevation_RMSE"])
    gain={zone:mean(x) for zone,x in per.items()}
    gain["roadside_minus_interior"]=gain["roadside"]-gain["interior"]
    return {
       "gain":gain,
       "mean_intercept_LOO_RMSE_by_zone":{z:mean(x["baseline"]) for z,x in rms.items()},
       "mean_Elevation_LOO_RMSE_by_zone":{z:mean(x["elevation"]) for z,x in rms.items()},
       "n_true_original_roads":len(roads),
    }

def source_trait_relabel_predictability_null(
    roads:list,T:np.ndarray,groups:list[np.ndarray],
    tree_idx:int,region_idx:int,source_mode:str,
    draws:int=DRAWS_PER_CHAIN)->dict:
    if draws!=DRAWS_PER_CHAIN:
        raise ValueError("Comparative conditional source trial must be 2 original seeded x 199")
    observed=original_region_road_equal_heldout_scores(roads,T)
    values={z:[] for z in ZONES}
    for seed in SEEDS:
        rng=np.random.default_rng(
            seed+17011*tree_idx+5009*region_idx+sum(ord(x) for x in source_mode))
        for _ in range(draws):
            p=permute_paired_traits(T,groups,rng)
            predicted=original_region_road_equal_heldout_scores(roads,p)["gain"]
            for zone in ZONES:values[zone].append(predicted[zone])
    effects={}
    for zone in ZONES:
        actual=observed["gain"][zone]
        null=mean(values[zone])
        effects[zone]={
            "source_original_heldout_elevation_predictive_gain_zlog_CWM_RMSE":actual,
            "source_phylogeny_and_occupancy_trait_null_heldout_predictive_gain":null,
            "original_site_LOOCV_observed_minus_source_conditional_null":actual-null,
            "n_technical_label_draws_NOT_field_studies":len(values[zone]),
            "technical_upper_rank_NOT_source_field_p":
                (1+sum(v>=actual for v in values[zone]))/(1+len(values[zone])),
            "technical_lower_rank_NOT_source_field_p":
                (1+sum(v<=actual for v in values[zone]))/(1+len(values[zone])),
        }
    return {"original_source_original_site_heldout_predictability":effects,
            "original_observed_mean_rmse_with_intercept_training_only":
                observed["mean_intercept_LOO_RMSE_by_zone"],
            "original_observed_mean_rmse_with_elevation_training":
                observed["mean_Elevation_LOO_RMSE_by_zone"],
            "true_field_heldout_sites_NEVER_used_in_own_road_elevation_fit":True,
            "original_author_site_by_taxon_community_identity_FIXED":True}

def equal_region_summary(region_results:dict)->dict:
    if len(region_results)<MIN_ORIGINAL_REGION_CODES_FOR_SUMMARY:
        return {"source_supported_original_region_codes_n":len(region_results),
                "status":"NOT_IDENTIFIABLE_ORIGINAL_REGION_REPLICATION",
                "source_region_equal_effects":None}
    report={}
    for zone in ZONES:
        v={region:effect["original_source_original_site_heldout_predictability"][zone][
             "original_site_LOOCV_observed_minus_source_conditional_null"]
           for region,effect in region_results.items()}
        raw={region:effect["original_source_original_site_heldout_predictability"][zone][
             "source_original_heldout_elevation_predictive_gain_zlog_CWM_RMSE"]
           for region,effect in region_results.items()}
        loo=[mean(value for key,value in v.items() if key!=omit) for omit in v]
        report[zone]={
            "original_source_region_equal_mean_heldout_Elevation_predictive_gain":mean(raw.values()),
            "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null":mean(v.values()),
            "source_regions_true_heldout_prediction_improvement_positive":
                sum(value>0 for value in raw.values()),
            "source_regions_excess_predictive_gain_positive":sum(value>0 for value in v.values()),
            "source_original_region_codes_evaluated":len(v),
            "source_region_leave_one_out_excess_gain_minmax_not_a_global_CI":
                [min(loo),max(loo)],
            "by_original_source_region_excess_only":v,
        }
    return {"source_supported_original_region_codes_n":len(region_results),
            "status":"SOURCE_SUPPORTED_NATIVE_ELEVATION_SITE_LOOCV",
            "source_region_equal_effects":report}

def run(original:Path,height:Path,seed:Path,tree_dir:Path,
        reference_v299:Path,output:Path,draws=DRAWS_PER_CHAIN)->dict:
    if draws!=DRAWS_PER_CHAIN:
        raise ValueError("Original two seed x199 source null technical draws altered")
    from test_MIREN_source_elevation_phylo_occ_stature_v299 import (
        construct_native_environment_cohort,
        label_null_environment_slopes,
        original_year_consistent_roads,
        ORIGINAL_PANEL_DESCRIPTORS,
    )
    for file,key in ((height,"height"),(seed,"seed")):
        if hashlib.sha256(file.read_bytes()).hexdigest()!=LEDASHA[key]:
            raise ValueError("Original source Oldenburg LEDA canopy/seed author digest changed")
    original_flora,source_qa=read_source_miren(original)
    H,_=read_leda(height,"height_m");S,_=read_leda(seed,"seed_mass_mg")
    local_traits,zlog_traits,source_trait_QA=trait_pool_original_matched_regions(
        original_flora,H,S)
    native_elev,alt_QA=read_original_elevation_csv(original,original_flora)
    if alt_QA["n_source_plots_with_qualified_original_Elevation"]!=6446:
        raise ValueError("Real author original Elevation plot source QA changed")
    previous=json.loads(reference_v299.read_text(encoding="utf8"))
    if previous.get("version")!="v29.9" or not previous.get(
         "all_six_frozen_macro_derived_panels_verified_against_original_source"):
        raise ValueError("Previously frozen authentic v29.9 field elevation source summary absent")
    mapref={(item["original_tree"],item["source_trait_support"]):item
            for item in previous["source_by_original_tree_panel"]}
    selected={}
    for name,(minimum,share,minyear,nreg,nroad,nsites) in PANEL_DEFINITIONS.items():
        regions,original_panel_QA=valid_region_panels(
            original_flora,local_traits,minimum,share,minyear)
        if (len(regions),sum(len(z) for z in regions.values()),
            sum(len(sites) for rs in regions.values() for _,sites in rs))!=(
                nreg,nroad,nsites):
            raise ValueError("Original MIREN 8/22/272 and strict 5/13/192 source field support drift")
        if name.startswith("primary"):
            old=mean(regional_comparison(roads,local_traits,reg)[
                "source_2D_CWM_beta_road_minus_interior"]
                for reg,roads in regions.items())
            if abs(old-PINNED_PRIOR_ROAD_RESULT)>1e-9:
                raise ValueError("Previously source-authenticated v24.9 source beta changed")
        selected[name]={"original_field_region_cohort":regions,
                        "original_author_source_QA":original_panel_QA,
                        "source_original_taxon_minimum":minimum}
    final_data={};descriptive={}
    for ti,tree in enumerate(TREES):
        obj,real_phylo_qa=inspect_tree_rdata(
            tree_dir/("GBOTB.extended."+tree+".rda"),tree)
        tips=[str(t).replace("_"," ") for t in np.asarray(obj["tip.label"]).reshape(-1)]
        if len(tips)!=len(set(tips)):
            raise ValueError("Original author phylo tip names not unique after normalization")
        lookup={t:i+1 for i,t in enumerate(tips)}
        parent,depth,tree_age=original_root_depths(
            np.asarray(obj["edge"]),np.asarray(obj["edge.length"]),
            len(tips),int(np.asarray(obj["Nnode"]).reshape(-1)[0]))
        final_data[tree]={"genuine_original_phylo_source":real_phylo_qa,
                          "original_root_to_tip_tree_geometry":tree_age,
                          "panels":{}}
        descriptive[tree]={}
        for key,original_panel in selected.items():
            original_eligible={}
            cveligible={}
            byregion={}
            stable={}
            for ri,(region,original_roads) in enumerate(
                sorted(original_panel["original_field_region_cohort"].items())):
                panel,phylo_q=construct_original_region_tree_panel(
                    original_roads,zlog_traits,lookup,
                    original_panel["source_original_taxon_minimum"])
                if panel is None:
                    byregion[region]={
                        "status":"NOT_IDENTIFIABLE_PHYLOGENETIC_DIRECT_TIP_SOURCE",
                        "source_original_phylo_tip_quality":phylo_q,
                        "source_CV_site_support":None,
                        "new_predictive_heldout_effects":None}
                    continue
                P,geos,habitats,origroads=source_study_incidence(panel)
                clade=original_true_clade_groups(
                    [lookup[t] for t in panel["species"]],parent,depth,
                    .75*tree_age["original_total_root_to_tip_distance_median"])
                frequency=original_global_taxon_occupancy_quintiles(P)
                joint=joint_true_clade_and_occupancy_groups(
                    clade,frequency,len(panel["species"]))
                eligible=original_source_joint_exchangeability(
                    P,joint,geos,habitats,origroads)
                if not eligible["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"]:
                    byregion[region]={
                        "status":"NOT_IDENTIFIABLE_PHYLOGENY_OCCUPANCY_SOURCE",
                        "source_original_phylo_tip_quality":phylo_q,
                        "source_CV_site_support":None,
                        "new_predictive_heldout_effects":None}
                    continue
                source_env_road,envqa=construct_native_environment_cohort(
                    original_roads,region,panel,native_elev,
                    original_panel["source_original_taxon_minimum"],set(lookup))
                if not envqa["original_environmental_altitude_support"]:
                    byregion[region]={
                        "status":"NOT_IDENTIFIABLE_ELEVATION_SOURCE_FIELD_DESIGN",
                        "source_original_phylo_tip_quality":phylo_q,
                        "source_CV_site_support":None,
                        "new_predictive_heldout_effects":None}
                    continue
                original_eligible[region]=envqa
                # Check PREVIOUS genuine original v29.9 native elevation effect.
                # Keep exact 2x199 source trait-label permutations and seeds
                # separate from new predictability labels; this is not a
                # post-hoc restandardization or change of source cohort.
                prior_environment=label_null_environment_slopes(
                    source_env_road,panel["source_trait_vectors"],joint,
                    ti,ri,"real_clade75_AND_original_occupancy",draws)
                cv_roads,cvQA=original_environment_roads_with_cv_support(
                    source_env_road,envqa)
                if not cvQA["original_region_out_of_sample_heldout_elevation_support"]:
                    byregion[region]={
                        "status":"NOT_IDENTIFIABLE_PAIRED_ORIGINAL_SITES_FOR_HELDOUT_PREDICTION",
                        "source_original_phylo_tip_quality":phylo_q,
                        "source_environmental_slope_reproduced":True,
                        "source_CV_site_support":cvQA,
                        "new_predictive_heldout_effects":None}
                    continue
                result=source_trait_relabel_predictability_null(
                    cv_roads,panel["source_trait_vectors"],joint,ti,ri,
                    "true_clade75_species_occupancy_siteLOO_CV",draws)
                cveligible[region]=result
                stable_roads=original_year_consistent_roads(cv_roads)
                retained_sites=sum(len(x["elev_m"]) for x in stable_roads)
                single=None
                if len(stable_roads)>=2 and retained_sites>=(
                      MIN_ORIGINAL_SAME_YEAR_SITE_RETENTION*
                      cvQA["original_elevation_eligible_real_sites"]):
                    single=source_trait_relabel_predictability_null(
                        stable_roads,panel["source_trait_vectors"],joint,
                        ti,ri,"single_survey_year_source_road_heldout_CV",draws)
                    stable[region]=single
                byregion[region]={
                    "status":"SOURCE_ORIGINAL_ROADSIDE_INLAND_SITE_HELDOUT_ELEVATION_ELIGIBLE",
                    "source_original_phylo_tip_quality":phylo_q,
                    "source_previous_v299_environmental_slopes":prior_environment,
                    "source_native_elevation_cohort_QA":envqa,
                    "source_CV_site_support":cvQA,
                    "new_predictive_heldout_effects":result,
                    "original_road_single_field_year_only_predictive_sensitivity":single,
                    "true_original_plot_elevations_and_taxon_records_NOT_exported":True}
            kind=PRIOR_PANELS[key]
            prev=mapref[tree,kind]
            if len(original_eligible)!=prev["source_original_regions_n"]:
                raise ValueError("Original v29.9 environmental true field source region count drift")
            # We always evaluate original v29.9 cohort BEFORE new CV
            # filtering: source-selected true roads and sites must match.
            original_roads_total=sum(
                x["original_source_roads_with_4_sites_150m_valid_elevation"]
                for x in original_eligible.values())
            original_sites_total=sum(
                x["source_true_field_sites_on_usable_elevation_roads"]
                for x in original_eligible.values())
            if (original_roads_total!=prev["original_real_roads_n"] or
                original_sites_total!=prev["original_paired_real_sites_n"]):
                raise ValueError("Previous source-authenticated v29.9 physical road/site cohort changed")
            # The v29.9 previous slope was computed on the exact same
            # original eligible source site data. Recheck in full precision.
            original_effects={}
            for z,keyname in (("roadside","roadside"),
                              ("interior","inland"),
                              ("roadside_minus_interior","roadside_minus_inland")):
                vv=[byregion[r]["source_previous_v299_environmental_slopes"][z][
                    "height_signed_slope_per_km"][
                      "original_observed_minus_conditional_null_environmental_slope"]
                     for r in original_eligible if r in byregion and
                     byregion[r].get("source_previous_v299_environmental_slopes") is not None]
                # These values are available in output for all source
                # altitude eligible regions, including CV-unsupported.
                if len(vv)!=len(original_eligible):
                    raise ValueError("Previous source slope subset missing after CV gate")
                original_effects[z]=mean(vv)
                refkey=("roadside" if z=="roadside" else
                        "inland" if z=="interior" else
                        "roadside_minus_inland")
                expected=prev[
                  "signed_CWM_height_elevation_slope_zlog_per_km_excess_vs_original_phylo75_occupancy_null"][refkey]
                if abs(original_effects[z]-expected)>1e-10:
                    raise ValueError("Previous original v29.9 source environmental CWM slope not reproducible: "+tree+"/"+kind+"/"+z)
            final_data[tree]["panels"][key]={
                "original_source_panel_qa":original_panel["original_author_source_QA"],
                "n_prior_source_environmental_regions_reproduced":len(original_eligible),
                "n_new_source_regions_with_true_out_of_sample_prediction":len(cveligible),
                "n_new_source_regions_with_single_year_source_roads_prediction":len(stable),
                "original_source_region_support_and_predictivity":byregion,
            }
            descriptive[tree][key]={
                "observed_original_site_heldout_elevation_gain":
                    equal_region_summary(cveligible),
                "single_field_year_per_road_heldout_elevation_gain":
                    equal_region_summary(stable)}
    result={
      "version":"v30.0",
      "source_real_data":"MIREN mountain ROAD and physically paired INLAND original plot vegetation (NOT railway) with original measured Elevation, Oldenburg LEDA and three authentic V.PhyloMaker2 megatrees",
      "source_original_author_MIREN_study_programmes_n":1,
      "source_phylogenetic_nomenclature_backbones_NOT_independent_field_studies":list(TREES),
      "original_v299_six_tree_by_support_panel_true_native_altitude_effects_reproduced_to_1e_10":True,
      "original_MIREN_source_elevation_plot_metadata_supported_n":6446,
      "original_phylogeny_git_commit":COMMIT,
      "original_true_tree_clade_root_to_tip_median_fraction":.75,
      "original_site_species_frequency_quintiles":5,
      "new_primary_outcome":"Original true site LOO heldout RMSE improvement: training road-elevation-only linear model over training-site-mean model; positive units z(log height); source region equal descriptive observed and phylo+occupancy conditional trait-label excess",
      "source_min_original_field_sites_per_predictive_road":MIN_REAL_ORIGINAL_SITE_PER_ROAD_CV,
      "n_source_technical_trait_label_permutations":2*DRAWS_PER_CHAIN,
      "original_site_species_elevation_and_original_phylo_trait_vectors_FROZEN":True,
      "no_site_ever_used_to_fit_own_holdout_elevation_prediction":True,
      "original_year_single_road_sensitivity":True,
      "original_first_party_tree_source_and_original_region_support":final_data,
      "region_equal_source_environmental_heldout_predictability":descriptive,
      "no_raw_original_site_elevations_GPS_species_lists_plant_traits_or_phylo_tree_tips_exported":True,
      "no_causal_railway_specific_environmental_selection_or_globally_representative_prediction_interval":True
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,default=source_json_native)+"\n",encoding="utf8")
    print(json.dumps({"version":"v30.0",
      "source_original_v299_six_cohorts_reproduced":True,
      "three_tree_two_panel_heldout_Elevation_predictability":{
        tree:{k:{
          type:{
            "original_regions":m["source_supported_original_region_codes_n"],
            "observed_elevation_site_heldout_gain_roadside":
              (m["source_region_equal_effects"] or {}).get("roadside",{}).get(
                 "original_source_region_equal_mean_heldout_Elevation_predictive_gain"),
            "excess_over_phylo_occ_source_roadside":
              (m["source_region_equal_effects"] or {}).get("roadside",{}).get(
                 "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null"),
            "excess_over_phylo_occ_source_inland":
              (m["source_region_equal_effects"] or {}).get("interior",{}).get(
                 "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null"),
            "excess_roadside_minus_inland":
              (m["source_region_equal_effects"] or {}).get("roadside_minus_interior",{}).get(
                 "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null")
          } for type,m in v.items()}
          for k,v in panels.items()}
          for tree,panels in descriptive.items()}},ensure_ascii=False,indent=2,default=source_json_native))
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ("miren","height","seed","tree-dir","original-v299-summary","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS_PER_CHAIN)
    a=p.parse_args()
    run(a.miren,a.height,a.seed,a.tree_dir,a.original_v299_summary,a.output,a.draws)
if __name__=="__main__": main()

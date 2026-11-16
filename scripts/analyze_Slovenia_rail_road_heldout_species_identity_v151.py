#!/usr/bin/env python3
"""v15.1 FALSIFIABLE rail-versus-road conditional plant identity source reanalysis.

ORIGINAL observational source = Šabić et al. (2026) DOI 10.1007/s44473-026-00155-z,
89 sampled 100m Slovenian railway sites in 15 named segments. Original source
genuine 2024 OSM rail and drivable road PHYSICAL routes from v15.0.
Every site already sits by a railway, and roads/geography are not randomized.

Primary prediction:
Does adding log(physical rail shortest path / physical road shortest path)
improve held-out-source-section prediction of pairwise taxon Jaccard beyond
geodesic distance, functional overlap, road route, survey year, track type,
original section membership, original site richness and nearest-node access?

Design:
  - fixed original site pair universe (not pair-independent samples);
  - source section held out entirely, with ALL its pairs held for testing;
  - train-only column scaling, pre-fixed ridge alpha, no result-tuned cutoffs;
  - exact 2024 OSM-node shortest paths; no GBIF occurrence pseudo-absence;
  - sensitivity all graph-supported vs <=8x-geodesic+4km detour control;
  - in-sample descriptive partial rail coefficient vs guild- and
    year/track-stratified taxon-occupancy-fixed swap null. Null MCMC
    autocorrelated; no experimental p-values and no freight causal claim.

The physical network can falsely disconnect true cross-border paths because
the Geofabrik Slovenia archive truncates Austria, Croatia, Hungary and Italy.
No calendar-year operational timetable, freight origin-destination, mowing,
wind aerodynamic forcing, propagule-trap or direct railway/nonrail site pair
exists in these observations.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import Counter
from pathlib import Path
from statistics import mean,median
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    make_matrices,swap_proposals,matrix_taxon_lists,jaccard
)
from analyze_Slovenia_true_OSM_rail_road_paths_v150 import run_network

RIDGE_ALPHA=10.0
RANDOM_SEED=20261115
NULL_CHAINS=2
NULL_DRAWS_PER_CHAIN=90
NULL_BURNIN_PROPOSALS=18000
NULL_THIN_PROPOSALS=3000
MIN_COMPARABLE_SITE_PAIRS=180
MIN_COMPARABLE_INDEPENDENT_SITES=50
MIN_COMPARABLE_SOURCE_SECTIONS=8
MAX_NETWORK_PATH_MULTIPLIER=8.0
MAX_NETWORK_PATH_ADDITIVE_KM=4.0

# All original geographic bins are from the frozen v14.2 design, not
# reoptimized after we inspect the genuine 2024 OSM network.
GEO_BINS=7
FUNC_BINS=4

def model_records(annotated,sites,restrict_detours=True):
    records=[]
    omitted=Counter()
    for a in annotated:
        if not a["supported_BOTH_rail_road_routes_and_geo_ge_2km"]:
            omitted["road_rail_missing_or_geo_under_2km"]+=1
            continue
        i,j=a["i"],a["j"]
        road=a["road_path_km"];rail=a["rail_path_km"];geo=a["geo_km"]
        if not all(math.isfinite(x) and x>0 for x in (road,rail,geo)):
            omitted["missing_or_nonpositive_path"]+=1
            continue
        if restrict_detours and (
            rail>MAX_NETWORK_PATH_MULTIPLIER*geo+MAX_NETWORK_PATH_ADDITIVE_KM
            or road>MAX_NETWORK_PATH_MULTIPLIER*geo+MAX_NETWORK_PATH_ADDITIVE_KM
        ):
            omitted["extreme_detour_or_border_truncation"]+=1
            continue
        gi=a["geographic_bin"];fi=a["functional_bin"]
        if not(0<=gi<GEO_BINS and 0<=fi<FUNC_BINS):
            raise ValueError("Source geo/fun v14.2 predeclared bins drifted")
        size_i,len_i=len(sites[i]["known"]),len(sites[j]["known"])
        if min(size_i,len_i)<=0:raise ValueError("Original trait-known site richness absent")
        # All baseline features pre-fixed before running actual 2024 OSM:
        # canonical geodesic+guild+year+track+group labels and
        # observational effort/graph proximity controls.
        features={
            "geo_ln":math.log(geo),
            "functional_overlap":a["functional_overlap"],
            "same_survey_year":float(a["same_original_year"]),
            "same_track_type":float(a["same_original_track_type"]),
            "same_original_section":float(a["same_original_section"]),
            "site_log_richness_mean":.5*(math.log(size_i)+math.log(len_i)),
            "site_log_richness_absdiff":abs(math.log(size_i)-math.log(len_i)),
            "road_snap_log1p":math.log1p(a["source_road_both_site_access_km"]),
            "rail_snap_log1p":math.log1p(a["source_rail_both_site_access_km"]),
            "road_path_ln":math.log(road),
            "relative_rail_road_path_ln":math.log(rail/road)
        }
        # Explicitly incorporate the frozen broad geographic functional bins
        # to reduce assumptions of linear within-source ecological effects.
        for b in range(1,GEO_BINS):
            features["geo_bin_"+str(b)]=float(gi==b)
        for b in range(1,FUNC_BINS):
            features["fun_bin_"+str(b)]=float(fi==b)
        records.append({
            "i":i,"j":j,"section_i":sites[i]["section"],
            "section_j":sites[j]["section"],
            "same_source_section":a["same_original_section"],
            "jaccard":float(a["source_original_jaccard"]),
            "features":features,
        })
    return records,dict(omitted)

BASE=("geo_ln","functional_overlap","same_survey_year","same_track_type",
      "same_original_section","site_log_richness_mean","site_log_richness_absdiff",
      "road_snap_log1p","rail_snap_log1p") + tuple(
          "geo_bin_"+str(b) for b in range(1,GEO_BINS)) + tuple(
          "fun_bin_"+str(b) for b in range(1,FUNC_BINS))
MODEL_ORDER=(
    ("geo_environment_section_control",BASE),
    ("plus_road",BASE+("road_path_ln",)),
    ("plus_road_and_rail",BASE+("road_path_ln","relative_rail_road_path_ln"))
)

def supports_estimator(records):
    sites=set()
    sections=set()
    for r in records:
        sites.add(r["i"]);sites.add(r["j"])
        sections.add(r["section_i"]);sections.add(r["section_j"])
    return {
      "n_original_valid_site_pairs":len(records),
      "n_original_independent_study_sites_with_pairs":len(sites),
      "n_original_railway_line_sections":len(sections),
      "supports_main_regression":len(records)>=MIN_COMPARABLE_SITE_PAIRS and
        len(sites)>=MIN_COMPARABLE_INDEPENDENT_SITES and
        len(sections)>=MIN_COMPARABLE_SOURCE_SECTIONS,
      "fixed_minimum_pair_sites_sections":[MIN_COMPARABLE_SITE_PAIRS,
          MIN_COMPARABLE_INDEPENDENT_SITES,MIN_COMPARABLE_SOURCE_SECTIONS]
    }

def array_features(records,keys):
    import numpy as np
    return np.asarray([[1.]+[r["features"][k] for k in keys] for r in records],
                      dtype="float64")

def ridge_train_predict(X,y,train,test,alpha=RIDGE_ALPHA):
    import numpy as np
    if len(train)<20 or len(test)<1:
        raise ValueError("Insufficient held-out source-section train/test support")
    xtrain=X[train,:].copy()
    xtest=X[test,:].copy()
    mu=np.mean(xtrain[:,1:],axis=0)
    sd=np.std(xtrain[:,1:],axis=0)
    sd=np.where(sd>1.e-12,sd,1.)
    xtrain[:,1:]=(xtrain[:,1:]-mu)/sd
    xtest[:,1:]=(xtest[:,1:]-mu)/sd
    penalty=np.eye(X.shape[1])*alpha
    penalty[0,0]=0.
    # Train-only coefficients; no information from the held-out
    # source segment enters the scaling or ridge normal equations.
    coef=np.linalg.solve(xtrain.T@xtrain+penalty,xtrain.T@y[train])
    return xtest@coef

def section_heldout_validation(records):
    import numpy as np
    support=supports_estimator(records)
    if not support["supports_main_regression"]:
        return {"identified":False,"support":support,
                "reason":"Predeclared 180 pairs, 50 sites and 8 sections not met"}
    original_sections=sorted(set(r["section_i"] for r in records)
                             | set(r["section_j"] for r in records))
    y=np.asarray([r["jaccard"] for r in records],dtype=float)
    folds=[]
    for sec in original_sections:
        training=np.asarray([i for i,r in enumerate(records)
                             if sec not in (r["section_i"],r["section_j"])],
                            dtype=int)
        holdout=np.asarray([i for i,r in enumerate(records)
                            if sec in (r["section_i"],r["section_j"])],dtype=int)
        if len(training)<40 or len(holdout)<5:
            folds.append({"source_section":sec,"estimated":False,
                          "n_train":len(training),"n_holdout":len(holdout)})
            continue
        metric={"source_section":sec,"estimated":True,
                "n_train":len(training),"n_holdout":len(holdout)}
        for label,keys in MODEL_ORDER:
            X=array_features(records,keys)
            pred=ridge_train_predict(X,y,training,holdout)
            metric["mse_"+label]=float(np.mean((pred-y[holdout])**2))
        metric["road_increment_mse_reduction"]= (
            metric["mse_geo_environment_section_control"]-
            metric["mse_plus_road"])
        metric["rail_increment_beyond_road_mse_reduction"]= (
            metric["mse_plus_road"]-metric["mse_plus_road_and_rail"])
        folds.append(metric)
    valid=[x for x in folds if x["estimated"]]
    if len(valid)<MIN_COMPARABLE_SOURCE_SECTIONS:
        return {"identified":False,"support":support,
                "reason":"Too few independent original sections with usable held-out folds",
                "heldout_source_section_folds":folds}
    pointmean={label:mean(x["mse_"+label] for x in valid) for label,_ in MODEL_ORDER}
    return {
       "identified":True,
       "source_pair_support":support,
       "heldout_original_rail_sections_total":len(original_sections),
       "n_heldout_sections_estimable":len(valid),
       "source_section_equal_weight_MSE_by_model":pointmean,
       "source_section_mean_incremental_road_MSE_reduction":mean(
           x["road_increment_mse_reduction"] for x in valid),
       "source_section_mean_incremental_rail_after_road_MSE_reduction":mean(
           x["rail_increment_beyond_road_mse_reduction"] for x in valid),
       "n_source_sections_rail_increment_positive":sum(
           x["rail_increment_beyond_road_mse_reduction"]>0 for x in valid),
       "n_source_sections_road_increment_positive":sum(
           x["road_increment_mse_reduction"]>0 for x in valid),
       "source_section_mean_rail_increment_pct_of_road_model_MSE":100*mean(
           x["rail_increment_beyond_road_mse_reduction"] for x in valid
       )/mean(x["mse_plus_road"] for x in valid),
       "details_per_original_section":folds,
       "cross_validation_sites_overlap_between_folds_NOT_independent_experiments":True,
       "no_per_pair_independent_p_values":True,
    }

def rail_partial_slope(records,values=None):
    import numpy as np
    keys=BASE+("road_path_ln",)
    xx=array_features(records,keys)
    x=np.asarray([r["features"]["relative_rail_road_path_ln"] for r in records])
    y=np.asarray([r["jaccard"] for r in records] if values is None else values)
    if y.shape!=(len(records),) or not all(np.isfinite(y)):
        raise ValueError("Original plant Jaccard vector invalid")
    beta=np.linalg.lstsq(xx,x,rcond=None)[0]
    xres=x-xx@beta
    den=xres@xres
    baseline_xvar=((x-x.mean())**2).sum()
    info=den/max(1.e-20,baseline_xvar)
    if den<1.e-10 or info<.0001:
        return {"identified":False,
                "source_rail_ratio_unique_variation_fraction":float(info),
                "reason":"Rail-to-road path ratio has too little independent variation"}
    beta_rail=float(xres@y/den)
    resid_y=y-xx@np.linalg.lstsq(xx,y,rcond=None)[0]
    cond_r2=((xres@resid_y)**2)/(den*(resid_y@resid_y)) if resid_y@resid_y>0 else None
    return {"identified":True,
            "partial_beta_jaccard_per_1unit_ln_rail_over_road":beta_rail,
            "partial_beta_per_doubling_rail_over_road":beta_rail*math.log(2.),
            "partial_association_signed_direction":"negative=more rail-direct paths linked to more original taxon similarity",
            "conditional_residual_R2_IN_SAMPLE_NOT_inference":float(cond_r2) if cond_r2 is not None else None,
            "source_rail_ratio_unique_variation_fraction":float(info),
            "OSM_rail_road_collinearity_weak_positivity_if_unique_fraction_under_0_02":bool(info<.02)}

def original_guild_species_swap_null(records,sites,guild,nchains=NULL_CHAINS,
                                      draws=NULL_DRAWS_PER_CHAIN,
                                      burnin=NULL_BURNIN_PROPOSALS,
                                      spacing=NULL_THIN_PROPOSALS,seed=RANDOM_SEED):
    import numpy as np
    support=supports_estimator(records)
    if not support["supports_main_regression"]:
        return {"identified":False,"reason":"Original site pairs or route sections insufficient","support":support}
    observed=rail_partial_slope(records)
    if not observed["identified"]:
        return {"identified":False,"reason":"Rail ratio is collinear with geography/road; cannot randomize a coefficient","source_beta":observed}
    # xres is fixed, because checkerboard swaps EXACTLY preserve
    # the functional guild count and each original site richness.
    design=array_features(records,BASE+("road_path_ln",))
    rail=np.asarray([r["features"]["relative_rail_road_path_ln"] for r in records])
    xres=rail-design@np.linalg.lstsq(design,rail,rcond=None)[0]
    den=float(xres@xres)
    null=[]
    chain=[]
    for cc in range(nchains):
        mat=make_matrices(sites,guild)
        rng=random.Random(seed+101*cc)
        accepted=swap_proposals(mat,sites,rng,burnin,True)
        draws_in_chain=[]
        for dd in range(draws):
            accepted+=swap_proposals(mat,sites,rng,spacing,True)
            species=matrix_taxon_lists(mat,len(sites))
            sims=np.asarray([jaccard(species[r["i"]],species[r["j"]])
                             for r in records],dtype=float)
            if not np.isfinite(sims).all():
                raise ValueError("Original fixed frequency/guild null generated missing Jaccard")
            val=float(xres@sims/den)
            draws_in_chain.append(val)
        null+=draws_in_chain
        chain.append({"chain":cc+1,"n_correlated_draws":len(draws_in_chain),
                      "accepted_swaps":accepted,
                      "proposed_swaps":burnin+spacing*draws})
    z=np.asarray(null)
    v=observed["partial_beta_jaccard_per_1unit_ln_rail_over_road"]
    return {
      "identified":True,
      "observed_original_conditional_rail_ratio_coefficient":v,
      "source_expected_rail_coefficient_under_year_track_guild_species_fixed_margin":float(np.mean(z)),
      "observed_minus_null_coefficient":float(v-np.mean(z)),
      "MC_null_percentile_2_5":float(np.percentile(z,2.5)),
      "MC_null_percentile_97_5":float(np.percentile(z,97.5)),
      "MC_tail_fraction_lower_plus_one_NOT_EXPERIMENTAL_P":float((1+sum(z<=v))/(1+len(z))),
      "n_autocorrelated_mcmc_draws":len(z),
      "chain_diagnostics":chain,
      "stratified_original_year_track_type_frequencies_preserved":True,
      "site_guild_richness_and_taxon_occupancy_frequency_preserved":True,
      "no_local_regional_species_pool_or_rail_freight_control":True,
      "NOT_per_pair_independent_significance":True,
    }

def run(original,pbf):
    topo,pairs,sites,guild=run_network(original,pbf,return_network_matrices=True)
    all_records,all_omitted=model_records(pairs,sites,False)
    restricted,restricted_omitted=model_records(pairs,sites,True)
    main=section_heldout_validation(restricted)
    sensitivity=section_heldout_validation(all_records)
    coef=rail_partial_slope(restricted) if main["identified"] else {"identified":False,"reason":"Source support inadequate"}
    null=(original_guild_species_swap_null(restricted,sites,guild)
          if main["identified"] and coef["identified"] else
          {"identified":False,"reason":"Predeclared source-support or rank gate failed"})
    return {
      "version":"v15.1_real_topology_road_vs_rail_plant_identity",
      "original_2026_field_and_2024_Geofabrik_source":topo,
      "predeclared_reference_design":"15 original route-section holdouts; pairs not independent",
      "source_network_supported_pairs_before_detour_filter":len(all_records),
      "source_main_path_detour_upper_bound_geographic":"8 x 2024 source geodesic km + 4 km",
      "source_main_eligible_pairs":len(restricted),
      "source_omitted_main":restricted_omitted,
      "main_section_heldout_rail_vs_road_predictive_result":main,
      "untrimmed_network_paths_heldout_sensitivity":sensitivity,
      "partial_rail_vs_road_source_coefficient":coef,
      "actual_original_site_species_guild_frequency_null_conditional_rail_result":null,
      "strict_limits":[
        "Every original sampled plot is on a railway; absence of independent nonrail plots prevents habitat-specific causal rail effects.",
        "OSM physical paths do not verify active trains, freight shipping direction, timetable or likely seed carriage.",
        "Geofabrik whole Slovenia source excludes routes through adjacent countries; detour trimming is predeclared not posthoc.",
        "Site coordinates are ORIGINAL July–September 2022–23 while physical OpenStreetMap extraction is 2024–01; historical map is approximate and not fully contemporaneous.",
        "The 4-guild composition is a coarse species code proxy; no measured aerodynamic dispersal of wild diaspores.",
        "Road-network distance can be a shared geographical collider or correlate of land use; adding it is a negative-control diagnostic and NOT rail randomization.",
        "Heldout source sections share partner sites and national source pool; 15 fold error values are descriptive not independent country studies.",
        "Null Monte Carlo checkerboard draws retain site guild richness and species prevalence by year/track type, not spatially varying regional pools or initial sampling selection.",
        "Cross-sectional taxon identity cannot distinguish arrival vs establishment vs local adaptation/plasticity; no mowing, measured airflow, fitness or heritability."
      ],
      "no_conclusion_of_rail_carried_seeds_or_evolved_stress_traits":True
    }

def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","osm-pbf","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    source=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    result=run(source,a.osm_pbf)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    report={k:v for k,v in result.items()
            if k in ("version","source_main_eligible_pairs",
                     "source_network_supported_pairs_before_detour_filter",
                     "source_omitted_main",
                     "main_section_heldout_rail_vs_road_predictive_result",
                     "partial_rail_vs_road_source_coefficient",
                     "actual_original_site_species_guild_frequency_null_conditional_rail_result")}
    print(json.dumps(report,ensure_ascii=False))

if __name__=="__main__":main()

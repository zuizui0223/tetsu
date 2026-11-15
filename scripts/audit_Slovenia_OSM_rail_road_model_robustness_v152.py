#!/usr/bin/env python3
"""v15.2 POSTHOC robustness of genuine Slovenian 89-site rail-vs-road result.

This is a sensitivity analysis AFTER seeing v15.1's null held-out result.
Do NOT label any improvement selected here as a pre-registered independent
positive test. 2024 original Geofabrik OSM PBF SHA pinned; all original
vegetation plots from Šabić et al 2026 SHA pinned by original reader.
Reuse exact v15.1 pair support; change ONLY known vs all author taxa,
model conditioning controls and fixed ridge regularization alpha.

Why this matters: same-original-section membership, functional-guild
composition and subset of taxa with complete trait codes can influence
the inference. Recompute every candidate, report full grid, not only
positive subsets. No cargo OD or mowing timing; all observational.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import jaccard
from analyze_Slovenia_true_OSM_rail_road_paths_v150 import run_network
from analyze_Slovenia_rail_road_heldout_species_identity_v151 import (
    model_records,array_features,ridge_train_predict,BASE
)

EXPECTED_PBF_SHA="b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39"
EXPECTED_MAIN_N_PAIRS=3896
EXPECTED_ORIGINAL_V151_ROAD_MSE=0.0037366620160146245
EXPECTED_ORIGINAL_V151_FULL_MSE=0.0037567314602607142

ALPHAS=(1.,10.,100.)
SPECIFICATIONS={
    "strict_v151_source_section_and_functional":BASE,
    "sensitivity_omit_same_original_section":
         tuple(k for k in BASE if k!="same_original_section"),
    "sensitivity_omit_all_functional_profile_controls":
         tuple(k for k in BASE if k!="functional_overlap" and not k.startswith("fun_bin_")),
    "sensitivity_omit_observed_site_species_richness":
         tuple(k for k in BASE if not k.startswith("site_log_richness_")),
}

def new_outcome(rows,sites,all_taxa):
    if not all_taxa:return rows
    result=[]
    for r in rows:
        i,j=r["i"],r["j"]
        x,y=sites[i]["all_original_taxa"],sites[j]["all_original_taxa"]
        if not x or not y:raise ValueError("Actual original source all-taxon floras incomplete")
        nr=dict(r)
        nr["jaccard"]=jaccard(x,y)
        nr["features"]=dict(r["features"])
        nr["features"]["site_log_richness_mean"]=(math.log(len(x))+math.log(len(y)))/2
        nr["features"]["site_log_richness_absdiff"]=abs(math.log(len(x))-math.log(len(y)))
        result.append(nr)
    return result

def heldout_variant(rows,basekeys,alpha):
    import numpy as np
    y=np.asarray([r["jaccard"] for r in rows],dtype=float)
    roadkeys=tuple(basekeys)+("road_path_ln",)
    fullkeys=roadkeys+("relative_rail_road_path_ln",)
    matrices={
        "road_control":array_features(rows,roadkeys),
        "road_plus_rail":array_features(rows,fullkeys)
    }
    sections=sorted({r["section_i"] for r in rows}|{r["section_j"] for r in rows})
    folds=[]
    for s in sections:
        train=np.asarray([i for i,r in enumerate(rows)
                         if s not in (r["section_i"],r["section_j"])],dtype=int)
        test=np.asarray([i for i,r in enumerate(rows)
                        if s in (r["section_i"],r["section_j"])],dtype=int)
        if len(train)<40 or len(test)<5:
            folds.append({"source_section":s,"estimated":False,
                          "n_train":len(train),"n_test":len(test)})
            continue
        scores={}
        for name,X in matrices.items():
            pred=ridge_train_predict(X,y,train,test,alpha=alpha)
            scores[name]=float(np.mean((pred-y[test])**2))
        folds.append({"source_section":s,"estimated":True,
                     "n_train":len(train),"n_test":len(test),
                     "road_mse":scores["road_control"],
                     "road_plus_rail_mse":scores["road_plus_rail"],
                     "rail_gain":scores["road_control"]-scores["road_plus_rail"]})
    valid=[v for v in folds if v["estimated"]]
    if len(valid)<8:
        return {"identified":False,"n_heldout_source_sections":len(valid)}
    road_mse=mean(v["road_mse"] for v in valid)
    full_mse=mean(v["road_plus_rail_mse"] for v in valid)
    return {
        "identified":True,
        "n_heldout_source_sections":len(valid),
        "source_section_equal_weight_road_MSE":road_mse,
        "source_section_equal_weight_road_plus_rail_MSE":full_mse,
        "source_section_equal_weight_rail_increment_MSE_improvement":road_mse-full_mse,
        "source_section_equal_weight_rail_increment_MSE_improvement_pct":100*(road_mse-full_mse)/road_mse,
        "n_source_sections_better_with_rail":sum(v["rail_gain"]>0 for v in valid),
        "n_source_sections_worse_with_rail":sum(v["rail_gain"]<0 for v in valid),
        "no_holdout_per_pair_independent_p_values":True,
        "n_original_source_sections_not_countries":len(valid),
    }

def sensitivity_grid(rows,sites):
    full_result=[]
    for outcome,all_taxa in (("261_trait_known_author_taxa",False),
                             ("ALL_original_author_observed_BS_BT_taxa",True)):
        data=new_outcome(rows,sites,all_taxa)
        for spec,basekeys in SPECIFICATIONS.items():
            for alpha in ALPHAS:
                stat=heldout_variant(data,basekeys,alpha)
                full_result.append({"taxon_outcome":outcome,
                                    "conditioning_specification":spec,
                                    "ridge_alpha_prespecified_grid":alpha,
                                    **stat})
    return full_result

def run(original,pbf):
    network,annotated,sites,guild=run_network(original,pbf,True)
    if network["original_OSM_pbf_SHA256"]!=EXPECTED_PBF_SHA:
        raise ValueError("v15.2 original Geofabrik PBF SHA does not match verified 2024 archive")
    rows,omitted=model_records(annotated,sites,True)
    if len(rows)!=EXPECTED_MAIN_N_PAIRS or len(sites)!=89:
        raise ValueError("v15.2 strict original v15.1 site/pair support changed")
    results=sensitivity_grid(rows,sites)
    frozen=[x for x in results if
            x["taxon_outcome"]=="261_trait_known_author_taxa" and
            x["conditioning_specification"]=="strict_v151_source_section_and_functional" and
            x["ridge_alpha_prespecified_grid"]==10.0]
    if len(frozen)!=1:
        raise ValueError("Frozen v15.1 reference model missing")
    x=frozen[0]
    if (not math.isclose(x["source_section_equal_weight_road_MSE"],
                         EXPECTED_ORIGINAL_V151_ROAD_MSE,rel_tol=1.e-11)
        or not math.isclose(x["source_section_equal_weight_road_plus_rail_MSE"],
                            EXPECTED_ORIGINAL_V151_FULL_MSE,rel_tol=1.e-11)):
        raise ValueError("Exact SOURCE v15.1 holdout result changed during sensitivity audit")
    positive=[x for x in results if x["identified"] and
              x["source_section_equal_weight_rail_increment_MSE_improvement"]>0]
    return {
      "version":"v15.2_POSTHOC_real_original_robustness_exhaustive_grid",
      "source_original_author_DOI":network["original_vegetation_DOI"],
      "source_original_author_supplement_SHA256":network["original_vegetation_source_SHA256"],
      "fixed_historical_2024_OSM_pbf_SHA256":EXPECTED_PBF_SHA,
      "source_same_v151_site_count":len(sites),
      "source_same_v151_strict_eligible_site_pairs":len(rows),
      "original_15_source_sections_not_15_country_studies":True,
      "exact_v151_main_reference_successfully_reproduced":True,
      "posthoc_grid_exhaustively_reported_not_fitted_to_positives":True,
      "n_sensitivity_combinations":len(results),
      "n_sensitivity_specs_where_rail_adds_section_average_prediction":len(positive),
      "n_sensitivity_specs_where_rail_worsens_or_equal_section_average_prediction":len(results)-len(positive),
      "actual_source_outcome_conditioning_ridge_sensitivity":results,
      "interpretation_limits":[
        "This is POSTHOC after original v15.1 source result; positive subsets cannot replace frozen main model.",
        "All site pair observations share only 89 original locations in 15 sections nested in ONE country.",
        "No measured local seed rainfall, propagation success, regional habitat pool, cutting exposure or independently confirmed freight itinerary.",
        "Rail ratio is based on physically mapped 2024 OSM graph, not biological travel, train number, active schedule, service time, direction or actual observed train-carried seeds.",
        "Covariates include original guild traits partly derived from the outcome site species lists; the CV therefore tests conditional additional source species identity sharing rather than purely out-of-sample prediction from external abiotic predictors.",
        "All-taxon source outcome retains published missing-trait taxa as true observed plants, not as nonwind control.",
        "No invented original unsampled station/track or rail vs nonrail experiment."
      ]
    }

def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","osm-pbf","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    source=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    d=run(source,a.osm_pbf)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:v for k,v in d.items() if
                      k in ("version","exact_v151_main_reference_successfully_reproduced",
                            "source_same_v151_site_count","source_same_v151_strict_eligible_site_pairs",
                            "n_sensitivity_combinations",
                            "n_sensitivity_specs_where_rail_adds_section_average_prediction",
                            "n_sensitivity_specs_where_rail_worsens_or_equal_section_average_prediction",
                            "actual_source_outcome_conditioning_ridge_sensitivity")},ensure_ascii=False))
if __name__=="__main__":main()

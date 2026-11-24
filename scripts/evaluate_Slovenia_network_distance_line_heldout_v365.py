#!/usr/bin/env python3
"""v365: can TRUE physical Overture railway-network distance predict genuine
Slovenian plant community Jaccard beyond geographical distance?

Full original section-heldout test, NOT pair-random CV, not proof of
train dispersal or globally replicated biology. True 2026 source GERS
normalized WGS84-geodesic lengths corrected in v364 before fitting.
No original botanical site names, taxon records, GPS or rail geometry exported.
"""
from __future__ import annotations
import argparse,json,math,statistics
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import (
    to_original_sites,jaccard_similarity
)
from analyze_Slovenia_Overture_connector_network_botany_v358 import (
    SOURCE_UTM,SOURCE_ORIGINAL_SHA,RELEASE,real_overture
)
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import (
    physical_network_distances
)
from recompute_Slovenia_GERS_true_geodesic_network_v364 import (
    corrected_physical_graph,corrected_source_attachments
)

LAM=0.05
HOLDOUT_MIN_TRAIN_PAIRS=200
HOLDOUT_MIN_TEST_PAIRS=10
HOLDOUT_MIN_SOURCE_SECTIONS=5
SOURCE_N_SITES=89
SOURCE_N_SECTIONS=15
BASE_FEATURES=("log_geodesic","log_geodesic_squared","same_year","same_layout")
RICH_FEATURES=("mean_log_site_richness","abs_diff_log_site_richness")
NETWORK_FEATURE="log_physical_network"
SECTION_FEATURE="same_author_section"
MODELS={
   "M0_geo_year_layout_richness":BASE_FEATURES+RICH_FEATURES,
   "M1_M0_plus_TRUE_network_distance":BASE_FEATURES+RICH_FEATURES+(NETWORK_FEATURE,),
   "M2_M0_plus_author_section_label":BASE_FEATURES+RICH_FEATURES+(SECTION_FEATURE,),
   "M3_M2_plus_TRUE_network_distance":BASE_FEATURES+RICH_FEATURES+(NETWORK_FEATURE,SECTION_FEATURE)
}
MODELS_NO_RICH={
   "M0_geo_year_layout":BASE_FEATURES,
   "M1_geo_plus_TRUE_network_distance":BASE_FEATURES+(NETWORK_FEATURE,),
   "M2_geo_plus_author_section_label":BASE_FEATURES+(SECTION_FEATURE,),
   "M3_geo_section_plus_TRUE_network_distance":BASE_FEATURES+(NETWORK_FEATURE,SECTION_FEATURE)
}
PANELS=("BS_and_BT","BS","BT")

def source_original_pairs(sites,network):
    """ONE undirected record per real botanical site pair; source missing paths excluded."""
    from pyproj import Geod
    g=Geod(ellps="WGS84")
    pairs=[];qa=Counter()
    for i,j in combinations(range(len(sites)),2):
        qa["source_original_site_pairs_total"]+=1
        net=network.get((i,j))
        if net is None:
            qa["source_no_verified_source_GERS_network_path_pairs"]+=1
            continue
        rev=network.get((j,i))
        if rev is None or abs(rev-net)>1e-8:
            raise ValueError("Source GERS railway shortest distance is not symmetric")
        a,b=sites[i],sites[j]
        geod=abs(float(g.inv(a["lon"],a["lat"],b["lon"],b["lat"])[2]))/1000.
        if not all(math.isfinite(v) and v>=0 for v in (geod,net)):
            raise ValueError("Source original field GPS geodesic or GERS route length malformed")
        qa["original_true_source_physical_network_pair_n"]+=1
        if net+0.03<geod: qa["source_real_rail_route_shorter_than_geo_by_more_than30m_n"]+=1
        row={"i":i,"j":j,"source_original_geodetic_km":geod,
           "source_true_rail_GERS_network_km":net,
           "same_original_survey_year":a["year"]==b["year"],
           "same_original_track_layout":a["layout"]==b["layout"],
           "same_original_author_section":a["section"]==b["section"],
           "site_source_field_richness":(len(a["BS"]|a["BT"]),len(b["BS"]|b["BT"])),
           "original_author_section_pair":(a["section"],b["section"]),
           "response":{}}
        for habitat in PANELS:
            def src(s):
                return s["BS"]|s["BT"] if habitat=="BS_and_BT" else s[habitat]
            row["response"][habitat]=jaccard_similarity(src(a),src(b))
        pairs.append(row)
    if len(pairs)<800 or len(pairs)+qa["source_no_verified_source_GERS_network_path_pairs"]!=qa["source_original_site_pairs_total"]:
        raise ValueError("Source genuine network-supported original botanical pairs too sparse/drifted")
    return pairs,dict(qa)

def features(row):
    geo=math.log1p(row["source_original_geodetic_km"])
    lng=math.log1p(row["source_true_rail_GERS_network_km"])
    ri,rj=row["site_source_field_richness"]
    if min(ri,rj)<1:raise ValueError("Original source field site richness absent")
    r1=math.log1p(ri);r2=math.log1p(rj)
    return {
        "log_geodesic":geo,
        "log_geodesic_squared":geo**2,
        "same_year":float(row["same_original_survey_year"]),
        "same_layout":float(row["same_original_track_layout"]),
        "mean_log_site_richness":.5*(r1+r2),
        "abs_diff_log_site_richness":abs(r1-r2),
        "log_physical_network":lng,
        "same_author_section":float(row["same_original_author_section"])
    }

def source_degree_pair_weights(pairs):
    """Same original site contributes through many dependent pairs; normalize."""
    degree=Counter()
    for r in pairs:
        degree[r["i"]]+=1;degree[r["j"]]+=1
    out=[.5*(1/degree[r["i"]]+1/degree[r["j"]]) for r in pairs]
    S=sum(out)
    if S<=0:raise ValueError("Original source pair weights zero")
    return [x/S for x in out]

def fit_predict(train,test,habitat,featureset,lam=LAM):
    """Same frozen weighted ridge penalty across models; no hyperparameter search."""
    import numpy as np
    if lam<0 or len(train)<2 or len(test)<1:raise ValueError("Cannot fit source original empty field species pairs")
    columns=tuple(featureset)
    x=np.asarray([[features(row)[k] for k in columns] for row in train],dtype=float)
    xt=np.asarray([[features(row)[k] for k in columns] for row in test],dtype=float)
    target=np.asarray([row["response"][habitat] for row in train],dtype=float)
    w=np.asarray(source_degree_pair_weights(train),dtype=float)
    mu=np.sum(w[:,None]*x,axis=0)
    sd=np.sqrt(np.sum(w[:,None]*(x-mu)**2,axis=0))
    sd=np.where(sd>1e-10,sd,1.0)
    X=np.concatenate([np.ones((len(train),1)),(x-mu)/sd],axis=1)
    T=np.concatenate([np.ones((len(test),1)),(xt-mu)/sd],axis=1)
    penalty=np.eye(X.shape[1]);penalty[0,0]=0.
    # Weights sum to one, giving fixed source-ridge penalty an invariant scale.
    lhs=X.T@(w[:,None]*X)+lam*penalty
    rhs=X.T@(w*target)
    coef=np.linalg.solve(lhs,rhs)
    pred=np.clip(T@coef,0.,1.)
    return pred

def line_holdout_fit(pairs,source_sites,habitat,models):
    from numpy import mean,asarray
    sections=sorted({s["section"] for s in source_sites})
    if len(sections)!=SOURCE_N_SECTIONS:raise ValueError("Source original 15 rail section groups changed")
    folds=[];source_excluded=Counter()
    for section in sections:
        train=[r for r in pairs if section not in r["original_author_section_pair"]]
        test=[r for r in pairs if section in r["original_author_section_pair"]]
        if len(train)<HOLDOUT_MIN_TRAIN_PAIRS or len(test)<HOLDOUT_MIN_TEST_PAIRS:
            source_excluded["source_author_section_with_low_pair_support"]+=1
            continue
        expected=asarray([r["response"][habitat] for r in test],dtype=float)
        out={}
        for model,keys in models.items():
            p=fit_predict(train,test,habitat,keys)
            out[model]={"test_mse":float(mean((expected-p)**2)),
                        "test_mae":float(mean(abs(expected-p)))}
        folds.append({"heldout_original_source_section_index":len(folds)+1,
          "source_heldout_site_n":sum(s["section"]==section for s in source_sites),
          "source_train_dependent_pair_records_not_independent":len(train),
          "source_test_dependent_pair_records_not_independent":len(test),
          "test_losses":out})
    if len(folds)<HOLDOUT_MIN_SOURCE_SECTIONS:
        return {"identified":False,"source_original_section_holdout_n":len(folds),
            "source_not_identified":True}
    average={model:{"source_equal_original_section_mean_test_MSE":
                   statistics.mean(f["test_losses"][model]["test_mse"] for f in folds),
                    "source_equal_original_section_mean_test_MAE":
                   statistics.mean(f["test_losses"][model]["test_mae"] for f in folds)}
             for model in models}
    names=list(models)
    if len(names)!=4:raise ValueError("Models not pre-frozen four family sequence")
    contrasts={
       "M0_minus_M1_TRUE_network_increment_MSE":(
           average[names[0]]["source_equal_original_section_mean_test_MSE"]-
           average[names[1]]["source_equal_original_section_mean_test_MSE"]),
       "M2_minus_M3_TRUE_network_after_line_indicator_MSE":(
           average[names[2]]["source_equal_original_section_mean_test_MSE"]-
           average[names[3]]["source_equal_original_section_mean_test_MSE"]),
       "M0_minus_M2_original_author_section_label_increment_MSE":(
           average[names[0]]["source_equal_original_section_mean_test_MSE"]-
           average[names[2]]["source_equal_original_section_mean_test_MSE"])
    }
    deltas=[f["test_losses"][names[0]]["test_mse"]-f["test_losses"][names[1]]["test_mse"]
            for f in folds]
    dl2=[f["test_losses"][names[2]]["test_mse"]-f["test_losses"][names[3]]["test_mse"]
            for f in folds]
    contrasts.update({
       "source_original_sections_where_network_improves_MSE_without_line_n":sum(d>0 for d in deltas),
       "source_original_sections_where_network_worsens_MSE_without_line_n":sum(d<0 for d in deltas),
       "source_original_sections_where_network_improves_MSE_after_line_n":sum(d>0 for d in dl2),
       "source_original_sections_where_network_worsens_MSE_after_line_n":sum(d<0 for d in dl2),
       "technical_per_section_M0_minus_M1_MSE_median":statistics.median(deltas),
       "technical_per_section_M2_minus_M3_MSE_median":statistics.median(dl2),
       "whole_original_line_heldout_pairs_unequal_and_NOT_independent":True
    })
    return {"identified":True,
      "source_heldout_author_original_sections_n":len(folds),
      "source_sections_below_source_minimum_support":source_excluded[
          "source_author_section_with_low_pair_support"],
      "source_no_training_field_site_shared_with_heldout_line":True,
      "source_fixed_ridge_penalty_each_model":LAM,
      "original_plant_site_richness_known_conditioned_on_in_primary":RICH_FEATURES[0] in models[names[0]],
      "source_n_original_survey_years":2,
      "source_models_test_losses":average,
      "source_predeclared_network_incremental_validation":contrasts,
      "source_fold_support_and_test_losses":[{"original_author_line_fold_index":f["heldout_original_source_section_index"],
        "n_source_heldout_botanical_sites":f["source_heldout_site_n"],
        "train_pairs_NOT_field_independent":f["source_train_dependent_pair_records_not_independent"],
        "test_pairs_NOT_field_independent":f["source_test_dependent_pair_records_not_independent"],
        "test_losses":f["test_losses"]} for f in folds],
      "not_trains_seed_transport_or_globally_replicated_evidence":True}

def run(original,rail_factory=real_overture):
    from pyproj import Transformer
    bysection,botqa=to_original_sites(original)
    sites=[s for sec in sorted(bysection) for s in bysection[sec]]
    if len(sites)!=SOURCE_N_SITES or len(bysection)!=SOURCE_N_SECTIONS:
        raise ValueError("Authentic original botanist 89 sites/15 sections missing")
    bbox=(min(s["lon"] for s in sites)-.12,min(s["lat"] for s in sites)-.12,
          max(s["lon"] for s in sites)+.12,max(s["lat"] for s in sites)+.12)
    real,source=rail_factory(bbox)
    T=Transformer.from_crs("EPSG:4326",SOURCE_UTM,always_xy=True)
    _,G,segments,qa,metricqa=corrected_physical_graph(real,T)
    if len(segments)!=6612 or G.number_of_edges()!=14181:
        raise ValueError("Original v364 real source rail topology differs")
    augmented,site,attachQA=corrected_source_attachments(G,segments,sites,T)
    status=Counter(x["status"] for x in site)
    if status["matched_to_real_rail_connector_network"]!=86 or status[
        "ambiguous_multiple_disconnected_source_tracks"]!=3:
        raise ValueError("Original source v364 botany site-to-track support drift")
    d=physical_network_distances(augmented,site)
    pairs,sourcepairqa=source_original_pairs(sites,d)
    summary={}
    for habitat in PANELS:
        summary[habitat]={
            "original_source_primary_effort_conditional":line_holdout_fit(
               pairs,sites,habitat,MODELS),
            "no_plant_richness_predictor_sensitivity":line_holdout_fit(
               pairs,sites,habitat,MODELS_NO_RICH)
        }
    return {"version":"v365-real-Slovenia-plant-beta-physical-rail-network-vs-geography-LOLO",
        "source_original_2026_Slovenia_publisher_DOI":"10.1007/s44473-026-00155-z",
        "source_original_publisher_401279B_SHA256":SOURCE_ORIGINAL_SHA,
        "source_real_Overture_release":RELEASE,
        "source_original_field_botanical_site_n":len(sites),
        "source_original_line_section_heldout_fold_n":len(bysection),
        "source_actual_original_track_side_GERS_matched_n":status["matched_to_real_rail_connector_network"],
        "source_ambiguous_actual_physical_track_match_n":status["ambiguous_multiple_disconnected_source_tracks"],
        "source_original_GERS_geodesic_graph_QA":metricqa,
        "source_original_botanical_site_geodetic_linear_ref":attachQA,
        "source_pair_opportunity_QA":sourcepairqa,
        "source_plant_similarity_source_habitats":summary,
        "source_rail_route_disconnected_site_pairs_are_unavailable_not_0":True,
        "source_only_one_country_author_sections_not_independent_global_replication":True,
        "source_true_RAIL_GERS_path_not_real_train_OD_or_verified_active_track":True,
        "source_line_folds_not_dependent_pair_level_p_values":True,
        "source_genetic_adaptation_seed_transport_causal_effect_unmeasured":True,
        "source_raw_botanical_GPS_taxon_plot_matrix_and_rail_geometry_not_republished":True}

def main():
    p=argparse.ArgumentParser()
    for k in ("original-xlsx","source-audit","output"):p.add_argument("--"+k,required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.original_xlsx,a.source_audit)
    obj=run(original)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    output={"version":obj["version"],
        "source_pair_opportunity_QA":obj["source_pair_opportunity_QA"],
        "source_original_line_section_heldout_fold_n":obj["source_original_line_section_heldout_fold_n"],
        "source_results":{k:{p:{
            "identified":v["identified"],
            "n_source_heldout_sections":v.get("source_heldout_author_original_sections_n"),
            "losses":v.get("source_models_test_losses"),
            "predictive_increments":v.get("source_predeclared_network_incremental_validation")}
            for p,v in panels.items()} for k,panels in obj["source_plant_similarity_source_habitats"].items()}}
    print(json.dumps(output,ensure_ascii=False),flush=True)

if __name__=="__main__":main()

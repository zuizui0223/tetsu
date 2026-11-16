#!/usr/bin/env python3
"""v15.6 genuine Czech source-35 railway station (2020–22 flora) comparison
of physical 2023 OSM railway route length vs physical 2023 road route length.

The botanical station coordinates are **2026 SR70 point proxies** independently
matched by ORIGINAL EXACT author station name; not source field quadrats.
Only 35/39 authors' original stations have strict exact-name correspondence.
Original source units: 5 geographical rail stretches in ONE country, not 595
independent survey replicates. Original Czech station scale differs from the
Slovenian 100m sites; do not pool common meta railway beta.

Outcome is 2024 original author presence-by-station taxonomic Jaccard of 544
Pladias original dual-known species (and 761 taxon all-original sensitivity).
Question: does adding physical shortest-rail/road ratio help predict heldout
ORIGINAL Czech stretch Jaccard after geodesic distance, functional guild
profile, station size, source species richness, and ACTUAL road route length?

Original 2023 OSM node-based paths do not represent operating 2020–22 train
freight, carriage of wild seeds, measured winds or cutting selection.
No non-rail habitat plots. Stop if route overlap or source independence fails.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter
from pathlib import Path
from statistics import mean
from analyze_Czech_2023_true_rail_road_osm_routes_v155 import run as network_run
from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import real_site_pairs
from analyze_Slovenia_rail_road_heldout_species_identity_v151 import ridge_train_predict
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    into_bin,DISTANCE_CUTS_KM,FUNCTIONAL_OVERLAP_CUTS
)
from analyze_Czech_2024_39station_functional_identity_v143 import function_sim,jaccard

ALPHA=10.
MIN_ROUTES=250
MIN_SITES=30
MIN_STRETCHES=4
DET_MUL=8.
DET_ADD=4.
MIN_GEO_KM=2.


def pair_records(topology_pairs,sites,guild,taxon_mode="both_Pladias_known"):
    if taxon_mode not in ("both_Pladias_known","all_original_S1_current"):
        raise ValueError("Czech original plant-taxonomic response group not preregistered")
    out=[];missing=Counter()
    for p in topology_pairs:
        i,j=p["i"],p["j"]
        a,b=sites[i],sites[j]
        if not p["network_route_comparison_eligible"]:
            missing["no_geography_and_both_true_rail_road_route_support"]+=1
            continue
        geo=p["source_station_geodesic_km"]
        rail=p["rail_graph_source_shortest_path_km"]
        road=p["road_graph_source_shortest_path_km"]
        if not all(v is not None and math.isfinite(v) and v>0 for v in (geo,rail,road)):
            missing["missing_or_nonpositive_real_route"]+=1
            continue
        if max(rail,road)>DET_MUL*geo+DET_ADD:
            missing["predeclared_extreme_route_detour"]+=1
            continue
        if not all(s["area_ha"] is not None and s["area_ha"]>0 for s in (a,b)):
            missing["original_publisher_station_area_missing"]+=1
            continue
        ga,gb=a["known"],b["known"]
        f=function_sim(ga,gb,guild)
        if f is None:raise ValueError("Original full Pladias matched site guild profile missing")
        geobin=into_bin(geo,DISTANCE_CUTS_KM)
        fbin=into_bin(f,FUNCTIONAL_OVERLAP_CUTS)
        if geobin is None or fbin is None:raise ValueError("Original SI fixed spatial+guild bins not supported")
        x=a["all_original_taxa"] if taxon_mode=="all_original_S1_current" else a["known"]
        y=b["all_original_taxa"] if taxon_mode=="all_original_S1_current" else b["known"]
        lna,lnb=math.log(len(x)),math.log(len(y))
        aa,ab=math.log(a["area_ha"]),math.log(b["area_ha"])
        feat={
          "geo_ln":math.log(geo),
          "functional_overlap":f,
          "same_original_stretch":float(p["source_same_stretch"]),
          "station_area_ln_mean":.5*(aa+ab),
          "station_area_ln_absdiff":abs(aa-ab),
          "site_ln_species_count_mean":.5*(lna+lnb),
          "site_ln_species_count_absdiff":abs(lna-lnb),
          "rail_snap_ln":math.log1p(p["source_both_rail_snap_km"]),
          "road_snap_ln":math.log1p(p["source_both_road_snap_km"]),
          "road_route_ln":math.log(road),
          "rail_to_road_route_ln":math.log(rail/road),
        }
        for k in range(1,7):feat["geobin_"+str(k)]=float(geobin==k)
        for k in range(1,4):feat["funbin_"+str(k)]=float(fbin==k)
        out.append({"i":i,"j":j,"source_section_i":a["section"],"source_section_j":b["section"],
                    "jaccard":jaccard(x,y),"features":feat,
                    "site_actual_original_botanical_station_geometry_is_proxy":True})
    return out,dict(missing)


BASE=("geo_ln","functional_overlap","same_original_stretch",
      "station_area_ln_mean","station_area_ln_absdiff",
      "site_ln_species_count_mean","site_ln_species_count_absdiff",
      "rail_snap_ln","road_snap_ln")+tuple("geobin_"+str(i) for i in range(1,7))+tuple(
          "funbin_"+str(i) for i in range(1,4))
MODELS=(("M0_geography_functional_area_station_size",BASE),
        ("M1_plus_actual_2023_road_route",BASE+("road_route_ln",)),
        ("M2_plus_actual_2023_road_and_rail_ratio",BASE+("road_route_ln","rail_to_road_route_ln")))


def cv_section_heldout(records):
    import numpy as np
    distinctsites=set()
    groups=set()
    for r in records:
        distinctsites.update((r["i"],r["j"]))
        groups.update((r["source_section_i"],r["source_section_j"]))
    support={
       "n_dependent_matched_source_station_pairs":len(records),
       "n_source_original_botanical_stations":len(distinctsites),
       "n_original_rail_stretches":len(groups),
       "minimum_required_route_pairs":MIN_ROUTES,
       "minimum_required_source_stations":MIN_SITES,
       "minimum_required_original_stretches":MIN_STRETCHES,
    }
    if len(records)<MIN_ROUTES or len(distinctsites)<MIN_SITES or len(groups)<MIN_STRETCHES:
        return {"identified":False,"source_overlap_support":support,
                "reason":"Insufficient authentic station location×physical rail/road pathway common support"}
    y=np.asarray([x["jaccard"] for x in records],dtype=float)
    folds=[]
    for group in sorted(groups):
        train=np.asarray([i for i,r in enumerate(records)
                          if group not in (r["source_section_i"],r["source_section_j"])],dtype=int)
        test=np.asarray([i for i,r in enumerate(records)
                         if group in (r["source_section_i"],r["source_section_j"])],dtype=int)
        if len(train)<40 or len(test)<10:
            folds.append({"heldout_original_source_stretch":group,"identified":False,
                          "n_training_pairs":len(train),"n_test_pairs":len(test)})
            continue
        row={"heldout_original_source_stretch":group,"identified":True,
             "n_training_pairs":len(train),"n_test_pairs":len(test)}
        for name,keys in MODELS:
            X=np.asarray([[1.]+[r["features"][k] for k in keys] for r in records],dtype=float)
            pred=ridge_train_predict(X,y,train,test,alpha=ALPHA)
            row["MSE_"+name]=float(np.mean((pred-y[test])**2))
        row["road_increment_MSE_gain"]=(row["MSE_"+MODELS[0][0]]-
                                         row["MSE_"+MODELS[1][0]])
        row["rail_after_road_increment_MSE_gain"]=(row["MSE_"+MODELS[1][0]]-
                                         row["MSE_"+MODELS[2][0]])
        folds.append(row)
    valid=[f for f in folds if f["identified"]]
    if len(valid)<MIN_STRETCHES:
        return {"identified":False,"source_overlap_support":support,
                "reason":"Under minimum 4 original section heldout folds; not identified",
                "folds":folds}
    avg={name:mean(f["MSE_"+name] for f in valid) for name,_ in MODELS}
    relative=100*(avg[MODELS[1][0]]-avg[MODELS[2][0]])/avg[MODELS[1][0]]
    return {"identified":True,"source_overlap_support":support,
            "ridge_alpha_fixed":ALPHA,
            "n_holdout_original_stretches":len(valid),
            "equal_original_stretch_weighted_mse":avg,
            "rail_added_MSE_improvement_percentage_after_road":relative,
            "road_added_MSE_improvement_percentage":100*(avg[MODELS[0][0]]-avg[MODELS[1][0]])/avg[MODELS[0][0]],
            "rail_improves_how_many_stretches":sum(f["rail_after_road_increment_MSE_gain"]>0 for f in valid),
            "road_improves_how_many_stretches":sum(f["road_increment_MSE_gain"]>0 for f in valid),
            "source_stretch_holdout":folds,
            "source_original_5_stretches_NOT_five_countries_or_5_independent_studies":True,
            "model_contains_original_site_species_richness_so_not_abiotic_only_prediction":True}


def rail_coefficient(records):
    import numpy as np
    x=np.asarray([r["features"]["rail_to_road_route_ln"] for r in records])
    y=np.asarray([r["jaccard"] for r in records])
    Z=np.asarray([[1.]+[r["features"][k] for k in BASE+("road_route_ln",)] for r in records])
    xx=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0]
    yy=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0]
    den=xx@xx
    if den<=1.e-11:return {"identified":False,"reason":"True rail and road path predictors collinear"}
    beta=float(xx@yy/den)
    r2=float((xx@yy)**2/(den*(yy@yy))) if yy@yy>0 else None
    return {"identified":True,"source_rail_over_road_path_logratio_partial_beta":beta,
            "source_delta_Jaccard_per_doubling_of_rail_road_ratio":beta*math.log(2),
            "source_conditional_residual_R2_DESCRIPTIVE":r2,
            "no_plant_seed_carriage_or_selection_causality_from_coefficient":True}


def run(s1,s2,disp,life,sr70,pbf):
    topo,pairs,sites,guild=network_run(s1,s2,disp,life,sr70,pbf,True)
    outcomes={}
    for key in ("both_Pladias_known","all_original_S1_current"):
        rows,missing=pair_records(pairs,sites,guild,key)
        cv=cv_section_heldout(rows)
        coef=(rail_coefficient(rows) if cv["identified"] else
              {"identified":False,"reason":"source common-support threshold not met"})
        outcomes[key]={
            "source_original_railroad_pair_support":len(rows),
            "excluded_pair_causes":missing,
            "original_stretch_heldout_CV":cv,
            "descriptive_partial_rail_after_road_source_association":coef,
        }
    return {
      "version":"v15.6_Czech2023_real_rail_road_conditional_station_flora_CV",
      "original_physical_network_provenance":topo,
      "source_2020_2022_Czech_2026_station_centroid_35_original_5_sections":True,
      "source_frozen_geographic_guild_control_inference":outcomes,
      "strict_original_assumptions":[
         "Source botanical plant station records from original 2020-2022, while station point centroids are 2026 SR70 derivative, historic 2023 OSM routes.",
         "Not true author field GPS and not standardized sample-area rail verge plots.",
         "Original plant 35 stations nested in 5 rail stretches; dependent station pairs/5 folds are NOT independent science experiments.",
         "There are no nonrail original vegetation sites, actual cargo OD, measured mowing/wind perturbations, trait heritability, seed catch-to-establishment or natural selection.",
         "Country cross-border rail/road shortest paths may be truncated; fixed extreme-detour gate excludes unsupported examples.",
         "Ridge alpha 10 and original source geodesic and functional bin cutoffs inherited from Slovenian earlier study, not tuned for Czech observed positive findings.",
         "No crosscountry pooled causal railway connectivity effect; preserve 2 independent source vegetation studies' heterogeneous units."
      ],
      "NO_claim_of_train_seed_disposal_or_adaptation":True
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("s1","s2","pladias-diaspore","pladias-lifeform","sr70","osm-pbf","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.pladias_diaspore,a.pladias_lifeform,a.sr70,a.osm_pbf)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k in (
        "version","source_frozen_geographic_guild_control_inference"
    )},ensure_ascii=False))


if __name__=="__main__":main()

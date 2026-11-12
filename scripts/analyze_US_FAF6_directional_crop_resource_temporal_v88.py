#!/usr/bin/env python3
"""v8.8 real source-only two-commodity 2022 FAF6 directionally estimated
RAIL OD vs 2018–21 origin recorded USA-GRIIS taxon candidate that appears
in 2023–25 destination GBIF indexed records.

Crop SCTG02/03 versus resource SCTG10/11/12; same 48 source-quality US
states and same true directed origin→destination pairs. Control BOTH
cargo TRUCK modes, other cargo RAIL OD, state freight hub exposure,
Census geography and original botanical recording scale.

NEVER call source first-record an actual biological introduction, state-
alien status, railroad seed transfer, mowing/wind filter or adaptation.
"""
from __future__ import annotations
import argparse,json,math
from collections import defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests,numpy as np

from analyze_US_FAF6_directional_source_taxa_v86 import (
    construct_directions,RESPONSES,COHORTS,PRIMARY_CONTROLS,
    partial_rank as original_partial_rank)
from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    CENSUS_TIGER_URL,parse_census_centroids,_rank)
from analyze_FAF6_crop_vs_resource_rail_floristic_v84 import (
    COMMODITIES,extract_OD,verify_matrix_marginals,assert_original_FAF_integrity)
from analyze_FAF6_crop_OD_rail_vs_truck_network_v80 import build_networks
from audit_BTS_FAF6_2022_agricultural_rail_freight_v79 import (
    URL,records_from_faf6_zip)

OUTCOMES=("expected_20record_source_linked_first_record_fraction",
          "raw_source_linked_first_record_fraction")
PREDICTORS=("crop_rail","resource_rail")
SOURCE_COHORTS=("all_qualified_directed_pairs",
                "crop_RAIL_positive","crop_RAIL_and_TRUCK_positive",
                "both_crop_and_resource_RAIL_positive")
MIN_STATES=20
MIN_PAIR_SUPPORT=30
JOINT_CONTROLS=tuple(PRIMARY_CONTROLS)+(
    "log1p_resource_truck_OD",
    "log1p_resource_origin_rail_out_strength",
    "log1p_resource_destination_rail_in_strength",
    "log1p_resource_origin_truck_out_strength",
    "log1p_resource_destination_truck_in_strength",
)


def source_matrices(raw,prior79,prior80):
    flows,states,trace=build_networks(records_from_faf6_zip(raw))
    assert_original_FAF_integrity(raw,prior79,prior80,trace["original_source_rows"],len(states))
    matrices={g:{m:extract_OD(flows,g,m) for m in ("RAIL","TRUCK")}
              for g in COMMODITIES}
    verify_matrix_marginals(matrices,prior79,prior80)
    return matrices,trace


def source_in_out(matrix):
    out=defaultdict(float);destination=defaultdict(float)
    for (o,d),tons in matrix.items():
        if o!=d:
            out[o]+=tons;destination[d]+=tons
    return out,destination


def augment_directed_rows(rows,matrix):
    crop=matrix["crop_sctg02_03"]
    resource=matrix["resource_sctg10_12"]
    strength={m:source_in_out(resource[m]) for m in ("RAIL","TRUCK")}
    built=[]
    for row in rows:
        o,d=row["o"],row["d"]
        crop_r=crop["RAIL"].get((o,d),0.0)
        crop_t=crop["TRUCK"].get((o,d),0.0)
        if (abs(crop_r-row["rail_2022_estimated_thousand_short_tons"])>1e-5 or
            abs(crop_t-row["truck_2022_estimated_thousand_short_tons"])>1e-5):
            raise ValueError("Original v8.6 actual direction source crop RAIL/TRUCK drift")
        rr=resource["RAIL"].get((o,d),0.0)
        rt=resource["TRUCK"].get((o,d),0.0)
        cov=dict(row["covars"])
        cov.update({
            "log1p_resource_truck_OD":math.log1p(rt),
            "log1p_resource_origin_rail_out_strength":math.log1p(strength["RAIL"][0][o]),
            "log1p_resource_destination_rail_in_strength":math.log1p(strength["RAIL"][1][d]),
            "log1p_resource_origin_truck_out_strength":math.log1p(strength["TRUCK"][0][o]),
            "log1p_resource_destination_truck_in_strength":math.log1p(strength["TRUCK"][1][d]),
        })
        built.append({**row,"crop_rail":crop_r,"crop_truck":crop_t,
                      "resource_rail":rr,"resource_truck":rt,
                      "covars":cov})
    return built


def source_cohort(rows,key):
    if key==SOURCE_COHORTS[0]:return rows
    if key==SOURCE_COHORTS[1]:
        return [r for r in rows if r["crop_rail"]>0]
    if key==SOURCE_COHORTS[2]:
        return [r for r in rows if min(r["crop_rail"],r["crop_truck"])>0]
    if key==SOURCE_COHORTS[3]:
        return [r for r in rows if min(r["crop_rail"],r["resource_rail"])>0]
    raise ValueError("Cargo cohort must be selected independent of first-recorded flora")


def conditional_partial_rank(rows,response,predictor,controls=JOINT_CONTROLS):
    if response not in OUTCOMES or predictor not in PREDICTORS:
        raise ValueError("Unpredeclared botanical source response or commodity RAIL exposure")
    if len(rows)<MIN_PAIR_SUPPORT or min(
        len({r["o"] for r in rows}),len({r["d"] for r in rows}))<MIN_STATES:
        return None
    alternate="resource_rail" if predictor=="crop_rail" else "crop_rail"
    x=np.asarray([math.log1p(r[predictor]) for r in rows],float)
    y=np.asarray([r["responses"][response] for r in rows],float)
    z=np.asarray([[r["covars"][name] for name in controls]+[math.log1p(r[alternate])]
                  for r in rows],float)
    if min(len(np.unique(x)),len(np.unique(y)))<5 or not all(
        np.all(np.isfinite(a)) for a in (x,y,z)):
        return None
    A=np.column_stack((np.ones(len(rows)),*(_rank(z[:,i]) for i in range(z.shape[1]))))
    if np.linalg.matrix_rank(A)<A.shape[1] or len(rows)-A.shape[1]<10:
        return None
    xr=_rank(x)-A@np.linalg.lstsq(A,_rank(x),rcond=None)[0]
    yr=_rank(y)-A@np.linalg.lstsq(A,_rank(y),rcond=None)[0]
    n=float(np.linalg.norm(xr)*np.linalg.norm(yr))
    if n<=1e-10:return None
    rho=float(xr@yr/n)
    if not -1.000001<=rho<=1.000001:raise ValueError("Source partial rho invalid")
    return round(rho,6)


def cohort_summary(rows):
    out={
        "n_directed_state_pairs_NOT_independent":len(rows),
        "n_source_origin_states":len({r["o"] for r in rows}),
        "n_source_destination_states":len({r["d"] for r in rows}),
        "source_quality_gate_at_least_20_origin_and_dest_states":False,
        "botanical_record_source_results":{},
    }
    supported=(len(rows)>=MIN_PAIR_SUPPORT and
               min(out["n_source_origin_states"],out["n_source_destination_states"])>=MIN_STATES)
    out["source_quality_gate_at_least_20_origin_and_dest_states"]=supported
    for response in OUTCOMES:
        x={}
        for p in PREDICTORS:
            x[p]=conditional_partial_rank(rows,response,p) if supported else None
        out["botanical_record_source_results"][response]={
            "crop_RAIL_full_adjusted_rho":x["crop_rail"],
            "resource_RAIL_full_adjusted_rho":x["resource_rail"],
            "crop_minus_resource_rho":round(x["crop_rail"]-x["resource_rail"],6)
              if all(v is not None for v in x.values()) else None}
    return out


def node_deletion(rows,predictor):
    estimates={}
    for var in ("o","d"):
        unique=sorted({r[var] for r in rows})
        scores=[conditional_partial_rank([r for r in rows if r[var]!=code],
                                         OUTCOMES[0],predictor) for code in unique]
        valid=[x for x in scores if x is not None]
        estimates["origin" if var=="o" else "destination"]={
            "n_deleted_source_state_nodes":len(unique),
            "n_estimable":len(valid),
            "n_positive":sum(z>0 for z in valid),
            "n_negative":sum(z<0 for z in valid),
            "min_rho":min(valid) if valid else None,
            "max_rho":max(valid) if valid else None,
            "NOT_independent_country_replication_or_valid_inference":True,
        }
    return estimates


def run(temporal,flora82,prior86,fa79,fa80,census_bytes,faf_raw):
    geo=parse_census_centroids(census_bytes.decode("utf-8"),flora82["states"])
    if prior86.get("complete") is not True:
        raise ValueError("Original temporal rail crop source v8.6 missing")
    pairs,sources,missing=construct_directions(temporal,fa80,geo)
    if len(pairs)!=2255 or len(sources)!=48 or (
          prior86.get("n_evaluable_directed_origin_destination_source_pairs")!=2255):
        raise ValueError("Previously verified 48-state source geography/botanical record frame drift")
    for cohort in COHORTS:
        src=([r for r in pairs] if cohort==COHORTS[0] else
             [r for r in pairs if r["rail_2022_estimated_thousand_short_tons"]>0
              and (cohort!=COHORTS[2] or r["truck_2022_estimated_thousand_short_tons"]>0)])
        for response in OUTCOMES:
            val=original_partial_rank(src,response) if len(src)>=30 else None
            original=prior86["all_results"][cohort]["results"][response]["rail_full"]
            if val!=original:
                raise ValueError("Verified directional temporal v8.6 cargo flora response did not reproduce")
    matrix,trace=source_matrices(faf_raw,fa79,fa80)
    joined=augment_directed_rows(pairs,matrix)
    result={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "source_US_48_quality_states_pre_2018_21_post_2023_25":True,
        "source_2022_FAF6_commodity_OD_sha_and_marginals_reconciled":True,
        "source_prior_directional_v86_partial_rho_reconciled":True,
        "source_2022_faf_rows":trace["original_source_rows"],
        "analysis_by_predefined_cargo_positive_cohort":{},
        "one_origin_and_destination_source_state_deleted":{},
        "2022_directional_cargo_overlap":{},
        "NO_observed_transport_propagules_or_true_naturalization_or_mowing_wind_adaptation":True,
    }
    for name in SOURCE_COHORTS:
        sub=source_cohort(joined,name)
        result["analysis_by_predefined_cargo_positive_cohort"][name]=cohort_summary(sub)
    for predictor in PREDICTORS:
        result["one_origin_and_destination_source_state_deleted"][predictor]=node_deletion(
            joined,predictor)
    joint=[r for r in joined if min(r["crop_rail"],r["resource_rail"])>0]
    result["2022_directional_cargo_overlap"]={
        "n_crop_and_resource_RAIL_both_positive_in_48_state_directional_cohort":len(joint),
        "n_eligible_source_origins":len({r["o"] for r in joint}),
        "n_eligible_source_destinations":len({r["d"] for r in joint}),
        "NOT_actually_seed_free_resource_material":True}
    result["limitations"]=[
        "Inferred 2023–25 newly recorded source speciesKey at a destination is NOT a documented plant immigration, local introduced-status change or known biological absence in 2018–21.",
        "Even a positive association of crop cargo with post-2022 indexed first records would not establish actual seed transport, railway route use or germination.",
        "SCTG10–12 resource cargo is not a proven viable-seed-free control. 2022 commodity state OD is estimated economic freight, not line-specific train traffic.",
        "All state directions repeatedly share origin/destination recording sources; no dyadic or leave-state-out pseudoindependent p-values are meaningful.",
        "Post-2022 three-year vs prior four-year stateProvince source lists have differing chance of recorded GBIF taxa; indexed record rarefaction does not correct provider/observer bias.",
        "No railway mowing, train gusts, ballast stress, morphological plasticity or heritable adaptive evolution is measured.",
        "Single country and freight base year, no country/time replication or counterfactual colonization history."
    ]
    return result


def main():
    ap=argparse.ArgumentParser()
    for x in ("temporal85","flora82","prior86","fa79","fa80","output"):
        ap.add_argument("--"+x,required=True)
    a=ap.parse_args()
    read=lambda field:json.loads(Path(getattr(a,field)).read_text(encoding="utf-8"))
    census=requests.get(CENSUS_TIGER_URL,timeout=(20,90))
    census.raise_for_status()
    faf=requests.get(URL,timeout=(30,180))
    faf.raise_for_status()
    data=run(read("temporal85"),read("flora82"),read("prior86"),
             read("fa79"),read("fa80"),census.content,faf.content)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "cohorts":data["analysis_by_predefined_cargo_positive_cohort"],
        "cargo_OD_overlap":data["2022_directional_cargo_overlap"],
        "state_leaveout":data["one_origin_and_destination_source_state_deleted"],
        "NO_causal_seed_transport_or_rail_mowing_wind_adaptation":True},
        ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

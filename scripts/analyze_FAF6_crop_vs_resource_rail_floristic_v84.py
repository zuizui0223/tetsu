#!/usr/bin/env python3
"""v8.4 observational macro source falsification: do rail agricultural
freight OD links predict recorded US-country-GRIIS plant similarity more
specifically than SAME YEAR non-crop resource freight OD?

BTS FAF6.0 2022 STATE estimated economic OD × Census 2025 51-state
geography × GBIF 2018-25 country-GRIIS recorded plant speciesKeys.
State-pairs reuse states; no pairwise independent biological p-values.

SCTG10-12 are DIFFERENT-COMMODITY comparison, *not* guaranteed sterile
or true seed-free controls. No observed viable seeds, operated rail
tracks, destination colonization, cutting, gusts or adaptation.
"""
from __future__ import annotations
import argparse,json,math
from collections import defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests
import numpy as np

from analyze_macro_50km_rail_road_aliens_v34 import _rank
from analyze_FAF6_crop_OD_rail_vs_truck_network_v80 import (
    build_networks,GROUPS,BASELINE_GROUPS,MODES)
from audit_BTS_FAF6_2022_agricultural_rail_freight_v79 import (
    records_from_faf6_zip,URL)
from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    CENSUS_TIGER_URL,PRIMARY_CONTROLS,METRICS,
    validated_botanical_states,parse_census_centroids,state_pair_rows,
    partial_rank_rail,node_strength_nonself)

COMMODITIES=("crop_sctg02_03","resource_sctg10_12")
PRIMARY_COHORTS=(
    "all_49_state_pairs",
    "crop_rail_positive",
    "crop_rail_and_truck_positive",
    "both_crop_and_resource_rail_positive",
)
MIN_STATE_PAIRS=30
MIN_DISTINCT_STATE_SOURCES=20
CONTROL_GEOGRAPHY=tuple(k for k in PRIMARY_CONTROLS
                        if k not in ("log1p_truck_OD",
                         "log1p_rail_node_strength_geomean",
                         "log1p_truck_node_strength_geomean"))
JOINT_CONTROLS=CONTROL_GEOGRAPHY+(
    "log1p_crop_truck_OD","log1p_resource_truck_OD",
    "log1p_crop_rail_strength_geomean","log1p_crop_truck_strength_geomean",
    "log1p_resource_rail_strength_geomean",
    "log1p_resource_truck_strength_geomean",
)
SINGLE_CONTROLS=CONTROL_GEOGRAPHY+(
    "log1p_rail_node_strength_geomean",
    "log1p_truck_node_strength_geomean","log1p_truck_OD")
LARGE_NO_INDEPENDENT_DYADIC_SIGNIFICANCE=True


def assert_original_FAF_integrity(raw,fa79,fa80,rows,nstates):
    expected=fa79.get("official_published_CSV_ZIP_sha256")
    if (fa79.get("complete") is not True or
            fa80.get("complete") is not True or
            fa80.get("source_hash_verified_vs_prior_v79") is not True or
            fa80.get("source_trace",{}).get("previous_verified_v79_total_and_modes_reconciled") is not True or
            nstates!=51 or
            sha256(raw).hexdigest()!=expected or
            rows!=fa79["data"]["n_original_official_FAF6_state_OD_commodity_mode_rows"] or
            rows!=fa80["source_trace"]["original_source_rows"]):
        raise ValueError("FAF6 official CSV SHA/2022 rows/states drift since verified v7.9/v8.0")
    return True


def extract_OD(flows,group,mode):
    if group not in COMMODITIES or mode not in ("RAIL","TRUCK"):
        raise ValueError("Only predeclared crop or resource freight comparisons allowed")
    ans={}
    for (g,m,o,d),w in flows.items():
        if g==group and m==mode:
            if (o,d) in ans or not math.isfinite(w) or w<=0:
                raise ValueError("FAF6 source economic OD duplicate or nonpositive")
            ans[o,d]=w
    return ans


def symmetric_OD(matrix,a,b):
    return matrix.get((a,b),0)+matrix.get((b,a),0)


def verify_matrix_marginals(by_group,fa79,fa80):
    for group in COMMODITIES:
        old=fa79["data"]["national_trade_type_x_commodity_group_2022"]["domestic"][BASELINE_GROUPS[group]]
        for mode in ("RAIL","TRUCK"):
            total=sum(by_group[group][mode].values())
            key="rail_thousand_short_tons" if mode=="RAIL" else "truck_thousand_short_tons"
            orig=fa80["mode_networks"][group][mode]["total_2022_estimated_thousand_short_tons"]
            if abs(total-old[key])>0.00001 or abs(total-orig)>0.00001:
                raise ValueError("2022 crop/resource actual OD matrix marginals disagree with original independent source")
    original=fa80["raw_official_state_OD_2022_crop_estimated_tonnage"]
    for mode in ("RAIL","TRUCK"):
        for row in original[mode]:
            o,d=row["o"],row["d"]
            if abs(by_group["crop_sctg02_03"][mode].get((o,d),0)-row["estimated_thousand_short_tons"])>1e-5:
                raise ValueError("v8.0 original crop FAF6 pair weights not reproduced")
    return True


def augment_base_pairs(source_pairs,by_group):
    strength={}
    for group in COMMODITIES:
        strength[group]={m:node_strength_nonself(by_group[group][m]) for m in ("RAIL","TRUCK")}
    output=[]
    for row in source_pairs:
        a,b=row["a"],row["b"]
        crop=symmetric_OD(by_group["crop_sctg02_03"]["RAIL"],a,b)
        crop_t=symmetric_OD(by_group["crop_sctg02_03"]["TRUCK"],a,b)
        if (abs(crop-row["rail_estimated_thousand_short_tons"])>1e-5 or
            abs(crop_t-row["truck_estimated_thousand_short_tons"])>1e-5):
            raise ValueError("Real crop source OD weights changed against v8.3 prior matched geography")
        rr=symmetric_OD(by_group["resource_sctg10_12"]["RAIL"],a,b)
        rt=symmetric_OD(by_group["resource_sctg10_12"]["TRUCK"],a,b)
        cov=dict(row["covars"])
        cov.update({
          "log1p_crop_truck_OD":math.log1p(crop_t),
          "log1p_resource_truck_OD":math.log1p(rt),
        })
        for group,short in (("crop_sctg02_03","crop"),("resource_sctg10_12","resource")):
            for mode in ("RAIL","TRUCK"):
                key=f"log1p_{short}_{mode.lower()}_strength_geomean"
                va=strength[group][mode]
                cov[key]=(math.log1p(va[a])+math.log1p(va[b]))/2
        output.append({
            "a":a,"b":b,
            "crop_rail":crop,
            "crop_truck":crop_t,
            "resource_rail":rr,
            "resource_truck":rt,
            "responses":row["responses"],
            "covars":cov,
        })
    if len(output)!=len(source_pairs):
        raise ValueError("FAF6 source comparison omitted original 49-state pairs")
    return output


def partial_rank_pair(rows,response,predictor,controls):
    if response not in METRICS or predictor not in ("crop_rail","resource_rail"):
        raise ValueError("Unprespecified flora response or railway commodity predictor")
    if len(rows)<MIN_STATE_PAIRS or len({s for z in rows for s in (z["a"],z["b"])})<MIN_DISTINCT_STATE_SOURCES:
        return None
    y=np.asarray([z["responses"][response] for z in rows],float)
    x=np.asarray([math.log1p(z[predictor]) for z in rows],float)
    control_series=[("resource_rail" if predictor=="crop_rail" else "crop_rail")]+list(controls)
    Z=np.asarray([[math.log1p(z[k]) if k in ("crop_rail","resource_rail") else z["covars"][k]
                   for k in control_series] for z in rows],float)
    if min(len(np.unique(x)),len(np.unique(y)))<5 or not all(np.isfinite(v).all() for v in (x,y,Z)):
        return None
    A=np.column_stack([np.ones(len(rows)),*[_rank(Z[:,j]) for j in range(Z.shape[1])]])
    if np.linalg.matrix_rank(A)<A.shape[1] or len(rows)-A.shape[1]<10:
        return None
    xr=_rank(x)-A@np.linalg.lstsq(A,_rank(x),rcond=None)[0]
    yr=_rank(y)-A@np.linalg.lstsq(A,_rank(y),rcond=None)[0]
    den=float(np.linalg.norm(xr)*np.linalg.norm(yr))
    if den<1e-8:return None
    rho=float(xr@yr/den)
    if not -1.0000001<=rho<=1.0000001:
        raise ValueError("Partial rank source outside [-1,1]")
    return round(rho,6)


def one_cohort(rows):
    states={s for row in rows for s in (row["a"],row["b"])}
    out={
        "n_source_US_state_pairs_NOT_INDEPENDENT":len(rows),
        "n_source_states_in_pair_cohort":len(states),
        "sufficient_independent_state_source_opportunity_not_independent_ecological_tests":
            len(states)>=MIN_DISTINCT_STATE_SOURCES and len(rows)>=MIN_STATE_PAIRS,
        "outcome_source":{},
    }
    for metric in METRICS:
        out["outcome_source"][metric]={}
        for cargo,control in (("crop_rail","resource_rail"),("resource_rail","crop_rail")):
            rho=partial_rank_pair(rows,metric,cargo,JOINT_CONTROLS)
            out["outcome_source"][metric][cargo]={
                "joint_adjusted_partial_rank_rho":rho,
                "other_commodity_rail_OD_adjusted":control,
                "same_state_pair_and_botanical_recording_geography":True,
                "hypothesis_direction_positive_more_rail_with_more_recorded_similarity":True,
            }
        a=out["outcome_source"][metric]
        x,y=a["crop_rail"]["joint_adjusted_partial_rank_rho"],a["resource_rail"]["joint_adjusted_partial_rank_rho"]
        out["outcome_source"][metric]["descriptive_crop_minus_resource_partial_rho"]=(round(x-y,6)
                                                                                    if x is not None and y is not None else None)
    return out


def same_pair_cohort(rows,label):
    if label=="all_49_state_pairs":
        return list(rows)
    if label=="crop_rail_positive":
        return [r for r in rows if r["crop_rail"]>0]
    if label=="crop_rail_and_truck_positive":
        return [r for r in rows if r["crop_rail"]>0 and r["crop_truck"]>0]
    if label=="both_crop_and_resource_rail_positive":
        return [r for r in rows if r["crop_rail"]>0 and r["resource_rail"]>0]
    raise ValueError("Unfrozen freight route cohort selection")


def refit_without_one_state(rows,metric,predictor):
    full=partial_rank_pair(rows,metric,predictor,JOINT_CONTROLS)
    estimates={}
    for state in sorted({s for r in rows for s in (r["a"],r["b"])}):
        remaining=[r for r in rows if state not in (r["a"],r["b"])]
        estimates[state]=partial_rank_pair(remaining,metric,predictor,JOINT_CONTROLS)
    good=[v for v in estimates.values() if v is not None]
    return {
        "full_rho":full,
        "n_withheld_states_not_independent_countries":len(estimates),
        "n_estimable_refits":len(good),
        "min_rho":min(good) if good else None,
        "max_rho":max(good) if good else None,
        "n_positive":sum(v>0 for v in good),
        "n_negative":sum(v<0 for v in good),
        "NOT_bootstrap_or_a_valid_statepair_p_value":True
    }


def commodity_specificity_stress_test(rows):
    """Posthoc sensitivity: same-OD crop/resource contrast under state
    deletion, independently defined cargo-tonnage exclusions and all
    predeclared covariate definitions. Do not select favorable cases."""
    metric="rarefied20_expected_inclusion_overlap"
    controls={
        "full_joint_adjustment":JOINT_CONTROLS,
        "omit_all_node_strength_proxies":CONTROL_GEOGRAPHY+(
            "log1p_crop_truck_OD","log1p_resource_truck_OD"),
        "distance_and_recording_effort_plus_both_truck_modes":(
            "log_state_centroid_geo_distance",
            "log_allplant_record_geomean","log_allplant_record_ratio",
            "log1p_crop_truck_OD","log1p_resource_truck_OD"),
    }
    out={
        "n_source_pairs_not_independent":len(rows),
        "n_source_states":len({s for r in rows for s in (r["a"],r["b"])}),
        "source_outcome":"2018-25 country-GRIIS-indexed GBIF K20",
        "all_predeclared_alternatives_reported_even_if_hypothesis_incompatible":True,
        "covariate_sensitivity":{},
        "one_state_deleted":{},
        "largest_economic_cargo_OD_pairs_removed":{},
    }
    for name,ctl in controls.items():
        rhos={p:partial_rank_pair(rows,metric,p,ctl)
              for p in ("crop_rail","resource_rail")}
        a,b=rhos["crop_rail"],rhos["resource_rail"]
        out["covariate_sensitivity"][name]={
            "source_partial_rhos":rhos,
            "source_crop_minus_resource_delta":round(a-b,6)
                if a is not None and b is not None else None,
            "n_rank_covariates_plus_other_rail_cargo":len(ctl)+1}
    for p in ("crop_rail","resource_rail"):
        out["one_state_deleted"][p]=refit_without_one_state(rows,metric,p)
    for n in (1,3,5):
        ordered=sorted(rows,key=lambda r:(-(r["crop_rail"]+r["resource_rail"]),r["a"],r["b"]))
        subset=ordered[n:]
        rhos={p:partial_rank_pair(subset,metric,p,JOINT_CONTROLS)
              for p in ("crop_rail","resource_rail")}
        a,b=rhos["crop_rail"],rhos["resource_rail"]
        out["largest_economic_cargo_OD_pairs_removed"][f"remove_top_{n}"]={
            "remaining_pairs":len(subset),
            "remaining_unique_states":len({s for r in subset for s in (r["a"],r["b"])}),
            "source_partial_rhos":rhos,
            "source_crop_minus_resource_delta":round(a-b,6)
                if a is not None and b is not None else None}
    out["NOT_evidence_for_seed_contamination_or_mowing_adaptation"]=True
    return out


def run(previous_fa79,previous_fa80,botanical,previous_v83,
        census_bytes,raw_faf_bytes):
    if (previous_v83.get("complete") is not True or
        previous_v83.get("n_eligible_country_GRIIS_indexed_US_state_sources")!=49 or
        previous_v83.get("Census_official_2025_TIGER_state_table_sha256")!=sha256(census_bytes).hexdigest()):
        raise ValueError("Census 2025 geographic source not the same independently verified v8.3 content")
    flows,states,trace=build_networks(records_from_faf6_zip(raw_faf_bytes))
    assert_original_FAF_integrity(raw_faf_bytes,previous_fa79,previous_fa80,
                                  trace["original_source_rows"],len(states))
    by_group={group:{mode:extract_OD(flows,group,mode) for mode in ("RAIL","TRUCK")}
              for group in COMMODITIES}
    verify_matrix_marginals(by_group,previous_fa79,previous_fa80)
    biol,excluded=validated_botanical_states(botanical)
    geo=parse_census_centroids(census_bytes.decode("utf-8"),botanical["states"])
    source_pairs=state_pair_rows(biol,geo,previous_fa80)
    if len(source_pairs)!=1176:
        raise ValueError("Original 49-state × 49-state cohort not preserved")
    for metric in METRICS:
        reference=previous_v83["analysis"]["all_eligible_state_pairs"]["responses"][metric]["rail_partial_rho"]
        actual=partial_rank_rail(source_pairs,metric)
        if reference!=actual:
            raise ValueError("Original v8.3 same-geography crop source coefficient did not exactly reproduce")
    joined=augment_base_pairs(source_pairs,by_group)
    out={
        "complete":True,
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "study":"v8.4 crops SCTG02/03 versus resources SCTG10/11/12 estimated RAIL freight associations to recorded country GRIIS USA flora",
        "genuine_original_2022_BTS_FAF6_zip_SHA256_reconciled":True,
        "genuine_original_2025_Census_geographic_table_SHA256_reconciled":True,
        "reproduced_independent_v83_crop_allpairs_source_coefficients":True,
        "same_preexisting_49_US_GBIF_quality_states":True,
        "cohorts":{},
        "leave_one_state_source_sensitivity":{},
        "source_flow_network_overlap":{},
        "conditional_crop_specificity_robustness":{},
        "NOT_a_negative_control_of_viable_seed_contamination":True,
        "NOT_observed_cargo_seed_flux_or_mowing_wind_evolution":True,
    }
    crop_set={p for p in by_group["crop_sctg02_03"]["RAIL"] if p[0]!=p[1]}
    resource_set={p for p in by_group["resource_sctg10_12"]["RAIL"] if p[0]!=p[1]}
    out["source_flow_network_overlap"]={
        "n_positive_directed_interstate_crop_RAIL_OD":len(crop_set),
        "n_positive_directed_interstate_resource_RAIL_OD":len(resource_set),
        "n_positive_directed_interstate_OD_both_cargo_groups":len(crop_set & resource_set),
        "n_positive_crop_RAIL_OD_not_resource_RAIL":len(crop_set-resource_set),
        "n_positive_resource_RAIL_OD_not_crop_RAIL":len(resource_set-crop_set),
        "NOT_verified_rail_segments_or_viable_propagules":True,
    }
    for cohort in PRIMARY_COHORTS:
        sub=same_pair_cohort(joined,cohort)
        out["cohorts"][cohort]=one_cohort(sub)
        if cohort in ("crop_rail_and_truck_positive","both_crop_and_resource_rail_positive"):
            out["conditional_crop_specificity_robustness"][cohort]=commodity_specificity_stress_test(sub)
    for predictor in ("crop_rail","resource_rail"):
        out["leave_one_state_source_sensitivity"][predictor]=refit_without_one_state(
            joined,"rarefied20_expected_inclusion_overlap",predictor)
    out["source_population"]={
        "original_faf2022_source_rows":trace["original_source_rows"],
        "US_Census_GBIF_matched_eligible_state_count":len(biol),
        "n_geographical_state_pair_rows_sharing_states":len(joined),
        "source_noneligible_GBIF_states":excluded,
        "source_freight_commodity_groups":{k:sorted(GROUPS[k]) for k in COMMODITIES}
    }
    out["interpretation_limits"]=[
        "Cross-sectional 2022 estimated economic cargo tonnage and 2018–25 GBIF source-recorded plants cannot determine train-mediated species dispersal or post-arrival colonization.",
        "SCTG10-12 resources are a different-commodity COMPARATOR, not sterile freight, randomized ecological negative control or verified absence of seeds.",
        "Rail and truck FAF region OD is NOT proven physical train route, port connection or exact freight track load.",
        "GBIF taxon country GRIIS listing is not locally checked alien status; stateProvince is text geography not validated occurrence-by-state geometry.",
        "Source rarefaction is based on GBIF indexed records and cannot eliminate provider/observer, climatic, agricultural-production and urban biases.",
        "All 1176 geographic state-pair outcomes reuse 49 states and share flora; no dyadic pseudo-independent p-values.",
        "Source transport comparators and train track stresses are conceptually different: none of mowing, gusts, thermal filtering, phenotype plasticity or genetic adaptation was measured.",
        "Single source USA/2022 does not establish an independent cross-national macroecological law."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for s in ("fa79","fa80","flora82","source83","output"):
        p.add_argument("--"+s,required=True)
    a=p.parse_args();load=lambda x:json.loads(Path(getattr(a,x)).read_text(encoding="utf-8"))
    geo=requests.get(CENSUS_TIGER_URL,timeout=(20,90),
        headers={"User-Agent":"tetsu-v84-same-census-freight-traits-control"})
    geo.raise_for_status()
    freight=requests.get(URL,timeout=(30,180),
        headers={"User-Agent":"tetsu-v84-same-FAF6-2022-source"})
    freight.raise_for_status()
    out=run(load("fa79"),load("fa80"),load("flora82"),load("source83"),
            geo.content,freight.content)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
      "complete":True,
      "source":out["source_population"],
      "cohorts":out["cohorts"],
      "leave_one_state_source_sensitivity":out["leave_one_state_source_sensitivity"],
      "economic_OD_overlap":out["source_flow_network_overlap"],
      "conditional_crop_specificity_robustness":out["conditional_crop_specificity_robustness"],
      "NO_biological_or_train_seed_mowing_wind_adaptation_claim":True,
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()

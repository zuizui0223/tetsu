#!/usr/bin/env python3
"""v8.9: origin/destination STATE fixed-effects plus constrained node-label
QAP-style SOURCE-NETWORK PERMUTATION for v8.6 directional 2022 FAF6
RAIL/TRUCK crop freight vs 2018-21 origin / 2023-25 destination
GRIIS-country-listed recorded plant taxon candidates.

No dyadic independent p-values. Node permutations, even constrained
by 2025 Census state REGION/DIVISION, do NOT guarantee exchangeability
or calibrated randomization inference when agriculture, rail routing,
environment and botanical reporting are historically structured.

Permutation exceedance is a source-network structural diagnostic,
NOT a causal ecological p-value or verified seed movement.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests
import numpy as np

from analyze_US_FAF6_directional_source_taxa_v86 import (
    COHORTS,RESPONSES,construct_directions,partial_rank as legacy_partial_rank)
from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    CENSUS_TIGER_URL,parse_census_centroids,original_crop_OD,_rank)

SEED=20261112
PERMUTATIONS={"Census_region_preserved":499,
              "Census_division_preserved":299,
              "unrestricted_state_codes_descriptive_only":199}
MODELS=("RAIL","TRUCK")
RESPONSE=RESPONSES[0]
MIN_STATES=40
MIN_ROWS=1000
PRIMARY_PAIR_CONTROLS=(
    "log_geographical_distance",
    "same_Census_region",
    "same_Census_division",
    "log_origin_candidate_taxon_categories",
)


def prepare_dyadic(rows,geometry,other="TRUCK"):
    """Two-way origin/destination state effects absorb state effort,
    freight origin/destination hub strength, land area, broad survey
    propensity. Pairwise truck/rail source mode must still be controlled.
    """
    ids=sorted({z["o"] for z in rows}|{z["d"] for z in rows})
    if len(ids)<MIN_STATES or len(rows)<MIN_ROWS:
        raise ValueError("Insufficient original 48-state directional source support")
    if any(z not in geometry for z in ids):
        raise ValueError("Missing original Census FIPS source geography")
    cat={key:i for i,key in enumerate(ids)}
    n=len(rows)
    # With first state's origin+destination coefficients zero as reference,
    # all other fixed effects are estimated on original 2255 direction rows.
    a=np.zeros((n,1+2*(len(ids)-1)+1+len(PRIMARY_PAIR_CONTROLS)),float)
    a[:,0]=1
    for i,r in enumerate(rows):
        o,d=cat[r["o"]],cat[r["d"]]
        if o:a[i,o]=1
        if d:a[i,(len(ids)-1)+d]=1
        j=1+2*(len(ids)-1)
        a[i,j]=r["x_truck" if other=="TRUCK" else "x_rail"]
        a[i,j+1:]=[r["covars"][s] for s in PRIMARY_PAIR_CONTROLS]
    # Check directly before rank transformation: full rank does not mean
    # origin/dest 48 observational state units are independent countries.
    if np.linalg.matrix_rank(a)<a.shape[1]:
        raise ValueError("Origin/destination fixed-effects control frame rank-deficient")
    a[:,1:]=np.column_stack([_rank(a[:,j]) for j in range(1,a.shape[1])])
    # Binary dummies rank to affine transforms; rank under complete
    # intercept remains full-rank in this separately matched source frame.
    q,r=np.linalg.qr(a,mode="reduced")
    if np.min(np.abs(np.diag(r)))<1e-8:
        raise ValueError("Two-way fixed-effects residualization rank deficient")
    return {
        "states":ids,
        "src_rows":rows,
        "Q":q,
        "Y":{metric:residual(_rank(np.array([r["responses"][metric] for r in rows])),q)
             for metric in RESPONSES},
        "n_observed_direction_pairs":n,
        "control_count":a.shape[1]-1,
        "Census_region_partition":{s:geometry[s]["region"] for s in ids},
        "Census_division_partition":{s:geometry[s]["division"] for s in ids}
    }


def residual(z,q):
    return z-q@(q.T@z)


def residual_rho(x,y,q):
    xr=residual(_rank(np.asarray(x,float)),q)
    den=float(np.linalg.norm(xr)*np.linalg.norm(y))
    if den<1e-10:return None
    result=float(xr@y/den)
    if not -1.0000001<=result<=1.0000001:
        raise ValueError("Invalid source FE partial rank correlation")
    return round(result,6)


def permutation_map(ids,groups,rng):
    """Simultaneously permute origins AND destinations, preserving
    the directed FAF6 network's internal topology and edge weights."""
    b=defaultdict(list)
    for state in ids:
        b[groups[state]].append(state)
    perm={}
    for val in sorted(b):
        members=sorted(b[val])
        target=members[:]
        rng.shuffle(target)
        perm.update(zip(members,target))
    if set(perm)!=set(ids) or len(set(perm.values()))!=len(ids):
        raise ValueError("Original source 48-state bijective permutation failed")
    return perm


def predictor_vector(rows,od,mapping):
    values=[math.log1p(od.get((mapping[z["o"]],mapping[z["d"]]),0.0))
            for z in rows]
    if any(not math.isfinite(x) or x<0 for x in values):
        raise ValueError("Invalid permuted real FAF state freight source OD")
    return values


def one_mode(mode,source,odmatrix):
    if mode not in MODELS:
        raise ValueError("Only original verified crop RAIL and TRUCK modes")
    rows=source["src_rows"];q=source["Q"];ids=source["states"]
    identity={code:code for code in ids}
    # Exactly same original paired source biological response, 48 states
    # and original FAF 2022 directional crop matrix for each mode.
    x=predictor_vector(rows,odmatrix,identity)
    observed={response:residual_rho(x,y,q) for response,y in source["Y"].items()}
    if observed[RESPONSE] is None:
        raise ValueError("No supported original QAP-style fixed-effect source coefficient")
    out={
        "mode":mode,
        "n_same_directed_source_pairs_not_independent":len(rows),
        "n_same_GBIF_state_sources":len(ids),
        "observed_origin_destination_FE_source_partial_rho":observed,
        "constrained_FA6_directed_network_node_label_null":{},
        "actual_observed_positive_economic_OD_links_in_same_state_frame":
            sum(odmatrix.get((r["o"],r["d"]),0)>0 for r in rows),
        "no_calibrated_causal_p_value_from_network_permutations":True,
    }
    for idx,(name,nperm) in enumerate(PERMUTATIONS.items()):
        if name=="Census_region_preserved":
            groups=source["Census_region_partition"]
        elif name=="Census_division_preserved":
            groups=source["Census_division_partition"]
        else:
            groups={s:"US" for s in ids}
        rng=np.random.default_rng(SEED+idx*100+MODELS.index(mode)*10000)
        null=[]
        same=0
        for _ in range(nperm):
            mapping=permutation_map(ids,groups,rng)
            if all(mapping[s]==s for s in ids):
                same+=1
                continue
            xx=predictor_vector(rows,odmatrix,mapping)
            rr=residual_rho(xx,source["Y"][RESPONSE],q)
            if rr is not None:null.append(rr)
        if len(null)<max(150,int(.8*nperm)):
            raise ValueError("Too many nonestimable constrained state-label permutation trials")
        arr=np.asarray(null)
        actual=observed[RESPONSE]
        out["constrained_FA6_directed_network_node_label_null"][name]={
            "original_node_label_network_topology_conserved":True,
            "n_requested_label_swaps":nperm,
            "n_same_identity_permutations_skipped":same,
            "n_valid_node_relabeling_trials":len(null),
            "null_median_rho":round(float(np.median(arr)),6),
            "null_025_and_975_quantiles":list(np.round(np.quantile(arr,[.025,.975]),6)),
            "observed_rho":actual,
            "fraction_null_with_rho_at_least_observed":
                round(float(np.mean(arr>=actual)),6),
            "fraction_null_with_abs_rho_at_least_observed_abs":
                round(float(np.mean(np.abs(arr)>=abs(actual))),6),
            "NOT_A_VALID_RANDOMIZED_PVALUE_UNLESS_NODE_EXCHANGEABILITY_HOLDS":True
        }
    return out


def run(temporal,fa,flora82,previous86,census_raw):
    if (previous86.get("complete") is not True or
        previous86.get("n_evaluable_directed_origin_destination_source_pairs")!=2255 or
        previous86.get("n_both_GBIF_source_quality_US_states")!=48 or
        previous86.get("Census_source_sha256")!=sha256(census_raw).hexdigest()):
        raise ValueError("Verified 2026 GBIF 48-state and Census source no longer identical to original v8.6")
    if flora82.get("complete") is not True or fa.get("complete") is not True:
        raise ValueError("2022 FAF economic OD or 2026 GBIF original source missing")
    geo=parse_census_centroids(census_raw.decode("utf-8"),flora82["states"])
    rows,eligible,missing=construct_directions(temporal,fa,geo)
    if len(eligible)!=48 or len(rows)!=2255:
        raise ValueError("Original country-GRIIS indexed source and 2022 FAF directiveness drifted")
    for mode in MODELS:
        r=legacy_partial_rank(rows,RESPONSE,mode)
        baseline=previous86["all_results"][COHORTS[0]]["results"][RESPONSE][
            "rail_full" if mode=="RAIL" else "truck_full"]
        if r!=baseline:
            raise ValueError("Actual original v8.6 recorded-source association coefficients not replicated")
    od=original_crop_OD(fa)
    out={
        "complete":True,"source_utc":datetime.now(timezone.utc).isoformat(),
        "FAF6_crop_economic_OD_year":2022,
        "recorded_source_gbif_year_bins":"origin 2018-2021; destination 2023-2025",
        "n_qualifying_US_Census_state_sources":48,
        "n_directed_state_record_source_pairs_sharing_nodes":2255,
        "precommitted_model":"nonparametric residual rank correlation with origin and destination state fixed effects, actual other-mode OD, Census geodesic geography and source candidate richness",
        "random_source_node_relabeling_SEED":SEED,
        "permutation_counts":PERMUTATIONS,
        "verified_original_v86_rail_and_truck_source_coefficients":True,
        "original_Census_TIGER_source_SHA256":sha256(census_raw).hexdigest(),
        "source_cargo_mode_comparison":{},
        "NO_actual_viable_seeds_causal_railway_plant_colonization_adaptation":True,
    }
    for mode in MODELS:
        source=prepare_dyadic(rows,geo,other=("TRUCK" if mode=="RAIL" else "RAIL"))
        out["source_cargo_mode_comparison"][mode]=one_mode(mode,source,od[mode])
        print(json.dumps({"mode":mode,
           "source_fixed_effects_rhos":out["source_cargo_mode_comparison"][mode][
               "observed_origin_destination_FE_source_partial_rho"],
           "region_preserved":out["source_cargo_mode_comparison"][mode][
               "constrained_FA6_directed_network_node_label_null"]["Census_region_preserved"],
           "no_causal_plant_effect":True},ensure_ascii=False),flush=True)
    out["limitations"]=[
        "GBIF indexed speciesKey 'first-recorded after' is not true new invasive establishment; country-GRIIS taxon may be native to a specific state.",
        "FAF6 2022 state-economic crop cargo OD is NOT observed use of physical rail tracks, actual seed shipment, contamination or railway release.",
        "Origin/destination fixed effects absorb all state-level variation including potentially genuine effects; association is strictly remaining between-state directed OD source residuals.",
        "Node-label source-network permutations within Census region/division do NOT make historical agriculture, climate, rail networks or GBIF reporting exchangeable. Exceedance fractions are descriptive structural diagnostics, NOT guaranteed calibrated p-values or causal confidence bounds.",
        "The same 48 origin and destination source lists recur in all 2,255 dyads; node-label permutation preserves one form of dependence but NOT all network/event/observer dependencies.",
        "One 2022 freight base-year and uneven pre/post GBIF source event-year window lengths, no source/provider-fixed independent surveys, no difference-in-differences intervention.",
        "No measured railway seed traps, freight-route loaded propagule assays, rail/road verge mowing, train-generated gusts, phenotype plasticity or population genetic adaptation.",
        "This analysis does not estimate specific railway colonization propensity or any global macroecological law."
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for s in ("temporal85","freight80","flora82","baseline86","output"):
        p.add_argument("--"+s,required=True)
    a=p.parse_args();read=lambda x:json.loads(Path(getattr(a,x)).read_text(encoding="utf-8"))
    r=requests.get(CENSUS_TIGER_URL,timeout=(20,90),
        headers={"User-Agent":"tetsu-v89-origin-destination-state-FE-node-label-FAF6-GBIF"})
    r.raise_for_status()
    result=run(read("temporal85"),read("freight80"),read("flora82"),read("baseline86"),r.content)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,
      "models":{m:{"rho":v["observed_origin_destination_FE_source_partial_rho"],
                 "source_node_label_null":v["constrained_FA6_directed_network_node_label_null"]}
                for m,v in result["source_cargo_mode_comparison"].items()},
      "not_causal_freight_seed_establishment":True},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

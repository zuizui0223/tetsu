#!/usr/bin/env python3
"""v8.0: real 2022 BTS FAF6 state×state freight OD network comparison.

Primary: crop-linked SCTG02+03 RAIL vs TRUCK, strictly same domestic
origin/destination state frame. Separate SCTG02 grain and SCTG03 other
agriculture, and SCTG10-12 resource comparison. Calculate observed-source
transport-OPPORTUNITY topology: within/inter-state tonnage shares, OD
support, origin-partner diversity and directional imbalance.

This is NOT a railway track routing assignment or evidence of train-borne
viable non-native seeds. No plant distribution, mowing, wind exposure,
evolutionary adaptation or causal ecological effect is inferred.
"""
from __future__ import annotations
import argparse,json,math
from collections import defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
from statistics import median
import requests

from audit_BTS_FAF6_2022_agricultural_rail_freight_v79 import (
    URL,records_from_faf6_zip,normalize_code,source_float)

MODES={"01":"TRUCK","02":"RAIL"}
GROUPS={
    "grain_sctg02":frozenset(("02",)),
    "other_agriculture_sctg03":frozenset(("03",)),
    "crop_sctg02_03":frozenset(("02","03")),
    "resource_sctg10_12":frozenset(("10","11","12")),
}
BASELINE_GROUPS={
    "crop_sctg02_03":"crop_codes_02_03",
    "resource_sctg10_12":"resource_codes_10_11_12",
}
MIN_EXPECTED_STATE_CODES=50
MIN_SOURCE_OD_ROWS=500
THOUSAND_SHORT_TON_THRESHOLDS=(0,10,100)
FREIGHT_DATA_ARE_ESTIMATES_NOT_ROUTE_OBSERVATIONS=True


def safe(x,places=6):
    return round(float(x),places)


def build_networks(rows):
    """Never choose states or OD links using species records or rail effects."""
    flows=defaultdict(float)
    nall=0
    nvalid=0
    states=set()
    excluded=defaultdict(int)
    for row in rows:
        nall+=1
        trade=normalize_code(row.get("trade_type"))
        o=normalize_code(row.get("dms_origst"))
        d=normalize_code(row.get("dms_destst"))
        mode=normalize_code(row.get("dms_mode"))
        sctg=normalize_code(row.get("sctg2"))
        tons=source_float(row.get("tons_2022"))
        if o is None or d is None or sctg is None:
            excluded["missing_OD_or_commodity"]+=1
            continue
        states.update((o,d))
        if trade!="01":
            excluded["import_or_export_not_domestic_primary"]+=1
            continue
        if mode not in MODES:
            excluded["other_transport_modes_not_primary_comparator"]+=1
            continue
        groups=[g for g,ss in GROUPS.items() if sctg in ss]
        if not groups:
            excluded["not_predeclared_commodity_group"]+=1
            continue
        nvalid+=1
        if tons==0:
            excluded["zero_estimated_tonnage_row"]+=1
            continue
        for group in groups:
            flows[(group,MODES[mode],o,d)]+=tons
    if nall<MIN_SOURCE_OD_ROWS or len(states)<MIN_EXPECTED_STATE_CODES:
        raise ValueError("FAF6 state-OD original source population unexpected or truncated")
    if not flows:
        raise ValueError("No genuine FAF state OD network source rows")
    if len(states)!=51:
        raise ValueError("FAF6 state FIPS code universe not original verified 51")
    return dict(flows),sorted(states),{
        "original_source_rows":nall,
        "domestic_targeted_commodity_rail_or_truck_source_rows_including_zero":nvalid,
        "source_exclusion_counts":dict(excluded),
        "n_source_state_codes":len(states),
    }


def mode_network(flows,group,mode,states):
    # These are economic shipment OD weights, NOT observed routes on rail.
    matrix={(o,d):v for (g,m,o,d),v in flows.items()
            if g==group and m==mode and v>0}
    if not matrix:
        raise ValueError("Missing mode/commodity source OD shipment network")
    within=sum(v for (o,d),v in matrix.items() if o==d)
    between=sum(v for (o,d),v in matrix.items() if o!=d)
    total=within+between
    if total<=0:
        raise ValueError("No positive tonnage in original source")
    ordered_nonself=[(o,d,v) for (o,d),v in matrix.items() if o!=d]
    origin=defaultdict(dict)
    destination=defaultdict(float)
    for o,d,v in ordered_nonself:
        origin[o][d]=v
        destination[d]+=v
    partner_stats={}
    for state in states:
        links=origin.get(state,{})
        vol=sum(links.values())
        if vol<=0:
            continue
        probs=[v/vol for v in links.values()]
        shannon=math.exp(-sum(p*math.log(p) for p in probs))
        inverse_hhi=1/sum(p*p for p in probs)
        partner_stats[state]={
            "n_interstate_destination_codes_with_positive_OD_tonnage":len(links),
            "estimated_interstate_outgoing_thousand_short_tons":safe(vol),
            "effective_destination_count_shannon":safe(shannon),
            "effective_destination_count_simpson_inverseHHI":safe(inverse_hhi),
            "largest_destination_share":safe(max(probs)),
        }
    existing_unordered={tuple(sorted((o,d))) for o,d,v in ordered_nonself}
    mismatch=sum(abs(matrix.get((a,b),0)-matrix.get((b,a),0))
                 for a,b in sorted(existing_unordered))
    pairtotal=sum(matrix.get((a,b),0)+matrix.get((b,a),0)
                  for a,b in existing_unordered)
    bidirectional=sum(matrix.get((a,b),0)>0 and matrix.get((b,a),0)>0
                      for a,b in existing_unordered)
    thresholds={
        f"thousand_short_ton_threshold_gt_{t}":{
            "n_positive_directed_interstate_OD_state_pairs":
                sum(v>t for o,d,v in ordered_nonself),
            "n_state_origins_with_positive_interstate_outflow":
                len({o for o,d,v in ordered_nonself if v>t}),
            "n_state_destinations_with_positive_interstate_inflow":
                len({d for o,d,v in ordered_nonself if v>t}),
        } for t in THOUSAND_SHORT_TON_THRESHOLDS
    }
    return {
        "n_positive_OD_state_pairs_including_within_state":len(matrix),
        "n_positive_interstate_directed_OD_state_pairs":len(ordered_nonself),
        "total_2022_estimated_thousand_short_tons":safe(total),
        "intra_state_2022_estimated_thousand_short_tons":safe(within),
        "inter_state_2022_estimated_thousand_short_tons":safe(between),
        "share_interstate_in_mode":safe(between/total),
        "interstate_directional_weight_imbalance_0reciprocal_to1oneway":safe(mismatch/pairtotal) if pairtotal else None,
        "n_undirected_state_pairs_with_any_estimated_interstate_flow":len(existing_unordered),
        "n_state_pairs_with_both_OD_directions_positive":bidirectional,
        "n_origin_states_with_positive_interstate_flow":len(partner_stats),
        "n_destination_states_with_positive_interstate_flow":len(destination),
        "median_effective_outgoing_destination_count_shannon":safe(median(
            z["effective_destination_count_shannon"] for z in partner_stats.values())) if partner_stats else None,
        "by_origin_state_code":partner_stats,
        "source_OD_density_threshold_sensitivity":thresholds,
        "ORIGINAL_ECONOMIC_OD_TONNAGE_NOT_TRACK_ROUTING":True,
    }


def paired_origin_contrasts(a,b):
    """Equal state-code support, no fake n=51 independent ecological replicates."""
    rail=a["by_origin_state_code"]
    truck=b["by_origin_state_code"]
    states=sorted(set(rail)&set(truck))
    if not states:return {"n_eligible_state_sources":0}
    comparisons={}
    for name in ("effective_destination_count_shannon",
                 "n_interstate_destination_codes_with_positive_OD_tonnage",
                 "largest_destination_share"):
        differences=[rail[s][name]-truck[s][name] for s in states]
        comparisons[name]={
            "n_state_origin_codes_with_both_modes_positive_interstate_outgoing":len(states),
            "n_states_rail_greater_than_truck":sum(v>0 for v in differences),
            "n_states_rail_less_than_truck":sum(v<0 for v in differences),
            "equal_state_median_rail_minus_truck":safe(median(differences)),
            "equal_state_median_rail":safe(median(rail[s][name] for s in states)),
            "equal_state_median_truck":safe(median(truck[s][name] for s in states)),
        }
    return {
        "n_origin_states_with_both_RAIL_and_TRUCK_interstate_flow":len(states),
        "matched_origin_state_codes":states,
        "comparisons":comparisons,
        "NOT_independent_species_or_route_based_replication":True,
    }


def audit_source_invariants(out,previous):
    if previous.get("complete") is not True or previous.get("FAF_release")!="FAF6.0 / 2022 BASE YEAR / STATE OD CSV":
        raise ValueError("Original validated FAF6 v7.9 source missing")
    a=previous["data"]
    if out["source_trace"]["original_source_rows"]!=a["n_original_official_FAF6_state_OD_commodity_mode_rows"]:
        raise ValueError("FAF6 original CSV source-row population changed from frozen source")
    if out["source_trace"]["n_source_state_codes"]!=a["n_distinct_state_code_categories_in_source"]:
        raise ValueError("FAF6 source state code universe changed")
    for group,oldgroup in BASELINE_GROUPS.items():
        x=out["mode_networks"][group]
        baseline=a["national_trade_type_x_commodity_group_2022"]["domestic"][oldgroup]
        for mode in ("RAIL","TRUCK"):
            key="rail_thousand_short_tons" if mode=="RAIL" else "truck_thousand_short_tons"
            if abs(x[mode]["total_2022_estimated_thousand_short_tons"]-baseline[key])>0.0001:
                raise ValueError("FAF6 OD matrix marginal tonnage not conserved with independent v7.9 source")
    out["source_trace"]["previous_verified_v79_total_and_modes_reconciled"]=True
    return out


def run(previous,fetch=requests.get):
    response=fetch(URL,timeout=(30,180),headers={"User-Agent":"tetsu-official-FAF6-2022-directed-state-OD-source-v80"})
    response.raise_for_status()
    raw=response.content
    if sha256(raw).hexdigest()!=previous.get("official_published_CSV_ZIP_sha256"):
        raise ValueError("Source FAF6 original ZIP SHA differs from previously verified v7.9: stop for review")
    flows,states,trace=build_networks(records_from_faf6_zip(raw))
    out={
        "complete":True,
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_source_url":URL,
        "source_hash_verified_vs_prior_v79":True,
        "country":"US",
        "source_year":2022,
        "source_measurement_unit":"thousands US short tons of estimated freight by state of origin and destination",
        "source_trace":trace,
        "predeclared_commodity_groups":{g:sorted(v) for g,v in GROUPS.items()},
        "mode_networks":{},
        "matched_state_origin_rail_vs_truck":{},
        "NO_plant_community_railway_seed_evolution_inference":True,
    }
    for group in GROUPS:
        out["mode_networks"][group]={
            mode:mode_network(flows,group,mode,states) for mode in ("RAIL","TRUCK")}
        out["matched_state_origin_rail_vs_truck"][group]=paired_origin_contrasts(
            out["mode_networks"][group]["RAIL"],
            out["mode_networks"][group]["TRUCK"])
    audit_source_invariants(out,previous)
    out["raw_official_state_OD_2022_crop_estimated_tonnage"]={
        mode:[{"o":o,"d":d,"estimated_thousand_short_tons":safe(v)}
              for (group,m,o,d),v in sorted(flows.items()) if group=="crop_sctg02_03" and m==mode]
        for mode in ("RAIL","TRUCK")}
    out["limitations"]=[
        "FAF6 state OD weights are estimated ECONOMIC freight flows, not actual routes, line segment volumes, observed carriage of invasive seeds or railway verge colonization.",
        "Nonzero modeled flow is not an observed shipment or propagule. FAF covers 2022 domestic freight; no 2010–2025 rail plant colonization chronology can follow.",
        "Rail vs truck OD graphs share 51 states but different reporting/commodity suitability; this is descriptive topology, not a statistical rail biological effect.",
        "Large state-specific OD weights cannot be allocated to specific shortest physical Overture 50km railway edges without an independent route assignment model.",
        "Inter-state does NOT always mean long-distance, and a state code is an economic OD location, not a detailed ecological study site.",
        "Shannon effective number of destination states depends on estimated tonnage composition; not observed species diversity or a seed-dispersal distance measure.",
        "No measured railway grass mowing, train slipstream gust, phenotype, plasticity, population genetics, or local evolutionary adaptation is contained here.",
    ]
    print(json.dumps({"complete":True,
        "source_OD_rows":trace["original_source_rows"],
        "crop_RAIL":{k:v for k,v in out["mode_networks"]["crop_sctg02_03"]["RAIL"].items()
                     if k not in ("by_origin_state_code",)},
        "crop_TRUCK":{k:v for k,v in out["mode_networks"]["crop_sctg02_03"]["TRUCK"].items()
                     if k not in ("by_origin_state_code",)},
        "crop_paired_origins":out["matched_state_origin_rail_vs_truck"]["crop_sctg02_03"],
        "source_trace":trace,
        "noncausal":True},ensure_ascii=False),flush=True)
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--previous-v79",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    old=json.loads(Path(a.previous_v79).read_text(encoding="utf-8"))
    out=run(old)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    if not out["source_trace"]["previous_verified_v79_total_and_modes_reconciled"]:
        raise SystemExit("Source marginal conservation failed")

if __name__=="__main__":
    main()

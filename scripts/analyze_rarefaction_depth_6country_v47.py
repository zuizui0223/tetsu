#!/usr/bin/env python3
"""v4.7 six-country fixed-source GBIF GRIIS plant rarefaction-depth sensitivity.

Research question: did 5/6 RAW-positive vs 3/6 K=20-rarefied-positive
source association depend on *the arbitrary 20 indexed RECORD depth*?

Run the exact same pre-frozen six country 50km source grids, same source
RAIL/ROAD routes, period restrictions and rank residualization used in v4.4.
Evaluate rarefaction K=20/40/80/160/320 at:
 (A) depth-specific source-eligible cells, with potential cohort drift
 (B) fixed K80 source cohort for K=20/40/80, ensuring identical cells.
Other unqualified cells are missing, NOT actual plant absence or zero.
No per-pair p values, no actual independent survey samples, no claims of
train-mediated seed dispersal or biological abundance. K is indexed GBIF
RECORDS not organisms. This is a POSTHOC diagnostic, not preregistered.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from itertools import combinations
import json,math
from pathlib import Path
from statistics import median

from analyze_macro_50km_rail_road_aliens_v34 import (
    _eligible_both_period_cells, _routes, partial_rank_r
)
from analyze_6country_50km_network_plants_v44 import (
    BLOCKS_SIX, STRATEGIES, PERIODS,
    combine_source_regions, MIN_SOURCE_CELLS, MIN_VALID_PAIRS
)
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS

DEPTHS=(20,40,80,160,320)
FIXED_COHORT_MIN_DEPTH=80
MIN_COUNTRY_BLOCKS=5


def stable_record_rarefaction_probability(n,c,k):
    """P(taxon seen at least once in k-record without replacement draw).

    Uses sum(log1p(-c/(N-j))) to avoid cancellation in large-N four-gamma
    differences, and -expm1 for small probabilities. Exact to FP accuracy
    under the source hypergeometric model, not a survey-aware model.
    """
    if any(isinstance(x,bool) or not isinstance(x,int) for x in (n,c,k)):
        raise ValueError("Integer GBIF source record counts required")
    if n<1 or not(1<=c<=n) or not(1<=k<=n):
        raise ValueError("Impossible GBIF source record sample")
    if n-c<k:
        return 1.0
    logzero=math.fsum(math.log1p(-c/(n-i)) for i in range(k))
    p=-math.expm1(logzero)
    if not(0<=p<=1 and math.isfinite(p)):
        raise ValueError("Invalid rarefied source taxon detection probability")
    return p


def stable_probabilities(counts,k):
    if not isinstance(counts,dict) or not counts:
        raise ValueError("Missing observed introduced species source")
    if any(not isinstance(v,int) or isinstance(v,bool) or v<=0 for v in counts.values()):
        raise ValueError("Invalid GBIF indexed speciesKey counts")
    n=sum(counts.values())
    if n<k:
        raise ValueError("Insufficient source indexed introduced records for rarefaction")
    return {name:stable_record_rarefaction_probability(n,v,k)
            for name,v in counts.items()}


def expected_overlap_ratio(left,right):
    keys=set(left)|set(right)
    if not keys:return None
    intersection=math.fsum(left.get(s,0)*right.get(s,0) for s in keys)
    union=math.fsum(left.get(s,0)+right.get(s,0)-
                    left.get(s,0)*right.get(s,0) for s in keys)
    return intersection/union if union>0 else None


def network_source_pairs(block,rail,road,strategy):
    rr=_routes(rail,block,strategy,"pair_routes")
    dd=_routes(road,block,strategy,"pair_paths")
    if set(rr)!=set(dd):
        raise ValueError("RAIL/ROAD geographic source cell pairs do not match")
    out={}
    for key in sorted(rr):
        a,b=rr[key],dd[key]
        railratio=a.get("rail_to_geodesic_ratio")
        roadratio=b.get("road_to_geodesic_ratio")
        rp=a.get("cell_centroid_snap_augmented_rail_path_km")
        dp=b.get("centroid_snap_augmented_road_path_km")
        geo=a.get("source_cell_centroid_geodesic_km")
        if not all(isinstance(x,(int,float)) and not isinstance(x,bool) and
                   math.isfinite(x) and x>0 for x in (railratio,roadratio,rp,dp,geo)):
            continue
        if railratio<.95 or roadratio<.95:
            continue
        out[key]={
            "X":math.log(rp/dp),
            "geo":math.log(geo),
            "road":math.log(dp),
        }
    return out


def source_cell_cohort(block,flora):
    if block not in BLOCKS_SIX or flora["regions"][block].get("country")!=BLOCKS[block][0]:
        raise ValueError("Country or block changed from prior source experiment")
    cells=_eligible_both_period_cells(block,flora)
    for cid,periods in cells.items():
        for period in PERIODS:
            row=periods[period]
            if sum(row["counts"].values())<20:
                raise ValueError("GBIF source eligibility violated")
    return cells


def calculate_block(block,flora,rail,road):
    cells=source_cell_cohort(block,flora)
    pair_universe={st:network_source_pairs(block,rail,road,st)
                   for st in STRATEGIES}
    original20={k for k,v in cells.items() if all(
        sum(v[p]["counts"].values())>=20 for p in PERIODS)}
    maxcounts={
        k:min(sum(v[p]["counts"].values()) for p in PERIODS)
        for k,v in cells.items()
    }
    fixed80={k for k in original20 if maxcounts[k]>=FIXED_COHORT_MIN_DEPTH}
    result={
        "country":BLOCKS[block][0],
        "n_original_qualifying_cells_K20":len(original20),
        "n_fixed_K80_qualifying_cells":len(fixed80),
        "both_period_min_indexed_introduced_records_by_cell":maxcounts,
        "source_rarefaction_depths":list(DEPTHS),
        "fixed80_cohort_is_selected_by_record_counts_only":True,
        "results":{}
    }
    for strategy in STRATEGIES:
        network=pair_universe[strategy]
        result["results"][strategy]={}
        for mode in ("depth_specific","fixed80"):
            result["results"][strategy][mode]={}
            for k in DEPTHS:
                cohort={cid for cid in cells if
                        maxcounts[cid]>=k} if mode=="depth_specific" else fixed80
                if mode=="fixed80" and k>FIXED_COHORT_MIN_DEPTH:
                    # Do not change the fixed80 cohort to enable a larger K.
                    result["results"][strategy][mode][str(k)]={
                        "unsupported_above_fixed80_source_depth":True,
                        "descriptive_partial_rank":{}
                    }
                    continue
                profiles={
                    cid:{p:stable_probabilities(cells[cid][p]["counts"],k)
                         for p in PERIODS}
                    for cid in cohort
                }
                byperiod={}
                for period in PERIODS:
                    observations=[]
                    for a,b in combinations(sorted(cohort),2):
                        src=network.get((a,b))
                        if src is None:continue
                        y=expected_overlap_ratio(profiles[a][period],profiles[b][period])
                        if y is None or not math.isfinite(y):
                            continue
                        n1=cells[a][period]["all_plant_indexed_records"]
                        n2=cells[b][period]["all_plant_indexed_records"]
                        if n1<=0 or n2<=0:
                            raise ValueError("Invalid source all-plant recording opportunity")
                        observations.append({
                            "a":a,"b":b,
                            "x":src["X"],"y":y,
                            "cov":[src["geo"],src["road"],
                                   abs(math.log(n1/n2)),
                                   (math.log(n1)+math.log(n2))/2.]
                        })
                    rho=partial_rank_r(
                        [o["x"] for o in observations],
                        [o["y"] for o in observations],
                        [o["cov"] for o in observations]
                    ) if observations else None
                    byperiod[period]={
                        "n_eligible_cells":len(cohort),
                        "n_nonindependent_paired_network_observations":len(observations),
                        "n_unique_cells_with_at_least_one_network_pair":len(set(
                            o["a"] for o in observations)|set(o["b"] for o in observations)),
                        "descriptive_partial_rank":rho,
                        "source_supported":len(cohort)>=MIN_SOURCE_CELLS
                          and len(observations)>=MIN_VALID_PAIRS and rho is not None,
                        "not_survey_occupancy_or_train_seed_dispersal":True
                    }
                result["results"][strategy][mode][str(k)]={"by_period":byperiod}
                if strategy=="nearest_in_cell":
                    print(json.dumps({
                        "block":block,"country":BLOCKS[block][0],
                        "cohort":mode,"depth":k,
                        "eligible_cells":len(cohort),
                        "early_partial_rho":byperiod["2010_2017"]["descriptive_partial_rank"],
                        "late_partial_rho":byperiod["2018_2025"]["descriptive_partial_rank"],
                        "late_paired_records":byperiod["2018_2025"]["n_nonindependent_paired_network_observations"],
                        "late_source_supported":byperiod["2018_2025"]["source_supported"]
                    },ensure_ascii=False),flush=True)
    return result


def aggregate_country(blocks,period,strategy,mode,k):
    bycountry={}
    source={}
    for block,v in blocks.items():
        r=v["results"][strategy][mode][str(k)]
        if "unsupported_above_fixed80_source_depth" in r:
            continue
        row=r["by_period"][period]
        source[block]={
            "country":v["country"],
            "eligible_cells":row["n_eligible_cells"],
            "pairwise_observations_sharing_cells":row["n_nonindependent_paired_network_observations"],
            "support":row["source_supported"]
        }
        if row["source_supported"]:
            bycountry[block]=float(row["descriptive_partial_rank"])
    vals=list(bycountry.values())
    loo={b:round(median(x for name,x in bycountry.items() if name!=b),6)
         for b in bycountry if len(bycountry)>1}
    return {
        "depth_indexed_GBIF_records":k,
        "cohort":mode,
        "period":period,
        "source_country_blocks_total":len(blocks),
        "n_independent_supported_country_blocks":len(vals),
        "rail_shorter_more_similar_hypothesis_negative_direction":sum(x<0 for x in vals),
        "opposite_positive_direction":sum(x>0 for x in vals),
        "country_equal_weight_median_rho":round(median(vals),6) if vals else None,
        "min_5_country_source_coverage_gate":len(vals)>=MIN_COUNTRY_BLOCKS,
        "leave_one_country_out_median_rho":loo,
        "leave_one_out_direction_stable":bool(loo) and (all(x<0 for x in loo.values()) or all(x>0 for x in loo.values())),
        "country_rhos":bycountry,
        "country_source_coverage":source,
        "no_p_value_or_train_propagule_inference":True,
    }


def analyze(flora,rail,road):
    if not flora.get("complete") or flora.get("source_errors"):
        raise ValueError("Original full country-flora source incomplete")
    for frame in (flora,rail,road):
        if not set(BLOCKS_SIX).issubset(frame.get("regions",{})):
            raise ValueError("Not all six originally selected independent country blocks available")
    if any(set(v.get("regions",{}))!=set(BLOCKS_SIX) for v in (rail,road)):
        raise ValueError("Network source selection changed")
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "question":"Is the six-country RAW-vs-rarefied sign reversal driven by fixed 20 indexed RECORD standardization depth?",
        "posthoc_sensitivity_not_independent_preregistered_test":True,
        "six_country_blocks":list(BLOCKS_SIX),
        "fixed_source_coverage_depth":FIXED_COHORT_MIN_DEPTH,
        "depths":list(DEPTHS),
        "country_results":{},"country_equal_weight":{},
        "not_a_causal_railway_seed_effect":True,
    }
    for block in BLOCKS_SIX:
        out["country_results"][block]=calculate_block(block,flora,rail,road)
    for period in PERIODS:
        out["country_equal_weight"][period]={}
        for strategy in STRATEGIES:
            out["country_equal_weight"][period][strategy]={}
            for mode in ("depth_specific","fixed80"):
                out["country_equal_weight"][period][strategy][mode]={
                    str(k):aggregate_country(out["country_results"],period,strategy,mode,k)
                    for k in DEPTHS
                }
    out["caveats"]=[
        "All observations are GBIF indexed positive source records (not plants) and country GRIIS-list matching, not systematically observed absence or source DOI downloaded events.",
        "Changing depth changes source cohort in depth_specific results; fixed80 results hold identical botanical cells at k20/40/80 only.",
        "Even fixed rarefaction of indexed records does NOT control observer/dataset geography, sample-site duplication or metadata quality.",
        "An n=20 source minimum does not imply sample representativeness; normalized depth k does not define expected true ecological richness.",
        "Networks from Overture 2026 are not historical 2010-2025 active freight flows or direct seed-trap transport records.",
        "Rho uses nonindependent cell-pair rank residualization; country blocks are the replication unit and no cell pair p-values can be used.",
        "This posthoc descriptive analysis cannot establish actual railway seed/spore transport or absence of such a mechanism.",
    ]
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--flora",required=True)
    a.add_argument("--rail-sources",nargs="+",required=True)
    a.add_argument("--road-sources",nargs="+",required=True)
    a.add_argument("--output",default="results/audit/6country_GRIIS_rarefaction_depth_v47.json")
    args=a.parse_args()
    read=lambda paths:[json.loads(Path(x).read_text(encoding="utf-8")) for x in paths]
    flora=json.loads(Path(args.flora).read_text(encoding="utf-8"))
    rail=combine_source_regions(read(args.rail_sources))
    road=combine_source_regions(read(args.road_sources))
    out=analyze(flora,rail,road)
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    for period in PERIODS:
        for mode in ("depth_specific","fixed80"):
            for k in DEPTHS:
                row=out["country_equal_weight"][period]["nearest_in_cell"][mode][str(k)]
                print(json.dumps({"period":period,"mode":mode,"k":k,
                    "n_country_blocks":row["n_independent_supported_country_blocks"],
                    "negative_direction":row["rail_shorter_more_similar_hypothesis_negative_direction"],
                    "median_rho":row["country_equal_weight_median_rho"],
                    "leave_one_out_stable":row["leave_one_out_direction_stable"]
                },ensure_ascii=False),flush=True)
    print(json.dumps({"complete":True,"n_blocks":len(out["country_results"]),
          "posthoc_no_biological_causality":True}),flush=True)


if __name__=="__main__":
    main()

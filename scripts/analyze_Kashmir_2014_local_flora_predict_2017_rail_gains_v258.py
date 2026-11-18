#!/usr/bin/env python3
"""v25.8: does the 2014 LOCAL non-track flora predict taxa newly *recorded*
on the same Kashmir railway track in 2017, above region-wide species frequency?

Original Rashid et al. (2021) Dryad Matrix_plot.xlsx, fixed 2014/2017 field
locations and four 100m2 zones. Reuse source-authenticated v25.4 XLSX reader.
This is a predictive LOCAL-POOL association, NEVER proof of seed movement,
true immigration, recruitment, plant fitness or local adaptation.

For target z, source pool p, spatial site i:
  G_iz = species(z,2017) - species(z,2014)  # observed gains, not colonizations
  h_i  = |G_iz intersection p(2014,i)| / |G_iz| if G_iz nonempty
  e_i  = average of the same fraction across all OTHER source sites in its
         original station / between-station stratum.
  local_excess = mean_i(h_i - e_i), equal original field site weight.
The expected value is EXACT under uniform derangements; site permutation
draws only visualize/rank the finite label-assignment null.

Controls: lost target taxa seen in 2014 but not 2017, and embankment
2014->2017 gains measured against the same interior 2014 pool. Both
remove the claim that a gain-only association uniquely implies transport.

No original site × taxon data, source coordinates or plant lists exported.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import (
    parse_xlsx_matrix, matching_sites
)
from analyze_rail_local_site_coupling_v257 import derangement

YEARS = (2014, 2017)
ZONES = ("Track", "Embankment", "Interior Plot 1", "Interior Plot 2")
GROUPS = ("station", "between")
N_DRAWS = 999
SEED = 258118
SCENARIOS = (
    ("Track", "interiors_2014"),
    ("Track", "embankment_2014"),
    ("Track", "all_nontrack_2014"),
    ("Embankment", "interiors_2014"),
)


def source_pool(plots: dict, site: str, definition: str) -> frozenset[int]:
    def get(z: str) -> set[int]:
        x=plots[(site,2014,z)]
        if x is None:raise ValueError("Source zone missing")
        return x
    if definition=="interiors_2014":
        return frozenset(get("Interior Plot 1")|get("Interior Plot 2"))
    if definition=="embankment_2014":
        return frozenset(get("Embankment"))
    if definition=="all_nontrack_2014":
        return frozenset(get("Embankment")|get("Interior Plot 1")|get("Interior Plot 2"))
    raise ValueError("Unknown source pool, cannot silently impute")


def time_records(plots: dict, sites: list[str], zone: str) -> dict[str, dict]:
    out={}
    for s in sites:
        before,after=plots[(s,2014,zone)], plots[(s,2017,zone)]
        if before is None or after is None:raise ValueError("Source temporal target not present")
        out[s]={"gained":frozenset(after-before),
                "lost":frozenset(before-after),
                "retained":frozenset(after&before),
                "before_n":len(before),"after_n":len(after)}
        if (len(out[s]["gained"])+len(out[s]["retained"]) != len(after)
            or len(out[s]["lost"])+len(out[s]["retained"]) != len(before)):
            raise AssertionError("Temporal gain/loss identity invalid")
    return out


def hit_fraction(target: frozenset[int], donor: frozenset[int]) -> float | None:
    return len(target & donor)/len(target) if target else None


def observed_and_analytic_expected(
    records: dict[str,dict], pools: dict[str,frozenset[int]], category: str
) -> dict:
    """Exact per-site derangement expectation (uniform over n-1 other donors)."""
    sites=sorted(records)
    if set(sites)!=set(pools) or len(sites)<3:
        raise ValueError("No common site/donor stratum")
    eligible=[s for s in sites if records[s][category]]
    if not eligible:
        return {"identified": False, "n_eligible_sites":0,
                "n_target_taxon_events":0}
    actual=[];expected=[];hit_events=0;event_n=0
    for s in eligible:
        target=records[s][category]
        h=hit_fraction(target,pools[s])
        e=mean(hit_fraction(target,pools[t]) for t in sites if t!=s)
        actual.append(h);expected.append(e)
        hit_events+=len(target&pools[s]); event_n+=len(target)
    return {
        "identified":True,
        "n_original_sites_all":len(sites),
        "n_eligible_sites":len(eligible),
        "n_target_taxon_events":event_n,
        "observed_local_source_hit_events":hit_events,
        "observed_event_weighted_hit_share":hit_events/event_n,
        "observed_site_equal_hit_share":mean(actual),
        "source_other_site_exact_expected_site_equal":mean(expected),
        "local_excess_exact":mean(actual)-mean(expected),
        "source_pool_size_site_mean":mean(len(pools[s]) for s in sites),
        "site_specific_excess_min":min(a-b for a,b in zip(actual,expected)),
        "site_specific_excess_max":max(a-b for a,b in zip(actual,expected)),
        "original_site_recorded_counts": {
            "mean_2014":mean(records[s]["before_n"] for s in sites),
            "mean_2017":mean(records[s]["after_n"] for s in sites)
        }
    }


def permuted_score(
    records: dict[str,dict],pools: dict[str,frozenset[int]],
    category: str,mapping: dict[str,str]
) -> float | None:
    use=[s for s in records if records[s][category]]
    if not use:return None
    return mean(hit_fraction(records[s][category],pools[mapping[s]]) for s in use)


def qtile(values: list[float],p:float) -> float:
    x=sorted(values)
    if not x:raise ValueError("No original-source derangement draws")
    k=p*(len(x)-1);i=int(k)
    return x[i] if i==len(x)-1 else x[i]*(1-(k-i))+x[i+1]*(k-i)


def technical_diagnostic(observed: float, null: list[float]) -> dict:
    if len(null)<99:raise ValueError("Insufficient original-source label permutations")
    n=len(null)
    return {
        "draws":n,
        "null_mean":mean(null),
        "null_q025_technical":qtile(null,.025),
        "null_q975_technical":qtile(null,.975),
        "n_permutations_ge_observed":sum(x>=observed-1e-12 for x in null),
        "upper_tail_rank_plusone_mapping_only_NOT_field_p":
            (1+sum(x>=observed-1e-12 for x in null))/(n+1),
        "not_a_causal_effect_or_independent_field_CI":True
    }


def group_analysis(plots: dict,sites: list[str],draws:int,seed:int) -> dict:
    sites=sorted(sites)
    if len(sites)<10:raise ValueError("Original 4-zone site support under 10")
    record={zone:time_records(plots,sites,zone) for zone in ("Track","Embankment")}
    pool={label:{s:source_pool(plots,s,label) for s in sites}
          for label in ("interiors_2014","embankment_2014","all_nontrack_2014")}
    rng=random.Random(seed)
    actual={}
    for zone,definition in SCENARIOS:
        key=f"{zone}__{definition}"
        actual[key]={
            "gain":observed_and_analytic_expected(record[zone],pool[definition],"gained"),
            "loss_diagnostic":observed_and_analytic_expected(record[zone],pool[definition],"lost"),
            "retained_diagnostic":observed_and_analytic_expected(record[zone],pool[definition],"retained"),
        }
    null={key:{"gain":[], "loss_diagnostic":[]} for key in actual}
    gaps=[]
    for _ in range(draws):
        # Exactly the same re-pairing for all target zones and source pools.
        mapping=derangement(sites,rng)
        vals={}
        for zone,definition in SCENARIOS:
            key=f"{zone}__{definition}"
            for category,name in (("gained","gain"),("lost","loss_diagnostic")):
                v=permuted_score(record[zone],pool[definition],category,mapping)
                if v is not None:null[key][name].append(v)
            vals[key]=permuted_score(record[zone],pool[definition],"gained",mapping)
        a,b=vals["Track__interiors_2014"],vals["Embankment__interiors_2014"]
        if a is not None and b is not None:gaps.append(a-b)
    for key,v in actual.items():
        for section in ("gain","loss_diagnostic"):
            p=v[section]
            if p["identified"]:
                p["permutation_reference"]=technical_diagnostic(
                    p["observed_site_equal_hit_share"],null[key][section])
                if abs(p["source_other_site_exact_expected_site_equal"]-
                       mean(null[key][section]))>0.065:
                    raise AssertionError("Monte Carlo derangement average unexpectedly far from exact expectation")
    k1="Track__interiors_2014"; k2="Embankment__interiors_2014"
    p=actual[k1]["gain"];q=actual[k2]["gain"]
    contrasted=None
    if p["identified"] and q["identified"]:
        contrasted={
            "observed_track_minus_embankment_local_hit_share":
                p["observed_site_equal_hit_share"]-q["observed_site_equal_hit_share"],
            "track_minus_embankment_local_excess_above_exact_own_null":
                p["local_excess_exact"]-q["local_excess_exact"],
            "same_permutation_null":technical_diagnostic(
                p["observed_site_equal_hit_share"]-q["observed_site_equal_hit_share"],gaps),
            "not_equivalent_target_richness_or_detection":True
        }
    return {
        "original_source_sites":len(sites),
        "all_4_zones_both_original_years_complete":True,
        "source_set_fixed_at_2014_not_2017":True,
        "source_matched_null":"Uniform fixed-point-free reassignment of the entire 2014 source flora among same group, one-to-one, holding all original source sets and all target gains fixed",
        "scenario_metrics":actual,
        "track_vs_embankment_using_identical_interior_2014_source_pool":contrasted,
    }


def run(path:Path,draws:int=N_DRAWS) -> dict:
    if draws<99:raise ValueError("99 or more mappings required")
    plots,qc=parse_xlsx_matrix(path)
    cohorts={g:matching_sites(plots,g,ZONES) for g in GROUPS}
    if len(cohorts["station"])!=15 or len(cohorts["between"])!=12:
        raise ValueError("Previously verified original Kashmir 4-zone common cohort changed")
    result={g:group_analysis(plots,sites,draws,SEED+i*119)
            for i,(g,sites) in enumerate(cohorts.items())}
    return {
        "version":"v25.8",
        "status":"ACTUAL ORIGINAL DRYAD CALCULATION_IF_WORKFLOW_SUCCEEDS",
        "source":"Rashid et al. 2021, official Dryad version 99735",
        "original_sha256_verified":qc["original_file_sha256"],
        "source_qc":qc,
        "results_by_stratum":result,
        "source_baseline_year":2014,
        "newly_recorded_2017_not_proven_arrivals":True,
        "absence_in_2014_is_not_a_verified_biological_absence":True,
        "no_individual_seed_carriage_or_fitness_or_mowing_or_train_volume_measured":True,
        "comparison_uses_same_27_original_spatial_sites_and_both_original_years":True,
        "not_an_independent_country_replication":True,
        "no_original_plant_site_matrix_or_taxon_IDs_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=N_DRAWS)
    a=p.parse_args()
    result=run(a.source,a.draws)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for g,val in result["results_by_stratum"].items():
        t=val["scenario_metrics"]["Track__interiors_2014"]["gain"]
        e=val["scenario_metrics"]["Embankment__interiors_2014"]["gain"]
        print(json.dumps({
            "version":"v25.8","group":g,"source_spatial_sites":val["original_source_sites"],
            "track_gain_taxon_events":t.get("n_target_taxon_events"),
            "track_gain_sites":t["n_eligible_sites"],
            "track_gain_local_excess":t.get("local_excess_exact"),
            "track_gain_local_share":t.get("observed_site_equal_hit_share"),
            "track_gain_pool_expectation":t.get("source_other_site_exact_expected_site_equal"),
            "embank_gain_local_excess":e.get("local_excess_exact"),
            "track_vs_embank_excess":(val["track_vs_embankment_using_identical_interior_2014_source_pool"] or {}).get("track_minus_embankment_local_excess_above_exact_own_null")
        },ensure_ascii=False))


if __name__=="__main__":
    main()

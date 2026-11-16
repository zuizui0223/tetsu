#!/usr/bin/env python3
"""v10.6 real Slovenian railway 2-trait Kitagawa decomposition.

Question: is the v10.4-v10.5 BS−BT deficit of species with strict
wind-diaspore strategy a compositional consequence of BS annual enrichment?
Separate (a) annual-versus-other-known-lifespan mixing from (b) wind
fraction shifts *within* these two source-classified lifespan categories.
Original 2026 Slovenia railway field data, DOI 10.1007/s44473-026-00155-z.
This is not an experiment or measurement of true train-borne seeds.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_wind_diaspore_vs_annual_v104 import (
    WIND_STRA, VALID_SOURCE_CATEGORIES, SOURCE_DOI, SOURCE_SHA256,
    extract_original_xlsx, validate_original_keys, pair_source,
    measure, location_level,
)
from analyze_Slovenia_rail_v105_cluster_missingness import (
    segment_cluster_diagnostic, verify_against_locked_v104
)

WIND_VALUES=frozenset(("Epilobium","Phragmites"))
COMPONENTS=(
    "joint_known_wind_BS_minus_BT",
    "annual_mix_component",
    "within_lifespan_wind_component",
    "within_annual_component",
    "within_other_known_lifespan_component",
    "source_strict_wind_delta_on_exact_same_eligible_sides",
)


def joint_trait_counts(codes, traits):
    """Retain ONLY original taxa with known wind strategy AND known lifespan."""
    bins={"strict_annual":{"n":0,"wind":0},
          "other_known_lifespan_including_mixed":{"n":0,"wind":0}}
    missing_strategy=missing_lifespan=0
    for c in codes:
        r=traits[c]
        strategy=r.get("Diaspore.strategy")
        lifespan=r.get("Lifespan.categories")
        if strategy is None or strategy=="":
            missing_strategy+=1
        if lifespan is None or lifespan=="":
            missing_lifespan+=1
        if strategy is None or strategy=="" or lifespan is None or lifespan=="":
            continue
        if strategy not in VALID_SOURCE_CATEGORIES:
            raise ValueError("Unexpected original 2026 Diaspore.strategy class")
        cls="strict_annual" if lifespan=="annual" else "other_known_lifespan_including_mixed"
        bins[cls]["n"]+=1
        bins[cls]["wind"]+=int(strategy in WIND_VALUES)
    return {"bins":bins,"n_original_zone_taxa":len(codes),
            "n_source_missing_strategy":missing_strategy,
            "n_source_missing_lifespan":missing_lifespan,
            "n_original_taxa_with_both_source_traits":
                sum(v["n"] for v in bins.values())}


def wind_decomposition(bs, bt):
    """Exact two-stratum symmetric Kitagawa identity for two trait-known zones.

    Wind fraction = p_annual * wind_given_annual +
                    (1-p_annual) * wind_given_other_known_lifespan.
    Only identifiable where both strata occur at BOTH original zones.
    """
    result={}
    for label,obj in (("BS",bs),("BT",bt)):
        annual=obj["bins"]["strict_annual"]
        other=obj["bins"]["other_known_lifespan_including_mixed"]
        na, no=annual["n"], other["n"]
        if na<=0 or no<=0:
            return None
        pa=na/(na+no)
        qa=annual["wind"]/na
        qo=other["wind"]/no
        result[label]=(pa,qa,qo,pa*qa+(1-pa)*qo)
    p1,q1,o1,w1=result["BS"]
    p0,q0,o0,w0=result["BT"]
    dp=p1-p0
    pq=(p1+p0)/2
    annual_mix=dp*((q1+q0)/2-(o1+o0)/2)
    within_annual=pq*(q1-q0)
    within_other=(1-pq)*(o1-o0)
    within=within_annual+within_other
    observed=w1-w0
    if not math.isclose(observed,annual_mix+within,abs_tol=1e-12):
        raise ValueError("Invalid exact composition+within-lifespan decomposition")
    return {
        "joint_known_wind_BS_minus_BT":observed,
        "annual_mix_component":annual_mix,
        "within_lifespan_wind_component":within,
        "within_annual_component":within_annual,
        "within_other_known_lifespan_component":within_other,
        "joint_known_strict_annual_share_BS_minus_BT":dp,
        "original_exact_algebra_identity":True,
    }


def composition_sensitivity(original, prev103, prev104, prev105):
    if (prev105.get("complete") is not True or
        prev105.get("version") != "v10.5" or
        prev105.get("original_supplement_sha256") != SOURCE_SHA256):
        raise ValueError("Original source section bootstrap v10.5 receipt missing")
    rebuilt=verify_against_locked_v104(original,prev103,prev104)
    strict=rebuilt["observed_within_location_BS_minus_BT_group_differences"][WIND_STRA]
    ref=prev105["section_cluster_bootstrap_diagnostics"]["strict_wind_diaspore"][
        "all_15_source_sections"]["equal_field_location_mean_BS_minus_BT"]
    if strict["mean_BS_minus_BT_fraction"]!=ref:
        raise ValueError("Actual source locked v10.5 strict wind estimate changed")
    geos,zones,species,traits,_,_=validate_original_keys(original)
    paired=pair_source(geos,zones,species,traits)
    measured=measure(paired,species,traits)
    prev_by={(r["location_id"],r["track_side"]):r for r in measured}
    original_trait_taxonomy=Counter()
    for code,tr in traits.items():
        lifespan=tr.get("Lifespan.categories")
        d=tr.get("Diaspore.strategy")
        ls="unknown" if not lifespan else ("strict_annual" if lifespan=="annual"
            else "other_known_lifespan_including_mixed")
        disp="unknown" if not d else ("strict_wind" if d in WIND_VALUES
            else "other_known_strategy")
        original_trait_taxonomy[(ls,disp)]+=1
    valid=[]
    exclusions=Counter()
    for pair in paired:
        lid,side=pair["location_id"],pair["track_side"]
        bs=joint_trait_counts(species[(lid,side+"BS")],traits)
        bt=joint_trait_counts(species[(lid,side+"BT")],traits)
        d=wind_decomposition(bs,bt)
        if d is None:
            exclusions["no_annual_or_no_other_in_at_least_one_zone"]+=1
            continue
        matched=prev_by[(lid,side)]
        # The prior strict wind comparison uses ALL strategy-known taxa.
        # This joint-coded analysis uses ONLY BOTH trait-known taxa.
        d["source_strict_wind_delta_on_exact_same_eligible_sides"]=matched[
            "wind_and_diaspore_group_shares_BS_minus_BT"][WIND_STRA]
        d["location_id"]=lid
        d["section_number"]=pair["section_number"]
        valid.append(d)
    by=defaultdict(list)
    for r in valid:
        by[r["location_id"]].append(r)
    summary=[]
    for lid, rows in sorted(by.items()):
        row={"location_id":lid, "section_number":rows[0]["section_number"],
             "n_original_track_side_pairs_contributing":len(rows)}
        if any(r["section_number"]!=row["section_number"] for r in rows):
            raise ValueError("Railway track segment composition error")
        for k in COMPONENTS:
            vals=[r[k] for r in rows if r[k] is not None]
            row[k]=mean(vals) if vals else None
        summary.append(row)
    if sum(r["n_original_track_side_pairs_contributing"] for r in summary)!=len(valid):
        raise ValueError("Pseudoreplicated original track side count")
    if any(not math.isclose(r["joint_known_wind_BS_minus_BT"],
                            r["annual_mix_component"]+
                            r["within_lifespan_wind_component"],abs_tol=1e-12)
           for r in summary):
        raise ValueError("Location-level decomposition identity violated")
    metrics={k:segment_cluster_diagnostic(summary,k) for k in COMPONENTS}
    sections=Counter(str(x["section_number"]) for x in summary)
    adequate={k for k,n in sections.items() if n>=3}
    metrics_adequate={k:segment_cluster_diagnostic(
        [r for r in summary if str(r["section_number"]) in adequate], k)
        for k in COMPONENTS}
    withincode=dict(sorted({ls+"__"+d:n for (ls,d),n
                            in original_trait_taxonomy.items()}.items()))
    return {
        "complete":True,
        "version":"v10.6",
        "source_original_DOI":SOURCE_DOI,
        "source_original_xlsx_sha256":SOURCE_SHA256,
        "all_original_railway_source_field_locations":89,
        "all_original_true_BS_BT_sides":173,
        "original_side_pairs_with_both_lifespan_strata_identifiable":len(valid),
        "original_side_pairs_excluded_for_stratum_sparsity":dict(exclusions),
        "original_distinct_field_locations_retained":len(summary),
        "original_segments_retained":len(sections),
        "original_rail_segments_with_3plus_retained_sites":len(adequate),
        "source_312_taxon_lifespan_x_diaspore_class_inventory":withincode,
        "original_same_site_means_15_section_bootstrap":metrics,
        "source_subset_sections_3plus_sites_sensitivity":metrics_adequate,
        "descriptive_interpretation_rule":(
            "annual_mix accounts for wind gap due to lifespan reshuffling; "
            "within_lifespan remains after controlling for two original jointly "
            "known lifespan strata. A within term is still observational and "
            "cannot identify arrival, train gusts, mowing, plasticity or selection."),
        "caveat":(
            "Other known lifespan includes mixed annual/biennial codes, not "
            "all are strict perennials. Excluding pairs without both groups "
            "may change site coverage; no zero-denominator imputations. "
            "One-country nonrandom railway plots, unequal zone sampled areas, "
            "and original field study authors credited; do not infer global "
            "railway causal effect or seed vector."),
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","verified103",
                 "verified104","verified105","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    src103=json.loads(Path(a.verified103).read_text(encoding="utf-8"))
    src104=json.loads(Path(a.verified104).read_text(encoding="utf-8"))
    src105=json.loads(Path(a.verified105).read_text(encoding="utf-8"))
    r=composition_sensitivity(original,src103,src104,src105)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

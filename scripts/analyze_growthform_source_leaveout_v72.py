#!/usr/bin/env python3
"""v7.2: REAL GIFT herbaceous/woody GRIIS source flora, full vs excluding the
largest indexed GBIF all-plant dataset on exactly the same 50km cell pairs.
RAIL/ROAD paths are physical 2026 geography, NOT observed freight/seed flow.
No wind, mowing, plasticity, evolution or train dispersal inference.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from datetime import datetime,timezone
from statistics import median
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,PERIODS,STRATEGIES,combine_source_regions
from analyze_macro_50km_rail_road_aliens_v34 import _eligible_both_period_cells
from analyze_GRIIS_topdataset_leaveout_rail_macro_v57 import merge_species_sources
from analyze_sixcountry_functional_rail_v71 import GROUPS,validate_traits,eligible_group_counts,make_pair_rows,rho_for_rows

COUNT_KEYS={"full":"full_source_current_country_GRIIS_aliens",
            "without_dominant_dataset":"remaining_datasets_country_GRIIS_aliens"}

def fixed_group_source_cells(block,period,flora,split,trait_map):
    old=_eligible_both_period_cells(block,flora)
    source=split["source_cells"]
    if len(source)!=16 or not set(old).issubset(source):
        raise ValueError("Dataset source cells do not match frozen 16-cell geography")
    good={}
    for cid,values in old.items():
        p=source[cid]
        if not p.get("original_two_period_country_GRIIS_eligible") is True:
            raise ValueError("Country original two-period GRIIS source gate missing")
        row=p.get("periods",{}).get(period)
        if not row or row.get("dataset_partition_is_count_of_indexed_speciesKey_records_NOT_independent_surveys") is not True:
            raise ValueError("GBIF country speciesKey dataset split not independently verified")
        for key in COUNT_KEYS.values():
            counts=row[key]
            if not isinstance(counts,dict) or any(type(n)!=int or n<=0 for n in counts.values()):
                raise ValueError("Malformed indexed source species records")
        if all(eligible_group_counts({"counts":row[key]},trait_map,g) is not None
               for key in COUNT_KEYS.values() for g in GROUPS):
            good[cid]=values
    return good

def compare_country(block,period,strategy,flora,split,rail,road,traits):
    if split.get("complete") is not True or split.get("source_errors"):
        raise ValueError("Country GBIF provider partition incomplete")
    if split.get("source_geo_block")!=block or split.get("source_country")!=flora["regions"][block]["country"]:
        raise ValueError("Country/source geography mismatch")
    shared=fixed_group_source_cells(block,period,flora,split,traits)
    records={}
    for variant,key in COUNT_KEYS.items():
        cells={}
        for cid,old in shared.items():
            row={p:dict(v) for p,v in old.items()}
            row[period]["counts"]=split["source_cells"][cid]["periods"][period][key]
            cells[cid]=row
        pairrows,_=make_pair_rows(block,period,strategy,cells,rail,road,traits)
        filtered={(r["a"],r["b"]):r for r in pairrows if all(r[g] is not None for g in GROUPS)}
        if len(filtered)!=sum(all(r[g] is not None for g in GROUPS) for r in pairrows):
            raise ValueError("Duplicate geographic pair")
        records[variant]=filtered
    pairs=sorted(set.intersection(*(set(v) for v in records.values())))
    for p in pairs:
        a,b=(records[key][p] for key in COUNT_KEYS)
        if a["x"]!=b["x"] or a["covars"]!=b["covars"]:
            raise ValueError("Actual RAIL/ROAD paths and source covariates differ across provider exclusion")
    n_unique_cells=len({c for p in pairs for c in p})
    estimates={}
    for variant,allrows in records.items():
        subset=[allrows[p] for p in pairs]
        estimates[variant]={g:rho_for_rows(subset,g)["partial_rho"] for g in GROUPS}
    deltas={}
    for variant,g in estimates.items():
        deltas[variant]=(round(g["herbaceous"]-g["woody"],6)
                         if all(g[x] is not None for x in GROUPS) else None)
    return {"country":flora["regions"][block]["country"],
            "n_same_provider_and_both_group_source_cells":len(shared),
            "n_common_physical_pair_observations_NOT_independent":len(pairs),
            "n_unique_source_cells_in_common_pairs":n_unique_cells,
            "full_vs_without_dominant_source_partial_rho":estimates,
            "herb_minus_woody_same_pair_rho_delta":deltas,
            "no_pairwise_p_values_or_seed_transport_claim":True}

def equal_country_summary(regions,period,strategy):
    output={}
    for category in (*GROUPS,"herb_minus_woody"):
        output[category]={}
        for variant in COUNT_KEYS:
            values={}
            for block,data in regions.items():
                r=data["periods"][period][strategy]
                x=(r["herb_minus_woody_same_pair_rho_delta"][variant] if category=="herb_minus_woody"
                   else r["full_vs_without_dominant_source_partial_rho"][variant][category])
                if x is not None:values[block]=x
            output[category][variant]={
                "country_coefficients":values,
                "n_independent_country_blocks":len(values),
                "n_negative_rho":sum(v<0 for v in values.values()),
                "n_positive_rho":sum(v>0 for v in values.values()),
                "median_equal_country":round(median(values.values()),6) if values else None,
                "leave_one_country_out_medians":{
                    c:round(median(v for d,v in values.items() if d!=c),6)
                    for c in values if len(values)>1},
                "replication_five_country_gate":len(values)>=5}
    return output

def run(flora,gift,split,rail,road):
    traits=validate_traits(gift)
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Frozen GBIF country-GRIIS flora source incomplete")
    if any(set(x["regions"])!=set(BLOCKS_SIX) for x in (split,rail,road)):
        raise ValueError("Country source block identity differs from six-country frozen baseline")
    result={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "study":"Six-country GIFT growth form vs physical rail/road recorded GRIIS species similarity with largest indexed dataset removed",
        "primary_period":"2018_2025","primary_anchor":"nearest_in_cell",
        "K_indexed_records_per_growth_form":20,
        "identical_country_cell_pairs_across_both_source_variants_and_growth_groups":True,
        "regions":{},"equal_country_summary":{},
        "NO_rail_freight_seed_wind_mowing_adaptation_identified":True}
    for block in BLOCKS_SIX:
        result["regions"][block]={"periods":{}}
        for period in PERIODS:
            result["regions"][block]["periods"][period]={
                strategy:compare_country(block,period,strategy,flora,
                                         split["regions"][block],rail,road,traits)
                for strategy in STRATEGIES}
    for period in PERIODS:
        result["equal_country_summary"][period]={
            strategy:equal_country_summary(result["regions"],period,strategy)
            for strategy in STRATEGIES}
    result["limitations"]=[
        "Source dataset UUID does not mean independent observers, plants, or repeated standardized surveys.",
        "Indexed species keys are positive recorded taxa; absence and true seed movement not measured.",
        "2026 physical RAIL/ROAD paths are not historical active freight or train operations.",
        "GIFT growth form does not measure mowing, gust intensity, plasticity or genetic evolution.",
        "Strict same-provider plus herbaceous/woody country-cell eligibility can reduce the number of geographic replicate blocks relative to original v7.1.",
        "Do not treat dependent within-region cell pairings as independent p-value samples or select favorable anchor/period."
    ]
    return result

def main():
    ap=argparse.ArgumentParser()
    for f in ("flora","gift","output"):ap.add_argument("--"+f,required=True)
    for f in ("split-sources","rail-sources","road-sources"):ap.add_argument("--"+f,nargs="+",required=True)
    a=ap.parse_args();read=lambda x:json.loads(Path(x).read_text(encoding="utf-8"))
    data=run(read(a.flora),read(a.gift),
             merge_species_sources([read(p) for p in a.split_sources]),
             combine_source_regions([read(p) for p in a.rail_sources]),
             combine_source_regions([read(p) for p in a.road_sources]))
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":True,"period":"2018_2025","K":20,
      "nearest_in_cell":data["equal_country_summary"]["2018_2025"]["nearest_in_cell"],
      "not_causal_seed_or_stress_inference":True},ensure_ascii=False),flush=True)

if __name__=="__main__":main()

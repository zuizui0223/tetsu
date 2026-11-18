#!/usr/bin/env python3
"""v25.6 rail microhabitat beta: species replacement vs nestedness.

Two independent genuine railway FIELD STUDIES: Kashmir track / embankment /
two inland plots (2014 and 2017), and Slovenia matched ballast shoulder /
ballast toe (2022/23). Never pool raw areas, taxonomic nomenclature or sites.

This computes pairwise Jaccard beta partition (Baselga 2012):
  beta_jac = (b+c)/(a+b+c)
  beta_jtu = 2*min(b,c)/(a+2*min(b,c))
  beta_jne = beta_jac - beta_jtu
where a=shared species, b=exclusive railway-side taxa, c=exclusive
reference taxa. Nesting index is a mathematical partition; it does NOT
establish causal railway habitat filtering, railway seed carriage or fitness.

Kashmir: same original site and same zone pair in both field years; sites
are retained for all 4 original zones to ensure exactly the same sample.
Station vs between is a field-design category, NOT trains per day.
Slovenia: original same-location same-side BS/BT paired species sets;
side averages -> original location -> original railway section. Survey
areas for BS/BT are NOT verified identical and are not standardized.

No original species-by-plot records or source GPS are published here.
"""
from __future__ import annotations
import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import (
    parse_xlsx_matrix, matching_sites
)
from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx, validate_original_keys
)

INDIA_GROUPS = ("station", "between")
INDIA_YEARS = (2014, 2017)
INDIA_ZONES = ("Track", "Embankment", "Interior Plot 1", "Interior Plot 2")
INDIA_REFERENCES = ("Embankment", "Interior Plot 1", "Interior Plot 2")
METRICS = (
    "rail_richness", "reference_richness", "richness_delta", "shared_taxa",
    "exclusive_rail_taxa", "exclusive_reference_taxa", "beta_jac",
    "beta_replacement_jtu", "beta_nestedness_jne", "replacement_fraction_of_beta",
    "strict_rail_subset", "fraction_of_rail_taxa_shared"
)
BOOT_METRICS = (
    "richness_delta", "exclusive_rail_taxa", "beta_jac",
    "beta_replacement_jtu", "beta_nestedness_jne",
    "replacement_fraction_of_beta", "strict_rail_subset"
)
SEED = 256118


def components(rail: set, reference: set) -> dict:
    if not rail or not reference:
        raise ValueError("Pairwise dissimilarity requires both actually recorded nonempty plots")
    a = len(rail & reference)
    b = len(rail - reference)
    c = len(reference - rail)
    union = a + b + c
    m = min(b, c)
    jac = (b+c)/union
    turnover = 2*m/(a+2*m) if a+2*m else 0.0
    nestedness = jac - turnover
    if nestedness < -1e-12 or turnover < -1e-12 or jac > 1+1e-12:
        raise ArithmeticError("Baselga partition must be in [0,1]")
    if abs((turnover + nestedness) - jac) > 1e-12:
        raise ArithmeticError("Pairwise Jaccard decomposition failed")
    return {
        "rail_richness": len(rail),
        "reference_richness": len(reference),
        "richness_delta": len(rail)-len(reference),
        "shared_taxa": a,
        "exclusive_rail_taxa": b,
        "exclusive_reference_taxa": c,
        "beta_jac": jac,
        "beta_replacement_jtu": turnover,
        "beta_nestedness_jne": max(0.0, nestedness),
        "replacement_fraction_of_beta": turnover/jac if jac > 0 else None,
        "strict_rail_subset": float(b == 0),
        "fraction_of_rail_taxa_shared": a/len(rail),
    }


def average(records: list[dict]) -> dict:
    if not records:
        raise ValueError("No actual independent sites to average")
    out = {}
    for k in METRICS:
        valid = [x[k] for x in records if x[k] is not None]
        out[k] = mean(valid) if valid else None
    out["observed_field_unit_count"] = len(records)
    out["positive_beta_unit_count"] = sum(x["replacement_fraction_of_beta"] is not None for x in records)
    return out


def qtile(values: list[float], p: float) -> float:
    v = sorted(values)
    if not v:
        raise ValueError("No valid resampled estimates")
    z = (len(v)-1)*p
    i = int(z)
    if i == len(v)-1:
        return v[i]
    return v[i]*(1-(z-i)) + v[i+1]*(z-i)


def bootstrap_summary(values: list[float]) -> dict:
    if not values:
        return {"identified": False}
    return {
        "draws": len(values),
        "q025_technical_not_field_CI": qtile(values, 0.025),
        "q975_technical_not_field_CI": qtile(values, 0.975),
        "fraction_negative": sum(x<0 for x in values)/len(values),
    }


def india_source(source: Path) -> tuple[dict, dict]:
    plots, provenance = parse_xlsx_matrix(source)
    cohort = {
        g: matching_sites(plots, g, INDIA_ZONES) for g in INDIA_GROUPS
    }
    if (len(cohort["station"]) < 10 or len(cohort["between"]) < 10
        or provenance["all_spatial_sites"] != 31
        or provenance["plot_year_zone_rows"] != 248
        or provenance["source_missing_species_cells"] != 4):
        raise ValueError("Original Dryad 31-site / all-zone repeated cohort failed source support gate")
    values = {}
    for group in INDIA_GROUPS:
        values[group] = {}
        for site in cohort[group]:
            values[group][site] = {}
            for year in INDIA_YEARS:
                values[group][site][year] = {}
                rail = plots[site,year,"Track"]
                assert rail is not None
                for ref in INDIA_REFERENCES:
                    refset = plots[site,year,ref]
                    assert refset is not None
                    values[group][site][year][ref] = components(rail,refset)
    return values, {
        "verified_original_Dryad_source": provenance,
        "original_31_locations_source_years": list(INDIA_YEARS),
        "all_four_zones_same_source_sites_retained_both_years": True,
        "station_independent_sites": len(cohort["station"]),
        "between_independent_sites": len(cohort["between"]),
        "source_area_per_zone_m2": 100,
        "station_is_not_train_traffic_exposure": True,
    }


def india_aggregate(data: dict, bootdraws=1999) -> dict:
    cohort = {g: sorted(data[g]) for g in INDIA_GROUPS}
    observed = {}
    for g in INDIA_GROUPS:
        observed[g] = {
            str(year): {
                ref: average([data[g][site][year][ref] for site in cohort[g]])
                for ref in INDIA_REFERENCES
            } for year in INDIA_YEARS
        }
    comparisons = {}
    for year in INDIA_YEARS:
        for ref in INDIA_REFERENCES:
            key = f"{year}:{ref}"
            comparisons[key] = {
                "station_minus_between":
                    {m: observed["station"][str(year)][ref][m] -
                        observed["between"][str(year)][ref][m]
                     if observed["station"][str(year)][ref][m] is not None
                        and observed["between"][str(year)][ref][m] is not None
                     else None
                     for m in BOOT_METRICS}
            }
    temporal = {}
    for ref in INDIA_REFERENCES:
        temporal[ref] = {
            m:(comparisons[f"2017:{ref}"]["station_minus_between"][m] -
               comparisons[f"2014:{ref}"]["station_minus_between"][m])
            if comparisons[f"2017:{ref}"]["station_minus_between"][m] is not None
            and comparisons[f"2014:{ref}"]["station_minus_between"][m] is not None
            else None
            for m in BOOT_METRICS
        }
    # Resample original spatial site as a whole: both years and all
    # comparisons travel together. Site × year is NOT a new site.
    rng = random.Random(SEED)
    boot = {k:{m:[] for m in BOOT_METRICS} for k in comparisons}
    boot_change = {ref:{m:[] for m in BOOT_METRICS} for ref in INDIA_REFERENCES}
    for _ in range(bootdraws):
        drawn={g:rng.choices(cohort[g],k=len(cohort[g])) for g in INDIA_GROUPS}
        averages = {g:{
            str(year): {
                ref:average([data[g][site][year][ref] for site in drawn[g]])
                for ref in INDIA_REFERENCES
            } for year in INDIA_YEARS
        } for g in INDIA_GROUPS}
        diffs={}
        for year in INDIA_YEARS:
            for ref in INDIA_REFERENCES:
                key=f"{year}:{ref}"
                diffs[key]={}
                for m in BOOT_METRICS:
                    a,b=averages["station"][str(year)][ref][m],averages["between"][str(year)][ref][m]
                    if a is None or b is None:
                        continue
                    d=a-b
                    diffs[key][m]=d
                    boot[key][m].append(d)
        for ref in INDIA_REFERENCES:
            for m in BOOT_METRICS:
                a,b=diffs.get(f"2017:{ref}",{}).get(m),diffs.get(f"2014:{ref}",{}).get(m)
                if a is not None and b is not None:
                    boot_change[ref][m].append(a-b)
    return {
        "source_unit": "actual spatial site repeated across original 2014 and 2017, all four 100m2 zones",
        "station_vs_between": observed,
        "year_by_reference_station_minus_between": {
            key: {**v, "bootstrap":{m:bootstrap_summary(x) for m,x in boot[key].items()}}
            for key,v in comparisons.items()
        },
        "paired_site_year_change_in_station_minus_between": {
            ref: {"observed":temporal[ref],
                  "bootstrap":{m:bootstrap_summary(x) for m,x in boot_change[ref].items()}}
            for ref in INDIA_REFERENCES
        },
        "bootstrap_unit":"original location (whole 2-year × 4-zone histories), not pair/site-year",
        "bootstrap_draws_requested":bootdraws,
        "bootstrap_percentiles_NOT_randomized_field_CIs":True,
    }


def slovenia_source(original_xlsx: Path, original_audit: Path) -> tuple[dict, dict]:
    source = extract_original_xlsx(original_xlsx,original_audit)
    loc, zones, species, traits, codebook, provenance = validate_original_keys(source)
    grouped=defaultdict(list)
    side_pairs=0
    for (lid,zid),zone in zones.items():
        if zone["area"]!="BS":
            continue
        if zid not in ("LBS","RBS"):
            raise ValueError("Published BS original zone-ID not in LBS/RBS")
        toe=zones.get((lid,zid[0]+"BT"))
        if toe is None:
            continue
        if toe["area"]!="BT" or zone.get("area_side") != toe.get("area_side"):
            raise ValueError("Source same-side matching failed")
        a=species.get((lid,zid),set())
        b=species.get((lid,zid[0]+"BT"),set())
        grouped[lid].append(components(a,b))
        side_pairs+=1
    if side_pairs!=173 or len(grouped)!=89 or len({str(loc[lid]["track_section_number"]) for lid in grouped})!=15:
        raise ValueError("Original Slovenia 89 locations / 173 original paired track sides / 15 sections unsupported")
    sections=defaultdict(list)
    for lid,side_records in grouped.items():
        section=str(loc[lid]["track_section_number"])
        sections[section].append(average(side_records))
    return sections,{
        "source": "Slovenia Šabić et al. 2026 original publisher supplement",
        "doi": "10.1007/s44473-026-00155-z",
        "original_source_qc":provenance,
        "original_locations":len(grouped),
        "same_side_pairs":side_pairs,
        "original_line_sections":len(sections),
        "original_segment_location_counts":sorted(len(v) for v in sections.values()),
        "within_location_left_right_sides_are_not_independent_field_sites":True,
        "BS_BT_original_survey_area_equivalence_not_verified":True,
    }


def slovenia_aggregate(sections: dict, bootdraws=1999) -> dict:
    sec=sorted(sections)
    bysection={s:average(sections[s]) for s in sec}
    observed=average(list(bysection.values()))
    rng=random.Random(SEED+761)
    draws={m:[] for m in METRICS}
    for _ in range(bootdraws):
        sample=rng.choices(sec,k=len(sec))
        b=average([bysection[s] for s in sample])
        for m in METRICS:
            if b[m] is not None:
                draws[m].append(b[m])
    loo={m:[] for m in METRICS}
    for excluded in sec:
        reduced=average([bysection[s] for s in sec if s!=excluded])
        for m in METRICS:
            if reduced[m] is not None:
                loo[m].append(reduced[m])
    return {
        "source_unit": "same-side BS and BT -> equal location -> equal original railway section (15)",
        "section_equal_observed":observed,
        "section_bootstrap":{m:bootstrap_summary(x) for m,x in draws.items()},
        "leave_one_section_out_min_max":{
            m:{"min":min(v),"max":max(v)} if v else None
            for m,v in loo.items()
        },
        "bootstrap_source_section_count":len(sec),
        "bootstrap_percentiles_NOT_randomized_field_CIs":True,
        "no_original_plot_area_identification":True,
    }


def run(kashmir_source: Path, slov_source: Path, slov_audit: Path,
        bootdraws: int =1999) -> dict:
    if bootdraws<100:
        raise ValueError("Field unit bootstrap requires >=100 technical draws")
    india,india_qc=india_source(kashmir_source)
    slov,slov_qc=slovenia_source(slov_source,slov_audit)
    india_out=india_aggregate(india,bootdraws)
    slov_out=slovenia_aggregate(slov,bootdraws)
    return {
        "version":"v25.6",
        "status":"TWO_AUTHENTIC_RAILWAY_STUDIES_SOURCE_VERIFIED_IF_WORKFLOW_PASSES",
        "source_studies":2,
        "models":"Baselga (2012) pairwise Jaccard replacement (jtu) and nestedness-resultant (jne)",
        "method_reference_DOI":"10.1111/j.1466-8238.2011.00756.x",
        "source_audits":{"Kashmir":india_qc,"Slovenia":slov_qc},
        "Kashmir_2014_2017":india_out,
        "Slovenia_2022_2023":slov_out,
        "cross_source_meta_effect_not_computed":True,
        "no_causal_train_dispersal_or_mowing_or_evolved_trait_claim":True,
        "Kashmir_true_equal_100m2_vs_Slovenia_BS_BT_area_unverified":True,
        "within_study_site_and_section_weights_not_to_be_read_as_two_independent_country_replicates":True,
        "per_plot_taxa_GPS_not_redistributed":True,
    }


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-source",required=True,type=Path)
    p.add_argument("--slovenia-xlsx",required=True,type=Path)
    p.add_argument("--slovenia-audit",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--bootstrap",default=1999,type=int)
    args=p.parse_args()
    result=run(args.kashmir_source,args.slovenia_xlsx,
               args.slovenia_audit,args.bootstrap)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                           encoding="utf8")
    print(json.dumps({
        "version":result["version"],"status":result["status"],
        "Kashmir_support":result["source_audits"]["Kashmir"],
        "Slovenia_support":result["source_audits"]["Slovenia"]["original_locations"],
        "Slovenia_equal_section_beta":result["Slovenia_2022_2023"]["section_equal_observed"]["beta_jac"],
        "Kashmir_track_ref2014":result["Kashmir_2014_2017"]["station_vs_between"]["station"]["2014"]["Interior Plot 2"]["beta_jac"]
    },ensure_ascii=False))


if __name__=="__main__":
    main()

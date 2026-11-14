#!/usr/bin/env python3
"""v10.4: rail microhabitat WIND DISPERSAL TRAIT vs ANNUAL LIFE FORM
cross-sectional field observation, exact published Slovenia 2026 source.

Does the within-location same-side Ballast Shoulder (BS) vs Ballast
Toe (BT) share of wind-dispersed flora shift in SAME direction as annual
and therophyte share? An independent dispersion-vs-establishment-filter
diagnostic, not actual measured seed carriage or wind velocity.

Source Šabić, Grašič & Rovanšek (2026), DOI 10.1007/s44473-026-00155-z.
Original supplementary 1 XLSX checksum fixed and validated by v10.2,
original source 92 RAIL sites/482 field zones/6331 taxon records
(6295 unique zone×EPPO incidence after 36 source duplications)/312 traits.
Use original 'Diaspore.strategy' EXACT codes:
  WIND ('Epilobium' + 'Phragmites') as used by authors;
  additional 'Lycopodium' spores an EXPLICIT unverified sensitivity,
  'Allium' predominantly autochory, 'Bidens' epizoochory, 'Cornus'
  endozoochory as separate source-category comparators.

Missing Diaspore.strategy is UNKNOWN, excluded from trait-known
denominator, NEVER treated as nonwind. Source trait species-level mean
is not actual local plasticity; zone area differences forbid a direct
species richness per 100m2 comparison. Both track sides average
within each one of 89 actual railway locations; source railway segments
15; only 11 have 3+ independent field locations for section sign.

No observed railway slipstream speed, wind seed movement, freight
propagule release, mowing, herbicide-before-after, trait selection
experiment, heritable evolution or counterfactual road sampling.
"""
from __future__ import annotations
import argparse,json,math
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
from statistics import mean,median

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_SHA256,SOURCE_DOI,SOURCE_BYTES,
    extract_original_xlsx,validate_original_keys,pair_source,
    locations_independent,location_summary,section_summary)

WIND_STRA="wind_Epilobium_Phragmites"
WIND_SENS="wind_plus_spore_Lycopodium"
GROUPS={
    WIND_STRA:frozenset(("Epilobium","Phragmites")),
    WIND_SENS:frozenset(("Epilobium","Phragmites","Lycopodium")),
    "source_Allium_autochory":frozenset(("Allium",)),
    "source_Bidens_epizoochory":frozenset(("Bidens",)),
    "source_Cornus_endozoochory":frozenset(("Cornus",)),
}
VALID_SOURCE_CATEGORIES=frozenset((
    "Allium","Epilobium","Cornus","Bidens","Sparganium",
    "Phragmites","Lycopodium","Zea"))
MIN_SEGMENT_FIELD_LOCATIONS=3
MIN_LOCATION_SOURCE_SAME_SIDE_PAIRS=1


def strategy_fraction(codes,traits,strategies):
    if not set(strategies).issubset(VALID_SOURCE_CATEGORIES):
        raise ValueError("Not a declared original Slovenian Diaspore.strategy code")
    vals=[]
    for code in codes:
        source=traits.get(code)
        if source is None:
            raise ValueError("Original EPPO code in railway species source lacking trait metadata")
        category=source.get("Diaspore.strategy")
        if category is None or category=="":
            continue
        if category not in VALID_SOURCE_CATEGORIES:
            raise ValueError("Unexpected 2026 publisher Diaspore.strategy category")
        vals.append(category)
    return {
        "n_taxa_with_source_known_Diaspore_strategy":len(vals),
        "n_taxa_missing_original_Diaspore_strategy":len(codes)-len(vals),
        "fraction_of_source_trait_known_taxa":sum(v in strategies for v in vals)/len(vals)
                 if vals else None,
        "NO_direct_field_seed_movement_or_mowing_wind_measurement":True,
    }


def measure(pairs,sp,traits):
    out=[]
    for item in pairs:
        lid=item["location_id"]; side=item["track_side"]
        src={}
        for micro,zone_id in (("BS",side+"BS"),("BT",side+"BT")):
            taxa=sp.get((lid,zone_id))
            if taxa is None:
                raise ValueError("Original same-side BS or BT field taxon sources missing")
            src[micro]={name:strategy_fraction(taxa,traits,strats)
                        for name,strats in GROUPS.items()}
        differences={}
        for group in GROUPS:
            b=src["BS"][group]["fraction_of_source_trait_known_taxa"]
            t=src["BT"][group]["fraction_of_source_trait_known_taxa"]
            differences[group]=b-t if b is not None and t is not None else None
        out.append({
            "location_id":lid,
            "section_number":item["section_number"],
            "track_type":item["track_type"],
            "survey_year":item["survey_year"],
            "track_side":side,
            "source_strict_annual_BS_minus_BT": (
                item["BS"].get("strict_annual_share_of_known_lifespans")-
                item["BT"].get("strict_annual_share_of_known_lifespans")
                if item["BS"].get("strict_annual_share_of_known_lifespans") is not None
                and item["BT"].get("strict_annual_share_of_known_lifespans") is not None
                else None),
            "source_strict_therophyte_BS_minus_BT": (
                item["BS"].get("strict_therophyte_share_of_known_lifeforms")-
                item["BT"].get("strict_therophyte_share_of_known_lifeforms")
                if item["BS"].get("strict_therophyte_share_of_known_lifeforms") is not None
                and item["BT"].get("strict_therophyte_share_of_known_lifeforms") is not None
                else None),
            "source_cover_BS_minus_BT":(
                item["BS_cover_source_numeric_not_guaranteed_percentage"]-
                item["BT_cover_source_numeric_not_guaranteed_percentage"]
                if item["BS_cover_source_numeric_not_guaranteed_percentage"] is not None
                and item["BT_cover_source_numeric_not_guaranteed_percentage"] is not None
                else None),
            "original_source_cover_GT100":item["source_cover_value_above_100_on_either_side"],
            "wind_and_diaspore_group_shares_BS_minus_BT":differences,
            "n_original_taxa_known_strategy_BS":src["BS"][WIND_STRA][
                "n_taxa_with_source_known_Diaspore_strategy"],
            "n_original_taxa_known_strategy_BT":src["BT"][WIND_STRA][
                "n_taxa_with_source_known_Diaspore_strategy"],
            "n_original_taxa_missing_strategy_BS":src["BS"][WIND_STRA][
                "n_taxa_missing_original_Diaspore_strategy"],
            "n_original_taxa_missing_strategy_BT":src["BT"][WIND_STRA][
                "n_taxa_missing_original_Diaspore_strategy"],
        })
    if len(out)!=173 or len({v["location_id"] for v in out})!=89:
        raise ValueError("Same true 89 railway sites / 173 side pairs source changed")
    return out


def location_level(rows):
    by=defaultdict(list)
    for p in rows:by[p["location_id"]].append(p)
    loc=[]
    for lid,items in sorted(by.items()):
        first=items[0]
        if (any(r["section_number"]!=first["section_number"] or
                r["survey_year"]!=first["survey_year"] or
                r["track_type"]!=first["track_type"] for r in items)
            or len(items) not in (1,2)):
            raise ValueError("Different source RAIL site sides not nested within same location")
        r={"location_id":lid,"section_number":first["section_number"],
           "survey_year":first["survey_year"],
           "track_type":first["track_type"],
           "n_original_BS_BT_source_track_side_pairs":len(items),
           "source_cover_GT100":any(z["original_source_cover_GT100"] for z in items)}
        for group in GROUPS:
            x=[p["wind_and_diaspore_group_shares_BS_minus_BT"][group] for p in items]
            v=[q for q in x if q is not None]
            r[group+"_delta"]=mean(v) if v else None
        for metric in ("source_strict_annual_BS_minus_BT","source_strict_therophyte_BS_minus_BT",
                       "source_cover_BS_minus_BT"):
            x=[p[metric] for p in items]
            v=[q for q in x if q is not None]
            r[metric]=mean(v) if v else None
        loc.append(r)
    if len(loc)!=89 or len({r["section_number"] for r in loc})!=15:
        raise ValueError("Actual published RAIL source field location / section count changed")
    return loc


def summarize_locations(rows,measure,report_parts=False):
    subset=[r for r in rows if r[measure] is not None]
    x=[r[measure] for r in subset]
    out={
        "n_independent_actual_RAIL_field_locations":len(subset),
        "mean_BS_minus_BT_fraction":round(mean(x),6) if x else None,
        "median_BS_minus_BT_fraction":round(median(x),6) if x else None,
        "n_locations_BS_larger":sum(v>0 for v in x),
        "n_locations_BT_larger":sum(v<0 for v in x),
        "n_locations_equal":sum(v==0 for v in x),
    }
    by=defaultdict(list)
    for r in subset:by[r["section_number"]].append(r[measure])
    sign=[mean(v) for v in by.values() if len(v)>=MIN_SEGMENT_FIELD_LOCATIONS]
    out["n_source_RAIL_segments_with_3plus_independent_locations"]=len(sign)
    out["n_segments_BS_larger"]=sum(v>0 for v in sign)
    out["n_segments_BT_larger"]=sum(v<0 for v in sign)
    out["median_segment_equal_mean_BS_minus_BT"]=(round(median(sign),6)
                                                  if sign else None)
    out["region_codes_are_rail_segments_not_independent_countries"]=True
    if report_parts:
        out["segments"]=[{"track_section_number":str(k),"n_field_locations":len(v),
                           "mean_BS_minus_BT":round(mean(v),6)}
                         for k,v in sorted(by.items(),key=lambda t:str(t[0]))
                         if len(v)>=MIN_SEGMENT_FIELD_LOCATIONS]
    return out


def joint_direction(loc):
    both=[v for v in loc if v[WIND_STRA+"_delta"] is not None
          and v["source_strict_annual_BS_minus_BT"] is not None]
    counter={
        "BS_annual_up_BS_wind_down":0,
        "BS_annual_up_BS_wind_up":0,
        "BS_annual_down_BS_wind_down":0,
        "BS_annual_down_BS_wind_up":0,
        "zero_or_tied_trait_differences":0,
    }
    for v in both:
        a=v["source_strict_annual_BS_minus_BT"]
        w=v[WIND_STRA+"_delta"]
        if a==0 or w==0:
            counter["zero_or_tied_trait_differences"]+=1
            continue
        counter[f"BS_annual_{'up' if a>0 else 'down'}_BS_wind_{'up' if w>0 else 'down'}"]+=1
    return {"n_original_locations_with_valid_both_trait_source_differences":len(both),
            "contrasting_trait_directions":counter,
            "source_species_lifespan_does_NOT_measure_selection_or_gust_causality":True}


def run(supplement,previous):
    if (previous.get("complete") is not True or
        previous.get("n_same_track_side_pairs")!=173 or
        previous.get("n_distinct_actual_railway_field_locations")!=89):
        raise ValueError("Independently verified RAIL v10.3 site source not available")
    loc,zones,sp,traits,cb,source=validate_original_keys(supplement)
    original_pairs=pair_source(loc,zones,sp,traits)
    paired=measure(original_pairs,sp,traits)
    location=location_level(paired)
    out={
        "complete":True,
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_source_DOI":SOURCE_DOI,
        "original_supplement_sha256":SOURCE_SHA256,
        "n_original_source_rail_field_locations":89,
        "n_original_BS_vs_BT_same_side_pairs":173,
        "n_original_line_segments":15,
        "n_original_tracked_species_trait_keys":len(traits),
        "source_original_2022_2023_cross_sectional_RAIL_not_ROAD":True,
        "known_source_diaspore_categories":sorted(VALID_SOURCE_CATEGORIES),
        "source_wind_classes_explicitly_Epilobium_and_Phragmites":True,
        "observed_within_location_BS_minus_BT_group_differences":{},
        "independent_rail_segment_field_location_checks":{},
        "within_same_89_field_locations_joint_annual_vs_wind_direction":joint_direction(location),
        "mandatory_source_track_and_year_sensitivity":{},
        "original_EPPO_codebook_discrepancy_not_silently_reconciled":source[
            "trait_EPPO_codes_absent_from_species_codebook"],
        "NO_measured_train_slipstream_wind_or_seed_delivery_mowing_adaptation":True,
    }
    for label in GROUPS:
        col=label+"_delta"
        out["observed_within_location_BS_minus_BT_group_differences"][label]=summarize_locations(
            location,col)
        out["independent_rail_segment_field_location_checks"][label]=summarize_locations(
            location,col,report_parts=True)["segments"]
    for yr in (2022,2023):
        source_l=[r for r in location if r["survey_year"]==yr]
        out["mandatory_source_track_and_year_sensitivity"]["year_"+str(yr)]={
            label:summarize_locations(source_l,label+"_delta") for label in GROUPS}
    for kind in ("SINGLE","DOUBLE"):
        source_l=[r for r in location if r["track_type"]==kind]
        out["mandatory_source_track_and_year_sensitivity"]["track_"+kind]={
            label:summarize_locations(source_l,label+"_delta") for label in GROUPS}
    keep=[r for r in location if not r["source_cover_GT100"]]
    out["mandatory_source_track_and_year_sensitivity"]["remove_original_cover_over100"]={
        label:summarize_locations(keep,label+"_delta") for label in GROUPS}
    out["source_integrity"]={
        "original_482_zone_richness_conserves_6295_unique_EPPO_zone_incidence":source[
            "all_482_original_zone_richness_equals_unique_species_codes"],
        "original_taxon_species_zone_rows":source["n_original_species_zone_occurrence_rows"],
        "source_EPPO_codebook_trait_discrepancy":source[
            "trait_EPPO_codes_absent_from_species_codebook"],
        "no_nonfield_plant_species_added":True,
    }
    out["limitations"]=[
        "Wind-dispersal categories Epilobium and Phragmites are original plant species-level classification, not observed railway-train wind, seed movement or seed viability.",
        "Higher annual-plant fraction in ballast shoulder despite lower relative wind-dispersed taxa suggests DIFFERENT potential functional axes but is not a factorial independent selection experiment.",
        "Only 89 field locations from 15 nonrandom railway segments in one country; left/right sides averaged within location, and 11 segments have at least 3 field-site observations.",
        "Published BS and BT zone dimensions/true areas differ or are uncertain; evaluate fractions among known strategy-coded taxa, not zone species counts as equivalent quadrats.",
        "Missing Diaspore.strategy classes are excluded from strategy denominator and cannot be classified as non-wind or native.",
        "2022–2023 field survey post-dates 2021 herbicide policy change; no before-policy baseline, site mowing frequency, measured gust velocity, or controlled train traffic exposure.",
        "A wind-diaspore strategy permits natural ambient wind dispersal, independent of moving trains; railroad microhabitat is not proof seeds were delivered there by freight.",
        "No heritability, genetic change, common-garden, plant size regrowth plasticity, actual seeds on train cars or live seed traps.",
        "Original dataset belongs to Slovenia study authors and is supplementary material, not original tetsu field collection.",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","verified103","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    prior=json.loads(Path(a.verified103).read_text(encoding="utf-8"))
    original=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    result=run(original,prior)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "source_original_RAIL_locations":result["n_original_source_rail_field_locations"],
        "paired_BS_BT_sides":result["n_original_BS_vs_BT_same_side_pairs"],
        "wind_vs_annual_direction":result["within_same_89_field_locations_joint_annual_vs_wind_direction"],
        "original_wind":result["observed_within_location_BS_minus_BT_group_differences"][WIND_STRA],
        "wind_plus_spore_sensitivity":result["observed_within_location_BS_minus_BT_group_differences"][WIND_SENS],
        "rail_segment_signed":result["independent_rail_segment_field_location_checks"][WIND_STRA],
        "no_observed_train_wind_or_seed_delivery_causation":True
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

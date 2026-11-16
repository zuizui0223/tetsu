#!/usr/bin/env python3
"""v7.8: source-only falsification diagnostic of reported TRAINS vs plants.

Exactly the two pre-frozen US 50km blocks from v7.6, SAME eligible cells
and SAME GBIF/freight-crossing data (not independently sampled rail verge
plots). Compare GRIIS woody/herbaceous richness against ALL-PLANT
recorded speciesKey richness, non-GRIIS taxon-category richness (UNKNOWN
native/alien status), and the woody share of GRIIS recorded taxa.
An additional CLOSED-crossing fraction infrastructure comparator tests
whether track/road geography also covaries with a botanical response.
NO actual train-borne seeds, active freight volume, stress or adaptation.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
from statistics import median

from analyze_sixcountry_functional_rail_v71 import validate_traits
from analyze_US_train_counts_GRIIS_traits_v76 import (
    BLOCKS,PERIOD,K,GROUPS,MIN_TRAIN_CROSSINGS_PER_CELL,
    MIN_GIFT_TRAIT_COVERAGE,COVARIATES,split_counts,species_metrics,
    partial_rank_cells,build_block)

CONTROL_OUTCOMES=(
    "all_Plantae_observed_speciesKey_categories",
    "non_GRIIS_speciesKey_categories_UNKNOWN_nativity",
    "country_GRIIS_observed_woody_fraction_of_known_growth_forms",
    "country_GRIIS_recorded_speciesKey_categories",
    "woody_GRIIS_rarefied20_expected_richness",
    "herbaceous_GRIIS_rarefied20_expected_richness")
EXPOSURES=("real_reported_train_count","closed_crossing_fraction_infrastructure_comparator")
MIN_CELLS=10

def source_rho(rows,exposure,response,control_fields):
    """Re-use unchanged v7.6 exact partial-rank model and gates."""
    key={
       "real_reported_train_count":"reported_train_exposure_log1p_median_day_night_through_trains",
       "closed_crossing_fraction_infrastructure_comparator":"_closed_crossing_fraction",
    }[exposure]
    # Adapt only predictor key; all source rank/covariates and cell gates
    # are identical to the previously independently verified v7.6.
    adjusted=[]
    for row in rows:
        a=dict(row)
        a["reported_train_exposure_log1p_median_day_night_through_trains"]=row[key]
        adjusted.append(a)
    return partial_rank_cells(adjusted,response,control_fields)

def source_rows(block,flora,fra,traitmap):
    if block not in BLOCKS or flora.get("country")!="US":
        raise ValueError("Only original frozen two US FRA-GRIIS macro blocks")
    if fra.get("source_open_closed_balance_verified") is not True:
        raise ValueError("Closed FRA crossings must not contribute train intensity")
    cells=flora.get("cells",{})
    crossings=fra.get("by_original_50km_cells",{})
    if len(cells)!=16 or set(cells)!=set(crossings):
        raise ValueError("Source 16 x 50km GBIF/FRA cell geography mismatch")
    rows=[]
    exclusions=Counter()
    for cid,old in cells.items():
        if old.get("source_pre_outcome_both_period_recording_eligible") is not True:
            exclusions["GBIF_two_period_original_gate"]+=1
            continue
        p=old.get("periods",{}).get(PERIOD,{})
        if p.get("full_facet_source_pagination_completed") is not True or (
            p.get("alien_observed_minimum_5_species_and_20_records") is not True):
            exclusions["GRIIS_facet_min_source_opportunity"]+=1
            continue
        crossing=crossings[cid]
        if crossing.get("train_opportunity_metrics_USE_confirmed_OPEN_crossings_ONLY") is not True:
            raise ValueError("FRA v7.5 current traffic source may include closed crossings")
        nobs=crossing.get("n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields")
        nopen=crossing.get("n_crossings_confirmed_open_in_current_source")
        nclosed=crossing.get("n_crossings_confirmed_closed_excluded_from_operating_exposure")
        nall=crossing.get("n_distinct_crossing_inventory_records")
        train=crossing.get("median_reported_total_daily_through_trains_at_eligible_crossings")
        if type(nobs)!=int or nobs<MIN_TRAIN_CROSSINGS_PER_CELL or train is None:
            exclusions["no_year_eligible_open_crossing_traffic"]+=1
            continue
        if (type(train) not in (int,float) or not math.isfinite(train)
            or train<0 or any(type(n)!=int or n<0 for n in (nopen,nclosed,nall))
            or nall<=0 or nopen+nclosed>nall):
            raise ValueError("Invalid open/closed crossing or train-count source")
        nrec=p.get("indexed_all_plantae_source_records")
        nalltax=p.get("n_all_observed_species_key_categories")
        counts=p.get("observed_alien_specieskey_counts")
        if (type(nrec)!=int or nrec<=0 or type(nalltax)!=int
            or not isinstance(counts,dict) or nalltax<len(counts)):
            raise ValueError("All-Plantae GBIF total taxon and GRIIS source counts inconsistent")
        groups=split_counts(counts,traitmap)
        known=len(groups["herbaceous"])+len(groups["woody"])
        cov=known/len(groups["all_GRIIS_introduced"])
        if cov<MIN_GIFT_TRAIT_COVERAGE:
            exclusions["below_half_GIFT_growth_form_coverage"]+=1
            continue
        coords=cid.split("_x")[-1].split("_y")
        if len(coords)!=2 or any(not v.isdigit() for v in coords):
            raise ValueError("Original EPSG:6933 50km source cell coordinate key invalid")
        wp=species_metrics(groups["woody"])
        hp=species_metrics(groups["herbaceous"])
        responses={
            "all_Plantae_observed_speciesKey_categories":nalltax,
            "non_GRIIS_speciesKey_categories_UNKNOWN_nativity":nalltax-len(counts),
            "country_GRIIS_observed_woody_fraction_of_known_growth_forms":len(groups["woody"])/known,
            "country_GRIIS_recorded_speciesKey_categories":len(counts),
            "woody_GRIIS_rarefied20_expected_richness":wp["rarefied20_expected_species_richness"] if wp else None,
            "herbaceous_GRIIS_rarefied20_expected_richness":hp["rarefied20_expected_species_richness"] if hp else None
        }
        years=crossing.get("reported_train_count_year_histogram",{})
        if not isinstance(years,dict) or any(type(n)!=int or n<0 for n in years.values()):
            raise ValueError("Malformed source crossing train count year histogram")
        rows.append({
            "cell_id":cid,
            "reported_train_exposure_log1p_median_day_night_through_trains":math.log1p(train),
            "_closed_crossing_fraction":nclosed/nall,
            "raw_reported_crossing_median_through_train_count":train,
            "log_all_Plantae_indexed_records":math.log(nrec),
            "log_confirmed_open_crossing_count":math.log1p(nopen),
            "cell_grid_x":int(coords[0]),"cell_grid_y":int(coords[1]),
            "responses":responses,"year_hist":years,
        })
    return rows,dict(exclusions)

def one_block(block,flora,fra,traits,original_v76):
    rows,excluded=source_rows(block,flora,fra,traits)
    # Fail closed if this 'new' test inadvertently changes its source
    # selection or baseline outcome from independently successful v7.6.
    previous=original_v76["regions"][block]
    if len(rows)!=previous["n_eligible_cells_for_source_observations"]:
        raise ValueError("Matched source cell cohort has drifted from frozen v7.6")
    original_woody=previous["descriptive_partial_rank_results"]["woody"][
        "rarefied20_expected_species_richness"]["gbif_and_crossing_opportunity"]["rho"]
    target="woody_GRIIS_rarefied20_expected_richness"
    verify=source_rho(rows,"real_reported_train_count",target,
                      COVARIATES["gbif_and_crossing_opportunity"])["rho"]
    if original_woody!=verify:
        raise ValueError("Original real v7.6 country woody coefficient did not replicate exactly")
    results={}
    for exposure in EXPOSURES:
        results[exposure]={}
        for response in CONTROL_OUTCOMES:
            results[exposure][response]={
                name:source_rho(rows,exposure,response,fields)
                for name,fields in COVARIATES.items()
            }
    positive=[r for r in rows if r["raw_reported_crossing_median_through_train_count"]>0]
    positive_results={
        response:source_rho(positive,"real_reported_train_count",response,
                            COVARIATES["gbif_and_crossing_opportunity"])
        for response in CONTROL_OUTCOMES}
    years=Counter()
    for row in rows:years.update(row["year_hist"])
    return {
        "n_same_eligible_50km_source_cells_as_v76":len(rows),
        "n_positive_reported_train_median_cells":len(positive),
        "excluded_original_source_cells_by_reason":excluded,
        "reproduced_original_v76_woody_rho":verify,
        "control_outcomes_and_infrastructure_comparator":results,
        "positive_through_train_only_same_source_cell_sensitivity":positive_results,
        "open_crossing_train_reported_source_year_category_counts_over_eligible_cells":dict(years),
        "n_cells_with_any_source_count_reported_2026":sum(
            r["year_hist"].get("2026",0)>0 for r in rows),
        "NO_inference_on_true_native_species_adaptation_or_seed_flux":True,
    }

def run(flora,fra,gift,original):
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("Original GBIF GRIIS source incomplete")
    if fra.get("complete") is not True or original.get("complete") is not True:
        raise ValueError("Previously audited FRA/GIFT 2-block source missing")
    if original.get("frozen_regions")!=list(BLOCKS):
        raise ValueError("Original source comparison country blocks changed")
    traits=validate_traits(gift)
    out={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "study":"FRA traffic versus GRIIS alien plants with all-plant recording/outcome comparator and CLOSED crossing infrastructure comparator",
        "frozen_50km_US_blocks":list(BLOCKS),
        "observational_period":"GBIF 2018_2025 vs mixed FRA 2018_2026 reported train count years",
        "country_macro_region_not_cell_as_independent_replicate":True,
        "not_plant_evolution_mowing_wind_or_rail_carried_seeds":True,
        "regions":{}}
    for block in BLOCKS:
        out["regions"][block]=one_block(block,flora["regions"][block],
                                       fra["regions"][block],traits,original)
        row=out["regions"][block]
        main=row["control_outcomes_and_infrastructure_comparator"]["real_reported_train_count"]
        print(json.dumps({"block":block,"n_cells":row["n_same_eligible_50km_source_cells_as_v76"],
            "negative_control_all_Plantae_rho":main["all_Plantae_observed_speciesKey_categories"]
                 ["gbif_and_crossing_opportunity"]["rho"],
            "woody_GRIIS_rho":row["reproduced_original_v76_woody_rho"],
            "woody_GRIIS_share_rho":main["country_GRIIS_observed_woody_fraction_of_known_growth_forms"]
                 ["gbif_and_crossing_opportunity"]["rho"]},ensure_ascii=False),flush=True)
    out["limits"]=[
        "The all-Plantae recorded-species response is a botanical observation comparator, not actual native-plant richness; GRIIS non-listed taxa have UNKNOWN naturalization status.",
        "Closed crossing density is an infrastructure/history comparator, NOT a randomized or validated negative-control exposure.",
        "All source comparisons remain observational in only two selected US geographic macro blocks.",
        "Reported FRA 2026 snapshot includes 2018-2026 train-count source years, does not measure 2018-2025 simultaneous traffic, and includes no freight-specific volume.",
        "Only 14/13 actual independent 50km source cells within two adjacent/selected geographic windows; no spatially valid significance tests.",
        "GIFT growth form is not measured seed mass, aerodynamic dispersal, mowing tolerance, gust wind or adaptation.",
        "Source record effort adjustment and fractions cannot remove GBIF site/dataset/observer/urban/ecoregion biases."
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for name in ("flora","fra","gift","baseline","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args();read=lambda x:json.loads(Path(x).read_text(encoding="utf-8"))
    out=run(read(a.flora),read(a.fra),read(a.gift),read(a.baseline))
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":True,"blocks":list(out["regions"]),
       "no_train_causal_mechanism_claim":True},ensure_ascii=False),flush=True)

if __name__=="__main__":main()

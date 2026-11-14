#!/usr/bin/env python3
"""v7.6: two PRESELECTED US 50km macro regions, FRA crossing 2026 open-line
train COUNT opportunity versus GBIF 2018–25 GRIIS-listed introduced flora
observational richness and GIFT herb/woody classification.

SOURCE-ONLY descriptive 2-block ecological association, not observed
freight movement, train-carried propagules, railway verge vegetation,
wind/mowing exposure, species adaptation or historical change.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import numpy as np
from analyze_macro_50km_rail_road_aliens_v34 import _rank
from analyze_macro_rail_assemblages_v22 import incidence_probabilities
from analyze_sixcountry_functional_rail_v71 import validate_traits

BLOCKS=("us_midwest","us_southeast")
PERIOD="2018_2025"
GROUPS=("all_GRIIS_introduced","herbaceous","woody")
MEASURES=("raw_recorded_species_richness","rarefied20_expected_species_richness")
COVARIATES={
  "gbif_and_crossing_opportunity":("log_all_Plantae_indexed_records","log_confirmed_open_crossing_count"),
  "gbif_only":("log_all_Plantae_indexed_records",),
  "crossing_only":("log_confirmed_open_crossing_count",),
  "plus_spatial_grid_position":("log_all_Plantae_indexed_records","log_confirmed_open_crossing_count","cell_grid_x","cell_grid_y")
}
MIN_COUNTRY_CELLS=10
MIN_TRAIN_CROSSINGS_PER_CELL=5
K=20
MIN_GROUP_SPECIES=3
MIN_GIFT_TRAIT_COVERAGE=.5

def partial_rank_cells(rows,response,controls):
    usable=[r for r in rows if r["responses"].get(response) is not None]
    if len(usable)<MIN_COUNTRY_CELLS:
        return {"rho":None,"n_50km_cells":len(usable),"support":False}
    x=np.asarray([r["reported_train_exposure_log1p_median_day_night_through_trains"] for r in usable],float)
    y=np.asarray([r["responses"][response] for r in usable],float)
    z=np.asarray([[r[c] for c in controls] for r in usable],float)
    if (len(set(x))<5 or len(set(y))<5 or not np.all(np.isfinite(x))
            or not np.all(np.isfinite(y)) or not np.all(np.isfinite(z))):
        return {"rho":None,"n_50km_cells":len(usable),"support":False}
    A=np.column_stack((np.ones(len(usable)),*(_rank(z[:,i]) for i in range(z.shape[1]))))
    if np.linalg.matrix_rank(A)<A.shape[1] or len(usable)-A.shape[1]<3:
        return {"rho":None,"n_50km_cells":len(usable),"support":False}
    rx=_rank(x)-A@np.linalg.lstsq(A,_rank(x),rcond=None)[0]
    ry=_rank(y)-A@np.linalg.lstsq(A,_rank(y),rcond=None)[0]
    norm=float(np.linalg.norm(rx)*np.linalg.norm(ry))
    rho=float(rx@ry/norm) if norm>1e-8 else None
    if rho is not None and not -1.000001<=rho<=1.000001:
        raise ValueError("Invalid partial rank regression")
    return {"rho":round(rho,6) if rho is not None else None,
            "n_50km_cells":len(usable),"support":rho is not None,
            "no_cell_level_p_value_due_to_spatial_autocorrelation":True}

def split_counts(source,traits):
    if any(type(n) is not int or n<=0 for n in source.values()):
        raise ValueError("Invalid GBIF country GRIIS source taxon count")
    return {
        "all_GRIIS_introduced":source,
        "herbaceous":{k:n for k,n in source.items() if traits.get(k)=="herbaceous"},
        "woody":{k:n for k,n in source.items() if traits.get(k)=="woody"},
    }

def species_metrics(group):
    if len(group)<MIN_GROUP_SPECIES or sum(group.values())<K:
        return None
    return {
        "raw_recorded_species_richness":len(group),
        "rarefied20_expected_species_richness":sum(incidence_probabilities(group,k=K).values())
    }

def leave_one_cell_out_diagnostic(rows,response,controls):
    """Posthoc source sensitivity. NOT spatial-block independent inference.

    Remove every original cell once from the already matched source rows,
    retaining the SAME response/controls and immutable exposure selection.
    """
    original=partial_rank_cells(rows,response,controls)
    estimates={}
    for item in sorted(rows,key=lambda v:v["cell_id"]):
        cid=item["cell_id"]
        value=partial_rank_cells([r for r in rows if r["cell_id"]!=cid],response,controls)
        estimates[cid]=value["rho"]
    valid=[v for v in estimates.values() if v is not None]
    base=original["rho"]
    return {
        "complete_n_withheld_cells":len(estimates),
        "n_eligible_withheld_50km_cells":len(valid),
        "n_unestimable_withheld_cells":len(estimates)-len(valid),
        "baseline_partial_rho":base,
        "min_and_max_partial_rho_excluding_one_cell":[min(valid),max(valid)] if valid else None,
        "n_negative_rho_jackknife":sum(v<0 for v in valid),
        "n_positive_rho_jackknife":sum(v>0 for v in valid),
        "same_sign_for_all_estimable_leaveone_cellout":bool(valid) and base is not None and all(
            (v<0)==(base<0) for v in valid if v!=0) and all(v!=0 for v in valid),
        "no_country_replication_or_statistical_significance_from_withheld_cells":True,
    }

def build_block(block,flora,fra,traits):
    if block not in BLOCKS or flora["country"]!="US":
        raise ValueError("Pre-frozen US botanical source block required")
    if (not isinstance(fra.get("by_original_50km_cells"),dict) or
            fra.get("source_open_closed_balance_verified") is not True):
        raise ValueError("Current FRA train crossing source excludes closed lines incorrectly")
    cells=flora["cells"]
    op=fra["by_original_50km_cells"]
    if set(cells)!=set(op) or len(cells)!=16:
        raise ValueError("FRA and GRIIS same 50km source cell identifiers required")
    rows=[]
    excluded={}
    for cid,src in cells.items():
        if src.get("source_pre_outcome_both_period_recording_eligible") is not True:
            excluded[cid]="original_GBIF_not_source_eligible"
            continue
        p=src.get("periods",{}).get(PERIOD)
        if not p or p.get("full_facet_source_pagination_completed") is not True or (
                p.get("alien_observed_minimum_5_species_and_20_records") is not True):
            excluded[cid]="GBIF_positive_GRIIS_source_unqualified"
            continue
        ft=op[cid]
        if ft.get("train_opportunity_metrics_USE_confirmed_OPEN_crossings_ONLY") is not True:
            raise ValueError("Actual crossing source must exclude CLOSED sites from traffic intensity")
        nc=ft["n_crossings_with_2018to2026_train_count_year_and_both_daily_count_fields"]
        train=ft["median_reported_total_daily_through_trains_at_eligible_crossings"]
        if nc<MIN_TRAIN_CROSSINGS_PER_CELL or train is None:
            excluded[cid]="no_adequate_reported_train_count"
            continue
        if type(train) not in (int,float) or not math.isfinite(train) or train<0:
            raise ValueError("FRA source measured median train count invalid")
        total=p["indexed_all_plantae_source_records"]
        if type(total)!=int or total<=0:
            raise ValueError("Missing independent GRIIS recording source effort")
        groups=split_counts(p["observed_alien_specieskey_counts"],traits)
        mapped=set(groups["herbaceous"])|set(groups["woody"])
        cov=len(mapped)/len(groups["all_GRIIS_introduced"])
        if cov<MIN_GIFT_TRAIT_COVERAGE:
            excluded[cid]="growth_form_source_coverage_below_half"
            continue
        grid=cid.split("_x")[-1]
        try:
            ix,iy=(int(q) for q in grid.split("_y"))
        except (ValueError,IndexError):
            raise ValueError("Original 50km source cell grid identity not recognized")
        response={}
        for g in GROUPS:
            stats=species_metrics(groups[g])
            for metric in MEASURES:
                response[f"{g}__{metric}"]=stats[metric] if stats else None
        rows.append({
            "cell_id":cid,
            "reported_train_exposure_log1p_median_day_night_through_trains":math.log1p(train),
            "raw_reported_crossing_median_through_train_count":train,
            "n_open_crossings_with_year_eligible_2018to2026_count":nc,
            "log_all_Plantae_indexed_records":math.log(total),
            "log_confirmed_open_crossing_count":math.log1p(ft["n_crossings_confirmed_open_in_current_source"]),
            "cell_grid_x":ix,"cell_grid_y":iy,
            "source_growth_form_taxon_match_fraction":round(cov,6),
            "responses":response})
    results={}
    for group in GROUPS:
        results[group]={}
        for metric in MEASURES:
            name=f"{group}__{metric}"
            results[group][metric]={
              controls:partial_rank_cells(rows,name,fields)
              for controls,fields in COVARIATES.items()
            }
    # A TRUE positive-through-train sensitivity is scientifically critical:
    # a presently OPEN crossing can nevertheless report 0 trains.
    # Exclude zero-train median cells BEFORE estimating sign; never relax
    # the prespecified >=10 distinct 50km cells when it removes a region.
    nonzero_rows=[r for r in rows if r["raw_reported_crossing_median_through_train_count"]>0]
    positive_only={}
    for group in GROUPS:
        positive_only[group]={}
        for metric in MEASURES:
            response=f"{group}__{metric}"
            positive_only[group][metric]={
                control:partial_rank_cells(nonzero_rows,response,fields)
                for control,fields in COVARIATES.items()}
    jackknife={}
    for group in GROUPS:
        jackknife[group]={}
        for metric in MEASURES:
            response=f"{group}__{metric}"
            jackknife[group][metric]=leave_one_cell_out_diagnostic(
                rows,response,COVARIATES["gbif_and_crossing_opportunity"])
    return {
      "leave_one_cell_out_sensitivity_gbif_and_crossing_adjusted":jackknife,
      "positive_reported_train_only_sensitivity":positive_only,
      "n_nonzero_reported_train_50km_cells":len(nonzero_rows),
      "nonzero_trains_cannot_be_assumed_operated_freight_only":True,
      "n_original_cells":len(cells),
      "n_eligible_cells_for_source_observations":len(rows),
      "n_GIFT_source_eligible_cells":sum(r["source_growth_form_taxon_match_fraction"]>=MIN_GIFT_TRAIT_COVERAGE for r in rows),
      "n_source_cells_zero_or_positive_reported_train_median":{
          "zero":sum(r["raw_reported_crossing_median_through_train_count"]==0 for r in rows),
          "positive":sum(r["raw_reported_crossing_median_through_train_count"]>0 for r in rows)},
      "n_source_cells_excluded_and_reasons":excluded,
      "descriptive_partial_rank_results":results,
      "within_region_cells_NOT_independent_country_replicates":True,
    }

def run(flora,fra,gift):
    if flora.get("complete") is not True or flora.get("source_errors"):
        raise ValueError("GBIF GRIIS original flora not successful")
    if fra.get("complete") is not True or set(BLOCKS)-set(fra.get("regions",{})):
        raise ValueError("FRA 2026 source blocks incomplete")
    traits=validate_traits(gift)
    out={"complete":True,"utc":datetime.now(timezone.utc).isoformat(),
         "study":"Exploratory two-US-macro-block observed crossing train counts versus GRIIS introduced plant recorded richness by GIFT growth form",
         "frozen_regions":list(BLOCKS),"recent_GBIF_positive_source_period":PERIOD,
         "rail_traffic_data_source_vintage":"FRA 2026 crossing inventory; reported train-count year 2018-2026, NOT a synchronized year panel",
         "n_independent_geographic_macro_blocks_at_most":2,
         "source_plant_records_NOT_true_railway_verge_establishments":True,
         "NOT_train_freight_seed_transport_or_mowing_wind_evolved_adaptation":True,
         "regions":{}}
    for block in BLOCKS:
        out["regions"][block]=build_block(block,flora["regions"][block],fra["regions"][block],traits)
        x=out["regions"][block]
        print(json.dumps({"block":block,"n_eligible_source_cells":x["n_eligible_cells_for_source_observations"],
          "main_train_vs_woody_rarefied_rho":x["descriptive_partial_rank_results"]["woody"]["rarefied20_expected_species_richness"]["gbif_and_crossing_opportunity"],
          "source_excluded_cells":x["n_source_cells_excluded_and_reasons"]}),flush=True)
    out["limitations"]=[
       "Only 2 prespecified US geographic blocks, 27 or fewer cells; insufficient spatially independent replication for a universal macro ecological claim.",
       "Reported train movements at road-rail crossing point do NOT measure freight-specific train volume, seed flow or slipstream gusts.",
       "Crossing sampling is biased by roads, towns and active crossing locations; log crossing count is an incomplete road/land-use control.",
       "The 2026 snapshot includes train-count source years 2018–2026 and cannot establish train exposure at 2018–2025 plant observation times.",
       "GBIF species richness is recorded species, not true biological species richness or independent occupancy, and record-level source provider exclusion was not done for US Southeast.",
       "GIFT species growth form is neither seed dispersal mode nor grass-cutting or wind resistance; unclassified taxa are kept unknown.",
       "Spatial gradients and phylogenetic/climatic/urban effects remain incompletely controlled; descriptive rank coefficients have no p-values or causal interpretation.",
       "Do not select a favorable trait, region, count rarefaction or covariate sensitivity for a causal adaptation claim."
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for name in ("flora","fra","gift","output"):p.add_argument("--"+name,required=True)
    a=p.parse_args();read=lambda filename:json.loads(Path(filename).read_text(encoding="utf-8"))
    out=run(read(a.flora),read(a.fra),read(a.gift))
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"complete":True,"regions":list(out["regions"]),"no_causal_railway_claim":True}),flush=True)

if __name__=="__main__":main()

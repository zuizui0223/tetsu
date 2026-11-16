#!/usr/bin/env python3
"""v9.7 source-only exploratory correlation matrix: do states where the
original NONRANDOM iNaturalist observer/site GBIF source sample has
more observer and coarse-location continuity show fewer apparently
"late-only" country-GRIIS speciesKey recorded categories?

Independent inputs: v9.5 non-identifying aggregate of exactly
9,600 original records before and after (200 per 48 Census states,
two deterministic occurrence/search source offsets), and v9.1
48-state GBIF original indexed source-record reporting-effort audit.
SAME 48 geographical state sources; no individual observer names,
record IDs, original occurrence coordinates or repeated surveys.

Report ALL 7 precommitted observer/site metrics × TWO responses,
raw Spearman and partial-rank conditional on late total ALL-Plantae
indexed record count and annualized post/early recording ratio,
plus leave-ONE-state-out signed coefficient range. This set is
POSTHOC EXPLORATORY, NOT pre-registered statistical validation,
not a causal measurement of detection bias or plant immigration.

No freight train seed carriage, newly established state-level alien
plant, railway mowing/gust/plasticity/evolution inference.
"""
from __future__ import annotations
import argparse,json,math
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import numpy as np
from analyze_macro_50km_rail_road_aliens_v34 import _rank

PREDICTORS=(
    "sampled_observer_label_recurrence",
    "sampled_site_recurrence_10km",
    "sampled_site_recurrence_50km",
    "sampled_observer_and_site_recurrence_10km",
    "sampled_observer_and_site_recurrence_50km",
    "sampled_GRIIS_record_at_previously_sampled_site_10km",
    "sampled_GRIIS_record_at_previously_sampled_site_50km",
)
RESPONSES=(
    "late_only_country_GRIIS_taxa_fraction_of_all_late_recorded_GRIIS_taxa",
    "late_only_country_GRIIS_taxa_count_NOT_immigrations",
)
CONTROLS=(
    "log_2023_25_all_Plantae_indexed_GBIF_record_count",
    "post_pre_annualized_all_Plantae_recording_ratio",
)
STATES=48
MIN_STATES=40


def source_state_frame(original95,original91):
    if (original95.get("complete") is not True or
        original91.get("complete") is not True or
        original95.get("nonrandom_source_probe_DO_NOT_EXTRAPOLATE_TO_ALL_INAT_USERS") is not True or
        original95.get("no_raw_observer_names_ids_record_ids_or_latlon_saved") is not True or
        original95.get("summary",{}).get("n_actual_GBIF_iNaturalist_original_sampled_rows_before")!=9600 or
        original95.get("summary",{}).get("n_actual_GBIF_iNaturalist_original_sampled_rows_after")!=9600):
        raise ValueError("Independently source-verified nonrandom original GBIF iNaturalist recording probe missing")
    if (original91.get("n_qualified_Census_US_states_for_both_year_windows")!=STATES or
        len(original95.get("states",{}))!=STATES or
        len(original91.get("state_FIPS_source_audit_rows",[]))!=STATES):
        raise ValueError("V9.5 and V9.1 original 48 Census-state source units disagree")
    original={r["Census_FIPS_state_code"]:r for r in original91[
        "state_FIPS_source_audit_rows"]}
    if set(original)!=set(original95["states"]):
        raise ValueError("Source state FIPS labels not identical")
    rows={}
    for code in sorted(original):
        s=original95["states"][code];r=original[code]
        if (s["source_pre_2018_21"]["n_original_rows"]!=200 or
            s["source_post_2023_25"]["n_original_rows"]!=200):
            raise ValueError("Fixed 2-by-100 original iNaturalist source sample missing")
        overlap=s["overlap"]; tiles=overlap["coarse_site_overlap"]
        known=overlap["n_late_source_rows_with_observer_identity"]
        if not 0<=known<=200:
            raise ValueError("Invalid deidentified source observer completeness")
        controls={
            CONTROLS[0]:math.log(r["n_post_all_Plantae_indexed_records"]),
            CONTROLS[1]:r["post_before_all_Plantae_annual_record_ratio"],
        }
        y={
            RESPONSES[0]:r["post_only_fraction_of_post_recorded_taxon_categories"],
            RESPONSES[1]:r["post_only_first_recorded_taxon_categories_NOT_immigrations"],
        }
        x={
            PREDICTORS[0]:s["fraction_late_original_rows_with_repeated_prior_observer_label"],
        }
        for km in ("10km","50km"):
            tile=tiles[km]
            denom=tile["n_late_records_with_original_valid_coordinates"]
            obsden=tile["n_late_records_with_both_observer_and_coarse_location"]
            grden=s["source_post_2023_25"]["n_source_country_GRIIS_listed_plant_records_in_sample"]
            x[f"sampled_site_recurrence_{km}"]=(
                tile["n_late_records_in_early_sampled_geographical_cells"]/denom
                if denom else None)
            x[f"sampled_observer_and_site_recurrence_{km}"]=(
                tile["n_late_records_with_same_observer_same_coarse_site_as_early"]/obsden
                if obsden else None)
            x[f"sampled_GRIIS_record_at_previously_sampled_site_{km}"]=(
                tile["n_late_GRIIS_records_in_repeated_geographical_cells"]/grden
                if grden else None)
        nums=[*controls.values(),*y.values(),*[a for a in x.values() if a is not None]]
        if any(type(n) not in (int,float) or not math.isfinite(n) for n in nums):
            raise ValueError("Source original GBIF state recording variable nonfinite")
        if any(z is not None and not 0<=z<=1 for z in x.values()):
            raise ValueError("Out-of-range recurrent source observer/location fraction")
        rows[code]={"predictors":x,"responses":y,"controls":controls,
                    "n_late_original_inat_sample_records":200}
    return rows


def rank_correlation(rows,xkey,ykey,partial=False):
    available=[r for r in rows if r["predictors"].get(xkey) is not None and
               r["responses"].get(ykey) is not None]
    if len(available)<MIN_STATES:return None
    x=np.asarray([r["predictors"][xkey] for r in available],float)
    y=np.asarray([r["responses"][ykey] for r in available],float)
    if len(np.unique(x))<5 or len(np.unique(y))<5:return None
    rx=_rank(x)
    ry=_rank(y)
    if partial:
        z=np.column_stack([np.asarray([r["controls"][k] for r in available],float)
                           for k in CONTROLS])
        A=np.column_stack([np.ones(len(x)),*(_rank(z[:,k]) for k in range(z.shape[1]))])
        if np.linalg.matrix_rank(A)<A.shape[1]:return None
        rx=rx-A@np.linalg.lstsq(A,rx,rcond=None)[0]
        ry=ry-A@np.linalg.lstsq(A,ry,rcond=None)[0]
    else:
        rx=rx-rx.mean();ry=ry-ry.mean()
    den=float(np.linalg.norm(rx)*np.linalg.norm(ry))
    return round(float(rx@ry/den),6) if den>1e-10 else None


def compute(source):
    v=list(source.values())
    out={}
    for predictor in PREDICTORS:
        out[predictor]={}
        for response in RESPONSES:
            raw=rank_correlation(v,predictor,response,False)
            adjusted=rank_correlation(v,predictor,response,True)
            leave={}
            for state in sorted(source):
                leave[state]=rank_correlation([
                    r for code,r in source.items() if code!=state],
                    predictor,response,True)
            valid=[x for x in leave.values() if x is not None]
            out[predictor][response]={
                "n_state_source_units":len(v),
                "raw_source_spearman_rho":raw,
                "record_count_and_annual_growth_partial_rank_rho":adjusted,
                "leave_one_Census_state_out_refits_not_independent_replicates":{
                    "n_deleted_states":len(leave),"n_estimable":len(valid),
                    "min_rho":min(valid) if valid else None,
                    "max_rho":max(valid) if valid else None,
                    "n_negative":sum(x<0 for x in valid),
                    "n_positive":sum(x>0 for x in valid)},
                "POSTHOC_hypothesis_no_p_value_multiple_predictors":True,
            }
    return out


def run(data95,data91):
    original=source_state_frame(data95,data91)
    out={
        "complete":True,"utc":datetime.now(timezone.utc).isoformat(),
        "study":"Posthoc 48-state GBIF recorded-source observation turnover vs later-only GRIIS source-taxonomy share",
        "n_same_Census_state_source_units":len(original),
        "n_original_NONRANDOM_iNat_GBIF_before_samples":9600,
        "n_original_NONRANDOM_iNat_GBIF_after_samples":9600,
        "posthoc_all_exploratory_source_predictors":list(PREDICTORS),
        "source_recording_control_fields":list(CONTROLS),
        "source_outcomes_not_colonization":list(RESPONSES),
        "full_exploratory_correlation_matrix":compute(original),
        "no_original_observer_labels_record_ids_or_geolocation_saved":True,
        "NO_rail_cargo_seed_causation_new_introduced_plants_mowing_wind_adaptation":True,
        "limitations":[
            "All 48 state observer overlap metrics originate from TWO nonrandom ordered GBIF occurrence/search source offset pages (0,1000) per state×period, not probability samples or all-user data.",
            "The study is POSTHOC exploratory (7 correlated observability metrics × 2 source first-record responses); no multiple-comparison or randomized design; coefficients are only descriptive.",
            "State FIPS source units vary in geography, climate, existing flora, agriculture, data platform sampling and number of uploads; residualizing total indexed records and annual growth does NOT control actual observer/site detection.",
            "Repeat source observer label may not identify same human; same 10/50km source grid does not mean an actual botanical plot or within-site revisits.",
            "The two indexed GBIF event bins are four and three years with retrospective submissions; later-only speciesKey category is NEVER proof of biological immigration, true absence or state-level nonnative status.",
            "Source GRIIS U.S.-country listing is not a state-specific alien checklist; there are no standardized flora surveys or independent country-level ecological experiments.",
            "No verified 2022 railway freight line operating history, viable seeds carried or shed, direct rail/road vegetation mowing/airflow exposure, plasticity or heritable evolution.",
        ]
    }
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("original95","source91","output"):p.add_argument("--"+name,required=True)
    args=p.parse_args()
    load=lambda k:json.loads(Path(getattr(args,k)).read_text(encoding="utf-8"))
    result=run(load("original95"),load("source91"))
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "number_state_source_units":result["n_same_Census_state_source_units"],
        "all_exploratory_observability_and_late_only_source_metrics":
            result["full_exploratory_correlation_matrix"],
        "no_obs_site_field_surveys_or_freight_seed_effect_claim":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

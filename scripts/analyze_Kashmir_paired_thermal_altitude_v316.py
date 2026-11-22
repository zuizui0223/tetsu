#!/usr/bin/env python3
"""v316 source-authentic Kashmir paired thermal-affinity/elevation interaction.

Main ecological independent unit is one original transect (2014, 31 sites);
2017 recovers the same sites and is NOT 31 independent new climatic regions.
Published Dryad original workbook, two publisher trait workbooks and v315
baseline are SHA verified. Results contain aggregates only (no site IDs, GPS,
taxon×site matrix, or raw elevations by site).

The study concerns thermal affinity indicators, NOT in situ temperatures.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,statistics
from collections import defaultdict
from pathlib import Path
import numpy as np

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix,DRYAD_SHA256
from analyze_thermal_affinity_reference_v315 import contrast,paired_summary
from analyze_matched_foreign_rail_road_v314 import normalize_name
from source_thermal_affinity_v315 import load_traits,kashmir_pairs,KASHMIR_ZONES
from source_matched_rail_road_v314 import audit_names,unresolved_keys

PRIMARY_YEAR=2014
YEARS=(2014,2017)
N_SITE_TOTAL=31
N_STATIONS=16
N_BETWEEN=15
TRAIT_SCENARIOS=("tichy_numeric","eive_all")
PERMUTATIONS=4999
SEED=316042
EXPECTED_LEGACY={
 "tichy_numeric":(16,0.002653787399820917),
 "eive_all":(16,0.06004162707505073),
}
# Original species × habitat matrix has 4 missing source cells across 2014/17.
# The XML/parser excludes any incomplete plot, never silently fills absence.


def min_support(rows,min_numeric=3):
    """Retain entire physical transects with adequately scored plant taxa."""
    return [r for r in rows
            if (r.get("difference") is not None and math.isfinite(r["difference"])
                and r["n_rail_numeric"]>=min_numeric
                and r["n_interior_numeric"]>=min_numeric
                and r["stratum"] in ("station","between")
                and math.isfinite(r["altitude"]))]


def _ols(rows):
    if not rows:raise ValueError("Empty ecological transects")
    elev=np.asarray([r["altitude"] for r in rows],float)
    strata=np.asarray([r["stratum"]=="station" for r in rows],float)
    yy=np.asarray([r["difference"] for r in rows],float)
    design=np.column_stack((np.ones(len(rows)),(elev-np.mean(elev))/100.,strata))
    if not np.isfinite(design).all() or not np.isfinite(yy).all():
        raise ValueError("Nonfinite original field altitude or thermal difference")
    if np.linalg.matrix_rank(design)<3:
        raise ValueError("Original 2014 transect altitude+station design lacks rank")
    coef=np.linalg.lstsq(design,yy,rcond=None)[0]
    return (float(coef[1]),float(coef[2]),float(np.mean(yy)),
            float(np.corrcoef(elev,yy)[0,1]) if np.std(yy)>0 else None)


def slope_model(rows,minimum_sites=20,minimum_per_stratum=6,minimum_altitude_span_m=40):
    """Original field-site-level model; pseudo-replicated quadrats forbidden."""
    all_rows=min_support(rows)
    alt=[r["altitude"] for r in all_rows]
    station=sum(r["stratum"]=="station" for r in all_rows)
    between=len(all_rows)-station
    span=max(alt)-min(alt) if alt else 0.
    reasons=[]
    if len(all_rows)<minimum_sites:reasons.append("insufficient_original_transects")
    if station<minimum_per_stratum or between<minimum_per_stratum:
        reasons.append("insufficient_station_or_between_source_stratum")
    if span<minimum_altitude_span_m:reasons.append("altitude_source_range_lt_40m")
    if len(all_rows)>=3 and len({round(v,8) for v in alt})>=2 and station and between:
        X=np.column_stack((np.ones(len(all_rows)),
                           (np.array(alt)-np.mean(alt))/100.,
                           np.array([r["stratum"]=="station" for r in all_rows],float)))
        if np.linalg.matrix_rank(X)<3:reasons.append("source_design_rank_deficient")
    else:reasons.append("source_design_rank_deficient")
    result={
        "estimable":not reasons,
        "reasons":reasons,
        "n_candidate_original_transects":len(rows),
        "n_transects":len(all_rows),
        "n_station":station,"n_between":between,
        "source_altitude_m_minmax":[float(min(alt)),float(max(alt))] if alt else None,
        "source_altitude_span_m":float(span),
        "altitude_slope_per_100m":None,
        "station_effect_at_mean_altitude":None,
        "mean_paired_thermal_difference":statistics.mean(r["difference"] for r in all_rows) if all_rows else None,
        "no_pseudoreplicated_quadrats":True,
    }
    if not reasons:
        beta,station_beta,mean,corr=_ols(all_rows)
        result.update(altitude_slope_per_100m=beta,
                      station_effect_at_mean_altitude=station_beta,
                      mean_paired_thermal_difference=mean,
                      source_sample_corr_altitude_vs_thermal_delta=corr)
    return result


def source_site_permutation(rows,seed=SEED,draws=PERMUTATIONS):
    """Permute source altitudes only within real station/between-site strata."""
    eligible=min_support(rows)
    fit=slope_model(eligible)
    if not fit["estimable"]:
        return {"estimable":False,"reasons":fit["reasons"],
                "n_original_transects":len(eligible)}
    observed=fit["altitude_slope_per_100m"]
    rng=random.Random(seed)
    strata={g:[i for i,r in enumerate(eligible) if r["stratum"]==g]
            for g in ("station","between")}
    values=[]
    for i in range(draws):
        changed=[dict(r) for r in eligible]
        for g,indices in strata.items():
            altitudes=[eligible[k]["altitude"] for k in indices]
            rng.shuffle(altitudes)
            for j,a in zip(indices,altitudes):
                changed[j]["altitude"]=a
        values.append(_ols(changed)[0])
    tail_up=(1+sum(v>=observed for v in values))/(draws+1)
    tail_both=(1+sum(abs(v)>=abs(observed) for v in values))/(draws+1)
    loo=[]
    for i in range(len(eligible)):
        after=eligible[:i]+eligible[i+1:]
        fit_lo=slope_model(after)
        if not fit_lo["estimable"]:
            raise ValueError("Leave-whole-site-out design failed support while original passes")
        loo.append(fit_lo["altitude_slope_per_100m"])
    return {
        "estimable":True,
        "n_original_transects":len(eligible),
        "n_station_original_stratum":len(strata["station"]),
        "n_between_original_stratum":len(strata["between"]),
        "observed_per_100m":observed,
        "null_mean_per_100m":float(np.mean(values)),
        "null_2p5_97p5_technical_NOT_field_CI":list(map(float,np.quantile(values,[.025,.975]))),
        "n_permutations":draws,
        "permutation_within_source_station_between_only":True,
        "upper_rank":float(tail_up),
        "two_sided_rank":float(tail_both),
        "tail_ranks_are_exploratory_not_causal_field_p_values":True,
        "whole_site_LOO_n":len(loo),
        "whole_site_LOO_n_positive":sum(v>0 for v in loo),
        "whole_site_LOO_minmax":[float(min(loo)),float(max(loo))],
    }


def build_normalized_plot_pairs(sets,normalize):
    """From real in-memory raw labels to same-year, same-site 2×100m² unions.

    Apply the SAME ambiguous-taxon identity mask to every year and habitat.
    No missing plot can be interpreted as an empty real plant community.
    """
    labels=set().union(*(x for x in sets.values() if x is not None))
    audit=audit_names(labels,normalize)
    blocked=unresolved_keys(audit,labels,normalize)
    all_sites=sorted({site for year,site,zone in sets})
    res={year:{} for year in YEARS}
    for year in YEARS:
        for site in all_sites:
            source=[sets.get((year,site,zone)) for zone in KASHMIR_ZONES]
            if any(x is None for x in source):
                continue
            rail=source[0]|source[1]
            interior=source[2]|source[3]
            norm_rail={audit[name] for name in rail}-{None}-blocked
            norm_int={audit[name] for name in interior}-{None}-blocked
            if not norm_rail or not norm_int:
                continue
            res[year][site]=(norm_rail,norm_int)
    return res


def make_pair(rail,interior,trait):
    """Index-score means with 3-taxa support, not °C or microclimate."""
    out=contrast({t:1. for t in rail},{t:1. for t in interior},
                 trait,min_species=3,species_equal=True)
    return {
        "estimable":out["estimable"],
        "difference":out["difference"],
        "n_rail_numeric":out["rail"]["n_numeric"],
        "n_interior_numeric":out["reference"]["n_numeric"],
        "rail_numeric_taxon_coverage":out["rail"]["numeric_taxon_coverage"],
        "interior_numeric_taxon_coverage":out["reference"]["numeric_taxon_coverage"],
        "rail_indifferent_frequency_fraction":out["rail"]["status_frequency_fraction"].get("indifferent",0.),
        "interior_indifferent_frequency_fraction":out["reference"]["status_frequency_fraction"].get("indifferent",0.),
    }


def genuine_source_site_metadata(xlsx):
    """Read exact original site heights, assert stable across years/plots."""
    from openpyxl import load_workbook
    plots,source_qa=parse_xlsx_matrix(xlsx)
    book=load_workbook(xlsx,read_only=True,data_only=True)
    try:
        ws=book.worksheets[0]
        rows=ws.iter_rows(values_only=True)
        header=next(rows)
        if (ws.max_row,ws.max_column)!=(249,129):
            raise ValueError("Source original field workbook shape drift")
        expected=("Name of the Site","Latitude ","Longitude ","Altitude","Year of sampling","Location")
        if tuple(header[:6])!=expected:
            raise ValueError("Original publisher altitude metadata headers changed")
        elevations=defaultdict(set)
        raw={}
        species_labels=header[6:]
        for row in rows:
            site=" ".join(str(row[0]).split())
            if not isinstance(row[3],(int,float)) or not math.isfinite(row[3]):
                raise ValueError("Original plot altitude is not a finite height")
            elevations[site].add(float(row[3]))
            year=int(row[4]);zone=" ".join(str(row[5]).split())
            key=(site,year,zone)
            observed=plots[key]
            raw[key]=({species_labels[i] for i in observed} if observed is not None else None)
        if len(elevations)!=31 or any(len(v)!=1 for v in elevations.values()):
            raise ValueError("Original same-site altitude not stable in two sampling years")
        altitude={site:next(iter(v)) for site,v in elevations.items()}
        if not (min(altitude.values())>=1500 and max(altitude.values())<=1800):
            raise ValueError("Original field altitude outside plausible Kashmir original range")
        if len({(y,s,z) for y,s,z in raw})!=248:
            raise ValueError("Source year×site×microhabitat records incomplete")
        return altitude,raw,source_qa
    finally:
        book.close()


def source_site_rows(allpairs,altitude,trait,year):
    rows={}
    for site,(rail,interior) in allpairs[year].items():
        site_out=make_pair(rail,interior,trait)
        rows[site]={
            "altitude":altitude[site],
            "stratum":"between" if site.startswith("Between ") else "station",
            **site_out
        }
    return rows


def replay_legacy_v315(data,traits,baseline):
    pairs=kashmir_pairs(data,normalize_name)
    out={}
    for scenario,trait in (("tichy_numeric",traits["tichy"]),("eive_all",traits["eive"])):
        replay=paired_summary(pairs,trait,min_species=3,min_pairs=8)
        original=next(x for x in baseline["paired_india"] if x["scenario"]==scenario)
        n,mean=EXPECTED_LEGACY[scenario]
        if (replay["n_pairs"]!=n or original["n_pairs"]!=n
                or abs(replay["mean_difference"]-mean)>1e-10
                or abs(replay["mean_difference"]-original["mean_difference"])>1e-10):
            raise ValueError("v315 actual real 16 station year-2014 paired thermal baseline NOT reproduced")
        out[scenario]={"n_original_stations":n,"mean_same_station_diff_verified":float(mean),
                       "strict_legacy_replay_passed":True}
    return out


def _bounded_stats(rows):
    values=min_support(rows)
    if not values:
        return {"n":0}
    return {
        "n":len(values),
        "source_station_strata":{g:sum(r["stratum"]==g for r in values)
                                 for g in ("station","between")},
        "mean_rail_numeric_species":statistics.mean(r["n_rail_numeric"] for r in values),
        "mean_interior_numeric_species":statistics.mean(r["n_interior_numeric"] for r in values),
        "min_rail_numeric_species":min(r["n_rail_numeric"] for r in values),
        "min_interior_numeric_species":min(r["n_interior_numeric"] for r in values),
        "mean_rail_numeric_species_fraction":statistics.mean(r["rail_numeric_taxon_coverage"] for r in values),
        "mean_interior_numeric_species_fraction":statistics.mean(r["interior_numeric_taxon_coverage"] for r in values),
        "mean_rail_Tichy_temperature_indifferent_fraction":statistics.mean(r["rail_indifferent_frequency_fraction"] for r in values),
        "mean_interior_Tichy_temperature_indifferent_fraction":statistics.mean(r["interior_indifferent_frequency_fraction"] for r in values),
    }


def analyze(data,traitdir,baseline,draws=PERMUTATIONS):
    original=Path(data)/"Matrix_plot.xlsx"
    if hashlib.sha256(original.read_bytes()).hexdigest()!=DRYAD_SHA256:
        raise ValueError("Original original publisher source hash changed")
    altitude,raw,origqa=genuine_source_site_metadata(original)
    traits=load_traits(Path(traitdir),normalize_name)
    pairs=build_normalized_plot_pairs(raw,normalize_name)
    legacy=replay_legacy_v315(data,traits,baseline)
    ret={
      "version":"v316",
      "design":"Exploratory after v315, equal-physical-site temperature-index rail vs interior gradient",
      "original_Dryad_source_sha256":DRYAD_SHA256,
      "original_source_metadata_QA":origqa,
      "original_site_altitude_range_m":[min(altitude.values()),max(altitude.values())],
      "n_publisher_original_transects":31,
      "original_station_n":16,"original_between_n":15,
      "same_original_source_altitude_all_2014_2017_plots":True,
      "source_trait_sha256":{k:traits["metadata"]["sources"][k]["sha256"] for k in ("tichy","eive")},
      "v315_original_16_station_baseline_replayed":legacy,
      "data_reuse_note":"Original source per-site species, names, site altitudes, coordinates never exported",
      "local_source_climate_temperature_measured":False,
      "inference_scope":"one Kashmir corridor with limited elevation range, no independent multi-region climate effect",
      "main_2014":{},"year_2017_sensitivity":{},"temporal_change_sensitivity":{}
    }
    records_by_scenario={}
    for scenario,trait in (("tichy_numeric",traits["tichy"]),("eive_all",traits["eive"])):
        records_by_scenario[scenario]={}
        for year in YEARS:
            site_rows=source_site_rows(pairs,altitude,trait,year)
            records_by_scenario[scenario][year]=site_rows
            values=list(site_rows.values())
            model=slope_model(values)
            perm=source_site_permutation(values,seed=SEED+(year-2014)*11+(0 if scenario=="tichy_numeric" else 5),
                                         draws=draws) if model["estimable"] else {"estimable":False,"reason":"source model ineligible"}
            out={"model":model,"permutation":perm,"trait_support":_bounded_stats(values),
                 "original_rows_after_require_all_4_source_plots":len(site_rows),
                 "original_unpaired_missing_or_empty_source_sites":31-len(site_rows),
                 "no_station_ids_or_raw_temperature_values_exported":True}
            where="main_2014" if year==2014 else "year_2017_sensitivity"
            ret[where][scenario]=out
        first,second=records_by_scenario[scenario][2014],records_by_scenario[scenario][2017]
        change=[]
        for site in sorted(first.keys()&second.keys()):
            a,b=first[site],second[site]
            if not (a["estimable"] and b["estimable"]):
                continue
            change.append({"altitude":a["altitude"],"stratum":a["stratum"],
                           "difference":b["difference"]-a["difference"],
                           "n_rail_numeric":min(a["n_rail_numeric"],b["n_rail_numeric"]),
                           "n_interior_numeric":min(a["n_interior_numeric"],b["n_interior_numeric"])})
        model=slope_model(change)
        ret["temporal_change_sensitivity"][scenario]={
            "n_original_same_2014_2017_site_pairs_with_both_trait_scores":len(change),
            "model":model,
            "permutation":source_site_permutation(change,seed=SEED+35+(0 if scenario=="tichy_numeric" else 5),draws=draws)
                          if model["estimable"] else {"estimable":False,"reason":"source temporal common-cohort ineligible"},
            "not_independent_2017_replicate":True,
        }
    # Core result is the 2014 source blocked within-transect slope only.
    return ret


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--traits",type=Path,required=True)
    parser.add_argument("--baseline",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    baseline=json.loads(args.baseline.read_text())
    result=analyze(args.data,args.traits,baseline)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print(json.dumps({
      "status":"v316_COMPLETED_ORIGINAL_PUBLISHER_SOURCE_KASHMIR_PAIRED_THERMAL_ELEVATION",
      "original_2014":{k:{
          "n":v["model"]["n_transects"],
          "beta_T_per_100m":v["model"]["altitude_slope_per_100m"],
          "permutation_positive_upper_rank":v["permutation"].get("upper_rank"),
          "whole_site_loo_positive":v["permutation"].get("whole_site_LOO_n_positive"),
          "mean_pair_difference":v["model"]["mean_paired_thermal_difference"]}
           for k,v in result["main_2014"].items()},
      "year_2017_sensitivity":{k:v["model"]["altitude_slope_per_100m"] for k,v in result["year_2017_sensitivity"].items()},
      "source_climate_hypothesis_between_regions_identified":False,
    },ensure_ascii=False,indent=2))


if __name__=="__main__":main()

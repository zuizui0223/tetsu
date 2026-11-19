#!/usr/bin/env python3
"""v28.6 genuine Czech railway 14 station-wide floras against SAME-LOCATION
nearby independent CzechVeg vegetation: isolate HEIGHT vs SEED-mass spatial
differentiation under the original v24.4 reproducible strict 100m-cell
maximum-coverage matched reference cohort.

Unlike internal four-study v28.5 null, this compares RAIL and NONTRANSPORT
context real observed plant communities in the same broader Czech region.
Source/authenticity:
- Kutlvašr et al (2024) Preslia original 39 station floras, original 2020-22,
  Preslia S1+corrected S2 SHA verified in GitHub Actions.
- third-party 2026 SR70 station centroids exact 35/39 name match (not actual
  2020-22 botanical plot GPS);
- CzechVeg-Open 2026 117739 original relevé headers and 2510340 source
  species-row observations (MD5 verified, 2015-23), with positive original
  coordinate uncertainty <=250m, no rail/transport words in source locality
  and habitat fields, geographically within 10km of the railroad station.
  A "nontransport text" field does NOT prove natural or never road-disturbed.
- LEDA Oldenburg original canopy height (m) and seed mass (mg) exact source
  species-median values, uniform zlog scaling over verified station+reference
  matched taxa, NOT restandardized separately for railway vs local source.

Use original v24.4 true binary MILP selection: 16 eligible rail stations
after original source plant/traits QA, select maximal 14 actual stations
across 4 author railway stretches. Each station exactly five DIFFERENT
CzechVeg source plot IDs from >=2 independently sampled CzechVeg studies,
nonoverlapping ORIGINAL latitude/longitude rounded 3 decimals (~100m)
across all 70 true source reference plots. GIS uncertainty <=250m,
so 100m grid cells do NOT certify physical polygon non-overlap.
Matching uses only distance and original study metadata, not outcomes.
The original 14 station rail community samples all station habitat areas
(3 original station zones), the reference comprises five individual
nontransport-text relevés; no physical equal-AREA control is established.

For each of 100 original v24.4 pre-fixed deterministic technical
species-rarefaction draws:
  - each SAME original rail station and pooled 5-reference source plots
    are random-subset to exact SAME source number of 2D-known taxa
    min(original site rail known species, reference known species);
  - preserve the SAME selected field plot IDs and station cohort across all
    draws, and use identical site pairs to compute pairwise absolute CWM
    difference separately for zlog canopy height, zlog seed mass, and
    Euclidean 2D CWM distance;
  - report both ALL pairs of 14 source stations (91 dependent pairs),
    and within each of four original rail stretches (equal stretch means).
  - site 2D trait variance trace alpha is distinct from pairwise beta.
  - leave one actual original rail study stretch out, no 91 independent n.

Source v24.4 canonical 100-draw 3-decimal reference must be reproduced
within 0.001 for the prior 2D means (rail 0.270484 vs local ref 0.697218,
rail minus local=-0.426733). If not, FAIL CLOSED: not silently a new
different control cohort.

The equal-species rarefaction fixes count per paired site but cannot make
the STATION-WIDE vs 5 SMALL PLOTS source field-grain, age, soil, climate,
management, local habitat classes or source coordinate error equivalent.
This is observational matched-reference evidence, NOT railway-specific
selection, functional evolution, mower/wind adaptation or seed arrival.
Output only aggregated Czech publisher study/field species coverage and
macro estimates, NEVER raw original site coords, plant taxon sets, source
plot IDs, individual source study authors, or exact station allocation.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_Czech_39station_true_local_nontransport_plant_LEDA_beta_v240 import (
    LEDA_SHA,source_station_and_reference_coords,source_species_to_traits,
    source_valid_local_control_plots,source_species,quart
)
from validate_Czech_rail_local_reference_alternative_14station_assignments_v243 import (
    station_reference_profiles
)
from validate_Czech_rail_station_local_reference_geo_cells_unique_v244 import (
    exact_max_stations_with_cell_uniqueness
)

SOURCE_STATIONS=14
SOURCE_REFERENCE_PLOTS=70
SOURCE_STRETCHES=4
SOURCE_ALL_SITE_PAIRS=91
SOURCE_RAREFACTION_DRAWS=100
SOURCE_SEED=20261117+244
SOURCE_SEED_STEP=291
TOL=0.001
OLD_100M_SOURCE_2D={
    "absolute_rail_CWM_distance":0.270484,
    "absolute_local_control_CWM_distance":0.697218,
    "rail_minus_local_nontransport_CWM_beta":-0.426733
}
METRICS=("CWM_2D_Euclidean","CWM_height_1D_absolute","CWM_seed_mass_1D_absolute")


def group_metric(vecs:dict[str,tuple[float,float]], original_stretches:dict[str,str],
                 metric:str, *, within_only:bool=False, equal_stretches:bool=False)->float|None:
    if metric not in METRICS:raise ValueError("Unknown real source plant trait")
    if set(vecs)!=set(original_stretches):
        raise ValueError("Railway station and original independent reference site sets differ")
    by=defaultdict(list)
    ids=sorted(vecs)
    for i,j in combinations(ids,2):
        x=vecs[i]; y=vecs[j]
        if metric=="CWM_2D_Euclidean":value=math.dist(x,y)
        elif metric=="CWM_height_1D_absolute":value=abs(x[0]-y[0])
        else:value=abs(x[1]-y[1])
        if not within_only or original_stretches[i]==original_stretches[j]:
            by[original_stretches[i] if original_stretches[i]==original_stretches[j] else "between_source_stretches"].append(value)
    if not by:return None
    if within_only and equal_stretches:
        return mean(mean(vals) for vals in by.values())
    return mean(value for group in by.values() for value in group)


def comparative_stats(rail:dict,control:dict,source_stretches:dict)->dict:
    if set(rail)!=set(control) or set(rail)!=set(source_stretches):
        raise ValueError("Original matched rail/reference station identities changed")
    rv={k:rail[k]["CWM"] for k in rail}
    cv={k:control[k]["CWM"] for k in control}
    metrics={}
    for m in METRICS:
        metrics[m]={}
        for scope in ("all_14_source_station_pairs","within_original_stretches_stationpair_equal","within_stretches_equal_original_stretch"):
            w=scope!="all_14_source_station_pairs"
            equal=scope=="within_stretches_equal_original_stretch"
            ra=group_metric(rv,source_stretches,m,within_only=w,equal_stretches=equal)
            re=group_metric(cv,source_stretches,m,within_only=w,equal_stretches=equal)
            if ra is None or re is None:
                raise ValueError("No original same rail stretch plots with independently published vegetation")
            metrics[m][scope]={
                "rail":ra,"local_control":re,"rail_minus_local_control":ra-re
            }
    metrics["site_2D_trait_alpha_variance_trace"]={
        "rail":mean(rail[k]["trace"] for k in rail),
        "local_control":mean(control[k]["trace"] for k in control),
        "rail_minus_local_control":mean(rail[k]["trace"]-control[k]["trace"] for k in rail)
    }
    return metrics


def study_stretch_leave_one_out(rail:dict,reference:dict,stretches:dict)->dict:
    if len(set(stretches.values()))!=SOURCE_STRETCHES:
        raise ValueError("The original 4 source Preslia railway stretches must be kept separate")
    output={}
    for g in sorted(set(stretches.values())):
        keep={k for k in rail if stretches[k]!=g}
        if len(keep)<8 or len({stretches[k] for k in keep})<3:
            raise ValueError("Original source 4 railway stretch LOO support collapsed")
        r={k:rail[k] for k in keep}
        n={k:reference[k] for k in keep}
        meta={k:stretches[k] for k in keep}
        results=comparative_stats(r,n,meta)
        output[g]={m:results[m]["all_14_source_station_pairs"]["rail_minus_local_control"]
                   for m in METRICS}
    return output


def summarize_technical_draws(rows:list[dict],stretches_loo:list[dict])->dict:
    if len(rows)<20 or len(rows)!=len(stretches_loo):raise ValueError("Source technical independent subset draws missing")
    metrics={}
    for trait in METRICS+("site_2D_trait_alpha_variance_trace",):
        metrics[trait]={}
        scopes=(("station_equal_2D_alpha" if trait=="site_2D_trait_alpha_variance_trace"
                 else "all_14_source_station_pairs"),)
        if trait in METRICS:
            scopes+=("within_original_stretches_stationpair_equal","within_stretches_equal_original_stretch")
        for scope in scopes:
            values=[x[trait] if trait=="site_2D_trait_alpha_variance_trace" else x[trait][scope] for x in rows]
            def fld(key):
                return [x[key] for x in values]
            diffs=fld("rail_minus_local_control")
            metrics[trait][scope]={
                "source_technical_draws":len(rows),
                "rail_original_study_mean":mean(fld("rail")),
                "original_geographic_local_reference_mean":mean(fld("local_control")),
                "rail_minus_local_reference_mean":mean(diffs),
                "original_technical_rarefaction_q025_NOT_field_CI":quart(diffs,.025),
                "original_technical_rarefaction_q975_NOT_field_CI":quart(diffs,.975),
                "technical_draws_rail_minus_ref_negative":sum(x<0 for x in diffs),
                "technical_draws_rail_minus_ref_positive":sum(x>0 for x in diffs),
                "source_station_pairs_depend_on_shared_real_stations":True
            }
    loo={}
    for trait in METRICS:
        current={s:[] for s in stretches_loo[0]}
        for draw in stretches_loo:
            for s,p in draw.items():current[s].append(p[trait])
        loo[trait]={
            "original_field_stretches_heldout":len(current),
            "mean_effects_by_heldout_stretch_no_source_names":sorted(mean(a) for a in current.values()),
            "each_source_stretch_holdout_rail_less_local_reference":
                all(mean(a)<0 for a in current.values())
        }
    return {"source_rarefied_original_rail_vs_local_reference_trait_components":metrics,
            "source_stretch_leave_one_out":loo}


def run(s1:Path,s2:Path,sr70:Path,headers:Path,species:Path,
        height:Path,seed:Path,draws:int=SOURCE_RAREFACTION_DRAWS)->dict:
    if draws!=SOURCE_RAREFACTION_DRAWS:
        raise ValueError("Frozen 100 original source v24.4 rarefaction draws required to verify original beta")
    h,hq=read_leda(height,"height_m")
    sm,sq=read_leda(seed,"seed_mass_mg")
    stations_original,geo,read,original_qa=source_station_and_reference_coords(s1,s2,sr70,headers)
    plants,plant_qa=source_species(species,read)
    traits,traitqa=source_species_to_traits(stations_original,plants,h,sm)
    stations,stationqa=source_valid_local_control_plots(stations_original,plants,geo,traits)
    if len(stations)!=16:
        raise ValueError("Authentic original Czech 16 true sites for exact reference MILP missing")
    choice,matching=exact_max_stations_with_cell_uniqueness(stations,geo,decimals=3)
    if choice is None or len(choice)!=SOURCE_STATIONS or matching["n_source_independent_releve_plot_IDs"]!=SOURCE_REFERENCE_PLOTS:
        raise ValueError("Genuine v24.4 source 14 railway + 70 nontransport unique 100m WGS84 cells cohort drift")
    if len({stations[i]["stretch"] for i in choice})!=SOURCE_STRETCHES:
        raise ValueError("Genuine source four independently sampled railway stretches not present")
    stretched={stations[i]["id"]:stations[i]["stretch"] for i in choice}
    all_draws=[];leave=[]
    original_2d_check=[]
    for d in range(draws):
        source_rail,source_ctrl,original_2d,qa=station_reference_profiles(
            stations,choice,geo,traits,SOURCE_SEED+SOURCE_SEED_STEP*d)
        if len(source_rail)!=SOURCE_STATIONS or set(source_rail)!=set(stretched):
            raise ValueError("Original paired source plants beyond selected field site cohort")
        measures=comparative_stats(source_rail,source_ctrl,stretched)
        if abs(measures["CWM_2D_Euclidean"]["all_14_source_station_pairs"]["rail_minus_local_control"]-
               original_2d["rail_minus_local_nontransport_CWM_beta"])>1e-10:
            raise ValueError("Our new height/seed original trait beta decomposition fails to reproduce genuine source v24.4 prior 2D metric on SAME draw")
        all_draws.append(measures)
        leave.append(study_stretch_leave_one_out(source_rail,source_ctrl,stretched))
        original_2d_check.append(original_2d)
    for old,value in OLD_100M_SOURCE_2D.items():
        actual=mean(q[old] for q in original_2d_check)
        if abs(actual-value)>TOL:
            raise ValueError(f"Original v24.4 source 100m centroid-unique 14site 70plot same rarefaction benchmark mismatch {old}: got {actual:.9f} expected {value:.6f}")
    full=summarize_technical_draws(all_draws,leave)
    return {
        "version":"v28.6",
        "original_railway_study":"CZ Kutlvašr et al 2024 Preslia original 2020–22 39 station flora, only 14 accepted",
        "original_local_vegetation_study":"CzechVeg-Open 2026 actual 2015–23 other-study original nontransport-text site floras",
        "source_real_CZ_preslia_rail_and_CzechVeg_original_QA":original_qa,
        "source_original_independently_published_CzechVeg_species_QA":plant_qa,
        "source_original_LEDA_union_axis_QA":traitqa,
        "source_16_stations_five_reference_plot_candidates_QA":stationqa,
        "source_strict_milp_100m_rounding_reference_assignment_QA":matching,
        "source_original_14_railway_stations_70_geo_unique_reference_plots":True,
        "original_4_preslia_stretches_rail_stations":SOURCE_STRETCHES,
        "true_source_matched_reference_species_rarefaction_technical_draws":draws,
        "verified_prior_original_v24_4_twoD_CWM_source_results":{
            k:mean(x[k] for x in original_2d_check) for k in OLD_100M_SOURCE_2D},
        "original_station_and_stretch_equal_trait_beta_1D_vs_local_controls":full,
        "source_reference_not_equal_surveyed_field_plot_area_or_habitat_classes":True,
        "source_nontransport_habitat_word_filter_not_independently_onsite_nonrail_verified":True,
        "source_reference_geocoordinate_uncertainty_max_250m_cannot_certify_true_field_nonoverlap":True,
        "original_railway_station_centroid_not_exact_original_field_plot_GPS":True,
        "no_evolution_mowing_air_gust_or_train_seed_transport_identified":True,
        "no_original_station_plot_ids_GPS_species_by_site_or_trait_table_exported":True
    }


def main()->None:
    p=argparse.ArgumentParser()
    for name in ("s1","s2","sr70","headers","species","height","seed","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--draws",type=int,default=SOURCE_RAREFACTION_DRAWS)
    a=p.parse_args()
    report=run(a.s1,a.s2,a.sr70,a.headers,a.species,a.height,a.seed,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
        "version":"v28.6",
        "source_verified_original_publisher_14_station_70_reference_100m_cell":True,
        "prior_v244":report["verified_prior_original_v24_4_twoD_CWM_source_results"],
        "original_new_height_and_seed_beta":report["original_station_and_stretch_equal_trait_beta_1D_vs_local_controls"],
    },ensure_ascii=False))


if __name__=="__main__":
    main()

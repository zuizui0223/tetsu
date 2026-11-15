#!/usr/bin/env python3
"""tetsu v16.0 POSTHOC group-specific plant identity test on the SAME original
Czech railway vs ROAD source network.

Actual botanical observations: Kutlvašr et al (2024), Preslia 96, 247-266,
DOI 10.23855/preslia.2024.247 (2020-22, 39 full-station source samples).
Original publishers' S1 exact current-station flora + corrected S2, Pladias
traits and 35/39 original full names joined to publicly sourced 2026 SR70
railway STATION CENTROIDS (not individual botanical plot GPS).
Real Geofabrik 2023-01-01 Czech PBF SHA
26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747.
Previously original-source verified v15.5 and v15.6.

Hypothesis-generating question AFTER seeing v15.7-v15.9:
Does the log(railway real shortest physical route / drivable road path)
predict original explicitly tagged ARCHAEOPHYTE / NEOPHYTE taxon identity
beyond basic geography, whole-assemblage functional profile, source
station area, study route-stretch membership and real ROAD distance?

Do not merge new/old alien classes with "native" (source untagged unknown),
and do not equate archaeophyte/neophyte with exact arrival year, train
seed-carriage rate, succession, mowing/airflow exposure or adaptation.
The 594 matched original station pairs are dependent; ONLY five original
study line stretches from ONE country. All 3 class effects are from the
SAME one original field study, NOT 3 meta study effects.

Base main v15.6 544 trait-known model is recomputed exactly as an
independent live-source regression checksum before the group analyses.
Compare source subgroups on the SAME 594 matched, geodesic>=2km,
predeclared max detour≤8*geodesic+4km original observation universe.
Models are fixed v15.6 ridge alpha=10 and 5 original-stretch holdouts;
no outcome-dependent hyperparameter or station subset selection.
Two sensitivities: subgroup-specific source station log1p richness
versus main model's original whole known-flora richness.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from collections import Counter
from statistics import mean

from analyze_Czech_2023_true_rail_road_osm_routes_v155 import (
    run as network_run,EXPECTED_SOURCE_PBF_BYTES
)
from analyze_Czech_2023_true_rail_road_heldout_taxon_identity_v156 import (
    pair_records,cv_section_heldout,rail_coefficient
)
from analyze_Czech_2024_real_archaeophyte_neophyte_geographic_identity_v157 import (
    read_source_status,group_sets
)
from analyze_Czech_2024_39station_functional_identity_v143 import jaccard

DOI="10.23855/preslia.2024.247"
OSM_SHA="26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
PUBLISHER_S1_SHA="df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06"
SR70_SHA="fc43a1969af64f17285c7c6e3a30b00c3994d0ba1cfda6db09d8e4baaa166019"
N_STATIONS=35
N_COMPARABLE_SOURCE_PAIRS=594
ORIGINAL_5_STRETCHES=5
FROZEN_V156_KNOWN_M1=0.0027606838643546936
FROZEN_V156_KNOWN_M2=0.0028333736388869778
FROZEN_V156_ALL_M1=0.00227047338324627
FROZEN_V156_ALL_M2=0.0023482748839512187
SUBGROUP_EXPECTED={"archaeophyte":(151,1763),
                   "neophyte":(150,929),
                   "explicit_aliens_union":(301,2692)}
RICHNESS_MODES=("subgroup_log1p_richness","whole_original_trait_known_log_richness")

def assert_original_reference(cv,name,road,rail):
    if not cv["identified"] or cv["n_holdout_original_stretches"]!=ORIGINAL_5_STRETCHES:
        raise ValueError("Authentic original Czech v15.6 reference not estimable")
    result=cv["equal_original_stretch_weighted_mse"]
    if (not math.isclose(result["M1_plus_actual_2023_road_route"],road,rel_tol=1e-9,abs_tol=1e-12)
        or not math.isclose(result["M2_plus_actual_2023_road_and_rail_ratio"],rail,rel_tol=1e-9,abs_tol=1e-12)):
        raise ValueError(name+" original verified Czech road+rail v15.6 reference MSE drifted")


def reassign_original_source_class_outcome(original_pair_records,source_sites,original_class_sets,
                                           richness_mode="subgroup_log1p_richness"):
    """Preserve SAME original train-road route, geography, functional controls and
    station-pair support; only species class Jaccard and station richness vary.

    Fail on missing source subgroup occupancy/double-zero: no fabricated 0
    Jaccard for a group absent at both site endpoints.
    """
    if richness_mode not in RICHNESS_MODES:
        raise ValueError("Undeclared original richness sensitivity")
    if len(original_class_sets)!=len(source_sites)!=N_STATIONS:
        raise ValueError("Original SR70 exact station count differs from 35")
    if any(not original_class_sets[i] for i in range(len(source_sites))):
        raise ValueError("Explicit original subgroup absent in actual source station; shared source support would drift")
    out=[]
    for r in original_pair_records:
        i,j=r["i"],r["j"]
        if i==j or not(0<=i<len(source_sites) and 0<=j<len(source_sites)):
            raise ValueError("Source comparison station indices invalid")
        x,y=original_class_sets[i],original_class_sets[j]
        val=jaccard(x,y)
        if val is None:
            raise ValueError("Never convert double-zero original class to fake presence similarity")
        if not(0<=val<=1):
            raise ValueError("Source original subgroup Jaccard outside range")
        nr=dict(r)
        nr["features"]=dict(r["features"])
        nr["jaccard"]=float(val)
        if richness_mode=="subgroup_log1p_richness":
            la,lb=math.log1p(len(x)),math.log1p(len(y))
            nr["features"]["site_ln_species_count_mean"]=(la+lb)/2
            nr["features"]["site_ln_species_count_absdiff"]=abs(la-lb)
        nr["source_original_explicit_alien_group_richness_mode"]=richness_mode
        out.append(nr)
    if len(out)!=len(original_pair_records):
        raise ValueError("Changed original authentic rail-road-geo pair support under plant class subsetting")
    return out


def group_cv_paired(original,sites,class_groups):
    if len(original)!=N_COMPARABLE_SOURCE_PAIRS:
        raise ValueError("Czech author-original rail/road/geo common support changed from 594")
    if len({r["i"] for r in original}|{r["j"] for r in original})!=N_STATIONS:
        raise ValueError("One original station lost all comparable physical network routes")
    results={}
    for label,(taxa,incidence) in SUBGROUP_EXPECTED.items():
        d=class_groups[label]
        if (d["n_source_distinct_observed_class_taxa"]!=taxa or
            d["n_source_site_taxon_incidence"]!=incidence or
            d["n_actual_original_site_stations_with_class"]!=N_STATIONS):
            raise ValueError("Original Czech author classification source class counts drifted for "+label)
        sets=d["source_site_sets"]
        variants={}
        for mode in RICHNESS_MODES:
            rows=reassign_original_source_class_outcome(original,sites,sets,mode)
            cv=cv_section_heldout(rows)
            if not cv["identified"] or cv["n_holdout_original_stretches"]!=ORIGINAL_5_STRETCHES:
                raise ValueError("Class "+label+" failed five section source holdout: do not infer a group effect")
            variants[mode]={
                "original_stretch_heldout_CV":cv,
                "in_sample_rail_route_association_NOT_transport":rail_coefficient(rows),
                "original_railroad_geographic_site_pairs_IDENTICAL_to_v156":len(rows),
            }
        results[label]={
          "n_distinct_source_author_tagged_taxa":taxa,
          "n_original_source_35station_taxon_incidence":incidence,
          "mean_source_taxa_per_author_study_station":incidence/N_STATIONS,
          "n_actual_stations_with_source_class":N_STATIONS,
          "class_specific_and_full_flora_richness_sensitivity":variants,
        }
    return results


def validate_source_topology_and_reference(topology,pairs,sites,guild,s1,sr70):
    if topology["original_botanical_source_DOI"]!=DOI:
        raise ValueError("Original field study DOI mismatch")
    if (topology["actual_download_source_SHA256"]!=OSM_SHA
        or topology["original_OSM_archive_bytes"]!=EXPECTED_SOURCE_PBF_BYTES):
        raise ValueError("2023 Czech historical OSM archive SHA/byte provenance drift")
    src=topology["original_author_plant_source_and_station_geographic_provenance"]
    if src["2026_SR70_source_provenance"]["downloaded_csv_sha256"]!=SR70_SHA:
        raise ValueError("Czech SR70 station point original derivative source sha mismatch")
    if len(sites)!=N_STATIONS or len({x["section"] for x in sites})!=ORIGINAL_5_STRETCHES:
        raise ValueError("Original five Czech botanical stretches or exact name matched stations missing")
    if len(pairs)!=595:
        raise ValueError("Czech true graph must yield 35 choose 2 dependent station pairs")
    reference={}
    for kind,road,rail in (
        ("both_Pladias_known",FROZEN_V156_KNOWN_M1,FROZEN_V156_KNOWN_M2),
        ("all_original_S1_current",FROZEN_V156_ALL_M1,FROZEN_V156_ALL_M2)):
        r,omitted=pair_records(pairs,sites,guild,kind)
        if len(r)!=N_COMPARABLE_SOURCE_PAIRS or omitted!={"no_geography_and_both_true_rail_road_route_support":1}:
            raise ValueError("Actual original 594 pair comparative cohort not recovered")
        c=cv_section_heldout(r)
        assert_original_reference(c,kind,road,rail)
        reference[kind]={"original_5_stretch_CV_verified_exactly":True,
                         "original_source_road_mse":road,
                         "original_source_road_plus_rail_mse":rail,
                         "n_original_station_pairs":len(r)}
        if kind=="both_Pladias_known":
            base_records=r
    return base_records,reference


def run(s1,s2,disp,life,sr70,pbf):
    topo,pairs,sites,guild=network_run(s1,s2,disp,life,sr70,pbf,True)
    base,reference=validate_source_topology_and_reference(topo,pairs,sites,guild,s1,sr70)
    labels,provenance=read_source_status(s1)
    if provenance["original_S1_publisher_SHA256"]!=PUBLISHER_S1_SHA:
        raise ValueError("Original publisher source explicit tags SHA drift")
    groups=group_sets(sites,labels)
    result=group_cv_paired(base,sites,groups)
    # Contrasts are predictive diagnostics on the same original pair/fold
    # opportunities, not species genetic or arrival-stage effect sizes.
    primary={}
    for label,v in result.items():
        cv=v["class_specific_and_full_flora_richness_sensitivity"]["subgroup_log1p_richness"]["original_stretch_heldout_CV"]
        primary[label]={
          "source_stretch_equal_rail_increment_MSE_gain_pct":
              cv["rail_added_MSE_improvement_percentage_after_road"],
          "n_original_stretches_positive_rail_increment":
              cv["rail_improves_how_many_stretches"],
          "n_original_stretches":cv["n_holdout_original_stretches"],
        }
    return {
      "complete":True,
      "version":"v16.0_POSTHOC_real_Czech_original_aliens_vs_physical_rail_road_path_holdout",
      "original_author_field_study_DOI":DOI,
      "historical_Geofabrik_Czech_2023_01_01_full_source_SHA256":OSM_SHA,
      "original_author_Czech_explicit_alien_status_source_SHA256":PUBLISHER_S1_SHA,
      "source_original_2026_SR70_station_centroid_SHA256":SR70_SHA,
      "source_2020_22_author_field_and_2026_station_centroid_provenance":topo,
      "independently_recomputed_original_v156_whole_flora_road_rail_reference":reference,
      "posthoc_explicit_alien_class_source_road_adjusted_predictive_results":result,
      "posthoc_source_class_rail_increment_comparison":primary,
      "same_594_source_pairs_35_original_stations_5_original_stretches":True,
      "biological_limits":[
        "POSTHOC after v15.7-v15.9: not preregistered, no corrected inference/independent confirmation for comparing group effects.",
        "Original S1 arch/neo tags are time-of-introduction categories, NOT year of entry into each railway station, seed passage rate, plant establishment stage or author-confirmed native status for untagged species.",
        "Original 2020-22 whole-station source flora is matched to 2026 official-derived SR70 STATION CENTROIDS (35/39 exact names), not original surveyed plant plots or station edges.",
        "The 2023 OSM true rail/road graph gives geometric shortest path only, NOT service traffic, directional freight OD, active cargo load, slipstream wind exposure or actual propagule deposition.",
        "Five route stretches, 35 stations, 594 dependent site pairs in ONE author field cohort in ONE country. Three class outcomes are NOT three independent study effect sizes.",
        "Station size, plant regional source pool, land use, local management, detection and colonization history remain confounded.",
        "Subgroup Jaccard scales differ with per-station species richness/prevalence; apparent class differences cannot be interpreted as different per-seed transport efficiencies.",
        "No same-seed trap catch→germination→same rail adult→reproductive fitness link, and no phenology or genetic common garden results from this station study.",
        "Any stable signal is conditional observational prediction only, not evidence of selection for wind/mowing adaptation.",
      ],
      "not_meta_pooled_with_Slovenia_100m_original_source":True
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("s1","s2","pladias-diaspore","pladias-lifeform","sr70","osm-pbf","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    result=run(a.s1,a.s2,a.pladias_diaspore,a.pladias_lifeform,a.sr70,a.osm_pbf)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    compact={
        "complete":True,"version":result["version"],
        "source_v156_reference_verified":result["independently_recomputed_original_v156_whole_flora_road_rail_reference"],
        "group_primary":result["posthoc_source_class_rail_increment_comparison"],
        "group_sensitivity":{
          label:{
            mode:{
              "MSE":variant["original_stretch_heldout_CV"].get("equal_original_stretch_weighted_mse"),
              "rail_gain_pct":variant["original_stretch_heldout_CV"].get("rail_added_MSE_improvement_percentage_after_road"),
              "n_source_stretch_rail_gain_positive":variant["original_stretch_heldout_CV"].get("rail_improves_how_many_stretches"),
              "source_original_folds":variant["original_stretch_heldout_CV"].get("source_stretch_holdout")
            } for mode,variant in dat["class_specific_and_full_flora_richness_sensitivity"].items()
          } for label,dat in result["posthoc_explicit_alien_class_source_road_adjusted_predictive_results"].items()
        },
        "not_train_seed_causal_or_3_study_meta":True}
    print(json.dumps(compact,ensure_ascii=False))

if __name__=="__main__":main()

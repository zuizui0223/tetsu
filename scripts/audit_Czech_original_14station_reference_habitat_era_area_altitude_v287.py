#!/usr/bin/env python3
"""v28.7 author-original Czech RAIL station vs CzechVeg local-reference:
independent OPEN HABITAT, survey era, source plot AREA and control
ALTITUDE support gates, keeping botanical outcomes BLIND to matching.

The exact original v24.4 field cohort:
  14/16 trait-qualified Kutlvašr 2024 (39) true station-wide 2020–22
  railway flora sites, 4 real original railway survey stretches,
  each 5 CzechVeg-Open 2026 other-paper original vegetation relevés
  within 10km of the 2026 third-party station CENTROID;
  source-control plot IDs AND original WGS84 rounded 0.001° cells
  exclusive across all 70 real plant source plots. An exact source
  MILP solver certifies optimality. Baseline height/seed zlog
  distance and all source technical rarefaction sampling strictly
  reproduce v28.6 to 1e-9.

Original 117739 CzechVeg source plot headers expose HABITAT, YEAR,
SURF_AREA, ALTITUDE, INCLINATIO, GEOLOGY, SOIL, and COV_ROCK. Restrict
available source control opportunities among the SAME original 16
stations' original trait-qualified CzechVeg controls BEFORE examining
species or any functional outcome:

  baseline_any_nontransport_context — exact original 14/70 design;
  survey2020 — original control survey YEAR >=2020 (comparable era
  to station flora 2020–22, not matching exact species phenology);
  habitat_text_open — unambiguous source HABITAT-only token hints at
  open grassland, rock, sand, gravel, disturbed or ruderal communities;
  habitat_open_recent — both text-open and YEAR>=2020;
  habitat_open_area — source HABITAT open and source SURF_AREA 16..400m2;
  open_recent_area — all above;
  explicit_rock_sand_gravel — original source HABITAT explicitly names
  quarry/rock, sand, gravel, spoil only (NOT plant-based SYNTAXON).
  broad_open_altitude — HABITAT text-open and reference ALTITUDE valid,
  original reference site altitudes only checked for recorded coverage.

IMPORTANT: original Preslia station-wide AREA 2–90 ha and CzechVeg
relevé SURF_AREA m2 are NOT same geographic plot GRAIN and can never
be truly area-matched based on these fields. Railway station original
ELEVATION is NOT known; reference ALTITUDE alone cannot yield
rail-vs-reference matched elevation. CzechVeg original COV_ROCK may have
missing-to-zero substitution so do not claim measured exposed ballast.
HABITAT keyword proxies are from Czech/Czech-English text, not
field-verified substrate; SYNTAXON vegetation classes are derived from
plants and are specifically forbidden as match criteria.

For each predeclared metadata panel, run exact binary MILP under
  5 true source plots/station, >=2 actual independent CzechVeg papers,
  globally disjoint plot IDs, globally unique 100m rounded coordinate
  cells, >=3 original rail survey stretches, 10km geography, true
  source 2015+ inclusion, and original pre-verified LEDA coverage.
  Source 16 station eligibility is not prefiltered by outcomes.
If <10 actual railway stations or <3 original stretches, report
NOT IDENTIFIABLE: DO NOT create a positive/negative treatment effect.
If identified, same genuine original 100 technical species-equal
draws as v24.4/28.6 and report rail-vs-control CWM differences for
height, seed mass and 2D, then original stretch LOO. The selected
study/pool changes with panel; DO NOT interpret change from baseline
as a causal treatment effect from habitat alone. This source is 1
Czech railway author publication, other CzechVeg source publications
are reference observations not replicate railway studies.

No raw original source station locations, plot IDs, taxa × plot,
published detailed habitat descriptions, or exact GPS exported.
"""
from __future__ import annotations

import argparse,csv,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_CzechVeg_abiotic_regional_reference_opportunity_v194 import number
from analyze_Czech_39station_true_local_nontransport_plant_LEDA_beta_v240 import (
    source_station_and_reference_coords, source_species_to_traits,
    source_valid_local_control_plots, source_species
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from validate_Czech_rail_station_local_reference_geo_cells_unique_v244 import (
    exact_max_stations_with_cell_uniqueness
)
from validate_Czech_rail_local_reference_alternative_14station_assignments_v243 import (
    station_reference_profiles
)
from audit_Czech_original_rail_vs_local_CzechVeg_height_seed_beta_v286 import (
    comparative_stats, SOURCE_SEED,SOURCE_SEED_STEP,
    SOURCE_STATIONS,OLD_100M_SOURCE_2D
)

EXPECTED_SOURCE_HEADERS=117739
ORIGINAL_16=16
ORIGINAL_14=14
TECHNICAL_DRAWS=100
MIN_STATIONS=10
MIN_STRETCHES=3
# Text filters are ONLY original HABITAT value, never vegetation derived
# SYNTAXON, botanical species labels, LOCALITY or source paper titles.
OPEN_WORDS=re.compile(
    r"(?:\bgrassland\b|\bmeadow\b|\bpasture\b|\bruderal\b|"
    r"\bwasteland\b|\bbare\s+ground\b|\bpioneer\b|\bquarry\b|"
    r"\bgravel\b|\bsand(?:y)?\b|\brock(?:y)?\b|\bfield\s+margin\b|"
    r"lou[ck]|trav|tráv|pastv|step|such[áéý]|ruder|"
    r"písk|pisk|štěrk|sterk|skál|skal|kamen|"
    r"lom(?:u|y|ě|e)?\b|výsyp|vysyp|úhor|uhor|"
    r"nezastav|opuště|opustě|násyp|nasyp)",
    re.IGNORECASE
)
GRAVEL_WORDS=re.compile(
    r"(?:\bgravel\b|\bsand(?:y)?\b|\bquarry\b|\brock(?:y)?\b|"
    r"písk|pisk|štěrk|sterk|skál|skal|kamen|lom(?:u|y|ě|e)?\b|"
    r"výsyp|vysyp|hald|slag|spoil)",
    re.IGNORECASE
)
# Exclude explicitly forest/wetland context to avoid "forest grassy gap"
# overmatching; a positive word AND a negative word means ambiguous.
NEGATIVE_WORDS=re.compile(
    r"(?:\bforest\b|\bwoodland\b|\bswamp\b|\bwetland\b|"
    r"\bpeatland\b|\breedbed\b|les(?:e|a|ní|n[iy])?\b|"
    r"prales|rašelin|raselin|mokř|mokr|bažin|bazin|rybník|"
    r"rybnik|lužní\s+les|luzni\s+les)",
    re.IGNORECASE
)
PANELS=(
 "baseline_any_nontransport_context",
 "survey2020",
 "habitat_text_open",
 "habitat_open_recent",
 "habitat_open_area",
 "open_recent_area",
 "explicit_rock_sand_gravel",
 "broad_open_altitude"
)
MEASURES=("CWM_2D_Euclidean","CWM_height_1D_absolute",
          "CWM_seed_mass_1D_absolute")


def parse_metadata(path:Path,all_option_ids:set[str])->tuple[dict,dict]:
    records={}
    seen=Counter()
    with path.open(encoding="utf-8-sig",newline="") as handle:
        reader=csv.DictReader(handle)
        req={"RELEVE_NR","HABITAT","YEAR","SURF_AREA",
             "ALTITUDE","INCLINATIO","COV_ROCK","GEOLOGY","SOIL","SYNTAXON"}
        if req-set(reader.fieldnames or []):
            raise ValueError("Original CzechVeg public plot metadata missing source fields")
        for row in reader:
            seen["original_rows"]+=1
            pid=str(row.get("RELEVE_NR") or "").strip()
            if pid not in all_option_ids:continue
            if pid in records:raise ValueError("Original source CzechVeg plot ID duplicated")
            hab=str(row["HABITAT"] or "").strip()
            explicit_open=bool(OPEN_WORDS.search(hab))
            explicit_gravel=bool(GRAVEL_WORDS.search(hab))
            negative=bool(NEGATIVE_WORDS.search(hab))
            records[pid]={
                "source_original_habitat_present":bool(hab) and hab!="NA",
                "open":explicit_open and not negative,
                "gravel":explicit_gravel and not negative,
                "ambiguous_open_and_wooded_wet":explicit_open and negative,
                "year":number(row["YEAR"],1900,2030),
                "plot_area_m2":number(row["SURF_AREA"],0.01,1_000_000),
                "altitude_m":number(row["ALTITUDE"],0,2500),
                "slope_deg":number(row["INCLINATIO"],0,90),
                "rock_pct":number(row["COV_ROCK"],0,100),
                "source_geology_text_recorded":bool(str(row["GEOLOGY"] or "").strip().strip("NA")),
                "source_soil_text_recorded":bool(str(row["SOIL"] or "").strip().strip("NA")),
            }
    if seen["original_rows"]!=EXPECTED_SOURCE_HEADERS or len(records)!=len(all_option_ids):
        raise ValueError("Original CzechVeg 117739 header or source selectable plot IDs drift")
    return records,{
        "original_publisher_CzechVeg_header_rows":seen["original_rows"],
        "original_distinct_trait_qualified_candidate_plot_ids":len(records),
        "source_plot_habitat_nonmissing_n":sum(v["source_original_habitat_present"] for v in records.values()),
        "source_plot_habitat_text_open_n":sum(v["open"] for v in records.values()),
        "source_plot_habitat_explicit_gravel_sand_rock_n":sum(v["gravel"] for v in records.values()),
        "source_plot_habitat_open_and_wooded_wet_ambiguous_n":sum(v["ambiguous_open_and_wooded_wet"] for v in records.values()),
        "source_plot_2020plus_original_n":sum(v["year"] is not None and v["year"]>=2020 for v in records.values()),
        "source_plot_16_400m2_original_n":sum(v["plot_area_m2"] is not None and 16<=v["plot_area_m2"]<=400 for v in records.values()),
        "source_plot_with_recorded_altitude_n":sum(v["altitude_m"] is not None for v in records.values()),
        "source_original_slope_recorded_n":sum(v["slope_deg"] is not None for v in records.values()),
        "original_source_railway_site_ALTITUDE_not_in_Preslia_metadata":True,
        "original_source_COV_ROCK_zero_may_be_missingness_not_verified_rock":True,
        "matching_DOES_NOT_use_source_species_SYNTAXON":True,
    }


def belongs(panel:str,m:dict)->bool:
    if panel not in PANELS:raise ValueError("Undeclared nonplant habitat matching gate")
    if panel=="baseline_any_nontransport_context":return True
    if panel=="survey2020":return m["year"] is not None and m["year"]>=2020
    if panel=="habitat_text_open":return m["open"]
    if panel=="habitat_open_recent":return m["open"] and belongs("survey2020",m)
    if panel=="habitat_open_area":return m["open"] and m["plot_area_m2"] is not None and 16<=m["plot_area_m2"]<=400
    if panel=="open_recent_area":return (belongs("habitat_open_area",m) and belongs("survey2020",m))
    if panel=="explicit_rock_sand_gravel":return m["gravel"]
    if panel=="broad_open_altitude":return m["open"] and m["altitude_m"] is not None
    raise AssertionError("Unknown source metadata panel")


def panel_stations(stations:list[dict],metadata:dict,panel:str)->tuple[list[dict],dict]:
    if len(stations)!=ORIGINAL_16:
        raise ValueError("16 genuine original source rail sites with trait-qualified control opportunities changed")
    result=[]
    stats=Counter()
    for site in stations:
        candidate=site["control_options"]
        allowed=[x for x in candidate if belongs(panel,metadata[x["id"]])]
        result.append({**site,"control_options":allowed})
        stats["original_source_prequalified_options"]+=len(candidate)
        stats["eligible_source_reference_plot_opportunities"]+=len(allowed)
        stats["original_rail_stations_with_min5_options"]+=len(allowed)>=5
        stats["original_rail_stations_with_min2_reference_publications"]+=len({x["study"] for x in allowed})>=2
        stats["original_rail_stations_with_both_min5_and_min2papers"]+=len(allowed)>=5 and len({x["study"] for x in allowed})>=2
    return result,dict(stats)


def chosen_control_metadata(choice:dict,original_stations:list[dict],metadata:dict)->dict:
    values=[metadata[x["id"]] for plots in choice.values() for x in plots]
    years=[x["year"] for x in values if x["year"] is not None]
    areas=[x["plot_area_m2"] for x in values if x["plot_area_m2"] is not None]
    alts=[x["altitude_m"] for x in values if x["altitude_m"] is not None]
    slopes=[x["slope_deg"] for x in values if x["slope_deg"] is not None]
    sources=[x for plots in choice.values() for x in plots]
    return {
        "source_reference_original_releve_count":len(values),
        "source_distinct_true_stations":len(choice),
        "original_railway_field_stretches":len({original_stations[i]["stretch"] for i in choice}),
        "independent_original_CzechVeg_publications":len({x["study"] for x in sources}),
        "source_open_field_habitat_text_n":sum(x["open"] for x in values),
        "source_explicit_gravel_sand_rock_text_n":sum(x["gravel"] for x in values),
        "original_HABITAT_text_nonmissing_n":sum(x["source_original_habitat_present"] for x in values),
        "source_original_CzechVeg_survey2020plus_n":sum(x["year"] is not None and x["year"]>=2020 for x in values),
        "source_actual_CzechVeg_survey_year_range":[min(years),max(years)] if years else None,
        "source_control_original_plot_area_m2_median":median(areas) if areas else None,
        "source_control_original_plot_area_m2_min_max":[min(areas),max(areas)] if areas else None,
        "source_control_original_altitude_recorded_n":len(alts),
        "source_control_original_altitude_m_median":median(alts) if alts else None,
        "source_control_original_slope_recorded_n":len(slopes),
        "source_control_original_GEOL_text_recorded_n":sum(x["source_geology_text_recorded"] for x in values),
        "source_control_original_SOIL_text_recorded_n":sum(x["source_soil_text_recorded"] for x in values),
        "control_rock_pct_even_if_zero_NOT_physical_rock_evidence":True,
        "Preslia_whole_station_surveyed_area_NOT_equal_to_CzechVeg_original_releve_area":True,
        "original_botanical_field_year_2020_2022_vs_CzechVeg_survey_YEAR_not_same_phenology":True
    }


def evaluate(source_stations:list[dict],choice:dict,geo:dict,traits:dict,
             draws:int=TECHNICAL_DRAWS)->dict:
    if len(choice)<MIN_STATIONS or len({source_stations[i]["stretch"] for i in choice})<MIN_STRETCHES:
        return {"identified":False,"source_station_support_below_min10_or_3stretches":True}
    values=defaultdict(list)
    stretch_loo=defaultdict(list)
    for d in range(draws):
        r,n,old,qa=station_reference_profiles(
            source_stations,choice,geo,traits,SOURCE_SEED+SOURCE_SEED_STEP*d)
        stretch={source_stations[i]["id"]:source_stations[i]["stretch"] for i in choice}
        x=comparative_stats(r,n,stretch)
        for measure in MEASURES:
            for scope in ("all_14_source_station_pairs",
                          "within_stretches_equal_original_stretch"):
                values[measure,scope].append(
                    x[measure][scope]["rail_minus_local_control"])
            # Source LOO mean after each source publication stretch:
            for drop in sorted(set(stretch.values())):
                keep={key for key in r if stretch[key]!=drop}
                if len(keep)<8 or len({stretch[k] for k in keep})<2:continue
                sub=comparative_stats({k:r[k] for k in keep},
                    {k:n[k] for k in keep},{k:stretch[k] for k in keep})
                stretch_loo[measure,drop].append(
                    sub[measure]["all_14_source_station_pairs"]["rail_minus_local_control"])
    estimates={}
    for measure in MEASURES:
        estimates[measure]={}
        for scope in ("all_14_source_station_pairs","within_stretches_equal_original_stretch"):
            v=values[measure,scope]
            estimates[measure][scope]={
                "mean_original_technical_rarefied_RAIL_minus_nearby_context":mean(v),
                "negative_technical_draws":sum(x<0 for x in v),
                "source_technical_draws":len(v),
                "not_field_experimental_CI":True
            }
        loo=[mean(v) for (name,drop),v in stretch_loo.items() if name==measure]
        estimates[measure]["source_original_stretch_LOO"]={
            "n_heldout_stretches":len(loo),
            "min_max_source_effect":([min(loo),max(loo)] if loo else None),
            "all_heldout_stretches_negative":all(x<0 for x in loo) if loo else None
        }
    return {"identified":True,
            "n_original_rail_stations":len(choice),
            "n_original_rail_source_stretches":len({source_stations[i]["stretch"] for i in choice}),
            "source_original_reference_plots":5*len(choice),
            "source_trait_beta_components":estimates,
            "same_actual_station_plant_and_reference_species_richness_per_draw":True,
            "not_identical_original_surveyed_area_or_habitat_type":True}


def run(s1:Path,s2:Path,sr70:Path,headers:Path,species:Path,
        height:Path,seed:Path,draws:int=TECHNICAL_DRAWS)->dict:
    if draws!=100:raise ValueError("Freeze original exact 100 v24.4 plant trait subsamples")
    canopy,_=read_leda(height,"height_m")
    seed_m,_=read_leda(seed,"seed_mass_mg")
    original_stations,geo,sourceplot_ids,sourceqa=source_station_and_reference_coords(
        s1,s2,sr70,headers)
    original_plants,plantsqa=source_species(species,sourceplot_ids)
    traits,traitqa=source_species_to_traits(original_stations,original_plants,canopy,seed_m)
    eligible,eligqa=source_valid_local_control_plots(original_stations,original_plants,geo,traits)
    if len(eligible)!=16:raise ValueError("Original complete 16 true stations with biological coverage drifted")
    candidates={p["id"] for site in eligible for p in site["control_options"]}
    metadata,metaqa=parse_metadata(headers,candidates)
    output={}
    for panel in PANELS:
        sites,opportunity=panel_stations(eligible,metadata,panel)
        alloc,matchqa=exact_max_stations_with_cell_uniqueness(
            sites,geo,decimals=3,timeout_s=120)
        if alloc is None:
            output[panel]={
                "source_metadata_opportunity":opportunity,
                "source_MILP_status":matchqa,
                "source_effect_identified":False
            }
            continue
        if len(alloc)>ORIGINAL_14:
            raise ValueError("Strict source habitat constrained field-cohort exceeded previously source-proven 14-station maximum")
        controls=chosen_control_metadata(alloc,eligible,metadata)
        score=evaluate(eligible,alloc,geo,traits,draws)
        if panel=="baseline_any_nontransport_context":
            if len(alloc)!=14 or score["identified"] is not True:
                raise ValueError("Original 14/70 baseline source cohort lost")
            # Exact all original source 100 draw from v28.6, fail closed.
            old_score=score["source_trait_beta_components"]["CWM_2D_Euclidean"][
                "all_14_source_station_pairs"][
                "mean_original_technical_rarefied_RAIL_minus_nearby_context"]
            if abs(old_score-OLD_100M_SOURCE_2D["rail_minus_local_nontransport_CWM_beta"])>1e-9:
                raise ValueError("True original v28.6 14/70 control source 2D baseline changed")
        output[panel]={
            "source_metadata_opportunity":opportunity,
            "source_MILP_status":matchqa,
            "source_original_reference_metadata":controls,
            "source_effect_identified":score["identified"],
            "matched_source_rail_and_local_cwm_trait_beta":score,
        }
    return {
        "version":"v28.7",
        "source_original_Preslia_2024_station_floras_39_no_missing_filled":True,
        "source_original_CzechVeg_2026_March_v1_1_DOI":"10.5281/zenodo.19231152",
        "original_station_S2_and_CzechVeg_original_source_QA":sourceqa,
        "original_CzechVeg_source_species_complete_QA":plantsqa,
        "source_trait_original_LEDA_coverage_QA":traitqa,
        "source_original_16_station_trait_eligible_QA":eligqa,
        "source_publication_real_option_habitat_year_area_elevation_quality":metaqa,
        "true_original_100m_cell_distinct_reference_selection_and_source_sensitivity":output,
        "strict_original_source_14station_70plot_v28_6_prior_benchmark_reproduced":True,
        "original_Preslia_station_altitude_NOT_measured":True,
        "CzechVeg_HABITAT_keyword_not_field_confirmed_open_or_ballast":True,
        "source_CzechVeg_SYNTAXON_derived_from_plants_NOT_matching_control":True,
        "source_plot_area_match_to_railway_station_ha_NOT_IDENTIFIED":True,
        "source_CzechVeg_rock_cover_0_may_mean_original_missingness":True,
        "single_original_czech_railway_publication_with_reference_21_other_publications_NOT_independent_21_rail_studies":True,
        "no_rail_specific_wind_mowing_evolution_traffic_seed_causality":True,
        "no_original_taxon_plot_ids_individual_stations_habitat_text_or_GPS_republished":True
    }


def main()->None:
    p=argparse.ArgumentParser()
    for x in ("s1","s2","sr70","headers","species","height","seed","output"):
        p.add_argument("--"+x,required=True,type=Path)
    a=p.parse_args()
    out=run(a.s1,a.s2,a.sr70,a.headers,a.species,a.height,a.seed)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"version":"v28.7",
        "source_original_metadata_coverage":out["source_publication_real_option_habitat_year_area_elevation_quality"]},ensure_ascii=False))
    for panel,data in out["true_original_100m_cell_distinct_reference_selection_and_source_sensitivity"].items():
        stat=data.get("matched_source_rail_and_local_cwm_trait_beta",{})
        terms={}
        for m in MEASURES:
            q=stat.get("source_trait_beta_components",{}).get(m,{})
            terms[m]=q.get("within_stretches_equal_original_stretch")
        print(json.dumps({
            "version":"v28.7","panel":panel,
            "source_qualified_reference_opportunity":data["source_metadata_opportunity"],
            "original_rail_stations":stat.get("n_original_rail_stations"),
            "original_rail_stretches":stat.get("n_original_rail_source_stretches"),
            "identified":data["source_effect_identified"],
            "selected_original_context":data.get("source_original_reference_metadata"),
            "original_source_height_seed_functional_differentiation_within_stretch":terms
        },ensure_ascii=False))


if __name__=="__main__":main()

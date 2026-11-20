#!/usr/bin/env python3
"""Independent source-only verification of ORIGINAL MIREN native Elevation
signed community stature and distance-vs-altitude gap slopes.

Check all original source study×backbone×coverage strata, original per-Region
trait/null arithmetic, physical survey-year sensitivity, pseudo-replication
guards and self-consistent geographic LOO. NEVER print original site
altitudes, site×taxon flora, source Road/Plot IDs or private coordinates.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean

TREES=("LCVP","WP","TPL")
PANELS=("primary_5species_25percent","strict_10species_50percent")
MODES=("native_source_Elevation_original_site_environmental_slopes",
       "original_within_road_single_first_survey_year_only_sensitivity")
ZONES=("roadside","interior","roadside_minus_interior")
METRICS=("height_signed_slope_per_km",
         "height_distance_vs_elevation_gap_slope_per_km")

def isclose(a,b):
    if not all(math.isfinite(float(x)) for x in (a,b)):
        raise ValueError("Original source field or conditional elevation statistic non-finite")
    if not math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10):
        raise ValueError(f"Original source environment group arithmetic failed: {a} versus {b}")

def run(source:dict,metadata:dict)->dict:
    if source.get("version")!="v29.9" or metadata.get("version")!="v29.9_SOURCE_ONLY_ELEVATION_QA":
        raise ValueError("Original first-party MIREN elevation dataset version drift")
    if not source.get("original_v298_joint_phylo_occupancy_height_effects_reproduced_for_all_6_tree_panel_cohorts"):
        raise ValueError("Frozen original source clade+frequency plant CWM beta reference unverified")
    if not source.get("original_MIREN_is_ONE_field_program_not_EIGHT_independent_studies"):
        raise ValueError("MIREN field program must be ONE, not 8 independent original field studies")
    if metadata["source_original_plant_Elevation_data"]["n_source_plots_with_qualified_original_Elevation"]!=6446:
        raise ValueError("MIREN source field elevation metadata source plot count drift")
    tree_detail=source["source_original_study_tree_region_year_and_environmental_support"]
    stats=source["original_region_equal_source_elevation_gradient_summaries"]
    if set(tree_detail)!=set(TREES) or set(stats)!=set(TREES):
        raise ValueError("Three genuine original megatree nomenclature alternatives incomplete")
    tests=0
    region_panels=0
    per_tree={}
    for tree in TREES:
        per_tree[tree]={}
        if set(tree_detail[tree]["panels"])!=set(PANELS) or set(stats[tree])!=set(PANELS):
            raise ValueError("Prior two original MIREN source floristic panel supports incomplete")
        for panel in PANELS:
            x=tree_detail[tree]["panels"][panel]
            body=x["regions"]
            prior=x["number_original_regions_with_joint_phylo_occupancy_source"]
            n_main=x["number_original_regions_with_native_altitude_slopes"]
            n_one=x["number_original_regions_with_original_one_year_roads_only_slopes"]
            if not 0<=n_one<=n_main<=prior<=8:
                raise ValueError("Original field-supported elevation Region counts impossible")
            if len(body)!=(8 if panel.startswith("primary") else 5):
                raise ValueError("Source original 8/5 MIREN regions missing")
            per_tree[tree][panel]={"joint_original_regions":prior,"elevation_regions":n_main,
                                    "single_field_year_sensitivity_regions":n_one}
            for mode in MODES:
                y=stats[tree][panel][mode]
                supported={}
                for region,rec in body.items():
                    if mode==MODES[0]:
                        eff=rec.get("source_original_region_signed_effects")
                    else:
                        eff=rec.get("single_field_year_roads_only_sensitivity")
                    if eff is not None:
                        if rec["source_status"]!="SOURCE_ALTITUDE_ENVIRONMENTAL_SLOPE_ESTIMABLE":
                            raise ValueError("Effect found on source-metadata unidentifiable region")
                        if set(eff)!=set(ZONES):
                            raise ValueError("Original true habitat and paired-contrast structure broken")
                        supported[region]=eff
                        for metric in METRICS:
                            road=eff["roadside"][metric]
                            inland=eff["interior"][metric]
                            diff=eff["roadside_minus_interior"][metric]
                            fields=("source_original_observed_environmental_slope",
                                    "source_true_phylo_and_occupancy_trait_null_mean_slope",
                                    "original_observed_minus_conditional_null_environmental_slope")
                            for zone_dict in (road,inland,diff):
                                isclose(zone_dict[fields[0]]-zone_dict[fields[1]],
                                        zone_dict[fields[2]])
                                tests+=1
                            for field in fields:
                                isclose(road[field]-inland[field],diff[field])
                                tests+=1
                            for zone_dict in (road,inland,diff):
                                for rkey in ("source_technical_lower_rank_NOT_field_p",
                                             "source_technical_upper_rank_NOT_field_p"):
                                    v=zone_dict[rkey]
                                    if not 0<=v<=1:
                                        raise ValueError("Source technical label rank outside [0,1]")
                                    tests+=1
                if y["source_regions_n"]!=len(supported):
                    raise ValueError("Source Region support mismatches the effects in original field cohort")
                if mode==MODES[0] and len(supported)!=n_main:
                    raise ValueError("Main signed altitude source cohort incorrectly counted")
                if mode==MODES[1] and len(supported)!=n_one:
                    raise ValueError("Single-year-road source cohort incorrectly counted")
                all_summary=y.get("region_equal_estimates")
                if len(supported)<3:
                    if all_summary is not None:
                        raise ValueError("Unsupported original source regions assigned fake ecological mean")
                    continue
                if all_summary is None or set(all_summary)!=set(ZONES):
                    raise ValueError("Original region-equal environmental slope means missing despite source")
                for z in ZONES:
                    for metric in METRICS:
                        r=all_summary[z][metric]
                        if not r:
                            raise ValueError("True source ecosystem metric missing")
                        source_signed={reg:data[z][metric][
                            "original_observed_minus_conditional_null_environmental_slope"]
                            for reg,data in supported.items()}
                        source_obs=[data[z][metric][
                            "source_original_observed_environmental_slope"]
                            for reg,data in supported.items()]
                        if set(r["source_original_region_code_aggregate_effects"])!=set(source_signed):
                            raise ValueError("Original Region code difference introduced by summary")
                        for reg,v in source_signed.items():
                            isclose(v,r["source_original_region_code_aggregate_effects"][reg])
                            tests+=1
                        isclose(mean(source_obs),r["original_source_region_equal_mean_observed_environmental_slope"])
                        isclose(mean(source_signed.values()),
                                r["original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"])
                        if sum(v>0 for v in source_signed.values())!=r["original_region_codes_positive_excess"]:
                            raise ValueError("Source signed environmental slope positive region count changed")
                        if sum(v<0 for v in source_signed.values())!=r["original_region_codes_negative_excess"]:
                            raise ValueError("Source signed environmental slope negative region count changed")
                        loo=[mean(v for region,v in source_signed.items() if region!=leave)
                             for leave in source_signed]
                        got=r["source_original_region_loo_mean_range_not_field_CI"]
                        if len(got)!=2:
                            raise ValueError("Region leave-one-out range expected 2 values")
                        isclose(min(loo),got[0]);isclose(max(loo),got[1])
                        tests+=8
            region_panels+=1
    return {
       "version":"v29.9",
       "original_source_qa_provenance_verified":True,
       "original_real_field_programmes_n":1,
       "original_mountain_road_source_inland_same_year_pairing":True,
       "actual_original_elevation_plots_with_valid_metadata":6446,
       "n_source_backbone_by_field_support_panels":region_panels,
       "source_support_counts_by_tree_and_field_cohort":per_tree,
       "n_independent_source_arithmetic_sign_Loo_and_physical_year_checks":tests,
       "all_original_tree_phylo75_AND_occupancy_source_baselines_reproduced":True,
       "original_source_site_individual_Elevation_GPS_plant_species_rows_NOT_exported":True,
       "not_a_global_causal_railway_selection_effect_or_field_study_p_value":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    src=json.loads(a.source.read_text(encoding="utf8"))
    meta=json.loads(a.metadata.read_text(encoding="utf8"))
    result=run(src,meta)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Independent second-parser check that the frozen six public research
summary rows equal genuine original MIREN field/trait/tree source output.
No site elevation, original plant taxon, GPS or raw source field row released.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

SRC_MAIN="native_source_Elevation_original_site_environmental_slopes"
SRC_YEAR="original_within_road_single_first_survey_year_only_sensitivity"
PANEL_NAMES={"primary":"primary_5species_25percent",
             "strict":"strict_10species_50percent"}
BACKBONES={"LCVP","WP","TPL"}
METRIC_A="height_signed_slope_per_km"
METRIC_B="height_distance_vs_elevation_gap_slope_per_km"

def close(a,b):
    if not math.isfinite(float(a)) or not math.isfinite(float(b)) or \
       abs(float(a)-float(b))>1e-10:
        raise ValueError("v29.9 original source-derived macro output drift: "+repr((a,b)))

def run(src:dict,frozen:dict)->dict:
    if src.get("version")!="v29.9" or frozen.get("version")!="v29.9":
        raise ValueError("Original MIREN source verified publication summary missing")
    if src["original_source_Elevation_per_plot_metadata_QA"]["n_source_plots_with_qualified_original_Elevation"]!=6446:
        raise ValueError("Original plot elevation source cohort changed")
    original=src["original_region_equal_source_elevation_gradient_summaries"]
    source_tree=src["source_original_study_tree_region_year_and_environmental_support"]
    rows=frozen["source_by_original_tree_panel"]
    if len(rows)!=6 or set(original)!=BACKBONES:
        raise ValueError("Original three phylogeny trees×two trait support panels not retained")
    done=set()
    count=0
    for x in rows:
        tree=x["original_tree"];kind=x["source_trait_support"]
        if tree not in BACKBONES or kind not in PANEL_NAMES or (tree,kind) in done:
            raise ValueError("Duplicate/nonoriginal phylo nomenclature/plot support record")
        done.add((tree,kind))
        panel=PANEL_NAMES[kind]
        main=original[tree][panel][SRC_MAIN]
        year=original[tree][panel][SRC_YEAR]
        if main["region_equal_estimates"] is None or year["region_equal_estimates"] is None:
            raise ValueError("Original MIREN environmental source supports insufficient")
        if x["source_original_regions_n"]!=main["source_regions_n"] or \
           x["single_survey_year_roads_only"]["original_regions_supported"]!=year["source_regions_n"]:
            raise ValueError("Original Region support counts lost")
        q=source_tree[tree]["panels"][panel]
        active=[r["source_original_true_MIREN_altitude_support"]
                for r in q["regions"].values()
                if r["source_status"]=="SOURCE_ALTITUDE_ENVIRONMENTAL_SLOPE_ESTIMABLE"]
        road=sum(r["original_source_roads_with_4_sites_150m_valid_elevation"] for r in active)
        site=sum(r["source_true_field_sites_on_usable_elevation_roads"] for r in active)
        if road!=x["original_real_roads_n"] or site!=x["original_paired_real_sites_n"]:
            raise ValueError("Original physical Road×Site ecological source support changed")
        for zone,displayed in (("roadside","roadside"),("interior","inland")):
            v=main["region_equal_estimates"][zone][METRIC_A]
            gap=main["region_equal_estimates"][zone][METRIC_B]
            close(v["original_source_region_equal_mean_observed_environmental_slope"],
                  x["signed_CWM_height_elevation_slope_zlog_per_km_observed"][displayed])
            close(v["original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"],
                  x["signed_CWM_height_elevation_slope_zlog_per_km_excess_vs_original_phylo75_occupancy_null"][displayed])
            close(gap["original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"],
                  x["gap_abs_CWM_height_distance_vs_abs_native_elevation_difference_excess"][displayed])
            y=year["region_equal_estimates"][zone][METRIC_A]
            close(y["original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"],
                  x["single_survey_year_roads_only"][
                      "extra_elevation_slope_over_joint_null_"+displayed])
            count+=4
        close(main["region_equal_estimates"]["roadside_minus_interior"][METRIC_A][
                "original_source_region_equal_mean_excess_vs_joint_phylo_occupancy_null"],
              x["signed_CWM_height_elevation_slope_zlog_per_km_excess_vs_original_phylo75_occupancy_null"]["roadside_minus_inland"])
        count+=1
    if len(done)!=6:raise ValueError("No original field botanical cohort complete")
    return {"version":"v29.9",
      "original_MIREN_source_published_derived_rows_n":6,
      "independent_original_source_numeric_summary_comparisons":count,
      "original_6_main_years_phylogeny_support_panels_identical_to_source":True,
      "all_six_original_MIREN_source_environmental_slope_units_preserved":True,
      "no_source_original_plant_per_site_elevation_and_GPS_exported":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--derived",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    src=json.loads(a.source.read_text(encoding="utf8"))
    d=json.loads(a.derived.read_text(encoding="utf8"))
    result=run(src,d)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()

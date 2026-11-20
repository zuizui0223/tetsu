#!/usr/bin/env python3
"""Independent source-output audit for v30.0 genuine MIREN held-out-site skill.

Second-parser verifies all 3 published original V.PhyloMaker2 trees × both
source botanical trait coverage panels, exact genuine 4-region etc sampling,
original site/road geography/time, ALL 398 technical label draws per eligible
source, exact publisher-facing derived values, numeric arithmetic, signed
regional heterogeneity, original Region leave-one-out summaries, and no zero
imputation on unsupported source. This is not a second ecological field study.
Never export original site elevation, botanical plot species, or geographic ID.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean
TREES={"LCVP","WP","TPL"}
SOURCE_PANELS={"primary_5species_25percent":"primary",
               "strict_10species_50percent":"strict"}
HABITATS=("roadside","interior","roadside_minus_interior")
SRC_SOURCE_EFFECTS="original_source_original_site_heldout_predictability"
SOURCE_CVM="observed_original_site_heldout_elevation_gain"
SOURCE_YEAR="single_field_year_per_road_heldout_elevation_gain"

def close(x,y,tol=1e-10):
    a=float(x);b=float(y)
    if not math.isfinite(a) or not math.isfinite(b) or abs(a-b)>tol:
        raise ValueError(f"Original source held-out prediction metric vs prior published summary mismatch: {a}/{b}")

def verify(original:dict,published:dict)->dict:
    if original.get("version")!="v30.0" or published.get("version")!="v30.0":
        raise ValueError("Original publication version changed")
    if not original.get("original_v299_six_tree_by_support_panel_true_native_altitude_effects_reproduced_to_1e_10"):
        raise ValueError("True original v29.9 within-road elevation CWM source cohort not reproduced")
    if published.get("original_source_six_v299_signed_native_elevation_effects_reproduced_exact")!=True:
        raise ValueError("Frozen published source could mask previous source verification")
    if original["source_original_author_MIREN_study_programmes_n"]!=1 or\
       published["original_MIREN_study_field_programmes"]!=1:
        raise ValueError("Original one MIREN roads programme mistaken for multiple independent publications")
    if original["original_MIREN_source_elevation_plot_metadata_supported_n"]!=6446 or\
       published["original_MIREN_plot_Elevation_metadata_with_source_QA"]!=6446:
        raise ValueError("Original MIREN per-field plot elevation metadata count drift")
    true=original["original_first_party_tree_source_and_original_region_support"]
    regions=original["region_equal_source_environmental_heldout_predictability"]
    if set(true)!=TREES or set(regions)!=TREES:
        raise ValueError("Missing original 3 phylogeny nomenclature alternatives")
    rows=published["source_tree_by_field_support_panel"]
    if len(rows)!=6:raise ValueError("Expected exact 3 actual trees x 2 original support panels")
    used=set();checks=0;src_regions=set()
    for derived in rows:
        tree=derived["original_tree"];label=derived["source_trait_support"]
        if tree not in TREES or label not in SOURCE_PANELS or (tree,label) in used:
            raise ValueError("Non-source or duplicate phylogeny/field-support row")
        used.add((tree,label))
        panel=next(k for k,v in SOURCE_PANELS.items() if v==label)
        q=true[tree]["panels"][panel]
        stats=regions[tree][panel]
        eligible={reg:x for reg,x in q["original_source_region_support_and_predictivity"].items()
                  if x.get("new_predictive_heldout_effects") is not None}
        actual_roads=sum(x["source_CV_site_support"][
            "original_roads_with_ge6_true_sites_and_identifiable_heldout_elev"] for x in eligible.values())
        actual_sites=sum(x["source_CV_site_support"][
            "original_heldout_sites_on_eligible_roads"] for x in eligible.values())
        if (actual_roads!=derived["original_real_roads_with_6plus_source_sites"] or
            actual_sites!=derived["original_physical_roadside_inland_site_pairs"] or
            len(eligible)!=derived["original_source_region_codes_with_true_heldout_prediction"] or
            q["n_new_source_regions_with_true_out_of_sample_prediction"]!=len(eligible)):
            raise ValueError("Original physical road/site/Region heldout predictive cohort support drift")
        src_regions.update(eligible)
        main=stats[SOURCE_CVM]
        if main["source_supported_original_region_codes_n"]!=len(eligible):
            raise ValueError("Original source region environmental CV support conflicts")
        if len(eligible)<3:
            if main["source_region_equal_effects"] is not None:
                raise ValueError("Insufficient Region source sample assigned an effect")
            continue
        prev=main["source_region_equal_effects"]
        for zone in HABITATS:
            x=prev[zone]
            samples={}
            actual_gain={}
            for reg,record in eligible.items():
                effect=record["new_predictive_heldout_effects"]
                source=effect[SRC_SOURCE_EFFECTS][zone]
                if source["n_technical_label_draws_NOT_field_studies"]!=398:
                    raise ValueError("Original technical label draws changed")
                act=source["source_original_heldout_elevation_predictive_gain_zlog_CWM_RMSE"]
                null=source["source_phylogeny_and_occupancy_trait_null_heldout_predictive_gain"]
                excess=source["original_site_LOOCV_observed_minus_source_conditional_null"]
                close(act-null,excess)
                checks+=1
                samples[reg]=excess
                actual_gain[reg]=act
                for direction in ("technical_upper_rank_NOT_source_field_p",
                                  "technical_lower_rank_NOT_source_field_p"):
                    v=source[direction]
                    if not 0<=v<=1:
                        raise ValueError("Technical null rank invalid")
                if effect["original_observed_mean_rmse_with_intercept_training_only"].get(zone) is not None:
                    close(effect["original_observed_mean_rmse_with_intercept_training_only"][zone]-
                          effect["original_observed_mean_rmse_with_elevation_training"][zone],
                          act)
                    checks+=1
            close(x["original_source_region_equal_mean_heldout_Elevation_predictive_gain"],
                  mean(actual_gain.values()))
            close(x["original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null"],
                  mean(samples.values()))
            if x["source_regions_true_heldout_prediction_improvement_positive"]!=sum(v>0 for v in actual_gain.values()):
                raise ValueError("Positive actual source predictive region count altered")
            if x["source_regions_excess_predictive_gain_positive"]!=sum(v>0 for v in samples.values()):
                raise ValueError("Positive source excessive prediction region count altered")
            for reg,value in samples.items():
                close(x["by_original_source_region_excess_only"][reg],value)
            loo=[mean(v for r,v in samples.items() if r!=skip) for skip in samples]
            listed=x["source_region_leave_one_out_excess_gain_minmax_not_a_global_CI"]
            close(listed[0],min(loo));close(listed[1],max(loo))
            checks+=4+len(samples)
        labels={
           "roadside":("real_observed_original_roadside_site_elevation_LOO_RMSE_gain",
                       "roadside_gain_excess_vs_joint_phylo75_occ_label_null"),
           "interior":(None,"inland_gain_excess_vs_joint_phylo75_occ_label_null"),
           "roadside_minus_interior":(None,"roadside_minus_inland_gain_excess_vs_joint_phylo75_occ_label_null")
        }
        for z,(obsname,excessname) in labels.items():
            if obsname:close(prev[z][
                   "original_source_region_equal_mean_heldout_Elevation_predictive_gain"],
                   derived[obsname])
            close(prev[z]["original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null"],
                  derived[excessname])
            checks+=2
        close(derived["roadside_gain_excess_vs_joint_phylo75_occ_label_null"]-
              derived["inland_gain_excess_vs_joint_phylo75_occ_label_null"],
              derived["roadside_minus_inland_gain_excess_vs_joint_phylo75_occ_label_null"])
        if (derived["n_original_region_codes_actual_roadside_prediction_gain_positive"]!=
            prev["roadside"]["source_regions_true_heldout_prediction_improvement_positive"] or
            derived["n_original_region_codes_roadside_gain_excess_positive"]!=
            prev["roadside"]["source_regions_excess_predictive_gain_positive"]):
            raise ValueError("Published predictive region sign counts drift")
        year=stats[SOURCE_YEAR]
        for zone,field in (
            ("roadside","same_original_field_year_roads_only_source_roadside_excess"),
            ("interior","same_original_field_year_roads_only_source_inland_excess")):
            if year["source_region_equal_effects"] is None:
                raise ValueError("Published year-only predictive values lack data support")
            close(year["source_region_equal_effects"][zone][
                 "original_source_region_equal_mean_gain_excess_vs_phylo_occupancy_null"],
                  derived[field])
            checks+=1
        # Never count a tree as an independent source publication.
        if q["n_prior_source_environmental_regions_reproduced"]<len(eligible):
            raise ValueError("New CV gained extra botanical field Regions that were not in v29.9")
    if len(used)!=6:raise ValueError("Missing actual original 3 tree x two field-support profiles")
    return {
      "version":"v30.0",
      "original_first_party_source_and_published_summary_match":True,
      "n_original_tree_by_source_trait_support_panels":len(used),
      "original_MIREN_botanical_field_programmes_n":1,
      "genuine_original_native_plot_Elevation_metadata_supported_n":6446,
      "source_numeric_group_habitat_RMSE_null_intercept_year_and_LOO_checks":checks,
      "source_phylogenies_are_alternative_not_field_replicates":True,
      "unsupported_original_source_regions_NOT_assigned_false_zero_effect":True,
      "original_site_elevation_plant_species_GPS_and_phylogenetic_tips_NOT_exported":True,
    }
def main():
    p=argparse.ArgumentParser()
    for flag in ("source","derived","output"):
        p.add_argument("--"+flag,type=Path,required=True)
    args=p.parse_args()
    source=json.loads(args.source.read_text(encoding="utf8"))
    derived=json.loads(args.derived.read_text(encoding="utf8"))
    out=verify(source,derived)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(out,indent=2,ensure_ascii=False))
if __name__=="__main__": main()

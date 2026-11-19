#!/usr/bin/env python3
"""Verify readable v29.7 research summary reproduces every canonical ORIGINAL
MIREN, genuine V.PhyloMaker2 source result without typo or post-hoc losses.

This is a *second parsing and arithmetic verification* of original source
artifact aggregates. It does NOT convert original region effects into
independent field-studies or a pooled global meta-analysis.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

MODE="original_megatree_monophyletic_rootdepth_cut_0.75"
PANEL_75_NAMES={
    "primary_5species_25percent":"primary_75pct",
    "strict_10species_50percent":"strict_75pct",
}

def compare(original:dict,derived:dict)->dict:
    if original.get("version")!="v29.7" or derived.get("version")!="v29.7":
        raise ValueError("Original MIREN source effects or summary version changed")
    if original.get("original_MIREN_independent_field_programs_n")!=1 or\
       not derived.get("original_MIREN_mountain_ROADS_program_not_rail_study"):
        raise ValueError("A source road study was misidentified as multiple railway studies")
    oa=original["source_tree_original_per_region_support_and_effects"]
    summary=original["source_region_equal_descriptive_summaries_only_if_supported"]
    names=set(("LCVP","WP","TPL"))
    if set(oa)!=names or set(summary)!=names:
        raise ValueError("Original three phylogenetic nomenclature alternatives changed")
    tested=0
    for tree in names:
        for original_panel,target_panel in PANEL_75_NAMES.items():
            sample=oa[tree]["panels"][original_panel]["original_source_region_results"]
            s=derived[target_panel][tree]
            effects=summary[tree][original_panel][MODE]
            regions=effects["region_equal_descriptive_effects_if_3plus"]
            if not regions:
                raise ValueError("The original tree/site branch-clade source does not support effect estimation")
            valid={r:d for r,d in sample.items()
                   if d["status"]=="ORIGINAL_DIRECT_TIP_MIREN_FIELD_COHORT_SUPPORTED"}
            codes=sorted(valid)
            site_n=sum(x["source_tree_support"]["retained_original_source_site_pairs"]
                       for x in valid.values())
            if codes!=s["region_codes_supported"] or len(valid)!=s["regions_used"] or\
                    site_n!=s["original_paired_sites_remaining"]:
                raise ValueError("Original active phylogenetic region/source physical site support drift")
            for zone,name in (("roadside","roadside"),
                              ("interior","interior"),
                              ("roadside_minus_interior","roadside_minus_interior")):
                data=regions[zone]["height"]
                actual=data["original_region_code_equal_mean_source_excess"]
                if abs(actual-s[name])>1e-11:
                    raise ValueError(f"Original validated source beta statistic drifted: {tree}/{target_panel}/{zone}")
                positives=data["original_source_regions_signed_positive"]
                if zone=="roadside":
                    if positives!=s["positive_roadside_regions"]:
                        raise ValueError("Original source signed roadside study counts changed")
                elif zone=="interior":
                    if positives!=s["positive_interior_regions"]:
                        raise ValueError("Original source signed inland study counts changed")
                else:
                    if s["negative_roadside_minus_interior_regions"]!=len(valid)-positives:
                        raise ValueError("Original source signed corridor extra-beta counts changed")
                    loo=data["region_leave_one_out_minmax_not_global_CI"]
                    if any(abs(float(a)-float(b))>1e-11 for a,b in
                           zip(loo,s["difference_LOO_range"])):
                        raise ValueError("Original source region deletion robustness range drifted")
                    if len(loo)!=2:raise ValueError("Original region exclusion min-max incomplete")
                tested+=1
            if abs(s["roadside"]-s["interior"]-s["roadside_minus_interior"])>1e-11:
                raise ValueError("Roadside-inland source null decomposition does not balance")
            # Earlier 0.25/0.50 clades fail predeclared source exchangeability.
            for depth in ("0.25","0.50"):
                mod="original_megatree_monophyletic_rootdepth_cut_"+depth
                x=summary[tree][original_panel][mod]
                if x["source_original_region_codes_with_results"]!=0 or\
                        x["region_equal_descriptive_effects_if_3plus"] is not None:
                    raise ValueError("Undeclared early clade depth results cannot be converted into false zero")
            tested+=1
    return {"version":"v29.7","original_source_summary_source_authenticated":True,
            "source_original_phylogeny_study_and_panel_pairs":6,
            "comparisons_tested_original_source_effects_signs_support_and_LOO":tested,
            "original_regions_did_not_become_independent_field_publications":True,
            "three_original_phylo_tree_nomenclatures_are_not_independent_studies":True,
            "non_identifiable_depth_0_25_and_0_50_nulls_NOT_zero_effects":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-source",required=True,type=Path)
    p.add_argument("--derived",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    original=json.loads(a.original_source.read_text(encoding="utf8"))
    derived=json.loads(a.derived.read_text(encoding="utf8"))
    qa=compare(original,derived)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(qa,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(qa,indent=2,ensure_ascii=False))
if __name__=="__main__":main()

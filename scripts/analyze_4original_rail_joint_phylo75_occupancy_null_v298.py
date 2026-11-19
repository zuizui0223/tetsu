#!/usr/bin/env python3
"""v29.8 4 genuine European RAIL field publications:
TRUE V.PhyloMaker2 root-depth 75% clades × same original species source
occupancy QUINTILES, simultaneously preserved before label swapping.

The original source site-by-species matrix, source field geography and paired
canopy height+seed mass per species NEVER change. No source taxon imputation.
All 12 v29.6 SOURCE exact-tip phylogeny-only original effects must be fully
reproduced to 1e-10 or new simultaneous null is NOT reported.

This analysis is NOT a fifth field paper and is NOT a causal railway effect.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from statistics import mean
import numpy as np
from audit_VPhyloMaker2_RData_python_no_R_environment_v296 import (
  inspect_tree_rdata,TREES,COMMIT
)
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import (
 STUDIES,SITES,EXPECTED_TIP_COUNTS,EXPECTED_SOURCE_DIRECT_TIP_HEIGHT,
 original_source_private_csv,source_study_matrix,original_root_depths,
 original_true_clade_groups,source_null,source_json_native
)
from phylogeny_occupancy_joint_label_null_v298 import (
  original_global_taxon_occupancy_quintiles,
  joint_true_clade_and_occupancy_groups,
  original_source_joint_exchangeability,
)
DRAWS_PER_CHAIN=199
REFERENCE_PHYLO_MODE_INDEX=3
NEW_JOINT_PHYLO_AND_OCCUPANCY_MODE_INDEX=4
REF_GENUINE_V296={
 "LCVP":{
   "CZ_Czechia":0.030777984357754737,
   "IT_Italy_Parma":0.03541281361924449,
   "HR_BA_2017_one_study":0.09194637808943475,
   "SI_Slovenia_100m_track":0.06864392609187583},
 "WP":{
   "CZ_Czechia":0.02267026897853258,
   "IT_Italy_Parma":0.031800961395502314,
   "HR_BA_2017_one_study":0.10519811633404898,
   "SI_Slovenia_100m_track":0.06204654920215141},
 "TPL":{
   "CZ_Czechia":0.034796332075753736,
   "IT_Italy_Parma":0.035237121740295196,
   "HR_BA_2017_one_study":0.10197649659671024,
   "SI_Slovenia_100m_track":0.06398028682522175},
}

def run(original_private_csv:Path, tree_dir:Path, output:Path,draws=DRAWS_PER_CHAIN)->dict:
    if draws!=DRAWS_PER_CHAIN:
        raise ValueError("Do not alter published-source v29.6 reference technical draws when comparing new null")
    source=original_source_private_csv(original_private_csv)
    full={}
    bytree={}
    for tree_ix,tree in enumerate(TREES):
        obj,sourceqa=inspect_tree_rdata(
            tree_dir/("GBOTB.extended."+tree+".rda"),tree)
        tips=[str(v).replace("_"," ") for v in np.asarray(obj["tip.label"]).reshape(-1)]
        if len(set(tips))!=len(tips):raise ValueError("Original source phylo duplicated normalized species tips")
        lookup={s:i+1 for i,s in enumerate(tips)}
        par,depth,tree_depth=original_root_depths(
            np.asarray(obj["edge"]),np.asarray(obj["edge.length"]),
            len(tips),int(np.asarray(obj["Nnode"]).reshape(-1)[0]))
        backbones={}
        for study_ix,study in enumerate(STUDIES):
            P,T,names,geo,info=source_study_matrix(source,study,lookup)
            sourceinfo=info["qa"]
            if (sourceinfo["source_original_exact_megatree_tip_species_n"]!=
                EXPECTED_TIP_COUNTS[tree][study_ix]):
                raise ValueError("Genuine original European rail study phylo source tip coverage drift")
            if abs(sourceinfo["source_CWM_height_observed_on_original_direct_tips"]-
                   EXPECTED_SOURCE_DIRECT_TIP_HEIGHT[tree][study_ix])>1e-10:
                raise ValueError("Genuine original European rail study observed height CWM drift")
            clades=original_true_clade_groups(
                [lookup[n] for n in names],par,depth,
                .75*tree_depth["original_total_root_to_tip_distance_median"])
            source_baseline=source_null(
                P,T,clades,geo,info["ii"],info["jj"],study_ix,tree_ix,
                REFERENCE_PHYLO_MODE_INDEX,draws)
            if source_baseline["metrics"] is None:
                raise ValueError("Original v29.6 source 75%-tree-clade null not identifiable anymore")
            previous=source_baseline["metrics"]["height"][
                "original_observed_minus_original_clade_label_null"]
            if abs(previous-REF_GENUINE_V296[tree][study])>1e-10:
                raise ValueError(f"Original source v29.6 75%-clade effect changed: {tree}/{study}")
            occupancy=original_global_taxon_occupancy_quintiles(P)
            groups=joint_true_clade_and_occupancy_groups(clades,occupancy,len(names))
            support=original_source_joint_exchangeability(P,groups,np.asarray(geo,dtype=object))
            joint=None
            if support["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"]:
                joint=source_null(P,T,groups,geo,info["ii"],info["jj"],study_ix,tree_ix,
                                  NEW_JOINT_PHYLO_AND_OCCUPANCY_MODE_INDEX,draws)
                if joint["status"]!="SOURCE_REAL_PHYLOGENETIC_CLADE_DEPTH_ESTIMABLE" or joint["metrics"] is None:
                    raise ValueError("New joint taxon source null passed stated support yet original engine rejected")
            backbones[study]={
                "original_source_site_tree_tip_and_geography_QA":sourceinfo,
                "original_phylo75_v296_reproduced_exact":True,
                "original_v296_height_phylo75_excess":previous,
                "original_joint_phylo75_occupancy_source_support":support,
                "joint_phylo75_and_occupancy_observed_minus_null":
                    {key:joint["metrics"][key] for key in ("height","seed","joint")}
                    if joint is not None else None,
                "original_physical_author_source_sites_and_species_UNCHANGED":True,
                "original_75pct_source_clade_trait_vectors_AND_occupancy_group_preserved":True,
            }
        full[tree]={
            "original_VPhyloMaker2_tree_source_QA":sourceqa,
            "original_tree_root_to_tip_branch_length_QA":tree_depth,
            "original_rail_publications":backbones}
        supported={
          study:rec["joint_phylo75_and_occupancy_observed_minus_null"]["height"][
              "original_observed_minus_original_clade_label_null"]
          for study,rec in backbones.items()
          if rec["joint_phylo75_and_occupancy_observed_minus_null"] is not None}
        loo=[mean(v for k,v in supported.items() if k!=leave)
             for leave in supported] if len(supported)>=3 else []
        bytree[tree]={
            "original_true_rail_publications_with_joint_source_support":len(supported),
            "original_unidentifiable_rail_publications":[k for k in STUDIES if k not in supported],
            "source_study_equal_height_excess_if_at_least_three":
                mean(supported.values()) if len(supported)>=3 else None,
            "source_studies_positive_if_at_least_three":
                sum(v>0 for v in supported.values()) if len(supported)>=3 else None,
            "source_study_leave_one_out_mean_range_NOT_global95CI":
                [min(loo),max(loo)] if loo else None,
            "original_study_individual_height_signed_effects":supported,
            "source_4_publications_are_not_171_independent_studies":True,
        }
    result={
        "version":"v29.8",
        "status":"TRUE_ORIGINAL_4RAIL_AUTHOR_SOURCE_REAL_PHYLOGENY_AND_OCCUPANCY_CONDITIONED",
        "original_first_party_rail_publications_n":len(STUDIES),
        "original_source_rail_physical_sites_n":sum(SITES),
        "original_plant_source_original_publisher_studies_are_not_alternative_tree_replicates":True,
        "original_mega_tree_repository_commit":COMMIT,
        "original_phylogenetic_cut_fraction_of_true_root_to_tip_median":.75,
        "original_source_species_occupancy_frequency_bins":5,
        "source_original_site_species_incidence_and_locality_FIXED":True,
        "original_taxon_paired_2D_height_seed_vectors_conserved_WITHIN_phylo75_OCCUPANCY_INTERSECTION":True,
        "original_v296_source_398_draw_label_null_reproduced_for_all_12_pairs":True,
        "source_genus_family_phylogenetic_ancestry_not_a_causal_rail_selection_estimate":True,
        "source_true_original_tree_and_study_outputs":full,
        "source_original_study_equal_descriptive_by_tree":bytree,
        "original_source_plant_taxon_or_gps_rows_NOT_exported":True,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,default=source_json_native)+"\n",encoding="utf8")
    print(json.dumps({"version":"v29.8","status":result["status"],
       "phylogeny_occupancy_joint_source_support_and_height":bytree},
       indent=2,ensure_ascii=False,default=source_json_native))
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-original-source",type=Path,required=True)
    p.add_argument("--tree-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS_PER_CHAIN)
    a=p.parse_args()
    run(a.private_original_source,a.tree_dir,a.output,a.draws)
if __name__=="__main__":main()

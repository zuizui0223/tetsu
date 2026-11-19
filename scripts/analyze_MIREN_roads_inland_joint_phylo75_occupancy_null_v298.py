#!/usr/bin/env python3
"""v29.8 MIREN original mountain ROADS+paired inland: phylo75 x occurrence bin null.

True original MIREN site/year/ROAD Plot1 and same-site INLAND Plot3 botany
is FROZEN; exact original published tree-tip taxa and Oldenburg LE DA original
paired height+seed vectors only. ORIGINAL v29.7 effects MUST reproduce
source-authenticated results before calculating the new JOINT null.
One global road field program != additional railway study or 8 papers.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from statistics import mean
import numpy as np
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import read_source_miren,regional_comparison
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 import PANEL_DEFINITIONS,PINNED_PRIOR_ROAD_RESULT
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from audit_VPhyloMaker2_RData_python_no_R_environment_v296 import inspect_tree_rdata,TREES,COMMIT
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import original_root_depths,original_true_clade_groups,source_json_native
from validate_MIREN_original_roads_true_phylogeny_clade_negative_control_v297 import (
  trait_pool_original_matched_regions,construct_original_region_tree_panel,
  run_region_null,MODES as PRIOR_MODES)
from phylogeny_occupancy_joint_label_null_v298 import (
  original_global_taxon_occupancy_quintiles,
  joint_true_clade_and_occupancy_groups,original_source_joint_exchangeability)
DRAWS=199
JOINT_MODE="real_phylo_clade75_and_original_source_occupancy_quintile_joint"
ZONES=("roadside","interior","roadside_minus_interior")
METRICS=("height","seed","joint")
PANEL_TYPES=("primary_5species_25percent","strict_10species_50percent")

def source_study_incidence(panel):
    """Recover true observed 0/1 source vegetation separately for both zones."""
    both=[];geos=[];habitat={"roadside":[],"interior":[]};roads={}
    for i,road in enumerate(panel["roads"]):
        a=(np.asarray(road["roadside"])>0).astype(np.uint8)
        b=(np.asarray(road["interior"])>0).astype(np.uint8)
        if a.shape!=b.shape or len(a)<4:raise ValueError("Original MIREN 1:1 roadside/inland source sites changed")
        label="road_ordinal_"+str(i+1)  # anonymize actual source road ID
        both.extend((a,b))
        geos.extend([label]*(2*len(a)))
        roads[label]=(a,b)
        habitat["roadside"].append(a)
        habitat["interior"].append(b)
    P=np.vstack(both)
    actual=np.asarray(geos,dtype=object)
    habitat={k:np.vstack(v) for k,v in habitat.items()}
    if P.shape[1]!=len(panel["species"]) or len(P)!=len(actual):
        raise ValueError("Original source MIREN site×taxon field matrix broken")
    if not np.array_equal(P.sum(axis=0),panel["original_occupancy"]):
        raise ValueError("Original roadside plus inland taxon occupancy did not remain fixed")
    return P,actual,habitat,roads

def run(miren:Path,height:Path,seed:Path,tree_dir:Path,reference:Path,out:Path,draws:int=DRAWS):
    if draws!=DRAWS:raise ValueError("Must use same two-seed 398 permutations as frozen original v29.7")
    for filename,key in ((height,"height"),(seed,"seed")):
        if hashlib.sha256(filename.read_bytes()).hexdigest()!=LEDASHA[key]:
            raise ValueError("Oldenburg original height/seed source SHA changed")
    raw,source_qa=read_source_miren(miren)
    h,_=read_leda(height,"height_m");s,_=read_leda(seed,"seed_mass_mg")
    original_lead,traits,trait_qa=trait_pool_original_matched_regions(raw,h,s)
    if source_qa["original_source_record_QA"]["original_source_species_observation_rows"]!=167804:
        raise ValueError("Original MIREN field observation row count drifted")
    old=json.loads(reference.read_text())
    if old.get("version")!="v29.7":raise ValueError("Original v29.7 source-backed aggregate must be provided")
    panels={}
    for label,(minimum,share,year,nreg,nroads,nsites) in PANEL_DEFINITIONS.items():
        regions,qa=valid_region_panels(raw,original_lead,minimum,share,year)
        counts=(len(regions),sum(len(r) for r in regions.values()),
                sum(len(sites) for roads in regions.values() for _,sites in roads))
        if counts!=(nreg,nroads,nsites):
            raise ValueError("Original true MIREN source 8/22/272 or 5/13/192 field cohorts changed")
        if label.startswith("primary"):
            m=mean(regional_comparison(roads,original_lead,reg)[
                "source_2D_CWM_beta_road_minus_interior"] for reg,roads in regions.items())
            if abs(m-PINNED_PRIOR_ROAD_RESULT)>1e-9:
                raise ValueError("Source-authenticated v24.9 original roadside−inland 2D beta drift")
        panels[label]=(regions,qa,minimum)
    alltrees={};descriptions={}
    for tree_ix,tree in enumerate(TREES):
        obj,treeqa=inspect_tree_rdata(tree_dir/("GBOTB.extended."+tree+".rda"),tree)
        tips=[str(x).replace("_"," ") for x in np.asarray(obj["tip.label"]).reshape(-1)]
        if len(set(tips))!=len(tips):raise ValueError("Real source tree normalized plant tips not unique")
        lookup={v:i+1 for i,v in enumerate(tips)}
        parent,depth,age=original_root_depths(
            np.asarray(obj["edge"]),np.asarray(obj["edge.length"]),
            len(tips),int(np.asarray(obj["Nnode"]).reshape(-1)[0]))
        alltrees[tree]={"source_tree":treeqa,"source_tree_depth":age,"original_panels":{}}
        descriptions[tree]={}
        for panel_name,(regions,qa,minimum) in panels.items():
            original_supported={};new_supported={};perregion={}
            for reg_ix,(reg,roads) in enumerate(sorted(regions.items())):
                panel,support=construct_original_region_tree_panel(roads,traits,lookup,minimum)
                if panel is None:
                    perregion[reg]={"status":"NOT_IDENTIFIABLE_ORIGINAL_PHYLOGENY_DIRECT_TIP_SUPPORT",
                                    "original_phylo_source_cohort_QA":support,
                                    "clade_and_occupancy_joint_source_QA":None,"new_joint_null":None}
                    continue
                original_clades=original_true_clade_groups(
                    [lookup[t] for t in panel["species"]],parent,depth,
                    .75*age["original_total_root_to_tip_distance_median"])
                source_baseline=run_region_null(panel,original_clades,PRIOR_MODES[3],tree_ix,reg_ix,draws)
                original_effects=source_baseline.get("source_original_region_two_zone_effects")
                if original_effects is None:raise ValueError("Previously source verified clade75 region became unidentifiable")
                original_supported[reg]=original_effects
                P,geo,habitat,road_study=source_study_incidence(panel)
                bins=original_global_taxon_occupancy_quintiles(P)
                joint=joint_true_clade_and_occupancy_groups(original_clades,bins,len(panel["species"]))
                eligibility=original_source_joint_exchangeability(P,joint,geo,habitat,road_study)
                new=None
                if eligibility["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"]:
                    test=run_region_null(panel,joint,JOINT_MODE,tree_ix,reg_ix,draws)
                    new=test.get("source_original_region_two_zone_effects")
                    if new is None:
                        raise ValueError("Original joint source support passed yet source numerical null returned unidentifiable")
                    new_supported[reg]=new
                perregion[reg]={"status":"ORIGINAL_MIREN_ROAD_INLAND_PAIRED_PHYLOTIP_SUPPORTED",
                    "original_phylo_source_cohort_QA":support,
                    "original_v297_75pct_clade_null_reproduced":original_effects,
                    "clade_and_occupancy_joint_source_QA":eligibility,
                    "new_joint_null":new}
            ref=old["primary_75pct" if panel_name.startswith("primary") else "strict_75pct"][tree]
            if sorted(original_supported)!=sorted(ref["region_codes_supported"]):
                raise ValueError("Original v29.7 source region and taxonomy support changed")
            site_n=sum(perregion[r]["original_phylo_source_cohort_QA"][
                "retained_original_source_site_pairs"] for r in original_supported)
            if site_n!=ref["original_paired_sites_remaining"]:
                raise ValueError("Original MIREN true paired roadside/inland site retention changed")
            for zone,key in (("roadside","roadside"),("interior","interior"),
                              ("roadside_minus_interior","roadside_minus_interior")):
                value=mean(d[zone]["height"]["observed_minus_trait_label_source_null"]
                           for d in original_supported.values())
                if abs(value-ref[key])>1e-10:
                    raise ValueError("Earlier original MIREN source exact 75%-clade effect failed reproduction")
            stats={}
            for z in ZONES:
                stats[z]={}
                for metric in METRICS:
                    v={reg:d[z][metric]["observed_minus_trait_label_source_null"]
                       for reg,d in new_supported.items()}
                    if len(v)>=3:
                        loo=[mean(val for k,val in v.items() if k!=leave) for leave in v]
                        stats[z][metric]={
                            "original_region_code_evaluable":len(v),
                            "original_region_equal_mean_excess":mean(v.values()),
                            "n_original_region_codes_positive":sum(a>0 for a in v.values()),
                            "source_region_leave_one_out_range_NOT_CI":[min(loo),max(loo)],
                            "by_source_original_region":v}
                    else:stats[z][metric]=None
            alltrees[tree]["original_panels"][panel_name]={
                "source_original_road_year_site_support":qa,
                "n_original_v297_phylo75_supported_source_regions":len(original_supported),
                "n_new_joint_phylo_occupancy_source_regions":len(new_supported),
                "source_original_region_records":perregion}
            descriptions[tree][panel_name]=stats
    ans={"version":"v29.8",
       "status":"ORIGINAL_NONRAIL_ROADS_AND_INLAND_JOINT_PHYLO75_OCCUPANCY",
       "original_MIREN_publisher_source_QA":source_qa,
       "original_LEDA_2D_source_taxa_QA":trait_qa,
       "original_v297_source_baseline_REPRODUCED_in_all_six_tree_panels":True,
       "original_publisher_railway_field_studies_n":0,
       "original_independent_mountain_road_field_programs_n":1,
       "source_original_plant_studies_not_tree_alternative_replicates":True,
       "source_original_first_eligible_year_same_road_site_roadside_inland_frozen":True,
       "original_source_phylo_genuine_tree_commit":COMMIT,
       "source_tree_cut_root_to_tip_fraction":.75,
       "source_region_taxon_occupancy_bin_count":5,
       "original_plant_taxon_2D_height_seed_vectors_permitted_SWAP_only_in_true_clade_AND_occupancy_BIN":True,
       "source_two_technical_label_chains":2,"technical_draws_per_chain":draws,
       "source_true_trees_and_regions":alltrees,
       "source_original_region_equal_joint_null_summaries_if_at_least_three":descriptions,
       "no_raw_original_site_taxon_or_road_ids_and_phylogeny_tip_leaf_names_exported":True,
       "no_rail_specific_adaptation_or_seed_transport_identified":True}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(ans,ensure_ascii=False,indent=2,default=source_json_native)+"\n")
    print(json.dumps({"version":"v29.8","source":"ONE genuine MIREN ROAD field survey",
       "tree_panel_road_inland_height":{tree:{panel:{
         "n_joint_source_regions":alltrees[tree]["original_panels"][panel][
             "n_new_joint_phylo_occupancy_source_regions"],
         "mean_height_excess":{
           z:(descriptions[tree][panel][z]["height"] or {}).get("original_region_equal_mean_excess")
           for z in ZONES}}
         for panel in PANEL_TYPES} for tree in TREES}},ensure_ascii=False,indent=2,default=source_json_native))
    return ans

def main():
    p=argparse.ArgumentParser()
    for name in ("miren","height","seed","tree-dir","original-v297-summary","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS)
    a=p.parse_args()
    run(a.miren,a.height,a.seed,a.tree_dir,a.original_v297_summary,a.output,a.draws)
if __name__=="__main__":main()

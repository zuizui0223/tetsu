#!/usr/bin/env python3
"""v29.6 fail-closed original 12 study-by-backbone source support and true clade null.

The original 3 tree nomenclatures are ALTERNATIVE backbones, never 12
independent botanical field studies. Source R files are original study level
only; no species/plot records or exact source GPS are exported.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

STUDIES=("CZ_Czechia","IT_Italy_Parma","HR_BA_2017_one_study","SI_Slovenia_100m_track")
BACKBONES=("LCVP","WP","TPL")
MODES=("original_species_occupancy_quintile",
       "true_tree_clade_depth_fraction_0.25",
       "true_tree_clade_depth_fraction_0.50",
       "true_tree_clade_depth_fraction_0.75")
SOURCE_COUNTS={
 "LCVP":(433,245,163,195),
 "WP":(416,242,161,189),
 "TPL":(444,248,170,202)
}
EXPECTED_ORIGINAL_PHYSICAL_STUDY_SITES=(39,28,15,89)
N_PERM=398
PHYS_ORIGINAL_GEO_MIN=2

def to_f(v,*,optional=False)->float|None:
    if v in ("",None):
        if optional:return None
        raise ValueError("Expected finite source empirical output missing")
    f=float(v)
    if not math.isfinite(f):raise ValueError("Nonfinite observed or source null value")
    return f

def validate_source_tip_qa(v295:dict)->dict:
    if v295.get("version")!="v29.5":raise ValueError("Real original v29.5 source QA missing")
    rows=v295.get("study_backbone_source_coverage")
    if not isinstance(rows,list) or len(rows)!=12:
        raise ValueError("12 original source paper tree matches absent")
    indexed={}
    for item in rows:
        key=(item["phylogeny"],item["original_study"])
        if key in indexed:raise ValueError("Repeated true author study×backbone source ID")
        indexed[key]=item
    for tree in BACKBONES:
        for i,study in enumerate(STUDIES):
            q=indexed.get((tree,study))
            if not q:raise ValueError("Missing authoritative v29.5 original field sample")
            if q["prespecified_exact_tip_no_grafting_inference_gate_status"]!=\
              "DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_PASS":
                raise ValueError("An original field cohort cannot support genuine branches")
            if int(q["exact_original_backbone_tip_taxon_n"])!=SOURCE_COUNTS[tree][i]:
                raise ValueError("True original source exact-tip species number changed")
            if int(q["source_original_station_or_ballast_site_n"])!=\
               EXPECTED_ORIGINAL_PHYSICAL_STUDY_SITES[i]:
                raise ValueError("Original field-study physical site grain changed")
            site_fraction=to_f(q["original_field_site_tip_match_fraction_median"])
            if not (0<=site_fraction<=1):
                raise ValueError("v29.5 source site-tip fraction invalid: must be actual rows/sites")
            if to_f(q["original_field_site_ge5_direct_tip_fraction"])<.8:
                raise ValueError("Too few original field sites with five real unaugmented tree tips")
    return indexed

def aggregate(source_rows:list[dict],v295:dict,tree_properties:str)->dict:
    v295_index=validate_source_tip_qa(v295)
    expected={(tree,study,mode) for tree in BACKBONES for study in STUDIES for mode in MODES}
    actual={(r["backbone"],r["study"],r["null_model"]) for r in source_rows}
    if len(source_rows)!=48 or actual!=expected:
        raise ValueError("Incomplete or duplicate actual tree×study×source model rows")
    unique={}
    for r in source_rows:
        key=(r["backbone"],r["study"],r["null_model"])
        if key in unique:raise ValueError("Duplicate source effect within a genuine field study")
        unique[key]=r
    # Check R runner original-source file was always valid and source
    # real tip-restriction independently reproduces v29.5 author data.
    for tree in BACKBONES:
        for i,study in enumerate(STUDIES):
            cohort=[unique[(tree,study,mode)] for mode in MODES]
            qa=v295_index[tree,study]
            h_expected=to_f(qa["same_study_retained_sites_direct_tip_taxa_CWM_height_between_site_distance"])
            v=[to_f(r["original_observed_exact_backbone_CWM_height"]) for r in cohort]
            if not all(abs(x-h_expected)<1e-9 for x in v):
                raise ValueError("True original height CWM beta mismatched original v29.5 authentic source")
            for r in cohort:
                if int(r["n_original_physical_field_sites"])!=EXPECTED_ORIGINAL_PHYSICAL_STUDY_SITES[i]:
                    raise ValueError("Source same physical field site population changed")
                if int(r["n_original_trait_known_exact_tree_tips"])!=SOURCE_COUNTS[tree][i]:
                    raise ValueError("True tree original matched-taxonomy source count changed")
                if int(r["n_orig_source_study_plant_geographical_groups"])<PHYS_ORIGINAL_GEO_MIN:
                    raise ValueError("Insufficient original geographic field group support")
                if r["original_site_species_geography_FIXED"]!="TRUE" or\
                   r["2D_LEDA_height_and_seed_mass_paired_vector_preserved"]!="TRUE":
                    raise ValueError("Original source trait/site invariants not preserved")
                status=r["prespecified_support_status"]
                n=int(r["n_technical_trait_label_shuffles_if_identifiable"])
                if status in ("IDENTIFIABLE_DIRECT_BACKBONE_CLADE_LABEL_NULL",
                             "ORIGINAL_TIP_SUBSET_OCCUPANCY_NULL_REFERENCE"):
                    if n!=N_PERM:raise ValueError("Original fixed two-chain null technical assignments changed")
                    obs=to_f(r["original_observed_exact_backbone_CWM_height"])
                    null=to_f(r["source_height_null_mean"])
                    delta=to_f(r["observed_minus_clade_null_height"])
                    if abs(obs-null-delta)>1e-10:
                        raise ValueError("Source distance null subtraction does not balance")
                    if not (0<=to_f(r["source_height_technical_upper_rank_NOT_field_p"])<=1):
                        raise ValueError("Technical rank outside [0,1]")
                else:
                    if status!="NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY":
                        raise ValueError("Unknown original phylogenetic support status")
                    if n!=0 or to_f(r["observed_minus_clade_null_height"],optional=True) is not None:
                        raise ValueError("An unidentified original tree clade was assigned a false biological effect")
    model_aggregates={}
    for tree in BACKBONES:
        model_aggregates[tree]={}
        for mode in MODES:
            samples=[unique[(tree,study,mode)] for study in STUDIES]
            done={r["study"]:to_f(r["observed_minus_clade_null_height"],optional=True)
                  for r in samples if r["prespecified_support_status"] in
                  ("IDENTIFIABLE_DIRECT_BACKBONE_CLADE_LABEL_NULL",
                   "ORIGINAL_TIP_SUBSET_OCCUPANCY_NULL_REFERENCE")}
            n=len(done)
            if n>=3:
                loo=[mean(v for k,v in done.items() if k!=omitted) for omitted in done]
                values={"original_source_pub_equal_height_excess":mean(done.values()),
                        "original_source_publications_positive":sum(v>0 for v in done.values()),
                        "publication_leave_one_out_range":[min(loo),max(loo)],
                        "all_original_independent_source_paper_signed_results":done}
            else:
                values=None
            model_aggregates[tree][mode]={
                "original_publications_supported":n,
                "source_studies_not_estimable":[x for x in STUDIES if x not in done],
                "descriptive_height":values,
                "not_a_population_prediction_interval":True,
            }
    return {
      "version":"v29.6",
      "status":"GENUINE_ORIGINAL_RAILWAY_PHYLOGENY_BRANCH_CLADE_NULL",
      "original_railway_true_field_publications":len(STUDIES),
      "true_sites_not_repeated_meta_studies":171,
      "original_tree_versions_three_nomenclatures_not_three_independent_publications":list(BACKBONES),
      "original_source_original_VPhyloMaker2_git_commit":"7af3fb5152f691af2e4ec9d5e2e467d1b50505e9",
      "original_tree_source_barcodes":{
        "LCVP":"0d761cb0bf36392cc13759a7ac9026db564217eb",
        "WP":"8d3777762cb32a80e8c6ea097a0456e2b63a94c6",
        "TPL":"d7bc9b69bf3cb11c3c3d090282723b6e73b1de55",
      },
      "original_tree_source_qa_dcf":tree_properties,
      "source_null_definitions":"paired zlog height+seed species trait vectors, true source original site plant membership and geographic grouping fixed, phylogeny clades at fractions of source root-to-tip median branch depth; actual tips ONLY, no graft",
      "predeclared_clade_branch_depth_fractions":[.25,.50,.75],
      "source_original_direct_backbone_phylogeny_all_12_study_tree_pairs_pass":True,
      "two_fixed_technical_trait_null_chains_draws_per_chain":N_PERM//2,
      "by_original_backbone_and_clade_null":model_aggregates,
      "by_original_publication_source_summary_only":source_rows,
      "no_molecular_DNA_sequence_evidence_asserted_for_every_original_tip":True,
      "no_original_source_species_list_raw_trait_site_matrix_or_location_exported":True,
      "no_direct_rail_ballast_vs_nonrail_causal_treatment_identified":True,
      "no_global_field_p_value_or_mega_tree_phylogenetic_variance_partition":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-r",type=Path,required=True)
    p.add_argument("--previous-v295",type=Path,required=True)
    p.add_argument("--tree-properties",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    with a.source_r.open(encoding="utf8",newline="") as f:
        rows=list(csv.DictReader(f))
    previous=json.loads(a.previous_v295.read_text(encoding="utf8"))
    out=aggregate(rows,previous,a.tree_properties.read_text(encoding="utf8"))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
       "version":out["version"],
       "source_null_summary_by_tree":{t:{
           m:{"n_supported_papers":r["original_publications_supported"],
              "effect":r["descriptive_height"]["original_source_pub_equal_height_excess"]
                       if r["descriptive_height"] else None,
              "study_positive":r["descriptive_height"]["original_source_publications_positive"]
                       if r["descriptive_height"] else None}
           for m,r in parts.items()} for t,parts in out["by_original_backbone_and_clade_null"].items()}
    },indent=2))
if __name__=="__main__":
    main()

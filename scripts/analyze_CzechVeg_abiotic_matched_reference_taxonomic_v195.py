#!/usr/bin/env python3
"""v19.5 source-original abiotic matched regional taxonomic railway test.

Original CzechVeg v1.1 117739 source plot headers + 2510340 vascular plot
species records, full MD5 gate before fitting. Source RAIL exposure is TEXT,
not GPS-confirmed trackbed. Reference plots are original non-transport-text
sites from *other publisher studies* in same region, decade and area.
Two contrasting choices on exactly same original RAIL source plots:
  (a) all eligible other-publication regional source plots
  (b) among (a), measured ALTITUDE within ±100 m, or altitude ±100m
      AND slope within ±10 degrees.
No syntaxonomic classification (TH/TF) or plant response is used to choose
the physical reference. 500 resampling draws from observed source metadata
are TECHNICAL uncertainty, not independent studies or p values.
"""
from __future__ import annotations
import argparse,json,random,statistics,sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    DOI,get_meta,matched_rail_source
)
from audit_CzechVeg_abiotic_regional_reference_opportunity_v194 import (
    abiotic_rows,source_reference_options,matching
)
from analyze_CzechVeg_original_independent_reference_taxonomic_v191 import (
    stable_sample,source_species,genus_set,cross_similarity,stats
)

SEED=20261116
DRAW_COUNT=500
LOO_DRAW_COUNT=80
MIN_CONTROLS=5
CAP_PER_SOURCE_RAIL=250
FRAME_TYPES=("elev100","elev100_slope10")
POOL_NAMES=("rail","abiotic","all")

def source_arms(metadata,blocks,measured):
    orig=source_reference_options(metadata,blocks)
    out={key:{} for key in FRAME_TYPES+("all",)}
    for rid,refs in orig.items():
        if len(refs)<MIN_CONTROLS:continue
        all_candidates=stable_sample(refs,CAP_PER_SOURCE_RAIL,
                                     rid+"|v195|all")
        out["all"][rid]=all_candidates
        for kind in FRAME_TYPES:
            match_ids=[pid for pid in refs if kind in matching(measured[rid],measured[pid])]
            if len(match_ids)>=MIN_CONTROLS:
                out[kind][rid]=stable_sample(match_ids,CAP_PER_SOURCE_RAIL,
                                            rid+"|v195|"+kind)
    needed=set()
    for panel in out.values():
        needed.update(panel)
        for candidates in panel.values():needed.update(candidates)
    return out,needed,{
       "source_publisher_rail_plots_with_min5_prematched_background":len(out["all"]),
       "source_publisher_rail_plots_with_min5_elev100":len(out["elev100"]),
       "source_publisher_rail_plots_with_min5_elev100_slope10":len(out["elev100_slope10"]),
       "source_plot_ids_chosen_before_plant_species_source_rows":True,
       "only_original_field_elevation_slope_selected":True,
       "source_original_environmental_ref_max250_per_rail":True,
       "original_matched_source_plot_ids_to_read":len(needed)
    }

def frame_for_type(arms,metadata,plants,match_type,richness=False):
    selected=defaultdict(list)
    refs_abiotic=arms[match_type]
    refs_all=arms["all"]
    for rid in sorted(set(refs_abiotic)&set(refs_all)):
        if rid not in plants or not plants[rid]:
            continue
        a=[p for p in refs_abiotic[rid] if plants.get(p)]
        b=[p for p in refs_all[rid] if plants.get(p)]
        if richness:
            richness_R=len(plants[rid])
            a=[p for p in a if 2/3 <=len(plants[p])/richness_R<=1.5]
            b=[p for p in b if 2/3 <=len(plants[p])/richness_R<=1.5]
        if len(a)<MIN_CONTROLS or len(b)<MIN_CONTROLS:
            continue
        selected[metadata[rid]["region"]].append((rid,{"abiotic":a,"all":b}))
    return dict(selected),{
        "selected_original_rail_plots":sum(len(v) for v in selected.values()),
        "distinct_source_phyto_regions":len(selected),
        "same_source_rail_sites_for_abiotic_and_all_background":True,
        "original_source_rail_sites_per_region_no_plot_IDs":{
            str(k):len(v) for k,v in sorted(selected.items())},
        "control_plot_richness_2thirds_to_1point5x_original_rail_richness":richness,
        "strict_original_study_year_region_area_then_abiotic_matching":True
    }

def sample_once(frame,plants,regions,rng):
    plot_sets={k:[] for k in POOL_NAMES}
    for reg in regions:
        rid,refs=rng.choice(frame[reg])
        ids={"rail":rid,"abiotic":rng.choice(refs["abiotic"]),
             "all":rng.choice(refs["all"])}
        for typ,pid in ids.items():plot_sets[typ].append(plants[pid])
    result={}
    for rank in ("species","genus"):
        sets={k:(v if rank=="species" else [genus_set(s) for s in v])
              for k,v in plot_sets.items()}
        for typ,val in sets.items():
            result[f"{rank}_{typ}"]=cross_similarity(val)
        result[f"{rank}_rail_minus_abiotic"]=result[f"{rank}_rail"]-result[f"{rank}_abiotic"]
        result[f"{rank}_rail_minus_all"]=result[f"{rank}_rail"]-result[f"{rank}_all"]
        result[f"{rank}_abiotic_minus_all"]=result[f"{rank}_abiotic"]-result[f"{rank}_all"]
    for typ,arr in plot_sets.items():
        result[f"mean_source_taxa_{typ}"]=statistics.mean(len(s) for s in arr)
    return result

FOCUS=("species_rail_minus_abiotic","species_rail_minus_all",
       "species_abiotic_minus_all","genus_rail_minus_abiotic",
       "genus_rail_minus_all","genus_abiotic_minus_all")

def evaluate(frame,plants,seed=SEED,draws=DRAW_COUNT):
    regs=sorted(frame)
    source_n=sum(len(v) for v in frame.values())
    quality={
        "original_source_phyto_region_count":len(regs),
        "original_rail_plot_n_on_common_abiotic_vs_all_reference":source_n,
        "one_study_original_source_CzechVeg_database":True,
        "same_original_rail_plot_reused_for_two_reference_types_per_draw":True,
        "source_region_pairs_nonindependent_per_draw":len(regs)*(len(regs)-1)//2,
        "reference_source_studies_distinct_from_rail_original_study":True
    }
    if len(regs)<3:
        return {**quality,"identified":False,"not_estimable_reason":"<3 region groups after exact source matching and species coverage"}
    rng=random.Random(seed)
    samples=[sample_once(frame,plants,regs,rng) for _ in range(draws)]
    result={k:stats([s[k] for s in samples]) for k in samples[0]}
    loo={}
    for i,region in enumerate(regs):
        rest=[x for x in regs if x!=region]
        if len(rest)<3:continue
        rr=random.Random(seed+10000+i)
        ss=[sample_once(frame,plants,rest,rr) for _ in range(LOO_DRAW_COUNT)]
        loo[str(region)]={k:round(statistics.mean(s[k] for s in ss),7) for k in FOCUS}
    return {**quality,"identified":True,
            "source_technical_draws_not_independent_field_n":draws,
            "summary":result,"region_leave_one_out":loo,
            "region_leave_one_out_sign_counts":{
                k:{"negative":sum(s[k]<0 for s in loo.values()),
                   "positive":sum(s[k]>0 for s in loo.values()),
                   "n":len(loo)} for k in FOCUS},
            "resampling_intervals_not_field_confidence_intervals":True,
            "does_not_measure_real_railway_transport_or_environment_causality":True}

def run(header_csv,species_csv):
    source_meta,original_header=get_meta(header_csv)
    blocks,rail_frame=matched_rail_source(source_meta)
    abiotic,abiotic_quality=abiotic_rows(header_csv,source_meta)
    arms,needed,armqa=source_arms(source_meta,blocks,abiotic)
    plant,species_qa=source_species(species_csv,needed)
    panels={}
    for kind in FRAME_TYPES:
        panels[kind]={}
        for richness in (False,True):
            key="richness_matched" if richness else "metadata_only"
            frame,quality=frame_for_type(arms,source_meta,plant,kind,richness)
            panels[kind][key]={"source_frame":quality,
                               "taxonomic_result":evaluate(frame,plant,seed=SEED+129*int(richness)+FRAME_TYPES.index(kind)*997)}
    return {
        "version":"v19.5",
        "original_CzechVeg_source_DOI":DOI,
        "source_original_headers_md5":original_header["publisher_original_headers_md5"],
        "source_original_species_md5":species_qa["publisher_species_csv_md5"],
        "source_original_2510340_rows_scanned":species_qa["all_original_species_rows"],
        "original_matched_rail_road_context_frame":rail_frame,
        "original_field_abiotic_availability":abiotic_quality,
        "preselected_abiotic_matched_control_opportunity":armqa,
        "original_selected_vascular_source_row_quality":species_qa,
        "taxonomic_results_from_abiotic_not_species_derived_syntaxon_controls":panels,
        "taxon_and_genus_source_phyto_region_contrasts_only":True,
        "no_rail_habitat_legal_location_verified":True,
        "control_not_verified_true_natural_land":True,
        "no_independent_study_p_value":True,
        "no_multi_country_global_generalization":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-headers",required=True)
    p.add_argument("--original-species",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    out=run(args.original_headers,args.original_species)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    quick={}
    for kind,cases in out["taxonomic_results_from_abiotic_not_species_derived_syntaxon_controls"].items():
        quick[kind]={}
        for scope,res in cases.items():
            val=res["taxonomic_result"]
            quick[kind][scope]={
                "frame":res["source_frame"],
                "identified":val["identified"],
                "key_responses":{k:val["summary"][k] for k in FOCUS if k in val.get("summary",{})},
                "loo_signs":val.get("region_leave_one_out_sign_counts")
            }
    print(json.dumps({"version":"v19.5","source_full_plant_rows":out["source_original_2510340_rows_scanned"],
                      "source_abiotic_reference":out["preselected_abiotic_matched_control_opportunity"],
                      "result":quick},ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()

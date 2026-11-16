#!/usr/bin/env python3
"""v19.6 SAME source rail sites for all vegetation-defined and abiotic controls.

Four independent original-source plot reference choices per same original
RAIL-text plot (first select RAIL before species contents inspected):
 (1) TH/TF/TG SYNTAXON-based vegetation class (species-composition derived)
 (2) nonanalog nontransport-text plots
 (3) ALTITUDE within ±100m original measured abiotic, irrespective SYNTAXON
 (4) all nontransport-text plots, all of a different original publication.
Same CzechVeg-Open source, same phytoregion, decade, plot-area band.
Only original species and genus Jaccard tested, not direct substrate identity
or railway causal intervention. Cross-region pairs are nonindependent;
Monte Carlo source resampling is NOT field CI or study p-value.
"""
from __future__ import annotations
import argparse,json,random,statistics,sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    DOI,get_meta,matched_rail_source
)
from audit_CzechVeg_abiotic_regional_reference_opportunity_v194 import abiotic_rows
from analyze_CzechVeg_same_rail_sites_analog_vs_background_v192 import candidate_arms
from analyze_CzechVeg_abiotic_matched_reference_taxonomic_v195 import source_arms
from analyze_CzechVeg_original_independent_reference_taxonomic_v191 import (
    source_species,genus_set,cross_similarity,stats
)

SEED=20261116
DRAWS=500
LOO_DRAWS=80
MIN_CONTROL=5
CONTROL_TYPES=("analog","nonanalog","abiotic","all")
TRACK_TYPES=("species","genus")
RESPONSE_TYPES=tuple(
    f"{t}_rail_minus_{c}"
    for t in TRACK_TYPES for c in CONTROL_TYPES
)+tuple(f"{t}_analog_minus_abiotic" for t in TRACK_TYPES)

def combined_arms(original_meta,railblocks,environment):
    code_arms,_=candidate_arms(original_meta,railblocks)
    abiotic_arms,_,abiotic_qa=source_arms(original_meta,railblocks,environment)
    sites=set(code_arms["analog"]) & set(code_arms["nonanalog"]) & set(code_arms["all"]) & set(abiotic_arms["elev100"])
    all_arms={}
    for rid in sorted(sites):
        per={"analog":code_arms["analog"][rid],
             "nonanalog":code_arms["nonanalog"][rid],
             "abiotic":abiotic_arms["elev100"][rid],
             "all":code_arms["all"][rid]}
        if any(len(v)<MIN_CONTROL for v in per.values()):
            continue
        all_arms[rid]=per
    read=set(all_arms)
    for group in all_arms.values():
        for refs in group.values():read.update(refs)
    return all_arms,read,{
        "source_eligible_rail_sites_pre_species":len(all_arms),
        "original_site_ids_to_read":len(read),
        "abiotic_previous_source_opportunity":abiotic_qa,
        "all_reference_types_formed_before_species_contents":True
    }

def common_frame(arms,plants,meta,richness=False):
    frame=defaultdict(list)
    overlaps=[]
    for rid in sorted(arms):
        if not plants.get(rid):continue
        groups={}
        for kind in CONTROL_TYPES:
            candidates=[x for x in arms[rid][kind] if plants.get(x)]
            if richness:
                n=len(plants[rid])
                candidates=[x for x in candidates if 2/3<=len(plants[x])/n<=1.5]
            groups[kind]=candidates
        if any(len(groups[k])<MIN_CONTROL for k in CONTROL_TYPES):
            continue
        overlaps.append(len(set(groups["analog"]) & set(groups["abiotic"])))
        frame[meta[rid]["region"]].append((rid,groups))
    return dict(frame),{
        "n_same_original_rail_text_sites":sum(map(len,frame.values())),
        "n_covered_source_phyto_regions":len(frame),
        "rail_original_sites_per_region_no_ID":{
            str(k):len(v) for k,v in sorted(frame.items())},
        "analog_vs_abiotic_original_control_site_overlap_total":sum(overlaps),
        "mean_overlap_same_region_rail_site":round(statistics.mean(overlaps),4) if overlaps else None,
        "same_rail_site_and_geographic_study_design_for_all_reference_choices":True,
        "matched_original_source_taxa_richness_within_1p5fold":richness
    }

def one_draw(frame,plant,regions,rng):
    lists={k:[] for k in ("rail",)+CONTROL_TYPES}
    for reg in regions:
        rail,controls=rng.choice(frame[reg])
        site_ids={"rail":rail}
        site_ids.update({kind:rng.choice(controls[kind]) for kind in CONTROL_TYPES})
        for typ,pid in site_ids.items():
            lists[typ].append(plant[pid])
    values={}
    for rank in TRACK_TYPES:
        transformed={
           k:([genus_set(s) for s in v] if rank=="genus" else v)
           for k,v in lists.items()}
        for k,v in transformed.items():
            values[f"{rank}_{k}"]=cross_similarity(v)
        for k in CONTROL_TYPES:
            values[f"{rank}_rail_minus_{k}"]=values[f"{rank}_rail"]-values[f"{rank}_{k}"]
        values[f"{rank}_analog_minus_abiotic"]=values[f"{rank}_analog"]-values[f"{rank}_abiotic"]
    for k,v in lists.items():
        values[f"mean_source_richness_{k}"]=statistics.mean(len(s) for s in v)
    return values

def analyze(frame,plant,seed=SEED,draws=DRAWS):
    regions=sorted(frame)
    info={"n_geographic_regions":len(regions),
          "n_source_rail_sites":sum(map(len,frame.values())),
          "nonindependent_source_region_pairs":len(regions)*(len(regions)-1)//2,
          "study_n_is_not_region_pair_count":True}
    if len(regions)<3:
        return {**info,"identified":False,"reason":"fewer than three phyto-regions with all reference cohorts; cannot compare"}
    rng=random.Random(seed)
    s=[one_draw(frame,plant,regions,rng) for _ in range(draws)]
    estimates={k:stats([r[k] for r in s]) for k in s[0]}
    leave={}
    for i,reg in enumerate(regions):
        smaller=[x for x in regions if x!=reg]
        if len(smaller)<3:continue
        rr=random.Random(seed+1007+i)
        sub=[one_draw(frame,plant,smaller,rr) for _ in range(LOO_DRAWS)]
        leave[str(reg)]={k:round(statistics.mean(x[k] for x in sub),7) for k in RESPONSE_TYPES}
    return {**info,"identified":True,
            "source_technical_resamples":draws,
            "summary":estimates,
            "leave_one_region_out":leave,
            "leave_region_sign_counts":{
                key:{"negative":sum(v[key]<0 for v in leave.values()),
                     "positive":sum(v[key]>0 for v in leave.values()),
                     "region_tests":len(leave)}
                for key in RESPONSE_TYPES},
            "syntaxonomic_reference_is_plant_composition_defined_not_independent_abiotic":True,
            "abiotic_reference_is_altitude_only_not_true_gravel_habitat":True,
            "no_field_study_p_values_or_CI":True}

def run(header,species):
    meta,head_qa=get_meta(header)
    original_blocks,matched=matched_rail_source(meta)
    environmental,env_qa=abiotic_rows(header,meta)
    arms,read,original_qa=combined_arms(meta,original_blocks,environmental)
    plants,species_qa=source_species(species,read)
    cases={}
    for richness in (False,True):
        name="species_richness_matched" if richness else "metadata_only"
        f,frameqa=common_frame(arms,plants,meta,richness)
        cases[name]={"source_common_frame":frameqa,
                     "source_comparative_result":analyze(f,plants,seed=SEED+2013*int(richness))}
    return {
        "version":"v19.6",
        "original_CzechVeg_DOI":DOI,
        "full_original_source_headers_md5":head_qa["publisher_original_headers_md5"],
        "full_original_source_species_md5":species_qa["publisher_species_csv_md5"],
        "source_original_species_rows":species_qa["all_original_species_rows"],
        "original_source_matched_study_rail_context":matched,
        "original_environment_header_qa":env_qa,
        "source_arm_pre_species_selection":original_qa,
        "source_species_qa":species_qa,
        "same_original_rail_plots_syntaxon_vs_elevation_vs_nonanalog_vs_all":cases,
        "one_country_CzechVeg_compilation":True,
        "absence_transport_text_does_not_prove_true_natural_reference":True,
        "unverified_source_rail_text_real_trackbed":True,
        "noncausal_exploratory_reference_definition_robustness":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-headers",required=True)
    p.add_argument("--original-species",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.original_headers,a.original_species)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    brief={}
    for case,val in out["same_original_rail_plots_syntaxon_vs_elevation_vs_nonanalog_vs_all"].items():
        r=val["source_comparative_result"]
        brief[case]={
            "frame":val["source_common_frame"],
            "identified":r["identified"],
            "means_and_technical_quantiles":{
                k:r["summary"][k] for k in RESPONSE_TYPES if k in r.get("summary",{})
            },
            "leave_region_signs":r.get("leave_region_sign_counts")
        }
    print(json.dumps({"version":out["version"],
                      "all_original_species_rows_verified":out["source_original_species_rows"],
                      "source_prematched":out["source_arm_pre_species_selection"],
                      "same_source_comparisons":brief},ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()

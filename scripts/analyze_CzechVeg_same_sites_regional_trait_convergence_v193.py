#!/usr/bin/env python3
"""v19.3: CzechVeg original source regional taxonomic versus functional convergence.

Same source rail-text plots, same regions, same selected reference plot IDs
for each taxonomic / Pladias-trait metric. Separate:
  lifeform (original public Pladias Raunkiaer categories)
  dispersal strategy (Sádlo nine archetypes, NOT observed train seed transport)
  traits jointly known (sensitivity for missingness; no imputation)

Source TH/TF/TG includes managed/secondary grassland. Original RAIL/ROAD
are TEXT labels, not verified railway trackbed. Nontransport-unlabeled
is not confirmed natural. A within-CzechVeg study-conditional descriptive
comparison, NOT railway specific cause, global transport or field p-value.
"""
from __future__ import annotations
import argparse,hashlib,json,random,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    DOI,get_meta,matched_rail_source
)
from analyze_CzechVeg_original_independent_reference_taxonomic_v191 import (
    source_species,genus_set,cross_similarity,stats
)
from analyze_CzechVeg_same_rail_sites_analog_vs_background_v192 import (
    candidate_arms,CATEGORIES
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import parse_traits

SEED=20261116
DRAW_COUNT=500
LOO_DRAW_COUNT=80
MIN_REFERENCES=5
MIN_KNOWN_TAXA=5
MIN_TRAIT_COVERAGE=0.50
TRAIT_KEYS=("lifeform","dispersal","both_known")
RESPONSE_KEYS=("taxonomic","genus","lifeform","dispersal")

def known_taxa(plant_set,trait,life,disp):
    if trait=="lifeform":
        return [t for t in plant_set if t in life]
    if trait=="dispersal":
        return [t for t in plant_set if t in disp]
    if trait=="both_known":
        return [t for t in plant_set if t in life and t in disp]
    raise ValueError("Unsupported trait coverage strategy")

def site_eligible(pid,plants,life,disp,trait):
    raw=plants.get(pid,set())
    if not raw:return False
    valid=known_taxa(raw,trait,life,disp)
    return len(valid)>=MIN_KNOWN_TAXA and len(valid)/len(raw)>=MIN_TRAIT_COVERAGE

def eligible_frame(arms,meta,plants,life,disp,trait,richness=False):
    good={pid:site_eligible(pid,plants,life,disp,trait) for pid in plants}
    choices={k:{} for k in CATEGORIES}
    for k in CATEGORIES:
        for rid,original_controls in arms[k].items():
            if not good.get(rid,False):
                continue
            eligible=[s for s in original_controls if good.get(s,False)]
            if richness:
                nr=len(plants[rid])
                eligible=[s for s in eligible
                          if 2/3 <=len(plants[s])/nr<=1.5]
            if len(eligible)>=MIN_REFERENCES:
                choices[k][rid]=eligible
    common=set.intersection(*(set(choices[k]) for k in CATEGORIES))
    frame=defaultdict(list)
    for rid in sorted(common):
        frame[meta[rid]["region"]].append((
            rid,{k:choices[k][rid] for k in CATEGORIES}
        ))
    return dict(frame),{
        "n_study_source_rail_sites_with_all_three_reference_types_and_trait_quality":len(common),
        "n_phyto_regions":len(frame),
        "rail_site_n_per_region":{str(reg):len(v) for reg,v in sorted(frame.items())},
        "n_total_sites_with_trait_eligibility_among_selected_original_source_sites":sum(good.values()),
        "min_original_observed_taxa_with_known_trait":MIN_KNOWN_TAXA,
        "min_original_trait_coverage_fraction":MIN_TRAIT_COVERAGE,
        "all_taxonomy_and_traits_evaluated_on_exact_same_selected_site_ids":True,
        "richness_matched_to_original_taxa_not_trait_subset":richness
    }

def profile(taxa,traits):
    obs=Counter(traits[t] for t in taxa if t in traits)
    n=sum(obs.values())
    if not n:return {}
    return {t:c/n for t,c in obs.items()}

def profile_overlap(a,b):
    return sum(min(a.get(t,0),b.get(t,0)) for t in set(a)|set(b))

def cross_profile_mean(profiles):
    vals=[profile_overlap(profiles[i],profiles[j])
          for i in range(len(profiles))
          for j in range(i+1,len(profiles))]
    if not vals:raise ValueError("Need >=2 source region profiles")
    return statistics.mean(vals)

def draw_once(regions,frame,plants,life,disp,rng):
    chosen={k:[] for k in ("rail",)+CATEGORIES}
    for region in regions:
        rail,reference_options=rng.choice(frame[region])
        sources={"rail":rail}
        sources.update({kind:rng.choice(reference_options[kind])
                        for kind in CATEGORIES})
        for label,pid in sources.items():
            chosen[label].append(plants[pid])
    metrics={}
    for group,taxon_sets in chosen.items():
        metrics[f"taxonomic_{group}"]=cross_similarity(taxon_sets)
        metrics[f"genus_{group}"]=cross_similarity(
            [genus_set(s) for s in taxon_sets])
        metrics[f"lifeform_{group}"]=cross_profile_mean(
            [profile(s,life) for s in taxon_sets])
        metrics[f"dispersal_{group}"]=cross_profile_mean(
            [profile(s,disp) for s in taxon_sets])
        metrics[f"original_taxa_richness_{group}"]=statistics.mean(
            len(s) for s in taxon_sets)
        metrics[f"lifeform_coverage_{group}"]=statistics.mean(
            sum(t in life for t in s)/len(s) for s in taxon_sets)
        metrics[f"dispersal_coverage_{group}"]=statistics.mean(
            sum(t in disp for t in s)/len(s) for s in taxon_sets)
    for response in RESPONSE_KEYS:
        for pool in CATEGORIES:
            metrics[f"{response}_rail_minus_{pool}"]=(
                metrics[f"{response}_rail"]-metrics[f"{response}_{pool}"])
        metrics[f"{response}_analog_minus_nonanalog"]=(
            metrics[f"{response}_analog"]-metrics[f"{response}_nonanalog"])
    return metrics

FOCUS_KEYS=tuple(
    f"{resp}_rail_minus_{pool}"
    for resp in RESPONSE_KEYS for pool in ("analog","nonanalog","all")
)+tuple(f"{resp}_analog_minus_nonanalog" for resp in RESPONSE_KEYS)

def result_for_frame(frame,plants,life,disp,seed=SEED,draws=DRAW_COUNT):
    regs=sorted(frame)
    info={
        "n_original_phyto_regions":len(regs),
        "n_source_rail_sites":len({rid for arr in frame.values() for rid,_ in arr}),
        "n_region_pairs_nonindependent":len(regs)*(len(regs)-1)//2,
        "source_one_rail_site_chosen_per_region_per_technical_draw":True,
        "source_same_selected_rail_and_reference_site_ids_for_all_response_metrics":True,
        "original_relevé_year_area_region_and_distinct_study_source_matched":True
    }
    if len(regs)<3:
        return {**info,"identified":False,"reason":"Fewer than three phyto-regions after trait and reference sampling gates"}
    rng=random.Random(seed)
    samples=[draw_once(regs,frame,plants,life,disp,rng) for _ in range(draws)]
    keys=list(samples[0])
    summaries={k:stats([x[k] for x in samples]) for k in keys}
    loo={}
    for i,reg in enumerate(regs):
        subset=[x for x in regs if x!=reg]
        if len(subset)<3:continue
        rr=random.Random(seed+9100+i)
        sample=[draw_once(subset,frame,plants,life,disp,rr)
                for _ in range(LOO_DRAW_COUNT)]
        loo[str(reg)]={k:round(statistics.mean(x[k] for x in sample),7)
                        for k in FOCUS_KEYS}
    return {**info,"identified":True,"n_technical_sampling_draws":draws,
            "summary_statistics":summaries,"region_leave_one_out_estimates":loo,
            "leave_one_region_out_sign_stability":{
                k:{"positive":sum(s[k]>0 for s in loo.values()),
                   "negative":sum(s[k]<0 for s in loo.values()),
                   "regions_tested":len(loo)} for k in FOCUS_KEYS},
            "intervals_are_technical_resampling_quantiles_not_field_CI":True,
            "independent_study_p_values_not_estimated":True}

def run(header_csv,species_csv,life_xlsx,dispersal_xlsx):
    meta,header_qa=get_meta(header_csv)
    blocks,source_frame=matched_rail_source(meta)
    arms,_=candidate_arms(meta,blocks)
    original_sites=set()
    for k in CATEGORIES:
        original_sites.update(arms[k])
        for options in arms[k].values():
            original_sites.update(options)
    plants,species_qa=source_species(species_csv,original_sites)
    disp,life,trait_qa=parse_traits(dispersal_xlsx,life_xlsx)
    # Original Pladias source has 3147 species-strategy and 4599 species-lifeform
    # rows. Only exact taxon-to-single unambiguous traits are used.
    if trait_qa["source_Pladias_original_diaspore_rows"]!=3147 or trait_qa[
        "source_Pladias_original_lifeform_rows"]!=4599:
        raise ValueError("Original Pladias trait source population changed")
    if len(set(life.values()))<4 or len(set(disp.values()))<7:
        raise ValueError("Original Pladias category coverage invalid")
    out={}
    for trait in TRAIT_KEYS:
        out[trait]={}
        for richness in (False,True):
            key="richness_matched" if richness else "metadata_only"
            frame,qa=eligible_frame(arms,meta,plants,life,disp,trait,richness)
            res=result_for_frame(frame,plants,life,disp,
                 seed=SEED+int(richness)*101+TRAIT_KEYS.index(trait)*203)
            out[trait][key]={"original_source_frame":qa,"responses":res}
    known={t for ss in plants.values() for t in ss}
    return {
        "version":"v19.3",
        "original_source_CzechVeg_DOI":DOI,
        "source_header_md5":header_qa["publisher_original_headers_md5"],
        "source_species_md5":species_qa["publisher_species_csv_md5"],
        "source_Pladias_DOI":"10.23855/preslia.2021.001",
        "original_species_row_count":species_qa["all_original_species_rows"],
        "original_trait_source_metadata":trait_qa,
        "n_source_original_distinct_vegetation_taxa_in_eligible_source_sites":len(known),
        "exact_original_taxon_match_lifeform":sum(t in life for t in known),
        "exact_original_taxon_match_diaspore":sum(t in disp for t in known),
        "Pladias_lifeform_original_values":dict(Counter(life[t] for t in known if t in life)),
        "Pladias_diaspore_original_values":dict(Counter(disp[t] for t in known if t in disp)),
        "source_sample_QC":species_qa,
        "reference_group_full_same_source_site_taxonomic_vs_functional_comparison":out,
        "this_is_observational_Czech_source_one_country":True,
        "natural_reference_land_verified":False,
        "natural_analogs_TH_TF_include_managed_secondary_vegetation":True,
        "rail_text_source_not_verified_active_ballast":True,
        "plot_independent_publication_does_not_imply_country_independence":True,
        "same_field_estimator_as_Croatia_Italy_rail_only_null":False,
        "Pladias_substantial_use_governance_contact_board_before_manuscript":True,
        "no_imputation_of_missing_traits":True,
        "no_macro_cross_country_or_causal_selection_claim":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--headers",required=True)
    p.add_argument("--species",required=True)
    p.add_argument("--lifeform-xlsx",required=True)
    p.add_argument("--diaspore-xlsx",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=run(a.headers,a.species,a.lifeform_xlsx,a.diaspore_xlsx)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    summary={}
    for trait,panel in result["reference_group_full_same_source_site_taxonomic_vs_functional_comparison"].items():
        summary[trait]={}
        for design,b in panel.items():
            x=b["responses"];values=x.get("summary_statistics",{})
            summary[trait][design]={
                "frame":b["original_source_frame"],
                "identified":x["identified"],
                "n_regions":x["n_original_phyto_regions"],
                "n_rail_source_plots":x["n_source_rail_sites"],
                "response":{k:values[k] for k in FOCUS_KEYS if k in values},
                "loo":{k:v for k,v in x.get("leave_one_region_out_sign_stability",{}).items()
                       if k.endswith(("rail_minus_analog","rail_minus_nonanalog"))}
            }
    print(json.dumps({
        "version":"v19.3","original_CzechVeg_source_rows":result["original_species_row_count"],
        "matched_lifeform_taxa":result["exact_original_taxon_match_lifeform"],
        "matched_diaspore_taxa":result["exact_original_taxon_match_diaspore"],
        "trait_source_row_counts":{
          "life":result["original_trait_source_metadata"]["source_Pladias_original_lifeform_rows"],
          "strategy":result["original_trait_source_metadata"]["source_Pladias_original_diaspore_rows"]},
        "all_response_summary":summary
    },ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()

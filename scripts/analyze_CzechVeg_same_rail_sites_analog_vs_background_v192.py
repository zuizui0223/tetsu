#!/usr/bin/env python3
"""v19.2: SAME RAIL SITES and SAME phyto-regions, only reference pool changes.

Validates whether a v19.1 sign difference survives strictly identical
source rail plots and regional covariate frames. Reference pools:
  (a) original source TH/TF/TG class analogue (incl secondary habitats);
  (b) original other nontransport-source class background EXCLUDING a;
  (c) all original nontransport-source background (a+b).
Each reference original plot is a different publisher study from its
paired rail-text original plot, region/decade/area-class identical.
Still one CzechVeg source country (text RAIL, not verified track habitat),
not experimental causality or meta-analytic field CI.
"""
from __future__ import annotations
import argparse,hashlib,json,random,statistics,sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    DOI,get_meta,matched_rail_source,ANALOG_CODES)
from analyze_CzechVeg_original_independent_reference_taxonomic_v191 import (
    source_species,stable_sample,jaccard,cross_similarity,genus_set,stats)
SEED=20261116
MIN_SOURCE_CONTROL=5
CAP_PER_RAIL=250
DRAWS=500
LOO_DRAWS=80
CATEGORIES=("analog","nonanalog","all")
METRICS=("species","genus")

def candidate_arms(meta,blocks):
    by_geo=defaultdict(list)
    for key,rec in meta.items():
        if rec["label"] is None and rec["strict"] is None:
            by_geo[(rec["region"],rec["decade"],rec["area"])].append(key)
    out={name:{} for name in CATEGORIES}
    railmeta={}
    for (paper,region,dec,area),v in sorted(blocks.items()):
        independent=[x for x in by_geo[(region,dec,area)]
                     if meta[x]["study"] != paper]
        cand={
            "analog":[x for x in independent if meta[x]["source_class"] in ANALOG_CODES],
            "nonanalog":[x for x in independent if meta[x]["source_class"] not in ANALOG_CODES],
            "all":independent
        }
        for rid in v["rail"]:
            for k in CATEGORIES:
                if len(cand[k])>=MIN_SOURCE_CONTROL:
                    out[k][rid]=stable_sample(cand[k],CAP_PER_RAIL,
                                              rid+"|v192|"+k)
                    railmeta[rid]=meta[rid]
    return out,railmeta

def common_frame(arms,meta,plants,richness=False):
    eligible={}
    for kind in CATEGORIES:
        eligible[kind]={}
        for rid,ids in arms[kind].items():
            if not plants.get(rid):continue
            options=[x for x in ids if plants.get(x)]
            if richness:
                n=len(plants[rid])
                options=[p for p in options if 2/3<=len(plants[p])/n<=1.5]
            if len(options)>=MIN_SOURCE_CONTROL:
                eligible[kind][rid]=options
    intersect=set.intersection(*(set(eligible[k]) for k in CATEGORIES))
    by_region=defaultdict(list)
    for rid in sorted(intersect):
        by_region[meta[rid]["region"]].append((rid,{k:eligible[k][rid] for k in CATEGORIES}))
    return dict(by_region),{
        "source_rail_sites_with_all_three_reference_kinds":len(intersect),
        "n_original_phyto_regions_with_all_reference_kinds":len(by_region),
        "source_rail_site_n_per_region_no_ids":{
            str(k):len(v) for k,v in sorted(by_region.items())},
        "richness_matched_independently_for_EACH_reference_kind":richness
    }

def sample_same_rail(by_region,plants,regions,rng):
    groups={"rail":[],"analog":[],"nonanalog":[],"all":[]}
    n_site={name:[] for name in groups}
    for region in regions:
        rail,cands=rng.choice(by_region[region])
        selected={"rail":rail}
        selected.update({kind:rng.choice(cands[kind]) for kind in CATEGORIES})
        for kind,pid in selected.items():
            groups[kind].append(plants[pid])
            n_site[kind].append(len(plants[pid]))
    out={}
    for taxrank in METRICS:
        transformed={k:(v if taxrank=="species" else [genus_set(s) for s in v])
                     for k,v in groups.items()}
        for k in transformed:out[taxrank+"_"+k]=cross_similarity(transformed[k])
        for k in CATEGORIES:
            out[taxrank+"_rail_minus_"+k]=out[taxrank+"_rail"]-out[taxrank+"_"+k]
        out[taxrank+"_analog_minus_nonanalog"]=out[taxrank+"_analog"]-out[taxrank+"_nonanalog"]
    out.update({"mean_richness_"+k:statistics.mean(vals) for k,vals in n_site.items()})
    return out

def run_pair_rail_reference(by_region,plants,rng_seed=SEED,draws=DRAWS):
    regions=sorted(by_region)
    result={"regions_count":len(regions),
            "source_rail_n":len({p for arr in by_region.values() for p,_ in arr}),
            "same_original_rail_plot_drawn_for_ALL_reference_choices":True,
            "same_original_phyto_region_between_all_pool_types":True,
            "n_pseudo_sampling_draws_not_field_replicates":draws,
            "nonrail_regional_reference_not_verified_natural":True}
    if len(regions)<3:
        return {**result,"identified":False,"reason":"Need three supported original phyto-regions"}
    rng=random.Random(rng_seed)
    sample=[sample_same_rail(by_region,plants,regions,rng) for _ in range(draws)]
    keys=sorted(sample[0])
    result["resampling_results"]={k:stats([x[k] for x in sample]) for k in keys}
    l=[]
    for i,region in enumerate(regions):
        smaller=[x for x in regions if x!=region]
        rr=random.Random(rng_seed+9000+i)
        sub=[sample_same_rail(by_region,plants,smaller,rr) for _ in range(LOO_DRAWS)]
        l.append({
            "removed_region":region,
            "species_rail_minus_analog":round(statistics.mean(x["species_rail_minus_analog"] for x in sub),7),
            "species_rail_minus_nonanalog":round(statistics.mean(x["species_rail_minus_nonanalog"] for x in sub),7),
            "genus_rail_minus_analog":round(statistics.mean(x["genus_rail_minus_analog"] for x in sub),7)
        })
    result["loo_n_region"]=len(l)
    result["loo_species_analog_positive"]=sum(x["species_rail_minus_analog"]>0 for x in l)
    result["loo_species_nonanalog_positive"]=sum(x["species_rail_minus_nonanalog"]>0 for x in l)
    result["leave_one_region_out_descriptive"]=l
    result["no_field_p_value_or_true_confidence_interval"]=True
    result["identified"]=True
    return result

def run(header,species):
    meta,hinfo=get_meta(header)
    blocks,binfo=matched_rail_source(meta)
    arms,_=candidate_arms(meta,blocks)
    read_ids=set()
    for pool in CATEGORIES:
        read_ids.update(arms[pool])
        for ids in arms[pool].values():read_ids.update(ids)
    plant,proof=source_species(species,read_ids)
    out={}
    for strict in (False,True):
        frame,qa=common_frame(arms,meta,plant,richness=strict)
        out["richness_matched" if strict else "metadata_only"]={
            "frame":qa,
            "comparison":run_pair_rail_reference(frame,plant,
                         rng_seed=SEED+(1001 if strict else 0))}
    return {
        "version":"v19.2",
        "original_DOI":DOI,
        "original_source_headers_md5":hinfo["publisher_original_headers_md5"],
        "original_source_species_md5":proof["publisher_species_csv_md5"],
        "full_original_species_records_scanned":proof["all_original_species_rows"],
        "source_plot_eligibility":binfo,
        "source_original_row_retrieval":proof,
        "strictly_identical_original_rail_plots_and_georegions_comparison":out,
        "main_source_rail_habitat_text_not_physically_verified":True,
        "natural_analog_TH_TF_TG_may_include_managed_secondary_vegetation":True,
        "all_three_reference_kinds_same_CzechVeg_original_compilation":True,
        "no_measured_functional_traits_or_CAUSAL_rail_effects":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--headers",required=True)
    p.add_argument("--species",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    out=run(a.headers,a.species)
    dest=Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    small={"version":"v19.2","original_rows":out["full_original_species_records_scanned"],
           "source_original_MD5":out["original_source_species_md5"],
           "comparisons":{}}
    for name,r in out["strictly_identical_original_rail_plots_and_georegions_comparison"].items():
        cm=r["comparison"]
        small["comparisons"][name]={"frame":r["frame"],"identified":cm["identified"],
          "loo_rail_vs_analog_positive":cm.get("loo_species_analog_positive"),
          "loo_rail_vs_nonanalog_positive":cm.get("loo_species_nonanalog_positive"),
          "metrics":{k:v for k,v in cm.get("resampling_results",{}).items()
                     if k.startswith(("species_rail_minus_","genus_rail_minus_",
                                       "species_analog_minus_nonanalog",
                                       "mean_richness_"))}}
    print(json.dumps(small,ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()

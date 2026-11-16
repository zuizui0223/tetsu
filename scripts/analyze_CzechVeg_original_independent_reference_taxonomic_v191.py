#!/usr/bin/env python3
"""v19.1: True CzechVeg plot SPECIES regional-reference first macro contrast.

Two independently sourced *plots / field publications*, all drawn from the
SAME published CzechVeg-Open dataset, DOI 10.5281/zenodo.19231152.
The reference is either (A) geographically matched nontransport-text source
plots, or (B) original Pladias syntaxon TH/TF/TG dry/pioneer analogous
vegetation. Reference context lacks direct verification of *natural* land
status. Neither is a proof of an unbiased local environmental species pool.
RAIL is original source TEXT rather than verified active trackbed habitat.

An observational CROSS-REGION Jaccard contrast, not train transport,
selection, convergent evolution, nor independent multi-country replication.
Original vegetation plot sources are not rehosted. Every original site x
species source row scanned. Controls are chosen without looking at species.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    DOI,get_meta,matched_rail_source,ANALOG_CODES
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import normalize_taxon

EXPECTED_SPECIES_BYTES=382939947
EXPECTED_SPECIES_MD5="b296e5e19f03fb3bad0a8347c317f8a3"
EXPECTED_ROWS=2510340
SEED=20261116
CAP_BACKGROUND_PER_RAIL=250
MIN_BACKGROUND=5
REPS=500
LEAVE_REGION_REPS=80

def stable_sample(ids,n,tag):
    return sorted(ids,key=lambda x:hashlib.sha256((str(tag)+"|"+x).encode()).digest())[:n]

def pick_source_candidates(meta,railblocks):
    bygeo=defaultdict(list)
    for sid,m in meta.items():
        if m["label"] is None and m["strict"] is None:
            bygeo[(m["region"],m["decade"],m["area"])].append(sid)
    arms={"nontransport":defaultdict(list),"analog":defaultdict(list)}
    for (study,reg,dec,area),cohort in sorted(railblocks.items()):
        source=[i for i in bygeo[(reg,dec,area)] if meta[i]["study"]!=study]
        for kind in ("nontransport","analog"):
            pool=source if kind=="nontransport" else [
                i for i in source if meta[i]["source_class"] in ANALOG_CODES]
            if len(pool)<MIN_BACKGROUND:continue
            for rail_id in cohort["rail"]:
                arms[kind][rail_id]=stable_sample(pool,CAP_BACKGROUND_PER_RAIL,rail_id+"|"+kind)
    to_read=set()
    for kind in arms:
        to_read.update(arms[kind])
        for pid_list in arms[kind].values():to_read.update(pid_list)
    return arms,to_read,{
        "candidate_rail_plots_after_geographical_sample_support":{k:len(v) for k,v in arms.items()},
        "original_sites_to_read_including_all_preselected_controls":len(to_read),
        "controls_preselected_by_original_site_id_hash_before_plant_species":True,
        "max_stable_reference_control_candidates_per_rail_plot":CAP_BACKGROUND_PER_RAIL
    }

def source_species(path,wanted):
    p=Path(path)
    if p.stat().st_size!=EXPECTED_SPECIES_BYTES:
        raise ValueError("Publisher original 2.51M species source byte count differs")
    md=hashlib.md5()
    with p.open("rb") as handle:
        for b in iter(lambda:handle.read(8*1024*1024),b""):
            md.update(b)
    if md.hexdigest()!=EXPECTED_SPECIES_MD5:
        raise ValueError("Published original CzechVeg species CSV MD5 differs")
    got=defaultdict(set);c=Counter()
    with p.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if {"RELEVE_NR","KAPLAN_PLUS_NAME_SIMPLE","PLANT_GROUP"}-set(reader.fieldnames or []):
            raise ValueError("Not original source taxon names or vascular group")
        for row in reader:
            c["all_original_species_rows"]+=1
            id=str(row.get("RELEVE_NR") or "").strip()
            if id not in wanted:continue
            if str(row.get("PLANT_GROUP") or "").strip().lower()!="vascular":continue
            name=normalize_taxon(row.get("KAPLAN_PLUS_NAME_SIMPLE"))
            if name and name.lower() not in {"na","nan","null"}:
                got[id].add(name)
                c["selected_original_source_vascular_rows"]+=1
    if c["all_original_species_rows"]!=EXPECTED_ROWS:
        raise ValueError("Source complete 2.51034M original species×plot rows not read")
    return got,{
        **c,"publisher_species_csv_md5":md.hexdigest(),
        "source_sites_with_real_vascular_taxa":len(got),
        "source_sites_with_no_real_plant_rows":len(wanted)-len(got)
    }

def regional_frame(arms,plant,meta,richness_gate=False):
    bykind={}
    for kind,control in arms.items():
        per_region=defaultdict(list)
        for rail_id,refs in control.items():
            if rail_id not in plant or not plant[rail_id]:continue
            cand=[pid for pid in refs if pid in plant and plant[pid]]
            if richness_gate:
                target=len(plant[rail_id])
                # Restrict to rough x1.5 richness match without inspecting
                # plant identities (higher-order sampling conditioning).
                cand=[pid for pid in cand if 2/3 <=len(plant[pid])/target<=1.5]
            if len(cand)<MIN_BACKGROUND:continue
            reg=meta[rail_id]["region"]
            per_region[reg].append((rail_id,cand))
        bykind[kind]=dict(per_region)
    return bykind

def genus_set(taxa):
    return {s.split(" ",1)[0] for s in taxa if s and s[0].isupper()}

def jaccard(a,b):
    union=a|b
    return (len(a&b)/len(union)) if union else float("nan")

def cross_similarity(sets):
    values=[]
    for i in range(len(sets)):
        for j in range(i+1,len(sets)):
            values.append(jaccard(sets[i],sets[j]))
    return statistics.mean(values)

def sample_once(source,plant,rng,regions):
    rail,ref=[],[]
    count_r,count_n=[],[]
    for reg in regions:
        railid,pool=rng.choice(source[reg])
        bgid=rng.choice(pool)
        rail.append(plant[railid]);ref.append(plant[bgid])
        count_r.append(len(plant[railid]));count_n.append(len(plant[bgid]))
    if len(rail)<3:raise ValueError("Expected >=3 geographic original regions")
    return {
        "rail_species_jaccard":cross_similarity(rail),
        "reference_species_jaccard":cross_similarity(ref),
        "rail_genus_jaccard":cross_similarity([genus_set(s) for s in rail]),
        "reference_genus_jaccard":cross_similarity([genus_set(s) for s in ref]),
        "rail_mean_n_source_taxa":statistics.mean(count_r),
        "reference_mean_n_source_taxa":statistics.mean(count_n)
    }

def stats(values):
    v=sorted(values)
    return {"mean":round(statistics.mean(v),7),
            "q025_resampling_not_field_CI":round(v[int(.025*len(v))],7),
            "q975_resampling_not_field_CI":round(v[min(len(v)-1,int(.975*len(v)))],7),
            "min":round(v[0],7),"max":round(v[-1],7)}

def run_reference_panel(by_region,plant,label,seed=SEED,reps=REPS):
    regions=sorted(by_region)
    sites=set(r for entries in by_region.values() for r,_ in entries)
    candidates={c for entries in by_region.values() for _,controls in entries for c in controls}
    out={"reference_label":label,
         "n_independent_phyto_regions":len(regions),
         "n_eligible_original_rail_text_plots":len(sites),
         "n_independent_original_reference_candidate_plots":len(candidates),
         "observational_source_only":True}
    if len(regions)<3:
        return {**out,"identified":False,"reason":"Fewer than 3 independent source phytogeographical regions"}
    rng=random.Random(seed)
    vals=[sample_once(by_region,plant,rng,regions) for _ in range(reps)]
    measures={
       "species_delta":[a["rail_species_jaccard"]-a["reference_species_jaccard"] for a in vals],
       "genus_delta":[a["rail_genus_jaccard"]-a["reference_genus_jaccard"] for a in vals],
       "rail_species_jaccard":[a["rail_species_jaccard"] for a in vals],
       "reference_species_jaccard":[a["reference_species_jaccard"] for a in vals],
       "rail_genus_jaccard":[a["rail_genus_jaccard"] for a in vals],
       "reference_genus_jaccard":[a["reference_genus_jaccard"] for a in vals],
       "rail_mean_source_richness":[a["rail_mean_n_source_taxa"] for a in vals],
       "reference_mean_source_richness":[a["reference_mean_n_source_taxa"] for a in vals],
    }
    summaries={k:stats(v) for k,v in measures.items()}
    leave={}
    for i,region in enumerate(regions):
        rest=[r for r in regions if r!=region]
        if len(rest)<3:continue
        rngloo=random.Random(seed+1030+i)
        samples=[sample_once(by_region,plant,rngloo,rest) for _ in range(LEAVE_REGION_REPS)]
        leave[region]={
           "species_delta":round(statistics.mean(x["rail_species_jaccard"]-
                                   x["reference_species_jaccard"] for x in samples),7),
           "genus_delta":round(statistics.mean(x["rail_genus_jaccard"]-
                                   x["reference_genus_jaccard"] for x in samples),7)
        }
    return {
       **out,"identified":True,"n_pseudo_resampling_draws":reps,
       "n_nonindependent_cross_region_pairs_each_draw":len(regions)*(len(regions)-1)//2,
       "source_sample_draws_not_independent_field_studies":True,
       "cross_region_species_or_genus_similarity_difference_positive_means_more_rail_similarity":True,
       "response_statistics":summaries,
       "region_leave_one_out_descriptive":leave,
       "n_leave_one_out_regions_with_positive_species_difference":
            sum(x["species_delta"]>0 for x in leave.values()),
       "n_leave_one_out_regions_with_positive_genus_difference":
            sum(x["genus_delta"]>0 for x in leave.values()),
       "source_regions_have_potential_nonuniform_railway_survey_type":True,
       "source_study_clusters_not_estimated_as_independent_sampling_distribution":True
    }

def run(header,species):
    meta,headqa=get_meta(header)
    blocks,frame=matched_rail_source(meta)
    arms,read_ids,selection=pick_source_candidates(meta,blocks)
    plants,plantqa=source_species(species,read_ids)
    results={}
    for kind in ("nontransport","analog"):
        results[kind]={}
        for richness in (False,True):
            tag="richness_1p5_matched" if richness else "area_year_region_only"
            frame_r=regional_frame({kind:arms[kind]},plants,meta,richness_gate=richness)[kind]
            results[kind][tag]=run_reference_panel(
                frame_r,plants,kind+"|"+tag,
                seed=SEED+int(richness)*20002+(0 if kind=="analog" else 1003))
    return {
       "version":"v19.1","original_CzechVeg_DOI":DOI,
       "source_original_headers_md5":headqa["publisher_original_headers_md5"],
       "source_original_species_md5":plantqa["publisher_species_csv_md5"],
       "data_original_117739_headers_and_2510340_species_rows_verified":True,
       "source_only_rail_vs_road_matched":frame,
       "preselected_controls":selection,
       "actual_source_species_matrix_quality":plantqa,
       "true_species_vs_genus_cross_region_matched_source_reference":results,
       "actual_natural_land_status_verified":False,
       "source_syntaxon_TH_TF_TG_analog_is_natural_or_secondary_open_vegetation":True,
       "single_country_one_common_compiled_database":True,
       "original_source_rail_text_plots_not_experimental_ballast_habitat":True,
       "functional_lifeform_traits_not_tested":True,
       "statistical_CI_p_from_independent_papers_not_estimated":True,
       "interpretation":"Source-conditional exploratory multi-region taxonomic comparison; study-dependent botanical survey plots and original vegetation classes do not identify rail causality or global convergence."}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--headers",required=True)
    ap.add_argument("--species",required=True)
    ap.add_argument("--output",required=True)
    p=ap.parse_args()
    out=run(p.headers,p.species)
    dest=Path(p.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    short={}
    for pool,panel in out["true_species_vs_genus_cross_region_matched_source_reference"].items():
        short[pool]={}
        for design,stat in panel.items():
            short[pool][design]=({
               "identified":stat["identified"],
               "n_rail":stat["n_eligible_original_rail_text_plots"],
               "n_regions":stat["n_independent_phyto_regions"],
               "n_reference":stat["n_independent_original_reference_candidate_plots"],
               "responses":{k:v for k,v in stat.get("response_statistics",{}).items()
                            if k in ("species_delta","genus_delta","rail_mean_source_richness",
                                     "reference_mean_source_richness")},
               "loo_positive":stat.get("n_leave_one_out_regions_with_positive_species_difference")})
    print(json.dumps({"version":"v19.1","real_source_original_full":True,
                      "source_plot_taxon_quality":out["actual_source_species_matrix_quality"],
                      "reference_analyses":short},ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()

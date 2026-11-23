#!/usr/bin/env python3
"""v343 independent EIVE thermal-indicator source study benchmark in CzechVeg.

Original RAIL-context versus ROAD-context historical vegetation labels.
Same original author/paper, phytoregion, decade, area-band matching, but NOT
field-matched nonrail nature or active ballast. No raw plot records or IDs
exported. Main design frozen in docs/tetsu_v343..._frozen_ja.md.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random,re,sys,unicodedata
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median
sys.path.insert(0,str(Path(__file__).resolve().parent))

from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file,plot_labels,plot_taxa,HEADERS_NAME,SPECIES_NAME)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import header_source_context
from analyze_CzechVeg_rail_road_same_study_trait_v123 import (
    source_study_meta,common_original_study_blocks,SPECS)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_C4_mechanism_v339 import C4

SOURCE_PLOTS=117739
SOURCE_SPECIES_ROWS=2510340
SPEC=next(x for x in SPECS if x[0]=="S2_study_region_decade_area4_min2")
PANEL=((3,0.5),(3,0.7),(1,0.5))
MIN_SOURCE_TAXA=5
MIN_EIVE_TAXA=3
MIN_PLOTS_PER_ARM=2
N_NULL=999
SEED=343
RE_BINOMIAL=re.compile(r"^[a-z-]+ [a-z-]+$")
REMOVAL=(("full",frozenset()),("exclude_Cynodon",frozenset(("cynodon dactylon",))),
         ("exclude_C4_candidates",C4))


def normalized_binomial(v):
    s=unicodedata.normalize("NFKC",str(v or "")).strip().casefold()
    s=" ".join(s.split())
    return s if RE_BINOMIAL.fullmatch(s) else None


def source_cohort(zenodo,headers,species,eive_path):
    meta=json.loads(Path(zenodo).read_text(encoding="utf8"))
    q1=verify_file(meta,HEADERS_NAME,headers)
    q2=verify_file(meta,SPECIES_NAME,species)
    eive,eive_qa=eive_table(eive_path)
    labelled,plot_qa=plot_labels(headers)
    year,year_qa=header_source_context(headers,labelled)
    original_meta,paper_qa=source_study_meta(headers,year)
    plants,plants_qa=plot_taxa(species,original_meta,{}, {})
    if (plot_qa["counts"]["original_header_rows"]!=SOURCE_PLOTS or
        plants_qa["n_original_source_species_csv_records_scanned"]!=SOURCE_SPECIES_ROWS):
        raise ValueError("CzechVeg source record counts not identical to prior source QA")
    blocks,block_qa=common_original_study_blocks(original_meta,plants,"broad",SPEC)
    counts=(len(blocks),len({k[0] for k in blocks}),
            sum(len(v["rail"]) for v in blocks.values()),
            sum(len(v["road"]) for v in blocks.values()))
    if counts!=(44,36,229,462):
        raise ValueError(f"Previously source frozen 44/36/229/462 original design drift: {counts}")
    traits={}
    for name,entry in eive.items():
        traits[name]=(float(entry["position"]),entry["systems"])
    site_names={}
    for pid in {p for arms in blocks.values() for v in arms.values() for p in v}:
        source_taxa=plants[pid]
        x=set(filter(None,(normalized_binomial(q[0]) for q in source_taxa)))
        site_names[pid]={"source_total_taxa":len(source_taxa),
                         "exact_binomial_taxa":x}
    qa={"CzechVeg_headers":q1,"CzechVeg_species":q2,
        "EIVE":eive_qa,
        "original_publisher_source_quality":paper_qa,
        "original_original_matched_blocks":block_qa,
        "original_source_counts":{
            "header_plots":SOURCE_PLOTS,
            "plant_rows":SOURCE_SPECIES_ROWS,
            "matched_blocks":44,"source_studies":36,"rail_text_plots":229,
            "road_text_plots":462},
        "original_label_area":plot_qa["counts"],
        "CzechVeg_source_original_geographic_and_year":year_qa,
        "source_EIVE_crosswalk_explicit_binomial_only":True}
    return blocks,site_names,traits,qa


def plot_score(site,trait_map,minsystems,minimum_cov):
    original_n=site["source_total_taxa"]
    eligible=[n for n in site["exact_binomial_taxa"]
              if n in trait_map and trait_map[n][1]>=minsystems]
    matched_n=len(eligible)
    coverage=matched_n/original_n if original_n else 0.0
    if original_n<MIN_SOURCE_TAXA or matched_n<MIN_EIVE_TAXA or coverage<minimum_cov:
        return None
    return {"eligible":tuple(sorted(eligible)),
            "baseline_EIVE":mean(trait_map[n][0] for n in eligible),
            "matched_n":matched_n,"source_n":original_n,"coverage":coverage}


def frozen_frame(blocks,site_names,trait_map,minsystems,minimum_cov):
    plot_s={p:plot_score(s,trait_map,minsystems,minimum_cov)
            for p,s in site_names.items()}
    selected={}
    for k,arms in blocks.items():
        rail=tuple(pid for pid in arms["rail"] if plot_s[pid] is not None)
        road=tuple(pid for pid in arms["road"] if plot_s[pid] is not None)
        if min(len(rail),len(road))>=MIN_PLOTS_PER_ARM:
            selected[k]={"rail":rail,"road":road}
    valid_plots={pid for arm in selected.values()
                for z in ("rail","road") for pid in arm[z]}
    pub={k[0] for k in selected}
    plot_cov=[plot_s[p]["coverage"] for p in valid_plots]
    return selected,plot_s,{
        "EIVE_independent_sources_threshold":minsystems,
        "min_species_coverage_each_source_plot":minimum_cov,
        "n_selected_original_field_studies":len(pub),
        "n_selected_original_matched_blocks":len(selected),
        "n_selected_original_RAIL_text_plots":sum(len(v["rail"]) for v in selected.values()),
        "n_selected_original_ROAD_text_plots":sum(len(v["road"]) for v in selected.values()),
        "mean_source_plot_trait_coverage":mean(plot_cov) if plot_cov else None,
        "min_source_plot_trait_coverage":min(plot_cov) if plot_cov else None,
        "same_source_study_region_decade_area4_matching_before_trait_support":True,
        "original_source_not_a_randomized_field_study":True}


def report_stats(selected,plot_s,traits,excluded=frozenset()):
    paper_blocks=defaultdict(list)
    missed=Counter()
    for k,v in selected.items():
        obs={}
        for side in ("rail","road"):
            current=[]
            for pid in v[side]:
                terms=[traits[name] for name in plot_s[pid]["eligible"]
                       if name not in excluded and name in traits]
                if terms:current.append(mean(terms))
                else:missed["plots_with_no_qualifying_traits_after_deletion"]+=1
            if len(current)<MIN_PLOTS_PER_ARM:
                missed["blocks_disqualified_due_to_unscorable_original_plot"]+=1
                obs=None
                break
            obs[side]=mean(current)
        if obs is not None:
            paper_blocks[k[0]].append(obs["rail"]-obs["road"])
    paper={p:mean(v) for p,v in paper_blocks.items()}
    vals=list(paper.values())
    if not vals:
        return {"n_original_studies":0,"n_original_blocks":0,
                "n_source_plots_missing_post_deletion":missed.get("plots_with_no_qualifying_traits_after_deletion",0),
                "no_confirmed_evaluable_external_benchmark":True},None
    out={"n_original_studies":len(vals),
         "n_original_matched_blocks":sum(len(x) for x in paper_blocks.values()),
         "n_original_studies_positive_RAIL_minus_ROAD":sum(x>0 for x in vals),
         "n_original_studies_negative_RAIL_minus_ROAD":sum(x<0 for x in vals),
         "mean_study_equal_RAIL_minus_ROAD_EIVE":mean(vals),
         "median_study_RAIL_minus_ROAD_EIVE":median(vals),
         "min_study_thermal_contrast":min(vals),"max_study_thermal_contrast":max(vals),
         "n_source_plots_missing_post_deletion":missed.get("plots_with_no_qualifying_traits_after_deletion",0),
         "n_original_blocks_dropped_post_deletion":missed.get("blocks_disqualified_due_to_unscorable_original_plot",0),
         "original_study_hashes_only_no_original_paper_plot_ids_in_results":True}
    return out,paper


def quantify(original,bs=None,ns=1999):
    vals=list(original.values())
    if len(vals)<3:return None
    n=len(vals)
    rng=random.Random(SEED)
    boot=[mean(vals[rng.randrange(n)] for _ in vals) for _ in range(ns)]
    loo=[mean(vals[:i]+vals[i+1:]) for i in range(n)]
    return {"n_source_studies":n,"study_cluster_resampling_draws":ns,
            "resample_q025_q975_technical_not_randomized_field_CI":
               [percentile(boot,.025),percentile(boot,.975)],
            "leave_one_study_out_min_max":[min(loo),max(loo)],
            "all_leave_one_paper_out_positive":all(x>0 for x in loo),
            "all_leave_one_paper_out_negative":all(x<0 for x in loo)}


def percentile(vals,p):
    v=sorted(vals);q=(len(v)-1)*p;a=math.floor(q);b=math.ceil(q)
    return v[a]+(v[b]-v[a])*(q-a)


def trait_null(selected,plot_s,trait_map,obs,mode,draws=N_NULL):
    active={n for v in selected.values() for side in ("rail","road") for pid in v[side]
            for n in plot_s[pid]["eligible"]}
    vals={name:trait_map[name] for name in active}
    genus=defaultdict(list)
    for name in active:
        genus[name.split()[0]].append(name)
    exchange=[sorted(names) for names in genus.values()
              if len(names)>1 and len({vals[x] for x in names})>1]
    if mode=="genus" and (sum(map(len,exchange))<20 or len(exchange)<8):
        return {"mode":mode,"n_active_taxa":len(vals),
                "n_movable_genera":len(exchange),"n_movable_taxa":sum(map(len,exchange)),
                "eligible_for_genus_shuffle":False}
    if mode not in ("global","genus"):raise ValueError("Unknown trait shuffle mode")
    rng=random.Random(SEED+ (0 if mode=="global" else 300))
    names=sorted(active)
    original_vals=[vals[x] for x in names]
    simulated=[]
    for _ in range(draws):
        perm=vals.copy()
        if mode=="global":
            y=original_vals.copy();rng.shuffle(y)
            perm=dict(zip(names,y))
        else:
            for gp in exchange:
                y=[vals[x] for x in gp];rng.shuffle(y)
                perm.update(zip(gp,y))
        outcome,paper=report_stats(selected,plot_s,perm)
        if outcome["n_original_studies"]!=obs["n_original_studies"]:
            raise AssertionError("Trait shuffle changed original plot completeness/matched studies")
        simulated.append(outcome["mean_study_equal_RAIL_minus_ROAD_EIVE"])
    mu=mean(simulated)
    target=obs["mean_study_equal_RAIL_minus_ROAD_EIVE"]
    count=sum(abs(x-mu)>=abs(target-mu)-1e-12 for x in simulated)
    return {"mode":mode,"draws":draws,
            "n_source_taxa_redistributed":len(vals) if mode=="global" else sum(map(len,exchange)),
            "original_source_incidence_and_trait_known_masks_fixed":True,
            "actual_mean_study_difference":target,"null_mean":mu,
            "null_95pct_technical":[percentile(simulated,.025),percentile(simulated,.975)],
            "two_sided_fraction_centered_plus_one":(count+1)/(draws+1),
            "technical_trait_null_not_original_field_p_value":True}


def analyze(zenodo,headers,species,eive):
    blocks,site_names,trait_data,qa=source_cohort(zenodo,headers,species,eive)
    out={"version":"v343-CzechVeg-independent-original-study-EIVE-RAILtext-vs-ROADtext",
         "source_QA":qa,"panels":[],
         "not_same_site_rail_natural_reference":True,
         "not_independent_country_or_climate_region_validation":True,
         "original_Czech_source_field_habitat_is_text_only":True,
         "no_original_source_data_rows_or_raw_GPS_republished":True}
    for systems,coverage in PANEL:
        selected,scores,frame=frozen_frame(blocks,site_names,trait_data,systems,coverage)
        panel={"min_EIVE_systems":systems,"coverage_min":coverage,
               "source_selected_frame":frame,"deletion_effects":[],
               "source_caveats":[
                  "Original rail and road source plots share paper/region/year/area class, not same GPS",
                  "Historical site descriptions are ambiguous, train operating/ballast conditions unverified",
                  "Field study pseudoreplication and nonindependence persist",
                  "EIVE trait is European ordinal preference, not measured thermal microclimate",
                  "Independent regional climate interaction and Australia/Kashmir original contrasts are not identified",
                  "A technical taxon label shuffle does not provide randomized field p values"
               ]}
        for label,exc in REMOVAL:
            metrics,paper=report_stats(selected,scores,trait_data_to_positions(trait_data),exc)
            if paper is not None and label=="full":
                metrics["source_study_spatial_resampling"]=quantify(paper)
                metrics["source_taxon_label_nulls"]=[
                    trait_null(selected,scores,trait_data_to_positions(trait_data),metrics,"global"),
                    trait_null(selected,scores,trait_data_to_positions(trait_data),metrics,"genus")
                ] if systems==3 and coverage==0.5 else []
            panel["deletion_effects"].append({"source_taxon_subset":label,**metrics})
        out["panels"].append(panel)
    return out

def trait_data_to_positions(traits):
    return {k:v[0] for k,v in traits.items()}

def main():
    p=argparse.ArgumentParser()
    for field in ("zenodo-record","headers","species","eive","output"):
        p.add_argument("--"+field,required=True)
    a=p.parse_args()
    r=analyze(a.zenodo_record,a.headers,a.species,a.eive)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(r,ensure_ascii=False),flush=True)

if __name__=="__main__":main()

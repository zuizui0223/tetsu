#!/usr/bin/env python3
"""v341: source-authentic 2014→2017 re-observation stratified by EIVE warm trait.

Temperature preference EIVE is not measured local heat tolerance.
Each original site/zone compares 2014 present taxa observed again vs not
recorded again in 2017. No source GPS, species×site records or site names
are written to reports. See frozen docs before outcome.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_cover_weighted_matched_temperature_v334 import read_original_plant_values
from analyze_Kashmir_C4_mechanism_v339 import C4
from analyze_Kashmir_Cynodon_transition_v340 import source_sites,quantile,CYNODON

PRIMARY=("Track","Interior Plot 2")
NEAR=("Track","Interior Plot 1")
REMOVALS={
  "all_EIVE_matched_2014_taxa":frozenset(),
  "exclude_Cynodon_dactylon":frozenset((CYNODON,)),
  "exclude_all_C4_candidates":C4,
}
DRAWS=1999


def cohort(source_abundance,names,zones):
    eligible,exclusion=source_sites(source_abundance,zones)
    rows=[]
    for site in eligible:
        byzone={}
        for z in zones:
            a=source_abundance[site,2014,z]
            b=source_abundance[site,2017,z]
            byzone[z]=[(names[i],i in b,float(v)) for i,v in sorted(a.items()) if names[i]]
        rows.append(byzone)
    return rows,{"n_same_source_sites":len(rows),
                 "n_dropped_original_spatial_sites":31-len(rows),
                 "sites_with_at_least_one_missing_source_cell":exclusion}


def site_trait_split(entries,trait_map,exclude=frozenset()):
    old=[(float(trait_map[sp]),kept,math.log(ab))
         for sp,kept,ab in entries if sp in trait_map and sp not in exclude]
    kept=[x for x,flag,ab in old if flag]
    lost=[x for x,flag,ab in old if not flag]
    amounts_kept=[ab for x,flag,ab in old if flag]
    amounts_lost=[ab for x,flag,ab in old if not flag]
    return {"n_2014_EIVE_matched":len(old),"n_2017_reobserved":len(kept),
            "n_2017_not_reobserved":len(lost),
            "observed_thermal_alignment_kept_minus_lost":
                mean(kept)-mean(lost) if kept and lost else None,
            "initial_log_source_abundance_kept_minus_lost":
                mean(amounts_kept)-mean(amounts_lost)
                  if amounts_kept and amounts_lost else None}


def panel_stats(cohorts,zones,trait_map,exclude=frozenset(),include_extra=True):
    byzone={z:[] for z in zones}
    matched=[]
    total=Counter()
    for sample in cohorts:
        row={}
        for z in zones:
            rec=site_trait_split(sample[z],trait_map,exclude)
            row[z]=rec
            byzone[z].append(rec)
            total["source_2014_EIVE_matched_species_site_occurrences_"+z]+=rec["n_2014_EIVE_matched"]
            total["source_2014_EIVE_matched_2017_reobserved_"+z]+=rec["n_2017_reobserved"]
        x=row[zones[0]]["observed_thermal_alignment_kept_minus_lost"]
        y=row[zones[1]]["observed_thermal_alignment_kept_minus_lost"]
        if x is not None and y is not None:
            matched.append(x-y)
    results={"source_complete_site_count":len(cohorts),
             "n_informative_temperature_alignment_site_pairs":len(matched),
             "n_source_sites_excluded_from_joint_stat_due_to_no_retained_or_no_lost_EIVE_species":
                 len(cohorts)-len(matched),
             "mean_matched_track_minus_interior_thermal_persistence_alignment":
                 mean(matched) if matched else None,
             "joint_site_signs":{
                 "track_more_negative":sum(v<0 for v in matched),
                 "track_more_positive":sum(v>0 for v in matched),
                 "same":sum(v==0 for v in matched)
             },
             "zone_stats":{},
             "total_2014_EIVE_species_occurrence_denominators":dict(total)}
    for z in zones:
        yes=[r["observed_thermal_alignment_kept_minus_lost"] for r in byzone[z]
             if r["observed_thermal_alignment_kept_minus_lost"] is not None]
        a=[r["initial_log_source_abundance_kept_minus_lost"] for r in byzone[z]
           if r["initial_log_source_abundance_kept_minus_lost"] is not None]
        nt=sum(r["n_2014_EIVE_matched"] for r in byzone[z])
        nr=sum(r["n_2017_reobserved"] for r in byzone[z])
        results["zone_stats"][z]={
             "n_informative_original_site_plots":len(yes),
             "mean_thermal_EIVE_of_reobserved_minus_not_reobserved":mean(yes) if yes else None,
             "n_plot_alignments_negative":sum(v<0 for v in yes),
             "n_plot_alignments_positive":sum(v>0 for v in yes),
             "n_plot_alignments_equal":sum(v==0 for v in yes),
             "2014_species_site_record_count":nt,
             "same_species_same_site_reobservation_count":nr,
             "reobservation_fraction_with_EIVE_2014_baseline":nr/nt if nt else None,
             "mean_initial_log_source_abundance_reobserved_minus_not_reobserved":
                 mean(a) if a else None,
             "n_plots_with_abundance_diagnostic":len(a),
        }
    if include_extra and len(matched)>1:
        rng=random.Random(342)
        b=[mean(matched[rng.randrange(len(matched))] for _ in matched)
             for _ in range(DRAWS)]
        result_boot={
           "n_original_spatial_site_pairs":len(matched),
           "site_resample_draws":DRAWS,
           "percentiles_025_975_technical_not_independent_field_CI":
               [quantile(b,.025),quantile(b,.975)],
           "leave_one_original_site_pair_out_min_max":
               [min(mean(matched[:i]+matched[i+1:]) for i in range(len(matched))),
                max(mean(matched[:i]+matched[i+1:]) for i in range(len(matched)))],
        }
        results["site_resampling_diagnostics"]=result_boot
    return results


def shuffle_traits(trait_map,mode,rng):
    nxt=dict(trait_map)
    if mode=="all_source_matched_taxa":
        names=sorted(trait_map)
        values=[trait_map[x] for x in names]
        rng.shuffle(values)
        nxt.update(dict(zip(names,values)))
    elif mode=="within_genus_only":
        gen=defaultdict(list)
        for name in trait_map:
            gen[name.split()[0]].append(name)
        for names in gen.values():
            names=sorted(names)
            if len(names)<2:continue
            vals=[trait_map[x] for x in names]
            rng.shuffle(vals)
            nxt.update(dict(zip(names,vals)))
    else:raise ValueError("Unknown frozen technical null")
    return nxt


def null_diagnostic(rows,zones,trait_map,exclude,observed,draws=DRAWS):
    if observed["n_informative_temperature_alignment_site_pairs"]<10:
        return {"status":"under10_matched_informative_original_sites_not_inferential"}
    modes=("all_source_matched_taxa","within_genus_only")
    group=defaultdict(list)
    for name in trait_map:group[name.split()[0]].append(name)
    movable={x for group_names in group.values() if len(group_names)>1
             and len({trait_map[n] for n in group_names})>1 for x in group_names}
    out={"n_original_taxa_with_EIVE":len(trait_map),
         "n_original_taxa_movable_in_genus_shuffle":len(movable),
         "n_original_genera_movable_in_genus_shuffle":
            sum(len(group_names)>1 and len({trait_map[n] for n in group_names})>1
                for group_names in group.values()),
         "nulls":[]}
    obs=observed["mean_matched_track_minus_interior_thermal_persistence_alignment"]
    for mode,seed in ((modes[0],341),(modes[1],1341)):
        rng=random.Random(seed)
        samples=[]
        for _ in range(draws):
            shuffled=shuffle_traits(trait_map,mode,rng)
            res=panel_stats(rows,zones,shuffled,exclude,False)
            assert res["n_informative_temperature_alignment_site_pairs"] == observed["n_informative_temperature_alignment_site_pairs"]
            samples.append(res["mean_matched_track_minus_interior_thermal_persistence_alignment"])
        mu=mean(samples)
        cnt=sum(abs(x-mu)>=abs(obs-mu)-1e-12 for x in samples)
        low=sum(x<=obs+1e-12 for x in samples)
        out["nulls"].append({
            "mode":mode,"draws":draws,"seed":seed,"observed":obs,
            "null_mean":mu,"null_95pct":[quantile(samples,.025),quantile(samples,.975)],
            "technical_two_sided_fraction_centered_plus_one":(cnt+1)/(draws+1),
            "negative_direction_one_sided_fraction_plus_one":(low+1)/(draws+1),
            "original_taxon_trait_shuffle_not_field_randomization":True
        })
    return out


def analyze(source_path,eive_path):
    presence,sourceqa=parse_xlsx_matrix(Path(source_path))
    names=species_names_from_publisher(source_path)
    source_ab,abqa=read_original_plant_values(source_path,presence)
    eive,eiveqa=eive_table(eive_path)
    if set(GLOBAL_SOURCE_NAMES for GLOBAL_SOURCE_NAMES in names if GLOBAL_SOURCE_NAMES) & set(eive)==set():
        raise ValueError("No exact source EIVE traits")
    result={"version":"v341-EIVE-temperature-affinity-observational-reobservation",
            "source_Dryad_SHA256":sourceqa["original_file_sha256"],
            "EIVE_independent_original_MD5":eiveqa["EIVE_original_md5"],
            "panels":[],
            "limitations":[
                "post-v340 exploratory temporal trait test, not prospectively confirmatory",
                "European EIVE transferred to Asian species; does not quantify physiological heat tolerance",
                "2014 presence to 2017 re-observation is NOT individual survival or true population establishment",
                "Raw original plot observations and original site/GPS are never reprinted",
                "No train seed transport, ballast temperature, in situ reproduction, competitive replacement or exclusion can be inferred",
                "Species and site are repeat-observed, neither are independent replicates; technical global/within-genus trait permutations not field experiments",
                "Source abundances from 1m² quadrats do not have fully verified physical units",
            ]}
    for zones in (PRIMARY,NEAR):
        rows,qa=cohort(source_ab,names,zones)
        for nsystems in (3,1):
            trait_map={s:float(eive[s]["position"]) for s in sorted(set(n for n in names if n))
                       if s in eive and eive[s]["systems"]>=nsystems}
            pan={"original_source_zone_comparison":list(zones),
                 "source_design_completeness":qa,
                 "EIVE_temperature_min_independent_systems":nsystems,
                 "EIVE_known_species_in_original_taxon_universe":len(trait_map),
                 "exclusion_results":[]}
            for label,excluded in REMOVALS.items():
                vals=panel_stats(rows,zones,trait_map,excluded,True)
                out={"original_species_subset":label,
                     "analysis":vals}
                # Frozen primary only; sensitivity descriptive to limit multiple testing.
                if zones==PRIMARY and nsystems==3:
                    out["frozen_technical_trait_label_nulls"]=null_diagnostic(
                        rows,zones,trait_map,excluded,vals)
                pan["exclusion_results"].append(out)
            result["panels"].append(pan)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=analyze(a.kashmir_original,a.eive_original)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

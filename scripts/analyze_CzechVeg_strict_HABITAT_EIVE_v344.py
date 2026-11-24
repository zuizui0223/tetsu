#!/usr/bin/env python3
"""v344: stricter original HABITAT (not LOCALITY) source labels, era bias.

Source genuine CzechVeg 1.1 and independent EIVE original MD5.
Never publish original plot coordinates, source species rows, or paper labels.
The output is about historical vegetation context, NOT causally sampled ballast.
Frozen design: docs/tetsu_v344_CzechVeg_strict_HABITAT_and_era_frozen_ja.md
"""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path

from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file,plot_labels,plot_taxa,HEADERS_NAME,SPECIES_NAME)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    header_source_context)
from analyze_CzechVeg_rail_road_same_study_trait_v123 import (
    source_study_meta,common_original_study_blocks)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_CzechVeg_original_study_EIVE_thermal_v343 import (
    SPEC,MIN_PLOTS_PER_ARM,normalized_binomial,frozen_frame,
    report_stats,quantify,trait_data_to_positions,REMOVAL
)

SOURCE_HEADERS=117739
SOURCE_ROWS=2510340
PANELS=((3,.5),(3,.7))
ERA=("before_2000","year_2000_onward","year_2010_onward")


def source_context(zenodo_path,headers_path,species_path,eive_path):
    zenodo=json.loads(Path(zenodo_path).read_text(encoding="utf8"))
    h=verify_file(zenodo,HEADERS_NAME,headers_path)
    s=verify_file(zenodo,SPECIES_NAME,species_path)
    eive,eiveqa=eive_table(eive_path)
    labels,qa=plot_labels(headers_path)
    time,qa_time=header_source_context(headers_path,labels)
    meta,qa_studies=source_study_meta(headers_path,time)
    plants,qa_plants=plot_taxa(species_path,meta,{}, {})
    if (qa["counts"]["original_header_rows"]!=SOURCE_HEADERS
          or qa_plants["n_original_source_species_csv_records_scanned"]!=SOURCE_ROWS):
        raise ValueError("Original CzechVeg source has unexpected plot/plant row counts")
    broad,broad_qa=common_original_study_blocks(meta,plants,"broad",SPEC)
    strict,strict_qa=common_original_study_blocks(meta,plants,"strict",SPEC)
    counts=(len(broad),len({x[0] for x in broad}),
            sum(len(b["rail"]) for b in broad.values()),
            sum(len(b["road"]) for b in broad.values()))
    if counts!=(44,36,229,462):
        raise AssertionError("v343 source broad study-block frame drift: "+str(counts))
    sourceplots={pid for k in (strict,broad) for group in k.values()
                 for arm in ("rail","road") for pid in group[arm]}
    plots={}
    for pid in sourceplots:
        plants_here=plants[pid]
        plots[pid]={
          "source_total_taxa":len(plants_here),
          "exact_binomial_taxa":{n for name,_,_ in plants_here
                                if (n:=normalized_binomial(name)) is not None}
        }
    traits={name:(float(v["position"]),v["systems"]) for name,v in eive.items()}
    provenance={"genuine_original_czechveg_headers":h,
         "genuine_original_czechveg_species":s,
         "genuine_original_EIVE":eiveqa,
         "full_source_plots":SOURCE_HEADERS,
         "full_original_species_records":SOURCE_ROWS,
         "source_original_study_metadata_verified":bool(qa_studies),
         "source_original_header_quality":qa["counts"],
         "source_original_time_quality":qa_time,
         "broad_original_before_traits":broad_qa,
         "HABITAT_only_original_before_traits":strict_qa,
         "original_raw_taxa_and_GPS_NOT_republished":True}
    return broad,strict,plots,traits,provenance


def subset_by_era(blocks,era):
    result={}
    for key,arms in blocks.items():
        year=key[2]  # source 10-year census bin, not actual single date
        if year is None:continue
        if era=="before_2000":keep=year<2000
        elif era=="year_2000_onward":keep=year>=2000
        elif era=="year_2010_onward":keep=year>=2010
        else:raise ValueError("Unknown frozen era")
        if keep:result[key]=arms
    return result


def safe_effect(frame,score,traits,removal):
    status,papers=report_stats(frame,score,trait_data_to_positions(traits),removal)
    if papers and len(papers)>=5:
        status["publication_bootstrap_and_leave_one_source_paper"]=quantify(papers)
    else:
        status["inference_gate"]="fewer_than_5_original_papers_or_no_eligible_plots"
    status["source_publications_are_not_independent_field_randomizations"]=True
    return status


def analyze(zenodo,headers,species,eive):
    broad,strict,plots,traits,sourceqa=source_context(zenodo,headers,species,eive)
    out={"version":"v344-CzechVeg-original-HABITAT-only-and-survey-era-contrast",
         "source_QA":sourceqa,
         "source_HABITAT_label_is_only_publisher_text_not_confirmed_ballast":True,
         "same_source_field_study_context_not_true_rail_natural_field_pair":True,
         "panels":[]}
    for nsystems,coverage in PANELS:
        sf,ss,sqa=frozen_frame(strict,plots,traits,nsystems,coverage)
        bf,bs,bqa=frozen_frame(broad,plots,traits,nsystems,coverage)
        panel={"EIVE_source_systems":nsystems,"minimum_source_trait_coverage":coverage,
            "strict_HABITAT_only_frame":sqa,
            "broad_HABITAT_plus_LOCALITY_frame":bqa,
            "strict_HABITAT_only_outcomes":{},
            "broad_outcomes_replay_only":{},
            "broad_matched_original_study_era":{}}
        for name,exclusion in REMOVAL:
            panel["strict_HABITAT_only_outcomes"][name]=safe_effect(sf,ss,traits,exclusion)
            panel["broad_outcomes_replay_only"][name]=safe_effect(bf,bs,traits,exclusion)
        for era in ERA:
            frozen=subset_by_era(bf,era)
            panel["broad_matched_original_study_era"][era]={
                "source_matched_original_blocks":len(frozen),
                "source_papers":len({k[0] for k in frozen}),
                "outcome":safe_effect(frozen,bs,traits,frozenset())
            }
        out["panels"].append(panel)
    out["limitations"]=[
        "This is a historical RAIL-text versus ROAD-text vegetation comparison in one country; label quality is original HABITAT fields only",
        "Survey era strata are temporal selection diagnoses, not independent climate regions",
        "Source publication clusters do not imply field treatment randomization; site locations and management are unmatched",
        "No ballast temperature, train frequency, propagule flux, phenology, functional adaptation measured",
        "C4 species names are exploratory trait labels from previous v339 not independently physiological confirmed here",
        "European EIVE temperature is ordinal plant ecology indicator and not measured physical temperature",
        "Original full species x plot rows, names of source locations and coordinates are not published",
    ]
    return out


def main():
    p=argparse.ArgumentParser()
    for name in ("zenodo-record","headers","species","eive","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    r=analyze(a.zenodo_record,a.headers,a.species,a.eive)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(r,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

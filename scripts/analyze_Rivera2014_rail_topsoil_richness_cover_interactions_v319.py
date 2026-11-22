#!/usr/bin/env python3
"""Reconstruct source-author treatment×year mean changes from true PLOS Table2.

This is original published aggregation, NOT individual-plot randomization
variance or independent causal treatment contrast p-values. All cells are
re-read from authentic publisher JATS-extracted source Table2.
"""
from __future__ import annotations
import argparse,json,math,re
from pathlib import Path
from audit_Rivera2014_Spain_rail_topsoil_original_JATS_v319 import SOURCE_DOI
MEASUREMENTS=("Plant cover","Richness","Equitability")
TREATMENTS=("Control","T10","T30")
YEARS=(2010,2011)
RECORDS=6
NUMBER=re.compile(r"^\s*(?P<mean>\d+(?:\.\d+)?)\s*±\s*(?P<se>\d+(?:\.\d+)?)\s*(?P<letter>[a-z]*)\s*$")

def value(text):
    m=NUMBER.fullmatch(str(text))
    if m is None:raise ValueError("Original author mean±SE has unresolved source formatting: "+repr(text))
    return {"mean":float(m.group("mean")),"SE_from_original_author":float(m.group("se")),
            "original_published_significance_letter":m.group("letter") or None}

def parse_original_table2(source):
    if source.get("author_paper_DOI")!=SOURCE_DOI or not source.get("original_source_doi_verified"):
        raise ValueError("Original Rivera 2014 JATS not authenticated")
    tables=[t for t in source["original_source_tables"] if t["label"].lower().replace(" ","")=="table2"]
    if len(tables)!=1 or "plant cover" not in tables[0]["source_author_caption"].lower():
        raise ValueError("Not one original author Table2 plant response source")
    rows=tables[0]["rows"]
    if not rows or rows[0]!=["Period","Treatment","Plant cover","Richness","Equitability"]:
        raise ValueError("Source exact treatment year and cover/richness metadata changed")
    out={}
    current=None
    for row in rows[1:]:
        if len(row)!=5:raise ValueError("Original treatment row shape changed")
        if row[0]:
            if row[0] not in ("2010","2011"):
                raise ValueError("Source period not original 2010/11")
            current=int(row[0])
        treat=row[1]
        if current is None or treat not in TREATMENTS:
            raise ValueError("Source treatment label could not be mapped")
        key=(current,treat)
        if key in out:raise ValueError("Repeated original randomized source group")
        out[key]={metric:value(row[n+2]) for n,metric in enumerate(MEASUREMENTS)}
    if len(out)!=RECORDS or set(out)!={(y,t) for y in YEARS for t in TREATMENTS}:
        raise ValueError("Incomplete true randomized paper Table2 2years×3 treatment groups")
    if not (out[2010,"Control"]["Plant cover"]["mean"]>out[2010,"T10"]["Plant cover"]["mean"] and
            out[2011,"T30"]["Richness"]["mean"]>out[2011,"Control"]["Richness"]["mean"]):
        raise ValueError("Independent original Table2 qualitative author results not reproduced")
    return out

def summarize(source):
    d=parse_original_table2(source)
    vals=[]
    for y in YEARS:
        for treat in TREATMENTS:
            obj={"year":y,"original_treatment":treat}
            obj.update(d[y,treat])
            vals.append(obj)
    change={}
    for metric in MEASUREMENTS:
        c=d[2011,"Control"][metric]["mean"]-d[2010,"Control"][metric]["mean"]
        effects={"original_control_group_2010_2011_change":float(c)}
        for treat in ("T10","T30"):
            pre=d[2010,treat][metric]["mean"]
            post=d[2011,treat][metric]["mean"]
            actual=float(post-pre)
            effects[treat]={
                "original_treatment_2010_2011_change":actual,
                "original_2010_treatment_minus_control_mean":float(pre-d[2010,"Control"][metric]["mean"]),
                "original_2011_treatment_minus_control_mean":float(post-d[2011,"Control"][metric]["mean"]),
                "descriptive_did_mean_change_minus_control_change":float(actual-c),
                "descriptive_difference_SE_p_unidentifiable_without_original_plot_level_repeated_pair_covariance":True,
            }
        effects["post_t30_minus_post_t10_mean"]=float(d[2011,"T30"][metric]["mean"]-d[2011,"T10"][metric]["mean"])
        change[metric]=effects
    expected={
      "Plant cover":{ "T10":9.3,"T30":23.2},
      "Richness":{"T10":3.5,"T30":5.9}
    }
    for key,values in expected.items():
        for treat,value0 in values.items():
            if abs(change[key][treat]["descriptive_did_mean_change_minus_control_change"]-value0)>1e-9:
                raise ValueError("Source published mean arithmetic unexpectedly changed")
    # The design is a 5 randomized complete-block sample, not 2x3 independent experiments.
    return {
        "version":"v319",
        "original_2014_PLOS_author_source_SHA256":source["true_original_article_sha256"],
        "source_original_DOI":SOURCE_DOI,
        "original_publication_context":"2014 PLOS original Railway near Madrid Spain; 5 randomized blocks, 3 treatments per block and repeated 2010/2011",
        "n_original_randomization_blocks":5,
        "n_original_treatment_plots":15,
        "plot_area_m2":54,
        "n_original_independent_second_year_studies":0,
        "source_original_table2_treatment_group_year_mean_and_SE":vals,
        "outcomes_descriptive_changes_not_independent_plot_effects":change,
        "source_additional_soil_mean_observations_not_randomized_mediation":{
            "original_substrate_organic_matter_g_kg":1.2,
            "topsoil_organic_matter_g_kg":10.4,
            "original_substrate_total_N_g_kg":0.2,
            "topsoil_total_N_g_kg":0.4,
            "source_origin":"Original PLOS Table1, unmanipulated soil material characterization, not independent plot-level mediator samples"
        },
        "published_author_interaction_results_REPORTED_not_recalculated":{
            "plant_cover_year_by_treatment_F":9.36,
            "plant_cover_year_by_treatment_P":"<0.001",
            "species_richness_year_by_treatment_F":10.46,
            "species_richness_year_by_treatment_P":"<0.001",
            "equitability_year_by_treatment_F":3.6,
            "equitability_year_by_treatment_P":0.03
        },
        "no_source_randomization_p_values_or_new_95_CI_reestimated":True,
        "original_treatment_block_repeated_year_covariances_not_reported":True,
        "does_not_measure_train_traffic_seed_transport_or_thermal_niche":True,
        "new_rail_vs_adjacent_nonrail_source_field_studies_added":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    source=json.loads(a.source.read_text())
    j=summarize(source)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(j,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps({
       "source_original_DOI":j["source_original_DOI"],
       "source_15_original_plots":j["n_original_treatment_plots"],
       "species_richness_post_vs_control":{t:j["outcomes_descriptive_changes_not_independent_plot_effects"]["Richness"][t]["original_2011_treatment_minus_control_mean"] for t in ("T10","T30")},
       "post_minus_pre_control_adjusted_mean_changes":{
           var:{t:j["outcomes_descriptive_changes_not_independent_plot_effects"][var][t]["descriptive_did_mean_change_minus_control_change"] for t in ("T10","T30")}
           for var in ("Plant cover","Richness")},
       "no_new_CI_or_p_values":j["no_source_randomization_p_values_or_new_95_CI_reestimated"]
    },indent=2,ensure_ascii=False))
if __name__=="__main__":main()

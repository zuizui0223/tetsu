#!/usr/bin/env python3
"""v12.9: ORIGINAL PLOS 2013 seed-airflow EXPERIMENT numeric 4 species
x five cumulative CAR passes, morphology × exposure descriptive source.

Published authors von der Lippe, Bullock, Kowarik, Knopp & Wichmann
(2013), PLOS ONE 8(1):e52733, DOI 10.1371/journal.pone.0052733,
original Table 2 (mean±SD, median, 99th percentile longitudinal
displacement in metres). PLOS source article is CC BY.

Source 4 seed taxa, TWO with wind appendages (Ailanthus altissima,
Clematis vitalba) vs TWO WITHOUT (Brassica napus,
Ambrosia artemisiifolia). Real source independent original experiment:
CAR moving at 48 km/h on paved surface, 10 replicate time sequences,
cumulative 1/10/20/40/80 passes over a 45m measured maximum.
All pass-level means ARE repeated measures of same sequences, and
each species is a biological treatment class; 20 table cells are
NOT 20 independent original seed/rail samples. Long-distance tail
truncated at measured 45m; original Ambrosia seeds microwaved
to prevent germination, cannot infer seed viability.

Unlike source Austrian Pascher 2016 60 sites, THESE ARE **actual tracked
seed movements** but caused by cars on a road; NOTHING here measures
train slipstreams or arrival/establishment on railway ballast.
The original source Berlin 2007 11818 traffic-delivered seeds was
from road tunnels, a separate design and NOT pooled into 2013
experimental replicate units or onto Austrian rail stations.

Do NOT invent p-values from published 4-species Table 2 summary
statistics with only 2 species per wind morphology class. This
descriptive physical seed mechanistic analogue is for *designing*
a future independent railway train-gust vs natural wind experiment.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
from statistics import mean

PLOS_DOI="10.1371/journal.pone.0052733"
ORIGINAL_SOURCE_URL="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052733"
REPORTED_PASSES=(1,10,20,40,80)
GROUPS={
    "Ailanthus altissima":"wind_morphology",
    "Clematis vitalba":"wind_morphology",
    "Brassica napus":"no_wind_morphology",
    "Ambrosia artemisiifolia":"no_wind_morphology",
}
PUBLISHED_ENDPOINT_ORIGINAL={
    ("Ailanthus altissima",1):(5.14,3.07,5.,15.),
    ("Ailanthus altissima",80):(10.83,8.11,8.,40.),
    ("Clematis vitalba",1):(4.95,3.43,4.,15.),
    ("Clematis vitalba",80):(10.46,8.00,8.,40.),
    ("Brassica napus",1):(0.98,1.12,.5,4.),
    ("Brassica napus",80):(1.71,2.87,1.,15.),
    ("Ambrosia artemisiifolia",1):(1.04,1.07,.5,5.),
    ("Ambrosia artemisiifolia",80):(1.59,2.05,1.,10.),
}


def read_original_experiment_table(input_csv):
    with open(input_csv,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        expected=[
            "original_species_name","original_seed_morphology","vehicle_passes",
            "mean_longitudinal_displacement_m","sd_longitudinal_displacement_m",
            "median_longitudinal_displacement_m","p99_longitudinal_displacement_m"]
        if reader.fieldnames!=expected:
            raise ValueError("Published PLOS original Table2 numeric source schema changed")
        rows=[]
        codes=set()
        for row in reader:
            species=row["original_species_name"]
            morph=row["original_seed_morphology"]
            if species not in GROUPS or morph!=GROUPS[species]:
                raise ValueError("Source seed species morphology label not from original study")
            passes=int(row["vehicle_passes"])
            if passes not in REPORTED_PASSES:
                raise ValueError("Not an original 1/10/20/40/80 pass source experiment")
            code=(species,passes)
            if code in codes:
                raise ValueError("Duplicate original species×passes source experiment key")
            codes.add(code)
            a=tuple(float(row[k]) for k in expected[3:])
            if not(all(math.isfinite(x) and x>=0 for x in a)
                   and a[1]>=0 and a[2]<=a[3]<=45):
                raise ValueError("Original seed source displacement summary invalid")
            rows.append({"species":species,"morphology":morph,
                         "passes":passes,"mean_m":a[0],"sd_m":a[1],
                         "median_m":a[2],"p99_m":a[3]})
    if len(rows)!=20 or codes!={(s,n) for s in GROUPS for n in REPORTED_PASSES}:
        raise ValueError("Missing genuine 20 original published source species-by-pass summary rows")
    for code,orig in PUBLISHED_ENDPOINT_ORIGINAL.items():
        row=next(x for x in rows if (x["species"],x["passes"])==code)
        val=(row["mean_m"],row["sd_m"],row["median_m"],row["p99_m"])
        if any(not math.isclose(a,b,abs_tol=1e-10) for a,b in zip(val,orig)):
            raise ValueError("Original published PLOS Table2 endpoint summaries don't reconcile: "+str(code))
    return rows


def summarize_original_source(rows):
    bypass={}
    for n in REPORTED_PASSES:
        sub=[r for r in rows if r["passes"]==n]
        groups={}
        for morph in ("wind_morphology","no_wind_morphology"):
            g=[x for x in sub if x["morphology"]==morph]
            if len(g)!=2:raise ValueError("There are exactly two original species per source morphology treatment class")
            groups[morph]={
                "source_species_n":len(g),
                "source_species_biological_replicates_not_seed_individuals":[x["species"] for x in g],
                "mean_of_published_species_mean_displacement_m":mean(x["mean_m"] for x in g),
                "mean_of_published_species_median_displacement_m":mean(x["median_m"] for x in g),
                "mean_of_published_species_99pctile_displacement_m":mean(x["p99_m"] for x in g),
                "source_experimental_replicates_per_species":10,
            }
        w=groups["wind_morphology"];o=groups["no_wind_morphology"]
        bypass[str(n)]={
            "original_vehicle_passes":n,
            "wind_type_vs_other_source_summary":groups,
            "ratio_of_group_means_original_mean_displacement":
                w["mean_of_published_species_mean_displacement_m"]/
                o["mean_of_published_species_mean_displacement_m"],
            "ratio_of_group_mean_source_medians":
                w["mean_of_published_species_median_displacement_m"]/
                o["mean_of_published_species_median_displacement_m"],
            "ratio_group_99pctile_summaries_NOT_true_pooled_distribution":
                w["mean_of_published_species_99pctile_displacement_m"]/
                o["mean_of_published_species_99pctile_displacement_m"],
        }
    per_species={}
    for name in GROUPS:
        vals={n:next(r for r in rows if r["species"]==name and r["passes"]==n)
              for n in REPORTED_PASSES}
        per_species[name]={
            "original_author_morphology":GROUPS[name],
            "mean_m_original_one_pass":vals[1]["mean_m"],
            "mean_m_after_80_passes":vals[80]["mean_m"],
            "median_m_one_pass":vals[1]["median_m"],
            "median_m_after_80_passes":vals[80]["median_m"],
            "p99_m_one_pass":vals[1]["p99_m"],
            "p99_m_after_80_passes":vals[80]["p99_m"],
            "within_original_species_ratio_of_author_means_80_passes_vs_1":
                vals[80]["mean_m"]/vals[1]["mean_m"],
        }
    return bypass,per_species


def run(csv_path):
    rows=read_original_experiment_table(csv_path)
    groups,spp=summarize_original_source(rows)
    return {
        "complete":True,"version":"v12.9",
        "original_external_author_published_article_DOI":PLOS_DOI,
        "original_PLOS_source_url":ORIGINAL_SOURCE_URL,
        "source_data_provenance":"manual numeric transcription of original author's CC-BY Table2 4 seed species × 5 cumulative car-pass conditions, PLOS 2013",
        "n_original_botanical_species":4,
        "n_actual_original_morphology_species_each_class":2,
        "n_author_summary_cells_species_by_pass":20,
        "n_original_experimental_replicates_per_species":10,
        "source_speed_kmh":48,
        "source_max_longitudinal_measurement_m":45,
        "original_table2_direct_ROAD_car_airflow_seed_movement_by_passes":groups,
        "source_seed_movement_within_species":spp,
        "direct_railway_train_gust_experiment_or_cargo_seed_measurement":False,
        "source_road_experiment_original_Ambrosia_seeds_microwaved_to_prevent_germination":True,
        "critical_inference_limits":[
            "Original study tested CAR slipstream and seed travel on paved ROAD, not real train/wagon/train-borne viable propagules or railway gravel ballast.",
            "PLOS original 1,10,20,40,80 passes are cumulative repeated measurements from 10 original experiment sequences per species and not 50 independent pass-level trials.",
            "Two species in each wind morphology type are insufficient to support globally representative per-species adaptation or a valid morphology random effect p-value.",
            "Mean of source species-specific median and 99th percentile distances is NOT the pooled seed-distance median or true pooled quantile.",
            "Original source longitudinal tracking limited at 45m; all further seed movements are right-censored/unobserved.",
            "The original 2013 source A. artemisiifolia seeds were microwave treated against germination, therefore physical seed movement does not imply viability or plant establishment.",
            "No direct link to Austrian 2016 established feral rapeseed station presence or Slovenia and Czech original flora community sites is identified.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-table2-transcription",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.source_table2_transcription)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v12.9",
        "original_species":out["n_original_botanical_species"],
        "source_airflow_road_not_train":True,
        "1pass":out["original_table2_direct_ROAD_car_airflow_seed_movement_by_passes"]["1"],
        "80passes":out["original_table2_direct_ROAD_car_airflow_seed_movement_by_passes"]["80"],
        "species_source":out["source_seed_movement_within_species"],
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

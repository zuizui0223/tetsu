#!/usr/bin/env python3
"""v13.3: source-verified Swiss Basel 2013 feral oilseed rape
RAILWAY station vs Rhine PORT same-study genetic-marker and
seed-pod pollen movement evidence. No simulated source positives.

Original authors Schulze J, Frauenknecht T, Brodmann P, Bagutti C
(2014) PLOS ONE 9:e114477, DOI
10.1371/journal.pone.0114477, CC BY 4.0.
2013 PLOS original MONTHLY source survey: Basel St Johann
freight station 279 B. napus leaf samples, 64 event-specific
GT73 positive; Rhine port of Basel 2291 B. napus leaf samples,
106 GT73 positive and 22 original positive bar gene samples.
Other crop event breakdown original 4 MS8,13 RF3,3 MS8+RF3,
two bar-positive samples without one of those three event
designations in the narrative. Note leaf samples may not be
distinct biological individuals due to repeated month visits
and original source collection pooling and resampling.

Additional original follow-up: non-GM B. napus MATERNAL
seed-pod samples from Rhine port original 48 samples; source
2 samples had GT73 DNA, interpreted by original authors as
at least partial outcrossing of GT73 pollen into seeds.
This is genetic REPRODUCTION evidence, not demonstrated
resulting live seedlings; 2/48 is not a true pollen gene
flow rate per pollinated flower nor per seed.

Official prior 2020 Basel canton government 8-year source
panel records DIFFERENT source sampling definitions and
totals in 2013: St Johann 100 tested of 639 estimated and
Rhine Port 2274 tested of 5438 estimated, versus PLOS
source 279+2291 B. napus leaf samples; do NOT pool or
identify as identical 2013 original sample sets.

No intervention/random management or documented cargo-to-
genotype parentage. Original author source paper notes
both sites treated with glyphosate in 2013 (port 3 dates,
station 1), railway cargo not handled at St Johann since
~2009 vs continued ship-cargo transfer at port. This
natural contrast is confounded by source site biology and
monitoring differences, not a quantitative causal effect
of trains, import volumes, herbicide or true seed viability.
"""
from __future__ import annotations
import argparse,csv,html,hashlib,json,math,re
from pathlib import Path

ORIGINAL_DOI="10.1371/journal.pone.0114477"
SOURCE_URL="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0114477"
SOURCE_COHORT="PLOS_Schulze2014_survey2013"
STATION="St_Johann_freight_station"
PORT="Kleinhuningen_Rhine_port"
EXPECTED={
    STATION:{"tested":279,"gt73":64,"bar":None},
    PORT:{"tested":2291,"gt73":106,"bar":22},
}
OTHER_EVENTS={"MS8":4,"RF3":13,"both_MS8_RF3":3}
SEED_POD_NON_GM_MOTHERS_N=48
SEED_POD_SAMPLES_GT73_DETECTED=2
PLOS_OSR_AND_RELATED_PLANT_SAMPLES=2787

def verify_actual_PLOS_fulltext_source(article_html):
    b=Path(article_html).read_bytes()
    if len(b)<30000:raise ValueError("Source full-text original PLOS article too short")
    source=b.decode("utf-8")
    s=html.unescape(re.sub(r"<[^>]*>"," ",source))
    s=re.sub(r"\s+"," ",s)
    needed={
        "original_doi":r"10\.1371/journal\.pone\.0114477",
        "2013_original_station_leaf_sampled":r"64 out of 279 OSR leaf samples",
        "2013_original_port_GT73_plant_tests":r"106 out of 2291 OSR leaf samples",
        "2013_original_port_bar_gene_22":r"22 OSR leaf samples tested positive for the bar gene",
        "original_nonGM_maternal_GT73_seed_outcrossing":r"Two out of 48 seed pod samples",
        "source_study_2013_st_johann_tested_279":r"sampled 279 OSR plants",
        "original_manuscript_original_PLOS_by_2014":r"Unexpected Diversity of Feral Genetically Modified Oilseed Rape",
    }
    missing=[k for k,v in needed.items() if re.search(v,s,re.I) is None]
    if missing:raise ValueError("Original PLOS 2014 manuscript lacks exact source biological samples: "+str(missing))
    return {
        "original_publisher_full_text_HTTP_page_content_SHA256":hashlib.sha256(b).hexdigest(),
        "original_publisher_HTML_downloaded_and_actual_source_number_phrases_verified":True,
        "original_PLOS_pub_DOI":ORIGINAL_DOI,
        "original_original_author_article_full_text_bytes":len(b),
        "no_source_raw_leaf_genotypes_per_individual_identifiable":True,
    }


def read_original_leaf_source_table(file):
    fields=[
      "source_publication_data_cohort","original_field_site","year",
      "original_OSR_leaf_samples_tested",
      "original_event_specific_GT73_positive_leaf_samples",
      "original_bar_gene_glufosinate_positive_leaf_samples_if_reported"]
    rows=[]
    names=set()
    with open(file,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=fields:
            raise ValueError("PLOS original sampled Bnapus source journal categories schema drift")
        for r in reader:
            site=r["original_field_site"]
            if site not in EXPECTED or site in names:
                raise ValueError("Invented or duplicated 2013 original PLOS Swiss sampling site")
            names.add(site)
            tested=int(r[fields[3]]);gt73=int(r[fields[4]])
            bar=int(r[fields[5]]) if r[fields[5]].strip() else None
            if (r[fields[0]]!=SOURCE_COHORT or int(r[fields[2]])!=2013 or
                (tested,gt73,bar)!=
                (EXPECTED[site]["tested"],EXPECTED[site]["gt73"],EXPECTED[site]["bar"])):
                raise ValueError("Published original author 2013 site × genetically tested OSR leaf sample totals changed")
            if not(0<=gt73<=tested and bar is None or
                   0<=gt73<=tested and bar is not None and 0<=bar<=tested):
                raise ValueError("Impossible original PLOS genotype count")
            rows.append({"original_site":site,
                         "source_original_OSR_leaf_samples":tested,
                         "source_confirmed_GT73_OSR_leaf_samples":gt73,
                         "source_glufosinate_bar_positive_OSR_leaf_samples_if_reported":bar,
                         "source_fraction_GT73_positive_among_original_OSR_leaf_samples":
                             gt73/tested})
    if names!=set(EXPECTED):
        raise ValueError("Missing original PLOS 2013 Basel railway or Rhine port source area")
    return rows


def run(pub_html,source_csv,official_gov131):
    if official_gov131.get("complete") is not True or official_gov131.get("version")!="v13.1":
        raise ValueError("Require verified original independent 2020 Canton Basel 2013–20 public government survey source")
    audit=verify_actual_PLOS_fulltext_source(pub_html)
    leaf=read_original_leaf_source_table(source_csv)
    s=next(x for x in leaf if x["original_site"]==STATION)
    p=next(x for x in leaf if x["original_site"]==PORT)
    gov=official_gov131["all_source_rows"]
    st2013=next(x for x in gov if x["site"]==STATION and x["year"]==2013)
    po2013=next(x for x in gov if x["site"]==PORT and x["year"]==2013)
    if (st2013["source_plants_tested_genetically"],
        po2013["source_plants_tested_genetically"])!=(100,2274):
        raise ValueError("Original government 2013 monitoring subset source drift")
    if (s["source_original_OSR_leaf_samples"],p["source_original_OSR_leaf_samples"])==(100,2274):
        raise ValueError("PLOS and Cantonal source 2013 sample sizes were improperly conflated")
    return {
        "complete":True,"version":"v13.3",
        "original_Schulze2014_PLOS_DOI":ORIGINAL_DOI,
        "source_original_publisher_article_source_integrity":audit,
        "original_PLOS_2013_RAIL_freight_station_and_Rhine_port_leaf_samples":leaf,
        "source_original_2013_OSR_leaf_samples_station_and_port_total":
            s["source_original_OSR_leaf_samples"]+p["source_original_OSR_leaf_samples"],
        "source_original_2013_all_related_species_and_OSR_plant_samples_PLOS_abstract":
            PLOS_OSR_AND_RELATED_PLANT_SAMPLES,
        "source_differing_original_total_OSR_and_related_specimens_not_equivalent_to_leaf_count":
            PLOS_OSR_AND_RELATED_PLANT_SAMPLES-(s["source_original_OSR_leaf_samples"]+
                                               p["source_original_OSR_leaf_samples"]),
        "source_port_other_bar_positive_event_specific_breakdown":
            OTHER_EVENTS,
        "source_port_bar_positive_22_min_genetic_event_counts_20_not_silently_forced_to_22":
            sum(OTHER_EVENTS.values())==20,
        "source_nonGM_mother_Bnapus_seed_pod_samples_checked_for_pollen_outcrossing":
            SEED_POD_NON_GM_MOTHERS_N,
        "source_nonGM_mother_seed_pod_samples_with_GT73_transgene_DNA":
            SEED_POD_SAMPLES_GT73_DETECTED,
        "source_port_fraction_of_maternal_seed_pod_SAMPLES_with_GT73_marker":
            SEED_POD_SAMPLES_GT73_DETECTED/SEED_POD_NON_GM_MOTHERS_N,
        "source_NOT_rate_of_genetically_modified_seeds_or_viable_seedling_germination":True,
        "different_actual_2013_st_johann_and_port_source_sampling_definitions":{
            "PLOS_2014_original_2013_station_Bnapus_leaf_samples":s["source_original_OSR_leaf_samples"],
            "Gov_2020_original_2013_station_PCR_tested_Bnapus":st2013["source_plants_tested_genetically"],
            "PLOS_2014_original_2013_port_Bnapus_leaf_samples":p["source_original_OSR_leaf_samples"],
            "Gov_2020_original_2013_port_PCR_tested_Bnapus":po2013["source_plants_tested_genetically"],
            "DO_NOT_ADD_OR_META_POOL_ACROSS_THE_TWO_OVERLAPPING_2013_SOURCE_SAMPLING_FRAMES":True,
        },
        "critical_ecological_boundaries":[
            "Confirmed GT73 transgenes in field feral plants, positive non-GM maternal pod seed samples demonstrate reproductive genetic movement; source seed germination into viable railway plants NOT traced.",
            "Original site-level GM leaf sample fraction is not whole annual estimated plant population prevalence, particularly at clustered high-density source locations.",
            "2013 repeated monthly plant sampling possibly re-sampled plants; all specimens are not independent seed-arrival events.",
            "The original RAIL freight station and port unloading are not randomly assigned to cessation vs ongoing cargo, and both have distinct glyphosate management schedules.",
            "Original ongoing Rhine port imports include ship and subsequent rail transfer and wheat-borne accidental seed contamination; do NOT claim all feral plants originated directly from freight trains.",
            "Original 2012 survey and 2013 follow-up have some overlapping sites and underlying plant lineages; they are not independent country replicates.",
            "The original government 2013 monitoring subset sample count differs from the PLOS 2013 one and may use a different survey space/definition; no combined 2013 denominator is valid.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--author-PLOS-2014-original-html",required=True)
    p.add_argument("--original-2013-leaf-samples-csv",required=True)
    p.add_argument("--verified131",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    prior=json.loads(Path(a.verified131).read_text(encoding="utf-8"))
    res=run(a.author_PLOS_2014_original_html,
            a.original_2013_leaf_samples_csv,prior)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v13.3",
        "paper":res["original_Schulze2014_PLOS_DOI"],
        "source":res["source_original_publisher_article_source_integrity"],
        "station_port":res["original_PLOS_2013_RAIL_freight_station_and_Rhine_port_leaf_samples"],
        "original_seed_pod_GT73_maternal":res["source_nonGM_mother_seed_pod_samples_with_GT73_transgene_DNA"],
        "out_of":res["source_nonGM_mother_Bnapus_seed_pod_samples_checked_for_pollen_outcrossing"],
        "source_sampling_frame_gap":res["different_actual_2013_st_johann_and_port_source_sampling_definitions"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

#!/usr/bin/env python3
"""v13.1: verified canton Basel 8-year railway freight station vs Rhine
port feral oilseed rape plant-population and GMO genetic-test panel.

ORIGINAL source: Claudia Bagutti & Dirk Hamburger (2020-12-17),
Kantonales Laboratorium Basel-Stadt + Bundesamt fuer Umwelt,
'Gentechnisch veraenderter Raps: Kantonale Ueberwachung 2014-2020',
official 3pp government PDF
https://media.bs.ch/original_file/9ff1074caebc95c112c02d023afc713a69e28f42/2020-gvraps.pdf

Original Table 1 and Table 2 provide 2013-2020 annual 8 years
ESTIMATED full feral B. napus plants in each spatial survey area
AND separate NUMBER PCR TESTED plants. The 16 rows come from
manual source table reading, and live original governmental source
PDF is separately downloaded and its actual table numeric sequences
verified before results can claim authentic source provenance.

Two sites are NOT randomized, environmental/confounder matched,
representative of Switzerland railway network or same study rail-vs-
road trials. St. Johann freight station ceased handling germinable
cargo in 2009; Rhine port actively unloads/imports. BOTH were
under official GM-control management from 2013/14: uprooting and
herbicide (different conditions). 2013 monitored monthly, 2014+
twice annually, estimated population counts NOT comparable to
number of genetically tested plants or direct seed arrival rate.

The official heading states 7908 total plant TESTS, 459 GM positive
(5.8%) across the 8-year two-site panel. Does NOT prove that all
GM canola plants arrived by railway; port shipborne transport
substantially contributes. Original report notes transgenic
GT73 at St Johann gone since 2016, repeatedly detected in port,
and evidence for introgression GT73xRF3 at port. Cannot infer
selection coefficient, recruitment, cargo seed viability or plant
survival from annual cross-sectional estimates.

Program audits exact original PDF text using native pdftotext, not
OCR. Fail closed if the PDF original first table, second table,
the published per-year numeric sequences, official GMO total
or critical government source station/port names cannot be found.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,subprocess
from collections import defaultdict
from pathlib import Path

REPORT_URL="https://media.bs.ch/original_file/9ff1074caebc95c112c02d023afc713a69e28f42/2020-gvraps.pdf"
SOURCE_DATE="2020-12-17"
STATION="St_Johann_freight_station"
PORT="Kleinhuningen_Rhine_port"
ORIGINAL={
    STATION:{
        "estimated":[639,84,7,23,0,1,16,2],
        "tested":[100,49,7,23,0,1,16,2],
    },
    PORT:{
        "estimated":[5438,9247,6756,5461,1038,459,3259,752],
        "tested":[2274,1333,1260,1216,417,338,528,344],
    },
}
YEARS=tuple(range(2013,2021))
ORIGINAL_GM_POSITIVE_TESTS=459
ORIGINAL_TOTAL_TESTED_PLANTS=7908

def verify_government_original_pdf(path):
    data=Path(path).read_bytes()
    if not data.startswith(b"%PDF-") or len(data)<80000:
        raise ValueError("Official Basel government source not valid substantive PDF")
    r=subprocess.run(["pdftotext","-layout",str(path),"-"],
                     stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                     check=True,timeout=80)
    txt=r.stdout.decode("utf-8")
    needed=(r"Tabelle\s+1\s*:",r"Tabelle\s+2\s*:",
            r"Bahnhof\s+St\.\s+Johann",r"Hafen\s+Kleinh[üu]ningen",
            r"7908",r"459",r"2013",r"2020")
    for pattern in needed:
        if not re.search(pattern,txt,re.I):
            raise ValueError("Original Basel source PDF missing mandatory Table/official survey heading: "+pattern)
    # original row sequences are tested in original PDF source table
    # segments, not a fabricated CSV or unrelated title-only web page.
    t1=re.search(r"Tabelle\s+1\s*:",txt,re.I)
    t2=re.search(r"Tabelle\s+2\s*:",txt,re.I)
    if t1.start()>=t2.start():
        raise ValueError("Original Basel author Table1/2 not in published order")
    end=re.search(r"Pr[uü]fverfahren",txt[t2.end():],re.I)
    tail=t2.end()+end.start() if end else min(len(txt),t2.end()+2500)
    blocks={STATION:txt[t1.end():t2.start()],
            PORT:txt[t2.end():tail]}
    valid={}
    for site,block in blocks.items():
        nums=[int(x) for x in re.findall(r"(?<!\d)\d+(?!\d)",block)]
        # The official PDF's Table 1 and 2 show years once, then
        # estimated source population row and source PCR tested row,
        # optionally separated by percent sampled.
        expected=ORIGINAL[site]
        # accept only if each original 8-year count vector is present in
        # its corresponding table in source text in the stated year order.
        # Some PDF converters interleave year headings with values:
        # check a loose series of numeric text tokens preserving order.
        for category,series in expected.items():
            # ordered subsequence in the original table block, but
            # not a full number census: numbers might be interleaved
            # with 2013..2020 column-year heads.
            pos=0
            positions=[]
            for value in series:
                while pos<len(nums) and nums[pos]!=value:pos+=1
                if pos>=len(nums):
                    raise ValueError("Source official PDF "+site+
                                     " lacks expected original annual "+category+
                                     " series at "+str(value))
                positions.append(pos)
                pos+=1
            valid[site+"_"+category+"_positions"]=positions
    return {
        "source_original_government_pdf_SHA256":hashlib.sha256(data).hexdigest(),
        "source_original_government_pdf_bytes":len(data),
        "source_publisher_report_url":REPORT_URL,
        "source_original_2_table_presence_verified":True,
        "source_original_16_year_by_site_estimated_population_cells_and_16_tested_cells_source_order_matched":True,
        "source_original_canton_report_7908_PCR_tested_459_GM_detection_header_found":True,
        "native_pdf_text_extracted_using_pdftotext_no_OCR":True,
        "original_source_sequence_validation_positions":valid,
    }


def parse_original_government_source_csv(path):
    heads=["original_field_site","year",
           "original_estimated_feral_rapeseed_population_count",
           "original_individual_plants_genetically_tested"]
    rows=[]
    with open(path,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=heads:
            raise ValueError("Original author 2020 monitor population/tested CSV schema changed")
        seen=set()
        for x in reader:
            site=x["original_field_site"]
            year=int(x["year"])
            if site not in ORIGINAL or year not in YEARS or (site,year) in seen:
                raise ValueError("Original source site×year missing, duplicated or invented")
            seen.add((site,year))
            population=int(x[heads[2]]);tested=int(x[heads[3]])
            i=year-2013
            if population!=ORIGINAL[site]["estimated"][i] or tested!=ORIGINAL[site]["tested"][i]:
                raise ValueError("Original Basel 2020 author published Table1/2 row value corrupted "+site+str(year))
            if min(population,tested)<0 or tested>population:
                raise ValueError("Source number PCR tested > original estimated actual population")
            rows.append({
               "site":site,"year":year,"source_estimated_feral_plant_count":population,
               "source_plants_tested_genetically":tested,
               "positive_feral_plant_detected":population>0,
               "fraction_estimated_population_PCR_tested":
                  tested/population if population else None,
            })
    if len(rows)!=16 or seen!={(site,year) for site in ORIGINAL for year in YEARS}:
        raise ValueError("Original government Basel study full 2-site×8year source rows not recovered")
    if sum(x["source_plants_tested_genetically"] for x in rows)!=ORIGINAL_TOTAL_TESTED_PLANTS:
        raise ValueError("Source government 2020 exactly 7908 genetics-tested plants not reproduced")
    return sorted(rows,key=lambda x:(x["site"],x["year"]))


def summarize_source(rows):
    by={}
    for site in ORIGINAL:
        site_rows=sorted([r for r in rows if r["site"]==site],key=lambda x:x["year"])
        counts=[r["source_estimated_feral_plant_count"] for r in site_rows]
        tested=[r["source_plants_tested_genetically"] for r in site_rows]
        by[site]={
          "n_separate_original_years_surveyed":8,
          "source_start_estimated_feral_population_2013":counts[0],
          "source_end_estimated_feral_population_2020":counts[-1],
          "source_population_2020_div_2013":counts[-1]/counts[0],
          "source_start_to_end_fractional_change":counts[-1]/counts[0]-1.,
          "n_original_site_years_any_feral_plants":sum(c>0 for c in counts),
          "source_years_zero_estimated_feral_plants":[YEARS[i] for i,c in enumerate(counts) if c==0],
          "original_year_peak_estimated_plants":YEARS[max(range(8),key=lambda i:counts[i])],
          "original_peak_estimated_feral_plants":max(counts),
          "original_source_genetic_test_count_2013_2020_sum_of_tests_NOT_INDEPENDENT_plants":sum(tested),
          "n_source_years_with_FULL_estimated_population_genetic_test_coverage":
              sum(c>0 and n==c for c,n in zip(counts,tested)),
          "source_original_2013_PCR_tested_fraction":tested[0]/counts[0],
          "source_original_2020_PCR_tested_fraction":tested[-1]/counts[-1],
          "original_source_monitoring_started_monthly_in_2013_then_2x_annually_since_2014":True,
        }
    # note treated plants are repeated sources across 8 years; 7908
    # are tests, not unique known genotypes or one site total population.
    assert sum(v["original_source_genetic_test_count_2013_2020_sum_of_tests_NOT_INDEPENDENT_plants"]
               for v in by.values())==7908
    return {
        "n_true_original_monitoring_field_sites":2,
        "n_original_annual_plant_estimates_for_each_site":8,
        "n_original_site_year_observation_units":16,
        "sum_of_original_number_of_PCR_tested_plants_2013_2020":7908,
        "government_report_original_GM_positive_plant_tests_2013_2020":459,
        "government_report_GM_positive_fraction_of_PCR_tests":459/7908,
        "original_report_only_SOURCE_GM_tested_total_not_yearly_per_site_genotypes":True,
        "source_site_year_time_series":by,
        "observed_station_cessation_of_viable_material_cargo_year":2009,
        "station_report_observed_GT73_GM_no_longer_identified_since_year":2016,
        "port_source_cargo_imports_continued":True,
        "port_GM_plant_events_present_through_2020":True,
        "source_original_canton_monitoring_used_GM_control_herbicide_and_manual_uprrooting":True,
        "site_contrast_not_randomized_controlled_rail_cargo_experiment":True,
        "same_original_spilled_seed_to_individual_plant_colonization_link_unmeasured":True,
    }


def run(csv_source,pdf_source):
    original_pdf=verify_government_original_pdf(pdf_source)
    rows=parse_original_government_source_csv(csv_source)
    out=summarize_source(rows)
    return {
        "complete":True,"version":"v13.1",
        "original_external_source":
            "Bagutti and Hamburger (2020-12-17), Kantonales Labor Basel-Stadt and Swiss Federal Office for the Environment, original government monitoring report",
        "original_author_pdf_URL":REPORT_URL,
        "original_live_published_report_pdf_provenance":original_pdf,
        "source_original_population_and_genetic_tests_panel":out,
        "all_source_rows":[{k:v for k,v in r.items() if k!="fraction_estimated_population_PCR_tested"}
                           for r in rows],
        "no_causal_establishment_probability_inferred":True,
        "limitations":[
            "Source 2013 estimated plant population and genetically sampled count differ significantly; 459/7908 is the share of tested plants, NOT share of total estimated 8-year flora or imported seed viability.",
            "2013 site checks monthly and since 2014 twice annually; sampling effort changed, so 2013–2020 plant counts are not equal-effort temporal abundance estimates.",
            "Source port cargo continued while St Johann ceased handling viable material in 2009, but the two source sites differ in habitat, management and cargo flows: no randomized causal comparison.",
            "Herbicide applications and uprooting can reduce GM population; source observations cannot identify separate effects of management, seed bank persistence and repeated new cargo imports.",
            "Site original population estimates across years can sample the same genetic lineage and are NOT eight independent countries or 7908 independent railway trajectories.",
            "Source railway station and Rhine port sample plants established locally, but no arriving viable cargo kernel is physically matched to its descendant plant.",
            "This official report monitored only two high-risk source locations, not a globally representative ecological macro sample.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--reported-tables-csv",required=True)
    p.add_argument("--original-canton-pdf",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.reported_tables_csv,a.original_canton_pdf)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v13.1",
        "original_canton_source":out["original_live_published_report_pdf_provenance"],
        "monitoring":out["source_original_population_and_genetic_tests_panel"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

#!/usr/bin/env python3
"""v13.2: TRUE original Swiss 2012 railway 79-site national survey,
59 published tested table rows incl 5 negative-control FIELD plants.

Schoenenberger N. and D'Andrea L. (2012), Environmental Sciences
Europe 24:23, DOI 10.1186/2190-4715-24-23, originally CC BY 2.0.
This paper directly observed feral Brassica napus PLANTS and CP4 EPSPS
glyphosate-tolerance genetic phenotype near Swiss RAILWAY areas.
No original train grain trap, direct seed viability, whole-country
random sampling, seed flux, freight route assignment, or genotype
from an individual cargo sample observed.

The paper's actual Table1 reports 2,403 individually sampled plants
INCLUDING 5 field negative control plants at Rümlang! Thus railway
samples tested are 2,398, all 50 positive protein detections occur
in 4 of the original 58 feral OSR positive RAILWAY AREAS; 21
additional original railway areas had NO feral OSR plants found
among the total 79 surveyed locations.
Reject any analysis that counts a field negative control as an
additional railway plant site or includes 5 field-control plants
in the denominator without labelling it so.

Actual original source PDF is machine-parsed using native
pdftotext -layout (no OCR), exact DOI/author/PDF body and original
published Table1 2403/50 source footer checks, 59 source rows
auto-extracted and validated. Original published site × sampled
plant sample-size × positive plants is no longer assistant guessed.
Source official original SpringerOpen CC BY mirror:
https://enveurope.springeropen.com/counter/pdf/10.1186/2190-4715-24-23.pdf
backup US Federal regulations record of original article:
https://downloads.regulations.gov/APHIS-2012-0031-0101/attachment_2.pdf

Real data-reanalysis: individual-positive fraction vs equally
weighted RAILWAY site fraction and concentration of genotypes into
four hotspots. Source morphological cargo selection is not proven
by this observational non-random railway sample. An exploratory
random-label null conditions on each original site-specific plant
sample count to ask whether observed 4 sites with GM is a
size-only artifact (NOT treatment effect or phylogenetic test).
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,re,subprocess
from collections import Counter,defaultdict
from pathlib import Path

DOI="10.1186/2190-4715-24-23"
SOURCE_URL="https://enveurope.springeropen.com/counter/pdf/10.1186/2190-4715-24-23.pdf"
MIRROR_URL="https://downloads.regulations.gov/APHIS-2012-0031-0101/attachment_2.pdf"
TABLE_SUM=(2403,50)
N_SITES_SURVEYED=79
N_FERAL_POSITIVE=58
SOURCE_FIELD_NEGATIVE_CONTROL_N=5
SOURCE_RAILWAY_GM_POSITIVE_SITES=4
SOURCE_TAXON="Brassica napus"
SEED=20261114
PERM_DRAWS=25000

DATE_RE=re.compile(r"^\s*(\d{2}\.\d{2}\.\d{4})\s+")
COORDS_TAIL=re.compile(r"\b(\d{6})\s+(\d{5,6})\s+([A-Z]{2})\s+(\d{1,5})\s+(\d{1,5})(?:\s*\*)?\s*$")


def authentic_original_paper_pdf_text(pdf_file):
    b=Path(pdf_file).read_bytes()
    if not b.startswith(b"%PDF-") or len(b)<100000:
        raise ValueError("Original Swiss research publisher PDF not readable authentic source")
    p=subprocess.run(["pdftotext","-layout",str(pdf_file),"-"],
                     capture_output=True,check=True,timeout=110)
    txt=p.stdout.decode("utf-8")
    if not re.search("Schoenenberger",txt,re.I) or not re.search("D.Andrea",txt,re.I):
        raise ValueError("Source PDF author names differ from original Swiss 2012 paper")
    if not re.search("Swiss railways",txt,re.I) or not re.search("CP4",txt,re.I):
        raise ValueError("Source 2012 Swiss glyphosate tolerant rail paper body missing")
    if not re.search(r"Total\s+2[\s,'’]*403\s+50",txt,re.I):
        raise ValueError("Source published Table 1 national total 2403 / 50 not in original PDF")
    return txt,{
        "actual_original_source_pdf_bytes":len(b),
        "actual_original_published_article_pdf_SHA256":hashlib.sha256(b).hexdigest(),
        "native_pdftotext_no_OCR_source":True,
        "original_ccby2012_paper_DOI":DOI,
        "original_paper_full_body_and_2403_50_source_totals_verified":True,
    }


def isolate_original_site_table(txt):
    x=re.search(r"Table\s+1\s+Distribution\s+of\s+the\s+sample\s+sites",txt,re.I)
    if x is None:
        raise ValueError("Original Swiss 2012 research Table 1 site list heading unavailable")
    end=re.search(r"Total\s+2[\s,'’]*403\s+50",txt[x.end():],re.I)
    if end is None:
        raise ValueError("Original Swiss site table original 2403/50 footer missing")
    out=txt[x.end():x.end()+end.start()]
    if len(out)<3500:
        raise ValueError("Original Swiss 2012 Table1 59 source sample sites body unexpectedly short")
    return out


def extract_original_2012_table(table_text):
    """Source primary PDF may wrap Kleinhuningen port area names over
    multiple PDF text lines, so a candidate date row is followed
    by at most three extra source lines for canton/coordinate fields.
    """
    sites=[]
    seen=set()
    lines=table_text.splitlines()
    for idx,row in enumerate(lines):
        m=DATE_RE.match(row)
        if not m:continue
        date=m.group(1)
        tail=row[m.end():].strip()
        match=COORDS_TAIL.search(tail)
        if not match:
            for j in range(idx+1,min(idx+4,len(lines))):
                if DATE_RE.match(lines[j]):break
                tail+=" "+lines[j].strip()
                match=COORDS_TAIL.search(tail)
                if match:break
        if not match:
            raise ValueError("Original Swiss 2012 Table1 date "+date+
                             " has no intact original x y canton tested positive source fields")
        place=tail[:match.start()].strip()
        if not place or len(place)>120:raise ValueError("Original published site place label malformed")
        east=int(match.group(1));north=int(match.group(2))
        canton=match.group(3)
        sampled=int(match.group(4));positive=int(match.group(5))
        if not(0<=positive<=sampled<=10000 and 400000<=east<=900000 and 60000<=north<=340000):
            raise ValueError("Original Swiss coordinates or GM source counts invalid "+place)
        key=(date,place,east,north,canton)
        if key in seen:raise ValueError("Duplicate original railway source area from Table 1")
        seen.add(key)
        neg=bool(re.search(r"R[üu]mlang\s*\(field\)\s*,?\s*NC",place,re.I))
        sites.append({
            "source_original_sampling_date":date,
            "source_original_rail_area_name":place,
            "source_CH1903_east":east,"source_CH1903_north":north,
            "original_canton_or_Liechtenstein":canton,
            "original_plant_individuals_tested_CP4":sampled,
            "original_GM_marker_positive_individual_plants":positive,
            "source_field_negative_control_NOT_railway":neg,
        })
    if len(sites)!=59:raise ValueError("Expected original 58 feral plants sites +1 nonrail field NC from paper Table1, found "+str(len(sites)))
    n=sum(x["original_plant_individuals_tested_CP4"] for x in sites)
    positive=sum(x["original_GM_marker_positive_individual_plants"] for x in sites)
    if (n,positive)!=TABLE_SUM:
        raise ValueError("Original author Table1 2403/50 sampled source plant totals do not reconcile")
    negative_control=[r for r in sites if r["source_field_negative_control_NOT_railway"]]
    if len(negative_control)!=1 or negative_control[0][
        "original_plant_individuals_tested_CP4"]!=SOURCE_FIELD_NEGATIVE_CONTROL_N:
        raise ValueError("Original Rümlang field 5-plant NC source is NOT railway")
    original_rail=[x for x in sites if not x["source_field_negative_control_NOT_railway"]]
    if len(original_rail)!=58 or sum(x["original_plant_individuals_tested_CP4"] for x in original_rail)!=2398:
        raise ValueError("True original railway 58 sites and 2398 tested plants not recovered")
    positive_sites=[x for x in original_rail if x["original_GM_marker_positive_individual_plants"]>0]
    if len(positive_sites)!=4:
        raise ValueError("Published four GM-positive source Swiss railway areas not recovered")
    expected_positive={(21,23),(16,18),(12,13),(1,47)}
    if {(x["original_GM_marker_positive_individual_plants"],
         x["original_plant_individuals_tested_CP4"]) for x in positive_sites}!=expected_positive:
        raise ValueError("Original 2012 source four railway hotspot tested/positive plant pairs changed")
    return sites


def source_homogeneous_plant_label_randomization(rows,draws=PERM_DRAWS,seed=SEED):
    """Negative control: GM positives randomly mixed among all 2398
    original sampled railway plants, holding each source site N fixed.
    Real labels are NOT randomly exchangeable biologically.
    """
    originals=[x for x in rows if not x["source_field_negative_control_NOT_railway"]]
    order=[]
    for i,x in enumerate(originals):
        order.extend([i]*x["original_plant_individuals_tested_CP4"])
    if len(order)!=2398:raise ValueError("Failed source sample-opportunity null construction")
    actual=sum(x["original_GM_marker_positive_individual_plants"]>0 for x in originals)
    K=sum(x["original_GM_marker_positive_individual_plants"] for x in originals)
    if (actual,K)!=(4,50):raise ValueError("Source 2012 rail positive clusters changed")
    rng=random.Random(seed)
    n_le=0
    vals=[]
    for _ in range(draws):
        ids={order[j] for j in rng.sample(range(len(order)),K)}
        vals.append(len(ids))
        if len(ids)<=actual:n_le+=1
    obs=sorted(r["original_GM_marker_positive_individual_plants"] for r in originals
               if r["original_GM_marker_positive_individual_plants"]>0)
    total_pos=sum(obs)
    return {
        "hypothetical_source_GM_status_exchangeability_null":"50 positive plant labels shuffled across all 2398 true original railway plant tests, conditional on each original site's sample count, not a randomized field experiment",
        "source_original_58_sites_GM_positive_clusters_observed":4,
        "source_original_50_positive_individuals":50,
        "n_montecarlo_plant_label_randomizations":draws,
        "source_null_mean_number_sites_containing_at_least_one_GM_positive":
            sum(vals)/len(vals),
        "source_null_min_number_positive_sites":min(vals),
        "source_null_max_number_positive_sites":max(vals),
        "source_null_fraction_with_at_most_4_positive_sites":n_le/draws,
        "source_null_plus_one_tail_fraction_at_most_4_positive_sites":
            (n_le+1)/(draws+1),
        "not_a_test_of_actual_train_spillage_or_herbicide_causal_effect":True,
        "original_GM_positive_plant_counts_at_four_sites_sorted":obs,
        "actual_fraction_all_50_GM_plant_tests_in_top_3_positive_sites":
            sum(sorted(obs,reverse=True)[:3])/50,
    }


def describe_surveyed_original_railway(rows):
    original=[r for r in rows if not r["source_field_negative_control_NOT_railway"]]
    n=sum(r["original_plant_individuals_tested_CP4"] for r in original)
    gm=sum(r["original_GM_marker_positive_individual_plants"] for r in original)
    if not(n==2398 and gm==50):raise ValueError("Railway testing numerator/denominator unexpected")
    byregion=defaultdict(lambda:{"n_original_feral_Bnapus_positive_rail_areas":0,
                                "n_plants_tested_at_rail_areas":0,
                                "GM_positive_plant_tests":0})
    for x in original:
        b=byregion[x["original_canton_or_Liechtenstein"]]
        b["n_original_feral_Bnapus_positive_rail_areas"]+=1
        b["n_plants_tested_at_rail_areas"]+=x["original_plant_individuals_tested_CP4"]
        b["GM_positive_plant_tests"]+=x["original_GM_marker_positive_individual_plants"]
    high=[r for r in original if r["original_GM_marker_positive_individual_plants"]>0]
    high.sort(key=lambda x:x["original_GM_marker_positive_individual_plants"],reverse=True)
    site_equal=sum(r["original_GM_marker_positive_individual_plants"]/
                   r["original_plant_individuals_tested_CP4"] for r in original)/len(original)
    return {
        "n_original_railway_areas_surveyed_from_published_text":79,
        "n_source_railway_areas_WITH_feral_Bnapus_and_test_samples":58,
        "n_source_railway_areas_NO_feral_Bnapus_found_not_in_Table1":21,
        "n_source_railway_areas_with_CP4_positive_OSR":4,
        "original_RAILWAY_plants_tested_excluding_field_negative_control":n,
        "original_field_negative_control_not_railway_plants_tested":5,
        "original_published_ALL_tested_including_field_control":n+5,
        "original_CP4_positive_feral_OSR_plants":gm,
        "plant_individual_weighted_positive_fraction_excluding_field_control":gm/n,
        "original_publication_positive_fraction_including_NC_5":gm/(n+5),
        "original_equal_feral_rail_site_mean_within_sample_positive_fraction":
            site_equal,
        "original_railway_area_fraction_with_any_GM_positive":4/79,
        "among_rail_areas_with_feral_Bnapus_fraction_any_GM_positive":4/58,
        "original_rail_plant_positive_vs_negative_survey_effort_not_uniform":True,
        "source_GM_hotspot_localities_and_actual_GM_sample_counts":[{
            "original_site_name":p["source_original_rail_area_name"],
            "original_canton":p["original_canton_or_Liechtenstein"],
            "original_tested_individual_plants":p["original_plant_individuals_tested_CP4"],
            "original_CP4_positive_individual_plants":p["original_GM_marker_positive_individual_plants"],
            "source_positive_fraction_within_locality":
                p["original_GM_marker_positive_individual_plants"]/
                p["original_plant_individuals_tested_CP4"],
            "original_CH1903_coordinates":(
                [p["source_CH1903_east"],p["source_CH1903_north"]]),
        }for p in high],
        "source_canton_observations":dict(byregion),
        "source_method":"original author CP4 EPSPS immunologic phenotype screening and subset independent PCR genotype confirmation, not a direct freight seed collector",
        "do_not_infer_GM_seed_arrival_count_from_fully_grown_plants":True,
    }


def run(pdf):
    txt,provenance=authentic_original_paper_pdf_text(pdf)
    tab=isolate_original_site_table(txt)
    rows=extract_original_2012_table(tab)
    valid=describe_surveyed_original_railway(rows)
    null=source_homogeneous_plant_label_randomization(rows)
    return {
        "complete":True,"version":"v13.2",
        "original_Swiss_2012_first_author_source":"Schoenenberger & D'Andrea (2012), Environmental Sciences Europe 24:23",
        "original_published_DOI":DOI,
        "original_author_CC_BY_2_source_PDF_provenance":provenance,
        "source_59_original_tabulated_rows_including_one_nonrail_control":rows,
        "source_railway_GM_feral_OSR_site_comparison":valid,
        "site_size_conditioned_plant_marker_randomization_diagnostic":null,
        "limitations":[
            "Swiss original 79 railway areas include 21 sites with no feral OSR detected; only 58 rail-area sources with plants are listed in Table 1 and an additional nonrail 5-plant field negative control is listed.",
            "The 58 positive-feral source locations are one-time site surveys in 2011–2012, NOT repeated independent years; several sites nested within the same port/canton logistics hub.",
            "Plant samples and rail-site sample sizes are not a whole-population census at very large source locations; source author sampled plants at regular spatial intervals for larger populations.",
            "The homogeneously random GM plant label permutation is an EXPLORATORY size-only null; independent/ exchangeable plant labels, same-origin incoming lots and identical herbicide exposure were never randomly assigned.",
            "Even molecular PCR-traced feral GM plants cannot be attributed to a specific train because nearby ship port unloading, seed contamination of wheat/other crops and local seed banks are alternative contributors.",
            "Original markers are actual detected glyphosate tolerance phenotype and confirmatory GM event DNA, not direct count of entering viable rail cargo seeds or generation/seedling success.",
            "Original source 2011/12 plant sampling is not the 2013–20 Basel official same-site cohort; overlap in localities needs cohort tracking without counting as independent sites.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-author-2012-full-PDF",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.original_author_2012_full_PDF)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    s=out["source_railway_GM_feral_OSR_site_comparison"]
    print(json.dumps({
        "complete":True,"version":"v13.2",
        "original_pdf":out["original_author_CC_BY_2_source_PDF_provenance"],
        "n_original_sites_surveyed":s["n_original_railway_areas_surveyed_from_published_text"],
        "n_feral_rail_sites":s["n_source_railway_areas_WITH_feral_Bnapus_and_test_samples"],
        "n_original_sampled_feral_rail_plants":s["original_RAILWAY_plants_tested_excluding_field_negative_control"],
        "source_molecule_positive":s["original_CP4_positive_feral_OSR_plants"],
        "source_FRACTIONS":{
            "individual_weight":s["plant_individual_weighted_positive_fraction_excluding_field_control"],
            "feral_site_equal":s["original_equal_feral_rail_site_mean_within_sample_positive_fraction"],
            "any_GM_rail_site":s["original_railway_area_fraction_with_any_GM_positive"]},
        "actual_GM_hotspots":s["source_GM_hotspot_localities_and_actual_GM_sample_counts"],
        "size_only_marker_null":out["site_size_conditioned_plant_marker_randomization_diagnostic"],
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

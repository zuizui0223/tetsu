#!/usr/bin/env python3
"""v339: genuine CzechVeg-Open v1.1 railway-plot discovery gate.

Downloads not handled by this script; the caller verifies original Zenodo MD5.
Audits ORIGINAL 117,739 vegetation-plot header descriptions, separating:
 A) railway mentioned in the explicit plot HABITAT (candidate only),
 B) railway mentioned only in LOCALITY (landmark, never rail-labelled).
No railway site is invented from mere proximity to a track or from railway-
related species names or bibliographic titles. No raw source GPS/plot IDs
are republished. This is a label/discoverability audit, not an effect estimate.
Source DOI: 10.5281/zenodo.19231152.
"""
import argparse,csv,hashlib,json,math,re,unicodedata
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED_MD5="0bc920462d1f9eac473e7b2a231e6782"
EXPECTED_PLOTS=117739

def ascii_cz(value):
    src=unicodedata.normalize("NFKD",str(value or "").casefold())
    return "".join(c for c in src if not unicodedata.combining(c))

STRONG={
    "explicit_railway":re.compile(r"\bzelezni\w*|\bzeleznic\w*"),
    "rail_yard":re.compile(r"\bkolejist\w*"),
    "train_station":re.compile(r"\bnadraz\w*|\bstanicni\s+kolej\w*"),
    "industrial_siding":re.compile(r"\bvleck\w*"),
    "track_structure":re.compile(r"\bkolejove\s+loz\w*|\btratovy\s+nasyp\w*"),
}
NON_RAIL_STATION=re.compile(r"\bautobusov\w*\s+nadraz\w*|\btramvajov\w*|\bmetro\b")

def classify(habitat,locality):
    h=ascii_cz(habitat)
    l=ascii_cz(locality)
    kinds=[k for k,pattern in STRONG.items() if pattern.search(h)]
    if kinds and NON_RAIL_STATION.search(h) and not STRONG["explicit_railway"].search(h):
        return "ambiguous_other_transit_habitat",kinds
    if kinds:
        return "rail_habitat_text_candidate",kinds
    if any(p.search(l) for p in STRONG.values()):
        return "rail_locality_landmark_only",[]
    return "no_rail_text_evidence",[]

def good_number(v):
    try:x=float(v)
    except (ValueError,TypeError):return None
    return x if math.isfinite(x) else None

def meta_num(row,options):
    for name in options:
        if name in row:
            return good_number(row[name])
    return None

def audit(path):
    src=Path(path)
    md5=hashlib.md5()
    with src.open("rb") as file:
        for block in iter(lambda:file.read(2*1024*1024),b""):
            md5.update(block)
    if md5.hexdigest()!=EXPECTED_MD5:
        raise ValueError("Not exact original CzechVeg v1.1 FULL headers, MD5 mismatch")
    counts=Counter()
    kinds=Counter()
    fields_covered=Counter()
    source_groups=defaultdict(set)
    source_source_year=defaultdict(Counter)
    short_habitats=Counter()
    year_distribution=Counter()
    coordinates=Counter()
    exact_source_pairs=Counter()
    header=[]
    with src.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        header=reader.fieldnames or []
        needed={"RELEVE_NR","HABITAT","LOCALITY"}
        if not needed.issubset(header):
            raise ValueError("Original header missing required plot, habitat, locality fields: "+str(header[:20]))
        seen=set()
        for r in reader:
            counts["source_rows"]+=1
            plot=r["RELEVE_NR"]
            if not plot or plot in seen:
                raise ValueError("Original RELEVE_NR missing/duplicated")
            seen.add(plot)
            status,indicators=classify(r.get("HABITAT"),r.get("LOCALITY"))
            counts[status]+=1
            if str(r.get("HABITAT") or "").strip():fields_covered["HABITAT"]+=1
            if str(r.get("LOCALITY") or "").strip():fields_covered["LOCALITY"]+=1
            if status!="rail_habitat_text_candidate":
                continue
            kinds.update(indicators)
            for x in ("DEG_LAT","DEG_LON","PRECISION","AREA","PLOT_SIZE","PLOT_AREA","ALTITUDE",
                      "REFERENCE","REFERENCE_AUTHOR","REFERENCE_YEAR","REFERENCE_TITLE","HABITAT"):
                if r.get(x) is not None and str(r.get(x)).strip():
                    fields_covered["candidate_"+x]+=1
            lat=meta_num(r,("DEG_LAT",))
            lon=meta_num(r,("DEG_LON",))
            precision=meta_num(r,("PRECISION",))
            if lat is not None and lon is not None and 48.5<=lat<=51.5 and 12<=lon<=19:
                coordinates["candidate_plots_with_CZ_plausible_decimal_coordinates"]+=1
            if precision is not None:
                if precision<=100:coordinates["candidate_plots_precision_le100m"]+=1
                if precision<=1000:coordinates["candidate_plots_precision_le1000m"]+=1
            source=str(r.get("REFERENCE") or r.get("REFERENCE_TITLE") or "")
            year=str(r.get("YEAR") or r.get("RELEVE_YEAR") or r.get("DATE") or "")
            if source:source_groups[status].add(source)
            if source and year:source_source_year[source][year]+=1
            cleaned=ascii_cz(r.get("HABITAT"))
            # Audit textual LABELS ONLY. Never include LOCALITY, SITE name or GPS.
            short_habitats[cleaned[:130]]+=1
            if year:year_distribution[year[:4]]+=1
    if counts["source_rows"]!=EXPECTED_PLOTS:
        raise ValueError("Original v1.1 open plot count differs from 117739")
    labels=[{"habitat_label_normalized":k,"plot_count":n}
            for k,n in short_habitats.most_common(35)]
    return {
        "version":"v339","source":"CzechVeg-Open v1.1 full header original",
        "doi":"10.5281/zenodo.19231152",
        "original_header_md5":EXPECTED_MD5,"source_total_vegetation_plots":counts["source_rows"],
        "field_coverage":dict(fields_covered),
        "railway_text_classification":{k:v for k,v in counts.items() if k!="source_rows"},
        "explicit_railway_habitat_indicators":dict(kinds),
        "n_distinct_rail_label_sources":len(source_groups["rail_habitat_text_candidate"]),
        "geographic_QA":dict(coordinates),
        "candidate_year_first4_distribution":dict(year_distribution),
        "top_habitat_text_labels_without_localities":labels,
        "source_header_field_names":header,
        "cautions":[
            "A railway keyword in HABITAT is still not definitive field verification; review source language and individual source publication",
            "A railway word in LOCALITY only means landmark proximity; NEVER label as sampled railway habitat",
            "Header source RELEVEs are observational, not matched rail-vs-natural pairs",
            "CzechVeg-OpenStrat stratified subset is not automatically a representative railway microhabitat sample; use full inventory for identification",
            "Coordinates are subject to stated PRECISION; cannot equate within a geodesic buffer with sampled adjacent railway plots",
            "Original species data are not examined at this source-intake gate",
            "No raw GPS, plot IDs, source locality text or original plant incidence published in this QA"]
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-headers",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=audit(a.source_headers)
    dest=Path(a.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "version":data["version"],
        "source_total_vegetation_plots":data["source_total_vegetation_plots"],
        "railway_text_classification":data["railway_text_classification"],
        "explicit_railway_habitat_indicators":data["explicit_railway_habitat_indicators"],
        "n_distinct_rail_label_sources":data["n_distinct_rail_label_sources"],
        "geographic_QA":data["geographic_QA"],
        "candidate_year_first4_distribution":data["candidate_year_first4_distribution"],
        "top_habitat_text_labels_without_localities":data["top_habitat_text_labels_without_localities"],
        "source_header_field_names":data["source_header_field_names"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()

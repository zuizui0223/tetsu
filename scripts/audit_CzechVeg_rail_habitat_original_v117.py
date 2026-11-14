#!/usr/bin/env python3
"""v11.7: Original CzechVeg-Open v1.1 railway habitat plot discovery.

Original 2026 Czech Vegetation Database v1.1, DOI 10.5281/zenodo.19231152.
Scientific authors Chytry et al. (2026) and Masaryk University.
Public source: 117739 plot header CSV records, species CSV separately
available. This first stage verifies Zenodo metadata SHA/MD5, original
plot header uniqueness, and strictly textual original railway habitat
evidence in HABITAT and LOCALITY. It does NOT infer rail location from
mere proximity to tracks or from a bibliographic paper title.

No biological effect fit; no Railway-vs-Road comparability presumed.
No full original Czech vegetation dataset redistribution in artifact.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

DOI="10.5281/zenodo.19231152"
VERSION="1.1"
HEADERS_NAME="CzechVeg-Open_headers_v1-1_2026-03-26.csv"
EXPECTED_PLOTS=117739
CHECK_FIELDS=("RELEVE_NR","SURF_AREA","HABITAT","LOCALITY","REFERENCE_TITLE")
# Conservative original site-context indicators only, NOT railway GIS proximity.
RAIL=re.compile(
    r"(?:železni(?:c|č|t|h)|železniční|nádraž|kolejiš|kolejov|"
    r"\brailway\b|\brailroad\b|\brail(?:\s|-)?track\b|\bbahnhof\b|"
    r"\bgleisanlag|\bgleisbett\b|\bgleisfeld\b)",
    re.I)
# Separate controls for transport road language; do not identify as rail.
ROAD=re.compile(
    r"(?:silnic|dálnic|cestn[íi]|vozov|roadside|roadside|"
    r"\broad\b|\bhighway\b|\bautobusov)",
    re.I)
STRONG_RAIL=re.compile(
    r"(?:železni(?:c|č|t|h)|železniční|kolejiš|kolejov|"
    r"\brailway\b|\brailroad\b|\brail(?:\s|-)?track\b|"
    r"\bgleisanlag|\bgleisbett\b|\bgleisfeld\b)",
    re.I)
CONFOUND_BUS=re.compile(r"(?:autobus|bus station|bus depot|autobusov)",re.I)


def original_zenodo_entry(meta,name=HEADERS_NAME):
    records=meta.get("files",[])
    if not isinstance(records,list) or len(records)==0:
        raise ValueError("Official Zenodo record did not expose files")
    files=[x for x in records if x.get("key")==name]
    if len(files)!=1:
        raise ValueError("Expected exactly one genuine original CzechVeg header CSV")
    rec=files[0]
    check=rec.get("checksum")
    if not isinstance(check,str) or ":" not in check:
        raise ValueError("Official Zenodo checksum unavailable")
    algo,digest=check.split(":",1)
    if algo not in ("md5","sha256"):
        raise ValueError("Unexpected Zenodo original checksum algorithm")
    if len(digest)!=(32 if algo=="md5" else 64):
        raise ValueError("Publisher original checksum malformed")
    size=rec.get("size")
    if not isinstance(size,int) or size<500000:
        raise ValueError("Published source size too small for 117k plot headers")
    url=(rec.get("links") or {}).get("self")
    if not isinstance(url,str) or not url.startswith("https://zenodo.org/"):
        raise ValueError("Original publisher file URL not official Zenodo")
    return {"source_key":name,"checksum_algorithm":algo,
            "checksum":digest,"source_bytes":size,"official_url":url}


def verify_original(meta,headers_file):
    if str(meta.get("id"))!="19231152":
        raise ValueError("Official dataset ID differs from original source DOI")
    # DOI may be a concept DOI but exact record is pinned by version+file hash.
    item=original_zenodo_entry(meta)
    raw=Path(headers_file).read_bytes()
    if len(raw)!=item["source_bytes"]:
        raise ValueError("Source original CzechVeg header file size mismatch")
    digest=hashlib.new(item["checksum_algorithm"],raw).hexdigest()
    if digest!=item["checksum"]:
        raise ValueError("Original CzechVeg header source checksum mismatch")
    item["independent_original_file_SHA256"]=hashlib.sha256(raw).hexdigest()
    return item


def text(row,key):
    return str(row.get(key) or "").strip()


def plot_source_scan(headers_file):
    counts=Counter()
    original_id=set()
    regions=Counter()
    years=Counter()
    area={"strong_rail":[],"other":[]}
    samples=[]
    examples_by_field=Counter()
    top_titles=Counter()
    strict_rail_plot_ids=[]
    with open(headers_file,encoding="utf-8-sig",newline="") as f:
        head=f.readline()
        f.seek(0)
        delimiter=";" if head.count(";")>head.count(",") else ","
        reader=csv.DictReader(f,delimiter=delimiter)
        cols=reader.fieldnames or []
        if set(CHECK_FIELDS)-set(cols):
            raise ValueError("No original CzechVeg declared habitat/location/area/header fields")
        for row in reader:
            counts["original_header_rows"]+=1
            rid=text(row,"RELEVE_NR")
            if not rid or rid in original_id:
                raise ValueError("Missing/duplicated original Czech vegetation plot id")
            original_id.add(rid)
            h=text(row,"HABITAT")
            loc=text(row,"LOCALITY")
            ref=text(row,"REFERENCE_TITLE")
            hay=h+" "+loc
            rr=bool(RAIL.search(hay)) and not bool(CONFOUND_BUS.search(hay))
            rs=bool(STRONG_RAIL.search(hay)) and not bool(CONFOUND_BUS.search(hay))
            road=bool(ROAD.search(hay))
            reference=bool(RAIL.search(ref))
            counts["original_plot_source_habitat_nonempty"]+=bool(h)
            counts["original_plot_source_locality_nonempty"]+=bool(loc)
            counts["original_plot_source_reference_title_nonempty"]+=bool(ref)
            counts["plot_rail_source_site_text_broad"]+=rr
            counts["plot_rail_source_site_text_strict"]+=rs
            counts["plot_original_road_source_site_text"]+=road
            counts["plot_rail_source_site_and_road_text"]+=rr and road
            counts["plot_rail_only_in_reference_title_not_site"]+=reference and not rr
            if rs:
                strict_rail_plot_ids.append(rid)
                if len(samples)<24:
                    samples.append({
                        "plot_id":rid,"year":text(row,"YEAR"),
                        "plot_area_m2":text(row,"SURF_AREA"),
                        "habitat_excerpt":h[:120],"locality_excerpt":loc[:120],
                        "source_reference_title_excerpt":ref[:120],
                    })
                top_titles[ref[:140] or "(none)"]+=1
                examples_by_field["habitat"]+=bool(STRONG_RAIL.search(h))
                examples_by_field["locality"]+=bool(STRONG_RAIL.search(loc))
            for label,active in (("strong_rail",rs),("other",not rs)):
                if not active:continue
                raw=text(row,"SURF_AREA")
                try:a=float(raw)
                except (TypeError,ValueError):a=float("nan")
                if math.isfinite(a) and a>0:
                    counts["positive_area_"+label]+=1
                    area[label].append(a)
            y=text(row,"YEAR")
            if rs:years[y or "(missing)"]+=1
            if rs:regions[text(row,"SKALICKY1988") or "(missing)"]+=1
    if counts["original_header_rows"]!=EXPECTED_PLOTS or len(original_id)!=EXPECTED_PLOTS:
        raise ValueError("Original CzechVeg v1.1 expected 117739 plots not present")
    return {
        "n_original_plot_rows":counts["original_header_rows"],
        "n_original_railway_site_text_strict_potential_not_independently_field_confirmed":
            counts["plot_rail_source_site_text_strict"],
        "n_original_railway_site_text_broad_including_station_ambiguous":
            counts["plot_rail_source_site_text_broad"],
        "n_original_road_source_site_text":counts["plot_original_road_source_site_text"],
        "n_original_railway_candidate_road_overlap_text":counts["plot_rail_source_site_and_road_text"],
        "n_reference_only_railway_mentioned_not_accepted_as_rail_plot":
            counts["plot_rail_only_in_reference_title_not_site"],
        "n_header_rows_with_habitat":counts["original_plot_source_habitat_nonempty"],
        "n_header_rows_with_locality":counts["original_plot_source_locality_nonempty"],
        "n_header_rows_with_reference_title":counts["original_plot_source_reference_title_nonempty"],
        "strict_candidate_with_measured_positive_source_plot_area_m2":
            counts["positive_area_strong_rail"],
        "other_source_plots_with_measured_positive_plot_area_m2":
            counts["positive_area_other"],
        "strict_candidate_original_area_median_m2":
            median(area["strong_rail"]) if area["strong_rail"] else None,
        "noncandidate_original_area_median_m2":
            median(area["other"]) if area["other"] else None,
        "strict_candidate_by_year":dict(years.most_common(20)),
        "strict_candidate_original_phytogeographic_region_code_top20":dict(regions.most_common(20)),
        "strict_candidate_source_title_top20":dict(top_titles.most_common(20)),
        "strict_candidate_source_evidence_fields":dict(examples_by_field),
        "first_24_strong_rail_keyword_source_plot_examples":samples,
        "csv_columns":cols,
        "n_unique_source_ids":len(original_id),
        "no_railway_plant_species_matrix_or_causal_model_fitted":True,
        "strict_rail_candidate_ids_preserved_only_as_count_not_all_ids":True,
    }


def run(meta,headers_file):
    proof=verify_original(meta,headers_file)
    scan=plot_source_scan(headers_file)
    return {
      "complete":True,"version":"v11.7",
      "original_source_title":"Czech Vegetation Database v1.1 (CzechVeg-Open)",
      "publisher_doi":DOI,
      "source_organization":"Masaryk University Czech Vegetation Database, original Chytry et al. (2026)",
      "original_source_data_url":"https://zenodo.org/records/19231152",
      "source_header_file":proof,
      "original_habitat_and_location_rail_keyword_screen":scan,
      "critical_limits":[
         "Keyword match in original site description is a source candidate, NOT proof a surveyed polygon is within an active railway right of way.",
         "Railway reference-title-only plots are NOT counted as on-rail plots.",
         "Rail, road and arbitrary nonrail vegetation plot designs differ. No matched rail/road distance or area comparison identified.",
         "Species CSV and taxon strategy harmonization must be provenance-verified before biological effect estimates.",
         "Authors and source original 2026 DOI must be cited; no raw 117739-plot data redistributed in this report.",
      ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zenodo-record",required=True)
    p.add_argument("--source-headers",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    meta=json.loads(Path(a.zenodo_record).read_text(encoding="utf-8"))
    result=run(meta,a.source_headers)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    summary=result["original_habitat_and_location_rail_keyword_screen"].copy()
    summary.pop("csv_columns",None)
    summary.pop("first_24_strong_rail_keyword_source_plot_examples",None)
    print(json.dumps({"complete":True,"version":"v11.7",
        "original_source_file":result["source_header_file"],
        "rail_habitat_source_audit":summary},ensure_ascii=False),flush=True)


if __name__=="__main__":main()

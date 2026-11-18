#!/usr/bin/env python3
"""v28.7: Source-authenticated 235,979-species life-cycle labels for real Kashmir rail plots.

Retrieves Figshare public article 20764777 (Poppenwimer et al. Nature 2023)
metadata, checks file checksums, and uses only unambiguous exact species labels.
No genus-only imputation, prediction from other traits or unpublished rail rows.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import re
import sys
import urllib.request
import zipfile
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_Kashmir_external_alpine_lifespan_v286 import (
    ALPINE_DOI, binomial, crosswalk, coverage_summary,
)
from analyze_Kashmir_annual_microsites_v285 import (
    life_span_crosswalk, original_taxon_headers, summarize,
)
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix

FIGSHARE_ITEM_ID=20764777
FIGSHARE_META="https://api.figshare.com/v2/articles/"+str(FIGSHARE_ITEM_ID)
FIGSHARE_DOI="https://doi.org/10.6084/m9.figshare.20764777"
PAPER_DOI="https://doi.org/10.1038/s41586-023-06644-x"

def fetch(url: str,max_bytes: int=120_000_000) -> bytes:
    req=urllib.request.Request(url,headers={"User-Agent":"tetsu academic open data verification (+https://github.com/zuizui0223/tetsu)"})
    with urllib.request.urlopen(req, timeout=150) as response:
        buf=io.BytesIO()
        while True:
            chunk=response.read(1024*1024)
            if not chunk: break
            buf.write(chunk)
            if buf.tell()>max_bytes: raise ValueError("Figshare source exceeded size gate")
        return buf.getvalue()


def source_metadata_and_file() -> tuple[dict,bytes]:
    metadata=json.loads(fetch(FIGSHARE_META,max_bytes=5_000_000))
    files=metadata.get("files") or []
    choices=[f for f in files if f.get("name","").lower().endswith((".csv",".tsv",".txt",".zip"))]
    if len(choices)!=1:
        raise ValueError("Figshare life-cycle study requires exactly one supported original data file: "
                         + json.dumps([{"name":x.get("name"),"size":x.get("size")} for x in files]))
    source=choices[0]
    if (not source.get("download_url") or source.get("size",0)>120_000_000):
        raise ValueError("Source file too large or missing actual public Figshare URL")
    blob=fetch(source["download_url"])
    if len(blob)!=source["size"]:
        raise ValueError("Downloaded Figshare file length differs from original metadata")
    if source.get("computed_md5") and hashlib.md5(blob).hexdigest()!=source["computed_md5"]:
        raise ValueError("Original public Figshare MD5 mismatch")
    if source.get("supplied_md5") and hashlib.md5(blob).hexdigest()!=source["supplied_md5"]:
        raise ValueError("Original supplied Figshare MD5 mismatch")
    fileinfo={"file_id":source.get("id"),"original_name":source["name"],
              "sha256":hashlib.sha256(blob).hexdigest(),
              "size_bytes":len(blob),"original_computed_md5":source.get("computed_md5"),
              "original_source_article_id":FIGSHARE_ITEM_ID,
              "figshare_article_doi":metadata.get("doi"),
              "paper":PAPER_DOI}
    return fileinfo,blob


def read_rows(data: bytes,name: str) -> tuple[list[dict],list[str]]:
    if name.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            files=[f for f in z.infolist() if f.filename.lower().endswith((".csv",".tsv",".txt"))]
            if len(files)!=1:
                raise ValueError("Expected exactly one text life-cycle data file in original ZIP")
            name=files[0].filename
            data=z.read(name)
    if not name.lower().endswith((".csv",".tsv",".txt")):
        raise ValueError("Source not CSV/TSV/TXT")
    try: text=data.decode("utf-8-sig")
    except UnicodeDecodeError:
        text=data.decode("latin-1")
    header=text[:40000]
    dialect=csv.Sniffer().sniff(header,delimiters=",;\t")
    reader=csv.DictReader(io.StringIO(text),dialect=dialect)
    result=list(reader)
    if not 100_000 <= len(result) <= 1_000_000:
        raise ValueError("Expected life-cycle database >100,000 unique species rows; found "+str(len(result)))
    return result,list(reader.fieldnames or [])


def infer_source_columns(columns: list[str]) -> tuple[str | None,str | None,dict]:
    cleaned={c: re.sub("[^a-z0-9]","",c.casefold()) for c in columns}
    names={"species","speciesname","acceptedname","taxon","taxonname",
           "standardizedname","scientificname","resolvedname","wfoname",
           "plantname","plant"}
    spans={"lifecycle","lifespan","annualperennial","annualperennialstate",
           "annualperennialclassification","lifehistory","lifecycledb",
           "plantlifecycle","lifeform","lifeformtype"}
    name_hits=[c for c,v in cleaned.items() if v in names]
    span_hits=[c for c,v in cleaned.items() if v in spans]
    # Explicit ambiguity triggers a schema audit rather than a guessed data join.
    return (name_hits[0] if len(name_hits)==1 else None,
            span_hits[0] if len(span_hits)==1 else None,
            {"original_columns":columns,"taxon_candidates":name_hits,"life_cycle_candidates":span_hits})


def parse_life_cycles(rows: list[dict], name_col: str, span_col: str) -> tuple[list[dict],dict]:
    value_counts=collections.Counter(str(r.get(span_col,"")).strip().casefold() for r in rows)
    allowed={"annual":"annual","perennial":"perennial"}
    # Authors' database may code life cycle numerically: do not guess numeric meaning.
    unknown=[v for v,n in value_counts.items() if v and v not in allowed]
    if unknown:
        raise ValueError("Original Figshare life-cycle coding needs manual protocol verification; values "
                         +str(value_counts.most_common(20)))
    species=[]
    for r in rows:
        name=(r.get(name_col) or "").strip()
        span=(r.get(span_col) or "").strip().casefold()
        if not name or span not in allowed: continue
        # Poppenwimer Figshare Accepted Name column is all lower case (WFO).
        # Restore the genus initial solely for strict two-token taxon parsing.
        normalized_name = name[:1].upper() + name[1:] if name and name[0].islower() else name
        species.append({"name":normalized_name,"binomial":binomial(normalized_name),
                        "synonym":"","span":span})
    return species,{"source_rows":len(rows),
                    "life_cycle_records":len(species),
                    "consensus_values":dict(value_counts)}


def original_meta(meta: dict):
    # no download URLs copied into report unless verified
    return {k:v for k,v in meta.items() if k in
            {"file_id","original_name","sha256","size_bytes","original_computed_md5",
             "original_source_article_id","figshare_article_doi","paper"}}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rail",type=Path,required=True)
    ap.add_argument("--result",type=Path,required=True)
    ap.add_argument("--crosswalk",type=Path,required=True)
    args=ap.parse_args()
    original=original_taxon_headers(args.rail)
    plots,raw_meta=parse_xlsx_matrix(args.rail)
    meta,blob=source_metadata_and_file()
    rows,columns=read_rows(blob,meta["original_name"])
    namecol,spancol,schema=infer_source_columns(columns)
    result={"version":"v28.7","source":"Poppenwimer et al. 2023, Nature",
            "paper":PAPER_DOI,"life_cycle_data":original_meta(meta),
            "rail_data_sha256":raw_meta["original_file_sha256"],
            "global_life_cycle_is_trait_source_not_independent_railway_study":True,
            "schema":schema,"field_effect_status":"NOT_EVALUATED"}
    if namecol and spancol:
        # unknown categorical encoding fails closed; the report records data schema instead.
        try:
            donor,sourceinfo=parse_life_cycles(rows,namecol,spancol)
            original_binomials={binomial(n) for n in original if binomial(n)}
            donor_binomials={r["binomial"] for r in donor if r["binomial"]}
            print("LIFEFORM_SCHEMA_DIAGNOSTIC "+json.dumps({
                "original_first_names":original[:12],
                "global_first_names":[r["name"] for r in donor[:12]],
                "global_life_cycle_count":len(donor),
                "original_binomial_count":len(original_binomials),
                "global_binomial_count":len(donor_binomials),
                "binomial_intersection_count":len(original_binomials & donor_binomials),
                "global_example_common_species":[r["name"] for r in donor if "poa annua" in r["name"].lower()][:2]
            },ensure_ascii=False),flush=True)
            mapping,mapping_meta=crosswalk(original,donor,allow_binomial=True)
            mapping_meta["trait_identification"]="Poppenwimer et al. 2023 Nature global plant life-cycle database"
            exact,_=crosswalk(original,donor,allow_binomial=False)
            args.crosswalk.parent.mkdir(parents=True,exist_ok=True)
            with args.crosswalk.open("w",encoding="utf-8",newline="") as fp:
                wr=csv.DictWriter(fp,fieldnames=["original_taxon","life_span","trait_source"])
                wr.writeheader()
                for record in mapping:
                    record=dict(record)
                    record["trait_source"]=PAPER_DOI+" | "+FIGSHARE_DOI+" | "+record["trait_source"].split(" | ")[-2]+" | "+record["trait_source"].split(" | ")[-1]
                    wr.writerow(record)
            classified,crosswalkmeta=life_span_crosswalk(args.crosswalk,original)
            result.update({"source_species_summary":sourceinfo,
                           "match":mapping_meta,
                           "strict_exact_fullname_count":len(exact),
                           "original_matched_site_trait_coverage":coverage_summary(plots,classified),
                           "strict_70pct":summarize(plots,classified,min_coverage=.7,bootstrap=1999),
                           "strict_50pct":summarize(plots,classified,min_coverage=.5,bootstrap=1999),
                           "field_effect_status":"EVALUATED_WITH_SITE_COVERAGE_GATE"})
        except ValueError as exc:
            result["field_effect_status"]="SCHEMA_OR_CATEGORY_FAIL_CLOSED"
            result["source_schema_issue"]=str(exc)[:400]
    else:
        result["field_effect_status"]="SCHEMA_AMBIGUOUS_FAIL_CLOSED"
    args.result.parent.mkdir(parents=True,exist_ok=True)
    args.result.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"version":result["version"],
                      "figshare_source":meta["original_name"],
                      "source_file_sha256":meta["sha256"],
                      "schema":schema,
                      "result_status":result["field_effect_status"],
                      "original_taxa_match":result.get("match",{}).get("classified_original_taxa"),
                      "strict_primary":result.get("strict_70pct",{}).get("comparisons",{}).get(
                              "Track__minus__Interior Plot 1",{}).get("all"),
                      "strict_50_percent_sensitivity":result.get("strict_50pct",{}).get("comparisons",{}),
                      "by_zone_trait_coverage":result.get("original_matched_site_trait_coverage",{}),
                      "matched_species":result.get("match",{}),
                      "output":str(args.result)},ensure_ascii=False))


if __name__=="__main__":
    main()

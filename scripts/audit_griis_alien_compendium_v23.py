#!/usr/bin/env python3
"""v2.3: validate the COUNTRY-LEVEL non-native species GRIIS compendium source.

2022 static 196-country checklist; historical baseline, NOT 2026 status.
Source GRIIS = *reported introduced taxa in countries*, not railway-seeded,
and absence from the checklist MUST NOT imply native status.
No raw third-party CSV, names or country-level sensitive sites committed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import requests

ZENODO_RECORD=6348164
API=f"https://zenodo.org/api/records/{ZENODO_RECORD}"
DOI="10.5281/zenodo.6348164"
EXPECTED_MD5="6cdbed57a49f49d696ece545908baf1f"
EXPECTED_FILENAME="GRIIS - Country Compendium V1_0.csv"
COUNTRY_NAMES={
    "JP":"Japan", "KR":"South Korea","TW":"Taiwan","US":"United States",
    "CA":"Canada","GB":"United Kingdom","DE":"Germany","FR":"France",
    "AU":"Australia","NZ":"New Zealand",
}


def pick_csv_source(metadata):
    if not isinstance(metadata,dict):
        raise ValueError("Zenodo metadata must be JSON")
    files=metadata.get("files")
    if not isinstance(files,list):
        raise ValueError("Zenodo files list missing")
    found=[f for f in files if str(f.get("key","")).strip()==EXPECTED_FILENAME]
    if len(found)!=1:
        raise ValueError("Expected exactly one original versioned compendium CSV")
    f=found[0]
    checksum=str(f.get("checksum","")).lower()
    if checksum!=f"md5:{EXPECTED_MD5}":
        raise ValueError(f"Zenodo archived CSV checksum changed: {checksum}")
    url=f.get("links",{}).get("self")
    if not isinstance(url,str) or not url.startswith("https://"):
        raise ValueError("Zenodo file content link not verified")
    return url


def parse_source(raw):
    if hashlib.md5(raw).hexdigest()!=EXPECTED_MD5:
        raise ValueError("GRIIS dataset bytes differ from archived Zenodo MD5")
    text=raw.decode("utf-8-sig")
    reader=csv.DictReader(io.StringIO(text))
    headers=reader.fieldnames
    if not headers or len(set(headers))!=len(headers):
        raise ValueError("Bad/duplicate CSV header")
    if not any("species" in h.lower() for h in headers):
        raise ValueError("No GRIIS species identifier column")
    countries=[h for h in headers if "country" in h.lower()]
    taxon=[h for h in headers if any(x in h.lower() for x in (
        "gbif","taxon","species","scientific","accepted","kingdom","rank"))]
    if not countries:
        raise ValueError("No source country column: cannot make country-specific alien list")
    frequencies={h:Counter() for h in countries}
    kingdoms={h:Counter() for h in headers if h.lower()=="kingdom"}
    rows=0;name_present=0
    names=[h for h in headers if h.lower()=="species"]
    for item in reader:
        rows+=1
        for col in countries:
            value=str(item.get(col) or "").strip()
            if value: frequencies[col][value]+=1
        for col in kingdoms:
            value=str(item.get(col) or "").strip()
            if value: kingdoms[col][value]+=1
        if names and str(item.get(names[0]) or "").strip():
            name_present+=1
    if rows<10000:
        raise ValueError(f"Unexpected too-small country compendium {rows}")
    results={}
    for col,count in frequencies.items():
        vals=dict(count)
        results[col]={
            "n_nonblank_source_country_labels":sum(count.values()),
            "n_distinct_country_labels":len(count),
            "source_country_label_examples":count.most_common(15),
            "target_country_exact_label_counts":{
                code:sum(n for label,n in vals.items() if label.lower() in {
                    code.lower(),COUNTRY_NAMES[code].lower()
                })
                for code in COUNTRY_NAMES
            }
        }
    return {
        "source_country_columns":countries,
        "source_taxon_name_and_key_columns":taxon,
        "full_header":headers,
        "n_original_country_taxon_rows":rows,
        "n_nonempty_species_binomial":name_present,
        "country_label_quality":results,
        "kingdom_column_source_counts":{k:dict(v) for k,v in kingdoms.items()},
        "no_native_status_can_be_inferred_from_list_nonmembership":True,
    }


def get_json(client,url):
    r=client.get(url,timeout=(15,90))
    r.raise_for_status()
    return r.json()


def fetch_source(client,url):
    r=client.get(url,timeout=(30,150))
    r.raise_for_status()
    if len(r.content)>30_000_000 or len(r.content)<100_000:
        raise ValueError("Unexpected dataset download size")
    return r.content


def run():
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-macro-biodiversity-introduced-country-list-v2.3",
                           "Accept":"application/json"})
    metadata=get_json(client,API)
    url=pick_csv_source(metadata)
    raw=fetch_source(client,url)
    report=parse_source(raw)
    report.update({
        "source_metadata_api":API,
        "source_citation_doi":DOI,
        "source_file":EXPECTED_FILENAME,
        "source_md5":hashlib.md5(raw).hexdigest(),
        "source_bytes":len(raw),
        "source_collection_year":"2022 static compendium; does not guarantee status in 2010 or 2026",
        "audit_utc":datetime.now(timezone.utc).isoformat(),
        "quality_gate":"metadata checksum and downloaded CSV verified; taxonomy crosswalk NOT YET performed",
        "not_plant_community_biological_effect":True,
    })
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/griis_country_compendium_schema_v23.json")
    args=p.parse_args()
    res=run()
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:res[k] for k in (
        "source_citation_doi","source_md5","source_bytes","n_original_country_taxon_rows",
        "source_country_columns","source_taxon_name_and_key_columns","kingdom_column_source_counts"
    )},ensure_ascii=False),flush=True)
    for col,detail in res["country_label_quality"].items():
        print(json.dumps({"country_column":col,"n_distinct_labels":detail["n_distinct_country_labels"],
                          "examples":detail["source_country_label_examples"],
                          "target_exact_names":detail["target_country_exact_label_counts"]},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

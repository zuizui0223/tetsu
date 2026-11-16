#!/usr/bin/env python3
"""v2.4: country-specific GRIIS introduced species keys via indexed GBIF checklists.

Use independently citable country checklist datasets and SOURCE name usages,
retaining only directly backbone-matched species-level Plantae usage keys.
Introduced = listed by GRIIS for that COUNTRY; unlisted = UNKNOWN, not native.
US study is Wisconsin, within 48 contiguous states; use US lower-48 list.
Taiwan has no entry in the audited 2022 GRIIS country compendium and remains
UNASSESSED; never classify Taiwan as zero introduced species.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import requests
import time

DATASETS={
    "JP":("5c5a6e45-d510-45ed-b7bf-6d0624fac056","Japan"),
    "KR":("72340f0c-a018-4c3e-8aa5-88c0ffb2386f","Republic of Korea"),
    "US":("32ad19ed-6b89-447a-9242-795c0897f345","United States (Contiguous)"),
    "CA":("b95e74e0-b772-430c-a729-9d56ce0182e2","Canada"),
    "GB":("1288ee7d-d67c-4e23-8d95-409973067383","United Kingdom"),
    "DE":("2a281675-c184-4d54-a47e-5d3a7ec6653f","Germany"),
    "FR":("692a0aec-70e2-4038-ac9c-9a1d7acc025f","France"),
    "AU":("15147db1-27c3-49b5-9c69-c78d55a4b8ff","Australia"),
    "NZ":("b286f8ee-1e29-4f8a-a041-9106bb1796f3","New Zealand"),
}
UNASSESSABLE_COUNTRIES={"TW":"no adequate country-specific GRIIS dataset has been verified"}
GBIF="https://api.gbif.org/v1"
PAGE=300
MAX_PAGES=60


def source_params(dataset,offset):
    if not isinstance(offset,int) or offset<0 or offset%PAGE:
        raise ValueError("Source rank pagination offset must align with PAGE")
    return {"datasetKey":dataset,"origin":"SOURCE","rank":"SPECIES",
            "limit":PAGE,"offset":offset}


def verify_dataset(metadata,dataset,country_hint):
    if not isinstance(metadata,dict) or metadata.get("key")!=dataset:
        raise ValueError("Wrong GBIF registry dataset source")
    title=str(metadata.get("title",""))
    if "Global Register of Introduced and Invasive Species" not in title:
        raise ValueError("GBIF checklist is not GRIIS")
    if country_hint.lower() not in title.lower():
        raise ValueError("Country/region source-title mismatch")
    if metadata.get("type") not in (None,"CHECKLIST"):
        raise ValueError("Incorrect GBIF dataset type")
    return {"title":title,"modified":metadata.get("modified"),
            "doi_or_reference":metadata.get("doi") or metadata.get("doiReference"),
            "license":metadata.get("license")}


def parse_species_page(payload,offset):
    if not isinstance(payload,dict):
        raise ValueError("Malformed checklist response")
    count=payload.get("count")
    if isinstance(count,bool) or not isinstance(count,int) or count<0:
        raise ValueError("No GBIF checklist count")
    rows=payload.get("results")
    if not isinstance(rows,list) or len(rows)>PAGE:
        raise ValueError("Invalid checklist page result")
    result=[]
    for r in rows:
        if not isinstance(r,dict):
            raise ValueError("Non-object checklist usage")
        srcid=r.get("key")
        if isinstance(srcid,bool) or not isinstance(srcid,int):
            raise ValueError("Missing source usage ID")
        result.append(r)
    end=(offset+len(rows)>=count) or (payload.get("endOfRecords") is True)
    if not end and len(rows)==0:
        raise ValueError("Empty nonfinal checklist page")
    return count,result,end


def get_source(client,url,params=None):
    for attempt in range(5):
        try:
            res=client.get(url,params=params,timeout=(15,90))
            res.raise_for_status()
            return res.json()
        except (requests.RequestException,ValueError):
            if attempt==4: raise
            time.sleep(min(12,2**attempt))
    raise RuntimeError("GBIF source retry exceeded")


def audit_country(client,country,get_json=get_source,sleep=time.sleep):
    if country not in DATASETS:
        raise ValueError("Country lacks a source-verified GRIIS country checklist")
    uuid,hint=DATASETS[country]
    registry=verify_dataset(get_json(client,f"{GBIF}/dataset/{uuid}"),uuid,hint)
    seen=set(); plantkeys=set(); unmatched=0; rows_seen=0
    rank_mismatch=0; nonplants=0; has_nub=0
    total=None;page_count=0
    for page in range(MAX_PAGES):
        offset=page*PAGE
        n,rows,end=parse_species_page(
            get_json(client,f"{GBIF}/species/search",source_params(uuid,offset)),offset)
        if total is None:total=n
        elif total!=n:raise ValueError("Checklist source population count changed across pages")
        for r in rows:
            srcid=r["key"]
            if srcid in seen:
                raise ValueError("Duplicated source usage ID across pages: pagination not stable")
            seen.add(srcid);rows_seen+=1
            if str(r.get("rank","")).upper()!="SPECIES":
                rank_mismatch+=1;continue
            if str(r.get("kingdom","")).lower()!="plantae":
                nonplants+=1;continue
            key=r.get("nubKey")
            if isinstance(key,bool) or not isinstance(key,int) or key<=0:
                unmatched+=1;continue
            has_nub+=1
            plantkeys.add(str(key))
        page_count+=1
        if end:
            if len(seen)!=total:
                raise ValueError("GBIF source count and pages disagree; incomplete alien list")
            break
        sleep(.13)
    else:
        raise ValueError("Too many checklist pages; incomplete country species source")
    if not rows_seen:
        raise ValueError("No country checklist source records")
    return {
        "country":country,"dataset_uuid":uuid,"source_title":registry["title"],
        "source_dataset_reference":f"https://www.gbif.org/dataset/{uuid}",
        "source_registry":registry,
        "n_indexed_source_species_rank_usage_rows":rows_seen,
        "n_source_plant_species_usage_rows_with_backbone_key":has_nub,
        "n_directly_matched_plant_backbone_species_keys":len(plantkeys),
        "n_unmatched_plant_name_usages":unmatched,
        "n_nonplantae":nonplants,"n_rank_mismatch":rank_mismatch,
        "page_count":page_count,"full_source_pages_verified":True,
        "gbif_backbone_plant_species_key_set":sorted(plantkeys),
        "not_in_list_means_unknown_not_native":True,
        "country_localized_introduced_checklist_status_not_railway_arrival":True,
    }


def run(countries):
    if not countries or len(countries)!=len(set(countries)) or any(c not in DATASETS for c in countries):
        raise ValueError("Invalid fixed country subset")
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-railway-plants-macro-country-GRIIS-keys-v2.4"})
    report={"utc":datetime.now(timezone.utc).isoformat(),
            "scope":"GBIF indexed GRIIS alien country list of source-backed species keys, NOT occurrence analysis",
            "requested_countries":countries,
            "unassessed":{"TW":"No verified GRIIS alien country list; cannot treat nonmembership as native"},
            "countries":{},"errors":{}}
    for c in countries:
        try:
            result=audit_country(client,c)
            report["countries"][c]=result
            print(json.dumps({"country":c,"source_rows":result["n_indexed_source_species_rank_usage_rows"],
                              "plant_source_nub_rows":result["n_source_plant_species_usage_rows_with_backbone_key"],
                              "alien_backbone_keys":result["n_directly_matched_plant_backbone_species_keys"],
                              "unmatched_plant_usages":result["n_unmatched_plant_name_usages"],
                              "pages":result["page_count"]},ensure_ascii=False),flush=True)
        except Exception as exc:
            report["errors"][c]=f"{type(exc).__name__}: {str(exc)[:220]}"
            print(json.dumps({"country":c,"source_error":report["errors"][c]},ensure_ascii=False),flush=True)
    report["complete"]=not report["errors"] and len(report["countries"])==len(countries)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--countries",nargs="+",choices=list(DATASETS),default=["JP","GB"])
    p.add_argument("--output",default="results/audit/griis_9country_alien_backbone_keys_v24.json")
    args=p.parse_args()
    res=run(args.countries)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":res["complete"],"n_countries":len(res["countries"]),
                      "errors":res["errors"]},ensure_ascii=False),flush=True)
    if not res["complete"]:
        raise SystemExit("GBIF GRIIS checklist key source incomplete")


if __name__=="__main__":
    main()

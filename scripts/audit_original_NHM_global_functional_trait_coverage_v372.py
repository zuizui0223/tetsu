#!/usr/bin/env python3
"""v372: full original 59,369-species NHM ADEPT global trait completeness audit.

Genuine Natural History Museum London CKAN 2026 functional angiosperm
dataset DOI 10.5519/p3dm31kc, 2026 combined.csv read-only DataStore API.
Direct NHM bulk CSV source can return 403; official paginated DataStore
metadata is not a synthetic replacement. Repeat source collection with a
pinned official resource UUID, check exact SOURCE total and unique species,
and aggregate trait availability for global functional convergence feasibility.

No original taxon rows, names, individual descriptions or textual trait
passages redistributed. These are worldwide angiosperm trait annotations,
NOT a worldwide railway-plant-trait analysis.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time,urllib.error,urllib.parse,urllib.request
from collections import Counter
from pathlib import Path

SOURCE_RESOURCE="2552d157-9519-476c-995d-906e2433cc19"
DOI="10.5519/p3dm31kc"
ENDPOINT="https://data.nhm.ac.uk/api/3/action/datastore_search"
PAGE=2500
EXPECTED_TOTAL=59369
COLUMNS=("taxon","wfo","life form","habit","dispersal mode",
        "plant min_ height [m]","plant max_ height [m]",
        "seed min_ length [mm]","seed max_ length [mm]",
        "seed min_ diameter [mm]","seed max_ diameter [mm]",
        "reproduction system","pollination")
LIFE=("life form","habit")
DISP=("dispersal mode",)
NUMERIC=("plant min_ height [m]","plant max_ height [m]",
          "seed min_ length [mm]","seed max_ length [mm]",
          "seed min_ diameter [mm]","seed max_ diameter [mm]")

def source_value(v):
    if v is None:return False
    if isinstance(v,(list,dict)):return bool(v)
    s=str(v).strip().lower()
    return s not in ("","na","n/a","null","none","nan",
                     "not available","not reported","no data","unknown","?","-")

def query(offset,limit=PAGE,tries=5):
    url=ENDPOINT+"?"+urllib.parse.urlencode(
       {"resource_id":SOURCE_RESOURCE,"limit":limit,"offset":offset})
    req=urllib.request.Request(url,headers={
       "User-Agent":"tetsu-v372-world-trait-source-coverage/1.0",
       "Accept":"application/json"})
    last=None
    for i in range(tries):
        try:
            with urllib.request.urlopen(req,timeout=100) as resp:
                if resp.status!=200:raise ValueError("Official NHM API returned non-200")
                raw=resp.read(25_000_000)
            out=json.loads(raw)
            if out.get("success") is not True:
                raise ValueError("Official NHM DataStore source API error")
            result=out.get("result",{})
            if not isinstance(result.get("records"),list):
                raise ValueError("Original source NHM did not return records")
            return result,hashlib.sha256(raw).hexdigest()
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            last=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (408,429,500,502,503,504):
                raise
            time.sleep(min(18,2**i))
    raise RuntimeError("Original NHM DataStore source permanently unavailable") from last

def source_traits():
    seen_ids=set();seen_names=set()
    summaries=Counter()
    fields=None
    page_hashes=[]
    masks=Counter()
    genus=Counter()
    first_total=None
    offset=0
    while True:
        page,h=query(offset)
        total=page.get("total")
        if not isinstance(total,int) or total!=EXPECTED_TOTAL:
            raise ValueError(f"Original global NHM source total changed: {total}")
        if first_total is None:first_total=total
        current=[x["id"] for x in page.get("fields",[])
                 if x.get("id")!="_id"]
        if fields is None:
            fields=tuple(current)
            if not set(COLUMNS).issubset(fields):
                raise ValueError("Original NHM global plant trait source schema drift")
        elif tuple(current)!=fields:
            raise ValueError("Original NHM source pagination columns differ")
        records=page["records"]
        if not records:break
        if len(records)>PAGE:
            raise ValueError("Original source datastore pagination ignored limit")
        page_hashes.append(h)
        for row in records:
            summaries["source_original_trait_records"]+=1
            rid=row.get("_id")
            if rid in seen_ids:raise ValueError("Original NHM dataset page repeats row ID")
            seen_ids.add(rid)
            taxon=str(row.get("taxon") or "").strip()
            if not taxon:raise ValueError("NHM original trait record has no taxon name")
            if taxon in seen_names:summaries["source_duplicate_taxon_names"]+=1
            seen_names.add(taxon)
            for col in COLUMNS:
                if source_value(row.get(col)):
                    summaries["has_"+col]+=1
            li=any(source_value(row.get(col)) for col in LIFE)
            disp=source_value(row.get("dispersal mode"))
            height=any(source_value(row.get(col)) for col in NUMERIC[:2])
            seed=any(source_value(row.get(col)) for col in NUMERIC[2:])
            if li:summaries["source_any_life_form_or_habit"]+=1
            if li and disp:summaries["source_lifeform_AND_dispersal_trait_same_row"]+=1
            if li and height:summaries["source_lifeform_AND_height_same_row"]+=1
            if disp and seed:summaries["source_dispersal_AND_seed_dimensions_same_row"]+=1
            if li and disp and height:summaries["source_lifeform_disp_height_same_row"]+=1
            mask=(int(li),int(disp),int(height),int(seed))
            masks["".join(map(str,mask))]+=1
        offset+=len(records)
        if summaries["source_original_trait_records"]>=EXPECTED_TOTAL:break
        if offset%10000< PAGE:print(json.dumps({
           "NHM_original_source_trait_rows_validated_so_far":offset}),flush=True)
    if summaries["source_original_trait_records"]!=EXPECTED_TOTAL:
        raise ValueError("Full world source 59369 NHM original taxon records not received")
    if len(seen_ids)!=EXPECTED_TOTAL:
        raise ValueError("Original NHM canonical source IDs incomplete")
    return {
       "version":"v372-original-fullworld-59k-NHM-ADEPT-angiosperm-functional-trait-completeness",
       "source_original_DOI":DOI,
       "original_source_official_CKAN_DataStore_resource_UUID":SOURCE_RESOURCE,
       "source_original_complete_plant_trait_records":EXPECTED_TOTAL,
       "source_original_distinct_taxon_names":len(seen_names),
       "source_original_distinct_source_row_ids":len(seen_ids),
       "source_original_page_count":len(page_hashes),
       "source_official_DataStore_JSON_page_SHA256_manifest":hashlib.sha256(
           "\n".join(page_hashes).encode()).hexdigest(),
       "source_84_original_traits_and_taxon_schema":list(fields),
       "original_source_trait_availability":dict(summaries),
       "original_source_joint_trait_availability_mask_life_dispersal_height_seed":{
           key:value for key,value in sorted(masks.items())},
       "original_source_lifeform_dispersal_height_joint_fraction":round(
           summaries["source_lifeform_disp_height_same_row"]/EXPECTED_TOTAL,7),
       "source_trait_missing_is_NOT_a_negative_plant_trait_assignment":True,
       "source_plant_trait_DISTRIBUTION_not_crossmatched_to_global_iNat_taxon_ids":True,
       "source_trait_values_may_be_textual_literature_extractions_NOT_HARMONIZED_ECOLOGICAL_GUILDS":True,
       "original_raw_taxon_names_descriptions_and_plant_traits_NOT_published":True,
       "source_study_scope":"worldwide original functional angiosperm source, not railway-flora field evidence"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=source_traits()
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:result[k] for k in (
       "version","source_original_complete_plant_trait_records",
       "source_original_distinct_taxon_names","original_source_trait_availability",
       "original_source_lifeform_dispersal_height_joint_fraction",
       "original_source_joint_trait_availability_mask_life_dispersal_height_seed")},
       ensure_ascii=False),flush=True)

if __name__=="__main__":main()

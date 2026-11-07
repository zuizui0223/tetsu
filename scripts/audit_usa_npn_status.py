#!/usr/bin/env python3
"""Small real USA-NPN phenophase yes/no feasibility audit (no railway effect).

The source has repeated individual-level status observations, unlike
opportunistic flower photos. POST application/x-www-form-urlencoded mirrors
the official rnpn implementation. Do not infer rail disturbance or
phenological shifts. Downloads constrained to named species, state and year.

See https://usa-npn.github.io/rnpn/reference/npn_download_status_data.html
"""
from __future__ import annotations
import argparse
import collections
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import requests

ENDPOINT = "https://services.usanpn.org/npn_portal/observations/getObservations.ndjson"
# Public NPN species identifiers documented by rnpn:
# 82 = Liriodendron tulipifera (tuliptree), 210 = Carnegiea gigantea (saguaro)
# Initial single-state pilot only; NOT all United States
PILOTS = (
    {"species_id":82, "expected_species":"Liriodendron tulipifera", "state":"MD", "year":2024},
    {"species_id":210, "expected_species":"Carnegiea gigantea", "state":"AZ", "year":2024},
)
MAX_DOWNLOAD_BYTES = 22_000_000
MAX_RECORDS = 60_000
MISSING = {None, "", "-9999", -9999, "NA", "null"}


def form_args(case:dict)->dict:
    year=int(case["year"]);species=int(case["species_id"])
    state=str(case["state"])
    if not 2009 <= year <= 2025 or species < 1 or len(state)!=2:
        raise ValueError("Invalid bounded USA-NPN query")
    return {
        "request_src":"tetsu_research_feasibility",
        "climate_data":"0",
        "species_id[1]":str(species),
        "state[1]":state,
        "start_date":f"{year}-01-01",
        "end_date":f"{year}-12-31",
    }


def parse_stream(response)->tuple[list[dict],str,int]:
    """NDJSON, hard max, no partial summary accepted on overflow."""
    hasher=hashlib.sha256()
    records=[]
    size=0
    for line in response.iter_lines():
        if not line:
            continue
        size+=len(line)+1
        if size>MAX_DOWNLOAD_BYTES:
            raise ValueError("NPN download exceeds safety cap")
        hasher.update(line+b"\n")
        obj=json.loads(line)
        if not isinstance(obj,dict):
            raise ValueError("NPN NDJSON line is not an observation object")
        records.append(obj)
        if len(records)>MAX_RECORDS:
            raise ValueError("NPN observation count exceeds pilot cap")
    return records,hasher.hexdigest(),size


def value_counts(records,columns):
    out={}
    for key in columns:
        vals=collections.Counter(str(r.get(key)) for r in records)
        out[key]=dict(vals.most_common(14))
    return out


def summarize(records,case,digest,size):
    all_keys=sorted({k for r in records for k in r})
    status_cols=[c for c in all_keys if "status" in c.lower() or "phenophase" in c.lower()]
    site_cols=[c for c in all_keys if "station" in c.lower() or "site" in c.lower()]
    indiv_cols=[c for c in all_keys if "individual" in c.lower() or "plant_id" in c.lower()]
    id_cols=site_cols+indiv_cols
    if not records:
        return {"case":case,"records":0,"source_sha256":digest,"source_bytes":size,"column_names":[]}
    # Never read a missing status as a biological no. These counts are raw
    # field values awaiting source-schema confirmation.
    dates=[str(r.get("observation_date"))[:10] for r in records if r.get("observation_date") not in MISSING]
    ids={}
    for key in sorted(set(id_cols)):
        ids[key]=len({str(r[key]) for r in records if r.get(key) not in MISSING})
    keys_to_count=sorted(set(status_cols))
    result={
        "case":case,"records":len(records),"source_sha256":digest,
        "source_bytes":size,"column_names":all_keys,
        "candidate_status_columns":status_cols,
        "status_raw_value_counts":value_counts(records,keys_to_count[:16]),
        "id_cardinalities":ids,"first_and_last_observation_date":
            [min(dates),max(dates)] if dates else None,
        "sample_non_identifier_field_types":{
            k:type(records[0].get(k)).__name__
            for k in all_keys if not any(x in k.lower() for x in ("latitude","longitude","id","name","comment","email"))
        },
        "notes":[
            "Raw 0/1/-1 status entries are not directly comparable to iNaturalist photo labels.",
            "Source reports one row per phenophase observation of monitored plant; dedupe only with verified schema.",
            "NPN site coordinates may lack precision for railway verge assignment; do not release exact individual/site coordinates in artifact.",
            "Phenology differences require matched same-species rail/road/background, climate/GDD and spatial observer structure.",
        ],
    }
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/usa_npn_status_feasibility.json")
    args=p.parse_args()
    result={"created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "api":ENDPOINT,"status":"API source/field feasibility only",
            "cases":[],"failures":{}}
    for case in PILOTS:
        key=f"{case['state']}_{case['year']}_spp{case['species_id']}"
        try:
            with requests.post(ENDPOINT,data=form_args(case),
                               timeout=(20,140),stream=True,
                               headers={"User-Agent":"tetsu-biodiversity-phenology-audit/0.6"}) as response:
                response.raise_for_status()
                records,digest,size=parse_stream(response)
            report=summarize(records,case,digest,size)
            result["cases"].append(report)
            print(json.dumps({"case":key,"n":report["records"],
                              "columns":report.get("candidate_status_columns",[])},
                             ensure_ascii=False),flush=True)
        except Exception as exc:
            result["failures"][key]=f"{type(exc).__name__}: {str(exc)[:300]}"
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"case_counts":[x["records"] for x in result["cases"]],
                      "failure_codes":result["failures"]},indent=2))
    if not result["cases"]:
        raise SystemExit("No usable USA-NPN API response; no observational inference is possible")


if __name__=="__main__":
    main()

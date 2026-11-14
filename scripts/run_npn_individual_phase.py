#!/usr/bin/env python3
"""Real repeated-plant Open-flowers diagnostics from the verified NPN API.

Reports phenophase yes/no counts, actual observations supporting interval-
censored first detected flowers, and monitored flowering interruption.
No railway geometry or environmental interpretation in this script.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path
import requests

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.npn_individuals import (
    normalize_open_flowers,summarize_individuals,coverage_summary
)
from audit_usa_npn_status import (ENDPOINT,PILOTS,form_args,parse_stream,
                                  MAX_RECORDS,MAX_DOWNLOAD_BYTES)


def run():
    report={"created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "source_api":ENDPOINT,
            "meaning":"real within-individual phenophase feasibility, not railway, mowing or adaptation",
            "cases":[],"failures":{}}
    for case in PILOTS:
        key=f"{case['state']}_{case['year']}_spp{case['species_id']}"
        try:
            with requests.post(ENDPOINT,data=form_args(case),stream=True,
                               timeout=(20,140),headers={
                                   "User-Agent":"tetsu-biodiversity-phenology-audit/0.6"
                               }) as response:
                response.raise_for_status()
                raw,digest,size=parse_stream(response)
            clean,excluded=normalize_open_flowers(raw)
            individual=summarize_individuals(clean)
            coverage=coverage_summary(raw,clean,excluded,individual)
            # Dataset includes repeated phenophases and other site series.
            # Do not publish precise site or individual identifiers.
            coverage["case"]=case
            coverage["source_sha256"]=digest
            coverage["api_content_size_bytes"]=size
            coverage["n_individual_years_with_at_least_four_status_visits"]=sum(
                p["n_dates"]>=4 for p in individual
            )
            coverage["valid_prior_no_bracket_width_days"]={
                "1-7":sum(p["first_yes_bracket_days"] is not None and
                          p["first_yes_bracket_days"]<=7 for p in individual),
                "8-14":sum(p["first_yes_bracket_days"] is not None and
                           8<=p["first_yes_bracket_days"]<=14 for p in individual),
                "15-21":sum(p["first_yes_bracket_days"] is not None and
                            15<=p["first_yes_bracket_days"]<=21 for p in individual),
            }
            report["cases"].append(coverage)
            print(json.dumps({"case":key,"open_flower_rows":coverage["open_flower_rows_raw"],
                              "valid_yes":coverage["yes_visit_rows"],
                              "valid_no":coverage["no_visit_rows"],
                              "sites":coverage["n_distinct_sites"],
                              "individuals":coverage["n_distinct_individuals"],
                              "first_yes_brackets":coverage["individual_years_with_bracketed_first_yes"],
                              "recurrent_candidate":coverage["individual_years_with_2no_interrupted_flowering"]},
                             ensure_ascii=False),flush=True)
        except Exception as exc:
            report["failures"][key]=f"{type(exc).__name__}: {str(exc)[:260]}"
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/usa_npn_individual_flowering.json")
    a=p.parse_args()
    summary=run()
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    if not summary["cases"]:
        raise SystemExit("No valid real NPN phase data; see audit errors")
    print(f"success_cases={len(summary['cases'])} failed_cases={len(summary['failures'])}")


if __name__=="__main__":
    main()

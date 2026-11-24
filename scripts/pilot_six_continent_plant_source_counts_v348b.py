#!/usr/bin/env python3
"""v348b fast independent source accessibility and completeness count.
Exactly same pre-frozen six geographic windows/2025 plants as v348; no
plant records or species identities exported. No railway effect inferred.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.parse
from pathlib import Path
from pilot_self_extracted_six_continent_rail_flora_v348 import (
    CITIES,YEAR,RATE_DELAY_S,INAT,inat_parameters,get_bytes,
)
def run():
    records=[]
    for i,(name,continent,country,lat,lon) in enumerate(CITIES):
        if i:time.sleep(RATE_DELAY_S)
        p=inat_parameters(lat,lon,1)
        p["per_page"]=1
        url=INAT+"?"+urllib.parse.urlencode(p)
        try:
            b,final,status=get_bytes(url,timeout=40)
            j=json.loads(b)
            n=j.get("total_results")
            if not isinstance(n,int) or n<0 or not isinstance(j.get("results"),list):
                raise ValueError("invalid original iNaturalist API observation count schema")
            records.append({"city":name,"country":country,"continent":continent,
                "year":YEAR,"source_http":status,"source_total_public_research_grade_plant_observations":n,
                "within_original_3600_census_cap":n<=3600,
                "official_iNaturalist_original_results_schema_valid":True,
                "source_response_sha256":hashlib.sha256(b).hexdigest()})
        except Exception as e:
            records.append({"city":name,"country":country,"continent":continent,
                "year":YEAR,"source_status":"unavailable_NOT_zero",
                "error":type(e).__name__,"http":getattr(e,"code",None),
                "message":str(e)[:160]})
    return {"version":"v348b-iNat-source-count-six-continent",
       "predeclared_geography_year_identical_to_v348":True,
       "worldwide_biological_effects_NOT_calculated":True,
       "regions":records,
       "source_complete_count_regions":sum("source_total_public_research_grade_plant_observations" in x for x in records),
       "records_or_individual_sites_not_exported":True}
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args()
    obj=run();out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(obj,ensure_ascii=False),flush=True)

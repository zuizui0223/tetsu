#!/usr/bin/env python3
"""v331 source-only geographical metadata readiness for Kashmir railway sites.

Uses exact original Dryad version (SHA pinned in v254); no original site GPS,
species or station identities are exported. Do not assume row geographical
columns are parseable without examining the actual publisher header.
"""
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path
from statistics import mean,median
from openpyxl import load_workbook

from analyze_Kashmir_same_sites_reference_year_v254 import DRYAD_SHA256

def num(x):
    try: y=float(x)
    except (ValueError,TypeError):return None
    return y if math.isfinite(y) else None

def audit(path):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=DRYAD_SHA256 or len(raw)!=122187:
        raise ValueError("Original Dryad SHA or bytes mismatch")
    book=load_workbook(path,read_only=True,data_only=True)
    try:
        rows=book.worksheets[0].values
        header=next(rows)
        if len(header)!=129:raise ValueError("Unexpected original columns")
        meta=defaultdict(list)
        for row in rows:
            s=" ".join(str(row[0] or "").split())
            if not s:raise ValueError("Missing source site")
            meta[s].append(tuple(row[:6]))
    finally:book.close()
    if len(meta)!=31 or sum(s.startswith("Between ") for s in meta)!=15:
        raise ValueError("Kashmir 31 and between15 expected")
    variables={}
    for idx in (1,2,3):
        observed=[]
        per_site=[]
        for site,rr in meta.items():
            vals=[num(r[idx]) for r in rr]
            if None in vals:
                variables[idx]={"numeric_complete":False,
                     "header":str(header[idx])};break
            if max(vals)-min(vals)>1e-7:
                variables[idx]={"numeric_complete":True,
                     "stable_by_source_site":False,
                     "max_within_site_spread":max(vals)-min(vals),
                     "header":str(header[idx])};break
            per_site.append(vals[0]);observed.extend(vals)
        else:
            variables[idx]={"numeric_complete":True,
                 "stable_by_source_site":True,
                 "header":str(header[idx]),
                 "min":min(per_site),"max":max(per_site),
                 "median":median(per_site),
                 "site_unique_n":len(set(per_site))}
    return {"version":"v331","original_SHA256":DRYAD_SHA256,
       "source_header_first_6":[str(s) for s in header[:6]],
       "source_site_n":len(meta),
       "station_sites":16,"between_sites":15,
       "source_location_field_checks":variables,
       "coordinates_at_specific_source_sites_NOT_exported":True,
       "genuinely_local_nonrail_control_status":"Interior Plots1/2 sampled in same site",
       "no_climate_mechanism_or_species_evidence_estimated":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=audit(a.source)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out),flush=True)
if __name__=="__main__":main()

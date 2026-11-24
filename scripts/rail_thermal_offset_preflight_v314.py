#!/usr/bin/env python3
"""Preflight railway thermal-offset data without fitting an ecological model.

Input CSV columns:
 study,region,site_id,habitat,species,incidence,temperature_indicator,
 regional_temperature,trait_source
habitat must be 'rail' or 'nonrail'. Each site_id must have both habitats.
Missing or nonnumeric temperature indicators remain missing; never coerce 'x' to zero.
Usage: python scripts/rail_thermal_offset_preflight_v314.py input.csv output.csv
"""
import csv
import math
import sys
from collections import defaultdict

REQUIRED = {"study","region","site_id","habitat","species","incidence",
            "temperature_indicator","regional_temperature","trait_source"}

def numeric(s):
    try:
        x=float(s)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None

def run(src,dst):
    with open(src,newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError("Missing columns: "+", ".join(sorted(REQUIRED-set(reader.fieldnames or []))))
        groups=defaultdict(list)
        for row in reader:
            if row["habitat"] not in ("rail","nonrail"):
                raise ValueError("Unknown habitat")
            w=numeric(row["incidence"])
            if w is None or w < 0: raise ValueError("Invalid incidence")
            groups[(row["study"],row["region"],row["site_id"],row["habitat"])].append(row)
    sites=defaultdict(dict)
    for (study,region,site,habitat), rows in groups.items():
        if len({r["species"] for r in rows}) != len(rows):
            raise ValueError("Duplicate species in site/habitat")
        temps={numeric(r["regional_temperature"]) for r in rows}
        if len(temps)!=1 or None in temps: raise ValueError("Missing/inconsistent regional temperature")
        total=sum(float(r["incidence"]) for r in rows)
        covered=sum(float(r["incidence"]) for r in rows if numeric(r["temperature_indicator"]) is not None)
        weighted=sum(float(r["incidence"])*numeric(r["temperature_indicator"])
                     for r in rows if numeric(r["temperature_indicator"]) is not None)
        sites[(study,region,site)][habitat]={
            "cwm":weighted/covered if covered else None,
            "coverage":covered/total if total else None,
            "regional_temperature":next(iter(temps)),
            "n_species":len(rows)}
    out=[]
    for (study,region,site), habitats in sorted(sites.items()):
        if set(habitats)!={"rail","nonrail"}: continue
        a,b=habitats["rail"],habitats["nonrail"]
        if a["regional_temperature"]!=b["regional_temperature"]:
            raise ValueError("Regional climate differs within matched site")
        out.append({"study":study,"region":region,"site_id":site,
          "regional_temperature":a["regional_temperature"],
          "rail_cwm":a["cwm"],"nonrail_cwm":b["cwm"],
          "delta_cwm":a["cwm"]-b["cwm"] if a["cwm"] is not None and b["cwm"] is not None else "",
          "rail_coverage":a["coverage"],"nonrail_coverage":b["coverage"],
          "rail_n_species":a["n_species"],"nonrail_n_species":b["n_species"]})
    if not out: raise ValueError("No matched rail/nonrail site pairs")
    with open(dst,"w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0]))
        writer.writeheader();writer.writerows(out)
    return len(out)

if __name__=="__main__":
    if len(sys.argv)!=3: raise SystemExit(__doc__)
    print("Matched pairs:",run(sys.argv[1],sys.argv[2]))

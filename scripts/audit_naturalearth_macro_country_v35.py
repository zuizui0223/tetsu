#!/usr/bin/env python3
"""v3.5: Natural Earth 10m COUNTRY polygon fraction for frozen 50km cells.

Non-biological validity check: clipping 2500km² equal-area geographic
source cells by independent source-country polygons before ecological
inference. COUNTRY POLYGON fraction is not a vegetation/terrestrial
cover map and may include inland water; it is a conservative country
jurisdiction filter. Geographic/country map quality is assessed without
using rail connections or alien plant assemblage response.

Natural Earth 10m admin0 countries, public domain. Actual ZIP fetched
temporarily and SHA256 recorded, not republished in GitHub artifacts.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime,timezone
from hashlib import sha256
import io,json,math
from pathlib import Path
from statistics import median
import zipfile

import requests
import shapefile
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform as shapely_transform, unary_union
from shapely.validation import make_valid

from audit_gbif_equalarea_macroblocks_v29 import (
    BLOCKS,PILOT,GRID_CRS
)

ZIP_URL="https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"
COUNTRY_FRAC_MIN=0.80
MAX_ZIP_BYTES=15_000_000


def source_countries(raw):
    if len(raw)<500000 or len(raw)>MAX_ZIP_BYTES or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("Natural Earth admin0 source is missing/truncated/unexpected")
    with zipfile.ZipFile(io.BytesIO(raw)) as bundle:
        names={p.lower():p for p in bundle.namelist()}
        candidates=[p for p in bundle.namelist() if p.lower().endswith("_admin_0_countries.shp")]
        if len(candidates)!=1:
            raise ValueError("Unexpected natural earth admin0 SHP source schema")
        stem=candidates[0][:-4]
        def part(ext):
            member=names.get((stem+ext).lower())
            if not member:raise ValueError("Missing Natural Earth shapefile component "+ext)
            return io.BytesIO(bundle.read(member))
        with shapefile.Reader(shp=part(".shp"),shx=part(".shx"),dbf=part(".dbf"),
                              encoding="latin1") as rdr:
            columns=[c[0] for c in rdr.fields[1:]]
            if not ("ISO_A2" in columns and "ADMIN" in columns):
                raise ValueError("Unexpected country ISO_A2/ADMIN labels")
            byiso=defaultdict(list)
            to_eq=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
            codes={country for country,_,_ in BLOCKS.values()}
            for sr in rdr.iterShapeRecords():
                props=dict(zip(columns,sr.record))
                code=str(props.get("ISO_A2") or "").strip()
                # Natural Earth sometimes marks France ISO_A2=-99 despite
                # an explicit ADMIN=France / ADM0_A3=FRA. Match that
                # documented sovereign-country identity narrowly; never
                # classify an unknown country by proximity.
                if (code not in codes and
                    str(props.get("ADMIN") or "").strip()=="France" and
                    str(props.get("ADM0_A3") or "").strip()=="FRA"):
                    code="FR"
                if code not in codes:continue
                src=shape(sr.shape.__geo_interface__)
                if src.is_empty:continue
                if not src.is_valid:src=make_valid(src)
                eq=shapely_transform(to_eq.transform,src)
                if not eq.is_valid:eq=make_valid(eq)
                byiso[code].append(eq)
    if set(byiso)!=codes:
        raise ValueError("Country polygons missing required source iso2 codes: "+str(codes-set(byiso)))
    result={c:unary_union(v) for c,v in byiso.items()}
    if any(p.is_empty or p.area<=0 for p in result.values()):
        raise ValueError("Empty country polygon")
    return result


def _cell_geometry(block):
    country,lon,lat=BLOCKS[block]
    tr=Transformer.from_crs(4326,GRID_CRS,always_xy=True)
    x,y=tr.transform(lon,lat)
    full=200000;cell=50000
    from shapely.geometry import box
    out={}
    for j in range(4):
        for i in range(4):
            out[f"{block}_x{i}_y{j}"]=box(x-full/2+i*cell,y-full/2+j*cell,
                 x-full/2+(i+1)*cell,y-full/2+(j+1)*cell)
    return out


def summarize_block(block,geom,coverage):
    if block not in BLOCKS or geom.is_empty:
        raise ValueError("Invalid macro country source block")
    country=BLOCKS[block][0]
    items=coverage["blocks"][block]["cells"]
    geos=_cell_geometry(block)
    if set(items)!=set(geos):
        raise ValueError("GBIF sample cells and Natural Earth geographic cells differ")
    results={}
    for cellid,poly in geos.items():
        ratio=min(1.,max(0.,poly.intersection(geom).area/poly.area))
        entry=items[cellid]
        if entry["country"]!=country:
            raise ValueError("Country mismatch in prior GBIF coverage")
        period_ok=all(v.get("source_eligible") is True for v in entry["period_coverage"].values())
        if len(entry["period_coverage"])!=2:
            raise ValueError("Incomplete two-period plant source")
        results[cellid]={
            "source_country_polygon_fraction":round(ratio,6),
            "country_fraction_min0p80":ratio>=COUNTRY_FRAC_MIN,
            "prior_GBIF_both_period_coverage":period_ok,
            "country_and_prior_GBIF_source_qualified":ratio>=COUNTRY_FRAC_MIN and period_ok,
        }
    vals=[v["source_country_polygon_fraction"] for v in results.values()]
    return {
        "country":country,
        "n_2500km2_source_cells":len(results),
        "n_source_country_fraction_ge0p80":sum(v["country_fraction_min0p80"] for v in results.values()),
        "n_prior_GBIF_both_period_eligible":sum(v["prior_GBIF_both_period_coverage"] for v in results.values()),
        "n_country_and_GBIF_eligible":sum(v["country_and_prior_GBIF_source_qualified"] for v in results.values()),
        "country_fraction_min":min(vals),
        "country_fraction_median":round(median(vals),6),
        "cells":results,
    }


def run(coverage,raw,blocks):
    if not coverage.get("complete") or coverage.get("source_errors"):
        raise ValueError("Prior GBIF 50km country source incomplete")
    if not blocks or any(b not in BLOCKS for b in blocks) or len(set(blocks))!=len(blocks):
        raise ValueError("Source blocks not predefined")
    if not set(blocks).issubset(coverage.get("blocks",{})):
        raise ValueError("Prior GBIF source region omitted")
    polys=source_countries(raw)
    out={
        "utc":datetime.now(timezone.utc).isoformat(),
        "source":ZIP_URL,
        "source_license":"Natural Earth public-domain map data",
        "source_sha256":sha256(raw).hexdigest(),
        "source_zip_bytes":len(raw),
        "country_polygon_fraction_threshold":COUNTRY_FRAC_MIN,
        "regions":{}, "errors":{},
    }
    for block in blocks:
        try:
            c=BLOCKS[block][0]
            summary=summarize_block(block,polys[c],coverage)
            out["regions"][block]=summary
            print(json.dumps({
                "block":block,"country":c,
                "country_cells_ge80pct":summary["n_source_country_fraction_ge0p80"],
                "both_period_GBIF_source":summary["n_prior_GBIF_both_period_eligible"],
                "source_GBIF_and_country_qualified":summary["n_country_and_GBIF_eligible"],
                "median_country_fraction":summary["country_fraction_median"]
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            out["errors"][block]=f"{type(exc).__name__}: {str(exc)[:200]}"
            print(json.dumps({"source_error":block,"detail":out["errors"][block]},ensure_ascii=False),flush=True)
    out["complete"]=len(out["regions"])==len(blocks) and not out["errors"]
    out["summary"]={
        "n_2500km2_source_cells":16*len(blocks),
        "n_country_fraction_eligible_cells":sum(r["n_source_country_fraction_ge0p80"] for r in out["regions"].values()),
        "n_country_and_plant_record_eligible_cells":sum(r["n_country_and_GBIF_eligible"] for r in out["regions"].values())
    }
    out["limits"]=[
        "Country polygon fraction is NOT equivalent to terrestrial vegetation fraction; large inland lakes may remain.",
        "Natural Earth portrays disputed national boundaries from one map viewpoint and has geometric simplifications at 1:10m.",
        "Some border cells may host genuine plant records in the designated country but are intentionally excluded by >=80 percent surface rule.",
        "This 80 percent gate is a source quality sensitivity, not a preregistered biological hypothesis.",
        "No train transport, rail network, introduced flora response, or individual source GBIF occurrence rows are used in the country polygon screen."
    ]
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--gbif-coverage",required=True)
    a.add_argument("--blocks",nargs="+",choices=list(BLOCKS),default=list(PILOT))
    a.add_argument("--output",default="results/audit/50km_naturalearth_country_fraction_v35.json")
    args=a.parse_args()
    cov=json.loads(Path(args.gbif_coverage).read_text(encoding="utf-8"))
    response=requests.get(ZIP_URL,timeout=(20,120))
    response.raise_for_status()
    r=run(cov,response.content,args.blocks)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"complete":r["complete"],"summary":r["summary"],
                      "errors":r["errors"]},ensure_ascii=False),flush=True)
    if not r["complete"]:
        raise SystemExit("Country source fraction incomplete")
if __name__=="__main__":
    main()

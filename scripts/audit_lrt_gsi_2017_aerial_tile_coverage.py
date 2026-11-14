#!/usr/bin/env python3
"""2016/17/18 historical GSI annual orthophoto coverage over LRT pixels.

Official aerial photography is the only source here potentially fine enough
to visually verify which 2017 pixels were roads versus vegetation at 20m.
This script checks actual official PNG tile pixels and their year-labelled
coverage, NOT auto-classified historical roads or plant species.

Does NOT publish licensed N05 rail geometry, third-party photo bytes,
route coordinates, private imagery or road classification. All sampled
annual-photo tiles are fetched transiently and aggregated.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path

import requests
from PIL import Image
from pyproj import Transformer

from audit_n05_historical_railway import bounded_zip_download
from measure_utsunomiya_s2_ndvi_pixels import (
    get_lrt_lines,build_fixed_pairs
)
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.rail_verge_pixel_eligibility import near_points_from_140m

TEMPLATE="https://cyberjapandata.gsi.go.jp/xyz/nendophoto{year}/{z}/{x}/{y}.png"
GSI_PAGE="https://maps.gsi.go.jp/development/ichiran.html"
SEAMLESS="https://cyberjapandata.gsi.go.jp/xyz/seamlessphoto/{z}/{x}/{y}.jpg"
YEARS=(2016,2017,2018,2024)
ZOOM=17
MAX_SAMPLE_POINTS=36
# Every fourth along-rail anchor in the exact 48-anchor source frame.
ANCHOR_MODULUS=4


def world_tile(lon,lat,zoom=ZOOM):
    if not 130<=lon<=150 or not 20<=lat<=47:
        raise ValueError("unexpected Japanese source location")
    n=2**zoom
    latrad=math.radians(lat)
    xtile=int((lon+180)/360*n)
    ytile=int((1-math.asinh(math.tan(latrad))/math.pi)/2*n)
    return zoom,xtile,ytile


def image_is_nonblank(raw,expected_format='PNG'):
    """Do not interpret pure transparent/solid tiles as valid aerial photos."""
    if len(raw)>2_000_000:
        raise ValueError("oversized GSI individual tile")
    with Image.open(io.BytesIO(raw)) as im:
        if im.format!=expected_format or im.size not in ((256,256),(512,512)):
            raise ValueError("unexpected GSI tile format/size")
        im=im.convert("RGBA")
        # Inspect a deterministic grid rather than cherry-pick vegetation
        # pixels. Aerial PHOTO coverage is never a road/no-road label.
        w,h=im.size
        grid=[im.getpixel((int((i+.5)*w/16),int((j+.5)*h/16)))
              for i in range(16) for j in range(16)]
        visible=[rgb[:3] for rgb in grid if rgb[3]>5]
        if len(visible)<200:
            return "mostly_transparent_or_missing"
        unique=len(set(visible))
        if unique<=2:
            return "solid_color_not_orthophoto"
        return "photo_pixels_present"


def choose_points():
    source=bounded_zip_download()
    lines,_=get_lrt_lines(source)
    originals,pairs,geo=build_fixed_pairs(lines,n_anchors=48)
    selected=[
        (near,far) for near,far in pairs
        if int(near.split("_",1)[0][1:])%ANCHOR_MODULUS==0
    ]
    selected=selected[:MAX_SAMPLE_POINTS//3]
    pointsets=[]
    for offset in (20,35):
        new,pairlist=near_points_from_140m(originals,selected,offset)
        pointsets.extend((key,xy,f"rail_{offset}m")
                         for key,xy in new.items()
                         if key.endswith(f"_near{offset}m"))
    for near,far in selected:
        pointsets.append((far,originals[far],"control_1600m"))
    assert len(pointsets)==3*len(selected)
    return pointsets,hashlib.sha256(source).hexdigest(),len(selected)


def run():
    sampled,n05_hash,n_pairs=choose_points()
    convert=Transformer.from_crs(32654,4326,always_xy=True)
    tiles=defaultdict(set)
    for key,(east,north),group in sampled:
        lon,lat=convert.transform(east,north)
        coord=world_tile(lon,lat)
        tiles[coord].add(group)
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-historical-GSI-aerial-LRT-source/1.4"})
    report={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "official_GSI_annual_photographs":GSI_PAGE,
        "GSI_tile_template":TEMPLATE,
        "years":YEARS,
        "tile_zoom":ZOOM,
        "n_fixed_geo_pairs_sampled":n_pairs,
        "n_sampled_location_points":len(sampled),
        "n_distinct_tile_addresses_per_year":len(tiles),
        "N05_historical_source_sha256":n05_hash,
        "N05_noncommercial_geometry_not_published":True,
        "result_by_year":{},"errors":{},
        "interpretation":"An extant year-named orthophoto tile is a candidate for manually reviewing 2017 road/track land cover, not proof that the photo was actually taken in that year at every underlying pixel.",
        "license":"GSI survey mapping images must be attributed under official GSI terms; no raw images redistributed."
    }
    # Source-control: if identical coordinates have a valid modern GSI
    # seamlessphoto JPEG, the all-404 year-photo result reflects absent
    # YEAR-LABELED tiles, not bad longitude/latitude XYZ or internet access.
    modern=Counter()
    for z,x,y in sorted(tiles):
        try:
            response=client.get(
                SEAMLESS.format(z=z,x=x,y=y),timeout=(12,30)
            )
            if response.status_code==404:
                label="http404_no_modern_photo"
            else:
                response.raise_for_status()
                label=image_is_nonblank(response.content,expected_format="JPEG")
            modern[label]+=1
        except (requests.RequestException,ValueError,OSError):
            modern["modern_source_error"]+=1
    report["modern_GSI_seamlessphoto_positive_control"]={
        "source":SEAMLESS,
        "tile_status_counts":dict(modern),
        "n_nonblank_modern_tiles":modern["photo_pixels_present"],
        "not_historical_cannot_establish_2017_road":True
    }
    print(json.dumps({"modern_photo_source_control":dict(modern)}),flush=True)

    for year in YEARS:
        counts=Counter()
        by_group=defaultdict(set)
        for i,(z,x,y) in enumerate(sorted(tiles)):
            url=TEMPLATE.format(year=year,z=z,x=x,y=y)
            try:
                r=client.get(url,timeout=(12,30))
                if r.status_code==404:
                    state="http404_no_tile"
                else:
                    r.raise_for_status()
                    state=image_is_nonblank(r.content)
            except (requests.RequestException,ValueError,OSError) as exc:
                state="source_error"
                if len(report["errors"])<10:
                    report["errors"][f"{year}_item_{i}"]=f"{type(exc).__name__}: {str(exc)[:140]}"
            counts[state]+=1
            for group in tiles[z,x,y]:
                by_group[group].add((z,x,y,state))
        report["result_by_year"][str(year)]={
            "n_source_tile_addresses":len(tiles),
            "tile_status_counts":dict(counts),
            "n_nonblank_photo_tiles":counts["photo_pixels_present"],
            "sample_group_tile_coverage_counts":{
                group:dict(Counter(state for zz,xx,yy,state in entries))
                for group,entries in by_group.items()
            }
        }
        print(json.dumps({
            "year":year,"tile_addresses":len(tiles),
            "status_counts":dict(counts)
        },ensure_ascii=False),flush=True)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/lrt_gsi_historical_aerial_tile_coverage_v14.json")
    a=p.parse_args()
    report=run()
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    if not report["result_by_year"]:
        raise SystemExit("No annual GSI photo inventory succeeded")


if __name__=="__main__":
    main()

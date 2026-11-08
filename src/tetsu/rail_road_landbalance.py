"""Independent MLIT N13 road distance and preconstruction S2 balance audit.

Historical N05 geometry is transient/noncommercial, 2024 N13 road GIS is CC-BY
4.0, Sentinel COG public. A road distance measured in 2024 is a
POST-CONSTRUCTION COVARIATE, not proof a road existed in 2017.

All sample points within a single urban LRT remain spatially correlated.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import math
import shutil
import tempfile
from pathlib import Path
from zipfile import ZipFile
from urllib.parse import urlparse
import requests
import shapefile
from pyproj import Transformer
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree

from typing import Mapping, Iterable

N13_URL = "https://nlftp.mlit.go.jp/ksj/gml/data/N13/N13-24/N13-24_{mesh}_SHP.zip"
N13_META = "https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-2024.html"
MESHES = ("5439","5440")  # official JIS grid tiles spanning lon 140°
BOX_LONLAT = (139.80,36.45,140.15,36.68)
EPSG = 32654
MAX_BYTES_PER_ZIP = 130_000_000
ROAD_DISTANCE_MAX_FOR_BALANCE_M = 100.0
ROAD_DISTANCE_DIFF_50M = 50.0
ROAD_DISTANCE_DIFF_100M = 100.0
PRE_NDVI_MAX_DIFF = 0.15  # QC-only, not outcome-optimizing fit


def verify_url(url: str):
    parts=urlparse(url)
    if (parts.scheme!="https" or parts.hostname!="nlftp.mlit.go.jp"
        or not parts.path.startswith("/ksj/gml/data/N13/N13-24/")
        or not parts.path.endswith("_SHP.zip")):
        raise ValueError("Untrusted N13 download URL")


def stream_zip(mesh, dest):
    if mesh not in MESHES:
        raise ValueError("Mesh not approved for LRT corridor")
    url=N13_URL.format(mesh=mesh)
    verify_url(url)
    client=requests.Session()
    client.headers.update({"User-Agent":"tetsu-jp-LRT-road-balance/1.3"})
    h=hashlib.sha256()
    n=0
    with client.get(url,stream=True,timeout=(25,180)) as resp:
        resp.raise_for_status()
        with open(dest,"wb") as output:
            for chunk in resp.iter_content(chunk_size=1024*1024):
                if not chunk:
                    continue
                n+=len(chunk)
                if n>MAX_BYTES_PER_ZIP:
                    raise ValueError("N13 mesh unexpectedly large")
                output.write(chunk)
                h.update(chunk)
    if not ZipFile:
        raise ValueError("ZIP reader unavailable")
    return {"mesh":mesh,"sha256":h.hexdigest(),"bytes":n,"url":url}


def projected_road_segments_from_zip(source_zip: Path, bounds=BOX_LONLAT):
    """Stream SHP/SHX locally; avoid loading nationwide ~million road geometries."""
    with ZipFile(source_zip) as z:
        shp_files=[n for n in z.namelist() if n.lower().endswith(".shp")]
        if not shp_files:
            raise ValueError("Missing N13 SHP file")
        shp_files.sort(key=lambda name:(
            0 if "utf-8" in name.lower() else 1,len(name)
        ))
        shp=shp_files[0]
        stem=shp[:-4]
        by_lower={s.lower():s for s in z.namelist()}
        if (stem+".shx").lower() not in by_lower:
            raise ValueError("Missing N13 SHX companion")
        with tempfile.TemporaryDirectory() as tmp:
            for ext in (".shp",".shx"):
                member=by_lower.get((stem+ext).lower())
                if member is None:
                    raise ValueError("Missing N13 source companion")
                with z.open(member) as source, open(Path(tmp)/("road"+ext),"wb") as dest:
                    shutil.copyfileobj(source,dest,length=1024*1024)
            reader=shapefile.Reader(
                shp=str(Path(tmp)/"road.shp"),
                shx=str(Path(tmp)/"road.shx")
            )
            if reader.shapeType not in (3,13,23):
                raise ValueError("N13 geometry not a polyline")
            total=len(reader)
            t=Transformer.from_crs(4326,EPSG,always_xy=True)
            segments=[]
            screened=0
            for shape in reader.iterShapes(bbox=bounds):
                screened+=1
                offsets=list(shape.parts)+[len(shape.points)]
                for start,stop in zip(offsets,offsets[1:]):
                    if stop-start<2:
                        continue
                    xy=[t.transform(lon,lat) for lon,lat in shape.points[start:stop]]
                    line=LineString(xy)
                    if line.length>=1:
                        segments.append(line)
            reader.close()
    return segments,{"all_source_lines":total,"lines_in_source_bbox":screened,
                     "projected_parts_in_bbox":len(segments)}


def build_n13_road_index(*, download_dir: Path):
    line_parts=[]
    sources=[]
    for mesh in MESHES:
        target=download_dir/f"N13-{mesh}.zip"
        source=stream_zip(mesh,target)
        lines,counts=projected_road_segments_from_zip(target)
        # Delete after parsing; do not redistribute N13 line features.
        target.unlink(missing_ok=True)
        if not lines:
            raise ValueError(f"No local N13 roads in mesh {mesh}")
        line_parts.extend(lines)
        sources.append({**source,**counts})
    tree=STRtree(line_parts)
    def nearest(point):
        p=Point(*point)
        idx=int(tree.nearest(p))
        return float(p.distance(line_parts[idx]))
    return nearest,{"metadata_url":N13_META,
                    "source_meshes":sources,
                    "n_line_parts_in_spatial_index":len(line_parts),
                    "road_is_2024_not_2017_preconstruction":True,
                    "JGD2011_geographic_to_UTM_EPSG":EPSG}


def pair_balance_flags(road_near:float,road_far:float,
                       scl_near:int|None,scl_far:int|None,
                       pre_ndvi_near:float|None,pre_ndvi_far:float|None):
    """Predefined *descriptive* pair balance gates. No p-values."""
    if any(not math.isfinite(float(x)) or float(x)<0 for x in
           (road_near,road_far)):
        raise ValueError("Invalid road distance")
    both_road100=(road_near<=ROAD_DISTANCE_MAX_FOR_BALANCE_M and
                  road_far<=ROAD_DISTANCE_MAX_FOR_BALANCE_M)
    road50=both_road100 and abs(road_near-road_far)<=ROAD_DISTANCE_DIFF_50M
    road100=both_road100 and abs(road_near-road_far)<=ROAD_DISTANCE_DIFF_100M
    both_clear=(scl_near in (4,5) and scl_far in (4,5))
    both_2017veg=(scl_near==4 and scl_far==4)
    same_scl=both_clear and scl_near==scl_far
    pre_similar=(
        pre_ndvi_near is not None and pre_ndvi_far is not None and
        math.isfinite(pre_ndvi_near) and math.isfinite(pre_ndvi_far) and
        abs(pre_ndvi_near-pre_ndvi_far)<=PRE_NDVI_MAX_DIFF
    )
    return {
        "near_road_le100_and_far_road_le100":both_road100,
        "road_abs_diff_le50_both_le100":road50,
        "road_abs_diff_le100_both_le100":road100,
        "2017_both_clear_scl4or5":both_clear,
        "2017_same_scl_class":same_scl,
        "2017_both_vegetation_scl4":both_2017veg,
        "2017_abs_NDVI_gap_le015":pre_similar,
        "road50_and_2017_same_class":road50 and same_scl,
        "road50_sameclass_preNDVI_similar":
            road50 and same_scl and pre_similar,
        # Additional conservative posthoc ecological screening: 2017
        # vegetation must be TRUE on BOTH sides, not merely both bare.
        "road50_both_2017_vegetation":road50 and both_2017veg,
        "road50_both2017veg_preNDVI_similar":
            road50 and both_2017veg and pre_similar,
    }


def count_balance(flags_by_pair: Mapping):
    cats=sorted({key for flags in flags_by_pair.values() for key in flags})
    return {"n_all_geo_pairs":len(flags_by_pair),
            "n_pairs_by_gate":{
                k:sum(bool(flags.get(k)) for flags in flags_by_pair.values())
                for k in cats},
            "cannot_estimate_causal_railway_effect":True}


def blocked_summary(pairs,vals,flags_by_pair,year,gate="all"):
    """Descriptive pair-of-fixed-pixel NDVI double difference, 4 corridor blocks.

    No pseudoreplicate p-values; 48 anchors create at most 4 longitudinal
    blocks in ONE corridor. Gates requiring 2017 true road exposure cannot
    be inferred from 2024 N13.
    """
    import statistics
    block=Counter()
    diffs=[]
    by_block={}
    for a,b in pairs:
        key=(a,b)
        if gate!="all" and not flags_by_pair.get(key,{}).get(gate,False):
            continue
        v=vals.get(2017,{});w=vals.get(year,{})
        z=[v.get(a),v.get(b),w.get(a),w.get(b)]
        if any(x is None for x in z):
            continue
        d=(z[2]-z[0])-(z[3]-z[1])
        if not math.isfinite(d):
            continue
        # Fixed spatial block by indexed anchor along corridor.
        anchor=int(a.split("_",1)[0][1:])
        block_id=anchor//12 # assumes 48 fixed anchor grid, not independent lines
        by_block.setdefault(block_id,[]).append(d)
        diffs.append(d)
    def sm(values):
        return {"n":len(values),
                "median":round(statistics.median(values),6) if values else None,
                "mean":round(statistics.mean(values),6) if values else None}
    return {
        "year":year,
        "pair_gate":gate,
        "paired_complete":sm(diffs),
        "n_longitudinal_blocks_with_data":len(by_block),
        "block_medians":[{"block":i,**sm(v)} for i,v in sorted(by_block.items())],
        "single_corridor_not_independent_treatments":True,
    }

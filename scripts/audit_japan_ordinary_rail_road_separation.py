#!/usr/bin/env python3
"""Find conventional, potentially road-separated Japanese railway sites.

REAL official public-source screening in THREE prespecified spatially
independent regions: Miyagi/Sendai, Niigata, and Nagano. Railway N02-2025
(CC BY 4.0) includes passenger rail, not freight-only. Road N13-2024
(CC BY 4.0) gives contemporary road centreline, not road presence in 2017.
All source GIS data transient. No railway grass, seed, flower/fruit/spore
ROI, cutting/wind exposure, site vegetation or treatment effect inferred.
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
import shutil
import tempfile
import zipfile

import requests
import shapefile
from pyproj import Transformer
from shapely.geometry import LineString,Point,box
from shapely.strtree import STRtree

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.ordinary_rail_corridors import (
    MIN_ROUTE_WITHIN_BBOX_M, is_ordinary_rail,
    route_normal_pairs, safe_lines_from_shape,nearest_distance,
    classify_candidate,corridor_aggregate
)

N02="https://nlftp.mlit.go.jp/ksj/gml/data/N02/N02-25/N02-25_GML.zip"
N13="https://nlftp.mlit.go.jp/ksj/gml/data/N13/N13-24/N13-24_{mesh}_SHP.zip"
OFFICIAL_RAIL="https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N02-2025.html"
OFFICIAL_ROAD="https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-2024.html"

# Study rectangles are deliberately INSIDE their 1-degree x 40-minute
# first-level N13 mesh; the 1.75km sampling guard avoids mesh edges.
# This is a *source feasibility frame* in JP cities, not the whole prefecture.
REGIONS={
    "sendai_miyagi":{
        "mesh":"5740","utm":32654,
        "bbox":(140.53,38.12,140.96,38.57)
    },
    "niigata":{
        "mesh":"5639","utm":32654,
        "bbox":(139.04,37.56,139.73,37.95)
    },
    "nagano":{
        "mesh":"5438","utm":32654,
        "bbox":(138.04,36.24,138.73,36.61)
    }
}
MAX_ZIP_ROAD=100_000_000
MAX_ZIP_RAIL=50_000_000
MIN_ROAD_SEGMENTS=200
USER_AGENT="tetsu-independent-ordinary-rail-road-verge-source/1.5"

def download_source(url,max_size,to_file=None):
    client=requests.Session()
    client.headers.update({"User-Agent":USER_AGENT})
    h=hashlib.sha256()
    n=0
    if to_file is None:
        chunks=[]
    with client.get(url,stream=True,timeout=(25,170)) as response:
        response.raise_for_status()
        if to_file is None:
            for part in response.iter_content(1024*1024):
                if not part:continue
                n+=len(part)
                if n>max_size:raise ValueError("official source exceeded cap")
                chunks.append(part);h.update(part)
        else:
            with open(to_file,"wb") as output:
                for part in response.iter_content(1024*1024):
                    if not part:continue
                    n+=len(part)
                    if n>max_size:raise ValueError("official source exceeded cap")
                    output.write(part);h.update(part)
    return (b"".join(chunks) if to_file is None else None),{"source_sha256":h.hexdigest(),"bytes":n,"source_url":url}


def choose_n02_member(archive,needle):
    candidates=[n for n in archive.namelist()
                if n.lower().endswith(needle.lower())]
    if not candidates:raise ValueError(f"N02 missing {needle}")
    candidates.sort(key=lambda name:(0 if "utf-8" in name.lower() else 1,len(name)))
    return candidates[0]


def read_n02_region(archive_raw,region):
    cfg=REGIONS[region]
    lon0,lat0,lon1,lat1=cfg["bbox"]
    trans=Transformer.from_crs(4326,cfg["utm"],always_xy=True)
    x0,y0=trans.transform(lon0,lat0)
    x1,y1=trans.transform(lon1,lat1)
    bbox_m=(min(x0,x1),min(y0,y1),max(x0,x1),max(y0,y1))
    extent=box(*bbox_m)
    segments=defaultdict(list)
    all_lines=[]
    station_lines=[]
    rail_type=Counter()
    with zipfile.ZipFile(io.BytesIO(archive_raw)) as archive:
        mem={n.lower():n for n in archive.namelist()}
        def reader(member):
            stem=member[:-4]
            def b(ext):
                name=mem.get((stem+ext).lower())
                if not name:raise ValueError(f"N02 incomplete {stem+ext}")
                return io.BytesIO(archive.read(name))
            return shapefile.Reader(shp=b(".shp"),shx=b(".shx"),
                dbf=b(".dbf"),encoding="utf-8",encodingErrors="replace")
        shp=choose_n02_member(archive,"_RailroadSection.shp")
        rail=reader(shp)
        names=[f[0] for f in rail.fields[1:]]
        needed={"N02_001","N02_002","N02_003","N02_004"}
        if not needed.issubset(names):
            raise ValueError(f"N02 field schema drift {needed-set(names)}")
        for feature in rail.iterShapeRecords(bbox=(
            lon0-.03,lat0-.03,lon1+.03,lat1+.03)):
            record=dict(zip(names,feature.record))
            rails=safe_lines_from_shape(feature.shape,trans)
            for line in rails:
                clipped=line.intersection(extent)
                if clipped.is_empty:continue
                parts=[clipped] if clipped.geom_type=="LineString" else (
                    list(clipped.geoms) if clipped.geom_type=="MultiLineString" else []
                )
                for piece in parts:
                    if piece.length<500:continue
                    all_lines.append(piece)
                    if is_ordinary_rail(record):
                        key=(str(record["N02_003"]).strip(),
                             str(record["N02_004"]).strip())
                        segments[key].append(piece)
                        rail_type[str(record["N02_001"]).strip()]+=1
        rail.close()
        stations=reader(choose_n02_member(archive,"_Station.shp"))
        for feature in stations.iterShapes(bbox=(
            lon0-.02,lat0-.02,lon1+.02,lat1+.02)):
            station_lines.extend(safe_lines_from_shape_station(feature,trans))
        stations.close()
    if not all_lines or not segments:
        raise ValueError(f"N02 no conventional rail in study area {region}")
    if not station_lines:
        raise ValueError(f"N02 lacks station geometry {region}")
    return segments,all_lines,station_lines,bbox_m,{
        "n_N02_rail_segments_clipped":len(all_lines),
        "n_N02_ordinary_route_IDs":len(segments),
        "n_N02_station_line_geometries":len(station_lines),
        "ordinary_rail_type_row_counts":dict(rail_type)
    }


def safe_lines_from_shape_station(src_shape,trans):
    """Stations are source *line* geometries, often <500m, not points."""
    if src_shape.shapeType not in (3,13,23):
        return []
    out=[]
    parts=list(src_shape.parts)+[len(src_shape.points)]
    for a,b in zip(parts,parts[1:]):
        if b-a<2:continue
        pts=[trans.transform(lon,lat) for lon,lat in src_shape.points[a:b]]
        line=LineString(pts)
        if line.length>0:out.append(line)
    return out


def n13_region_index(region,workdir):
    cfg=REGIONS[region]
    mesh=cfg["mesh"]
    url=N13.format(mesh=mesh)
    filename=Path(workdir)/f"n13_{mesh}.zip"
    _,source=download_source(url,MAX_ZIP_ROAD,to_file=filename)
    lines=[]
    with zipfile.ZipFile(filename) as z:
        shp=[n for n in z.namelist() if n.lower().endswith(".shp")]
        if not shp:raise ValueError("N13 missing road SHP")
        shp.sort(key=lambda n:(0 if "utf-8" in n.lower() else 1,len(n)))
        name=shp[0]
        stem=name[:-4]
        names={n.lower():n for n in z.namelist()}
        # Stream SHP/SHX to temporary local files, not to repository.
        with tempfile.TemporaryDirectory(dir=workdir) as local:
            for ext in (".shp",".shx"):
                member=names.get((stem+ext).lower())
                if not member:raise ValueError("N13 source sidecar absent")
                with z.open(member) as src, open(Path(local)/("road"+ext),"wb") as dest:
                    shutil.copyfileobj(src,dest,length=1024*1024)
            reader=shapefile.Reader(
                shp=str(Path(local)/"road.shp"),
                shx=str(Path(local)/"road.shx")
            )
            if reader.shapeType not in (3,13,23):
                raise ValueError("N13 SHP is not polyline")
            raw_total=len(reader)
            trans=Transformer.from_crs(4326,cfg["utm"],always_xy=True)
            region_bbox=cfg["bbox"]
            # Extend by 0.02 degrees to include roads nearest an edge.
            target=(region_bbox[0]-.02,region_bbox[1]-.02,
                    region_bbox[2]+.02,region_bbox[3]+.02)
            features=0
            for shape in reader.iterShapes(bbox=target):
                features+=1
                parts=list(shape.parts)+[len(shape.points)]
                for a,b in zip(parts,parts[1:]):
                    if b-a<2:continue
                    ls=LineString([
                        trans.transform(lon,lat) for lon,lat in shape.points[a:b]
                    ])
                    if ls.length>=1:lines.append(ls)
            reader.close()
    filename.unlink(missing_ok=True)
    if len(lines)<MIN_ROAD_SEGMENTS:
        raise ValueError(f"N13 region {region} has only {len(lines)} line parts")
    return lines,{**source,
                  "n_total_road_source_shapes":raw_total,
                  "n_road_features_in_study_bbox":features,
                  "n_local_projected_road_line_parts":len(lines)}


def private_scope_includes(flags, scope):
    if scope=="strict":
        return flags["strict_near_road75_far_road50_far_rail1000_station300"]
    if scope=="all_fixed_geo":
        return True
    raise ValueError("Unexpected private candidate scope")


def analyze_region(region,rail_archive,*,private_site_positions=False,
                   private_candidate_scope="strict"):
    groups,rail_all,stations,bbox_m,rail_info=read_n02_region(
        rail_archive,region
    )
    with tempfile.TemporaryDirectory() as scratch:
        road,road_info=n13_region_index(region,scratch)
        road_tree=STRtree(road)
        station_tree=STRtree(stations)
        rail_tree=STRtree(rail_all)
        total=Counter()
        routes=[]
        private_candidates=[]  # in-memory only, NEVER publish coordinates
        for (name,company),segments in sorted(groups.items()):
            length_km=round(sum(x.length for x in segments)/1000.,3)
            if length_km*1000<MIN_ROUTE_WITHIN_BBOX_M:
                continue
            side_pairs,n_anchors=route_normal_pairs(
                segments,bbox_m=bbox_m,n_cap=42
            )
            measurements=[]
            for near,far in side_pairs:
                road_near=nearest_distance(road_tree,road,near)
                road_far=nearest_distance(road_tree,road,far)
                far_rail=nearest_distance(rail_tree,rail_all,far)
                nearest_station=nearest_distance(station_tree,stations,near)
                if None in (road_near,road_far,far_rail,nearest_station):
                    raise ValueError("missing required official line geometry")
                flags=classify_candidate(
                    near,far,near_road=road_near,far_road=road_far,
                    nearest_any_rail_far=far_rail,near_station=nearest_station
                )
                measurements.append({
                    "near_road_m":road_near,
                    "far_road_m":road_far,
                    "rail_far_m":far_rail,
                    "near_station_m":nearest_station,
                    "eligibility":flags
                })
                if private_site_positions and private_scope_includes(
                    flags, private_candidate_scope
                ):
                    # Only for transient source-matching inside an
                    # authorized source audit: never serialize coordinates.
                    private_candidates.append({
                        "route_name":name,"operator_name":company,
                        "near_xy":(near.x,near.y),
                        "road_only_xy":(far.x,far.y),
                        "near_road_m":road_near,
                        "far_road_m":road_far,
                        "far_rail_m":far_rail,
                        "near_station_m":nearest_station,
                        "source_gates":flags
                    })
            summary=corridor_aggregate(measurements)
            summary["rail_length_clipped_km"]=length_km
            summary["n_anchors_preselected_without_outcomes"]=n_anchors
            summary["route_name"]=name
            summary["operator_name"]=company
            summary["source"]="N02-2025 ordinary railway line, N13-2024 roads"
            summary["needs_ground_check_for_at_grade_rail_and_actual_verge"]=True
            routes.append(summary)
            total.update({
                "geometric_pairs":summary["n_geographically_sampled_rail_sides"],
                "n_strict_pair_opportunities":
                    summary["n_eligible_by_source_gate"].get(
                        "strict_near_road75_far_road50_far_rail1000_station300",0
                    )
            })
    routes.sort(key=lambda x:(
        -x["n_eligible_by_source_gate"].get(
            "strict_near_road75_far_road50_far_rail1000_station300",0
        ),
        -x["n_eligible_by_source_gate"].get("near_road_ge75",0),
        -x["rail_length_clipped_km"],
        x["route_name"],x["operator_name"]
    ))
    summary={
        "region":region,
        "source_rail":rail_info,
        "source_roads":road_info,
        "n_ordinary_route_groups_over_2p5km":len(routes),
        "total_rail_side_source_opportunities":total["geometric_pairs"],
        "total_strict_road_control_opportunities":
            total["n_strict_pair_opportunities"],
        "n_route_groups_with_any_road_isolated_and_matched_control":sum(
            v["n_eligible_by_source_gate"].get(
                "strict_near_road75_far_road50_far_rail1000_station300",0
            )>0 for v in routes
        ),
        "ranked_route_candidates_without_plant_outcomes":routes[:30],
        "rank_rule":"number of fixed geographically eligible rail-versus-road control samples, NOT plant cover, NDVI, species frequency or inferred biological effect",
    }
    # When requested by another authorized source audit, coordinates stay
    # in its transient runtime; public JSON must contain only 'summary'.
    return (summary,private_candidates) if private_site_positions else summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="results/audit/japan_ordinary_rail_road_separation_v15.json")
    p=parser.parse_args()
    raw,source=download_source(N02,MAX_ZIP_RAIL)
    if raw is None or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise SystemExit("N02 required source not a valid zip")
    report={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "study_type":"independently sourced geography only; NOT railway biodiversity or seed transport effects",
        "rail_source":{**source,"reference_year":2025,"license":"CC BY 4.0"},
        "road_source_page":OFFICIAL_ROAD,
        "rail_source_page":OFFICIAL_RAIL,
        "fixed_regions":REGIONS,
        "rail_type_codes":"11=JR ordinary, 12=other ordinary rail; 13–17 incl funicular/monorail and 21–25 tram excluded; high-speed operator class1 excluded",
        "min_nominal_line_length_km":2.5,
        "sample_per_route":"fixed 1.2km anchors along clipped N02 polyline; 35m normal near line vs 1600m perpendicular control",
        "control_source":"N13-2024 road centreline <=50m and any N02 passenger rail >=1000m",
        "source_thresholds":{
            "near_rail_road_separation_m":[50,75,100],
            "near_station_distance_m":300,
            "far_min_rail_distance_m":1000,
            "far_max_road_distance_m":50
        },
        "reports":{},
        "errors":{},
        "critical_warnings":[
            "Geometric centreline 35m offset does not establish vegetation, true railway land ownership, embankment/ballast substrate, exposed train gust or railway management.",
            "Some N02 ordinary rail segments could run through tunnels, over bridges or viaducts; NO at-grade rail source attribute is available here. Each candidate needs high-resolution source/human ground verification before biology.",
            "N13 2024 road centreline is a contemporaneous proxy; road pavement width, service traffic, past road changes and urban use not controlled.",
            "N02 2025 passenger rail omits freight-only railways; no active train freight traffic measured.",
            "More than one normal transect along the SAME route does not create independent railway treatments or license causal regression.",
            "Near-road >=75m and far-road <=50m are only GEOMETRIC source feasibility, not biological presence/absence.",
            "Route names/companies are not independent true track corridors without track-collision, same operator and topology verification.",
            "Biological NDVI, botany, ecological selection, phenology, photo ROI and railway propagule dispersal remain NOT measured by this stage."
        ]
    }
    for name in REGIONS:
        try:
            outcome=analyze_region(name,raw)
            report["reports"][name]=outcome
            print(json.dumps({
                "region":name,
                "n_routes":outcome["n_ordinary_route_groups_over_2p5km"],
                "road_parts":outcome["source_roads"]["n_local_projected_road_line_parts"],
                "fixed_side_pairs":outcome["total_rail_side_source_opportunities"],
                "strict_road_only_vs_rail_opportunities":outcome["total_strict_road_control_opportunities"],
                "n_routes_with_any":outcome["n_route_groups_with_any_road_isolated_and_matched_control"],
                "ranked_top3":[{
                    "route":v["route_name"],"operator":v["operator_name"],
                    "length_km":v["rail_length_clipped_km"],
                    "strict_source_pairs":v["n_eligible_by_source_gate"].get(
                        "strict_near_road75_far_road50_far_rail1000_station300",0
                    ),
                    "near_road_ge75":v["n_eligible_by_source_gate"].get("near_road_ge75",0)
                } for v in outcome["ranked_route_candidates_without_plant_outcomes"][:3]]
            },ensure_ascii=False),flush=True)
        except Exception as exc:
            report["errors"][name]=f"{type(exc).__name__}: {str(exc)[:500]}"
            print(json.dumps({"region":name,"source_error":report["errors"][name]}),flush=True)
    dest=Path(p.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    if len(report["reports"])<2:
        raise SystemExit("Fewer than 2 geographic regions verified, no cross-region candidate inference")


if __name__=="__main__":
    main()

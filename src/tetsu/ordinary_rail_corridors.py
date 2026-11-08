"""Source-independent sampling/eligibility rules for JP ordinary railway verges.

Official N02-2025 railway line (2025-12-31) and N13-2024 road line
(2024-09). Road proximity alone does NOT prove trackside plants,
railway operations, train wind, mowing or actual freight. N02 omits
freight-only rail lines and N13 does not date road opening.
"""
from __future__ import annotations
import math
from collections import Counter
from shapely.geometry import LineString, Point, box
from shapely.ops import unary_union

ORDINARY_TYPES={"11","12"}     # rail, NOT streetcar, subway, monorail etc.
EXCLUDE_ROUTE_TOKENS=("新幹線","地下鉄","リニア","メトロ","地下鉄線")
MIN_SECTION_M=500.
SAMPLE_INTERVAL_M=1200.
MIN_ROUTE_WITHIN_BBOX_M=2500.
NEAR_OFFSET_M=35.
ROAD_SEPARATED_M=(50.,75.,100.)
CONTROL_OFFSET_M=1600.
FAR_MIN_RAIL_M=1000.
FAR_MAX_ROAD_M=50.
MIN_STATION_DISTANCE_M=300.
MIN_ACCEPTED_SPACING_M=900.


def is_ordinary_rail(record):
    code=str(record.get("N02_001") or "").strip()
    if code not in ORDINARY_TYPES:
        return False
    operator_type=str(record.get("N02_002") or "").strip()
    route=str(record.get("N02_003") or "").strip()
    company=str(record.get("N02_004") or "").strip()
    if not route or not company:
        return False
    # High-speed rail has separate operator category 1. Ordinary subway
    # lines may be category '12' so exclude independently by labels.
    if operator_type=="1" or any(
        token in route or token in company for token in EXCLUDE_ROUTE_TOKENS
    ):
        return False
    return True


def safe_lines_from_shape(src_shape,transformer):
    """Part-wise projected railway geometry, including MultipartLineString."""
    if src_shape.shapeType not in (3,13,23):
        return []
    vertices=src_shape.points
    boundaries=list(src_shape.parts)+[len(vertices)]
    out=[]
    for a,b in zip(boundaries,boundaries[1:]):
        if b-a<2:
            continue
        points=[transformer.transform(lon,lat) for lon,lat in vertices[a:b]]
        line=LineString(points)
        if line.is_valid and line.length>=MIN_SECTION_M:
            out.append(line)
    return out


def route_normal_pairs(lines, *, bbox_m, n_cap=50):
    """Deterministic 1.2km route anchors with both normal directions.

    No road, plant, photograph, NDVI or visitation data used to choose
    anchor locations. Excludes outside region with 1.7km guard so
    comparator coordinates remain inside its 2024 road mesh.
    """
    guard=CONTROL_OFFSET_M+150.
    inner=box(
        bbox_m[0]+guard,bbox_m[1]+guard,
        bbox_m[2]-guard,bbox_m[3]-guard
    )
    if inner.is_empty:
        raise ValueError("region too small for controls")
    available=[line for line in lines if line.length>=MIN_SECTION_M]
    available.sort(key=lambda line:(
        -round(line.length,4),
        round(line.bounds[0],1),round(line.bounds[1],1)
    ))
    anchors=[]
    opportunities=[]
    for line in available:
        n=int(line.length/SAMPLE_INTERVAL_M)
        for i in range(n):
            d=(i+.5)*line.length/n
            center=line.interpolate(d)
            if not inner.covers(center):
                continue
            if any(center.distance(existing)<MIN_ACCEPTED_SPACING_M
                   for existing in anchors):
                continue
            earlier=line.interpolate(max(0.,d-30))
            later=line.interpolate(min(line.length,d+30))
            dx,dy=later.x-earlier.x,later.y-earlier.y
            size=math.hypot(dx,dy)
            if size<1:
                continue
            nx,ny=-dy/size,dx/size
            anchors.append(center)
            for side in (-1,1):
                near=Point(center.x+side*nx*NEAR_OFFSET_M,
                           center.y+side*ny*NEAR_OFFSET_M)
                far=Point(center.x+side*nx*CONTROL_OFFSET_M,
                          center.y+side*ny*CONTROL_OFFSET_M)
                if inner.covers(near) and inner.covers(far):
                    opportunities.append((near,far))
            if len(anchors)>=n_cap:
                break
        if len(anchors)>=n_cap:
            break
    return opportunities,len(anchors)


def nearest_distance(tree, geometries, point):
    if tree is None or not geometries:
        return None
    index=int(tree.nearest(point))
    return float(point.distance(geometries[index]))


def classify_candidate(near,far,*,near_road,far_road,
                       nearest_any_rail_far,near_station):
    """Noncausal source eligibility, each type of evidence kept distinct."""
    if not all(math.isfinite(v) and v>=0 for v in (
        near_road,far_road,nearest_any_rail_far,near_station
    )):
        raise ValueError("invalid metric distances")
    return {
        "near_road_ge50":near_road>=50.,
        "near_road_ge75":near_road>=75.,
        "near_road_ge100":near_road>=100.,
        "near_station_ge300":near_station>=MIN_STATION_DISTANCE_M,
        "far_rail_ge1000":nearest_any_rail_far>=FAR_MIN_RAIL_M,
        "far_road_le50":far_road<=FAR_MAX_ROAD_M,
        "strict_near_road75_far_road50_far_rail1000_station300":(
            near_road>=75. and far_road<=50. and
            nearest_any_rail_far>=1000. and
            near_station>=300.
        ),
    }


def bin_distance_m(d):
    if d is None:return "unknown"
    if d<25:return "lt25m"
    if d<50:return "25_50m"
    if d<75:return "50_75m"
    if d<100:return "75_100m"
    return "ge100m"


def corridor_aggregate(rows):
    """Aggregate repeated spatial samples without pseudo-effect p-values."""
    out=Counter()
    near_road_bins=Counter()
    far_road_bins=Counter()
    for record in rows:
        for key,val in record["eligibility"].items():
            out[key]+=int(val)
        near_road_bins[bin_distance_m(record["near_road_m"])]+=1
        far_road_bins[bin_distance_m(record["far_road_m"])]+=1
    return {
        "n_geographically_sampled_rail_sides":len(rows),
        "n_eligible_by_source_gate":dict(out),
        "road_distance_near_bins":dict(near_road_bins),
        "road_distance_far_bins":dict(far_road_bins),
        "each_side_is_correlated_with_rail_segment":True,
        "none_of_source_geometries_establishes_vegetation_or_management":True,
    }

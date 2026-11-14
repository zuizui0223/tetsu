"""Rail right-of-way pixel eligibility: NOT proof of actual trackside vegetation.

The official N05 rail centreline supports perpendicular transects. MLIT N13
is 2024 road centreline only, and Sentinel-2 red/NIR 10m plus SCL 20m
represent mixed urban land cover. All metrics are source diagnostics.

No inference on train seed transport, mowing, wind, adaptation or plant species.
"""
from __future__ import annotations
from collections import Counter
import math
import statistics

OFFSETS=(20.,35.)
ROAD_BUFFER_M=25.
ROAD_MATCH_MAX_DIST_M=100.
ROAD_MATCH_MAX_ABS_DIFF_M=50.


def near_points_from_140m(points,pairs,near_offset):
    """Interpolate true perpendicular near point from fixed 140m/1600m pair.

    Since both existing points were computed on the SAME outward normal
    at 140m and 1600m, linear interpolation at a specified smaller
    offset is geometrically exact (up to floating-point rounding).
    This does not use satellite reflectance or plant observations.
    """
    if not 0<near_offset<140:
        raise ValueError("near_offset must be strictly 0–140m")
    out={}
    mapped=[]
    fraction=(near_offset-140.)/(1600.-140.)
    for old_near,far in pairs:
        if old_near not in points or far not in points:
            raise ValueError("missing source near/far geometry")
        a=points[old_near]
        b=points[far]
        nid=old_near.replace("_near",f"_near{int(near_offset)}m")
        if nid==old_near:
            raise ValueError("unknown original near-ID convention")
        out[nid]=(a[0]+fraction*(b[0]-a[0]),
                  a[1]+fraction*(b[1]-a[1]))
        out[far]=b
        mapped.append((nid,far))
    if len(out)!=2*len(mapped):
        raise ValueError("duplicate same-rail comparison geometry")
    return out,mapped


def class_group(scl):
    if scl==4:return "vegetation"
    if scl==5:return "not_vegetated"
    return "missing_or_unsafe"


def preconstruction_comparability(road_near,road_far,
                                 scl_near,scl_far,
                                 ndvi_near,ndvi_far):
    """Outcomes independent tests on 2017 class + 2024 road context.

    Availability of a 'non-road centerline' sample does NOT imply the
    2024 road surface is farther than road width or that track vegetation
    existed in 2017. The thresholds are source eligibility, posthoc.
    """
    if not all(math.isfinite(x) and x>=0 for x in (road_near,road_far)):
        raise ValueError("invalid current N13 road distance")
    if not all(x is None or -1<=x<=1 for x in (ndvi_near,ndvi_far)):
        raise ValueError("invalid NDVI")
    observed=(ndvi_near is not None and ndvi_far is not None and
              scl_near in (4,5) and scl_far in (4,5))
    both_vegetation=(observed and scl_near==4 and scl_far==4)
    similar_2017=(both_vegetation and
                  abs(ndvi_near-ndvi_far)<=.15)
    road_similar=(road_near<=ROAD_MATCH_MAX_DIST_M and
                  road_far<=ROAD_MATCH_MAX_DIST_M and
                  abs(road_near-road_far)<=ROAD_MATCH_MAX_ABS_DIFF_M)
    away_from_road=(road_near>=ROAD_BUFFER_M)
    return {
        "2017_scl_both_clear":observed,
        "2017_both_vegetation":both_vegetation,
        "2017_both_vegetation_near_far_NDVI_similar":similar_2017,
        "2024_near_road_lt25m":road_near<ROAD_BUFFER_M,
        "2024_near_road_ge25m":away_from_road,
        "2024_road_distances_matched":road_similar,
        "2017_both_vegetation_and_2024_near_not_road25":both_vegetation and away_from_road,
        "strict_veg_road_isolated_and_road_matched":(
            similar_2017 and away_from_road and road_similar
        )
    }


def analyze_pairs(pairs,year_data,pre_data,pre_scl,road_by_id):
    gates={}
    classes=Counter()
    road_n=[]
    road_f=[]
    for near,far in pairs:
        rn=road_by_id[near];rf=road_by_id[far]
        road_n.append(rn);road_f.append(rf)
        sn=pre_scl.get(near)
        sf=pre_scl.get(far)
        classes[(class_group(sn),class_group(sf))]+=1
        gates[(near,far)]=preconstruction_comparability(
            rn,rf,sn,sf,pre_data.get(near),pre_data.get(far)
        )
    details={
        "n_geographic_pairs":len(pairs),
        "2024_road_distance_near_median_m":round(statistics.median(road_n),2) if road_n else None,
        "2024_road_distance_far_median_m":round(statistics.median(road_f),2) if road_f else None,
        "pair_2017_SCL_category_counts":{
            f"{a}_near__{b}_control":v for (a,b),v in classes.items()
        },
        "eligibility_counts":{
            name:sum(row[name] for row in gates.values())
            for name in sorted(next(iter(gates.values()),{}))
        }
    }
    for year,values in year_data.items():
        for gate_name in ("all",
                          "2017_both_vegetation",
                          "2017_both_vegetation_and_2024_near_not_road25",
                          "strict_veg_road_isolated_and_road_matched"):
            contrasts=[]
            block=Counter()
            for near,far in pairs:
                if gate_name!="all" and not gates[(near,far)][gate_name]:
                    continue
                pre_a=pre_data.get(near);pre_b=pre_data.get(far)
                post_a=values.get(near);post_b=values.get(far)
                if any(v is None for v in (pre_a,pre_b,post_a,post_b)):
                    continue
                contrasts.append(
                    (post_a-pre_a)-(post_b-pre_b)
                )
                anchor=int(near.split("_",1)[0][1:])
                block[anchor//12]+=1
            details.setdefault("pixel_NDVI_contrasts",{}).setdefault(str(year),{})[gate_name]={
                "n_complete":len(contrasts),
                "median":round(statistics.median(contrasts),6) if contrasts else None,
                "mean":round(statistics.mean(contrasts),6) if contrasts else None,
                "n_contiguous_along_rail_blocks":len(block),
                "block_pair_counts":dict(sorted(block.items())),
                "one_railway_corridor_not_independent_experiments":True,
            }
    return details

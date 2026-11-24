#!/usr/bin/env python3
"""v368b REAL worldwide length-uniform physical railway point probability sample.

Read ALL 16 SHA-verified Overture 2026-11-18 rail-only GeoParquet shards,
census 1,536,226 original standard/narrow/unknown non-abandoned,
non-disused, non-construction candidate rail GERS segments. Sample
EXACTLY 2048 independent uniform SOURCE WGS84 GEODESIC railway-km
locations, WITH replacement, by 2-pass full original rail centreline
cumulative length, not 2048 arbitrary handpicked cities/stations.

This is an original *rail-km* point population, not operated train
service, named commercial lines, botanical plot vegetation, or native
seed dispersal. Full point GEOJSON is written only as an original CI
artifact for future source-local plant matching, not GitHub main.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,time
from collections import Counter
from pathlib import Path

from build_global_Overture_physical_connector_components_v361 import (
    validate_files,eligible_panels,SOURCE_RELEASE,EXPECTED_INVENTORY_SHA
)
from audit_global_GERS_connector_position_consistency_v363 import (
    true_WGS84_geodesic_at
)

SAMPLE_K=2048
SEED="tetsu-v368-20261124-rail-km-global-random-seed"
EXPECTED_CANDIDATE_SEGMENTS=1536226
EXPECTED_CANDIDATE_KM=2237043.481
EXPECTED_WORLD_SOURCE_PARQUET_N=16
MAX_ABSOLUTE_GLOBAL_CANDIDATE_CENSUS_KM_DIFFERENCE=.035
BATCH=20000

def fixed_sample_targets(total_metres):
    if not math.isfinite(total_metres) or total_metres<=1e9:
        raise ValueError("Original all-world rail source length not complete")
    rng=random.Random(hashlib.sha256(SEED.encode("utf8")).hexdigest())
    # With-replacement independent rail-km draws. No after-outcome
    # rejection, duplicates retained as distinct independent targets.
    values=[(total_metres*rng.random(),j) for j in range(SAMPLE_K)]
    values.sort()
    if len(values)!=SAMPLE_K or any(v<0 or v>=total_metres for v,_ in values):
        raise AssertionError("Original global cumulative rail km targets invalid")
    return values

def stream_source_rails(paths):
    import pyarrow.parquet as pq
    from shapely import wkb
    for path in paths:
        pf=pq.ParquetFile(path)
        required=("id","subtype","class","rail_flags","geometry")
        if not set(required)<=set(pf.schema_arrow.names):
            raise ValueError("Original rail-only source physical geometry/flags column missing")
        for batch in pf.iter_batches(batch_size=BATCH,columns=list(required),use_threads=True):
            for row in batch.to_pylist():
                if row["subtype"]!="rail":
                    raise ValueError("Original source global rail-only shard contains roads")
                yield row

def candidate_segment(row):
    return eligible_panels(row)["conventional_active_candidate_source_rails"]

def seg_m(row,geod):
    from shapely import wkb
    shape=wkb.loads(bytes(row["geometry"]))
    if shape.geom_type!="LineString" or shape.is_empty:
        raise ValueError("Original GERS line malformed in global rail-only source")
    length=abs(float(geod.geometry_length(shape)))
    if not math.isfinite(length) or length<0:
        raise ValueError("Original global rail geometry geodesic negative/nonfinite")
    return shape,length

def census_global_source(paths):
    from pyproj import Geod
    geod=Geod(ellps="WGS84")
    counts=Counter();distance=0.
    t=time.monotonic()
    for row in stream_source_rails(paths):
        counts["original_all_rail_subtype_segments_processed"]+=1
        if not candidate_segment(row):continue
        source_geom,L=seg_m(row,geod)
        counts["original_candidate_normal_active_map_rail_segments"]+=1
        counts["source_normal_rail_class_"+str(row["class"])]+=1
        distance+=L
        if counts["original_candidate_normal_active_map_rail_segments"]%200000==0:
            print(json.dumps({"global_original_candidate_rail_segments_sofar":
                counts["original_candidate_normal_active_map_rail_segments"],
                "original_source_geodesic_candidate_rail_km_sofar":round(distance/1000,3),
                "firstpass_seconds":round(time.monotonic()-t,1)}),flush=True)
    if counts["original_candidate_normal_active_map_rail_segments"]!=EXPECTED_CANDIDATE_SEGMENTS:
        raise ValueError("Full original source conventional rail n != v361 1536226")
    if counts["original_all_rail_subtype_segments_processed"]!=2054855:
        raise ValueError("Original all world 2,054,855 GERS rail segments not scanned")
    if abs(distance/1000-EXPECTED_CANDIDATE_KM)>MAX_ABSOLUTE_GLOBAL_CANDIDATE_CENSUS_KM_DIFFERENCE:
        raise ValueError("Original global candidate WGS84 rail km drift >35m from true v361")
    return distance,dict(counts)

def original_rail_flags_at_fraction(rules,at):
    """Return only true Overture rail flags applying at exact SOURCE at.

    Overture rail_flags are geometrically SCOPED [start,end] along
    original geodetic rail segment. Any-flag-on-segment is NOT proof
    of bridge/tunnel at an independently sampled line location.
    """
    result=set()
    for rule in rules or ():
        if not isinstance(rule,dict):continue
        rng=rule.get("between")
        if rng is not None:
            if not isinstance(rng,(list,tuple)) or len(rng)!=2:
                raise ValueError("Original source rail flag scope not length-two")
            a,b=map(float,rng)
            if not (0<=a<b<=1):
                raise ValueError("Original source rail flag WGS84 interval invalid")
            if not (a<=at<b or (at==1.0 and b==1.0)):
                continue
        result.update(rule.get("values") or ())
    return result


def select_source_global_rail_km(paths,distance_m):
    from pyproj import Geod
    geod=Geod(ellps="WGS84")
    target=fixed_sample_targets(distance_m)
    current=0;acc=0.0
    out=[];counts=Counter()
    seen=0
    start=time.monotonic()
    for row in stream_source_rails(paths):
        seen+=1
        if not candidate_segment(row):continue
        geom,L=seg_m(row,geod)
        lo,hi=acc,acc+L
        if L>0:
            while current<len(target) and target[current][0]<hi:
                t,idx=target[current]
                if t+1e-8<lo:
                    raise ValueError("Source cumulative rail target misplaced across original sorted shard")
                fraction=min(1.,max(0.,(t-lo)/L))
                lon,lat=true_WGS84_geodesic_at(geom,fraction,geod)
                if not(-180<=lon<=180 and -90<=lat<=90):
                    raise ValueError("Original global source rail point coordinate invalid")
                flagval=original_rail_flags_at_fraction(
                    row.get("rail_flags"),fraction)
                any_segment_flags=set()
                for f in row.get("rail_flags") or []:
                    if isinstance(f,dict):
                        any_segment_flags.update(f.get("values") or [])
                out.append({"original_independent_rail_km_point_draw_index":idx,
                    "sampled_original_GERS_segment_SHA256":hashlib.sha256(
                        str(row["id"]).encode("utf8")).hexdigest(),
                    "source_original_2026_WGS84_track_point_lon":round(lon,7),
                    "source_original_2026_WGS84_track_point_lat":round(lat,7),
                    "original_S2_globally_uniform_rail_km_fraction":round(fraction,12),
                    "source_original_segment_length_m":round(L,3),
                    "source_original_rail_class":row["class"],
                    "source_original_known_rail_flags_APPLYING_AT_SAMPLED_POINT":sorted(flagval),
                    "source_original_any_flags_SOMEWHERE_ON_SEGMENT_NOT_AT_POINT":sorted(any_segment_flags),
                    "source_overpass_or_tunnel_candidate":bool(flagval&{"is_bridge","is_tunnel","is_covered"}),
                    "source_site_geocode_is_synthetic_rail_km_point_NOT_plant_or_observer_GPS":True})
                current+=1
        acc=hi
        if len(out) and len(out)%500==0 and counts.get("last_progress_n")!=len(out):
            counts["last_progress_n"]=len(out)
            print(json.dumps({"global_original_rail_km_drawn":len(out),
                "secondpass_seconds":round(time.monotonic()-start,1)}),flush=True)
    if (len(out)!=SAMPLE_K or current!=SAMPLE_K or seen!=2054855 or
        abs(acc-distance_m)>.05):
        raise ValueError("Original global 16-source rail-km 2048 points not complete/consistent")
    out.sort(key=lambda x:x["original_independent_rail_km_point_draw_index"])
    if [r["original_independent_rail_km_point_draw_index"] for r in out]!=list(range(SAMPLE_K)):
        raise ValueError("Original 2048 independent world source rail-km random draw indices missing")
    return out

def aggregate_by_continent_using_original_NaturalEarth(points,geometry):
    from shapely.geometry import shape,Point
    from shapely.ops import unary_union
    from shapely.prepared import prep
    from collections import defaultdict
    # This is only for DESCRIPTIVE regions of already sampled world source points.
    # Do not retrospectively replace a rail-km draw based on the continent.
    j=json.loads(Path(geometry).read_text(encoding="utf8"))
    polys=defaultdict(list)
    for feature in j["features"]:
        cont=feature["properties"].get("CONTINENT")
        if cont not in ("Africa","Asia","Europe","North America","South America","Oceania"):continue
        g=shape(feature["geometry"])
        if not g.is_valid:
            from shapely.validation import make_valid
            g=make_valid(g)
        if g.geom_type=="GeometryCollection":
            from shapely.geometry import MultiPolygon,Polygon
            z=[part for part in g.geoms if part.geom_type in ("Polygon","MultiPolygon")]
            if not z:continue
            g=unary_union(z)
        if feature["properties"].get("GEOUNIT")=="Russia":
            from shapely.geometry import box
            polys["Europe"].append(g.intersection(box(-170,-90,60,90)))
            polys["Asia"].append(g.intersection(box(60,-90,180,90)).union(g.intersection(box(-180,-90,-170,90))))
        else:polys[cont].append(g)
    preps={x:prep(unary_union(vals)) for x,vals in polys.items()}
    count=Counter();outside=0
    for x in points:
        p=Point(x["source_original_2026_WGS84_track_point_lon"],
                x["source_original_2026_WGS84_track_point_lat"])
        matches=[c for c,poly in preps.items() if poly.covers(p)]
        count[matches[0] if len(matches)==1 else "ocean_or_ambiguous_original_map_units"]+=1
    return dict(count)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--rail-dir",required=True)
    p.add_argument("--original-v357-census",required=True)
    p.add_argument("--natural-earth")
    p.add_argument("--output-points",required=True)
    p.add_argument("--output-summary",required=True)
    args=p.parse_args()
    source=json.loads(Path(args.original_v357_census).read_text(encoding="utf8"))
    paths=validate_files(args.rail_dir,source)
    if len(paths)!=EXPECTED_WORLD_SOURCE_PARQUET_N:
        raise ValueError("Not original 16 complete rail Parquet shards")
    total,qa=census_global_source(paths)
    sample=select_source_global_rail_km(paths,total)
    group=aggregate_by_continent_using_original_NaturalEarth(sample,args.natural_earth) if args.natural_earth else None
    rec={"version":"v368b-original-fullworld-physical-rail-2048-WGS84-length-uniform-random-points",
         "original_source_Overture_release":SOURCE_RELEASE,
         "original_source_GERS_S3_inventory_sha256":EXPECTED_INVENTORY_SHA,
         "original_128_global_source_parquet_complete_verified":True,
         "sampling_frame_original_fullworld_candidate_normal_rail_metres":round(total,3),
         "source_rail_all_classes_incidence_census":qa,
         "n_independent_uniform_world_rail_km_draws_with_replacement":len(sample),
         "sample_seed_SHA256":hashlib.sha256(SEED.encode()).hexdigest(),
         "sampling_strategy":"Independent U(0,sum eligible WGS84 source rail segment geodesic length) cumulative-target draws WITH replacement, not uniform rail segments or selected cities",
         "source_sampled_rail_flag_counts":dict(Counter(
            flag for row in sample for flag in row["source_original_known_rail_flags_APPLYING_AT_SAMPLED_POINT"])),
         "source_sampled_original_rail_class_counts":dict(Counter(
            row["source_original_rail_class"] for row in sample)),
         "source_draws_by_original_continent_map_units":group,
         "source_sampled_bridge_tunnel_covered_candidates":sum(row["source_overpass_or_tunnel_candidate"] for row in sample),
         "source_rail_flags_geometrically_scoped_at_actual_random_source_location":True,
         "original_normal_candidate_population_excludes_ENTIRE_SEGMENT_if_any_inactive_flag_even_partially_scoped":True,
         "source_rail_sample_full_point_GPS_only_CI_artifact_NOT_original_observer_GPS":True,
         "source_none_of_sample_selected_using_plant_data":True,
         "original_rail_line_length_sampling_NOT_botanical_plot_or_named_train_route":True,
         "original_rail_sampled_GERS_line_source_ODbL":"Overture Maps 2026-11-18.0 / OpenStreetMap contributors",
         "original_plant_GPS_or_community_Jaccard_NOT_processed":True}
    path=Path(args.output_points);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps({"source_global_rail_sample_meta":rec,
         "original_public_rail_point_draws_no_botanist_coordinates":sample},ensure_ascii=False,indent=2))
    target=Path(args.output_summary);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(rec,ensure_ascii=False,indent=2))
    print(json.dumps(rec,ensure_ascii=False),flush=True)

if __name__=="__main__":main()

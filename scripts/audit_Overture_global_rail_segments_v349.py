#!/usr/bin/env python3
"""v349: Extract independent, worldwide-coverage Overture rail SEGMENTS in
the SIX v348 predeclared sampling windows from official Overture GeoParquet.

This is an original rail geometry source gate and a possible replacement for
burdening public Overpass/OSM full map APIs at continental sampling scale.
It does not download iNaturalist observations or assert plant ecology.
"""
import argparse,hashlib,importlib.metadata,json,math,sys
from collections import Counter
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from pilot_self_extracted_six_continent_rail_flora_v348 import CITIES,bounds,project

ALLOWED_RAIL_CLASS={"standard_gauge","narrow_gauge","unknown"}
EXCLUDED_FLAGS={"is_tunnel","is_bridge","is_abandoned","is_disused","is_under_construction","is_covered"}
MAX_SOURCE_SEGMENTS_PER_WINDOW=300000

def eligible_segment(row):
    if row.get("subtype")!="rail":return False,"not_rail"
    if row.get("class") not in ALLOWED_RAIL_CLASS:return False,"wrong_class"
    if row.get("geometry") is None:return False,"no_geometry"
    flags=row.get("rail_flags") or []
    for f in flags:
        if not isinstance(f,dict):continue
        if EXCLUDED_FLAGS.intersection(f.get("values") or []):
            return False,"rail_flag_excluded"
    return True,None

def line_length_m(lines,lat0,lon0):
    length=0.0
    for line in lines:
        coords=list(line.coords)
        for (x,y), (x2,y2) in zip(coords,coords[1:]):
            a=project(y,x,lat0,lon0);b=project(y2,x2,lat0,lon0)
            length+=math.hypot(a[0]-b[0],a[1]-b[1])
    return length

def inspect(city):
    from overturemaps import record_batch_reader
    from shapely import wkb
    name,continent,country,lat,lon=city
    a,b,c,d=bounds(lat,lon)
    # Overture bbox is west, south, east, north, unlike our (S,W,N,E).
    bbox=(b,a,d,c)
    out={"city":name,"continent":continent,"country":country,
         "window_halfspan_degrees":0.012,
         "source_type":"Overture global transportation subtype rail"}
    try:
        reader=record_batch_reader("segment",bbox=bbox)
        if reader is None:raise ValueError("No original Overture source reader")
        tab=reader.read_all()
        if tab.num_rows>MAX_SOURCE_SEGMENTS_PER_WINDOW:
            raise ValueError("Overture source row bound exceeded")
        rows=tab.to_pylist()
        tally=Counter()
        eligible_lines=[]
        for row in rows:
            good,why=eligible_segment(row)
            if not good:
                tally[why]+=1;continue
            geom=wkb.loads(bytes(row["geometry"]))
            if geom.geom_type=="LineString":
                eligible_lines.append(geom)
            elif geom.geom_type=="MultiLineString":
                eligible_lines.extend(geom.geoms)
            else:tally["other_geometry"]+=1
        out.update({"source_total_transportation_segments":len(rows),
                    "source_rail_features_eligible":len(eligible_lines),
                    "source_rail_line_length_km":round(line_length_m(eligible_lines,lat,lon)/1000,4),
                    "source_filter_outcomes":dict(sorted(tally.items())),
                    "source_success":True,
                    "source_original_coordinates_or_segments_not_exported":True})
    except Exception as exc:
        out.update({"source_success":False,
                    "source_error":type(exc).__name__,
                    "source_error_description":str(exc)[:220]})
    return out

def run():
    import overturemaps
    version=importlib.metadata.version("overturemaps")
    out=[inspect(c) for c in CITIES]
    return {"version":"v349-six-continent-Overture-original-rail-source",
      "source_url":"https://docs.overturemaps.org/guides/transportation/",
      "source_cloud_public":"s3://overturemaps-us-west-2/release/2026-11-18.0/theme=transportation/type=segment/",
      "source_python_client_version":version,
      "temporal_mismatch_2026_rail_with_2025_plants":True,
      "source_rail_flags_conservatively_excluded":sorted(EXCLUDED_FLAGS),
      "source_rail_class_included":sorted(ALLOWED_RAIL_CLASS),
      "predeclared_fixed_six_city_windows_not_global_probability_sample":True,
      "individual_railway_segments_not_published":True,
      "regions":out,"n_regions_source_success":sum(x.get("source_success",False) for x in out),
      "no_rail_thermal_or_biological_effect_in_this_source_qa":True}

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);a=p.parse_args()
    r=run();dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False),flush=True)

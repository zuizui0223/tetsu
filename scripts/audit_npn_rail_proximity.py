#!/usr/bin/env python3
"""Join 2024 USA-NPN monitored sites to public 2026 FRA railway GIS.

Exploratory DATA COVERAGE audit: no disturbance, seed-transport, rail wind or
phenology causal effects. Source NARN was updated 2026-10-19, AFTER NPN 2024.
NPN site coordinates do not establish a plant growing on a railway verge.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import sys
import requests
from pyproj import Transformer
from shapely.geometry import shape,Point
from shapely.ops import transform
from shapely.strtree import STRtree

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.npn_individuals import normalize_open_flowers
from audit_usa_npn_status import PILOTS,ENDPOINT,form_args,parse_stream

LAYER="https://services.arcgis.com/xOi1kZaI0eWDREZv/ArcGIS/rest/services/NTAD_North_American_Rail_Network_Lines/FeatureServer/0"
ATTRS="OBJECTID,STATEAB,PASSNGR,STRACNET,NET,TRACKS,KM,RROWNER1,YARDNAME"
UTM_EPSG={"MD":32618,"AZ":32612}
BIN_THRESHOLDS=(250.,1000.)
MAX_FEATURES_PER_STATE=15000
SESSION=requests.Session()
SESSION.headers.update({"User-Agent":"tetsu-phenology-rail-overlay/0.6 (https://github.com/zuizui0223/tetsu)"})


def arc_get(path,params):
    r=SESSION.get(path,params=params,timeout=(15,120))
    r.raise_for_status()
    data=r.json()
    if not isinstance(data,dict) or "error" in data:
        raise RuntimeError(f"FRA ArcGIS query error: {str(data)[:220]}")
    return data


def get_rail_features(state):
    if state not in UTM_EPSG: raise ValueError("state not predeclared")
    metadata=arc_get(LAYER,{"f":"json"})
    names={x.get("name") for x in metadata.get("fields",[])}
    needed=set(ATTRS.split(","))
    missing=needed-names
    if missing: raise ValueError(f"missing expected FRA fields: {sorted(missing)}")
    query=LAYER+"/query"
    where=f"STATEAB='{state}'"
    args={"where":where,"f":"json","returnCountOnly":"true"}
    total=int(arc_get(query,args).get("count",-1))
    if not 0<total<=MAX_FEATURES_PER_STATE:
        raise ValueError(f"unexpected FRA {state} feature count {total}")
    lines=[]
    flags=Counter()
    ids=set()
    offsets=range(0,total,1000)
    for offset in offsets:
        d=arc_get(query,{
            "f":"geojson","where":where,
            "outFields":ATTRS,"returnGeometry":"true","outSR":4326,
            "orderByFields":"OBJECTID ASC","resultOffset":offset,"resultRecordCount":1000
        })
        features=d.get("features",[])
        if not features:
            raise RuntimeError(f"no page {offset} for {state}; truncate forbidden")
        for feat in features:
            properties=feat.get("properties") or {}
            gid=properties.get("OBJECTID")
            if gid is None or gid in ids:
                raise ValueError(f"Duplicate/missing FRA ID at offset {offset}")
            ids.add(gid)
            geom=feat.get("geometry")
            if not geom:
                continue
            g=shape(geom)
            if g.is_empty or g.geom_type not in ("LineString","MultiLineString"):
                continue
            lines.append(g)
            passenger=str(properties.get("PASSNGR") or "missing")
            flags[f"passenger_{passenger}"]+=1
            flags[f"network_{str(properties.get('NET') or 'missing')}"]+=1
    if len(ids)!=total or not lines:
        raise ValueError(f"FRA pages incomplete for {state}: unique={len(ids)} expected={total}")
    transform_xy=Transformer.from_crs(4326,UTM_EPSG[state],always_xy=True)
    metric_lines=[transform(transform_xy.transform,g) for g in lines]
    tree=STRtree(metric_lines)
    def proximity(lon,lat):
        pt=Point(*transform_xy.transform(lon,lat))
        nearest=int(tree.nearest(pt))
        return float(pt.distance(metric_lines[nearest]))
    return proximity,{
        "source":LAYER,"n_fra_features":total,
        "n_nonempty_rail_geometries":len(lines),
        "passenger_and_network_field_counts":dict(flags),
        "utm_epsg":UTM_EPSG[state],
        "source_version_note":"FRA NARN 2026-10-19; geography may differ from NPN 2024 exposure"
    }


def bin_for_distance(distance):
    if distance<=BIN_THRESHOLDS[0]: return "0-250m"
    if distance<BIN_THRESHOLDS[1]: return "250-1000m"
    return "1000m_plus"


def run():
    report={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "description":"2024 monitored NPN phenology site proximity to 2026 FRA railway geometry; no rail effect",
        "source_npn":ENDPOINT,"source_fra":LAYER,
        "cases":[],"failures":{}
    }
    for case in PILOTS:
        state=case["state"];key=f"{state}_{case['year']}_spp{case['species_id']}"
        try:
            with requests.post(ENDPOINT,data=form_args(case),stream=True,
                               timeout=(20,140),headers={"User-Agent":"tetsu-phenology-rail-overlay/0.6"}) as response:
                response.raise_for_status()
                raw,npn_sha,npn_bytes=parse_stream(response)
            observed,excluded=normalize_open_flowers(raw)
            sites=defaultdict(list)
            all_sites={}
            for r in raw:
                site=r.get("site_id")
                try:
                    lat=float(r.get("latitude"));lon=float(r.get("longitude"))
                    if not(math.isfinite(lat) and math.isfinite(lon)):
                        continue
                    if not (-125<=lon<=-65 and 24<=lat<=50):
                        continue
                    if site is None: continue
                except (TypeError,ValueError):
                    continue
                sites[str(site)].append((lon,lat))
            # A site's coords must agree across its repeated visits or be
            # rejected, rather than quietly assigning them to a rail band.
            for site,points in sites.items():
                if (max(p[0] for p in points)-min(p[0] for p in points)>0.001 or
                    max(p[1] for p in points)-min(p[1] for p in points)>0.001):
                    continue
                all_sites[site]=points[0]
            nearest,rail_source=get_rail_features(state)
            bins={}
            for site,(lon,lat) in all_sites.items():
                bins[site]=bin_for_distance(nearest(lon,lat))
            site_stats=defaultdict(Counter)
            for row in observed:
                b=bins.get(row["site_id"])
                if b is None:
                    continue
                site_stats[b][row["status"]]+=1
            counts=Counter(bins.values())
            result={
                "state":state,"species_id":case["species_id"],
                "expected_species":case["expected_species"],"year":case["year"],
                "npn_source_sha256":npn_sha,
                "npn_observations":len(raw),
                "npn_openflower_status_records":len(observed),
                "n_npnp_sites_with_coherent_usable_coordinates":len(all_sites),
                "n_npnp_sites_excluded_or_missing_coordinates":len(sites)-len(all_sites),
                "rail_site_proximity_bins":dict(counts),
                "openflower_yes_no_by_proximity":{
                    b:{"yes":c[1],"no":c[0]} for b,c in site_stats.items()
                },
                "rail_source":rail_source,
                "interpretation":"Source coverage only; site coordinate is NOT proof of railway habitat, train wind/cutting, or plant rail-mediated flowering adaptation."
            }
            report["cases"].append(result)
            print(json.dumps({"case":key,"rail_sites":result["rail_site_proximity_bins"],
                              "openflower_visits":result["openflower_yes_no_by_proximity"]},
                             ensure_ascii=False),flush=True)
        except Exception as exc:
            report["failures"][key]=f"{type(exc).__name__}: {str(exc)[:280]}"
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/usa_npn_rail_proximity.json")
    args=p.parse_args()
    result=run()
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"n_source_successes":len(result["cases"]),"errors":result["failures"]},ensure_ascii=False))
    if not result["cases"]: raise SystemExit("No NPN-to-FRA site overlay was completed")


if __name__=="__main__":
    main()

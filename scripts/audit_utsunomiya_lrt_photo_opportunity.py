#!/usr/bin/env python3
"""2023 Utsunomiya LRT as a verified new surface-rail operation event pilot.

Official operator: began regular passenger service August 26, 2023.
Rail construction began in 2018 and trials in 2022, SO 'pre-operation'
is NOT an untreated no-disturbance baseline.

Source iNaturalist 2018–2025 photos; compare all plant PHOTO opportunity
near 2023 N05 surface LRT vs 1–3 km away in a geometry-defined corridor.
No road, local land-use, plant microhabitat or observer effort causal
adjustment: not a rail dispersal, phenology, mowing or evolution effect.
No raw coordinates, IDs, route geometry or photo files are published.
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
import zipfile

from pyproj import Transformer
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree

from audit_n05_historical_railway import (
    SOURCE, known_year, is_current_end, shape_reader,bounded_zip_download
)
from audit_inat_rail_phenophases import (
    RateLimitedClient, URL, acceptable_observation
)

ROUTE="宇都宮芳賀ライトレール線"
COMPANY="宇都宮ライトレール"
OPENING_DATE="2023-08-26"
YEARS=tuple(range(2018,2026))
UTM_EPSG=32654
FOCUS={
    "Equisetum arvense":"horsetail_spore_plant",
    "Solidago altissima":"goldenrod",
    "Setaria":"annual_foxtail_genus",
    "Imperata cylindrica":"cogon_grass",
    "Miscanthus sinensis":"silvergrass",
    "Pueraria montana lobata":"kudzu_variety",
}
QUERY_PADDING_DEGREES=0.055


def rail_distance_function(source_zip):
    t=Transformer.from_crs(4326,UTM_EPSG,always_xy=True)
    with zipfile.ZipFile(io.BytesIO(source_zip)) as archive:
        reader,shp=shape_reader(archive,"RailroadSection2")
        fields=[f[0] for f in reader.fields[1:]]
        lines=[]
        bounds=[]
        source_n=0
        for feature in reader.iterShapeRecords():
            record=dict(zip(fields,feature.record))
            if (str(record.get("N05_002") or "").strip()!=ROUTE or
                str(record.get("N05_003") or "").strip()!=COMPANY):
                continue
            if known_year(record.get("N05_004"))!=2023:
                continue
            if not is_current_end(record.get("N05_005e")):
                continue
            shape=feature.shape
            if shape.shapeType not in (3,13,23):
                raise ValueError("Verified N05 LRT geometry is not a polyline")
            source_n+=1
            bounds.append(shape.bbox)
            breaks=list(shape.parts)+[len(shape.points)]
            for a,b in zip(breaks,breaks[1:]):
                if b-a<2:continue
                vertices=[t.transform(x,y) for x,y in shape.points[a:b]]
                line=LineString(vertices)
                if line.length>1:lines.append(line)
        reader.close()
    if not lines or source_n<1:
        raise ValueError("No independent N05 2023 Utsunomiya LRT records; do NOT create railway exposure proxy")
    xmin=min(x[0] for x in bounds)
    ymin=min(x[1] for x in bounds)
    xmax=max(x[2] for x in bounds)
    ymax=max(x[3] for x in bounds)
    envelope=(ymin-QUERY_PADDING_DEGREES,xmin-QUERY_PADDING_DEGREES,
              ymax+QUERY_PADDING_DEGREES,xmax+QUERY_PADDING_DEGREES)
    tree=STRtree(lines)
    def distance(x,y):
        pt=Point(*t.transform(x,y))
        near_idx=int(tree.nearest(pt))
        d=float(pt.distance(lines[near_idx]))
        grid=(math.floor(pt.x/1000), math.floor(pt.y/1000))
        return d,grid
    return distance,envelope,{
        "source_n05_2023_qualified_line_records":source_n,
        "indexed_line_parts":len(lines),
        "source_shapefile_basename":Path(shp).name,
        "official_opening_date":OPENING_DATE,
        "real_rail_construction_and_test_before_opening":True,
        "source_geom_cannot_be_used_as_causal_dose":True,
        "approx_line_length_km":round(sum(g.length for g in lines)/1000,2),
    }


def query_year(client,year,envelope):
    lat0,lon0,lat1,lon1=envelope
    params={
        "taxon_id":47126, # iNaturalist Plantae
        "quality_grade":"research","photos":"true",
        "d1":f"{year}-01-01","d2":f"{year}-12-31",
        "swlat":round(lat0,5),"swlng":round(lon0,5),
        "nelat":round(lat1,5),"nelng":round(lon1,5),
        "per_page":200,"order_by":"id","order":"asc"
    }
    first=client.get(URL,params={**params,"page":1}).json()
    total=int(first.get("total_results") or 0)
    if total>10000:
        raise RuntimeError(f"Year {year}: {total} exceeds iNaturalist 10k API search window; date-split is required")
    rows=list(first.get("results") or [])
    for page in range(2,math.ceil(total/200)+1):
        rows.extend(client.get(URL,params={**params,"page":page}).json().get("results") or [])
    if len(rows)!=total or len({x["id"] for x in rows})!=total:
        raise RuntimeError(f"Paging inconsistency year={year}")
    return rows,total


def exposure_bin(m):
    if not math.isfinite(m) or m<0:raise ValueError("invalid distance")
    if m<=250:return "near_0_250m"
    if 1000<=m<=3000:return "control_1_3km"
    return "excluded_mid_or_outer"


def organism_group(ob):
    taxon=(ob.get("taxon") or {}).get("name") or ""
    if taxon.startswith("Setaria "):
        return "Setaria"
    for name in FOCUS:
        if taxon==name or taxon.startswith(name+" "):
            return name
    return None


def summarize_photo_year(obs,dist_func):
    counts=Counter()
    observers=defaultdict(set)
    grids=defaultdict(set)
    species=defaultdict(lambda:defaultdict(set))
    focus=defaultdict(Counter)
    grid_phase=defaultdict(set)
    rejected=Counter()
    for row in obs:
        okay,reason=acceptable_observation(row)
        rejected[reason]+=1
        if not okay:
            continue
        lon,lat=row["geojson"]["coordinates"]
        d,grid=dist_func(lon,lat)
        group=exposure_bin(d)
        counts[group]+=1
        user=(row.get("user") or {}).get("id")
        if user is not None:observers[group].add(str(user))
        grids[group].add(grid)
        taxon=(row.get("taxon") or {}).get("name") or ""
        if taxon:
            species[taxon][group].add(grid)
        selected=organism_group(row)
        if selected:
            focus[selected][group]+=1
            grid_phase[(selected,group)].add(grid)
    near="near_0_250m"
    far="control_1_3km"
    common_species=sorted([
        name for name,g in species.items()
        if g.get(near) and g.get(far)
    ])
    return {
        "api_geo_quality":dict(rejected),
        "photo_counts_by_rail_proximity":dict(counts),
        "distinct_photographers_by_bin":{
            k:len(v) for k,v in observers.items()
        },
        "distinct_1km_cells_by_bin":{
            k:len(v) for k,v in grids.items()
        },
        "distinct_taxa_with_photo_near":sum(bool(x.get(near)) for x in species.values()),
        "distinct_taxa_with_photo_control":sum(bool(x.get(far)) for x in species.values()),
        "n_taxa_with_both_exposure_bins":len(common_species),
        "focal_source_taxa_annual_photo_counts":{
            name:dict(entries) for name,entries in focus.items()
        },
        "focal_taxa_spatial_cell_counts":{
            name:{
                group:len(grid_phase[(name,group)]) for group in (near,far)
            }
            for name in FOCUS
        },
        "quality_interpretation":"No photos were independently verified as trackside plants, nor is a photo non-detection an absence."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/utsunomiya_lrt_railway_opening_photo_opportunity.json")
    args=p.parse_args()
    source=bounded_zip_download()
    dist,envelope,route=rail_distance_function(source)
    client=RateLimitedClient(delay_seconds=1.1)
    reports={}
    for year in YEARS:
        rows,n=query_year(client,year,envelope)
        m=summarize_photo_year(rows,dist)
        m["api_total_query"] = n
        m["period"]="pre_regular_operation_INCLUDES_CONSTRUCTION" if year<2023 else (
            "opening_transition_EXCLUDED_from_before_after" if year==2023
            else "post_regular_operation"
        )
        reports[str(year)]=m
        print(json.dumps({
            "year":year,"all_photo_plant_records":n,
            "quality_near":m["photo_counts_by_rail_proximity"].get("near_0_250m",0),
            "quality_control":m["photo_counts_by_rail_proximity"].get("control_1_3km",0),
            "near_1km_cells":m["distinct_1km_cells_by_bin"].get("near_0_250m",0),
            "control_1km_cells":m["distinct_1km_cells_by_bin"].get("control_1_3km",0),
            "same_species_both_bins":m["n_taxa_with_both_exposure_bins"]
        },ensure_ascii=False),flush=True)
    result={
        "computed_utc":datetime.now(timezone.utc).isoformat(),
        "source":URL,
        "mlit_N05_source":SOURCE,
        "mlit_N05_sha256":hashlib.sha256(source).hexdigest(),
        "mlit_N05_license":"noncommercial source; only derived aggregated metadata, not actual source geometry",
        "line":ROUTE,
        "operator":COMPANY,
        "rail_geometry":route,
        "source_box_round_2decimal_degrees":[round(v,2) for v in envelope],
        "calendar":reports,
        "key_qualification":(
            "This is a pre/post OPERATION RECORDING-OPPORTUNITY audit only; "
            "2018-22 construction already changed habitat, 2023 transition "
            "excluded, photograph effort not randomized."
        ),
        "limits":[
            "A new tram line overlaps roads and formerly built land; proximity to route is not independently verified railway verge vegetation or train-carried propagule.",
            "City-specific iNaturalist observer effort, photo licensing, species ascertainment and urban redevelopment may differ across pre/post years.",
            "LRT service commenced 2023-08-26 but construction started 2018; even a change near the line cannot identify the effect of train operation separate from construction.",
            "No railway grass-cutting regimen, track-side biomass/flower/spore ROI or train aerodynamic gust is measured.",
            "Taxon-specific source opportunities are NOT local occupancy, establishment, seed dispersal or plant reproduction effect estimates.",
            "Each source observation is one public photo, not a complete standardized vegetation census; repeated photographers and grids are correlated.",
            "No inferential p values, first appearance dates as species colonization, or heritable adaptation conclusions are justified."
        ],
        "operator_event_source":"https://www.miyarail.co.jp/info/1143",
        "mlit_event_source":"https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N05-2025.html"
    }
    output=Path(args.output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    n_near=sum(r["photo_counts_by_rail_proximity"].get("near_0_250m",0)
               for r in reports.values())
    n_control=sum(r["photo_counts_by_rail_proximity"].get("control_1_3km",0)
               for r in reports.values())
    print(json.dumps({
        "photo_source_count_near_8years":n_near,
        "photo_source_count_control_8years":n_control,
        "pre_near":sum(reports[str(y)]["photo_counts_by_rail_proximity"].get("near_0_250m",0)
                       for y in range(2018,2023)),
        "post_near":sum(reports[str(y)]["photo_counts_by_rail_proximity"].get("near_0_250m",0)
                        for y in (2024,2025)),
        "report":str(output)
    },ensure_ascii=False))
    if n_near==0 or n_control==0:
        raise SystemExit("N05 LRT corridor lacks within-region photograph opportunity on both near and controls")


if __name__=="__main__":
    main()

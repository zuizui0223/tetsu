#!/usr/bin/env python3
"""v1.7: geography-frozen, *outcome-unseen* Japan replication of the v1.6
rail/road land-cover selection bias. This is a source-selection audit,
NOT a railway botanical, mowing, wind or propagule-transport effect.

Discovery: Sendai/Niigata/Nagano 370 transects (v1.6).
Replication: eight spatially separate first-level N13 map meshes selected
on geography/rail viability before their 2021 ESA classes are inspected.
Published summaries never contain source geometries or exact locations.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
from datetime import datetime, timezone
import io
import json
import math
from pathlib import Path
import sys
import zipfile

import rasterio
from rasterio.windows import Window
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.worldcover_landcover import (
    CLASS_LABELS, DOI, POTENTIAL_GROUND_HERBS,
    SOURCE_BASE, SOURCE_ALT, normalize_class, summarize_neighborhood,
    tile_name, tile_url,
)
from audit_japan_ordinary_rail_road_separation import (
    REGIONS, N02, MAX_ZIP_RAIL, analyze_region, download_source,
)
from audit_japan_ordinary_rail_worldcover_expanded import (
    ROAD_NEAR_THRESHOLDS, source_comparator_eligibility,
)

# All 8 new geographic windows are separate from the 3 v1.6 discovery
# regions. The N13 first-level mesh is derived from lat and lon, not
# selected from the WorldCover vegetation outcomes.
HOLDOUT = {
    "sapporo_hokkaido": {"mesh": "6441", "utm": 32654, "bbox": (141.12, 42.92, 141.71, 43.28)},
    "akita": {"mesh": "5940", "utm": 32654, "bbox": (140.05, 39.55, 140.72, 39.98)},
    "tokyo_west": {"mesh": "5339", "utm": 32654, "bbox": (139.17, 35.39, 139.78, 35.88)},
    "nagoya_aichi": {"mesh": "5236", "utm": 32653, "bbox": (136.58, 34.92, 136.97, 35.28)},
    "osaka": {"mesh": "5235", "utm": 32653, "bbox": (135.10, 34.69, 135.75, 35.04)},
    "hiroshima": {"mesh": "5132", "utm": 32653, "bbox": (132.19, 34.14, 132.86, 34.57)},
    "fukuoka": {"mesh": "5030", "utm": 32652, "bbox": (130.17, 33.39, 130.83, 33.91)},
    "kumamoto": {"mesh": "4930", "utm": 32652, "bbox": (130.34, 32.75, 130.93, 33.19)},
}


def validate_holdout_frame(regions=HOLDOUT):
    """Fail closed on N13 mesh crossing, bad projection or ESA tile split."""
    if len(regions) != len(set(regions)) or not regions:
        raise ValueError("empty/duplicate region frame")
    for name, cfg in regions.items():
        lon0, lat0, lon1, lat1 = cfg["bbox"]
        if not (lon0 < lon1 and lat0 < lat1):
            raise ValueError(f"bad bounds {name}")
        for lon, lat in ((lon0,lat0),(lon0,lat1),(lon1,lat0),(lon1,lat1)):
            mesh = f"{math.floor(lat * 1.5):02d}{math.floor(lon - 100):02d}"
            if mesh != cfg["mesh"]:
                raise ValueError(f"first-level N13 mesh crossing {name}: {mesh}")
            if tile_name(lon, lat) != tile_name(lon0, lat0):
                raise ValueError(f"WorldCover COG tile crossing {name}")
            zone = math.floor((lon + 180) / 6) + 1
            if cfg["utm"] != 32600 + zone:
                raise ValueError(f"bad projected UTM zone {name}")
        if cfg["mesh"] in {x["mesh"] for x in REGIONS.values()}:
            raise ValueError(f"holdout reused a discovery N13 mesh: {name}")


def read_esa_point(dataset, x, y, epsg, expected_tile):
    """Observed 2021 class + 3x3 map, with region-specific UTM conversion."""
    lon, lat = Transformer.from_crs(epsg, 4326, always_xy=True).transform(x,y)
    if tile_name(lon, lat) != expected_tile:
        raise ValueError("outside prespecified WorldCover tile")
    if dataset.count != 1 or dataset.crs.to_epsg() != 4326:
        raise ValueError("invalid ESA COG band / CRS")
    if not (dataset.bounds.left <= lon < dataset.bounds.right
            and dataset.bounds.bottom < lat <= dataset.bounds.top):
        raise ValueError("source point outside ESA COG extent")
    row, col = dataset.index(lon,lat)
    if not (1 <= row < dataset.height - 1 and 1 <= col < dataset.width - 1):
        raise ValueError("3x3 window crosses ESA COG edge")
    arr = dataset.read(1, window=Window(col-1,row-1,3,3))
    if arr.shape != (3,3):
        raise ValueError("incomplete source pixel window")
    return {
        "class": normalize_class(int(arr[1,1])),
        "context": summarize_neighborhood(arr.flatten().tolist()),
    }


def open_esa_tile(tile, stack):
    errors=[]
    for source in (SOURCE_BASE,SOURCE_ALT):
        url=tile_url(tile,source)
        try:
            dataset=rasterio.open(url)
            if (dataset.count!=1 or dataset.crs.to_epsg()!=4326
                    or dataset.width<10000 or dataset.height<10000
                    or dataset.nodata not in (0,None)):
                dataset.close()
                raise ValueError("bad WorldCover 2021 source raster schema")
            stack.callback(dataset.close)
            return dataset, url
        except Exception as exc:
            errors.append(f"{type(exc).__name__}: {str(exc)[:100]}")
    raise RuntimeError(f"no verified ESA WorldCover COG for {tile}: {errors}")


def summarize_region(site_rows):
    """Aggregate real fixed transects without turning them into independent plots."""
    classes=Counter()
    route_counts=Counter()
    thresholds={
        str(thr): {
            "n":0,"near_tree_n":0,"near_grass_or_wetland_n":0,
            "both_near_and_road_grass_n":0,"same_near_road_class_n":0,
            "route_names":set(),
        } for thr in ROAD_NEAR_THRESHOLDS
    }
    for site in site_rows:
        nc=site["near_class"]
        rc=site["road_class"]
        classes[CLASS_LABELS.get(nc,"unknown")]+=1
        route=(site["route_name"],site["operator_name"])
        route_counts[route]+=1
        for thr in ROAD_NEAR_THRESHOLDS:
            if not source_comparator_eligibility(site,int(thr)):
                continue
            item=thresholds[str(thr)]
            item["n"]+=1
            item["near_tree_n"]+=int(nc==10)
            item["near_grass_or_wetland_n"]+=int(nc in POTENTIAL_GROUND_HERBS)
            item["both_near_and_road_grass_n"]+=int(
                nc in POTENTIAL_GROUND_HERBS and rc in POTENTIAL_GROUND_HERBS
            )
            item["same_near_road_class_n"]+=int(
                nc is not None and rc is not None and nc==rc
            )
            item["route_names"].add(route)
    base_n=len(site_rows)
    base_tree_n=classes["tree_cover"]
    for value in thresholds.values():
        value["n_distinct_route_labels"]=len(value.pop("route_names"))
        value["near_tree_fraction"]=(round(value["near_tree_n"]/value["n"],6)
                                     if value["n"] else None)
        value["near_tree_minus_all_fixed_fraction"]=(
            round(value["near_tree_n"]/value["n"]-base_tree_n/base_n,6)
            if value["n"] and base_n else None
        )
    return {
        "n_geo_selected_rail_normals":base_n,
        "n_route_labels_not_independent_treatments":len(route_counts),
        "near_2021_WorldCover_class_counts":dict(sorted(classes.items())),
        "all_fixed_tree_fraction":round(base_tree_n/base_n,6) if base_n else None,
        "n_classified":sum(v for k,v in classes.items() if k!="unknown"),
        "road_control_eligible_by_near_road_min_m":thresholds,
        "no_pseudoreplicated_p_value":True,
        "not_a_vegetation_or_adaptation_effect":True,
    }


def holdout_summary(regions):
    """Region is the descriptive replication unit; missing eligibility is NA."""
    keys=("25","50","75","100")
    stats={}
    for key in keys:
        paired=[]
        for name,row in regions.items():
            a=row["road_control_eligible_by_near_road_min_m"][key]
            if row["n_geo_selected_rail_normals"]>=20 and a["n"]>=3:
                paired.append((name,a["near_tree_minus_all_fixed_fraction"]))
        stats[key]={
            "minimum_fixed_normals_per_region":20,
            "minimum_road_control_eligible_normals_per_region":3,
            "n_comparable_geographically_separate_regions":len(paired),
            "n_regions_with_higher_tree_fraction":sum(d>0 for _,d in paired),
            "n_regions_with_lower_tree_fraction":sum(d<0 for _,d in paired),
            "n_regions_equal":sum(d==0 for _,d in paired),
            "median_delta_tree_fraction":(
                round(sorted(d for _,d in paired)[len(paired)//2],6)
                if len(paired)%2 else (
                    round(sum(sorted(d for _,d in paired)[len(paired)//2-1:len(paired)//2+1])/2,6)
                    if paired else None
                )
            ),
            "no_independent_plant_effect_test":True,
        }
    return stats


def run():
    validate_holdout_frame()
    if any(name in REGIONS for name in HOLDOUT):
        raise ValueError("discovery/holdout region names overlap")
    REGIONS.update(HOLDOUT)  # in-memory only: source modules are not edited
    raw,metadata=download_source(N02,MAX_ZIP_RAIL)
    if not raw or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("invalid official N02 rail source")
    result={
        "run_utc":datetime.now(timezone.utc).isoformat(),
        "stage":"v1.7 unseen-area source-class replication of v1.6 GIS sampling confound",
        "discovery_regions_not_reused":["sendai_miyagi","niigata","nagano"],
        "holdout_region_names_frozen_before_ESA_source_reads":list(HOLDOUT),
        "2025_N02_source_sha256":metadata["source_sha256"],
        "2024_N13_source_not_historical":True,
        "ESA_WorldCover_2021_DOI":DOI,
        "strict_threshold_75m_fixed_before_holdout_outcomes":True,
        "posthoc_25_50_100_thresholds_are_sensitivities":True,
        "regions":{},"source_errors":{},
        "not_a_global_biological_effect":True,
    }
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="3",
        GDAL_HTTP_RETRY_DELAY="2",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="40000000",
        GDAL_CACHEMAX=256,
    ), ExitStack() as stack:
        rasters={}
        for name,cfg in HOLDOUT.items():
            try:
                summary,locations=analyze_region(
                    name,raw,private_site_positions=True,
                    private_candidate_scope="all_fixed_geo"
                )
                if not locations:
                    raise ValueError("no source GIS transects in heldout region")
                tile=tile_name(cfg["bbox"][0],cfg["bbox"][1])
                if tile not in rasters:
                    rasters[tile]=open_esa_tile(tile,stack)
                raster,url=rasters[tile]
                rows=[]
                for loc in locations:
                    near=read_esa_point(raster,*loc["near_xy"],cfg["utm"],tile)
                    far=read_esa_point(raster,*loc["road_only_xy"],cfg["utm"],tile)
                    rows.append({
                        "route_name":loc["route_name"],
                        "operator_name":loc["operator_name"],
                        "near_class":near["class"],
                        "road_class":far["class"],
                        "near_road_m":loc["near_road_m"],
                        "far_road_m":loc["far_road_m"],
                        "far_rail_m":loc["far_rail_m"],
                        "near_station_m":loc["near_station_m"],
                    })
                if len(rows)!=summary["total_rail_side_source_opportunities"]:
                    raise ValueError("source GIS side-count mismatch")
                outcome=summarize_region(rows)
                outcome["worldcover_2021_tile"]=tile
                outcome["worldcover_source_url"]=url
                result["regions"][name]=outcome
                print(json.dumps({
                    "region":name, "n":outcome["n_geo_selected_rail_normals"],
                    "classes":outcome["near_2021_WorldCover_class_counts"],
                    "road75":outcome["road_control_eligible_by_near_road_min_m"]["75"],
                },ensure_ascii=False),flush=True)
            except Exception as exc:
                result["source_errors"][name]=f"{type(exc).__name__}: {str(exc)[:250]}"
                print(json.dumps({"source_error":name,"type":type(exc).__name__}),flush=True)
    result["heldout_city_level_selection_bias_descriptors"]=holdout_summary(result["regions"])
    result["n_complete_heldout_regions"]=len(result["regions"])
    result["source_complete"]=(not result["source_errors"] and len(result["regions"])==len(HOLDOUT))
    result["interpretation_firewall"]=[
        "8 source windows are geographically preselected after discovery but before inspecting holdout 2021 landcover.",
        "A forest-enrichment replication is a selection/confounding diagnostic, NOT evidence that trains cause vegetation selection.",
        "35m GIS point is not verified at-grade railway property or an individual grass patch.",
        "N02 passenger rail cannot identify active freight rail or actual train propagule deposition.",
        "WorldCover 2021 and roads 2024 are different vintage; no mowing/wind/flowers/plant species measured.",
        "Overlapping transects on one route are not independent replicates; no individual-transect p-values.",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="results/audit/japan_heldout_8region_worldcover_2021_v17.json")
    args=p.parse_args()
    result=run()
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":result["source_complete"],
        "regions":result["n_complete_heldout_regions"],
        "source_errors":result["source_errors"],
        "road75":result["heldout_city_level_selection_bias_descriptors"]["75"]
    },ensure_ascii=False),flush=True)
    if not result["source_complete"]:
        raise SystemExit("Incomplete 8-region holdout sources; no replication claim")


if __name__=="__main__":
    main()

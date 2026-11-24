#!/usr/bin/env python3
"""v361: Build physical source railway connector components from ALL 16 source-
complete Overture rail-only shards AFTER globally authenticated 128-Parquet
census. Never equate physical components with named operated line routes.

Two graphs: all rail classes (including metro/tram), and normal active-
CANDIDATE classes (standard/narrow/unknown, excluding disused, abandoned,
construction but keeping bridges/tunnels). Exact GERS connector UUIDs are
the only topology keys; spatial intersections alone do not connect.
Outputs only aggregate component ranks and source hashes, no source IDs,
full network geometry, train OD or botanical similarity.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED_INVENTORY_SHA="6237da737bb783f23f352638cd5aef287bbb3d7aa8f0fe0bed8c1cc4558f6b34"
SOURCE_RELEASE="2026-11-18.0"
SHARDS=16
NORMAL_CLASSES=frozenset(("standard_gauge","narrow_gauge","unknown"))
EXCLUDED_INACTIVE=frozenset(("is_abandoned","is_disused","is_under_construction"))
PANELS=("all_original_rail_subtypes","conventional_active_candidate_source_rails")
MAX_CONNECTORS=30_000_000

class Disjoint:
    """Union-find with exact source GERS connector string identity."""
    def __init__(self):
        self.parent={}
        self.sz={}
        self.components=0
    def add(self,v):
        if v not in self.parent:
            if len(self.parent)>=MAX_CONNECTORS:
                raise MemoryError("Too many real source connector IDs for safe source graph")
            self.parent[v]=v
            self.sz[v]=1
            self.components+=1
    def find(self,v):
        parent=self.parent
        if v not in parent:
            raise KeyError("Original source physical connector graph node missing")
        while parent[v]!=v:
            parent[v]=parent[parent[v]]
            v=parent[v]
        return v
    def union(self,a,b):
        self.add(a);self.add(b)
        x,y=self.find(a),self.find(b)
        if x==y:return False
        if self.sz[x]<self.sz[y] or (self.sz[x]==self.sz[y] and x>y):
            x,y=y,x
        self.parent[y]=x
        self.sz[x]+=self.sz[y]
        del self.sz[y]
        self.components-=1
        return True

def original_connectors(row):
    pairs=[]
    for entry in row.get("connectors") or ():
        if not isinstance(entry,dict):continue
        connector=str(entry.get("connector_id") or "")
        try:at=float(entry.get("at"))
        except (ValueError,TypeError):continue
        if not connector or not math.isfinite(at) or not 0<=at<=1:continue
        pairs.append((at,connector))
    # Duplicate exact source connector location pairs do not create fake nodes.
    out=sorted(set(pairs))
    return [x[1] for x in out]

def eligible_panels(row):
    if row.get("subtype")!="rail":raise ValueError("Original source v357 rail-only Parquet contains nonrail data")
    flags=set()
    for f in row.get("rail_flags") or ():
        if isinstance(f,dict):flags.update(f.get("values") or ())
    normal=(row.get("class") in NORMAL_CLASSES
            and not (flags & EXCLUDED_INACTIVE))
    return {"all_original_rail_subtypes":True,
        "conventional_active_candidate_source_rails":normal}

def read_rows(paths,columns):
    import pyarrow.parquet as pq
    for path in paths:
        pf=pq.ParquetFile(path)
        if not set(columns).issubset(pf.schema_arrow.names):
            raise ValueError("Original derived source railway Parquet lacks requested graph fields")
        for batch in pf.iter_batches(batch_size=30000,columns=list(columns),use_threads=True):
            for r in batch.to_pylist():
                yield r

def validate_files(parent,census):
    if (census.get("version")!="v357-original-entire-world-Overture-2026-all-rail-scan-complete"
            or census.get("source_release")!=SOURCE_RELEASE
            or census.get("original_S3_source_Parquet_inventory_sha256")!=EXPECTED_INVENTORY_SHA
            or census.get("source_all_segment_source_objects_complete")!=128
            or len(census.get("source_shards",[]))!=SHARDS):
        raise ValueError("Not a source-complete original 128-file Overture global railway census")
    paths=[]
    for i in range(SHARDS):
        shard=census["source_shards"][i]
        if shard["shard_index"]!=i or shard["original_transport_source_files_n"]!=8:
            raise ValueError("Original source 128 Parquet files mispartitioned/incomplete")
        name=shard["source_rail_only_Parquet_artifact_filename"]
        expected=shard["source_rail_Parquet_sha256"]
        if not name or not expected:
            raise ValueError("Original full railway segment source missing derived Parquet")
        p=Path(parent,name)
        if not p.is_file() or p.stat().st_size!=shard["source_rail_Parquet_bytes"]:
            raise ValueError("Original full-source rail-only parquet file missing/size mismatch: "+name)
        h=hashlib.sha256()
        with p.open("rb") as f:
            for block in iter(lambda:f.read(2**20),b""):h.update(block)
        if h.hexdigest()!=expected:
            raise ValueError("Original rail-only parquet sha256 source artifact mismatch: "+name)
        paths.append(p)
    return paths

def topology(paths):
    d={k:Disjoint() for k in PANELS}
    source=Counter()
    seen_ids=set()
    for row in read_rows(paths,("id","subtype","class","connectors","rail_flags")):
        gid=row.get("id")
        if not isinstance(gid,str) or not gid or gid in seen_ids:
            raise ValueError("Original global GERS source segment id absent/duplicated")
        seen_ids.add(gid)
        source["source_all_rail_segments"]+=1
        connectors=original_connectors(row)
        source["source_segment_without_connectors"]+=not bool(connectors)
        source["source_segment_one_valid_connector"]+=len(connectors)==1
        for name,eligible in eligible_panels(row).items():
            if not eligible:continue
            source["source_segment_in_"+name]+=1
            if not connectors:
                source["source_unplaceable_segment_"+name]+=1
                continue
            uf=d[name]
            for node in connectors:uf.add(node)
            for a,b in zip(connectors,connectors[1:]):
                if a!=b:uf.union(a,b)
    return d,dict(source)

def aggregate_components(paths,census,d,firstpass):
    from shapely import wkb
    from pyproj import Geod
    geod=Geod(ellps="WGS84")
    comps={name:defaultdict(lambda:[0,0.0,0]) for name in PANELS}
    counts=Counter()
    for row in read_rows(paths,("id","subtype","class","connectors","rail_flags","geometry","names")):
        connectors=original_connectors(row)
        countable=eligible_panels(row)
        for name,eligible in countable.items():
            if not eligible:continue
            counts["n_source_"+name]+=1
            if not connectors:continue
            root=d[name].find(connectors[0])
            component=comps[name][root]
            component[0]+=1
            # Count all physical centerline length at current Overture 2026 source,
            # NOT active operating rail line length or track occupancy.
            geom=wkb.loads(bytes(row["geometry"]))
            if geom.geom_type!="LineString" or geom.is_empty:
                raise ValueError("Unexpected physical original source railway geometry")
            km=abs(geod.geometry_length(geom))/1000.
            if not math.isfinite(km):
                raise ValueError("Original source railway geodesic length not finite")
            component[1]+=km
            component[2]+=bool((row.get("names") or {}).get("primary"))
    out={}
    all_original_km=census["source_real_geodesic_all_rail_km"]
    for name in PANELS:
        uf=d[name];arr=list(comps[name].values())
        arr.sort(key=lambda x:(-x[0],-x[1]))
        n_segments=firstpass["source_segment_in_"+name]
        n_placed=sum(a[0] for a in arr)
        if n_placed+firstpass.get("source_unplaceable_segment_"+name,0)!=n_segments:
            raise AssertionError("Cannot account for all genuine source rail segments")
        total_length=sum(c[1] for c in arr)
        top=arr[:10]
        out[name]={"source_original_segments":n_segments,
            "original_source_rail_segments_not_in_connector_graph":firstpass.get(
                "source_unplaceable_segment_"+name,0),
            "physical_connector_nodes_n":len(uf.parent),
            "physical_connected_components_n":uf.components,
            "components_with_at_least_one_source_rail_segment":len(arr),
            "source_connected_components_with_one_original_rail_segment_n":sum(a[0]==1 for a in arr),
            "source_rail_segments_placed_on_physical_connector_components_n":n_placed,
            "source_total_mapped_rail_centerline_length_in_components_km":round(total_length,3),
            "source_top_component_segment_fraction":top[0][0]/n_placed if n_placed else None,
            "source_top_10_component_segment_fraction":sum(a[0] for a in top)/n_placed if n_placed else None,
            "source_top_component_centreline_km":round(top[0][1],3) if top else None,
            "source_top_10_component_centreline_km":round(sum(a[1] for a in top),3),
            "top_10_anonymous_component_ranks":[{
                "rank":i+1,"source_rail_segments":c[0],
                "source_rail_centerline_km":round(c[1],3),
                "source_segments_with_primary_name":c[2]}
                for i,c in enumerate(top)],
            "source_physical_components_NOT_named_railway_infrastructure_lines":True}
        if name=="all_original_rail_subtypes":
            # Lost segments having no original source connector are excluded
            # from component length, so total length cannot match all-source
            # geodesic length; validate placed length <= total original.
            if total_length>all_original_km*1.005:
                raise ValueError("Original rail graph has inflated more km than full Overture source")
    return out

def run(source_dir,census_path):
    census=json.loads(Path(census_path).read_text(encoding="utf-8"))
    paths=validate_files(source_dir,census)
    graph,first=topology(paths)
    if first["source_all_rail_segments"]!=census["source_original_rail_features_entire_release_n"]:
        raise ValueError("Input original global railway segments no longer match source-complete v357 census")
    components=aggregate_components(paths,census,graph,first)
    return {"version":"v361-full-authenticated-global-rail-physical-connector-components",
      "original_Overture_release":SOURCE_RELEASE,
      "original_128_source_parquet_file_census_sha256":EXPECTED_INVENTORY_SHA,
      "original_geodata_shards_n":len(paths),
      "source_full_original_rail_count_confirmed_against_v357":True,
      "original_v357_total_all_rail_subtype_segments":first["source_all_rail_segments"],
      "original_v357_centerline_km_all_rail_subtype":census["source_real_geodesic_all_rail_km"],
      "source_segment_topology_qa":first,
      "rail_physical_connector_graph_panels":components,
      "connected_components_do_not_equal_original_named_operating_railway_lines":True,
      "original_route_railway_relation_membership_NOT_INCLUDED_in_source_Oveture":True,
      "no_individual_original_source_GERS_connector_id_way_id_railgeometry_or_species_record_republished":True,
      "no_global_railway_plant_beta_diversity_or_seed_transport_causal_result":True,
      "original_source_ODbL":"Overture Maps 2026-11-18.0 / OpenStreetMap contributors"}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",required=True)
    p.add_argument("--census",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=run(a.source_dir,a.census)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!="rail_physical_connector_graph_panels"},
                     ensure_ascii=False,indent=2),flush=True)
    print(json.dumps(result["rail_physical_connector_graph_panels"],ensure_ascii=False),flush=True)

if __name__=="__main__":main()

#!/usr/bin/env python3
"""All seven checksum-pinned GRIP4 regional road vectors, same 100km clip and 200km rollup as rail; frozen plant candidate cells only."""
import argparse,hashlib,json,math,subprocess,zipfile
from collections import Counter
from pathlib import Path
import fiona
from pyproj import Transformer,Geod
from shapely.geometry import LineString
from scan_global_grip4_road_crossings import MD5S,digest,key_sha
from aggregate_entire_world_GERS_rail_100km_200km_grid_v370 import clip_source_line_to_tiles
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START,X_END,Y_END
W=200000
def candidates(path):
    h=json.loads(Path(path).read_text())
    if h["version"]!="railstatus-annual-verified-plant-size-v1" or h["width"]!=W or not h["exact_annual_counts_and_2025_moments_replay"]:raise ValueError("Wrong annual source")
    rows=h["grids"]["2025"];keys={tuple(r["key"]) for r in rows}
    if len(keys)!=len(rows) or len(keys)!=1545:raise ValueError("Frozen 2025 candidate geography changed")
    return keys
def possible_cells(xx,yy,keys):
    i,j=math.floor((min(xx)-X_START)/W),math.floor((min(yy)-Y_START)/W)
    a,b=math.floor((max(xx)-X_START)/W),math.floor((max(yy)-Y_START)/W)
    if (a-i+1)*(b-j+1)>len(keys):
        return any(i<=x<=a and j<=y<=b for x,y in keys)
    return any((x,y) in keys for x in range(i,a+1) for y in range(j,b+1))
def part_km(part,keys,to_ll,forward,inverse,geod,qa):
    if len(part)<2:qa["short_parts"]+=1;return {}
    xx,yy=zip(*[(v[0],v[1]) for v in part])
    lon,lat=to_ll.transform(xx,yy)
    if not all(math.isfinite(v) for v in (*lon,*lat)):raise ValueError("Nonfinite WGS84 source geometry")
    if max(lon)-min(lon)>180:qa["dateline_parts_omitted_like_rail"]+=1;return {}
    x,y=forward.transform(lon,lat)
    if not all(math.isfinite(v) for v in (*x,*y)):raise ValueError("Nonfinite road projection")
    if not possible_cells(x,y,keys):return {}
    qa["candidate_bbox_parts"]+=1
    geom=LineString(list(zip(lon,lat)));L=abs(float(geod.geometry_length(geom)))
    if L==0:qa["zero_length_parts"]+=1;return {}
    def clipped(g,l):
        z,omitted,status=clip_source_line_to_tiles(g,l,forward,inverse,geod)
        if omitted:
            # Rail also omits out-of-frame lines. If it intersects the
            # candidate bounding box, explicitly count the omitted amount.
            qa["outside_frame_candidate_bbox_parts"]+=1
            qa["outside_frame_candidate_bbox_m"]+=omitted
        return z
    try:z=clipped(geom,L)
    except ValueError:
        # Long or multipart source lines can be unlike short rail segments.
        # Deterministic subdivision at original vertices, never by outcome.
        qa["original_vertex_subdivision_parts"]+=1;z=Counter()
        for a,b in zip(geom.coords,list(geom.coords)[1:]):
            seg=LineString((a,b));l=abs(float(geod.geometry_length(seg)))
            if l>0:z.update(clipped(seg,l))
    ans=Counter()
    for (i,j),m in z.items():
        k=i//2,j//2
        if k in keys:ans[k]+=m/1000
    return ans
def run(region,keys,work):
    work.mkdir(parents=True,exist_ok=True);name=f"GRIP4_Region{region}_vector_fgdb.zip";z=work/name
    url=f"https://zenodo.org/records/6420961/files/{name}?download=1"
    subprocess.run(["curl","-LfsS","--retry","3","--retry-all-errors","--connect-timeout","30","--max-time","3600",url,"-o",str(z)],check=True)
    if digest(z,"md5")!=MD5S[region]:raise ValueError("Source archive checksum mismatch")
    sha=digest(z,"sha256");nbytes=z.stat().st_size
    with zipfile.ZipFile(z) as f:
        for member in f.namelist():
            p=Path(member)
            if p.is_absolute() or ".." in p.parts:raise ValueError("Unsafe archive path")
        f.extractall(work/"raw")
    dbs=list((work/"raw").rglob("*.gdb"))
    if len(dbs)!=1:raise ValueError("Expected one GDB")
    layers=[]
    for layer in fiona.listlayers(dbs[0]):
        with fiona.open(dbs[0],layer=layer,allow_unsupported_drivers=True) as src:
            for k in src.schema["properties"]:
                if k.upper()=="GP_RTP":layers.append((layer,k))
    if len(layers)!=1:raise ValueError("Ambiguous original road layer")
    layer,field=layers[0];qa=Counter();km=Counter();class_n=Counter()
    forward=Transformer.from_crs(4326,6933,always_xy=True);inverse=Transformer.from_crs(6933,4326,always_xy=True);geod=Geod(ellps="WGS84")
    with fiona.open(dbs[0],layer=layer,allow_unsupported_drivers=True) as src:
        total=len(src)
        if not src.crs_wkt or total<=0:raise ValueError("Missing original road CRS/count")
        to_ll=Transformer.from_crs(src.crs_wkt,4326,always_xy=True)
        for f in src:
            qa["features"]+=1;class_n[str(f["properties"][field])]+=1
            g=f["geometry"]
            if g is None:qa["null_geometry"]+=1;continue
            if g["type"]=="LineString":parts=[g["coordinates"]]
            elif g["type"]=="MultiLineString":parts=g["coordinates"]
            else:raise ValueError("Unexpected source road geometry type "+g["type"])
            for part in parts:km.update(part_km(part,keys,to_ll,forward,inverse,geod,qa))
            if qa["features"]%500000==0:print(json.dumps({"region":region,"qa":dict(qa),"candidate_km":sum(km.values())}),flush=True)
        if qa["features"]!=total or sum(class_n.values())!=total:raise ValueError("Original regional feature census incomplete")
    return {"version":"grip4-actual-road-vector-length-v1","region":region,"source_url":url,"source_md5":MD5S[region],"source_sha256":sha,"source_bytes":nbytes,"layer":layer,"source_features":total,"classes":dict(class_n),"candidate_sha256":key_sha(keys),"candidate_n":len(keys),"qa":dict(qa),"grids":[{"key":k,"km":float(km.get(k,0.))} for k in sorted(keys)]}
def combine(folder,keys):
    km=Counter();receipts=[]
    for region in range(1,8):
        paths=list(Path(folder).rglob(f"road_lengths_{region}.json"))
        if len(paths)!=1:raise ValueError("Missing/duplicate regional receipt")
        d=json.loads(paths[0].read_text())
        if d["version"]!="grip4-actual-road-vector-length-v1" or d["region"]!=region or d["source_md5"]!=MD5S[region] or d["candidate_sha256"]!=key_sha(keys) or d["candidate_n"]!=len(keys) or d["qa"]["features"]!=d["source_features"]:raise ValueError("Regional source receipt mismatch")
        grid={tuple(r["key"]):float(r["km"]) for r in d["grids"]}
        if len(grid)!=len(d["grids"]) or set(grid)!=keys or any(not math.isfinite(v) or v<0 for v in grid.values()):raise ValueError("Invalid regional road grid")
        km.update(grid);receipts.append({k:v for k,v in d.items() if k!="grids"})
    return dict(km),receipts
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--region",type=int,choices=range(1,8),required=True)
    for n in ("annual-size","work","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();d=run(a.region,candidates(a.annual_size),a.work)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,allow_nan=False))
    print(json.dumps({**{k:v for k,v in d.items() if k!="grids"},"candidate_km":sum(x["km"] for x in d["grids"])},allow_nan=False),flush=True)
if __name__=="__main__":main()

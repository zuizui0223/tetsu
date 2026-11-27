#!/usr/bin/env python3
"""Original GRIP4 regional road geometries on frozen bird candidate boundaries."""
import argparse,hashlib,json,math,subprocess,zipfile
from pathlib import Path
from collections import Counter
import fiona
from pyproj import Transformer
from aggregate_world_rail_boundary_crossings import crossing_edges,tile
MD5S={1:"18cb99e5e42ee224b6d6a761427c3fc1",2:"c9bed54134074a5818e2177eab8727ea",3:"6e56b32671e2945ee18deeba9da6596b",4:"2ae428508708db7a60ca192994186857",5:"426c6da68039e3f23cdce91b68803244",6:"a07ab65dfca03170619b86fc4b65620f",7:"5f99a4b3824937277fe663ad3fa1282e"}
W=200000
def candidates(path):
    h=json.loads(Path(path).read_text())
    if h["version"]!="railstatus-avonet-heldback-2018-22-hwi-pairs-v1" or not h["exact_2025_full_and_joint_guild_replay"] or h["width"]!=W:raise ValueError("Wrong frozen source")
    keys={(tuple(p["a"]),tuple(p["b"])) for p in h["pairs"]["2025"]}
    if len(keys)!=1645 or any(sum(abs(x-y) for x,y in zip(a,b))!=1 for a,b in keys):raise ValueError("Candidate geometry mismatch")
    return keys
def key_sha(keys):
    return hashlib.sha256(json.dumps(sorted(keys),separators=(",",":")).encode()).hexdigest()
def digest(path,algorithm):
    h=hashlib.new(algorithm)
    with open(path,"rb") as f:
        for x in iter(lambda:f.read(4*1024*1024),b""):h.update(x)
    return h.hexdigest()
def feature_edges(g,to_ll,to_eq,qa):
    if g["type"]=="LineString":parts=[g["coordinates"]]
    elif g["type"]=="MultiLineString":parts=g["coordinates"]
    else:raise ValueError("Unexpected road geometry "+g["type"])
    result=set()
    for part in parts:
        if len(part)<2:continue
        xx,yy=zip(*[(p[0],p[1]) for p in part])
        lon,lat=to_ll.transform(xx,yy)
        if max(lon)-min(lon)>180:qa["dateline_parts_skipped"]+=1;continue
        x,y=to_eq.transform(lon,lat)
        if not all(math.isfinite(v) for v in (*x,*y)):raise ValueError("Nonfinite geometry")
        if tile(min(x),min(y),W)==tile(max(x),max(y),W) and tile(min(x),min(y),W) is not None:continue
        result.update(crossing_edges(list(zip(x,y)),W))
    return result
def run(region,keys,work):
    work.mkdir(parents=True,exist_ok=True)
    name=f"GRIP4_Region{region}_vector_fgdb.zip";z=work/name
    url=f"https://zenodo.org/records/6420961/files/{name}?download=1"
    subprocess.run(["curl","-LfsS","--retry","3","--retry-all-errors","--connect-timeout","30","--max-time","3600",url,"-o",str(z)],check=True)
    if digest(z,"md5")!=MD5S[region]:raise ValueError("Original archive MD5 mismatch")
    sha=digest(z,"sha256");nbytes=z.stat().st_size
    with zipfile.ZipFile(z) as f:
        for member in f.namelist():
            p=Path(member)
            if p.is_absolute() or ".." in p.parts:raise ValueError("Unsafe archive path")
        f.extractall(work/"raw")
    dbs=list((work/"raw").rglob("*.gdb"))
    if len(dbs)!=1:raise ValueError("Expected one original GDB")
    layers=[]
    for layer in fiona.listlayers(dbs[0]):
        with fiona.open(dbs[0],layer=layer,allow_unsupported_drivers=True) as src:
            for k,v in src.schema["properties"].items():
                if k.upper()=="GP_RTP":layers.append((layer,k,v))
    if len(layers)!=1:raise ValueError(f"Ambiguous road layer {layers}")
    layer,field,typ=layers[0];qa=Counter();edges={1:Counter(),2:Counter()}
    with fiona.open(dbs[0],layer=layer,allow_unsupported_drivers=True) as src:
        total=len(src)
        if not src.crs_wkt or total<=0:raise ValueError("Missing source CRS/count")
        to_ll=Transformer.from_crs(src.crs_wkt,4326,always_xy=True);to_eq=Transformer.from_crs(4326,6933,always_xy=True)
        clause=f"{field} IN ('1','2')" if typ.startswith("str") else f"{field} IN (1,2)"
        for feature in src.filter(where=clause):
            t=int(feature["properties"][field])
            if t not in edges:raise ValueError("Source attribute filter failed")
            qa[f"type_{t}_features"]+=1
            if feature["geometry"] is None:qa[f"type_{t}_null"]+=1;continue
            edges[t].update(feature_edges(feature["geometry"],to_ll,to_eq,qa)&keys)
        if any(qa[f"type_{t}_features"]<=0 for t in edges) or sum(qa[f"type_{t}_features"] for t in edges)>total:raise ValueError("Implausible road source census")
    return {"version":"grip4-original-vector-candidate-crossings-v1","region":region,"url":url,"source_md5":MD5S[region],"source_sha256":sha,"source_bytes":nbytes,"layer":layer,"total_source_features":total,"candidate_sha256":key_sha(keys),"n_candidates":len(keys),"qa":dict(qa),"edges":{str(t):[{"a":a,"b":b,"features":n} for (a,b),n in sorted(v.items())] for t,v in edges.items()}}
def combine(folder,keys):
    by={1:Counter(),2:Counter()};receipts=[]
    for region in range(1,8):
        paths=list(Path(folder).rglob(f"road_crossings_{region}.json"))
        if len(paths)!=1:raise ValueError("Missing/duplicate regional artifact")
        d=json.loads(paths[0].read_text())
        if d["version"]!="grip4-original-vector-candidate-crossings-v1" or d["region"]!=region or d["source_md5"]!=MD5S[region] or d["candidate_sha256"]!=key_sha(keys) or d["n_candidates"]!=len(keys):raise ValueError("Regional receipt mismatch")
        receipts.append({k:v for k,v in d.items() if k!="edges"})
        for t in by:
            seen=set()
            for z in d["edges"][str(t)]:
                k=tuple(z["a"]),tuple(z["b"])
                if k in seen or k not in keys or type(z["features"])!=int or z["features"]<=0:raise ValueError("Invalid road crossing row")
                seen.add(k);by[t][k]+=z["features"]
    return by,receipts
def main():
    p=argparse.ArgumentParser();p.add_argument("--region",type=int,choices=range(1,8),required=True);p.add_argument("--historical",required=True);p.add_argument("--work",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();d=run(a.region,candidates(a.historical),a.work);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,allow_nan=False));print(json.dumps({**{k:v for k,v in d.items() if k!="edges"},"crossed_edges":{k:len(v) for k,v in d["edges"].items()}}),flush=True)
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Independent satellite test of mapped railway-status habitat structure.

Draw inactive railway-km points, pair with nearby mapped active surface
lines using fixed geometry-only rules. Analyze ESA WorldCover 2021 v200
10m class fractions, not species height, biomass or actual rail closure.
No iNaturalist queries. No site GPS or raw rasters exported.
"""
import argparse,json,hashlib,math,random
from pathlib import Path
from collections import Counter,defaultdict,OrderedDict
import numpy as np
from pyproj import Geod,Transformer
from shapely.ops import transform,unary_union
from shapely import points,distance
from shapely.geometry import box
import rasterio
from rasterio.windows import from_bounds,Window
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import audit_world_railstatus_local_pair_support as a

DRAWS=8192
SEED="tetsu-railstatus-worldcover-20261127-v1"
CAP=256
NEAR=80.
BG_INNER=150.
BG_OUTER=300.
BG_RAIL_EXCLUDE=80.
MIN_BG_AREA_FRACTION=.5
MIN_VALID_FRACTION=.95
MIN_PAIRS=50
MIN_BLOCKS=20
BUCKET="esa-worldcover"
PREFIX="v200/2021/map/"
CLASS={10:"tree",20:"shrub",30:"grass",40:"crop",50:"built",60:"bare",70:"snow",80:"water",90:"wetland",95:"mangrove",100:"moss_lichen"}
def frozen_points(paths,status):
    # Reuse exact verified geometry rules, but new independently seeded
    # length draw for the distinct satellite endpoint. No observation gate.
    lengths,cq=a.census(paths,status)
    a.N_DRAWS=DRAWS;a.SEED=SEED
    sites=a.sample(paths,lengths[a.INACTIVE]);local,lq=a.local_geometry(paths,sites)
    qa=Counter();pairs=[]
    for s in sites:
        if a.window(s) is None:qa["no_surface_or_geographic_window"]+=1;continue
        p,reason=a.pair(s,local.get(s["draw"],[]));qa[reason]+=1
        if p is not None:pairs.append(p)
    tf=Transformer.from_crs(4326,6933,always_xy=True);geod=Geod(ellps="WGS84")
    selected=[];seen=set()
    for p in pairs:
        k=a.source_coarse_cell(*p["inactive"],tf,100000)
        if k is None or k in seen:continue
        if any(min(abs(geod.inv(*u,*v)[2]) for u in (p["active"],p["inactive"]) for v in (q["active"],q["inactive"]))<1000 for q in selected):continue
        seen.add(k);selected.append(dict(p,block100=k))
    used=selected[:CAP]
    return used,{p["draw"]:local[p["draw"]] for p in used},{"draws":DRAWS,"census":cq,"length_km":{k:v/1000 for k,v in lengths.items()},"local_pass":lq,"pair_reasons":dict(qa),"n_spatially_separated":len(selected),"n_selected":len(used)}

def tile_name(south,west):
    if south%3 or west%3 or not -90<=south<90 or not -180<=west<180:raise ValueError("Invalid WorldCover tile origin")
    return f"{'N' if south>=0 else 'S'}{abs(south):02d}{'E' if west>=0 else 'W'}{abs(west):03d}"
def tile_keys(bounds):
    x0,y0,x1,y1=bounds
    for south in range(math.floor(y0/3)*3,math.floor(y1/3)*3+1,3):
        for west in range(math.floor(x0/3)*3,math.floor(x1/3)*3+1,3):
            yield south,west,PREFIX+"ESA_WorldCover_10m_2021_v200_"+tile_name(south,west)+"_Map.tif"

class Tiles:
    def __init__(self):
        self.client=boto3.client("s3",region_name="eu-central-1",config=Config(signature_version=UNSIGNED,connect_timeout=15,read_timeout=30,retries={"max_attempts":2}))
        self.cache=OrderedDict();self.receipts={}
    def close(self):
        for ds in self.cache.values():ds.close()
        self.cache.clear()
    def open(self,south,west,key):
        if key in self.cache:
            self.cache.move_to_end(key);return self.cache[key]
        h=self.client.head_object(Bucket=BUCKET,Key=key)
        if h["ContentLength"]<=0 or not h.get("ETag"):raise ValueError("Missing source receipt")
        receipt={"bytes":h["ContentLength"],"etag":h["ETag"],"last_modified":h["LastModified"].isoformat()}
        if key in self.receipts and any(self.receipts[key][k]!=v for k,v in receipt.items()):raise ValueError("Source changed during run")
        url="https://esa-worldcover.s3.eu-central-1.amazonaws.com/"+key
        ds=rasterio.open(url)
        try:
            if ds.crs.to_epsg()!=4326 or ds.count!=1 or ds.dtypes[0]!="uint8" or ds.width!=36000 or ds.height!=36000:raise ValueError("Unexpected WorldCover tile schema")
            if max(abs(a-b) for a,b in zip(tuple(ds.bounds),(west,south,west+3,south+3)))>1e-7 or abs(ds.transform.a-1/12000)>1e-10 or abs(ds.transform.e+1/12000)>1e-10:raise ValueError("Unexpected source tile bounds/resolution")
            receipt.update(shape=list(ds.shape),transform=list(ds.transform),nodata=ds.nodata)
        except Exception:
            ds.close();raise
        self.receipts[key]=receipt;self.cache[key]=ds
        if len(self.cache)>8:self.cache.popitem(last=False)[1].close()
        return ds

def fractions(v,weights,mask):
    n=int(mask.sum())
    if not n:return None
    selected=v[mask];w=weights[mask]
    valid=np.isin(selected,list(CLASS))
    validfrac=float(w[valid].sum()/w.sum())
    if validfrac<MIN_VALID_FRACTION:return None
    total=w[valid].sum()
    out={name:float(w[(selected==k)].sum()/total) for k,name in CLASS.items()}
    if abs(sum(out.values())-1)>1e-10:raise ValueError("Class fractions do not sum to one")
    return {"n_pixel_centers":n,"valid_fraction":validfrac,"fractions":out}

def site(tiles,xy,lines):
    lon,lat=xy;fwd=Transformer.from_crs(4326,f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs",always_xy=True)
    dy=.004;dx=dy/math.cos(math.radians(abs(lat)+dy))
    bounds=(lon-dx,lat-dy,lon+dx,lat+dy)
    # Exclude any mapped railway from the local background, including
    # disused, active, bridge/tunnel and unclassified segments.
    local=[]
    clip=box(*bounds)
    for _,g in lines:
        if not g.intersects(clip):continue
        g=g.intersection(clip)
        if not g.is_empty:local.append(transform(fwd.transform,g))
    if not local:raise ValueError("Rail point has no local geometry")
    rail=unary_union(local)
    vv=[];ww=[];nn=[];bb=[];base=[]
    for south,west,key in tile_keys(bounds):
        ds=tiles.open(south,west,key)
        left=max(bounds[0],ds.bounds.left);bottom=max(bounds[1],ds.bounds.bottom)
        right=min(bounds[2],ds.bounds.right);top=min(bounds[3],ds.bounds.top)
        if not left<right or not bottom<top:continue
        raw=from_bounds(left,bottom,right,top,ds.transform)
        i0=max(0,math.floor(raw.col_off));j0=max(0,math.floor(raw.row_off))
        i1=min(ds.width,math.ceil(raw.col_off+raw.width));j1=min(ds.height,math.ceil(raw.row_off+raw.height))
        if i1<=i0 or j1<=j0:continue
        win=Window(i0,j0,i1-i0,j1-j0)
        v=ds.read(1,window=win);tr=ds.window_transform(win)
        cy,cx=np.indices(v.shape)
        lo=tr.c+(cx+.5)*tr.a;la=tr.f+(cy+.5)*tr.e
        x,y=fwd.transform(lo,la);rad=np.hypot(x,y)
        near=rad<=NEAR;ann=(rad>=BG_INNER)&(rad<=BG_OUTER)
        bg=ann.copy()
        if ann.any():bg[ann]=distance(points(x[ann],y[ann]),rail)>=BG_RAIL_EXCLUDE
        # EPSG:4326 pixel area is proportional to cos(latitude) over
        # these small windows; consistent area weights, not resampled classes.
        weights=np.cos(np.radians(la))
        vv.append(v.ravel());ww.append(weights.ravel());nn.append(near.ravel());bb.append(bg.ravel());base.append(ann.ravel())
    if not vv:return {"status":"no_source_pixels"}
    v,w,n,bg,ann=map(np.concatenate,(vv,ww,nn,bb,base))
    if n.sum()<30 or ann.sum()<100:return {"status":"insufficient_geometry_pixels"}
    area=float(w[bg].sum()/w[ann].sum())
    if area<MIN_BG_AREA_FRACTION:return {"status":"insufficient_nonrail_background","background_area_fraction":area}
    near=fractions(v,w,n);background=fractions(v,w,bg)
    if near is None or background is None:return {"status":"insufficient_valid_classes","background_area_fraction":area}
    return {"status":"complete","near":near,"background":background,"background_area_fraction":area}

def interval(values):
    # Equal 1000km blocks; approximate normal cluster interval for the
    # finite conditional sample. No pixel-level pseudoreplication.
    blocks=defaultdict(list)
    for k,v in values:blocks[k].append(v)
    n=len(blocks);means=np.asarray([np.mean(v) for v in blocks.values()])
    if n<MIN_BLOCKS:return None
    mu=float(means.mean());se=float(means.std(ddof=1)/math.sqrt(n))
    return {"estimate":mu,"interval95_approx":[mu-1.96*se,mu+1.96*se],"n_pairs":len(values),"n_blocks":n}

def evaluate(rows):
    rr=[r for r in rows if all(r[k]["status"]=="complete" for k in ("active","inactive"))]
    nblocks=len({tuple(r["block1000"]) for r in rr})
    qa={"n_complete_pairs":len(rr),"n_blocks":nblocks}
    if len(rr)<MIN_PAIRS or nblocks<MIN_BLOCKS:return dict(qa,status="insufficient support")
    result={}
    for name in CLASS.values():
        z={}
        for label in ("near_difference","background_difference","difference_in_differences"):
            values=[]
            for r in rr:
                a0=r["active"]["near"]["fractions"][name]-r["inactive"]["near"]["fractions"][name]
                c0=r["active"]["background"]["fractions"][name]-r["inactive"]["background"]["fractions"][name]
                v=a0 if label=="near_difference" else c0 if label=="background_difference" else a0-c0
                values.append((tuple(r["block1000"]),v))
            z[label]=interval(values)
        result[name]=z
    return dict(qa,status="estimated",effects=result)

def main():
    p=argparse.ArgumentParser()
    for k in ("rail-dir","census","status","out"):p.add_argument("--"+k,type=Path,required=True)
    args=p.parse_args();status=json.loads(args.status.read_text())
    if status["source_Overture_release"]!=a.SOURCE_RELEASE or not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong rail source")
    paths=a.validate_files(args.rail_dir,json.loads(args.census.read_text()))
    pairs,local,geo=frozen_points(paths,status)
    print(json.dumps({"worldcover_geometry":geo},allow_nan=False),flush=True)
    tiles=Tiles();rows=[]
    try:
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_HTTP_TIMEOUT="45",GDAL_HTTP_MAX_RETRY="1",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
            for p in pairs:
                r={"draw":p["draw"],"block100":p["block100"],"block1000":[k//10 for k in p["block100"]],"separation_m":round(p["distance_m"],1)}
                for k in ("active","inactive"):
                    try:r[k]=site(tiles,p[k],local[p["draw"]])
                    except Exception as e:r[k]={"status":"source_or_geometry_error","error_type":type(e).__name__}
                rows.append(r)
                if len(rows)%16==0:print(json.dumps({"worldcover_progress":{"pairs":len(rows),"complete":sum(all(r[k]["status"]=="complete" for k in ("active","inactive")) for r in rows),"tiles":len(tiles.receipts)}}),flush=True)
    finally:tiles.close()
    ans=evaluate(rows)
    out={"version":"railstatus-local-worldcover-2021-v1","source":"ESA WorldCover 2021 v200","DOI":"10.5281/zenodo.7254221","geometry":geo,"tile_receipts":tiles.receipts,
        "status_counts":{k:dict(Counter(r[k]["status"] for r in rows)) for k in ("active","inactive")},
        "primary":ans,"pairs":rows,"limits":["2021 satellite product vs 2026 mapped rail flags; no verified operation/closure dates or changes through time",
        "near disk 80m and nearby 150–300m background excluding mapped rail within 80m; different habitat and built use may remain",
        "mapped tree cover is not canopy height, biomass, understory height, actual grass mowing or species composition",
        "conditional sample of geographically isolated nearby status pairs; not representative of every railway kilometre",
        "satellite land-cover classification error and mapping completeness; source class fractions not ground vegetation plots",
        "fixed minimum support, no iNat or eBird, no location or outcome-driven replacement"]}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({"worldcover_result":ans,"status_counts":out["status_counts"],"n_tiles":len(tiles.receipts)},allow_nan=False),flush=True)
if __name__=="__main__":main()

"""Independent GLC_FCS30D v2 before/after candidate-closure temporal pilot.

No outcome-dependent site selection, no actual freight/operation cessation claim.
Reuse frozen N05/Overture geography; no raw coordinates or pixels exported.
"""
import argparse
import collections
import io
import json
import math
import pathlib
import re
import urllib.parse
from functools import lru_cache
import numpy as np
import pyarrow.parquet as pq
import rasterio
from rasterio.windows import Window,from_bounds
from scipy.stats import t
from shapely.geometry import box
from shapely import points,distance
from shapely.ops import transform,unary_union
from shapely.strtree import STRtree
from pyproj import Transformer
import audit_japan_terminal_rail_geometry as j
import audit_japan_terminal_overture as o
import analyze_japan_terminal_worldcover as old
from audit_dated_rail_closure_sources import get
from build_global_Overture_physical_connector_components_v361 import validate_files

BASE="https://data.source.coop/portolan-mirrors/glc-fcs30d/glc_fcs30d/"
YEARS={1985,1990,1995,*range(2000,2023)}
FOREST={51,52,61,62,71,72,81,82,91,92}
BUILT={190}
VALID={10,11,12,20,*FOREST,120,121,122,130,140,150,152,153,181,182,183,184,185,186,187,190,200,201,202,210,220}
MIN_GROUPS=25
MIN_BLOCKS=8
NEAR=120.
INNER=240.
OUTER=480.
EXCLUDE=120.
OUTCOMES={"forest":FOREST,"impervious":BUILT}
def epochs(y):
    pre=[x for x in YEARS if x<y and x-5 in YEARS]
    post=[x for x in YEARS if x>y and x+5 in YEARS]
    if not pre or not post:
        return None
    p,q=max(pre),min(post)
    return p-5,p,q,q+5
def catalog():
    raw,receipt=get(BASE+"items.parquet",80_000_000)
    table=pq.read_table(io.BytesIO(raw),columns=["id","assets"])
    idx={}
    for r in table.to_pylist():
        m=re.fullmatch(r"(\d{4})-([EW])(\d+)([NS])(\d+)",r["id"])
        if not m:
            raise ValueError("unexpected catalog item ID")
        yr,ew,lon,ns,lat=m.groups()
        west=int(lon)*(1 if ew=="E" else -1)
        north=int(lat)*(1 if ns=="N" else -1)
        if not (120<=west<=155 and 20<=north<=50):
            continue
        assets=r["assets"]
        if isinstance(assets,str):assets=json.loads(assets)
        asset=assets["data"]
        url=asset["href"]
        u=urllib.parse.urlsplit(url)
        if u.scheme!="https" or u.hostname!="data.source.coop" or not url.startswith(BASE) or not u.path.endswith(".tif"):
            raise ValueError("unexpected COG href")
        k=(int(yr),west,north)
        if k in idx:
            raise ValueError("duplicate item")
        idx[k]={"href":url,"item_id":r["id"],"file_checksum":asset.get("file:checksum"),"file_size":asset.get("file:size")}
    if len(table)!=24986 or not idx:
        raise ValueError("catalog census mismatch")
    return idx,receipt
class Tiles:
    def __init__(self,idx):
        self.idx=idx
        self.opened=collections.OrderedDict()
        self.receipts={}
    def open(self,yr,w,n):
        k=(yr,w,n)
        if k in self.opened:
            self.opened.move_to_end(k)
            return self.opened[k]
        item=self.idx[k]
        ds=rasterio.open(item["href"])
        try:
            if ds.crs.to_epsg()!=4326 or ds.count!=1 or ds.dtypes[0]!="uint8" or ds.nodata!=0:
                raise ValueError("unexpected raster schema")
            if ds.width!=18554 or ds.height!=18554 or abs(ds.transform.a-.000269494585236)>1e-10 or abs(ds.transform.e+.000269494585236)>1e-10:
                raise ValueError("unexpected raster grid")
            if abs(ds.bounds.left-w)>.001 or abs(ds.bounds.top-n)>.001:
                raise ValueError("tile origin mismatch")
            self.receipts[item["item_id"]]={**item,"shape":list(ds.shape),"transform":list(ds.transform)}
            self.opened[k]=ds
            if len(self.opened)>32:self.opened.popitem(last=False)[1].close()
        except Exception:
            ds.close()
            raise
        return ds
    def close(self):
        for ds in self.opened.values():ds.close()
        self.opened.clear()

def windows(tiles,ep,lon,lat,lines):
    tf=Transformer.from_crs(4326,f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs",always_xy=True)
    dy=.006;dx=dy/math.cos(math.radians(abs(lat)+dy))
    bounds=(lon-dx,lat-dy,lon+dx,lat+dy)
    rail=unary_union([transform(tf.transform,g) for _,g in lines])
    if rail.is_empty:raise ValueError("missing historical line")
    data=[];ww=[];nn=[];bb=[];aa=[]
    for w in range(math.floor(bounds[0]/5)*5,math.floor(bounds[2]/5)*5+1,5):
        for n in range((math.floor(bounds[1]/5)+1)*5,(math.floor(bounds[3]/5)+1)*5+1,5):
            # Nominal 5-degree tile mask avoids duplicate edge pixels.
            left=max(bounds[0],w);right=min(bounds[2],w+5)
            bottom=max(bounds[1],n-5);top=min(bounds[3],n)
            if not left<right or not bottom<top:continue
            ds=tiles.open(ep[0],w,n)
            raw=from_bounds(left,bottom,right,top,ds.transform)
            i=max(0,math.floor(raw.col_off));k=max(0,math.floor(raw.row_off))
            i1=min(ds.width,math.ceil(raw.col_off+raw.width));k1=min(ds.height,math.ceil(raw.row_off+raw.height))
            win=Window(i,k,i1-i,k1-k)
            tr=ds.window_transform(win)
            v=[]
            for yr in ep:
                dd=tiles.open(yr,w,n)
                if dd.transform!=ds.transform or dd.shape!=ds.shape:
                    raise ValueError("temporal pixel alignment changed")
                v.append(dd.read(1,window=win))
            iy,ix=np.indices(v[0].shape)
            lo=tr.c+(ix+.5)*tr.a;la=tr.f+(iy+.5)*tr.e
            tile=(lo>=w)&(lo<w+5)&(la>n-5)&(la<=n)
            x,y=tf.transform(lo,la)
            rad=np.hypot(x,y)
            near=tile&(rad<=NEAR)
            ann=tile&(rad>=INNER)&(rad<=OUTER)
            bg=ann.copy()
            if ann.any():bg[ann]=distance(points(x[ann],y[ann]),rail)>=EXCLUDE
            data.append(np.stack(v).reshape(4,-1))
            ww.append(np.cos(np.radians(la)).ravel())
            nn.append(near.ravel());bb.append(bg.ravel());aa.append(ann.ravel())
    if not data:return None,"no_pixels"
    z=np.concatenate(data,axis=1);weight=np.concatenate(ww)
    near,bg,ann=map(np.concatenate,(nn,bb,aa))
    if near.sum()<30 or ann.sum()<100:return None,"insufficient_pixel_geometry"
    if weight[bg].sum()<.5*weight[ann].sum():return None,"insufficient_background"
    valid=np.all(np.isin(z,list(VALID)),axis=0)
    d={}
    for label,mask in (("near",near),("background",bg)):
        if weight[valid&mask].sum()<.95*weight[mask].sum():return None,"invalid_temporal_pixels"
        zz=z[:,valid&mask];w=weight[valid&mask];den=w.sum()
        d[label]={}
        for name,codes in OUTCOMES.items():
            b=np.isin(zz,list(codes))
            def freq(a):return float(w[a].sum()/den)
            pre=freq(b[1])-freq(b[0])
            post=freq(b[3])-freq(b[2])
            span=(freq(b[3])-freq(b[1]))*5/(ep[3]-ep[1])
            gap=(freq(b[2])-freq(b[1]))*5/(ep[2]-ep[1])
            d[label][name]={"pre5":pre,"post5":post,"span_per5":span,
                "span_minus_pre":span-pre,"post_minus_pre":post-pre,"gap_per5":gap,
                "span_gross_gain_per5":freq(~b[1]&b[3])*5/(ep[3]-ep[1]),
                "span_gross_loss_per5":freq(b[1]&~b[3])*5/(ep[3]-ep[1]),
                "pre1_cover":freq(b[1])}
    return d,"complete"

def block_interval(rows,endpoint,kind,name):
    by=collections.defaultdict(list)
    for r in rows:
        if endpoint=="delta":
            v=r["near"][name][kind]-r["background"][name][kind]
        else:v=r[endpoint][name][kind]
        by[r["block"]].append(v)
    G=len(by)
    if len(rows)<MIN_GROUPS or G<MIN_BLOCKS:
        return {"status":"insufficient support","n_groups":len(rows),"n_blocks":G}
    v=np.array([np.mean(a) for a in by.values()])
    mu=float(v.mean());se=float(v.std(ddof=1)/math.sqrt(G))
    crit=float(t.ppf(.975,G-1))
    return {"status":"estimated","n_groups":len(rows),"n_blocks":G,
        "estimate":mu,"interval95":[mu-crit*se,mu+crit*se]}
def summarize(rows):
    ans={}
    for name in OUTCOMES:
        ans[name]={}
        for e in ("delta","near","background"):
            ans[name][e]={k:block_interval(rows,e,k,name) for k in ("pre5","post5","span_per5","span_minus_pre","post_minus_pre","gap_per5","span_gross_gain_per5","span_gross_loss_per5")}
    # Descriptive divergent pathways on exact same support, not a significance rescue.
    cats=collections.Counter()
    for r in rows:
        f=r["near"]["forest"]["span_minus_pre"]-r["background"]["forest"]["span_minus_pre"]
        b=r["near"]["impervious"]["span_minus_pre"]-r["background"]["impervious"]["span_minus_pre"]
        label=("forest_only" if f>=.10 and b<.10 else "impervious_only" if b>=.10 and f<.10 else "both" if min(f,b)>=.10 else "neither")
        cats[label]+=1
    return {"effects":ans,"descriptive_groups_excess_at_least_10pp_per5":dict(cats)}
def selftest():
    assert epochs(1990) is None
    assert epochs(1991)==(1985,1990,1995,2000)
    assert epochs(2006)==(2000,2005,2007,2012)
    assert epochs(2015)==(2009,2014,2016,2021)
    assert epochs(2018) is None
    assert not (FOREST&BUILT) and 250 not in VALID
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    for k in ("rail-dir","census","out"):p.add_argument("--"+k,type=pathlib.Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    rr,geoms,receipt,_=j.read_n05()
    cand,attrs=j.classification(rr)
    pp=[];prior=j.audit_geometry(rr,geoms,cand,pp)
    if len(pp)!=4534 or len(cand)!=576 or prior["all"]["screened_length_proxy_km"]!=4388.978:
        raise ValueError("prior geometry replay failed")
    paths=validate_files(a.rail_dir,json.loads(a.census.read_text()))
    inactive,other,qa=o.split_current(paths)
    chosen,sel=old.select(old.clean_points(pp,inactive,other))
    if len(chosen)!=182 or sum(len(s) for k,y,s in chosen)!=541:
        raise ValueError("prior selected geometry replay failed")
    cases=[(k,y,ss,epochs(y)) for k,y,ss in chosen if epochs(y)]
    blocks={old.grid(ss[0][0][2],ss[0][0][3],200000) for k,y,ss,e in cases}
    result={"version":"japan-terminal-glcfcs30d-temporal-pilot-v1","n05_source_sha":receipt["sha256"],
        "previous_geometry_replay":True,"preoutcome_support":{"n_groups":len(cases),"n_blocks":len(blocks),
        "n_sites":sum(len(ss) for k,y,ss,e in cases),"epoch_patterns":dict(collections.Counter(str(e) for k,y,ss,e in cases))},
        "limits":["exploratory before/after temporal comparison, not verified railway operation cessation",
        "geometry selected using 2026 map; survivor/cohort selection and land-use conversion remain",
        "historical corridor position and remote-sensing class error; shared pixel mask; no raw source geodata exported",
        "near/background comparison and prior trend diagnostic do not replace an external active-line control",
        "net and gross land-cover changes not plant height, mowing, biomass, or species diversity",
        "two primary outcomes, forest and impervious; no selective threshold rescue; 30m pixel center area weights",
        "GLC_FCS30D v2, DOI 10.5281/zenodo.15063683, DOI 10.5194/essd-16-1353-2024, CC BY 4.0, Source Cooperative Portolan mirror"]}
    print(json.dumps({"preoutcome_support":result["preoutcome_support"]}),flush=True)
    if len(cases)<MIN_GROUPS or len(blocks)<MIN_BLOCKS:
        result["status"]="insufficient temporal support; no raster extraction"
    else:
        try:
            idx,source=catalog()
        except Exception as e:
            result.update(status="blocked catalog",error_type=type(e).__name__,reason=str(e)[:300])
            idx=None
        if idx:
            result["catalog_receipt"]=source
            historical=[transform(old.FWD.transform,g) for g in geoms if not g.is_empty and g.is_valid]
            geom=[historical,inactive,other];trees=[STRtree(v) for v in geom]
            ts=Tiles(idx);done=[];status=collections.Counter()
            try:
                with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_HTTP_TIMEOUT="30",GDAL_HTTP_MAX_RETRY="0",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_CACHEMAX=128_000_000):
                    for key,y,ss,ep in cases:
                        vals=[]
                        for p,z in ss:
                            x0,y0=p[2:4];lon,lat=old.INV.transform(x0,y0)
                            try:
                                v,s=windows(ts,ep,lon,lat,old.local_lines(trees,geom,x0,y0))
                            except Exception as e:
                                v,s=None,"source_or_geometry_error_"+type(e).__name__
                                if sum(a for k,a in status.items() if k.startswith("source_or_geometry_error"))<2:
                                    print(json.dumps({"source_error_type":type(e).__name__,"reason":str(e)[:240]}),flush=True)
                            status[s]+=1
                            if v is not None:vals.append(v)
                        if len(vals)==len(ss):
                            row={"block":old.grid(ss[0][0][2],ss[0][0][3],200000)}
                            for e in ("near","background"):
                                row[e]={name:{k:float(np.mean([v[e][name][k] for v in vals])) for k in vals[0][e][name]} for name in OUTCOMES}
                            done.append(row)
                        if len(done)%10==0 and len(done):
                            print(json.dumps({"sites":sum(status.values()),"groups":len(done)}),flush=True)
                        if sum(a for k,a in status.items() if k.startswith("source_or_geometry_error"))>=12 and status["complete"]==0:
                            result["early_stop"]="12 source/geometry errors, no complete sites"
                            break
            finally:ts.close()
            result.update(status="analyzed",site_status=dict(status),n_complete_groups=len(done),tile_receipts=ts.receipts,models=summarize(done))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps({k:v for k,v in result.items() if k!="tile_receipts"},ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()

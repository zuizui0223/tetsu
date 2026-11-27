"""Frozen US STB administrative completion land-cover change and pathway comparison."""
import argparse
import collections
import hashlib
import json
import math
import pathlib
import numpy as np
import rasterio
from scipy.stats import t
from scipy.optimize import linear_sum_assignment
from shapely.geometry import box
from shapely.ops import transform
from shapely.strtree import STRtree
from pyproj import Transformer
import audit_stb_dated_rail_cohort as s
import analyze_japan_terminal_glcfcs30d as glc
from build_global_Overture_physical_connector_components_v361 import validate_files

SEED="stb-glcfcs30d-20261127-v1"
EXPECTED={"abandoned":"e8eb0b33ba1f0be6c0f8176831f03648459b468891c83bc4699a85f8d734e73c",
          "railbanked":"8f9a7c70052fbcd60213bde0a872766c73dea56288337c6201ba47992b0eebc0"}
INV=Transformer.from_crs(5070,4326,always_xy=True)
KINDS=("abandoned","railbanked")
NAMES=("forest","impervious")
MIN_PAIRS=20
MIN_PAIR_BLOCKS=10
def block(p,width=200000):
    x,y=s.EQ.transform(p.x,p.y)
    return math.floor(x/width),math.floor(y/width)
def digest(*args):
    return hashlib.sha256("|".join(map(str,(SEED,*args))).encode()).digest()
def select(events,raw,active,inactive):
    at=STRtree(active)
    histories=[(key,g) for rr in raw.values() for key,y,g in rr if g is not None]
    ht=STRtree([g for k,g in histories])
    qa=collections.Counter()
    chosen=[];occupied=[]
    for kind,key,y,g in sorted(events,key=lambda r:digest(r[0],r[1])):
        ep=glc.epochs(y)
        if ep is None:continue
        good=[]
        for line in (list(g.geoms) if g.geom_type=="MultiLineString" else [g]):
            L=line.length;n=max(1,math.ceil(L/1000))
            for i in range(n):
                p=line.interpolate((i+.5)*L/n)
                if len(at.query(p,predicate="dwithin",distance=500)):continue
                good.append((p,L/n))
        total=sum(w for p,w in good)
        if total<5000:continue
        qa[kind+"_replay_eligible"]+=1
        screened=[]
        for p,w in good:
            if any(histories[int(i)][0]!=key for i in ht.query(p,predicate="dwithin",distance=500)):
                continue
            screened.append((p,w))
        km=sum(w for p,w in screened)/1000
        if km<5:
            qa[kind+"_cross_docket_geometry_excluded"]+=1
            continue
        pp=[];b=None
        for p,w in sorted(screened,key=lambda v:digest(key,round(v[0].x,3),round(v[0].y,3))):
            k=block(p)
            if b is not None and k!=b:continue
            if any(p.distance(q)<1000 for q in occupied):continue
            b=k
            pp.append(p);occupied.append(p)
            if len(pp)==3:break
        if pp:
            chosen.append({"kind":kind,"year":y,"km":km,"block":b,"epoch":ep,"points":pp})
        else:qa[kind+"_separation_excluded"]+=1
    if qa["abandoned_replay_eligible"]!=96 or qa["railbanked_replay_eligible"]!=50:
        raise ValueError("original 96/50 temporal census failed replay")
    return chosen,dict(qa)

def local_lines(trees,geo,p):
    region=box(p.x-1000,p.y-1000,p.x+1000,p.y+1000)
    lines=[]
    for tree,gg in zip(trees,geo):
        for i in tree.query(region,predicate="intersects"):
            g=gg[int(i)].intersection(region)
            if not g.is_empty:lines.append(("rail",transform(INV.transform,g)))
    return lines

def support(rr):
    return {k:{"n_dockets":sum(r["kind"]==k for r in rr),
               "n_blocks":len({r["block"] for r in rr if r["kind"]==k}),
               "n_sites":sum(len(r.get("points",[])) for r in rr if r["kind"]==k)} for k in KINDS}

def aggregate(r,vals):
    row={k:r[k] for k in ("kind","year","km","block")}
    for e in ("near","background"):
        row[e]={name:{k:float(np.mean([v[e][name][k] for v in vals])) for k in vals[0][e][name]} for name in NAMES}
    return row

def delta(r,name,metric):
    return r["near"][name][metric]-r["background"][name][metric]
def baseline_features(r):
    # Outcome-free pre-legal-completion covariates; all are fractions.
    return np.array([r[e][n]["pre1_cover"] for e in ("near","background") for n in NAMES]+
                    [delta(r,n,"pre5") for n in NAMES],float)
def match(rows,baseline=True):
    pairs=[]
    for b in sorted({r["block"] for r in rows}):
        A=[r for r in rows if r["kind"]=="abandoned" and r["block"]==b]
        B=[r for r in rows if r["kind"]=="railbanked" and r["block"]==b]
        if not A or not B:continue
        n,m=len(A),len(B)
        cost=np.full((n,m+n),1000.)
        cost[:,:m]=1e6
        for i,a in enumerate(A):
            for j,c in enumerate(B):
                dy=abs(a["year"]-c["year"])
                dl=abs(math.log(a["km"]/c["km"]))
                if dy>3 or dl>math.log(4):continue
                z=abs(baseline_features(a)-baseline_features(c)) if baseline else np.zeros(6)
                cap=np.array([.2]*4+[.1]*2)
                if (z>cap).any():continue
                cost[i,j]=dy/3+dl/math.log(4)+float((z/cap).mean())
        ii,jj=linear_sum_assignment(cost)
        for i,j in zip(ii,jj):
            if j<m and cost[i,j]<1000:pairs.append((A[i],B[j]))
    return pairs
def pair_result(pairs):
    blocks=collections.defaultdict(list)
    for a,b in pairs:blocks[a["block"]].append((a,b))
    n,G=len(pairs),len(blocks)
    qa={"n_pairs":n,"n_blocks":G}
    if n<MIN_PAIRS or G<MIN_PAIR_BLOCKS:return dict(qa,status="insufficient matching overlap")
    ans={}
    for name in NAMES:
        ans[name]={}
        for metric in ("span_minus_pre","pre5","span_per5","post_minus_pre"):
            xx=np.array([np.mean([delta(b,name,metric)-delta(a,name,metric) for a,b in vv]) for vv in blocks.values()])
            m=float(xx.mean());h=float(t.ppf(.975,G-1)*xx.std(ddof=1)/math.sqrt(G))
            ans[name][metric]={"estimate":m,"interval95":[m-h,m+h]}
    return dict(qa,status="estimated",contrast="railbanked minus abandoned",effects=ans,
                mean_abs_pre_covariate_difference=[float(x) for x in np.mean([abs(baseline_features(a)-baseline_features(b)) for a,b in pairs],axis=0)],
                mean_abs_year_difference=float(np.mean([abs(a["year"]-b["year"]) for a,b in pairs])))

def selftest():
    def row(k,y=2008,x=.2):
        return {"kind":k,"year":y,"km":10.,"block":(1,2),
                **{e:{n:{"pre1_cover":x,"pre5":0.,"span_minus_pre":.1,"span_per5":.1,"post_minus_pre":.1} for n in NAMES} for e in ("near","background")}}
    assert len(match([row("abandoned"),row("railbanked")]))==1
    assert not match([row("abandoned"),row("railbanked",2015)])
    assert not match([row("abandoned"),row("railbanked",x=.6)])
    assert pair_result(match([row("abandoned"),row("railbanked")]))["status"]=="insufficient matching overlap"
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    for name in ("rail-dir","census","out"):p.add_argument("--"+name,type=pathlib.Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    raw={};receipts={}
    for kind,(name,n) in s.LAYERS.items():
        raw[kind],receipts[kind]=s.fetch_layer(name,n)
        if receipts[kind]["receipt"]["sha256"]!=EXPECTED[kind]:raise ValueError("STB source digest changed")
    events,groupqa=s.group_events(raw)
    if groupqa.get("abandoned_complete_docket")!=423 or groupqa.get("railbanked_complete_docket")!=138:
        raise ValueError("previous STB docket grouping failed replay")
    paths=validate_files(a.rail_dir,json.loads(a.census.read_text()))
    active,inactive,mapqa=s.current_map(paths)
    prior=s.metric(events,active,inactive)
    for k,n,km in (("abandoned",96,2805.725),("railbanked",50,1661.985)):
        if prior[k]["n_temporal_ge5km"]!=n or prior[k]["screened_km_proxy"]!=km:
            raise ValueError("previous map screen failed replay")
    chosen,qa=select(events,raw,active,inactive)
    pre=support(chosen)
    # Geometry and date only overlap upper bound, no post outcome selection.
    upper=match(chosen,baseline=False)
    result={"version":"stb-glcfcs30d-temporal-pathways-v1","source":receipts,"map_replay":True,
        "map_qa":mapqa,"selection_qa":qa,"preoutcome_support":pre,
        "geometry_only_matching_upper_bound":{"pairs":len(upper),"blocks":len({a["block"] for a,b in upper})},
        "limits":["legal completion is not actual last train or trail opening; railbanking is not verified trail conversion",
                  "STB approximate applicant geometry; current 2026 status screens select on future geography",
                  "all other docket alignments and active/unknown rails excluded within 500m; up to three separated sites per group",
                  "same pixel masks before/after; local background and pretrend not external active-line controls",
                  "separate within-pathway change and pre-only matching are observational estimates, not causal effects",
                  "forest and impervious land-cover classes not actual plant height, biomass or species",
                  "GLC_FCS30D v2 CC BY 4.0 doi:10.5281/zenodo.15063683 and doi:10.5194/essd-16-1353-2024; no raw location/docket/pixels exported"]}
    print(json.dumps({"stb_preoutcome_support":pre,"geometry_matching":result["geometry_only_matching_upper_bound"],"qa":qa}),flush=True)
    if not all(v["n_dockets"]>=glc.MIN_GROUPS and v["n_blocks"]>=glc.MIN_BLOCKS for v in pre.values()):
        result["status"]="insufficient geometry support; no raster requests"
    else:
        try:
            idx,rc=glc.catalog(region=(-130,-65,20,55))
        except Exception as e:
            idx=None;result.update(status="blocked catalog",error_type=type(e).__name__,reason=str(e)[:300])
        if idx:
            result["catalog_receipt"]=rc
            histories=[g for rr in raw.values() for k,y,g in rr if g is not None]
            geo=[active,inactive,histories];trees=[STRtree(v) for v in geo]
            ts=glc.Tiles(idx);rows=[];status=collections.Counter()
            try:
                with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_HTTP_TIMEOUT="30",GDAL_HTTP_MAX_RETRY="0",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_CACHEMAX=128_000_000):
                    for r in chosen:
                        vals=[]
                        for point in r["points"]:
                            lon,lat=INV.transform(point.x,point.y)
                            try:
                                v,z=glc.windows(ts,r["epoch"],lon,lat,local_lines(trees,geo,point))
                            except Exception as e:
                                v,z=None,"source_or_geometry_error_"+type(e).__name__
                                if sum(n for k,n in status.items() if "error" in k)<2:
                                    print(json.dumps({"source_error_type":type(e).__name__,"reason":str(e)[:250]}),flush=True)
                            status[z]+=1
                            if v is not None:vals.append(v)
                        if len(vals)==len(r["points"]):rows.append(aggregate(r,vals))
                        if len(rows)>0 and len(rows)%20==0:
                            print(json.dumps({"completed_dockets":len(rows),"site_status":dict(status)}),flush=True)
                        if sum(n for k,n in status.items() if "error" in k)>=12 and status["complete"]==0:
                            result["early_stop"]="12 errors and zero complete sites";break
            finally:ts.close()
            result.update(status="analyzed",site_status=dict(status),tile_receipts=ts.receipts,
                complete_support=support(rows),within_pathways={k:glc.summarize([r for r in rows if r["kind"]==k]) for k in KINDS},
                primary_matched_comparison=pair_result(match(rows)))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps({k:v for k,v in result.items() if k!="tile_receipts"},ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()

"""Frozen 2021 land-cover chronosequence of historic passenger-rail terminal candidates."""
import argparse
import collections
import hashlib
import json
import math
import pathlib
import numpy as np
import rasterio
from scipy.stats import t
from shapely.geometry import Point,box
from shapely.ops import transform
from shapely.strtree import STRtree
from pyproj import CRS,Transformer
import audit_japan_terminal_rail_geometry as j
import audit_japan_terminal_overture as o
import analyze_world_railstatus_local_worldcover as wc
from build_global_Overture_physical_connector_components_v361 import validate_files

SEED="japan-terminal-worldcover-20261127-frozen-v1"
CR=CRS.from_proj4("+proj=aeqd +lat_0=36 +lon_0=137 +datum=WGS84 +units=m +no_defs")
INV=Transformer.from_crs(CR,4326,always_xy=True)
EQ=Transformer.from_crs(CR,6933,always_xy=True)
FWD=Transformer.from_crs(6668,CR,always_xy=True)
OUTCOMES=("tree","grass","built","shrub","crop","bare")
def grid(x,y,width):
    a,b=EQ.transform(x,y)
    return (math.floor(a/width),math.floor(b/width))
def score(p):
    return hashlib.sha256((SEED+"|"+str(p[0])+"|"+str(p[1])+"|"+str(round(p[2],2))+"|"+str(round(p[3],2))).encode()).digest()
def clean_points(points,inactive,other):
    ti,to=STRtree(inactive),STRtree(other)
    groups=collections.defaultdict(list)
    for p in points:
        key,yr,x,y,w=p
        q=Point(x,y)
        if len(to.query(q,predicate="dwithin",distance=500)):
            continue
        trace=bool(len(ti.query(q,predicate="dwithin",distance=80)))
        groups[(key,yr)].append((p,trace))
    assert abs(sum(p[4] for v in groups.values() for p,_ in v)/1000-3817.617)<.00051
    assert sum(sum(p[4] for p,_ in v)>=5000 for v in groups.values())==183
    assert sum(sum(p[4] for p,z in v if z)>=5000 for v in groups.values())==46
    return groups

def select(groups):
    # One group weight regardless of length. Geometry-only random ordering.
    # All points are in the same 200km block as the first selected point.
    # Cross-group separation is at least 1km in fixed projection; no replacements.
    chosen=[]; occupied=[]; qa=collections.Counter()
    order=sorted(groups,key=lambda z:hashlib.sha256((SEED+"|"+z[0]).encode()).digest())
    for key,yr in order:
        rows=groups[(key,yr)]
        if sum(p[4] for p,_ in rows)<5000:
            qa["below_5km"]+=1
            continue
        if yr>2015:
            qa["less_than_6_years_before_image"]+=1
            continue
        ss=[]; anchor=None
        for p,z in sorted(rows,key=lambda r:score(r[0])):
            x,y=p[2:4]
            b=grid(x,y,200000)
            if anchor is not None and b!=anchor:
                continue
            if any(math.hypot(x-a,y-bb)<1000 for a,bb in occupied):
                continue
            if anchor is None:
                anchor=b
            occupied.append((x,y))
            ss.append((p,z))
            if len(ss)==3:
                break
        if not ss:
            qa["no_spatially_separated_point"]+=1
            continue
        qa["selected_components"]+=1
        chosen.append((key,yr,ss))
    return chosen,dict(qa)

def local_lines(trees,geometries,x,y):
    # Retain the historical railway itself for the background exclusion, even
    # if erased from modern mapping. All passenger histories + all modern rails.
    region=box(x-800,y-800,x+800,y+800)
    lines=[]
    for tree,gg in zip(trees,geometries):
        for i in tree.query(region,predicate="intersects"):
            g=gg[int(i)].intersection(region)
            if not g.is_empty:
                lines.append(("any",transform(INV.transform,g)))
    return lines

def support(rows,width,min_n=50,min_blocks=15):
    by=collections.defaultdict(list)
    for r in rows:
        by[grid(r["x"],r["y"],width)].append(r)
    # Only within-block age variation is identified; no extrapolation across regions.
    eligible={b:v for b,v in by.items() if len(v)>=2 and max(r["age"] for r in v)-min(r["age"] for r in v)>=10}
    n=sum(map(len,eligible.values()))
    return eligible,{"n_groups":n,"n_blocks":len(eligible),"all_groups":len(rows),
        "all_blocks":len(by),"status":"estimated" if n>=min_n and len(eligible)>=min_blocks else "insufficient age overlap"}

def fit(rows,width=200000,endpoint="delta",name="tree",adjust_trace=False,min_n=50,min_blocks=15):
    blocks,qa=support(rows,width,min_n,min_blocks)
    if qa["status"]!="estimated":
        return qa
    X=[];Y=[];labels=[]
    for b,rr in blocks.items():
        x=np.array([[r["age"]/10,r["trace"]] if adjust_trace else [r["age"]/10] for r in rr],float)
        y=np.array([r[endpoint][name] for r in rr],float)
        X.extend(x-x.mean(axis=0));Y.extend(y-y.mean());labels.extend([b]*len(rr))
    X=np.asarray(X);Y=np.asarray(Y)
    n,p=X.shape;G=len(blocks)
    if np.linalg.matrix_rank(X)<p or n<=G+p:
        return dict(qa,status="rank deficient")
    inv=np.linalg.inv(X.T@X)
    beta=inv@X.T@Y
    e=Y-X@beta
    meat=np.zeros((p,p))
    for b in blocks:
        ix=np.array([c==b for c in labels]);s=X[ix].T@e[ix]
        meat+=np.outer(s,s)
    # CR1 including absorbed block intercepts, t critical df=G-1.
    cov=(G/(G-1))*((n-1)/(n-G-p))*inv@meat@inv
    se=math.sqrt(max(0.,float(cov[0,0])))
    crit=float(t.ppf(.975,G-1))
    return dict(qa,estimate=float(beta[0]),interval95=[float(beta[0]-crit*se),float(beta[0]+crit*se)],
        age_units="10 years",adjust_trace=adjust_trace)
def summarize(rows):
    out={"primary_tree_age_near_minus_background":fit(rows)}
    out["primary_near_background_components"]={e:fit(rows,endpoint=e) for e in ("near","background")}
    out["secondary_landcover_deltas"]={name:fit(rows,name=name) for name in OUTCOMES if name!="tree"}
    out["sensitivity"]={
        "100km":fit(rows,100000),"300km":fit(rows,300000),
        "adjust_current_map_trace_not_causal":fit(rows,adjust_trace=True),
        "no_trace_all_sampled_points":fit([r for r in rows if r["trace"]==0]),
        "trace_all_sampled_points":fit([r for r in rows if r["trace"]==1])}
    return out
def selftest():
    rr=[]
    # Two closure ages per spatial block; artificial exact +0.01 per decade.
    for b in range(4):
        for age in (10,30):
            x=(b*500000)
            rr.append({"x":x,"y":0,"age":age,"trace":0.,"delta":{"tree":age*.001+b*.03}})
    ans=fit(rr,min_n=8,min_blocks=3)
    assert ans["status"]=="estimated" and abs(ans["estimate"]-.01)<1e-9
    assert fit(rr)["status"]=="insufficient age overlap"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    for k in ("rail-dir","census","out"):
        p.add_argument("--"+k,type=pathlib.Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    rows,geoms,receipt,_=j.read_n05()
    cand,attrs=j.classification(rows)
    pp=[];old=j.audit_geometry(rows,geoms,cand,pp)
    if len(pp)!=4534 or old["all"]["screened_length_proxy_km"]!=4388.978 or len(cand)!=576:
        raise ValueError("prior N05 replay failed")
    paths=validate_files(a.rail_dir,json.loads(a.census.read_text()))
    inactive,other,qa=o.split_current(paths)
    replay=o.summarize(pp,inactive,other)
    if replay["all"]["no_competing_map_rail_km"]!=3817.617 or replay["all"]["components_ge_5km_without_competing_map_rail"]!=183:
        raise ValueError("prior mapped status replay failed")
    groups=clean_points(pp,inactive,other)
    chosen,sel=select(groups)
    potential=[{"x":ss[0][0][2],"y":ss[0][0][3],"age":2021-yr} for k,yr,ss in chosen]
    _,pre=support(potential,200000)
    print(json.dumps({"geometry_selection":sel,"preoutcome_age_support":pre,"points":sum(len(ss) for k,y,ss in chosen)}),flush=True)
    result={"version":"japan-terminal-worldcover-chronosequence-v1","source":receipt,"map_replay":True,
        "geometry_selection":sel,"preoutcome_age_support":pre,"n_selected_sites":sum(len(ss) for k,y,ss in chosen),
        "limits":["cross-sectional 2021 habitat, NOT a before-after or temporal recovery estimate",
        "candidate passenger terminal years, actual operations and freight/land conversion not independently verified",
        "2026 current map screening can select by survival or mapping; trace sensitivity not causal adjustment",
        "route history geometry is approximate, current near-background environment may differ from original habitat",
        "one group weight, fixed spatial blocks, conditional selection and support, no independent closure events guaranteed",
        "tree cover is not plant height, biomass, richness, or mowing; no iNat/eBird; no raw coordinates or group names exported"]}
    # No raster access at all if the frozen geometry overlap cannot identify slope.
    if pre["status"]!="estimated":
        result["status"]="blocked before raster extraction by frozen age-overlap gate"
    else:
        historical=[transform(FWD.transform,g) for g in geoms if not g.is_empty and g.is_valid]
        geom=[historical,inactive,other];trees=[STRtree(v) for v in geom]
        tiles=wc.Tiles();done=[];status=collections.Counter()
        try:
            with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_HTTP_TIMEOUT="45",GDAL_HTTP_MAX_RETRY="1",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
                for k,yr,ss in chosen:
                    vals=[]
                    for p,trace in ss:
                        x,y=p[2:4];xy=INV.transform(x,y)
                        try:
                            q=wc.site(tiles,xy,local_lines(trees,geom,x,y))
                        except Exception as e:
                            q={"status":"source_or_geometry_error_"+type(e).__name__}
                        status[q["status"]]+=1
                        if q["status"]=="complete":
                            vals.append((q,trace))
                    # No site replacement; require all preselected points complete.
                    if len(vals)==len(ss):
                        r={"x":ss[0][0][2],"y":ss[0][0][3],"age":2021-yr,"trace":float(np.mean([z for q,z in vals]))}
                        for endpoint in ("near","background"):
                            r[endpoint]={name:float(np.mean([q[endpoint]["fractions"][name] for q,z in vals])) for name in OUTCOMES}
                        r["delta"]={name:r["near"][name]-r["background"][name] for name in OUTCOMES}
                        done.append(r)
                    if sum(status.values())%25==0:
                        print(json.dumps({"sites_processed":sum(status.values()),"groups_complete":len(done)}),flush=True)
                    # Stop a clearly blocked source rather than cycling requests.
                    if sum(v for z,v in status.items() if z.startswith("source_or_geometry_error"))>=20 and status["complete"]==0:
                        result["early_stop"]="20 source/geometry errors and zero complete sites"
                        break
        finally:tiles.close()
        result.update(status="analyzed",site_status=dict(status),n_complete_groups=len(done),
            tile_receipts=tiles.receipts,models=summarize(done),
            age_bins_complete=[{"years":[lo,hi],"n":sum(lo<=2021-r["age"]<=hi for r in done)} for lo,hi in j.BINS])
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps({k:v for k,v in result.items() if k!="tile_receipts"},ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":
    main()

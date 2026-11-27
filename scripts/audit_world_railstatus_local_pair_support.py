#!/usr/bin/env python3
"""Global rail-length sample -> local mapped status pairs -> public API support.
Feasibility only, no plant-height outcomes, no observation/taxon/user IDs or GPS
written. No claims of active service, closed lines, or exact field vegetation.
"""
import argparse,hashlib,json,math,random,time,urllib.parse,urllib.request,urllib.error
from collections import Counter,defaultdict
from pathlib import Path
from shapely.geometry import Point,box
from shapely.ops import transform,nearest_points
from shapely.strtree import STRtree
from pyproj import Geod,Transformer
from build_global_Overture_physical_connector_components_v361 import validate_files,SOURCE_RELEASE
from sample_entire_global_GERS_railkm_WGS84_v368b import (
    stream_source_rails,seg_m,original_rail_flags_at_fraction,true_WGS84_geodesic_at)
from aggregate_world_rail_by_status_railstatus import category,flags_of,CATEGORIES,INACTIVE_STATUS
from scan_world_railstatus_cross_taxa import taxonomy,source_client
from scan_world_railstatus_disturbance_traits import taxonomy_metadata,load_life
from scan_world_railstatus_verified_plant_size import fetch_fig,parse_scaled,add_size,SIZE_SHA
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_coarse_cell
N_DRAWS=1024
MAX_API_PAIRS=128
SEED="tetsu-railstatus-local-pair-support-20261127-v1"
MIN_SEPARATION=600.
MAX_SEPARATION=5000.
EXCLUSION=350.
RADIUS=80.
MAX_UNCERT=30.
MAX_PAGES=3
PER_PAGE=200
BAD_SURFACE={"is_bridge","is_tunnel","is_covered"}
ACTIVE="active_normal";INACTIVE="inactive_normal"
ROOT=211194

def surface(row):
    return not bool(flags_of(row)&BAD_SURFACE)

def census(paths,status):
    geod=Geod(ellps="WGS84");q=Counter();total=Counter()
    for row in stream_source_rails(paths):
        q["rows"]+=1;c=category(row)
        if c is None:continue
        g,l=seg_m(row,geod);q[c]+=1;total[c]+=l
    if q["rows"]!=2054855:raise ValueError("Incomplete geometry census")
    for c in CATEGORIES:
        if q[c]!=status["segments"][c] or abs(total[c]/1000-status["total_km"][c])>.04:raise ValueError("Status length census changed "+c)
    return total,dict(q)

def sample(paths,total):
    geod=Geod(ellps="WGS84");rng=random.Random(SEED)
    targets=sorted((rng.random()*total,j) for j in range(N_DRAWS))
    out=[];k=0;acc=0.;seen=0
    for row in stream_source_rails(paths):
        seen+=1
        if category(row)!=INACTIVE:continue
        g,l=seg_m(row,geod)
        while k<len(targets) and targets[k][0]<acc+l:
            t,j=targets[k];f=(t-acc)/l
            x,y=true_WGS84_geodesic_at(g,f,geod)
            flags=original_rail_flags_at_fraction(row.get("rail_flags"),f)
            eligible=surface(row) and bool(flags&INACTIVE_STATUS)
            out.append({"draw":j,"lon":x,"lat":y,"eligible_surface":eligible})
            k+=1
        acc+=l
    if seen!=2054855 or k!=N_DRAWS or abs(acc-total)>.05:raise ValueError("Incomplete sampling replay")
    return sorted(out,key=lambda r:r["draw"])

def window(site):
    x,y=site["lon"],site["lat"]
    # At least 7km in each direction for |lat|<=75; exclude dateline
    # and high latitudes explicitly, without replacement.
    if not site["eligible_surface"] or abs(y)>75:return None
    dy=.07;dx=.07/math.cos(math.radians(abs(y)+dy))
    if x-dx<=-180 or x+dx>=180:return None
    return box(x-dx,y-dy,x+dx,y+dy)

def local_geometry(paths,sites):
    windows=[(s,window(s)) for s in sites];windows=[(s,g) for s,g in windows if g is not None]
    if not windows:return {},{"eligible_windows":0}
    tree=STRtree([g for s,g in windows]);local={s["draw"]:[] for s,g in windows};q=Counter()
    from shapely import wkb
    for row in stream_source_rails(paths):
        q["rows"]+=1
        g=wkb.loads(bytes(row["geometry"]))
        if g.geom_type!="LineString" or g.is_empty:raise ValueError("Invalid rail geometry")
        xmin,ymin,xmax,ymax=g.bounds
        if xmax-xmin>180:
            q["dateline_segments_skipped"]+=1
            continue
        for i in tree.query(g,predicate="intersects"):
            s,w=windows[int(i)]
            clipped=g.intersection(w)
            if clipped.is_empty:continue
            cat=category(row)
            # Do not treat bridge/tunnel/unknown/partly scoped inactive lines
            # as candidate sites. They still count as competing mapped rail.
            f=flags_of(row)
            whole_inactive=any(isinstance(r,dict) and set(r.get("values") or ())&INACTIVE_STATUS and
                (r.get("between") is None or list(r["between"])==[0,1]) for r in row.get("rail_flags") or ())
            ok=surface(row) and (cat==ACTIVE or (cat==INACTIVE and whole_inactive))
            local[s["draw"]].append((cat if ok else "other",clipped))
    if q["rows"]!=2054855:raise ValueError("Incomplete local geometry pass")
    return local,dict(q,eligible_windows=len(windows))

def pair(site,lines):
    if not lines:return None,"no_window_geometry"
    lon,lat=site["lon"],site["lat"]
    fwd=Transformer.from_crs(4326,f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs",always_xy=True)
    rev=Transformer.from_crs(fwd.target_crs,4326,always_xy=True)
    p=Point(0,0)
    xx=[(c,transform(fwd.transform,g)) for c,g in lines]
    def dist(x,cat):return min((x.distance(g) for c,g in xx if c==cat),default=math.inf)
    def competing(x,cat):return min((x.distance(g) for c,g in xx if c!=cat),default=math.inf)
    if dist(p,INACTIVE)>2:return None,"inactive_point_not_whole_surface"
    if competing(p,INACTIVE)<EXCLUSION:return None,"inactive_mixed_status"
    a=[(p.distance(g),g) for c,g in xx if c==ACTIVE]
    if not a:return None,"no_local_active"
    d,g=min(a,key=lambda z:z[0])
    if not MIN_SEPARATION<=d<=MAX_SEPARATION:return None,"separation_outside_gate"
    ap=nearest_points(p,g)[1]
    if competing(ap,ACTIVE)<EXCLUSION:return None,"active_mixed_status"
    # Nearest line point, not a station or measured trackside plot.
    x,y=rev.transform(ap.x,ap.y)
    if not all(map(math.isfinite,(x,y,d))):raise ValueError("Invalid local projection")
    return {"draw":site["draw"],"inactive":[lon,lat],"active":[x,y],"distance_m":d},"eligible"

def geometry(paths,status):
    lengths,cq=census(paths,status);sites=sample(paths,lengths[INACTIVE]);local,lq=local_geometry(paths,sites)
    qa=Counter();pairs=[]
    for s in sites:
        if window(s) is None:qa["no_surface_or_geographic_window"]+=1;continue
        p,reason=pair(s,local.get(s["draw"],[]));qa[reason]+=1
        if p is not None:pairs.append(p)
    # Avoid overlapping neighborhood queries: deterministic draw-order
    # selection, at most one pair per 100km equal-area cell. No biology used.
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    selected=[];seen=set();geod=Geod(ellps="WGS84")
    for p in pairs:
        k=source_coarse_cell(*p["inactive"],tf,100000)
        if k is None or k in seen:continue
        if any(min(abs(geod.inv(*a,*c)[2]) for a in (p["active"],p["inactive"]) for c in (q["active"],q["inactive"]))<1000 for q in selected):continue
        seen.add(k);selected.append(dict(p,block100=k))
    return selected[:MAX_API_PAIRS],{"length_km":{k:v/1000 for k,v in lengths.items()},"census":cq,"local_pass":lq,"draws":N_DRAWS,"pair_reasons":dict(qa),"n_spatially_separated_pairs":len(selected),"n_api_pairs":min(len(selected),MAX_API_PAIRS)}

class API:
    def __init__(self):self.n=0;self.last=0.
    def get(self,params):
        # Single client, at least 1.2 seconds between requests; bounded pages.
        # No media downloads, no authentication, no private coordinates.
        if self.n>=2*MAX_API_PAIRS*MAX_PAGES:raise RuntimeError("Fixed request budget exhausted")
        time.sleep(max(0.,1.2-(time.monotonic()-self.last)));self.last=time.monotonic();self.n+=1
        url="https://api.inaturalist.org/v1/observations?"+urllib.parse.urlencode(params)
        req=urllib.request.Request(url,headers={"User-Agent":"tetsu-railstatus-public-research-feasibility/1.0","Accept":"application/json"})
        with urllib.request.urlopen(req,timeout=45) as r:
            data=r.read(16000001)
        if len(data)>16000000:raise ValueError("Response cap")
        return json.loads(data)
    def records(self,site):
        x,y=site;params={"lat":y,"lng":x,"radius":.15,"taxon_id":ROOT,"quality_grade":"research",
          "d1":"2018-01-01","d2":"2025-12-31","per_page":PER_PAGE,"order_by":"id","order":"asc"}
        rows=[];total=None;ids=set()
        for p in range(1,MAX_PAGES+1):
            z=self.get(dict(params,page=p));n=z.get("total_results");r=z.get("results")
            if not isinstance(n,int) or n<0 or not isinstance(r,list) or len(r)>PER_PAGE:raise ValueError("Unexpected API response")
            if total is None:total=n
            if total!=n:raise ValueError("Changing query population")
            if total>MAX_PAGES*PER_PAGE:return None,"page_cap"
            for v in r:
                if not isinstance(v.get("id"),int) or v["id"] in ids:raise ValueError("Missing or repeated observation ID")
                ids.add(v["id"]);rows.append(v)
            if len(rows)==total:return rows,"complete"
            if len(rows)>total or not r:raise ValueError("Incomplete query")
        return None,"incomplete"

def accept(r,site,canonical,allowed):
    # Missing privacy flag is not confirmation. Do not inspect private_geojson
    # or use obscured coordinates. Require explicit unobscured public record.
    if r.get("obscured") is not False:return None,"obscured_or_unknown"
    if any(k not in r or r[k] not in (None,"open") for k in ("geoprivacy","taxon_geoprivacy")):return None,"privacy_status_not_confirmed"
    if r.get("quality_grade")!="research" or not r.get("photos"):return None,"quality_or_no_photo"
    date=r.get("observed_on")
    if not isinstance(date,str) or not "2018-01-01"<=date<="2025-12-31":return None,"date"
    u=r.get("positional_accuracy")
    if not isinstance(u,(float,int)) or isinstance(u,bool) or not math.isfinite(u) or not 0<=u<=MAX_UNCERT:return None,"accuracy"
    p=r.get("geojson")
    if not isinstance(p,dict) or p.get("type")!="Point":return None,"public_geometry"
    xy=p.get("coordinates")
    if not isinstance(xy,list) or len(xy)!=2 or not all(isinstance(x,(float,int)) and math.isfinite(x) for x in xy):return None,"public_geometry"
    if not -180<=xy[0]<=180 or not -90<=xy[1]<=90:return None,"public_geometry"
    d=abs(Geod(ellps="WGS84").inv(*site,*xy)[2])
    if d+u>RADIUS:return None,"outside_conservative_disk"
    t=r.get("taxon") or {};hit=canonical.get(t.get("id"))
    if hit is None or hit[0]!="plants" or hit[1] not in allowed:return None,"frozen_taxonomy_unmatched"
    o=(r.get("user") or {}).get("id")
    if not isinstance(o,int):return None,"observer_missing"
    return (hit[1],o),"eligible"

def support(api,site,canonical,allowed,life):
    try:rows,reason=api.records(site)
    except urllib.error.HTTPError as e:
        if e.code in (403,429):raise RuntimeError("API access/rate limit; stop without retry") from e
        return {"status":"http_"+str(e.code)}
    except (ValueError,TimeoutError,urllib.error.URLError):return {"status":"query_failed"}
    if rows is None:return {"status":reason}
    qa=Counter();species=set();observers=set()
    for r in rows:
        v,k=accept(r,site,canonical,allowed);qa[k]+=1
        if v is not None:species.add(v[0]);observers.add(v[1])
    known={s for s in species if life.get(s,{}).get("height_nonwoody") is not None}
    nw={s for s in species if life.get(s,{}).get("woody")==0.}
    return {"status":"complete","query_total":len(rows),"qa":dict(qa),"species":len(species),"observers":len(observers),
        "known_nonwoody":len(nw),"height_known":len(known),"height_coverage":len(known)/len(nw) if nw else None,
        "gate":qa["eligible"]>=20 and len(observers)>=5 and len(known)>=10 and len(known)>=.5*len(nw)}

def main():
    p=argparse.ArgumentParser()
    for k in ("rail-dir","census","status","taxa","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();status=json.loads(a.status.read_text());prior=json.loads(a.taxa.read_text())
    if status["source_Overture_release"]!=SOURCE_RELEASE or not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong map")
    paths=validate_files(a.rail_dir,json.loads(a.census.read_text()))
    pairs,qa=geometry(paths,status)
    print(json.dumps({"local_rail_geometry_support":qa},allow_nan=False),flush=True)
    client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=prior["taxa_sha256"]:raise ValueError("Frozen taxonomy changed")
    names,allowed,fam,roots=taxonomy_metadata(client,canonical,sha)
    if roots["tracheophyta"]!=ROOT:raise ValueError("Target root changed")
    life=load_life(names["plants"],allowed["plants"])
    rows,trait_audit=parse_scaled(fetch_fig(25600650,45643398,SIZE_SHA))
    trait_audit.update(add_size(life,rows,names["plants"],allowed["plants"]))
    if trait_audit["known_nonwoody"]["height_nonwoody"]!=5506:raise ValueError("Prior height coverage not replayed")
    del rows,names
    api=API();results=[];blocked=None
    for q in pairs:
        try:
            a0=support(api,q["active"],canonical,allowed["plants"],life)
            i0=support(api,q["inactive"],canonical,allowed["plants"],life)
        except RuntimeError as e:blocked=str(e);break
        results.append({"draw":q["draw"],"block100":q["block100"],"distance_m":round(q["distance_m"],1),"active":a0,"inactive":i0})
        if len(results)%16==0:print(json.dumps({"local_pair_api_progress":{"pairs":len(results),"requests":api.n,"both_gate":sum(r["active"].get("gate",False) and r["inactive"].get("gate",False) for r in results)}}),flush=True)
    done=[r for r in results if all(r[k]["status"]=="complete" for k in ("active","inactive"))]
    supported=[r for r in done if all(r[k]["gate"] for k in ("active","inactive"))]
    out={"version":"railstatus-global-local-pair-support-v1","geometry":qa,"trait_source":trait_audit,"source_taxonomy_sha256":sha,
         "api":{"requests":api.n,"blocked":blocked,"n_complete_pairs":len(done),"n_both_height_supported":len(supported),"n_attempted_pairs":len(results),
                "site_status":{k:dict(Counter(r[k]["status"] for r in results)) for k in ("active","inactive")}},
         "pairs":results,"limits":["feasibility only; no trait means or effects computed","live public API differs from frozen bulk snapshot, historical observations under current map classification",
          "80m conservative disks centered on sampled rail point, not whole corridor or matched habitat","same neighborhood does not equal same vegetation or sampling effort",
          "mapped inactive whole segment, not actual closure date or management","unobscured public coordinates only, no raw observation data exported","fixed selection cap and page cap; incomplete queries are not absence","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({"local_pair_api_support":out["api"]}),flush=True)
if __name__=="__main__":main()

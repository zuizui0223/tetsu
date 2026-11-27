#!/usr/bin/env python3
"""Historical map-state corroboration on the exact frozen WorldCover sample.

Query historical public OSM geometries before looking at satellite outcomes.
No claim of actual operation or closure dates. No individual mapper data or
source geometries exported. Failed/missing historical mapping is not absence.
"""
import argparse,hashlib,json,math,time,urllib.parse,urllib.request,urllib.error
from collections import Counter
from datetime import datetime
from pathlib import Path
from shapely.geometry import Point,LineString
from shapely.ops import transform,unary_union
from pyproj import Transformer
import analyze_world_railstatus_local_worldcover as wc
import audit_world_railstatus_local_pair_support as a

URL="https://overpass-api.de/api/interpreter"
DATES=("2021-01-01T00:00:00Z","2021-12-31T23:59:59Z")
BATCH=8
QUERY_RADIUS=400
CONFIRM_DISTANCE=30.
MIN_NEAR_LENGTH=40.
NEAR=80.
EXCLUSION=350.
RAIL_TYPES={"rail","narrow_gauge","disused","abandoned","razed","dismantled","construction","light_rail","tram","subway","monorail","funicular","preserved","miniature"}
NORMAL={"rail","narrow_gauge"}
INACTIVE={"abandoned","disused"}

def dateval(s):
    return datetime.fromisoformat(s.replace("Z","+00:00"))

def classify(tags):
    value=tags.get("railway")
    lifecycle={k.split(":")[0] for k in ("abandoned:railway","disused:railway") if tags.get(k) in NORMAL}
    if any(tags.get(k) not in (None,"no","0","false") for k in ("bridge","tunnel","covered")):return "other",None
    if tags.get("railway:preserved")=="yes":return "other",None
    if value in NORMAL and not lifecycle and tags.get("abandoned")!="yes" and tags.get("disused")!="yes":
        return a.ACTIVE,None
    if value in INACTIVE:
        other=tags.get(value+":railway")
        if other is not None and other not in NORMAL:return "other",None
        if lifecycle-{value}:return "other",None
        return a.INACTIVE,value
    if value is None and len(lifecycle)==1:return a.INACTIVE,next(iter(lifecycle))
    return "other",None

def query_string(sites,date):
    if date not in DATES:raise ValueError("Unfrozen date")
    q=[]
    for lon,lat in sites:
        if not all(map(math.isfinite,(lon,lat))) or not -180<=lon<=180 or not -75<=lat<=75:raise ValueError("Invalid public railway point")
        around=f"(around:{QUERY_RADIUS},{lat:.7f},{lon:.7f})"
        q.append('way["railway"~"^('+"|".join(sorted(RAIL_TYPES))+')$"]'+around+";")
        for k in ("abandoned:railway","disused:railway"):
            q.append('way["'+k+'"~"^(rail|narrow_gauge)$"]'+around+";")
    return '[out:json][timeout:60][maxsize:134217728][date:"'+date+'"];('+''.join(q)+');out meta geom;'

class HistoryAPI:
    def __init__(self):
        self.n=0;self.last=0.
        self.receipts=[]
    def get(self,sites,date):
        # A small number of batched local geometry requests, single client,
        # no automatic retry or rotation across public servers.
        if self.n>=2+len(DATES)*math.ceil(wc.CAP/BATCH):raise RuntimeError("Fixed public API budget exhausted")
        q=query_string(sites,date)
        time.sleep(max(0.,8.-(time.monotonic()-self.last)));self.last=time.monotonic();self.n+=1
        data=urllib.parse.urlencode({"data":q}).encode()
        req=urllib.request.Request(URL,data=data,headers={"User-Agent":"tetsu-railstatus-historical-map-academic/1.0","Accept":"application/json","Content-Type":"application/x-www-form-urlencoded"})
        try:
            with urllib.request.urlopen(req,timeout=90) as r:raw=r.read(16000001)
        except urllib.error.HTTPError as e:raise RuntimeError("Historical API HTTP "+str(e.code)) from e
        except (urllib.error.URLError,TimeoutError) as e:raise RuntimeError("Historical API unavailable") from e
        if len(raw)>16000000:raise RuntimeError("Historical API response cap")
        z=json.loads(raw)
        if not isinstance(z.get("elements"),list) or z.get("remark") or not z.get("osm3s",{}).get("timestamp_osm_base"):raise RuntimeError("Incomplete historical API response")
        rows=[];seen=set()
        for x in z["elements"]:
            if x.get("type")!="way" or not isinstance(x.get("id"),int) or x["id"] in seen or not isinstance(x.get("tags"),dict):raise ValueError("Unexpected OSM way")
            seen.add(x["id"])
            if not isinstance(x.get("timestamp"),str) or dateval(x["timestamp"])>dateval(date):raise ValueError("Historical way timestamp failed")
            v=x.get("geometry")
            if not isinstance(v,list) or len(v)<2 or any(p is None or "lat" not in p or "lon" not in p for p in v):raise ValueError("Incomplete historical way geometry")
            xy=[(float(p["lon"]),float(p["lat"])) for p in v]
            if any(not -180<=lo<=180 or not -90<=la<=90 for lo,la in xy) or max(lo for lo,la in xy)-min(lo for lo,la in xy)>180:raise ValueError("Ambiguous historical way geometry")
            g=LineString(xy)
            if not g.is_valid or g.is_empty:raise ValueError("Invalid historical way")
            c,sub=classify(x["tags"]);rows.append((c,sub,g))
        self.receipts.append({"date":date,"n_sites":len(sites),"n_ways":len(rows),"response_sha256":hashlib.sha256(raw).hexdigest(),
                              "bytes":len(raw),"query_sha256":hashlib.sha256(q.encode()).hexdigest(),
                              "osm_base":z["osm3s"]["timestamp_osm_base"]})
        return rows

def site_state(rows,xy,expected):
    lon,lat=xy
    tf=Transformer.from_crs(4326,f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs",always_xy=True)
    p=Point(0,0);disk=p.buffer(NEAR)
    matching=[];competing=[];subs=set()
    # Batch responses include geographically distant ways; bbox filter
    # avoids projecting unrelated objects near the antipode.
    dy=.005;dx=dy/math.cos(math.radians(abs(lat)+dy))
    from shapely.geometry import box
    local=box(lon-dx,lat-dy,lon+dx,lat+dy)
    for c,sub,g in rows:
        if not g.intersects(local):continue
        gg=transform(tf.transform,g.intersection(local))
        if gg.is_empty:continue
        if c==expected:
            matching.append(gg)
            if gg.distance(p)<=CONFIRM_DISTANCE and sub:subs.add(sub)
        else:competing.append(gg)
    md=min((g.distance(p) for g in matching),default=math.inf)
    cd=min((g.distance(p) for g in competing),default=math.inf)
    L=float(unary_union(matching).intersection(disk).length) if matching else 0.
    status="confirmed" if md<=CONFIRM_DISTANCE and L>=MIN_NEAR_LENGTH and cd>=EXCLUSION else (
        "competing_or_mixed" if cd<EXCLUSION else "expected_line_not_confirmed")
    return {"status":status,"expected_distance_m":round(md,2) if math.isfinite(md) else None,
            "competing_distance_m":round(cd,2) if math.isfinite(cd) else None,
            "expected_near_length_m":round(L,2),"inactive_subtypes":sorted(subs)}

def summarize(records):
    done=[r for r in records if all(y in r["dates"] for y in DATES)]
    stable=[r for r in done if all(r["dates"][y][k]["status"]=="confirmed" for y in DATES for k in ("active","inactive"))]
    sub=Counter()
    for r in stable:
        z=[r["dates"][y]["inactive"]["inactive_subtypes"] for y in DATES]
        label=z[0][0] if z[0]==z[1] and len(z[0])==1 else "mixed_or_changed"
        sub[label]+=1
    return {"complete_two_date_pairs":len(done),"both_status_confirmed_at_both_dates":len(stable),"stable_inactive_subtypes":dict(sub),
            "status_by_date":{y:{k:dict(Counter(r["dates"][y][k]["status"] for r in records if y in r["dates"])) for k in ("active","inactive")} for y in DATES}}

def main():
    p=argparse.ArgumentParser()
    for k in ("rail-dir","census","status","out"):p.add_argument("--"+k,type=Path,required=True)
    args=p.parse_args();api=HistoryAPI()
    out={"version":"railstatus-worldcover-historical-map-v1","source_url":URL,"dates":DATES,"source":"OpenStreetMap contributors, ODbL; historical Overpass extraction",
         "geometry":None,"pairs":[],"blocked":None,"preflight":None}
    try:
        # Fixed public railway example independent of the sampled outcomes.
        pre=api.get([(8.675,49.403)],DATES[0])
        if not any(c==a.ACTIVE for c,_,g in pre):raise RuntimeError("Historical rail preflight did not confirm source content")
        out["preflight"]={"status":"passed","n_ways":len(pre)}
        status=json.loads(args.status.read_text())
        if status["source_Overture_release"]!=a.SOURCE_RELEASE or not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong rail source")
        paths=a.validate_files(args.rail_dir,json.loads(args.census.read_text()))
        pairs,local,geo=wc.frozen_points(paths,status);del local
        out["geometry"]=geo
        print(json.dumps({"historical_map_geometry":geo},allow_nan=False),flush=True)
        # Identical geometry-only draw order and cap as the satellite run.
        for i in range(0,len(pairs),BATCH):
            batch=pairs[i:i+BATCH]
            rows=[{"draw":r["draw"],"block100":r["block100"],"separation_m":round(r["distance_m"],1),"dates":{}} for r in batch]
            out["pairs"].extend(rows)
            sites=[r[k] for r in batch for k in ("active","inactive")]
            for y in DATES:
                got=api.get(sites,y)
                for q,r in zip(batch,rows):
                    r["dates"][y]={k:site_state(got,q[k],c) for k,c in (("active",a.ACTIVE),("inactive",a.INACTIVE))}
            print(json.dumps({"historical_map_progress":summarize(out["pairs"]),"requests":api.n},allow_nan=False),flush=True)
    except RuntimeError as e:out["blocked"]=str(e)
    out.update(receipts=api.receipts,requests=api.n,support=summarize(out["pairs"]),
               limits=["map snapshots at start/end of 2021, not actual historical train service or maintenance, no proof of continuous status",
                       "mapped disused (tracks remain) and abandoned (removed) distinct historical tag proxies",
                       "unconfirmed, failed or incomplete queries not physical absence",
                       "historical OSM and Overture are not independent mapping systems",
                       "same geometry-only satellite sample; no outcome read, no raster effect estimates",
                       "no raw OSM IDs, user data or source GPS exported"])
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({"historical_map_result":{"preflight":out["preflight"],"blocked":out["blocked"],"requests":api.n,"support":out["support"]}},allow_nan=False),flush=True)
if __name__=="__main__":main()

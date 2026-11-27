"""Academic-use only. Source-only N05 terminal-route geometry audit; no ecology."""
import argparse
import collections
import hashlib
import io
import json
import math
import pathlib
import re
import zipfile

import shapefile
from pyproj import CRS, Transformer
from shapely.geometry import shape
from shapely.ops import transform, unary_union, linemerge
from shapely.strtree import STRtree

from audit_dated_rail_closure_sources import get, year

URL = "https://nlftp.mlit.go.jp/ksj/gml/data/N05/N05-25/N05-25_GML.zip"
SHA = "b8971863850d6d52b8916a6aec1de29e88f99063a4a8c3aacdc65724d3c41a3c"
RID = re.compile(r"EB02_\d{5}", re.I)
BINS = [(1950,1984),(1985,2000),(2001,2010),(2011,2015),(2016,2020),(2021,2025)]
def val(x):
    return str(x if x is not None else "").strip()
def decode_links(x):
    s = val(x).upper()
    if not s:
        return [], False
    matches = RID.findall(s)
    rem = RID.sub("",s)
    return matches, bool(re.sub(r"[\s,;、，；/|・]+","",rem)) or not matches

def classification(rows):
    groups = collections.defaultdict(list)
    bad = 0
    for i,r in enumerate(rows):
        rid = val(r["N05_006"]).upper()
        if not RID.fullmatch(rid):
            bad += 1
            continue
        groups[rid].append(i)
    parent = {k:k for k in groups}
    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k
    def union(a,b):
        a,b = find(a),find(b)
        if a != b:
            parent[max(a,b)] = min(a,b)
    unresolved = set()
    linked = 0
    for k,inds in groups.items():
        for i in inds:
            links,bad_link = decode_links(rows[i]["N05_010"])
            if bad_link:
                unresolved.add(k)
            for link in links:
                if link in groups:
                    union(k,link)
                    linked += 1
                else:
                    unresolved.add(k)
    components = collections.defaultdict(list)
    for k,inds in groups.items():
        components[find(k)].extend(inds)
    uncertain = {find(k) for k in unresolved}
    counts = collections.Counter()
    candidates = []
    for k,inds in sorted(components.items()):
        ends = [val(rows[i]["N05_005e"]) for i in inds]
        ys = [year(e) for e in ends]
        if "9999" in ends:
            counts["includes_current_record"] += 1
        elif k in uncertain:
            counts["unresolved_relationship"] += 1
        elif any(y is None or y < 1949 or y > 2024 for y in ys):
            counts["unknown_or_out_of_window_end"] += 1
        else:
            latest = max(ys)
            final = [i for i,y in zip(inds,ys) if y == latest]
            starts = [year(rows[i]["N05_005b"]) for i in final]
            if any(y is None or y < 1950 or y > latest for y in starts):
                counts["terminal_start_invalid"] += 1
            elif any(val(rows[i].get("N05_009")) for i in final):
                counts["terminal_geometry_note_requires_review"] += 1
            else:
                counts["terminal_component_candidate"] += 1
                candidates.append((k,latest+1,final))
    return candidates, {"rows":len(rows),"rows_invalid_primary_id":bad,
        "primary_groups":len(groups),"components":len(components),
        "other_id_link_references":linked,"groups_with_unresolved_link":len(unresolved),
        "component_classes":dict(counts)}

def parts(g):
    if g.geom_type == "LineString":
        if g.length > 0:
            yield g
    elif g.geom_type in ("MultiLineString","GeometryCollection"):
        for p in g.geoms:
            yield from parts(p)

def audit_geometry(rows, geoms, candidates, point_sink=None):
    # Japan-centred AEQD approximation, metre units. The 500m gate is intentionally
    # wider than the 300m future background buffer; final exact local audit required.
    crs=CRS.from_proj4("+proj=aeqd +lat_0=36 +lon_0=137 +datum=WGS84 +units=m +no_defs")
    project=Transformer.from_crs("EPSG:6668",crs,always_xy=True).transform
    projected=[]
    for g in geoms:
        if g.is_empty or g.geom_type not in ("LineString","MultiLineString") or not g.is_valid:
            projected.append(None)
        else:
            projected.append(transform(project,g))
    # Current or unknown status as possible competitor; no assumption of freight absence.
    competing=[g for r,g in zip(rows,projected) if g is not None and
        (val(r["N05_005e"])=="9999" or year(r["N05_005e"]) is None)]
    if not competing:
        raise ValueError("no current/unknown passenger geometry")
    tree=STRtree(competing)
    # Also remove alignments that appear in a later historical record, including
    # discontinued successor routes under different relationship IDs.
    later_inds=[i for i,g in enumerate(projected) if g is not None]
    hist_tree=STRtree([projected[i] for i in later_inds])
    summaries=[]
    for key,y,inds in candidates:
        if any(projected[i] is None for i in inds):
            continue
        u=unary_union([projected[i] for i in inds])
        lines=list(parts(u))
        if len(lines)>1:
            lines=list(parts(linemerge(lines)))
        total=clean=0.
        n=0
        for line in lines:
            L=line.length
            if not math.isfinite(L):
                raise ValueError("invalid metric length")
            total+=L
            # Exact equal-weight midpoints within each merged line, max step 1km.
            m=max(1,math.ceil(L/1000))
            for j in range(m):
                p=line.interpolate((j+.5)*L/m)
                if tree.query(p,predicate="dwithin",distance=500).size:
                    continue
                nearby=hist_tree.query(p,predicate="dwithin",distance=500)
                later=False
                for q in nearby:
                    i=later_inds[int(q)]
                    start=year(rows[i]["N05_005b"])
                    end=year(rows[i]["N05_005e"])
                    # A continuation in the same or another ID, possibly later closed.
                    if end is not None and end>=y and (start is None or start<=end):
                        later=True
                        break
                if not later:
                    clean+=L/m
                    n+=1
                    if point_sink is not None:
                        point_sink.append((key,y,p.x,p.y,L/m))
        summaries.append({"year":y,"total_km":total/1000,"screened_km":clean/1000,"screened_points":n})
    def sum_group(s):
        return {"candidate_components":len(s),"terminal_geometry_km":round(sum(r["total_km"] for r in s),3),
            "screened_length_proxy_km":round(sum(r["screened_km"] for r in s),3),
            "screened_points":sum(r["screened_points"] for r in s),
            "components_ge_5_screened_km":sum(r["screened_km"]>=5 for r in s)}
    return {"valid_current_or_unknown_geometries":len(competing),"components_with_geometry":len(summaries),
        "all":sum_group(summaries),"year_bins":[{"years":[a,b],**sum_group([r for r in summaries if a<=r["year"]<=b])} for a,b in BINS],
        "post_2000_candidates_ge_5km":sum(r["year"]>=2001 and r["screened_km"]>=5 for r in summaries),
        "post_2015_candidates_ge_5km":sum(r["year"]>=2016 and r["screened_km"]>=5 for r in summaries)}

def selftest():
    assert decode_links("EB02_00002、EB02_00003")==(["EB02_00002","EB02_00003"],False)
    assert decode_links("nonsense")[1]
    def row(k,end,link=""):
        return {"N05_006":k,"N05_005e":end,"N05_005b":"1950","N05_010":link,"N05_009":""}
    rows=[row("EB02_00001","2009","EB02_00002"),row("EB02_00002","9999"),row("EB02_00003","2019"),row("EB02_00004","2011","EB02_99999")]
    c,a=classification(rows)
    assert len(c)==1 and c[0][1]==2020 and a["component_classes"]["includes_current_record"]==1
    assert a["component_classes"]["unresolved_relationship"]==1

def read_n05():
    raw,receipt=get(URL,70_000_000)
    if receipt["sha256"]!=SHA:
        raise ValueError("pinned source digest mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=[n for n in z.namelist() if pathlib.PurePosixPath(n).name=="N05-25_RailroadSection2.dbf" and "UTF-8" in pathlib.PurePosixPath(n).parts]
        if len(members)!=1:
            raise ValueError("unique UTF-8 route DBF not found")
        prefix=members[0][:-4]
        src={e:z.read(prefix+"."+e) for e in ("dbf","shp","shx","prj")}
        reader=shapefile.Reader(dbf=io.BytesIO(src["dbf"]),shp=io.BytesIO(src["shp"]),shx=io.BytesIO(src["shx"]),encoding="utf-8")
        rows=[];geoms=[]
        for record in reader.iterShapeRecords():
            rows.append(record.record.as_dict())
            geoms.append(shape(record.shape.__geo_interface__))
    return rows,geoms,receipt,src

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--out",default="results/railstatus/japan_terminal_geometry_audit.json")
    args=p.parse_args()
    selftest()
    if args.self_test:
        return
    rows,geoms,receipt,src=read_n05()
    candidates,attrs=classification(rows)
    result={"version":"japan-n05-terminal-route-audit-v1","archive":receipt,
        "source_member_sha256":{e:hashlib.sha256(b).hexdigest() for e,b in src.items()},
        "row_end_class":dict(collections.Counter("current" if val(r["N05_005e"])=="9999" else "finite" if year(r["N05_005e"]) is not None else "unknown" for r in rows)),
        "attributes":attrs,"geometry":audit_geometry(rows,geoms,candidates),
        "limits":["terminal passenger relationship components, not verified independent closures","unknown and current passenger rail plus later historic geometry within 500m excluded","screened km approximate midpoint-weighted length, not observed abandoned land","freight, trail conversion, actual operations, map errors unverified","no raw coordinates, route IDs, route names or source files exported","no biological or satellite outcomes used"]}
    out=pathlib.Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False))

if __name__=="__main__":
    main()

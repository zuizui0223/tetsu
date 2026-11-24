#!/usr/bin/env python3
"""v378a: GloNAF 2.02 original source support on fixed all-world 2048 rail draws."""
import argparse, csv, hashlib, json, re
from pathlib import Path
from collections import Counter, defaultdict
from shapely.geometry import Point, Polygon
from pyproj import Geod
import geopandas as gpd
from analyze_true_world_railkm_pps_point_climate_v374 import valid_source_sample, land_source_continent
from generate_global_equal_area_climate_land_grid_v351 import source_continents

MD5={
 "glonaf_2024_nol.zip":"e081bda4b0768136a1f88f6133d26637",
 "glonaf_flora2.csv":"8e1a5d4e2c08a815cc415a0a9b57e31d",
 "glonaf_list.csv":"efbe00842ea4f6425721106e5a1d6dfb",
 "glonaf_region.csv":"6bf6da4a226cacddcd3c037299348d69",
 "glonaf_taxon_wcvp.csv":"02b6bb3d34305ed38a47c794c542159b",
}
def verified(root):
    for name,h in MD5.items():
        with open(root/name,"rb") as f:
            if hashlib.file_digest(f,"md5").hexdigest()!=h:
                raise ValueError("GloNAF original MD5 mismatch: "+name)

def load(path, required):
    with path.open(encoding="utf-8-sig",newline="") as f:
        d=csv.DictReader(f)
        if not required<=set(d.fieldnames or []):
            raise ValueError("Original source fields missing: "+str(path))
        return list(d)

def norm(value):
    return str(value or "").strip()
def valid_name(s):
    return re.fullmatch(r"[A-Z][a-z]+ [a-z][a-z-]+",s) is not None

def complete_vascular(row):
    return (norm(row["taxon_group_id"])=="1" and norm(row["completeness_id"])=="3"
        and norm(row["invasives_only"]).lower() in ("0","false","f"))

def unique_map(rows,key):
    result={}
    for r in rows:
        k=norm(r[key])
        if not k or k in result:raise ValueError("Missing or duplicate source ID: "+key)
        result[k]=r
    return result

def region_key(s, by_code, by_obj):
    vals={x.lower():norm(v) for x,v in s.items()}
    a=by_code.get(vals.get("code",""))
    b=by_obj.get(vals.get("objidsic",""))
    if a and b and a!=b:raise ValueError("Source region geometry identifiers disagree")
    if not (a or b):raise ValueError("Unresolved original GloNAF region polygon code/OBJIDsic; inspect schema")
    return a or b

def classify_window(geo, lon, lat, km):
    # Conservative local polygon approx; near antimeridian do not infer.
    coords=[geo.fwd(lon,lat,a,km*1000)[:2] for a in range(0,360,5)]
    if any(abs(x-lon)>180 for x,y in coords):return None
    return Polygon(coords)

def analyze(root, sample, ne):
    verified(root)
    regions=load(root/"glonaf_region.csv",{"id","code"})
    # Publisher's actual case for OBJIDsic is audited, not guessed.
    reg=unique_map(regions,"id")
    code=unique_map(regions,"code")
    by_code={k:norm(v["id"]) for k,v in code.items()}
    objfield=next((k for k in regions[0] if k.lower()=="objidsic"),None)
    by_obj={}
    if objfield:
        for row in regions:
            v=norm(row[objfield])
            if v and v.lower()!="na":
                if v in by_obj:raise ValueError("Ambiguous source OBJIDsic")
                by_obj[v]=norm(row["id"])
    lists=load(root/"glonaf_list.csv",{"id","region_id","taxon_group_id","completeness_id","invasives_only"})
    li=unique_map(lists,"id")
    taxa=unique_map(load(root/"glonaf_taxon_wcvp.csv",{"id","accepted_plant_name_id","taxa_accepted","wcvp_taxon_rank"}),"id")
    flora=load(root/"glonaf_flora2.csv",{"taxon_wcvp_id","region_id","list_id"})
    if (len(reg),len(lists),len(taxa),len(flora))!=(1343,2274,95635,349271):
        raise ValueError("Source row counts disagree with original 2.02 documentation")
    complete={norm(x["region_id"]) for x in lists if complete_vascular(x)}
    if not complete<=reg.keys():raise ValueError("List region not in original regions")
    name_ids=defaultdict(set); region_names=defaultdict(set); q=Counter()
    for r in flora:
        rid=norm(r["region_id"]); lid=norm(r["list_id"]); tid=norm(r["taxon_wcvp_id"])
        if rid not in reg or lid not in li or tid not in taxa:
            raise ValueError("Original flora join key missing")
        if norm(li[lid]["region_id"])!=rid:
            raise ValueError("Original flora/list region disagree")
        t=taxa[tid]; n=norm(t["taxa_accepted"]); aid=norm(t["accepted_plant_name_id"])
        if (not valid_name(n) or norm(t["wcvp_taxon_rank"]).lower()!="species"
              or aid.lower() in ("","na","nan")):
            q["flora_rows_excluded_non_exact_species_or_unaccepted"]+=1
            continue
        name_ids[n].add(aid); region_names[rid].add(n)
    unambig={n for n,ids in name_ids.items() if len(ids)==1}
    q["n_documented_exact_unambiguous_species_names"]=len(unambig)
    q["ambiguous_species_names_excluded"]=len(name_ids)-len(unambig)
    q["n_regions_with_complete_vascular_non_invasive_only_list"]=len(complete)
    g=gpd.read_file("zip://"+str((root/"glonaf_2024_nol.zip").resolve()))
    if g.crs is None:raise ValueError("Original region CRS missing")
    g=g.to_crs("EPSG:4326")
    geomname=g.geometry.name
    rids=[region_key(row.drop(labels=[geomname]).to_dict(),by_code,by_obj) for _,row in g.iterrows()]
    g["source_region_id"]=rids
    g=g[["source_region_id",geomname]].copy()
    g[geomname]=g.geometry.make_valid()
    g=g.dissolve(by="source_region_id").reset_index()
    if len(g)<900 or g.geometry.is_empty.any():raise ValueError("Original nonoverlap polygon coverage unexpected")
    points,meta=valid_source_sample(json.loads(Path(sample).read_text()))
    land,_=source_continents(ne)
    geo=Geod(ellps="WGS84")
    by=defaultdict(Counter)
    for r in points:
        lon=r["source_original_2026_WGS84_track_point_lon"];lat=r["source_original_2026_WGS84_track_point_lat"]
        cont=land_source_continent(lat,lon,land)
        counters=[by["all_world"],by[cont]]
        for c in counters:c["all_original_draws"]+=1
        bridge=bool(r["source_overpass_or_tunnel_candidate"])
        if bridge:
            for c in counters:c["bridge_tunnel_covered_candidate"]+=1
        p=Point(lon,lat)
        hits=[int(i) for i in g.sindex.query(p,predicate="covered_by")]
        if len(hits)!=1:
            for c in counters:c["no_unique_region"]+=1
            continue
        idx=hits[0];rid=g.iloc[idx]["source_region_id"];geom=g.geometry.iloc[idx]
        n=len(region_names[rid]&unambig)
        for c in counters:
            c["unique_GloNAF_region"]+=1
            if n:c["region_has_documented_exact_species"]+=1
            if rid in complete:c["complete_vascular_list_region"]+=1
        for km in (2,5):
            w=classify_window(geo,lon,lat,km)
            if w is not None and geom.covers(w):
                for c in counters:
                    c[f"unique_region_covers_approx_{km}km_window"]+=1
                    if not bridge and n and rid in complete:
                        c[f"strict_source_support_{km}km_window"]+=1
    assert by["all_world"]["all_original_draws"]==2048
    assert by["all_world"]["unique_GloNAF_region"]+by["all_world"]["no_unique_region"]==2048
    return {"version":"v378a","source_DOI":"10.5281/zenodo.17105725",
        "original_source_md5":MD5,"source_rail_sample_seed":meta["sample_seed_SHA256"],
        "source_nonoverlap_regions_n":len(g),"original_taxon_list_QA":dict(q),
        "source_counts_by_continent_and_world":{k:dict(v) for k,v in by.items()},
        "GloNAF_nonlisted_NOT_native":True,"approximate_geodesic_window_polygon_not_exact_circle":True,
        "only_source_coverage_NOT_ecological_effect_or_plant_record_presence":True}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--glonaf",type=Path,required=True);p.add_argument("--sample",required=True)
    p.add_argument("--natural-earth",required=True);p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();v=analyze(a.glonaf,a.sample,a.natural_earth)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(v,indent=2));print(json.dumps(v,indent=2))
if __name__=="__main__":main()

#!/usr/bin/env python3
"""v379: full original 2025 plant replay for regional documented naturalized fraction."""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path
import geopandas as gpd
import numpy as np
from shapely.geometry import box
from pyproj import Transformer
from audit_global_GloNAF_rail2048_v378a import verified,load,unique_map,norm,valid_name,region_key,MD5,audit_documented_rows,read_original_nol
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
 source_client, open_original_source, taxonomy_source, source_ancestry_ids,
 source_coarse_cell, OBS_ETAG,OBS_BYTES,TAXA_ETAG,TAXA_BYTES,MAX_UNCERT_M,
 YEAR,OBS_QUALITY,X_START,Y_START)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
 verify_original_sources,source_climate_lookup,original_cell_rows,continent_holdout,
 RICH_BASE,MIN_CONTINENT_RICH_CELLS,MIN_CONTINENTS,GRID_WIDTHS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import (
 read_source_builtup,inject_built,RICH_BUILT,ORIGINAL_NE_SHA)
from generate_global_equal_area_climate_land_grid_v351 import (
 WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file,CONTINENTS)

PADS=(0.3,0.0)
VASCULAR_ID="211194"
MIN_VASCULAR_SPECIES=20
MIN_REGION_CELLS=5
MIN_REGION_BLOCKS=2
REGIONAL_BASE=RICH_BASE+(RICH_BUILT,"log_species","log_documented_regional_species_pool")
REGIONAL_RAIL=REGIONAL_BASE+("log_rail_km",)

def reconcile_nol_geometry(g,regions,bycode,byobj,header):
    """Outcome-blind region reconciliation. Never guess an ambiguous polygon."""
    txids={bycode[k] for k in header if k in bycode}
    obj=[]; direct=[]
    for _,r in g.iterrows():
        oid=byobj.get(norm(r["OBJIDsic"]))
        raw=norm(r["IDregion"])
        # IDregion may store the original numeric region.id or its code.
        opts={x for x in (raw if raw in regions else None,bycode.get(raw)) if x is not None}
        direct.append(next(iter(opts)) if len(opts)==1 else None)
        obj.append(oid)
    p=set(direct)-{None};q=set(obj)-{None}
    audit={"header_columns_n":len(header),"matrix_region_codes_n":len(txids),
        "header_nonregion_columns":[x for x in header if x not in bycode],
        "original_shape_rows":len(g),"OBJIDsic_distinct_matches":len(q),
        "IDregion_distinct_matches":len(p),
        "IDregion_OBJIDsic_disagreement_rows":sum(x!=y for x,y in zip(direct,obj)),
        "IDregion_matrix_only_n":len(txids-p),"IDregion_geometry_only_n":len(p-txids),
        "OBJIDsic_matrix_only_n":len(txids-q),"OBJIDsic_geometry_only_n":len(q-txids)}
    print(json.dumps({"source_region_key_audit":audit}),flush=True)
    if len(header)!=987 or len(txids)!=984 or len(g)!=984:
        raise ValueError("Pinned original v2.02 header/geometry row counts changed")
    if p==txids and len(p)==len(g) and all(direct):
        # Direct field gives a unique bijection to the publisher matrix.
        ids=direct
        audit["method"]="IDregion_bijection_to_publisher_matrix"
    else:
        # If the direct field is unavailable or fails the independent
        # bijection check, use the previously source-verified OBJIDsic
        # mapping ONLY for uniquely assigned matrix-region rows, never
        # duplicate or contradictory rows. Do not merge duplicate polygons.
        counts=Counter(obj)
        ids=[o if o in txids and counts[o]==1 and (d is None or d==o) else None
             for d,o in zip(direct,obj)]
        audit["method"]="unique_nonconflicting_OBJIDsic_intersection"
        if sum(x is None for x in ids)>2:
            raise ValueError("More than two original region polygons unresolved; inspect source mapping")
    audit["excluded_geometry_rows"]=sum(x is None for x in ids)
    audit["matrix_regions_without_assigned_geometry"]=len(txids-(set(ids)-{None}))
    g=g.copy()
    # Preserve unresolved polygons as blockers when testing full-window
    # containment, with no biological region ID assigned.
    g["region_id"]=ids
    print(json.dumps({"source_region_reconciliation":audit}),flush=True)
    return g,audit

def source_lists(root):
    verified(root)
    regions=unique_map(load(root/"glonaf_region.csv",{"id","code"}),"id")
    lists=unique_map(load(root/"glonaf_list.csv",{"id","region_id"}),"id")
    taxa=unique_map(load(root/"glonaf_taxon_wcvp.csv",{"id","accepted_plant_name_id","taxa_accepted","wcvp_taxon_rank"}),"id")
    flora=load(root/"glonaf_flora2.csv",{"taxon_wcvp_id","region_id","list_id"})
    tableqa=audit_documented_rows(regions,lists,taxa,flora)
    names=defaultdict(set); mapped=defaultdict(set)
    # Detect ambiguous exact names across the whole original standardized archive,
    # rather than across only the subset of naturalization entries.
    for t in taxa.values():
        n=norm(t["taxa_accepted"]); aid=norm(t["accepted_plant_name_id"])
        if valid_name(n) and norm(t["wcvp_taxon_rank"]).lower()=="species" and aid.lower() not in ("","na","nan"):
            names[n].add(aid)
    for r in flora:
        rid,lid,tid=(norm(r[k]) for k in ("region_id","list_id","taxon_wcvp_id"))
        if tid not in taxa and tid.lower() in ("","na","nan"):continue
        if rid not in regions or lid not in lists or tid not in taxa or norm(lists[lid]["region_id"])!=rid:
            raise ValueError("Original region-list-taxon IDs not consistent")
        t=taxa[tid]; n=norm(t["taxa_accepted"]); aid=norm(t["accepted_plant_name_id"])
        if valid_name(n) and norm(t["wcvp_taxon_rank"]).lower()=="species" and len(names[n])==1 and aid in names[n]:
            mapped[rid].add(n)
    if len(set().union(*mapped.values()))<10000:raise ValueError("Exact original GloNAF species coverage unexpectedly small")
    codes=unique_map(list(regions.values()),"code")
    bycode={k:norm(v["id"]) for k,v in codes.items()}
    objfield=next((k for k in next(iter(regions.values())) if k.lower()=="objidsic"),None)
    byobj={}
    if objfield:
        for k,v in regions.items():
            obj=norm(v[objfield])
            if obj and obj.lower() not in ("na","nan"):
                if obj in byobj:raise ValueError("Original OBJIDsic ambiguity")
                byobj[obj]=k
    g=read_original_nol(root)
    print(json.dumps({"original_nonoverlap_geometry_rows":len(g),"columns":list(g.columns),"crs":str(g.crs)}),flush=True)
    if g.crs is None:raise ValueError("GloNAF source missing CRS")
    g=g.to_crs("EPSG:4326")
    geom=g.geometry.name
    # Cross-check original polygon field identity against the publisher
    # matrix header rather than requiring a guessed 985-region count.
    txr=root/"glonaf_TxR.csv"
    with txr.open("rb") as f:
        if hashlib.file_digest(f,"md5").hexdigest()!="af2b5ea09f3fe6c677c7a2ceadd2b543":
            raise ValueError("Original GloNAF TxR MD5 mismatch")
    import csv
    with txr.open(encoding="utf-8-sig",newline="") as f: header=next(csv.reader(f))
    g,regionqa=reconcile_nol_geometry(g,regions,bycode,byobj,header)
    g=g[["region_id",geom]].copy()
    g[geom]=g.geometry.make_valid()
    if g.geometry.is_empty.any():
        raise ValueError("Empty source nonoverlap region geometry")
    # Source-supported geometry alone does not prove there are accepted
    # binomial species matches in the direct-region lists. Such regions
    # remain in source geometry / denominator, but cannot enter a fraction.
    no_direct_names=sum(r is None or not mapped[r] for r in g["region_id"])
    print(json.dumps({"nonoverlap_regions_without_direct_exact_species_pool":no_direct_names}),flush=True)
    return g,mapped,{"publisher_original_regions":len(regions),"nonoverlap_regions":len(g),
        "all_known_exact_species":len(set().union(*mapped.values())),
        "ambiguous_exact_names_archive":sum(len(x)>1 for x in names.values()),
        "publisher_nonoverlap_matrix_header_reconciliation":regionqa,"publisher_table_audit":tableqa,
        "regions_without_direct_exact_species_pool":no_direct_names}

def cell_region(key,width,g,inverse,pad):
    x=X_START+key[0]*width;y=Y_START+key[1]*width
    # EPSG:6933 grid parallels/meridians map monotonically to WGS84.
    a,b=inverse.transform(x,y);c,d=inverse.transform(x+width,y+width)
    if not all(map(math.isfinite,(a,b,c,d))) or not(-180<=a<c<=180 and -90<=b<d<=90):
        return None
    a-=pad;b-=pad;c+=pad;d+=pad
    if a<=-180 or c>=180 or b<=-90 or d>=90:return None
    area=box(a,b,c,d)
    hits=g.sindex.query(area,predicate="covered_by")
    if len(hits)!=1:return None
    # A full-cell region hit must not also intersect an unresolved
    # original polygon or a second region with positive area.
    i=int(hits[0]);rid=g.iloc[i]["region_id"]
    if rid is None:return None
    for j in g.sindex.query(area,predicate="intersects"):
        if int(j)!=i and g.geometry.iloc[int(j)].intersection(area).area>0:
            return None
    return rid

def grid_counts(plant,width):
    return plant["source_grids"][str(width)]["source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts"]
def grid_key(r):
    return r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]

def design(plant,rail,built,ne,climate,g):
    if sha_file(ne)!=ORIGINAL_NE_SHA or sha_file(climate)!=WORLDCLIM_SHA:
        raise ValueError("Original independent global map or climate changed")
    land,_=source_continents(ne)
    raster,mem,_=worldclim_tif(climate)
    inv=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
    ghsl,_=read_source_builtup(built)
    panels={};qa={}
    try:
        scale,_=climate_scale(raster)
        for w in GRID_WIDTHS:
            ps=grid_counts(plant,w)
            ck={grid_key(r) for r in ps}
            geo,_=source_climate_lookup(land,raster,scale,ck,w)
            rs=rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)]
            elig,_,_=original_cell_rows(ps,rs,geo)
            original_n={100000:3455,200000:1442}[w]
            if len(elig)!=original_n:raise ValueError("Not the original v376c eligible plant-cell denominator")
            eligible,_=inject_built(elig,ghsl[str(w)])
            if len(eligible)!=original_n:raise ValueError("Not the original complete GHSL eligible panel")
            for pad in PADS:
                rows=[]
                for a in eligible:
                    rid=cell_region(a["key"],w,g,inv,pad)
                    if rid is not None:rows.append(dict(a,region=rid))
                panels[w,pad]=rows
                qa[f"{w}/{pad}"]={"original_global_denominator":original_n,
                    "original_continental_denominator":dict(Counter(x["continent"] for x in eligible)),
                    "assigned_single_region_by_continent":dict(Counter(x["continent"] for x in rows)),
                    "assigned_cells":len(rows),"assigned_regions":len({x["region"] for x in rows})}
    finally:raster.close();mem.close()
    print(json.dumps({"pre_plant_replay_region_assignment":qa}),flush=True)
    return panels,qa

def vascular_names(client,canonical):
    rows,raw,stream=open_original_source(client,"taxa.csv.gz",TAXA_ETAG,TAXA_BYTES)
    names={};root_ok=False
    try:
        for row in rows:
            tid=int(row["taxon_id"])
            if str(tid)==VASCULAR_ID:
                root_ok=(row["name"]=="Tracheophyta" and row["rank"].strip().lower()=="phylum")
            if canonical.get(tid)==tid and VASCULAR_ID in source_ancestry_ids(row["ancestry"]):
                names[tid]=norm(row["name"])
    finally:stream.close()
    if raw.bytes!=TAXA_BYTES or not root_ok or len(names)<10000:
        raise ValueError("Original source vascular root or taxonomy census invalid")
    return names,raw.hash.hexdigest()

def replay(client,canonical,vasc,plant,panels):
    selected={w:set() for w in GRID_WIDTHS}
    for (w,pad),entries in panels.items():selected[w].update(x["key"] for x in entries)
    counts={w:Counter() for w in GRID_WIDTHS}
    sets={w:defaultdict(set) for w in GRID_WIDTHS}
    # No individual coordinates/observer identities are retained or uploaded.
    transform=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    qa=Counter()
    try:
        for row in rows:
            qa["all_rows"]+=1
            if qa["all_rows"]%10000000==0:print(json.dumps(dict(qa)),flush=True)
            if not str(row["observed_on"]).startswith(YEAR+"-"):continue
            if str(row["quality_grade"]).strip().lower().replace(" ","_").replace("-","_")!=OBS_QUALITY:continue
            try:tid=canonical.get(int(row["taxon_id"]))
            except (ValueError,TypeError):tid=None
            if tid is None:continue
            try:unc=float(row["positional_accuracy"]);lon=float(row["longitude"]);lat=float(row["latitude"])
            except (TypeError,ValueError):continue
            if not math.isfinite(unc) or unc<0 or unc>MAX_UNCERT_M:continue
            keys={w:source_coarse_cell(lon,lat,transform,w) for w in GRID_WIDTHS}
            if any(k is None for k in keys.values()):continue
            qa["all_original_eligible_plant_observations"]+=1
            if tid in vasc:qa["all_vascular_observations"]+=1
            for w,k in keys.items():
                counts[w][k]+=1
                if tid in vasc and k in selected[w]:sets[w][k].add(tid)
    finally:stream.close()
    h=plant["source_summary"]
    if raw.bytes!=OBS_BYTES or raw.hash.hexdigest()!=h["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"]:
        raise ValueError("Full raw observation source fingerprint differs from v369")
    old=h["source_global_photo_observations_QA"]
    if qa["all_rows"]!=old["original_global_iNat_photo_linked_rows_read"] or qa["all_original_eligible_plant_observations"]!=old["source_year_2025_public_research_grade_georeferenced_plant_obs"]:
        raise ValueError("v369 all original record counts not reproduced")
    for w in GRID_WIDTHS:
        expected=Counter({grid_key(r):r["source_year_2025_photo_linked_plants_observation_count"] for r in grid_counts(plant,w)})
        if counts[w]!=expected:raise ValueError("v369 exact per-cell original observation counts not reproduced")
    return sets,dict(qa)

def within_region(rows):
    groups=defaultdict(list)
    for r in rows:groups[r["region"]].append(r)
    group_ok={k for k,a in groups.items() if len(a)>=MIN_REGION_CELLS and
              len({(r["continent"],r["block"]) for r in a})>=MIN_REGION_BLOCKS and
              len({r["log_rail_km"] for r in a})>1}
    a=[r for r in rows if r["region"] in group_ok]
    result={"n_regions":len(group_ok),"n_cells":len(a),
            "n_blocks":len({(r["continent"],r["block"]) for r in a}),"coefficient_pp_per_log1p_rail_km":None}
    if len(group_ok)<2 or len(a)<30 or result["n_blocks"]<5:return result
    counts=Counter((r["continent"],r["block"]) for r in a)
    w=np.array([1/counts[r["continent"],r["block"]] for r in a]);w=w/w.sum()
    fields=RICH_BASE+(RICH_BUILT,"log_species","log_rail_km","known_fraction")
    X=np.array([[r[k] for k in fields] for r in a],dtype=float)
    ids=defaultdict(list)
    for j,r in enumerate(a):ids[r["region"]].append(j)
    for ix in ids.values():
        X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    # weighted projection of the rail predictor and response off nuisance
    # columns; a rank-deficient nuisance design is not by itself fatal.
    Z=X[:,:-2]*np.sqrt(w)[:,None]
    rail=X[:,-2]*np.sqrt(w);y=X[:,-1]*np.sqrt(w)
    res=lambda u:u-Z@np.linalg.lstsq(Z,u,rcond=1e-10)[0]
    rr,yy=res(rail),res(y)
    support=float(rr@rr)
    result["rail_residual_weighted_variance"]=support
    if support<=1e-10:return result
    result["coefficient_pp_per_log1p_rail_km"]=float(100*(rr@yy)/support)
    return result

def score(panels,sets,names,listed):
    out=[]
    for (w,pad),items in panels.items():
        rows=[];qa=Counter()
        for r in items:
            ids=sets[w][r["key"]]
            if len(ids)>r["species"]:raise ValueError("Vascular richness exceeds original all-plant richness")
            if len(ids)<MIN_VASCULAR_SPECIES:
                qa["fails_min_vascular_species"]+=1;continue
            documented=listed[r["region"]]
            if not documented:
                qa["no_direct_exact_documented_species_pool"]+=1
                continue
            # Numerator is the count of source canonical species IDs
            # exactly matching an accepted, region-documented name.
            n=sum(names[i] in documented for i in ids)
            rows.append(dict(r,known_fraction=n/len(ids),
                         log_documented_regional_species_pool=math.log1p(len(documented))))
        cv=continent_holdout(rows,REGIONAL_BASE,REGIONAL_RAIL,"known_fraction",True,
                             MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,
                             "whole_landscape_recorded_vascular_species_fraction_documented_GloNAF")
        out.append({"width_m":w,"region_guard_degrees":pad,"after_vascular_gate":len(rows),
            "after_vascular_gate_by_continent":dict(Counter(r["continent"] for r in rows)),
            "excluded":dict(qa),"prediction":cv,"within_region_world":within_region(rows),
            "within_region_by_continent":{c:within_region([r for r in rows if r["continent"]==c]) for c in CONTINENTS},
            "mean_cell_fraction":float(np.mean([r["known_fraction"] for r in rows])) if rows else None})
    return out
def main():
    p=argparse.ArgumentParser()
    for key in ("plants","rail","built","ne","climate","glonaf","preflight","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    plant=json.loads(a.plants.read_text());rail=json.loads(a.rail.read_text())
    verify_original_sources(plant,rail)
    g,listed,lqa=source_lists(a.glonaf)
    panels,qa=design(plant,rail,a.built,a.ne,a.climate,g)
    a.preflight.parent.mkdir(parents=True,exist_ok=True)
    a.preflight.write_text(json.dumps({"version":"v379-preflight","source_lists":lqa,"design":qa,"ecological_effect_estimated":False},indent=2))
    # Require at least three eligible continental folds on an original
    # panel before paying for the full source replay. Report unsupported
    # panels rather than changing site selection.
    if not any(sum(n>=MIN_CONTINENT_RICH_CELLS for n in x["assigned_single_region_by_continent"].values())>=MIN_CONTINENTS for x in qa.values()):
        raise ValueError("Insufficient pre-outcome regional source assignment for frozen continental gates")
    client=source_client()
    canonical,_,_,_,sha=taxonomy_source(client)
    names,sha2=vascular_names(client,canonical)
    if sha!=sha2 or sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:
        raise ValueError("Original taxonomy SHA differs from v369")
    sets,sourceqa=replay(client,canonical,names,plant,panels)
    out={"version":"v379-original-global-coarse-documented-naturalized-species-fraction",
         "source_glonaf_DOI":"10.5281/zenodo.17105725","source_GloNAF_MD5":MD5,
         "source_list_QA":lqa,"original_v369_exact_replay_QA":sourceqa,"design":qa,
         "panels":score(panels,sets,names,listed),
         "limits":["GloNAF nonlisted does not mean native","Whole landscape photo-linked species detections, not railway flora",
          "No direct precise observation geoprivacy audit; region guard is a sensitivity only",
          "Road, population, forest and true survey access not controlled",
          "No causal transport, establishment, or local distance-decay inference",
          "No observer IDs, species IDs or individual original coordinates published"]}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False))
    print(json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()

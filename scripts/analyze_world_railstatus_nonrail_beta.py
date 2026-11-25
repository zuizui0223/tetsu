#!/usr/bin/env python3
"""Compare plant-bird composition coupling in mapped-rail and rail-free landscapes.
Same source cells, no raw scan; descriptive common-block contrasts, not causal.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import defaultdict,Counter
import analyze_world_railstatus_observer_split_connectivity as s
from analyze_world_railstatus_direct_boundary_beta import CELL_FIELDS,CATS
from analyze_world_railstatus_plant_insect_bird_beta import (
    attach_agb,join,setup,fields,score_fit,model,combine,interval,EXPECTED,ORIENT)
GROUPS=("plants","birds","insects","fungi")
PB=("plants","birds")

def classify(x,y):
    if min(x["active_normal_km"],y["active_normal_km"])>0:return "both_active"
    if all(r[c+"_km"]==0 for r in (x,y) for c in CATS):return "both_four_status_zero"
    return "mixed_or_other"

def unrestricted_pairs(pairs,rows,edges):
    by={r["key"]:r for r in rows};out=[]
    for p in pairs:
        a=tuple(p["source_cell_A_original_grid_ix_iy"]);b=tuple(p["source_cell_B_original_grid_ix_iy"])
        if b[0]-a[0]+b[1]-a[1]!=1 or sum(abs(x-y) for x,y in zip(a,b))!=1:raise ValueError("Nonadjacent pair")
        x=by.get(a);y=by.get(b)
        if x is None or y is None or x["continent"]!=y["continent"] or x["block"]!=y["block"]:continue
        if x["species"]!=p["source_grid_A_original_plant_species_n"] or y["species"]!=p["source_grid_B_original_plant_species_n"]:raise ValueError("Original species count mismatch")
        shared=p["source_shared_plants_species_n"];union=p["source_union_plants_species_n"]
        if union!=x["species"]+y["species"]-shared or union<=0:raise ValueError("Jaccard components invalid")
        j=float(p["source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells"])
        if not 0<=j<=1 or abs(j-shared/union)>1e-7:raise ValueError("Jaccard mismatch")
        r={"a":a,"b":b,"continent":x["continent"],"block":x["block"],"jaccard":j,
           "simpson":shared/min(x["species"],y["species"]),"horizontal":float(b[0]!=a[0]),
           "both_inactive":min(x["inactive_normal_km"],y["inactive_normal_km"])>0}
        for k in CELL_FIELDS:r["mean_"+k]=(x[k]+y[k])/2;r["diff_"+k]=abs(x[k]-y[k])
        for c in CATS:
            n=edges[c].get((a,b),0);r["cross_"+c]=float(n>0);r["log_cross_"+c]=math.log1p(n)
            if n>0 and min(x[c+"_km"],y[c+"_km"])<=0:raise ValueError("Geometry crossing inconsistent with prior rail km")
        r["rail_class"]=classify(x,y)
        out.append(r)
    return out

def key(r):return r["continent"],r["block"]

def shared_support(a,b):
    ca=Counter(key(r) for r in a);cb=Counter(key(r) for r in b)
    common={k for k in ca.keys()&cb.keys() if min(ca[k],cb[k])>=2}
    return [r for r in a if key(r) in common],[r for r in b if key(r) in common]

def coupling(rows,epoch=2024,metric="jaccard",flex=True,endpoint=False):
    rr=setup(rows,PB,epoch,flex);base=fields(PB,endpoint,flex)
    qa={"n_pairs":len(rr),"n_blocks":len({key(r) for r in rr}),"continents":dict(Counter(r["continent"] for r in rr))}
    f={}
    for i,j in ORIENT:
        p=f"plants_{i}_{metric}";y=f"birds_{j}_{metric}"
        q=score_fit(rr,base+(p,),y,(p,),endpoint)
        if q is None:return dict(qa,status="insufficient support or rank"),None
        f[i,j]=q[p]
    pooled=combine(f,{ij:.5 for ij in ORIENT})
    return dict(qa,status="estimated",effect=interval(pooled)),pooled

def comparison(a,b,**kw):
    aa,bb=shared_support(a,b)
    x,fx=coupling(aa,**kw);y,fy=coupling(bb,**kw)
    q={"rail":x,"four_status_zero":y,"n_common_blocks":len({key(r) for r in aa})}
    if fx is not None and fy is not None:
        q["rail_minus_zero"]=interval(combine({"rail":fx,"zero":fy},{"rail":1.,"zero":-1.}))
    return q

def panels(rows):
    a=[r for r in rows if r["rail_class"]=="both_active"]
    b=[r for r in rows if r["rail_class"]=="both_four_status_zero"]
    return {"class_counts":dict(Counter(r["rail_class"] for r in rows)),
        "native":{"rail_flexible":coupling(a)[0],"zero_flexible":coupling(b)[0],
                  "rail_linear":coupling(a,flex=False)[0],"zero_linear":coupling(b,flex=False)[0],
                  "rail_simpson":coupling(a,metric="simpson",flex=False)[0],
                  "zero_simpson":coupling(b,metric="simpson",flex=False)[0]},
        "common_blocks_flexible":comparison(a,b),"common_blocks_linear":comparison(a,b,flex=False),
        "common_blocks_simpson":comparison(a,b,metric="simpson",flex=False),
        "common_blocks_AGB2020":comparison(a,b,epoch=2020),
        "common_blocks_endpoint":comparison(a,b,endpoint=True,flex=False)}

def replay(rr,old):
    if len(rr)!=old["n_pairs"]:raise ValueError("Rail-positive comparison pair support changed")
    for label,kw in (("primary",{}),("flexible_effort",{"flex":True}),("simpson",{"metric":"simpson"})):
        z=model(rr,PB,**kw);v=old[label]
        if z["status"]!=v["status"] or z["n_pairs"]!=v["n_pairs"]:raise ValueError("Prior status changed")
        if z["status"]=="estimated":
            for name,e in z["effects"].items():
                q=v["effects"][name]
                if abs(e["estimate"]-q["estimate"])>1e-9 or max(abs(x-y) for x,y in zip(e["interval95_approx"],q["interval95_approx"]))>1e-9:
                    raise ValueError("Prior rail-positive estimate or interval changed")

def main():
    p=argparse.ArgumentParser()
    for k in ("beta_model","split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    src,prior,edge,agb,am,plant,rail,status,tree,road,pop=map(s.read,(a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    beta=s.read(a.beta_model)
    if beta["version"]!="railstatus-disjoint-community-composition-coupling-v1" or beta["biomass_source_SHA256"]!=agb["source_sha256"] or any(beta["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","edges","agb","agb_model")):raise ValueError("Wrong prior composition model")
    old={(r["width_m"],r["land"]):r for r in beta["panels"]}
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_gate_observations_observers_species"]!=[25,5,10] or src["split_seed"]!="tetsu-railstatus-observer-split-v1":raise ValueError("Wrong frozen split")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(src["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Frozen source mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):raise ValueError("Wrong geometry")
    if tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong road/tree")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong climate/map/pop")
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    results=[]
    try:
        scale,_=s.climate_scale(ras)
        for w in s.GRID_WIDTHS:
            grid=prior["grids"][str(w)];sp=src["pairs"][str(w)]
            keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for g in GROUPS for r in grid[g]}
            clim,_=s.source_climate_lookup(land,ras,scale,keys,w)
            tr=s.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=s.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=s.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=s.make_lookup(agb["grids"]["2020"][str(w)]);ag24=s.make_lookup(agb["grids"]["2024"][str(w)])
            el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(w)][c]} for c in s.CATS}
            raw={g:{(tuple(p["a"]),tuple(p["b"])):p for p in sp[g]} for g in GROUPS};geo={}
            for g in GROUPS:
                rows,_,_=s.original_cell_rows([s.convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=s.inject_built(rows,built[str(w)]);rows,_=s.inject_status(rows,status,w);geo[g]=rows
            for gate in s.LAND_THRESHOLDS:
                groups={};ag={}
                for g in GROUPS:
                    rr,_=s.enrich_cells(geo[g],rd,tr,w,gate);rr,_=s.join_population(rr,pp,rd)
                    full=unrestricted_pairs([s.adapt_pair(p) for p in sp[g]],rr,el)
                    groups[g]=[z for r in full if (z:=s.split_pair(r,raw[g][r["a"],r["b"]])) is not None]
                    for k,v in attach_agb(rr,ag20,ag24,rd).items():
                        if k in ag and ag[k]!=v:raise ValueError("AGB grid mismatch")
                        ag[k]=v
                original={g:[r for r in groups[g] if r["rail_class"]=="both_active"] for g in GROUPS}
                if tuple(len(original[g]) for g in GROUPS)!=EXPECTED[w,gate]:raise ValueError("Exact original pair support changed")
                rr=join(groups,raw,ag,PB)
                cls={(r["a"],r["b"]):r["rail_class"] for r in groups["plants"]}
                birdcls={(r["a"],r["b"]):r["rail_class"] for r in groups["birds"]}
                for r in rr:
                    k=(r["a"],r["b"])
                    if cls[k]!=birdcls[k]:raise ValueError("Pair classification mismatch")
                    r["rail_class"]=cls[k]
                replay([r for r in rr if r["rail_class"]=="both_active"],old[w,gate]["plants_birds"])
                z={"width_m":w,"land":gate,"source_pair_counts":{g:len(groups[g]) for g in GROUPS},**panels(rr)}
                results.append(z);print(json.dumps({"nonrail_beta_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-disjoint-nonrail-composition-control-v1","source_replay":checks,
        "input_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("beta_model","split","edges","agb","agb_model")},
        "panels":results,"primary":"100km land80 2025, direct both-active minus both-zero four-status plant-bird Jaccard slope contrast on common geographical blocks with >=2 pairs per category, flexible effort, observer halves disjoint across taxa. Native slopes and support diagnostics retained.",
        "limits":["zero mapped length of four included rail categories, not complete absence of any railway infrastructure or historical rail","geographically selected, shared blocks are not matched sites or randomized treatments","source gates favor well-sampled landscapes","same platform, observed similarity, no causal rail effect or mediation","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

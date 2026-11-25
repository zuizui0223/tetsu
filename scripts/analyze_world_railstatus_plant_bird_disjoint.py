#!/usr/bin/env python3
"""Observer-disjoint plant/bird richness coupling on frozen biomass geography.
All four orientations use identical cells and joint effort covariates.
No raw occurrence rescan. Cross-side observers cannot share identities.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_observer_annual_connectivity import cells_from_pairs
from scan_world_railstatus_observer_annual_beta import candidate_support,pair_key
from scan_world_railstatus_observer_split_beta import MIN_FULL,MIN_SPLIT,SALT
YEARS=("2023","2024","2025")
GROUPS=("plants","birds")
GEO=tuple(k for k in b.FULL+(b.POP,) if k not in ("log_observations","log_observers"))
EFF=tuple(f"{g}_{i}_{k}" for g in GROUPS for i in (0,1) for k in ("observations","observers"))
FLEX=tuple(f"{k}_{p}" for k in EFF for p in ("sq","cu"))+tuple(f"{g}_{i}_effort_interaction" for g in GROUPS for i in (0,1))
ORIENTATIONS=((0,0),(0,1),(1,0),(1,1))
TERMS={"cross":{(0,1):.5,(1,0):.5},"same":{(0,0):.5,(1,1):.5},
       "same_minus_cross":{(0,0):.5,(1,1):.5,(0,1):-.5,(1,0):-.5}}
EXPECTED={100000:{.5:(2636,.29977),.8:(2286,.30797)},200000:{.5:(1191,.45569),.8:(990,.43593)}}

def extract(pairs):
    out={}
    for p in pairs:
        if "split_similarity" not in p:continue
        for s in ("a","b"):
            k=tuple(p[s]);f=tuple(p["full_"+s]);h=tuple(tuple(x) for x in p["split_counts_"+s])
            if len(f)!=3 or len(h)!=2 or any(len(x)!=3 for x in h):raise ValueError("Wrong observer split schema")
            if any(x<y for x,y in zip(f,MIN_FULL)) or any(x<y for z in h for x,y in zip(z,MIN_SPLIT)):raise ValueError("Frozen source gate violated")
            if h[0][1]+h[1][1]!=f[1] or h[0][0]+h[1][0]>f[0] or not max(h[0][2],h[1][2])<=f[2] or any(x[2]>x[0] or x[1]>x[0] for x in h):raise ValueError("Split count inconsistency")
            if k in out and out[k]!=(f,h):raise ValueError("Repeated endpoint counts differ")
            out[k]=(f,h)
    return out

def attach(rows,maps):
    out=[]
    for r in rows:
        k=r["key"]
        if not all(k in maps[g] for g in GROUPS):continue
        z=dict(r)
        for g in GROUPS:
            full,hs=maps[g][k]
            expected=(r["observations"],r["observers"],r["species"]) if g=="plants" else (r["bird_observations"],r["bird_observers"],r["bird_species"])
            if tuple(expected)!=full:raise ValueError("Exact full count replay failed")
            for i,counts in enumerate(hs):
                for name,n in zip(("observations","observers","species"),counts):
                    f=f"{g}_{i}_{name}";z[f]=math.log1p(n)
                    if name!="species":z[f+"_sq"]=z[f]**2;z[f+"_cu"]=z[f]**3
                z[f"{g}_{i}_effort_interaction"]=z[f"{g}_{i}_observations"]*z[f"{g}_{i}_observers"]
        out.append(z)
    return out

def linear(rows,fields,target,term):
    groups=defaultdict(list)
    for i,r in enumerate(rows):groups[r["continent"],r["block"]].append(i)
    if len(rows)<100 or len(groups)<20 or len(rows)-len(groups)-len(fields)<30:return None
    w=b.block_weights(rows)
    X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite split design")
    for ii in groups.values():X[ii]-=np.average(X[ii],axis=0,weights=w[ii])
    X*=np.sqrt(w)[:,None];A=X[:,:-1];y=X[:,-1]
    scale=np.sqrt(np.sum(A*A,axis=0));keep=scale>1e-10;names=[k for k,v in zip(fields,keep) if v]
    if term not in names:return None
    Z=A[:,keep]/scale[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return None
    inv=np.linalg.inv(Z.T@Z);coef=inv@(Z.T@y);e=y-Z@coef;j=names.index(term)
    scores={k:float((inv@(Z[ii].T@e[ii]))[j]/scale[keep][j]) for k,ii in groups.items()}
    return float(coef[j]/scale[keep][j]),scores

def interval(x):
    beta,s=x;n=len(s);se=math.sqrt(n/(n-1)*sum(v*v for v in s.values()))
    return {"estimate":beta,"interval95_approx":[beta-1.96*se,beta+1.96*se],"n_blocks":n}

def combine(fits,weights):
    keys=set(next(iter(fits.values()))[1])
    if any(set(v[1])!=keys for v in fits.values()):raise ValueError("Paired blocks changed")
    return (sum(weights[k]*fits[k][0] for k in weights),
            {j:sum(weights[k]*fits[k][1][j] for k in weights) for j in sorted(keys)})

def one(rows,year=2024,flex=False,keep_scores=False):
    rr=b.with_epoch(b.local_support(rows),year)
    fields=GEO+(b.BIO,)+EFF+(FLEX if flex else ())
    f={}
    for i,j in ORIENTATIONS:
        term=f"plants_{i}_species"
        f[i,j]=linear(rr,fields+(term,),f"birds_{j}_species",term)
    out={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),"continents":dict(Counter(r["continent"] for r in rr))}
    if any(v is None for v in f.values()):return dict(out,status="insufficient support or rank"),None
    c={k:combine(f,w) for k,w in TERMS.items()}
    if abs(c["same"][0]-c["cross"][0]-c["same_minus_cross"][0])>1e-10:raise ValueError("Contrast identity")
    return dict(out,status="estimated",orientations={f"{i}{j}":interval(v) for (i,j),v in f.items()},
                contrasts={k:interval(v) for k,v in c.items()}),c if keep_scores else None

def evaluate(rows):
    rr=b.local_support(rows)
    return {"n_split_supported":len(rr),"primary":one(rr)[0],
            "flexible_effort":one(rr,flex=True)[0],"AGB2020_same_cells":one(rr,2020)[0],
            "full_counts_same_cells":b.models(rr,2024)["bird_plant_richness_association"],
            "high_effort":one([r for r in rr if min(r["observations"],r["bird_observations"])>=200 and min(r["observers"],r["bird_observers"])>=20])[0],
            "leave_continent_out":{c:one([r for r in rr if r["continent"]!=c])[0] for c in sorted({r["continent"] for r in rr})}}

def temporal(by):
    common=set.intersection(*({r["key"] for r in by[y]} for y in YEARS))
    by={y:b.local_support([r for r in by[y] if r["key"] in common]) for y in YEARS}
    if any({r["key"] for r in by[y]}!={r["key"] for r in by["2025"]} for y in YEARS):raise ValueError("Annual geography changed")
    output={"n_common_cells":len(by["2025"])}
    for name,flex in (("linear",False),("flexible",True)):
        z={};f={}
        for y in YEARS:z[y],f[y]=one(by[y],flex=flex,keep_scores=True)
        result={"years":z}
        if all(v is not None for v in f.values()):
            # 2023/24 mean excludes 2025 discovery; paired by geography.
            result["earlier_mean"]={t:interval(combine({y:f[y][t] for y in YEARS[:2]},{y:.5 for y in YEARS[:2]})) for t in TERMS}
            result["all_year_mean"]={t:interval(combine({y:f[y][t] for y in YEARS},{y:1/3 for y in YEARS})) for t in TERMS}
        output[name]=result
    return output

def main():
    p=argparse.ArgumentParser()
    for n in ("agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or old["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or old["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass replay mismatch")
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Frozen observer split invalid")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=200000 or annual["qa"]["rows"]!=288721938 or any(len({src[k],annual[k],taxa[k]})!=1 for k in ("obs_sha256","taxa_sha256")) or annual["split_seed"]!=src["split_seed"]:raise ValueError("Annual source mismatch")
    cand,_=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Frozen annual candidates changed")
    for y in YEARS[:2]:
        for g in GROUPS:
            if not {pair_key(p) for p in annual["pairs"][y][g]}<=cand[g]:raise ValueError("Annual candidates changed")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen map/climate/pop mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);raster,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(raster)
        for w in b.GRID_WIDTHS:
            grid=taxa["grids"][str(w)];yr=YEARS if w==200000 else ("2025",)
            ps={y:{g:src["pairs"][str(w)][g] if y=="2025" else annual["pairs"][y][g] for g in GROUPS} for y in yr}
            maps={y:{g:extract(ps[y][g]) for g in GROUPS} for y in yr}
            prior_counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in GROUPS}
            for g in GROUPS:
                if any(prior_counts[g].get(k)!=v[0] for k,v in maps["2025"][g].items()):raise ValueError("2025 exact source replay failed")
            keys={b.key(z) for g in GROUPS for z in grid[g]}|{k for y in yr for g in GROUPS for k in maps[y][g]}
            clim,_=b.source_climate_lookup(land,raster,scale,keys,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            geo={}
            for y in yr:
                geo[y]={}
                for g in GROUPS:
                    raw=[b.convert(z) for z in grid[g]] if y=="2025" else cells_from_pairs([p for p in ps[y][g] if "split_similarity" in p])
                    rows,_,_=b.original_cell_rows(raw,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                    rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w);geo[y][g]=rows
            for gate in b.LAND_THRESHOLDS:
                by={}
                for y in yr:
                    rows={}
                    for g in GROUPS:
                        m,_=b.enrich_cells(geo[y][g],rd,tr,w,gate);m,_=b.join_population(m,pp,rd)
                        rows[g]=[b.flexible(r) for r in m]
                    joined,_=b.join(rows,ag20,ag24,w,rd)
                    bc={r["key"]:r["species"] for r in rows["birds"]}
                    joined=[dict(r,bird_species=bc[r["key"]]) for r in joined]
                    if y=="2025":
                        n,expected=EXPECTED[w][gate];m=b.models(joined,2024)
                        q=m["bird_plant_richness_association"]["effects"]["plant_richness"]["bird_log_richness"]["estimate"]
                        if m["n"]!=n or abs(q-expected)>5.1e-6:raise ValueError("Prior full-count plant/bird coefficient not replayed")
                    by[y]=attach(joined,maps[y])
                z={"width_m":w,"land":gate,"2025":evaluate(by["2025"])}
                if w==200000:z["common_annual"]=temporal(by)
                panels.append(z);print(json.dumps({"disjoint_plant_bird_biomass_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-observer-disjoint-plant-bird-biomass-v1","source_replay":checks,"biomass_source_sha256":agb["source_sha256"],"source_artifact_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")},"panels":panels,
      "primary":"100km land80 2025, mean of two globally disjoint observer-side slopes of log(1+bird species) on log(1+plant species), all four half efforts and woody AGB + existing geography controls, local block FE; paired same-minus-cross slope diagnostic",
      "limits":["same platform and shared access/detection/environment despite disjoint identities","split-supported adjacent-pair endpoint selection; no representative global claim","same vs cross is diagnostic, not causal decomposition","2020/2024 AGB is woody, including zero ocean pixels; land80 primary","earlier annual panels same geography and static map/AGB, retrospective check","observed richness, not trophic interactions or bird biomass"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()

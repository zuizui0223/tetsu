#!/usr/bin/env python3
"""Within-cell annual observed richness coupling, disjoint observers.
Balanced 2023-25 panel; no static geographical effect estimated here.
Cell and year (or spatial block-by-year) FE remove time-invariant habitat
and common temporal trends. Paired score covariance across orientations.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,one,temporal,combine,interval)
from analyze_world_railstatus_plant_bird_coupling_modifiers import replay_check

def common_rows(by):
    maps={y:{r["key"]:r for r in by[y]} for y in YEARS}
    if any(len(maps[y])!=len(by[y]) for y in YEARS):raise ValueError("Repeated cell")
    keys=set.intersection(*(set(v) for v in maps.values()))
    rr={y:b.local_support([maps[y][k] for k in sorted(keys)]) for y in YEARS}
    kk={r["key"] for r in rr["2025"]}
    if any({r["key"] for r in rr[y]}!=kk for y in YEARS):raise ValueError("Annual geography differs")
    for k in kk:
        if len({(maps[y][k]["continent"],maps[y][k]["block"]) for y in YEARS})!=1:raise ValueError("Annual block differs")
    return rr

def fit_panel(by,fields,target,term,trend="year"):
    """Exact weighted two-way demeaning for balanced panel.
    All cell weights constant over time; each block has weight one.
    Clustered by original spatial block, including all years.
    """
    yy=tuple(by)
    if len(yy)<2:raise ValueError("Need two years")
    keys=sorted(r["key"] for r in by[yy[0]])
    if len(set(keys))!=len(keys):raise ValueError("Repeated cells")
    n=len(keys);maps={y:{r["key"]:r for r in by[y]} for y in yy}
    if any(set(m)!=set(keys) for m in maps.values()):raise ValueError("Unbalanced panel")
    groups=defaultdict(list)
    for i,k in enumerate(keys):
        r=maps[yy[0]][k];block=(r["continent"],r["block"]);groups[block].append(i)
        if any((maps[y][k]["continent"],maps[y][k]["block"])!=block for y in yy):raise ValueError("Moved cell")
    B=len(groups);T=len(yy)
    fe_rank=n+(T-1)*(1 if trend=="year" else B)
    if n<100 or B<20 or n*T-fe_rank-len(fields)<30:return None
    if trend not in ("year","block_year"):raise ValueError("Unknown trend")
    w=np.zeros(n)
    for ix in groups.values():w[ix]=1./(B*len(ix))
    X=np.array([[[maps[y][k][f] for f in (*fields,target)] for k in keys] for y in yy],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite panel")
    raw_norm=np.sqrt(np.sum(X[:,:,:-1]**2*w[None,:,None]/T,axis=(0,1)))
    # These projections commute on a balanced, time-invariant-weight panel.
    X-=X.mean(axis=0,keepdims=True)
    if trend=="year":X-=np.sum(X*w[None,:,None],axis=1,keepdims=True)
    else:
        for ix in groups.values():X[:,ix,:]-=X[:,ix,:].mean(axis=1,keepdims=True)
    if np.max(np.abs(X.mean(axis=0)))>1e-9:raise ValueError("Cell FE projection failed")
    if trend=="year":
        if np.max(np.abs(np.sum(X*w[None,:,None],axis=1)))>1e-9:raise ValueError("Year FE projection failed")
    else:
        if max(np.max(np.abs(X[:,ix,:].mean(axis=1))) for ix in groups.values())>1e-9:raise ValueError("Block-year FE projection failed")
    X*=np.sqrt(w[None,:,None]/T)
    X=X.reshape(n*T,len(fields)+1);A=X[:,:-1];y=X[:,-1]
    norm=np.sqrt(np.sum(A*A,axis=0));keep=norm>np.maximum(1e-10,raw_norm*1e-8)
    names=[k for k,v in zip(fields,keep) if v]
    if term not in names:return None
    Z=A[:,keep]/norm[keep];sv=np.linalg.svd(Z,compute_uv=False)
    if sv[-1]<1e-8*sv[0]:return None
    inv=np.linalg.inv(Z.T@Z);coef=inv@(Z.T@y);res=y-Z@coef;j=names.index(term)
    scores={}
    for k,ix in groups.items():
        ii=np.concatenate([np.asarray(ix)+t*n for t in range(T)])
        scores[k]=float((inv@(Z[ii].T@res[ii]))[j]/norm[keep][j])
    if sum(np.linalg.norm(Z[np.concatenate([np.asarray(ix)+t*n for t in range(T)]),j])>1e-8 for ix in groups.values())<20:return None
    return float(coef[j]/norm[keep][j]),scores

def model(by,trend="year",flex=False):
    n=len(next(iter(by.values())))
    out={"n_cells":n,"n_years":len(by),"n_blocks":len({(r["continent"],r["block"]) for r in next(iter(by.values()))})}
    fs={}
    for i,j in ORIENTATIONS:
        term=f"plants_{i}_species"
        fs[i,j]=fit_panel(by,EFF+(FLEX if flex else ())+(term,),f"birds_{j}_species",term,trend)
    if any(v is None for v in fs.values()):return dict(out,status="insufficient support or rank")
    return dict(out,status="estimated",contrasts={k:interval(combine(fs,w)) for k,w in TERMS.items()},
                orientations={f"{i}{j}":interval(v) for (i,j),v in fs.items()})

def analyze(by,old):
    rr=common_rows(by)
    if len(rr["2025"])!=old["n_common_cells"]:raise ValueError("Prior common panel changed")
    past=temporal(by)
    for lab in ("linear","flexible"):
        for y in YEARS:replay_check(past[lab]["years"][y],old[lab]["years"][y])
    specs={"cell_year":("year",False),"cell_year_flex":("year",True),
           "cell_block_year":("block_year",False),"cell_block_year_flex":("block_year",True)}
    out={"n_common_cells":len(rr["2025"]),"models":{k:model(rr,*v) for k,v in specs.items()}}
    # Two adjacent intervals: all 3-year common cells, no sample changes.
    out["periods"]={f"{ys[0]}_{ys[1]}":{k:model({y:rr[y] for y in ys},*v) for k,v in specs.items()} for ys in (YEARS[:2],YEARS[1:])}
    out["leave_continent_out"]={c:{k:model({y:[r for r in rr[y] if r["continent"]!=c] for y in YEARS},*v) for k,v in specs.items()} for c in sorted({r["continent"] for r in rr["2025"]})}
    return out

def main():
    p=argparse.ArgumentParser()
    for n in ("agb","agb_model","disjoint","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();agb,old,previous,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.agb,a.agb_model,a.disjoint,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if previous["version"]!="railstatus-observer-disjoint-plant-bird-biomass-v1" or previous["biomass_source_sha256"]!=agb["source_sha256"] or any(previous["source_artifact_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")):raise ValueError("Prior disjoint artifact mismatch")
    prior_panels={p["land"]:p for p in previous["panels"] if p["width_m"]==200000}
    if set(prior_panels)!={.5,.8}:raise ValueError("Prior annual support missing")
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
        for w in (200000,):
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
                oldp=prior_panels[gate]["common_annual"]
                z={"land":gate,"width_m":w,**analyze(by,oldp)}
                panels.append(z);print(json.dumps({"plant_bird_within_cell_panel":z},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    out={"version":"railstatus-observer-disjoint-plant-bird-within-cell-v1","source_replay":checks,
         "prior_disjoint_SHA256":hashlib.sha256(a.disjoint.read_bytes()).hexdigest(),"panels":panels,
         "primary":"200km land80 same 594 cells, 2023-25 balanced panel, cell and year FE, flexible all-four-half effort, mean two opposite-side slopes with paired spatial-block covariance; block-year FE and adjacent periods sensitivity.",
         "limits":["static geography, AGB, population and rail status absorbed, not time-varying biomass effect",
                   "annual recorded richness changes can reflect sampling, access, detectability and weather rather than community change",
                   "opposite-side observers cannot share identities but same platform/access remain",
                   "three years only; temporal correlation within blocks retained, across blocks may remain",
                   "no causal trophic mediation, no bird biomass, no rail causation","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Cross-taxon landscape community similarity coupling on fixed railway pairs.
Plant/insect within-half similarities predict bird opposite-half similarity.
Observers cannot overlap between predictors and bird; plant/insect may overlap.
No occurrence rescan. Same-pair conditional contrasts, not causal mediation.
"""
import argparse,json,math,hashlib
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
import analyze_world_railstatus_observer_split_connectivity as s
from analyze_world_railstatus_plant_bird_disjoint import combine,interval
from analyze_world_railstatus_plant_bird_biomass import BIO
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import fit_predict,metrics,CONTINENTS,block_weights
GROUPS=("plants","birds","insects","fungi")
EXPECTED={(100000,.5):(2696,2131,2243,876),(100000,.8):(2331,1807,1915,712),
          (200000,.5):(1077,1041,950,407),(200000,.8):(869,837,753,308)}
ORIENT=((0,1),(1,0))
GEO=s.FIXED_PAIR+("mean_"+BIO,"diff_"+BIO)
def attach_agb(rows,ag20,ag24,roads):
    out={}
    for r in rows:
        k=r["key"]
        a,c=ag20.get(k),ag24.get(k)
        if a is None or c is None:continue
        if min(a["valid_area_m2"],c["valid_area_m2"])<.9*roads[k]["source_GRIP4_total_land_area_m2"]:continue
        if any(not math.isfinite(z["agb_Mg_ha"]) or z["agb_Mg_ha"]<0 for z in (a,c)):raise ValueError("Invalid AGB")
        out[k]={2020:math.log1p(a["agb_Mg_ha"])-math.log(51),2024:math.log1p(c["agb_Mg_ha"])-math.log(51)}
    return out

def join(groups,raw,agb,gs):
    maps={g:{(r["a"],r["b"]):r for r in groups[g]} for g in gs}
    keys=set.intersection(*(set(m) for m in maps.values()))
    out=[]
    for a,b in sorted(keys):
        if a not in agb or b not in agb:continue
        k=a,b;p=maps["plants"][k];r={f:p[f] for f in (*s.FIXED_PAIR,"a","b","continent","block")}
        for yr in (2020,2024):
            x,y=agb[a][yr],agb[b][yr]
            r[f"mean_AGB_{yr}"]=(x+y)/2;r[f"diff_AGB_{yr}"]=abs(x-y)
        for g in gs:
            q=maps[g][k]
            if any(q[f]!=p[f] for f in s.FIXED_PAIR):raise ValueError("Geographic pair design mismatch")
            for f in s.EFF_PAIR:r[g+"_"+f]=q[f]
            for h in (0,1):
                for metric in ("jaccard","simpson"):
                    v=q[f"{metric}_{h}{h}"]
                    if not 0<=v<=1:raise ValueError("Invalid same-side similarity")
                    r[f"{g}_{h}_{metric}"]=v
                ov=raw[g][k]["same_side_observer_overlap"][str(h)]["jaccard"]
                if ov is None or not 0<=ov<=1:raise ValueError("Invalid same-half observer overlap")
                r[f"{g}_{h}_observer_overlap"]=ov
        out.append(r)
    return out

def fields(gs,endpoint=False,flex=False):
    effort=tuple(g+"_"+f for g in gs for f in s.EFF_PAIR if not endpoint or f.startswith("diff_"))
    geo=tuple(f for f in GEO if not endpoint or not f.startswith("mean_"))
    ov=tuple(f"{g}_{i}_observer_overlap" for g in gs for i in (0,1))
    ff=tuple(f+"_sq" for f in effort if "log_species" not in f) if flex else ()
    return geo+effort+ov+ff

def setup(rows,gs,epoch,flex):
    out=[]
    for r in rows:
        z=dict(r,**{"mean_"+BIO:r[f"mean_AGB_{epoch}"],"diff_"+BIO:r[f"diff_AGB_{epoch}"]})
        for f in fields(gs,False,flex):
            if f.endswith("_sq") and f not in z:z[f]=z[f[:-3]]**2
        out.append(z)
    return out

def score_fit(rows,pred,target,terms,endpoint):
    blocks=defaultdict(list)
    for i,r in enumerate(rows):blocks[r["continent"],r["block"]].append(i)
    if len(rows)<100 or len(blocks)<20:return None
    w=block_weights(rows);X=np.asarray([[r[f] for f in pred]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite pair regression")
    orig=np.sqrt(np.sum(w[:,None]*X[:,:-1]**2,axis=0))
    if endpoint:
        nodes={k:i for i,k in enumerate(sorted({r[q] for r in rows for q in ("a","b")}))}
        owners={};parent=list(range(len(nodes)))
        def root(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for r in rows:
            a,b=r["a"],r["b"];block=r["continent"],r["block"]
            for k in (a,b):
                if k in owners and owners[k]!=block:raise ValueError("Endpoint across spatial blocks")
                owners[k]=block
            # Adjacent regular-grid graph is bipartite.
            if (sum(a)-sum(b))%2!=1:raise ValueError("Non-bipartite grid edge")
            parent[root(nodes[a])]=root(nodes[b])
        rank_fe=len(nodes)-len({root(i) for i in range(len(nodes))})
        ii=np.repeat(np.arange(len(rows)),2);jj=np.array([nodes[r[k]] for r in rows for k in ("a","b")])
        D=csr_matrix((np.repeat(np.sqrt(w),2),(ii,jj)),shape=(len(rows),len(nodes)))
        X*=np.sqrt(w)[:,None]
        for j in range(X.shape[1]):
            sol=lsmr(D,X[:,j],atol=1e-12,btol=1e-12,maxiter=4000)
            if sol[1] not in (0,1,2,4,5):raise ValueError("Endpoint projection failed")
            X[:,j]-=D@sol[0]
        if np.max(np.abs(D.T@X))>1e-7:raise ValueError("Endpoint residual check")
    else:
        rank_fe=len(blocks)
        for ix in blocks.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
        X*=np.sqrt(w)[:,None]
    A=X[:,:-1];y=X[:,-1];norm=np.sqrt(np.sum(A*A,axis=0));keep=norm>np.maximum(1e-9,orig*1e-8)
    names=[f for f,k in zip(pred,keep) if k]
    if not set(terms)<=set(names) or len(rows)-rank_fe-len(names)<30:return None
    Z=A[:,keep]/norm[keep];sv=np.linalg.svd(Z,compute_uv=False)
    if sv[-1]<1e-8*sv[0]:return None
    inv=np.linalg.inv(Z.T@Z);beta=inv@(Z.T@y);e=y-Z@beta
    u=np.asarray([inv@(Z[ix].T@e[ix])/norm[keep] for ix in blocks.values()])
    effects={}
    for t in terms:
        j=names.index(t)
        if sum(np.linalg.norm(Z[ix,j])>1e-8 for ix in blocks.values())<20:return None
        effects[t]=(float(beta[j]/norm[keep][j]),{k:float(v) for k,v in zip(blocks,u[:,j])})
    return effects

def model(rows,gs,epoch=2024,metric="jaccard",endpoint=False,flex=False,do_cv=False):
    rr=setup(rows,gs,epoch,flex);base=fields(gs,endpoint,flex)
    qa={"n_pairs":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
        "n_endpoints":len({r[k] for r in rr for k in ("a","b")})}
    fits={}
    others=tuple(g for g in gs if g not in ("plants","birds"))
    for i,j in ORIENT:
        plant=f"plants_{i}_{metric}";target=f"birds_{j}_{metric}"
        oth=tuple(f"{g}_{i}_{metric}" for g in others)
        a=score_fit(rr,base+(plant,),target,(plant,),endpoint)
        c=score_fit(rr,base+(plant,)+oth,target,(plant,)+oth,endpoint)
        if a is None or c is None:return dict(qa,status="insufficient support or rank")
        z={"plant_only":a[plant],"plant_joint":c[plant]}
        z.update({g:c[k] for g,k in zip(others,oth)})
        z["plant_attenuation"]=combine(z,{"plant_only":1.,"plant_joint":-1.})
        if len(others)==2:z["insects_minus_fungi"]=combine(z,{"insects":1.,"fungi":-1.})
        fits[i,j]=z
    effects={k:interval(combine({ij:z[k] for ij,z in fits.items()},{ij:.5 for ij in ORIENT})) for k in fits[ORIENT[0]]}
    if abs(effects["plant_only"]["estimate"]-effects["plant_joint"]["estimate"]-effects["plant_attenuation"]["estimate"])>1e-9:raise ValueError("Paired effect identity failed")
    out=dict(qa,status="estimated",effects=effects)
    if do_cv and not endpoint:out["continent_holdout"]=cv(rr,base,others,metric)
    return out

def cv(rows,base,others,metric):
    folds=[]
    for c in CONTINENTS:
        te=[r for r in rows if r["continent"]==c];tr=[r for r in rows if r["continent"]!=c]
        if len(te)<50 or len(tr)<max(300,4*(len(base)+len(others)+1)):continue
        v=[]
        for i,j in ORIENT:
            target=f"birds_{j}_{metric}";full=base+tuple(f"{g}_{i}_{metric}" for g in ("plants",*others))
            mse=metrics(te,fit_predict(tr,te,full,target,True),target)[0]
            v.append({g:metrics(te,fit_predict(tr,te,tuple(k for k in full if k!=f"{g}_{i}_{metric}"),target,True),target)[0]-mse for g in ("plants",*others)})
        folds.append({"continent":c,"n_test":len(te),"n_train":len(tr),"dMSE":{g:sum(q[g] for q in v)/2 for g in v[0]}})
    return {"n_folds":len(folds),"adequate_folds":len(folds)>=4,"folds":folds,
       "mean_dMSE":{g:sum(q["dMSE"][g] for q in folds)/len(folds) for g in ("plants",*others)} if folds else {},
       "positive_folds":{g:sum(q["dMSE"][g]>0 for q in folds) for g in ("plants",*others)}}

def panel(rows,gs):
    return {"n_pairs":len(rows),"primary":model(rows,gs,do_cv=True),
      "flexible_effort":model(rows,gs,flex=True,do_cv=True),"simpson":model(rows,gs,metric="simpson"),
      "AGB2020":model(rows,gs,epoch=2020),"endpoint_FE":model(rows,gs,endpoint=True),
      "leave_continent_out":{c:model([r for r in rows if r["continent"]!=c],gs) for c in sorted({r["continent"] for r in rows})}}

def main():
    p=argparse.ArgumentParser()
    for k in ("split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    src,prior,edge,agb,am,plant,rail,status,tree,road,pop=map(s.read,(a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_gate_observations_observers_species"]!=[25,5,10] or src["split_seed"]!="tetsu-railstatus-observer-split-v1":raise ValueError("Wrong frozen split")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(src["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Frozen source mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):raise ValueError("Wrong geometry")
    if tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong road/tree")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong climate/map/pop")
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    panels=[]
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
                    full=s.build_pairs([s.adapt_pair(p) for p in sp[g]],rr,el)
                    groups[g]=[z for r in full if (z:=s.split_pair(r,raw[g][r["a"],r["b"]])) is not None]
                    for k,v in attach_agb(rr,ag20,ag24,rd).items():
                        if k in ag and ag[k]!=v:raise ValueError("AGB grid mismatch")
                        ag[k]=v
                if tuple(len(groups[g]) for g in GROUPS)!=EXPECTED[w,gate]:raise ValueError("Exact previous pair counts not replayed")
                z={"width_m":w,"land":gate,"source_pair_counts":{g:len(groups[g]) for g in GROUPS}}
                for label,gs in (("plants_birds",("plants","birds")),("plants_insects_birds",("plants","birds","insects")),("all_four",GROUPS)):
                    z[label]=panel(join(groups,raw,ag,gs),gs)
                panels.append(z);print(json.dumps({"beta_coupling_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-disjoint-community-composition-coupling-v1","source_replay":checks,"biomass_source_SHA256":agb["source_sha256"],"source_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","edges","agb","agb_model")},"panels":panels,
         "primary":"100km land80 2025 exact frozen rail-positive adjacent same-block pair geography; plant/insect same-half Jaccard, bird opposite-half Jaccard. Conditional insect and plant slopes, paired attenuation; all-taxa symmetric sampling and richness terms, each within-half observer overlap, biomass mean/difference, climate, population, road, mapped rail/crossing controls, block FE.",
         "limits":["same observers may occur in two endpoints within a taxon; measured overlap adjusted","plant and insect predictors may share observers; bird opposite side cannot share identity","whole landscape occurrence lists, no realized food interactions or causal mediation","Jaccard includes richness/nestedness; Simpson overlap sensitivity","endpoint FE may lack cycles; no threshold relaxation","conditional prediction not independent ecosystem validation","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

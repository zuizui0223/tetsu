#!/usr/bin/env python3
"""Selection-weighted sensitivity of cross-taxon recorded-richness coupling.
No imputed biodiversity in unsampled cells; no causal or population identification."""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from analyze_world_railstatus_coupling_climate import climate_extra, EXTRA as CLIMATE_FIELDS, fit_model as climate_model
import analyze_world_railstatus_observation_stages as stage
from analyze_world_railstatus_observation_support import (
    b,GEO,extract,linear,interval,combine,attach_agb,MIN_FULL,MIN_SPLIT,SALT,CATS,LOG,GROUPS,PB,make_frame)
from analyze_world_railstatus_plant_bird_disjoint import EFF,FLEX,attach,one
from analyze_world_railstatus_plant_bird_coupling_modifiers import replay_check
from analyze_world_railstatus_plant_insect_bird import extra_attach,fields_for,model
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import CONTINENTS
N_FOLDS=5
LAMBDA=.01
P_MIN=.1
CROSS=((0,1),(1,0))
PROP_FIELDS=GEO+(b.BIO,)
EXPECTED_PB={(100000,.5):1911,(100000,.8):1640,(200000,.5):966,(200000,.8):787}
EXPECTED_TRI={(100000,.5):1716,(100000,.8):1485,(200000,.5):849,(200000,.8):690}

def group(r):return r["continent"],r["block"]
def fold(r):
    z=json.dumps(group(r),separators=(",",":")).encode()
    return int.from_bytes(hashlib.sha256(b"tetsu-selection-fold-v1:"+z).digest()[:8],"big")%N_FOLDS

def design(rows):
    return np.array([[r[k] for k in PROP_FIELDS]+[float(r["continent"]==c) for c in CONTINENTS[1:]] for r in rows],float)

def propensity(frame,selected):
    if len({r["key"] for r in frame})!=len(frame) or not selected<={r["key"] for r in frame}:raise ValueError("Invalid selection frame")
    y=np.array([float(r["key"] in selected) for r in frame])
    X=design(frame);w=b.block_weights(frame)
    if not np.isfinite(X).all():raise ValueError("Nonfinite frame design")
    p=np.full(len(frame),np.nan);diag=[]
    folds=np.array([fold(r) for r in frame])
    for f in range(N_FOLDS):
        it=np.flatnonzero(folds!=f);iv=np.flatnonzero(folds==f)
        if not len(iv) or len(it)<100 or min(y[it].sum(),len(it)-y[it].sum())<20:raise ValueError("Insufficient propensity training support")
        wt=w[it]/w[it].sum();mu=np.average(X[it],axis=0,weights=wt)
        sd=np.sqrt(np.average((X[it]-mu)**2,axis=0,weights=wt));sd[sd<1e-10]=1.
        A=np.column_stack((np.ones(len(it)),(X[it]-mu)/sd))
        V=np.column_stack((np.ones(len(iv)),(X[iv]-mu)/sd))
        yy=y[it]
        def obj(q):
            z=A@q
            loss=float(wt@(np.logaddexp(0,z)-yy*z)+LAMBDA/2*(q[1:]@q[1:]))
            grad=A.T@(wt*(expit(z)-yy));grad[1:]+=LAMBDA*q[1:]
            return loss,grad
        opt=minimize(obj,np.zeros(A.shape[1]),jac=True,method="L-BFGS-B",options={"maxiter":1000,"ftol":1e-12,"gtol":1e-8})
        if not opt.success or not np.isfinite(opt.x).all():raise ValueError("Propensity optimization failed: "+str(opt.message))
        p[iv]=expit(V@opt.x)
        diag.append({"fold":f,"n_test":len(iv),"n_selected_test":int(y[iv].sum()),"n_blocks_test":len({group(frame[i]) for i in iv})})
    if not np.isfinite(p).all() or np.any((p<=0)|(p>=1)):raise ValueError("Invalid cross-fitted probabilities")
    bins=[]
    for lo,hi in zip((0,.1,.2,.4,.6,.8),(.1,.2,.4,.6,.8,1.)):
        ix=np.flatnonzero((p>=lo)&(p<hi))
        bins.append({"lo":lo,"hi":hi,"n":len(ix),"selected":int(y[ix].sum()),
            "predicted":float(np.average(p[ix],weights=w[ix])) if len(ix) else None,
            "observed":float(np.average(y[ix],weights=w[ix])) if len(ix) else None})
    out={"n_frame":len(frame),"n_selected":int(y.sum()),"weighted_brier":float(np.average((y-p)**2,weights=w)),
        "constant_brier":float(np.average((y-np.average(y,weights=w))**2,weights=w)),
        "calibration":bins,"folds":diag,"frame_fraction_p_ge_0_1":float(w[p>=P_MIN].sum()/w.sum()),
        "selected_fraction_p_ge_0_1":float(np.mean(p[y==1]>=P_MIN)) if y.sum() else None,
        "selection_probability_quantiles_selected":[float(z) for z in np.quantile(p[y==1],[0,.1,.5,.9,1])] if y.sum() else [],
        "by_class":{c:{"n":int(sum(r["rail_class"]==c for r in frame)),
            "n_selected":int(sum(r["rail_class"]==c and r["key"] in selected for r in frame)),
            "n_overlap":int(sum(r["rail_class"]==c and p[i]>=P_MIN for i,r in enumerate(frame)))}
            for c in ("active_positive","all_four_zero","other_only")}}
    return {r["key"]:float(v) for r,v in zip(frame,p)},out

def weights(rows,mode,frame_n,p):
    if mode=="native":w=np.array(b.block_weights(rows),float)
    elif mode=="frame":w=np.array([1/frame_n[group(r)] for r in rows],float)
    elif mode=="ipw":
        w=np.array([1/(frame_n[group(r)]*p[r["key"]]) for r in rows],float)
    else:raise ValueError("Unknown weight mode")
    if not np.isfinite(w).all() or np.any(w<=0):raise ValueError("Invalid analysis weights")
    return w/w.sum()

def weighted_fit(rows,fields,target,term,w):
    groups=defaultdict(list)
    for i,r in enumerate(rows):groups[group(r)].append(i)
    if len(rows)<100 or len(groups)<20 or 1/np.sum(w*w)<100 or len(rows)-len(groups)-len(fields)<30:return None
    X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite association")
    for ix in groups.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    X*=np.sqrt(w)[:,None];A=X[:,:-1];y=X[:,-1]
    scale=np.sqrt(np.sum(A*A,axis=0));keep=scale>1e-10;names=[k for k,v in zip(fields,keep) if v]
    if term not in names:return None
    Z=A[:,keep]/scale[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return None
    inv=np.linalg.inv(Z.T@Z);coef=inv@(Z.T@y);e=y-Z@coef;j=names.index(term)
    return float(coef[j]/scale[keep][j]),{k:float((inv@(Z[ix].T@e[ix]))[j]/scale[keep][j]) for k,ix in groups.items()}

def association(rows,groups,mode,frame_n,p,flex=True,extra=()):
    rr=b.with_epoch(b.local_support(rows),2024)
    f=fields_for(groups,flex)+extra
    if len(set(f))!=len(f):raise ValueError("Duplicate association predictors")
    qa={"n_cells":len(rr),"n_blocks":len({group(r) for r in rr}),"continents":dict(Counter(r["continent"] for r in rr))}
    if not rr:return dict(qa,status="insufficient support"),None
    w=weights(rr,mode,frame_n,p);qa.update(ESS=float(1/np.sum(w*w)),max_weight=float(w.max()))
    fits={}
    for i,j in CROSS:
        plant=f"plants_{i}_species";target=f"birds_{j}_species"
        z={};base=f+(plant,);full=base+tuple(f"{g}_{i}_species" for g in groups)
        for name,fields,term in (("plant",full,plant),("plant_without_insects",base,plant)):
            z[name]=weighted_fit(rr,fields,target,term,w)
        for g in groups:z[g]=weighted_fit(rr,full,target,f"{g}_{i}_species",w)
        if any(v is None for v in z.values()):return dict(qa,status="insufficient support or rank"),None
        z["attenuation"]=combine(z,{"plant_without_insects":1.,"plant":-1.})
        fits[i,j]=z
    effects={k:combine({ij:v[k] for ij,v in fits.items()},{ij:.5 for ij in CROSS}) for k in fits[CROSS[0]]}
    return dict(qa,status="estimated",effects={k:interval(v) for k,v in effects.items()}),effects

def diagnostic(frame,rows,groups):
    p,qa=propensity(frame,{r["key"] for r in rows})
    frame_n=Counter(group(r) for r in frame)
    subset=b.local_support([r for r in rows if p[r["key"]]>=P_MIN])
    out={"propensity":qa,"native_all":association(rows,groups,"native",frame_n,p)[0],
         "overlap_class_counts":dict(Counter(next(x["rail_class"] for x in frame if x["key"]==r["key"]) for r in subset))}
    scores={}
    for m in ("native","frame","ipw"):
        out[m],scores[m]=association(subset,groups,m,frame_n,p)
    if all(scores[m] is not None for m in ("native","ipw")):
        out["ipw_minus_native"]={k:interval(combine({m:scores[m][k] for m in ("ipw","native")},{"ipw":1.,"native":-1.})) for k in scores["native"]}
    return out

def check_effect(x,y):
    if abs(x["estimate"]-y["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(x["interval95_approx"],y["interval95_approx"]))>1e-9:raise ValueError("Unweighted effect replay failed")

def match_rows(frame,grid,sp,rail,clim,built,status,w,rd,tr,pp,gate,ag20,ag24):
    geo={};maps={g:extract(sp[g]) for g in ("plants","birds","insects")}
    full={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in maps}
    for g in maps:
        if any(full[g].get(k)!=v[0] for k,v in maps[g].items()):raise ValueError("Exact split counts changed")
    for g in PB:
        rr,_,_=b.original_cell_rows([b.convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
        rr,_=b.inject_built(rr,built[str(w)]);rr,_=b.inject_status(rr,status,w)
        rr,_=b.enrich_cells(rr,rd,tr,w,gate);rr,_=b.join_population(rr,pp,rd)
        geo[g]=[b.flexible(r) for r in rr]
    joined,_=b.join(geo,ag20,ag24,w,rd)
    bc={r["key"]:r["species"] for r in geo["birds"]}
    pb=b.local_support(attach([dict(r,bird_species=bc[r["key"]]) for r in joined],maps))
    tri=b.local_support(extra_attach(pb,maps,full,("insects",)))
    return pb,tri

def check_native(new,old):
    if new["status"]!=old["status"] or new["n_cells"]!=old["n_cells"]:raise ValueError("Prior weighted support changed")
    if new["status"]=="estimated":
        if abs(new["ESS"]-old["ESS"])>1e-7:raise ValueError("Prior ESS changed")
        for k,v in new["effects"].items():check_effect(v,old["effects"][k])

def joint(frame,rows,groups,ext,prior,climate_prior):
    p,qa=propensity(frame,{r["key"] for r in rows})
    frame_n=Counter(group(r) for r in frame)
    if qa["n_frame"]!=prior["propensity"]["n_frame"] or qa["n_selected"]!=prior["propensity"]["n_selected"] or abs(qa["weighted_brier"]-prior["propensity"]["weighted_brier"])>1e-9:
        raise ValueError("Prior selection model replay failed")
    overlap=b.local_support([r for r in rows if p[r["key"]]>=P_MIN])
    for m in ("native","frame","ipw"):
        old,_=association(overlap,groups,m,frame_n,p);check_native(old,prior[m])
    allcl=b.local_support([dict(r,**ext[r["key"]]) for r in rows if r["key"] in ext])
    q,_=climate_model(allcl,groups,True,CLIMATE_FIELDS)
    z=climate_prior["models"]["flexible_effort"]["extra_climate"]
    if q["n_cells"]!=z["n_cells"] or q["status"]!=z["status"]:raise ValueError("Prior climate support changed")
    if q["status"]=="estimated":
        for k,v in q["effects"].items():check_effect(v,z["effects"][k])
    common=b.local_support([dict(r,**ext[r["key"]]) for r in overlap if r["key"] in ext])
    out={"n_original_overlap":len(overlap),"n_climate_overlap":len(common),
         "weighted_frame_fraction_overlap":qa["frame_fraction_p_ge_0_1"],
         "zero_rail_cells":sum(all(r[c+"_km"]==0 for c in CATS) for r in common)}
    scores={}
    for m in ("native","ipw"):
        for label,fields in (("baseline",()),("climate",CLIMATE_FIELDS)):
            key=m+"_"+label
            out[key],scores[key]=association(common,groups,m,frame_n,p,extra=fields)
    if all(v is not None for v in scores.values()):
        out["climate_change_under_ipw"]={k:interval(combine({"extended":scores["ipw_climate"][k],"base":scores["ipw_baseline"][k]},{"extended":1.,"base":-1.})) for k in scores["ipw_climate"]}
        out["weight_change_under_climate"]={k:interval(combine({"ipw":scores["ipw_climate"][k],"native":scores["native_climate"][k]},{"ipw":1.,"native":-1.})) for k in scores["ipw_climate"]}
    return out

# Observation-only threshold does not contain richness or alpha response.
# Fixed extra climate and squares of non-status environmental terms that lack a square.
SQUARE_BASE=tuple(k for k in GEO if k not in tuple(LOG.values()) and not k.endswith("_sq") and k+"_sq" not in GEO)
def effort_fit(rows,fields):
    rr=b.local_support(rows)
    f={k:linear(rr,fields,"both_effort",k) for k in (b.ACTIVE,b.INACTIVE)}
    out={"n_cells":len(rr),"n_blocks":len({group(r) for r in rr})}
    if any(v is None for v in f.values()):return dict(out,status="insufficient support or rank"),None
    z={"active":f[b.ACTIVE],"inactive":f[b.INACTIVE],"active_minus_inactive":combine(f,{b.ACTIVE:1.,b.INACTIVE:-1.})}
    return dict(out,status="estimated",effects={k:interval(v) for k,v in z.items()}),z

def effort_sensitivity(frame,ext,stage_prior):
    old,_=effort_fit(frame,GEO)
    p=stage_prior["conditional"]["effects"]["both_effort"]
    if old["status"]!="estimated" or old["n_cells"]!=stage_prior["conditional"]["n_cells"]:raise ValueError("Prior effort support changed")
    for k,v in old["effects"].items():check_effect(v,p[k])
    rr=b.local_support([dict(r,**ext[r["key"]]) for r in frame if r["key"] in ext])
    for r in rr:
        for k in SQUARE_BASE:r["extra_sq_"+k]=r[k]**2
    specs={"baseline":GEO,"climate":GEO+CLIMATE_FIELDS,
           "flexible_environment":GEO+CLIMATE_FIELDS+tuple("extra_sq_"+k for k in SQUARE_BASE)}
    out={"n_frame":len(frame),"n_climate_common":len(rr),"models":{}}
    scores={}
    for lab,fields in specs.items():
        out["models"][lab],scores[lab]=effort_fit(rr,fields)
    if scores["baseline"] is not None and scores["flexible_environment"] is not None:
        out["flexible_minus_baseline"]={k:interval(combine({"new":scores["flexible_environment"][k],"old":scores["baseline"][k]},{"new":1.,"old":-1.})) for k in scores["baseline"]}
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("weighted","climate_model","support","stages","disjoint","insect_model","annual","split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    src,prior,edge,agb,am,plant,rail,status,tree,road,pop=map(b.read,(a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    previous=b.read(a.support)
    if previous["version"]!="railstatus-observation-support-audit-v1" or any(previous["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")):raise ValueError("Prior support audit provenance mismatch")
    replay={(p["width_m"],p["land"]):p for p in previous["panels"]}
    st=b.read(a.stages);dis=b.read(a.disjoint);ins=b.read(a.insect_model)
    if st["version"]!="railstatus-observation-stage-decomposition-v1" or st["support_model_SHA256"]!=hashlib.sha256(a.support.read_bytes()).hexdigest():raise ValueError("Wrong prior stage model")
    if dis["version"]!="railstatus-observer-disjoint-plant-bird-biomass-v1" or dis["biomass_source_sha256"]!=agb["source_sha256"] or any(dis["source_artifact_sha256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("agb","agb_model","split","annual")):raise ValueError("Prior disjoint source mismatch")
    if ins["version"]!="railstatus-observer-disjoint-plant-insect-bird-biomass-v1" or ins["prior_disjoint_SHA256"]!=hashlib.sha256(a.disjoint.read_bytes()).hexdigest() or ins["biomass_source_sha256"]!=agb["source_sha256"]:raise ValueError("Prior insect source mismatch")
    we=b.read(a.weighted);ce=b.read(a.climate_model)
    if we["version"]!="railstatus-selection-weighted-coupling-v1" or we["stage_model_SHA256"]!=hashlib.sha256(a.stages.read_bytes()).hexdigest() or any(we["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")):raise ValueError("Wrong weighted result source")
    if ce["version"]!="railstatus-observer-disjoint-richness-extra-climate-v1" or ce["triad_source_sha256"]!=hashlib.sha256(a.insect_model.read_bytes()).hexdigest() or ce["worldclim_sha256"]!=b.sha_file(a.climate):raise ValueError("Wrong climate result source")
    we_panels={(p["width_m"],p["land"]):p for p in we["panels"]}
    ce_panels={(p["width_m"],p["land"]):p for p in ce["panels"]}
    stage_panels={(p["width_m"],p["land"]):p for p in st["panels"]}
    dis_panels={(p["width_m"],p["land"]):p for p in dis["panels"]}
    ins_panels={(p["width_m"],p["land"]):p for p in ins["panels"]}
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Frozen source gate mismatch")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(src["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Source replay mismatch")
    if not status["version"].startswith("railstatus-v1") or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Rail frame mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass mismatch")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen environmental mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in b.GRID_WIDTHS:
            grid=prior["grids"][str(w)];sp=src["pairs"][str(w)]
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            # Candidate frame is independent of recorded species. Restrict to genuine road
            # source cell keys, then land/climate and common environmental support.
            clim,cqa=b.source_climate_lookup(land,ras,scale,set(rd),w)
            ext,cmeta=climate_extra(a.climate,ras,clim,w)
            base=make_frame(clim,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)])
            base,_=b.inject_built(base,built[str(w)]);base,_=b.inject_status(base,status,w)
            counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in PB}
            splits={g:extract(sp[g]) for g in PB}
            for g in PB:
                if any(counts[g].get(k)!=v[0] for k,v in splits[g].items()):raise ValueError("Exact split count mismatch")
            for gate in b.LAND_THRESHOLDS:
                rr,lqa=b.enrich_cells(base,rd,tr,w,gate);rr,pqa=b.join_population(rr,pp,rd)
                # Main denominator independent of AGB; repeat on AGB-supported cells.
                rows=[stage.outcomes(r,counts,splits) for r in rr]
                ag=attach_agb(rr,ag20,ag24,rd)
                sub=[r for r in rows if r["key"] in ag]
                z={"width_m":w,"land":gate,"environmental_frame_qa":{"road_keys":len(rd),"climate":cqa,"land":lqa,"population":pqa},
                   "funnel":stage.funnel(rows),"conditional":stage.fit(rows),
                   "unadjusted_local":stage.fit(rows,False),"AGB_supported_funnel":stage.funnel(sub)}
                stage.verify_replay(z,replay[w,gate])
                if z["funnel"]["all"]["counts"]!=stage_panels[w,gate]["funnel"]["all"]["counts"]:raise ValueError("Stage counts changed")
                # Recreate independently supported selected rows from exact previous alpha model.
                # climate lookup of biological records equals complete road frame on final matched rows.
                keys={b.key(x) for g in PB for x in grid[g]}
                bclim,_=b.source_climate_lookup(land,ras,scale,keys,w)
                pb,tri=match_rows(sub,grid,sp,rail,bclim,built,status,w,rd,tr,pp,gate,ag20,ag24)
                if len(pb)!=EXPECTED_PB[w,gate] or len(tri)!=EXPECTED_TRI[w,gate]:raise ValueError("Previous selected support changed")
                oldpb=dis_panels[w,gate]["2025"]["flexible_effort"]
                replay_check(one(pb,flex=True)[0],oldpb)
                oldtri=ins_panels[w,gate]["insect"]["flexible_effort"]
                q=model(tri,("insects",),flex=True)
                if q["status"]!=oldtri["status"] or q["n_cells"]!=oldtri["n_cells"]:raise ValueError("Previous insect support changed")
                for k,v in q["effects"].items():check_effect(v,oldtri["effects"][k])
                frame=[dict(r,**{b.BIO:ag[r["key"]][2024]}) for r in sub]
                pbd=joint(frame,pb,(),ext,we_panels[w,gate]["plant_bird"],ce_panels[w,gate]["plant_bird"])
                trid=joint(frame,tri,("insects",),ext,we_panels[w,gate]["plant_insect_bird"],ce_panels[w,gate]["insect"])
                z={"width_m":w,"land":gate,"plant_bird":pbd,"plant_insect_bird":trid,"effort":effort_sensitivity(rows,ext,stage_panels[w,gate]),"climate_source":cmeta}
                panels.append(z);print(json.dumps({"joint_climate_selection_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-joint-climate-selection-v1","weighted_source_SHA256":hashlib.sha256(a.weighted.read_bytes()).hexdigest(),"climate_source_SHA256":hashlib.sha256(a.climate_model.read_bytes()).hexdigest(),"stage_model_SHA256":hashlib.sha256(a.stages.read_bytes()).hexdigest(),"source_replay":checks,"source_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")},"panels":panels,
        "primary":"100km land80 fixed p>=0.1 overlap; same prior cross-fitted selection weights plus BIO4/BIO12/BIO15 and squares in outcome model; separate full-frame effort-only active-minus-inactive contrast with climate and nonlinear environmental controls.",
        "limits":["no record is not species absence","selection depends on detection and true richness; missing-not-at-random not solved","intervals conditional on fixed prior propensity weights; no identification outside overlap or of true biodiversity; same data sensitivities not independent replications","environmental complete-case frame is not all land","map status is not measured traffic, no causal observation bias or biodiversity effect","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

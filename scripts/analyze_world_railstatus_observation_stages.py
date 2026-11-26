#!/usr/bin/env python3
"""Decompose observation support using unchanged global environmental frame.
No new occurrence scan; exact additive stage-loss identities on common rows."""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
import analyze_world_railstatus_observation_support as old
from analyze_world_railstatus_observation_support import (
    b,GEO,extract,linear,interval,combine,attach_agb,MIN_FULL,MIN_SPLIT,SALT,CATS,LOG,GROUPS,PB,make_frame)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import fit_predict,metrics,CONTINENTS

STAGES=("both_record","both_effort","both_full","both_split_endpoint")
LOSSES=("effort_loss","richness_loss","post_full_loss")
OUTCOMES=old.OUTCOMES+("plant_effort","bird_effort","both_effort")+LOSSES
STATUS=tuple(LOG[c] for c in CATS)
FIELDS=GEO

def outcomes(row,counts,splits):
    q=old.outcomes(row,counts,splits)
    for g in PB:
        n=counts[g].get(row["key"],(0,0,0))
        q[("plant" if g=="plants" else "bird")+"_effort"]=float(n[0]>=MIN_FULL[0] and n[1]>=MIN_FULL[1])
    q["both_effort"]=q["plant_effort"]*q["bird_effort"]
    for x,y,k in zip(STAGES,STAGES[1:],LOSSES):
        q[k]=q[x]-q[y]
        if q[k] not in (0.,1.):raise ValueError("Non-nested eligibility stages")
    return q

def rates(rows):
    n=len(rows)
    counts={k:int(sum(r[k] for r in rows)) for k in OUTCOMES}
    return {"n":n,"counts":counts,"fractions":{k:v/n for k,v in counts.items()} if n else {},
        "conditional_pass":{y:(counts[y]/counts[x] if counts[x] else None) for x,y in zip(STAGES,STAGES[1:])}}

def funnel(rows):
    return {"all":rates(rows),"by_class":{c:rates([r for r in rows if r["rail_class"]==c]) for c in ("active_positive","all_four_zero","other_only")},
            "by_continent":{c:{cl:rates([r for r in rows if r["continent"]==c and r["rail_class"]==cl]) for cl in ("active_positive","all_four_zero","other_only")} for c in sorted({r["continent"] for r in rows})}}

def fit(rows,adjusted=True):
    rr=b.local_support(rows); fields=FIELDS if adjusted else STATUS
    fits={};effects={}
    for target in OUTCOMES:
        f={t:linear(rr,fields,target,t) for t in (b.ACTIVE,b.INACTIVE)}
        if any(v is None for v in f.values()):
            effects[target]={"status":"insufficient support or rank"};continue
        f["difference"]=combine(f,{b.ACTIVE:1.,b.INACTIVE:-1.})
        fits[target]=f
        effects[target]={"status":"estimated","active":interval(f[b.ACTIVE]),"inactive":interval(f[b.INACTIVE]),"active_minus_inactive":interval(f["difference"])}
    for x,y,k in zip(STAGES,STAGES[1:],LOSSES):
        if not all(t in fits for t in (x,y,k)):continue
        for t in (b.ACTIVE,b.INACTIVE,"difference"):
            v=combine({a:fits[a][t] for a in (x,y,k)},{x:1.,y:-1.,k:-1.})
            if abs(v[0])>1e-9 or max(abs(s) for s in v[1].values())>1e-9:raise ValueError("Stage coefficient or influence identity failed")
    return {"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),"effects":effects}

def predict(rows):
    base=tuple(k for k in FIELDS if k not in STATUS)
    model_fields={"environment":base,"all_status":FIELDS,
        "without_active":tuple(k for k in FIELDS if k!=b.ACTIVE),
        "without_inactive":tuple(k for k in FIELDS if k!=b.INACTIVE)}
    out={}
    for t in STAGES:
        folds=[]
        for c in CONTINENTS:
            train=[r for r in rows if r["continent"]!=c];test=[r for r in rows if r["continent"]==c]
            if len(test)<100 or len(train)<max(300,4*len(FIELDS)):continue
            m={k:metrics(test,fit_predict(train,test,p,t,True),t)[0] for k,p in model_fields.items()}
            folds.append({"continent":c,"n_test":len(test),"n_train":len(train),"brier":m,
                "all_status_gain":m["environment"]-m["all_status"],
                "active_increment":m["without_active"]-m["all_status"],
                "inactive_increment":m["without_inactive"]-m["all_status"]})
        terms=("all_status_gain","active_increment","inactive_increment")
        out[t]={"n_folds":len(folds),"adequate_folds":len(folds)>=4,"folds":folds,
                "mean_gains":{k:sum(z[k] for z in folds)/len(folds) for k in terms} if folds else {},
                "positive_folds":{k:sum(z[k]>0 for z in folds) for k in terms}}
    return out

def verify_replay(new,prior):
    for lab in ("all","by_class"):
        v=new["funnel"][lab];p=prior["funnel"][lab]
        pairs=[(v,p)] if lab=="all" else [(v[k],p[k]) for k in p]
        for x,y in pairs:
            if x["n"]!=y["n"] or any(x["counts"][k]!=y["counts"][k] for k in old.OUTCOMES):raise ValueError("Prior counts changed")
    for lab in ("conditional","unadjusted_local"):
        x=new[lab];y=prior[lab]
        if (x["n_cells"],x["n_blocks"])!=(y["n_cells"],y["n_blocks"]):raise ValueError("Prior support changed")
        for k in old.OUTCOMES:
            a=x["effects"][k];d=y["effects"][k]
            if a["status"]!=d["status"]:raise ValueError("Prior status changed")
            if a["status"]=="estimated":
                for t in ("active","inactive","active_minus_inactive"):
                    aa,dd=a[t],d[t]
                    if abs(aa["estimate"]-dd["estimate"])>1e-9 or max(abs(i-j) for i,j in zip(aa["interval95_approx"],dd["interval95_approx"]))>1e-9:raise ValueError("Prior coefficient/interval changed")

def main():
    p=argparse.ArgumentParser()
    for k in ("support","split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    src,prior,edge,agb,am,plant,rail,status,tree,road,pop=map(b.read,(a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    previous=b.read(a.support)
    if previous["version"]!="railstatus-observation-support-audit-v1" or any(previous["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")):raise ValueError("Prior support audit provenance mismatch")
    replay={(p["width_m"],p["land"]):p for p in previous["panels"]}
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
            base=make_frame(clim,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)])
            base,_=b.inject_built(base,built[str(w)]);base,_=b.inject_status(base,status,w)
            counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in PB}
            splits={g:extract(sp[g]) for g in PB}
            for g in PB:
                if any(counts[g].get(k)!=v[0] for k,v in splits[g].items()):raise ValueError("Exact split count mismatch")
            for gate in b.LAND_THRESHOLDS:
                rr,lqa=b.enrich_cells(base,rd,tr,w,gate);rr,pqa=b.join_population(rr,pp,rd)
                # Main denominator independent of AGB; repeat on AGB-supported cells.
                rows=[outcomes(r,counts,splits) for r in rr]
                ag=attach_agb(rr,ag20,ag24,rd)
                sub=[r for r in rows if r["key"] in ag]
                z={"width_m":w,"land":gate,"environmental_frame_qa":{"road_keys":len(rd),"climate":cqa,"land":lqa,"population":pqa},
                   "funnel":funnel(rows),"conditional":fit(rows),
                   "unadjusted_local":fit(rows,False),"AGB_supported_funnel":funnel(sub)}
                verify_replay(z,replay[w,gate])
                z["heldout"]=predict(rows)
                panels.append(z);print(json.dumps({"observation_stage_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-observation-stage-decomposition-v1","support_model_SHA256":hashlib.sha256(a.support.read_bytes()).hexdigest(),"source_replay":checks,"source_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")},"panels":panels,
        "primary":"100km land80, both taxa >=50 eligible observations and >=5 observers regardless of richness; active-minus-inactive local block coefficient and exact cumulative eligibility-loss decomposition. Fixed continent holdout as predictive diagnostic.",
        "limits":["no record is not species absence","observation-only gate still depends on detection and participation; split endpoints depend jointly on richness, observers and adjacency; no causal decomposition","environmental complete-case frame is not all land","map status is not measured traffic, no causal observation bias or biodiversity effect","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Outcome-blind overlap diagnostic for existing plant-height status result.

No new occurrence scan. Fixed extreme status shares, exact spatial block,
explicit length ratio and standardized covariate calipers, maximum-cardinality
one-to-one assignment. Estimate only if support and balance both pass.
"""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from scipy.optimize import linear_sum_assignment
from analyze_world_railstatus_verified_plant_size import (
    b,GEO,CORE,CLIMATE,EXTRA,SHARE,TOTAL,linear,combine,interval,extract,climate_extra,
    GROUPS,WIDTHS,TRAITS,LIFE_SHA,LEP_SHA,SIZE_SHA,SOURCE_SHA,attach,estimate)
from analyze_world_railstatus_plant_size_temporal import check_effects
HEIGHT="height_nonwoody"
MATCH=tuple(k for k in GEO if k not in (b.ACTIVE,b.INACTIVE) and not k.endswith("_sq"))+tuple(k for k in CLIMATE if not k.endswith("_sq"))+CORE+(b.BIO,TOTAL,"log_target_species")
HIGH=.75
LOW=.25
MAX_LENGTH_RATIO=2.
MAX_ABS_Z=1.
MAX_RMS_Z=.5
MIN_PAIRS=100
MIN_BLOCKS=20
MAX_SMD=.10

def groups(rows):
    out=defaultdict(list)
    for r in rows:out[r["continent"],r["block"]].append(r)
    return out

def balance(pairs,scale,fields=MATCH):
    if not pairs:return None
    g=defaultdict(list)
    for a,c in pairs:g[a["continent"],a["block"]].append(np.array([a[k]-c[k] for k in fields])/scale)
    v=np.mean([np.mean(x,axis=0) for x in g.values()],axis=0)
    return {k:float(x) for k,x in zip(fields,v)}

def select(rows,fields=MATCH):
    if len(set(fields))!=len(fields):raise ValueError("Duplicate matching covariates")
    # Outcome-blind matching: neither trait_value nor any family/height moment
    # is read in this function. Coverage-based eligibility is predetermined.
    rr=[r for r in rows if r["active_normal_km"]+r["inactive_normal_km"]>0 and (r[SHARE]>=HIGH or r[SHARE]<=LOW)]
    a=[r for r in rr if r[SHARE]>=HIGH];c=[r for r in rr if r[SHARE]<=LOW]
    qa={"n_high":len(a),"n_low":len(c),"n_middle_or_no_rail":len(rows)-len(rr),
        "n_blocks_high":len(groups(a)),"n_blocks_low":len(groups(c)),
        "n_shared_blocks":len(set(groups(a))&set(groups(c))),
        "n_pairs":0,"n_blocks_matched":0}
    if not a or not c:return [],qa,None
    # Frozen pooled extreme-arm scale; constant fields have zero differences.
    scale=np.std(np.asarray([[r[f] for f in fields] for r in rr],float),axis=0)
    scale=np.where(scale>1e-10,scale,1.)
    ag,cg=groups(a),groups(c);pairs=[];n_feasible=0;can_high=set();can_low=set()
    for key in sorted(set(ag)&set(cg)):
        aa=sorted(ag[key],key=lambda r:r["key"]);cc=sorted(cg[key],key=lambda r:r["key"])
        x=np.asarray([[r[f] for f in fields] for r in aa],float)/scale
        y=np.asarray([[r[f] for f in fields] for r in cc],float)/scale
        d=np.abs(x[:,None,:]-y[None,:,:])
        rms=np.sqrt(np.mean(d*d,axis=2))
        ax=np.log([r["active_normal_km"]+r["inactive_normal_km"] for r in aa])
        cy=np.log([r["active_normal_km"]+r["inactive_normal_km"] for r in cc])
        ok=(np.abs(ax[:,None]-cy[None,:])<=math.log(MAX_LENGTH_RATIO)+1e-12)&(d.max(axis=2)<=MAX_ABS_Z+1e-12)&(rms<=MAX_RMS_Z+1e-12)
        n_feasible+=int(ok.sum())
        can_high.update(aa[i]["key"] for i in np.flatnonzero(ok.any(axis=1)))
        can_low.update(cc[j]["key"] for j in np.flatnonzero(ok.any(axis=0)))
        # One dummy per high row. Dummy penalty exceeds any sum of eligible
        # distances, ensuring cardinality first, minimum total distance second.
        penalty=len(aa)+1.
        cost=np.concatenate((np.where(ok,rms,penalty*3),np.full((len(aa),len(aa)),penalty)),axis=1)
        ii,jj=linear_sum_assignment(cost)
        for i,j in zip(ii,jj):
            if j>=len(cc):continue
            if not ok[i,j]:raise ValueError("Invalid match selected")
            pairs.append((aa[i],cc[j]))
    if len({r["key"] for p in pairs for r in p})!=2*len(pairs):raise ValueError("Match reused cell")
    qa.update(n_pairs=len(pairs),n_blocks_matched=len({(a["continent"],a["block"]) for a,c in pairs}),
              n_feasible_edges=n_feasible,n_high_with_candidate=len(can_high),n_low_with_candidate=len(can_low),
              matched_by_continent=dict(Counter(a["continent"] for a,c in pairs)))
    bb=balance(pairs,scale,fields)
    qa["balance_SMD"]=bb;qa["max_abs_SMD"]=max(map(abs,bb.values())) if bb else None
    if pairs:
        ratios=[max((a["active_normal_km"]+a["inactive_normal_km"])/(c["active_normal_km"]+c["inactive_normal_km"]),(c["active_normal_km"]+c["inactive_normal_km"])/(a["active_normal_km"]+a["inactive_normal_km"])) for a,c in pairs]
        qa["length_ratio_max"]=max(ratios);qa["length_ratio_median"]=float(np.median(ratios))
        qa["mean_share_gap"]=float(np.mean([a[SHARE]-c[SHARE] for a,c in pairs]))
    return pairs,qa,scale

def outcome(pairs,qa,target="trait_value"):
    if qa["n_pairs"]<MIN_PAIRS or qa["n_blocks_matched"]<MIN_BLOCKS:return {"status":"insufficient matched support"}
    if qa["max_abs_SMD"] is None or qa["max_abs_SMD"]>MAX_SMD:return {"status":"insufficient matched covariate balance"}
    g=defaultdict(list)
    for a,c in pairs:g[a["continent"],a["block"]].append(a[target]-c[target])
    # Equal spatial block weights; conservative common block influence.
    z={k:float(np.mean(v)) for k,v in g.items()}
    mu=float(np.mean(list(z.values())))
    scores={k:(v-mu)/len(z) for k,v in z.items()}
    return {"status":"estimated","active_dominant_minus_inactive_dominant":interval((mu,scores))}

def evaluate(rows,tm,old,gate):
    rr=attach(rows,tm,HEIGHT,gate)
    for label,flex,share in (("primary",True,False),("linear_effort",False,False),("fixed_total_active_share",True,True)):
        z,_=estimate(rr,flex,share);check_effects(z,old[label])
    pairs,qa,_=select(b.with_epoch(rr,2024))
    return {"support":qa,"height":outcome(pairs,qa)}

def main():
    p=argparse.ArgumentParser()
    for k in ("traits","previous_cells","previous","size_model","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();src,oc,prev,pm,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(a.traits,a.previous_cells,a.previous,a.size_model,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if pm["version"]!="railstatus-verified-scaled-plant-size-association-v1" or pm["trait_cells_sha256"]!=hashlib.sha256(a.traits.read_bytes()).hexdigest() or pm["source"]!=src["source"]:raise ValueError("Prior size model differs")
    height_panels={(p["width"],p["land"],p["coverage"]):p for p in pm["panels"]}
    if src["version"]!="railstatus-verified-scaled-plant-size-cells-v1" or not src["exact_full_split_target_and_life_replay"] or any(src[k]!=split[k] for k in ("obs_sha256","taxa_sha256","split_seed")) or src["source_SHA256"]!={"LifeForm":LIFE_SHA,"size_archive":SIZE_SHA} or src["prior_cells_sha256"]!=hashlib.sha256(a.previous_cells.read_bytes()).hexdigest():raise ValueError("Wrong verified size source")
    if prev["version"]!="railstatus-disturbance-trait-association-v1" or prev["trait_cells_sha256"]!=src["prior_cells_sha256"] or oc["version"]!="railstatus-disturbance-trait-cells-v1":raise ValueError("Prior model source mismatch")
    old_panels={(p["width"],p["land"],p["coverage"]):p for p in prev["panels"] if p["taxon"]=="plants"}
    if any(prior[k]!=src[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Original census mismatch")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1":raise ValueError("Wrong AGB")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen map/climate/pop mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in WIDTHS:
            grid=prior["grids"][str(w)];keys={b.key(r) for g in GROUPS for r in grid[g]}
            clim,_=b.source_climate_lookup(land,ras,scale,keys,w);ext,qa=climate_extra(a.climate,ras,clim,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            for g in ('plants',):
                tm={tuple(r["key"]):r for r in src["grids"][str(w)][g]}
                old=extract(split["pairs"][str(w)][g])
                if {k:(tuple(v["full"]),tuple(tuple(x) for x in v["halves"])) for k,v in tm.items()}!=old:raise ValueError("Exact split counts not replayed")
                rows,_,_=b.original_cell_rows([b.convert(v) for v in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w)
                for landgate in (.8,.5):
                    rr,_=b.enrich_cells(rows,rd,tr,w,landgate);rr,_=b.join_population(rr,pp,rd)
                    ready=[]
                    for r in rr:
                        k=r["key"]
                        if k not in tm or k not in ext or k not in ag20 or k not in ag24:continue
                        if tuple(tm[k]["full"])!=(r["observations"],r["observers"],r["species"]):raise ValueError("Exact full census mismatch")
                        av,cv=ag20[k],ag24[k]
                        if min(av["valid_area_m2"],cv["valid_area_m2"])<.9*rd[k]["source_GRIP4_total_land_area_m2"]:continue
                        if any(not math.isfinite(z["agb_Mg_ha"]) or z["agb_Mg_ha"]<0 for z in (av,cv)):raise ValueError("Invalid AGB")
                        ready.append(dict(r,agb2020=av["agb_Mg_ha"],agb2024=cv["agb_Mg_ha"],**ext[k]))
                    if len(ready)!=height_panels[w,landgate,.5]["n_base"]:raise ValueError("Prior base geography drift")
                    # Fixed coverage50; no searching coverage thresholds.
                    p0=height_panels[w,landgate,.5]["traits"][HEIGHT]
                    z={"width":w,"land":landgate,"coverage":.5,**evaluate(ready,tm,p0,.5)}
                    panels.append(z);print(json.dumps({"height_overlap_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-height-status-overlap-v1","source_replay":checks,
         "prior_height_model_sha256":hashlib.sha256(a.size_model.read_bytes()).hexdigest(),
         "trait_cells_sha256":hashlib.sha256(a.traits.read_bytes()).hexdigest(),"matching_fields":MATCH,
         "design":{"high_share":HIGH,"low_share":LOW,"length_ratio":MAX_LENGTH_RATIO,"max_abs_z":MAX_ABS_Z,"max_rms_z":MAX_RMS_Z,"min_pairs":MIN_PAIRS,"min_blocks":MIN_BLOCKS,"max_abs_SMD":MAX_SMD},
         "panels":panels,"primary":"100km land80 coverage50, 2025 outcome-blind local matched extreme active share comparison. 200km and land50 fixed secondary.",
         "limits":["conditional overlap diagnostic on previously analyzed geography; not independent validation",
                   "estimated effect only if fixed support and balance gates pass; failure is lack of identification, not null effect",
                   "matched extremes are not all active/inactive railways; static 2026 mapping, no closure or management intervention",
                   "residual confounding and missing trait/sampling selection remain; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

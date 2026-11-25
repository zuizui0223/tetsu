#!/usr/bin/env python3
"""Paired temporal cluster covariance and exact block influence of bird beta.
No raw observation scan. Replay previous annual endpoint FE coefficients.
"""
import argparse,json,math
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
from analyze_world_railstatus_observer_annual_connectivity import (
    cells_from_pairs,subset,YEARS,W,GROUPS,read,adapt_pair,split_pair,
    build_pairs,fit,PRED,ENDPOINT,ACTIVE,CATS,verify_original_sources,verify_replay,
    EXPECTED_SHA,MD5,POP_SHA,ORIGINAL_NE_SHA,WORLDCLIM_SHA,sha_file,
    read_source_builtup,source_continents,worldclim_tif,climate_scale,source_climate_lookup,
    make_lookup,original_cell_rows,inject_built,inject_status,enrich_cells,join_population,
    candidate_support,pair_key)
from analyze_world_railstatus_observer_split_connectivity import EFF_PAIR,FIXED_PAIR
from analyze_world_railstatus_direct_boundary_beta import INACTIVE,block_weights
EXPECTED={.5:{"2023":.012075710500332027,"2024":.021451930146742047,"2025":.018907744583935054},
          .8:{"2023":.006277147167015711,"2024":.014775610654230303,"2025":.014471121131399969}}
TERMS={"active":{ACTIVE:1.},"inactive":{INACTIVE:1.},"active_minus_inactive":{ACTIVE:1.,INACTIVE:-1.}}

def prepare(rows,fields,target,endpoint):
    if len(rows)<100:return None
    blocks=defaultdict(list)
    for i,r in enumerate(rows):blocks[r["continent"],r["block"]].append(i)
    if len(blocks)<20:return None
    keys=sorted(blocks)
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields]+[r[target]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Invalid annual design")
    if endpoint:
        nodes={k:i for i,k in enumerate(sorted({r[q] for r in rows for q in ("a","b")}))}
        # Projections are separable by spatial block only if endpoints
        # belong to a unique block. Explicitly verify this before deletion.
        nodeblocks={}
        for r in rows:
            key=r["continent"],r["block"]
            for k in (r["a"],r["b"]):
                if k in nodeblocks and nodeblocks[k]!=key:raise ValueError("Endpoint spans blocks")
                nodeblocks[k]=key
        ii=np.repeat(np.arange(len(rows)),2);jj=np.array([nodes[r[k]] for r in rows for k in ("a","b")])
        D=csr_matrix((np.repeat(np.sqrt(w),2),(ii,jj)),shape=(len(rows),len(nodes)))
        X*=np.sqrt(w)[:,None]
        for j in range(X.shape[1]):
            sol=lsmr(D,X[:,j],atol=1e-12,btol=1e-12,maxiter=4000)
            if sol[1] not in (0,1,2,4,5):raise ValueError("Endpoint FE projection failed")
            X[:,j]-=D@sol[0]
        if np.max(np.abs(D.T@X))>1e-7:raise ValueError("Endpoint FE residual check")
    else:
        for ix in blocks.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
        X*=np.sqrt(w)[:,None]
    A=X[:,:-1];y=X[:,-1]
    norm=np.sqrt(np.sum(A*A,axis=0))
    orig=np.sqrt(np.sum(w[:,None]*np.array([[r[k] for k in fields] for r in rows],float)**2,axis=0))
    keep=norm>np.maximum(1e-9,orig*1e-8);names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return None
    Z=A[:,keep]/norm[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return None
    residual_blocks=int(sum(np.linalg.norm(Z[ix,names.index(ACTIVE)])>1e-8 for ix in blocks.values()))
    if residual_blocks<20:return None
    grams=np.array([Z[blocks[k]].T@Z[blocks[k]] for k in keys])
    xy=np.array([Z[blocks[k]].T@y[blocks[k]] for k in keys])
    G=grams.sum(axis=0);h=xy.sum(axis=0);inv=np.linalg.inv(G);b=inv@h
    scores=(xy-np.einsum("bij,j->bi",grams,b))@inv.T
    terms={k:np.array([t.get(f,0.) for f in names])/norm[keep] for k,t in TERMS.items()}
    effects={}
    for k,c in terms.items():
        v=float(c@b);score=scores@c
        se=math.sqrt(len(keys)/(len(keys)-1)*float(np.sum(score**2)))
        lever=np.einsum("i,bij,j->b",inv@c,grams,inv@c)/float(c@inv@c)
        if np.min(lever)<-1e-9 or abs(float(lever.sum())-1)>1e-8:raise ValueError("Invalid contrast information shares")
        effects[k]={"estimate":v,"scores":score,"ci":[v-1.96*se,v+1.96*se],
                    "effective_information_blocks":float(1/np.sum(lever**2)),
                    "max_information_share":float(np.max(lever)),
                    "continent_information_shares":{a:float(sum(lever[i] for i,k in enumerate(keys) if k[0]==a)) for a in sorted({k[0] for k in keys})}}
    prior=fit(rows,target,endpoint,fields)
    if prior["status"]!="estimated":raise ValueError("Prior fit not reproduced")
    for k in TERMS:
        old=prior["active_minus_inactive"] if k=="active_minus_inactive" else prior["effects"][ACTIVE if k=="active" else INACTIVE]
        if old is None or abs(effects[k]["estimate"]-old["estimate"])>1e-9 or max(abs(x-y) for x,y in zip(effects[k]["ci"],old["interval95_approx"]))>1e-9:
            raise ValueError("Prior estimates or intervals not reproduced")
    return {"keys":keys,"G":G,"h":h,"grams":grams,"xy":xy,"terms":terms,"effects":effects,
            "residual_blocks":residual_blocks,"n_pairs":len(rows)}

def deletion(m,drop):
    keep=[i for i,k in enumerate(m["keys"]) if k not in drop]
    if len(keep)<20:return None
    G=m["grams"][keep].sum(axis=0);h=m["xy"][keep].sum(axis=0)
    s=np.linalg.eigvalsh(G)
    if s[-1]<=0 or s[0]<1e-16*s[-1]:return None
    b=np.linalg.solve(G,h)
    return {k:float(c@b) for k,c in m["terms"].items()}

def average(m,years,term):
    keys=m[years[0]]["keys"]
    b=np.array([m[y]["effects"][term]["estimate"] for y in years])
    scores=np.stack([m[y]["effects"][term]["scores"] for y in years])
    cov=len(keys)/(len(keys)-1)*scores@scores.T
    se=float(np.sqrt(np.sum(cov)/len(years)**2));v=float(b.mean())
    return {"estimate":v,"interval95_approx":[v-1.96*se,v+1.96*se],
            "annual_covariance":cov.tolist(),"years":list(years)}

def summarize(vals):
    if not vals:return {"n_valid":0}
    return {"n_valid":len(vals),"min":min(vals),"max":max(vals),"positive":sum(v>0 for v in vals)}

def annual_audit(rows,fields,target,endpoint):
    m={y:prepare(rows[y],fields,target,endpoint) for y in YEARS}
    if any(v is None for v in m.values()):return {"status":"insufficient residual or rank support","n_pairs":{y:len(rows[y]) for y in YEARS}}
    keys=m["2025"]["keys"]
    if any(v["keys"]!=keys for v in m.values()):raise ValueError("Annual cluster keys differ")
    ids=lambda v:{(r["a"],r["b"]) for r in v}
    if any(ids(rows[y])!=ids(rows["2025"]) for y in YEARS):raise ValueError("Annual pair support differs")
    out={"status":"estimated","n_blocks":len(keys),"n_pairs":len(rows["2025"]),
         "residual_blocks":{y:v["residual_blocks"] for y,v in m.items()},
         "years":{y:{k:{j:z for j,z in e.items() if j!="scores"} for k,e in v["effects"].items()} for y,v in m.items()},
         "validation_2023_2024":{k:average(m,YEARS[:2],k) for k in TERMS},
         "descriptive_all_years":{k:average(m,YEARS,k) for k in TERMS}}
    blocks={k:[] for k in TERMS};cont={}
    for key in keys:
        a=[deletion(m[y],{key}) for y in YEARS[:2]]
        if all(x is not None for x in a):
            for k in TERMS:blocks[k].append(sum(x[k] for x in a)/2)
    for c in sorted({k[0] for k in keys}):
        a=[deletion(m[y],{k for k in keys if k[0]==c}) for y in YEARS[:2]]
        cont[c]={k:sum(x[k] for x in a)/2 for k in TERMS} if all(x is not None for x in a) else None
    out["validation_leave_one_block_out"]={k:summarize(v) for k,v in blocks.items()}
    out["validation_leave_continent_out"]=cont
    return out
def main():
    p=argparse.ArgumentParser()
    for n in ("annual","split","prior","edges","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args()
    annual,src,prior,edge,plant,rail,status,tree,road,pop=map(read,(a.annual,a.split,a.prior,a.edges,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail);verify_replay(prior,plant)
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"]:raise ValueError("Wrong 2025 split source")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=W or annual["qa"]["rows"]!=288721938 or any(annual[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong annual source")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Prior SHA mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"]:raise ValueError("Invalid edge replay")
    if any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Status km mismatch")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    cand,cells=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Frozen candidates changed")
    source={y:{g:annual["pairs"][y][g] for g in GROUPS} for y in YEARS[:-1]}
    source["2025"]={g:[r for r in src["pairs"][str(W)][g] if "split_similarity" in r] for g in GROUPS}
    for y in YEARS:
        for g in GROUPS:
            keys={pair_key(r) for r in source[y][g]}
            if len(keys)!=len(source[y][g]) or not keys<=cand[g] or any("split_similarity" not in r for r in source[y][g]):raise ValueError("Annual pair candidates changed")
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);clim,_=source_climate_lookup(land,raster,scale,set.union(*cells.values()),W)
    finally:raster.close();mem.close()
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(W)][c]} for c in CATS}
    geo={}
    for y in YEARS:
        geo[y]={}
        for g in GROUPS:
            raw=cells_from_pairs(source[y][g])
            if y=="2025":
                pmap={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]):(r["observations"],r["observers"],r["species"]) for r in prior["grids"][str(W)][g]}
                for r in raw:
                    key=(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"])
                    v=(r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_distinct_observers_with_photo_linked_plant_records"],r["source_2025_plant_species_richness"])
                    if pmap.get(key)!=v:raise ValueError("2025 exact cell replay failed")
            rows,_,_=original_cell_rows(raw,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
            rows,_=inject_built(rows,built[str(W)]);rows,_=inject_status(rows,status,W);geo[y][g]=rows

    panels=[]
    for threshold in (.5,.8):
        rr={}
        for y in YEARS:
            rr[y]={}
            for g in GROUPS:
                m,_=enrich_cells(geo[y][g],rd,tr,W,threshold);m,_=join_population(m,pp,rd)
                full=build_pairs([adapt_pair(p) for p in source[y][g]],m,el)
                lookup={pair_key(p):p for p in source[y][g]}
                rr[y][g]=[split_pair(r,lookup[r["a"],r["b"]]) for r in full]
        keys={g:set.intersection(*({(r["a"],r["b"]) for r in rr[y][g]} for y in YEARS)) for g in GROUPS}
        common=keys["plants"]&keys["birds"]
        if len(keys["birds"])!={.5:819,.8:649}[threshold] or len(common)!={.5:735,.8:575}[threshold]:raise ValueError("Annual common-pair denominator changed")
        birds={y:subset(rr[y]["birds"],keys["birds"]) for y in YEARS}
        bird=annual_audit(birds,ENDPOINT,"cross_jaccard",True)
        if bird["status"]!="estimated":raise ValueError("Original annual bird estimator no longer identifiable")
        for y,v in EXPECTED[threshold].items():
            if abs(bird["years"][y]["active"]["estimate"]-v)>1e-9:raise ValueError("Prior annual bird coefficient not reproduced")
        # Same geography and joint covariate design as prior plant-bird paired contrast.
        paired={}
        for y in YEARS:
            plants={(r["a"],r["b"]):r for r in rr[y]["plants"]}
            birds2={(r["a"],r["b"]):r for r in rr[y]["birds"]}
            z=[]
            for k in sorted(common):
                p,b=plants[k],birds2[k]
                if any(p[f]!=b[f] for f in FIXED_PAIR):raise ValueError("Common geographic fields changed")
                x=dict(p,**{"other_"+f:b[f] for f in EFF_PAIR})
                x["plant_minus_bird"]=p["cross_jaccard"]-b["cross_jaccard"]
                z.append(x)
            paired[y]=z
        fields=ENDPOINT+tuple("other_"+f for f in EFF_PAIR if f.startswith("diff_"))
        panel={"land":threshold,"n_bird_pairs":len(keys["birds"]),"n_plant_bird_pairs":len(common),
               "bird_endpoint":bird,"bird_block":annual_audit(birds,PRED,"cross_jaccard",False),
               "plant_minus_bird_endpoint":annual_audit(paired,fields,"plant_minus_bird",True)}
        panels.append(panel);print(json.dumps({"bird_temporal_influence_panel":panel},allow_nan=False),flush=True)
    out={"version":"railstatus-bird-temporal-joint-influence-v1","replay":checks,"panels":panels,
         "primary":"200km land50 exact 819 common bird pairs; mean of 2023 and 2024 active crossing slopes, endpoint FE with block-matched covariance, prior discovery 2025 excluded from validation mean",
         "limits":["posthoc audit, retrospective holdback not prospective or independent ecosystems","equal-year mean of year-specific regressions, not a temporal trend","cluster approximation, blocks not proven independent","endpoint FE variation from cycles; leave-block and leave-continent diagnostics not independent tests","fixed 2026 mapped rail status and whole-landscape occurrence data, no train or movement mechanism"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

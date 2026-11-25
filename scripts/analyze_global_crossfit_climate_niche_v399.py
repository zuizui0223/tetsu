#!/usr/bin/env python3
"""v399: temporally and spatially cross-fitted realized plant climate envelopes."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from pyproj import Transformer
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
 source_client,open_original_source,taxonomy_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,MAX_UNCERT_M,OBS_QUALITY)
from analyze_global_multiyear_naturalized_v397 import entry,actual_annual_grid
from analyze_global_rail_documented_naturalized_v379 import grid_counts,grid_key
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
 verify_original_sources,source_climate_lookup,original_cell_rows,block_weights,CONTINENTS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import read_source_builtup,inject_built,ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from analyze_global_turnover_climate_block_v387 import load_extra_climate
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_global_naturalized_population_v392 import add_population
from analyze_global_naturalized_joint_population_local_v394 import FULL,POP_SHA
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
W=100000
REFERENCE_YEARS=("2018","2019","2020","2021","2022")
MIN_REF_CELLS=20
MIN_REF_BLOCKS=5
MIN_KNOWN=20
MIN_COVERAGE=.5
BOOT=399
SEED=39920261125
OUTCOMES=("warm_offset_C","cold_tail_offset_C","below_reference_q05_fraction","temp_width_C",
          "among_species_temperature_sd_C","dry_offset_log_precip","precip_width_log","trait_coverage")
FIELDS=FULL+("bio4_raw","log1p_bio12","bio15_raw")
TARGET="log_rail_km"
def scan(client,canonical,plant):
    ref=defaultdict(set);test=defaultdict(entry);counts=Counter();total=0
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    try:
        required={"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}
        if not required<=set(rows.fieldnames or []):raise ValueError("Original source schema drift")
        for r in rows:
            total+=1
            if total%10000000==0:print(json.dumps({"v399_original_rows":total,"eligible_records":counts}),flush=True)
            y=str(r["observed_on"])[:4]
            if y not in (*REFERENCE_YEARS,"2025") or not str(r["observed_on"]).startswith(y+"-"):continue
            if str(r["quality_grade"]).strip().lower().replace(" ","_").replace("-","_")!=OBS_QUALITY:continue
            try:
                tid=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"])
                lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if tid is None or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k is None:continue
            counts[y]+=1
            if y=="2025":
                a=test[k];a[0]+=1;a[1].add(tid)
                if r["observer_id"]:a[2].add(int(r["observer_id"]))
            else:ref[tid].add(k)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or raw.hash.hexdigest()!=plant["source_summary"]["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"] or total!=288721938 or counts["2025"]!=7827936:
        raise ValueError("Original snapshot full replay mismatch")
    old={grid_key(r):(r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_plant_species_richness"],r["source_2025_distinct_observers_with_photo_linked_plant_records"]) for r in grid_counts(plant,W)}
    if old!={k:(a[0],len(a[1]),len(a[2])) for k,a in test.items()}:raise ValueError("2025 original exact cell replay failed")
    return ref,test,{"all_source_rows":total,"eligible_by_year":dict(counts),"source_observation_SHA256":raw.hash.hexdigest(),
                    "exact_2025_all_cell_effort_richness_replay":True}
def envelope(vals,blocks,continents,target_block,target_continent,mode):
    # EPSG6933 block coordinates (1000 km projected), not physical buffers.
    if mode=="exclude_3x3_blocks":
        keep=np.max(np.abs(blocks-np.asarray(target_block)),axis=1)>1
    elif mode=="exclude_continent":
        keep=continents!=target_continent
    else:raise ValueError("Unknown exclusion")
    x=vals[keep];b=blocks[keep]
    if len(x)<MIN_REF_CELLS or len(np.unique(b,axis=0))<MIN_REF_BLOCKS:return None
    q=np.quantile(x,[.05,.95],axis=0);mu=x.mean(axis=0)
    return (float(mu[0]),float(q[0,0]),float(q[1,0]-q[0,0]),float(mu[1]),float(q[1,1]-q[0,1]))
def encode_refs(ref,geo,extra):
    keys=sorted(set(geo)&set(extra))
    idx={k:i for i,k in enumerate(keys)}
    vals=np.array([[geo[k]["bio1"],extra[k]["log1p_bio12"]] for k in keys],float)
    blocks=np.array([geo[k]["block"] for k in keys],int)
    cont=np.array([geo[k]["continent"] for k in keys])
    by={}
    for tid,ks in ref.items():
        if len(ks)<MIN_REF_CELLS:continue
        a=np.asarray(sorted(idx[k] for k in ks if k in idx),dtype=np.int32)
        if len(a)>=MIN_REF_CELLS:by[tid]=a
    return by,vals,blocks,cont
def community(rows,test,encoded,mode):
    by,vals,blocks,cont=encoded
    groups=defaultdict(list)
    for r in rows:
        key=(r["continent"],r["block"]) if mode=="exclude_3x3_blocks" else (r["continent"],)
        groups[key].append(r)
    o=[]
    for key,rr in sorted(groups.items()):
        cache={}
        for r in rr:
            ids=test[r["key"]][1];known=[]
            for tid in sorted(ids):
                if tid not in cache:
                    a=by.get(tid)
                    cache[tid]=None if a is None else envelope(vals[a],blocks[a],cont[a],r["block"],r["continent"],mode)
                if cache[tid] is not None:known.append(cache[tid])
            n=len(known);coverage=n/len(ids)
            x=dict(r,known_niche_species=n,trait_coverage=coverage)
            if n:
                A=np.asarray(known);t=r["bio1"];p=r["log1p_bio12"]
                x.update(warm_offset_C=float(A[:,0].mean()-t),cold_tail_offset_C=float(A[:,1].mean()-t),
                    below_reference_q05_fraction=float(np.mean(t<A[:,1])),temp_width_C=float(A[:,2].mean()),
                    among_species_temperature_sd_C=float(A[:,0].std()),
                    dry_offset_log_precip=float(p-A[:,3].mean()),precip_width_log=float(A[:,4].mean()))
            o.append(x)
    return o
def group_indices(rows,group_mode):
    g=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        b=(r["continent"],r["block"]);blocks[b].append(i)
        g[b if group_mode=="block" else r["continent"]].append(i)
    return [np.asarray(v,int) for v in g.values()],blocks
def slopes(X,ids,w,fields,targets):
    if w.sum()<=0:return None
    w=w/w.sum();Z=X.copy()
    for ix in ids:
        if w[ix].sum()>0:Z[ix]-=np.average(Z[ix],axis=0,weights=w[ix])
    Z*=np.sqrt(w)[:,None];n=len(fields);A=Z[:,:n];Y=Z[:,n:]
    if not np.isfinite(Z).all():raise ValueError("Nonfinite source design")
    out={}
    for k in targets:
        j=fields.index(k);x=A[:,j];other=np.delete(A,j,axis=1)
        norm=np.sqrt((other**2).sum(axis=0));ok=norm>1e-12;other=other[:,ok]/norm[ok]
        rx=x-other@np.linalg.lstsq(other,x,rcond=1e-10)[0]
        var=float(rx@rx)
        if var<=1e-10:return None
        ry=Y-other@np.linalg.lstsq(other,Y,rcond=1e-10)[0]
        out[k]=list(map(float,rx@ry/var))
    return out
def regression(rows,mode,fields,targets,outcomes,bootstrap=True):
    g,b=group_indices(rows,mode)
    o={"n_cells":len(rows),"n_blocks":len(b),"n_groups":len(g),"cells_by_continent":dict(Counter(r["continent"] for r in rows))}
    if len(rows)<100 or len(b)<20 or len(rows)-len(g)-len(fields)<30:
        return dict(o,status="insufficient support")
    X=np.asarray([[r[k] for k in fields+outcomes] for r in rows],float)
    w=block_weights(rows);c=slopes(X,g,w,fields,targets)
    if c is None:return dict(o,status="insufficient residual rail variation")
    o.update(status="estimated",coef={k:dict(zip(outcomes,v)) for k,v in c.items()})
    if not bootstrap:return o
    bc=defaultdict(list)
    for k in sorted(b):bc[k[0]].append(k)
    rng=np.random.default_rng(SEED);draws={k:[] for k in targets}
    for _ in range(BOOT):
        bw=np.zeros(len(rows))
        for bs in bc.values():
            m=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,n in zip(bs,m):bw[b[k]]=n/len(b[k])
        q=slopes(X,g,bw,fields,targets)
        if q is not None:
            for k in targets:draws[k].append(q[k])
    o["block_percentile_95_descriptive"]={k:{t:{
        "valid":len(v),"low":float(np.quantile(np.asarray(v)[:,j],.025)) if len(v)>=.95*BOOT else None,
        "high":float(np.quantile(np.asarray(v)[:,j],.975)) if len(v)>=.95*BOOT else None} for j,t in enumerate(outcomes)} for k,v in draws.items()}
    return o
def panel(rows,mode,gate):
    result={"reference_exclusion":mode,"min_land_fraction":gate,"n_environment_cells":len(rows),"fits":[]}
    for label,base in (("all",rows),("high_effort",[r for r in rows if r["observations"]>=200 and r["observers"]>=20])):
        # Coverage model includes cells without identifiable species niches.
        coverage=regression(base,"block",FIELDS,(TARGET,),("trait_coverage",))
        r=[dict(a,rail_cold=a[TARGET]*(15.-a["bio1"])/10.) for a in base if a["known_niche_species"]>=MIN_KNOWN and a["trait_coverage"]>=MIN_COVERAGE]
        old=regression(r,"continent",FIELDS,(TARGET,),OUTCOMES)
        local=regression(r,"block",FIELDS,(TARGET,),OUTCOMES)
        interaction=regression(r,"block",FIELDS+("rail_cold",),(TARGET,"rail_cold"),OUTCOMES)
        continents={c:regression([a for a in r if a["continent"]==c],"block",FIELDS,(TARGET,),OUTCOMES,False) for c in CONTINENTS}
        result["fits"].append({"support":label,"n_before_niche_gate":len(base),"n_after_niche_gate":len(r),
             "coverage_mean":float(np.mean([a["trait_coverage"] for a in base])) if base else None,
             "coverage_diagnostic":coverage,"continent_FE":old,"block_FE":local,
             "cold_interaction_block_FE":interaction,"by_continent":continents})
    print(json.dumps({"v399_panel":result},allow_nan=False),flush=True)
    return result
def main():
    p=argparse.ArgumentParser()
    for k in ("plants","rail","built","tree","road","population","ne","climate","out","aggregates"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text())
    plant,rail,tree,road,pop=map(read,(a.plants,a.rail,a.tree,a.road,a.population))
    verify_original_sources(plant,rail)
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong road")
    if pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong population")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong climate/map")
    client=source_client();canonical,_,_,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Original taxa SHA mismatch")
    ref,test,qa=scan(client,canonical,plant)
    keys=set(test)
    for ks in ref.values():keys.update(ks)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);geo,gqa=source_climate_lookup(land,raster,scale,keys,W)
    finally:raster.close();mem.close()
    extra,_,cqa=load_extra_climate(a.climate,set(geo),W)
    built,_=read_source_builtup(a.built)
    rows,_,_=original_cell_rows(actual_annual_grid(test),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],geo)
    rows,_=inject_built(rows,built[str(W)])
    if len(rows)!=3455:raise ValueError("Original 100km v376c denominator mismatch")
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    encoded=encode_refs(ref,geo,extra);del ref
    qa.update(reference_taxa_with_at_least_20_climate_cells=len(encoded[0]),climate_source_cells=len(geo),climate_complete_cells=len(extra))
    # No individual GPS, taxon IDs or observer IDs persisted.
    aggregate=[];out=[]
    for gate in LAND_THRESHOLDS:
        r,_=enrich_cells(rows,rd,tr,W,gate)
        if len(r)!={.5:3350,.8:2959}[gate]:raise ValueError("Original v382 denominator mismatch")
        r,pqa=add_population(r,pp,rd)
        r=[dict(x,**extra[x["key"]]) for x in r if x["key"] in extra]
        for mode in ("exclude_3x3_blocks","exclude_continent"):
            c=community(r,test,encoded,mode)
            aggregate.append({"reference_exclusion":mode,"min_land_fraction":gate,"rows":c})
            out.append(panel(c,mode,gate))
    a.aggregates.parent.mkdir(parents=True,exist_ok=True)
    a.aggregates.write_text(json.dumps({"version":"v399-crossfit-niche-coarse-aggregate","source_QA":qa,"panels":aggregate},allow_nan=False))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({"version":"v399-original-crossfit-global-plant-climate-niche","source_QA":qa,"climate_QA":cqa,"panels":out,
      "limits":["realized occurrence climatology, not physiological tolerance or experimentally measured rail microclimate",
       "2018-22 reference and 2025 evaluation from the same observation platform, spatial blocks exclude target and neighbors or continent, not independent field survey",
       "reference unique occupied coarse cells equally weighted, reference rail exposure not excluded",
       "100km center climate and historical macroclimate, not fine scale actual railway flora",
       "coverage threshold selects widespread, sufficiently recorded taxa; unknowns not assigned means",
       "cold-tail excursion outside empirical 5th percentile is not niche boundary violation",
       "temperature offset slope is equivalent to mean species temperature slope conditional on local BIO1",
       "no causal inference, block percentile intervals descriptive, no annual temporal warming trend inference"]},allow_nan=False))
if __name__=="__main__":main()

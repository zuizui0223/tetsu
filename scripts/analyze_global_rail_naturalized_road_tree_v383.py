#!/usr/bin/env python3
"""v383: exact v379d source replay; road/tree-adjusted regional naturalization."""
import argparse, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from analyze_global_rail_documented_naturalized_v379 import (
    source_lists,design,vascular_names,replay,within_region,score,
    REGIONAL_BASE,REGIONAL_RAIL,MIN_VASCULAR_SPECIES,MIN_REGION_CELLS,MIN_REGION_BLOCKS)
from audit_global_GloNAF_rail2048_v378a import load,MD5
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_client,taxonomy_source
from analyze_global_rail_road_tree_richness_v382 import (
    enrich_cells,make_lookup,LAND_THRESHOLDS)
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources,continent_holdout,block_weights,RICH_BASE,CONTINENTS,
    MIN_CONTINENT_RICH_CELLS)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import RICH_BUILT
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5 as ROAD_MD5
ORIGINAL_REPLAY={(100000,.3):(725,.5968714407556076),
    (100000,0.):(1182,.6692331525448875),
    (200000,.3):(130,.17987006014758986),
    (200000,0.):(240,.7100798976312451)}
BASE=RICH_BASE+(RICH_BUILT,"log_species")
ENV=("land_fraction","tree_consensus_fullcell","log_road_density")
PRE=BASE+("log_rail_km",)
ADJ=BASE+ENV+("log_rail_km",)
FULL=ADJ+("road_tree","rail_tree")
BOOT=399
SEED=38320261124

def rows_from_source(panels,sets,names,listed):
    out={}
    for key,items in panels.items():
        w,_=key
        rows=[]
        for r in items:
            ids=sets[w][r["key"]]; documented=listed[r["region"]]
            if len(ids)>r["species"]:raise ValueError("Vascular richness exceeds all-plant richness")
            if len(ids)<MIN_VASCULAR_SPECIES or not documented:continue
            n=sum(names[i] in documented for i in ids)
            rows.append(dict(r,known_fraction=n/len(ids),
                vascular_richness=len(ids),documented_richness=n,
                log_documented_regional_species_pool=math.log1p(len(documented))))
        expected_n,expected_beta=ORIGINAL_REPLAY[key]
        old=within_region(rows)
        if len(rows)!=expected_n or abs(old["coefficient_pp_per_log1p_rail_km"]-expected_beta)>1e-9:
            raise ValueError("Exact v379d fraction/cell/within-region reproduction failed")
        out[key]=rows
    return out

def eligible_region_rows(rows):
    g=defaultdict(list)
    for r in rows:g[r["region"]].append(r)
    good={k for k,a in g.items() if len(a)>=MIN_REGION_CELLS and
          len({(r["continent"],r["block"]) for r in a})>=MIN_REGION_BLOCKS and
          len({r["log_rail_km"] for r in a})>1}
    return [r for r in rows if r["region"] in good]

def design_arrays(rows,fields):
    X=np.array([[r[k] for k in fields]+[100*r["known_fraction"]] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Nonfinite original model values")
    groups=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        groups[r["region"]].append(i);blocks[(r["continent"],r["block"])].append(i)
    return X,[np.array(a,dtype=int) for a in groups.values()],blocks

def conditional_coefs(X,groups,w,fields,targets):
    """WLS region FE; individual coefficients via weighted residualization."""
    if not np.isfinite(w).all() or np.any(w<0) or w.sum()<=0:return {k:None for k in targets}
    w=w/w.sum();Z=X.copy()
    for ix in groups:
        if w[ix].sum()>0:Z[ix]-=np.average(Z[ix],axis=0,weights=w[ix])
    Z*=np.sqrt(w)[:,None]
    result={}
    for target in targets:
        j=fields.index(target);x=Z[:,j];y=Z[:,-1]
        other=np.delete(Z[:,:-1],j,axis=1)
        # Scale nuisance columns before pseudoinverse to avoid unit artifacts.
        scale=np.sqrt((other**2).sum(axis=0));active=scale>1e-12
        other=other[:,active]/scale[active]
        def resid(u):return u-other@np.linalg.lstsq(other,u,rcond=1e-10)[0]
        rx,ry=resid(x),resid(y);var=float(rx@rx)
        result[target]=float(rx@ry/var) if var>1e-10 else None
    return result

def regression(rows,fields,targets,bootstrap=False):
    regions={r["region"] for r in rows};blocks0={(r["continent"],r["block"]) for r in rows}
    out={"n_cells":len(rows),"n_regions":len(regions),"n_1000km_blocks":len(blocks0),
         "coefficients_pp":{k:None for k in targets}}
    if len(regions)<2 or len(rows)<30 or len(blocks0)<5:return out
    X,groups,blocks=design_arrays(rows,fields)
    w=block_weights(rows)
    out["coefficients_pp"]=conditional_coefs(X,groups,w,fields,targets)
    if not bootstrap:return out
    bycont=defaultdict(list)
    for b in sorted(blocks):bycont[b[0]].append(b)
    rng=np.random.default_rng(SEED)
    samples={k:[] for k in targets}
    for _ in range(BOOT):
        bw=np.zeros(len(rows))
        for c,bs in bycont.items():
            mult=rng.multinomial(len(bs),np.ones(len(bs))/len(bs))
            for b,n in zip(bs,mult):
                ix=blocks[b];bw[ix]=n/len(ix)
        co=conditional_coefs(X,groups,bw,fields,targets)
        for k in targets:
            if co[k] is not None:samples[k].append(co[k])
    out["spatial_block_bootstrap_percentile_95_descriptive"]={
        k:{"valid":len(s),"low":float(np.quantile(s,.025)) if len(s)>=.95*BOOT else None,
            "high":float(np.quantile(s,.975)) if len(s)>=.95*BOOT else None}
        for k,s in samples.items()}
    # Whole-region deletion: re-evaluate WLS weights on the remaining rows.
    deletions={k:[] for k in targets}
    for region in sorted(regions):
        sub=[r for r in rows if r["region"]!=region]
        x=regression(sub,fields,targets,False)["coefficients_pp"]
        for k,v in x.items():
            if v is not None:deletions[k].append(v)
    out["leave_one_region_out"]={k:{"valid":len(s),"min":min(s) if s else None,
        "max":max(s) if s else None,"positive":sum(x>0 for x in s),"negative":sum(x<0 for x in s)}
        for k,s in deletions.items()}
    return out

def model_panel(rows,roads,trees,w,pad,minland):
    matched,qa=enrich_cells(rows,roads,trees,w,minland)
    supported=eligible_region_rows(matched)
    # SAME region-supported rows for pre-adjustment and adjusted direction.
    pre=regression(supported,PRE,("log_rail_km",),True)
    adj=regression(supported,ADJ,("log_rail_km","log_road_density"),True)
    full=regression(supported,FULL,("log_rail_km","log_road_density","rail_tree","road_tree"),True)
    if pre["coefficients_pp"]["log_rail_km"] is not None:
        old=within_region(supported)["coefficient_pp_per_log1p_rail_km"]
        if old is None or abs(old-pre["coefficients_pp"]["log_rail_km"])>1e-8:
            raise ValueError("Independent within-region estimator disagreement")
    base=REGIONAL_BASE+("land_fraction","tree_consensus_fullcell","log_road_density")
    both=base+("log_rail_km",)
    rx=both+("road_tree",)
    fx=rx+("rail_tree",)
    cv={}
    for name,a,b in (("rail_after_road",base,both),("rail_tree_after_road_tree",rx,fx)):
        cv[name]=continent_holdout(matched,a,b,"known_fraction",True,
            MIN_CONTINENT_RICH_CELLS,MIN_CONTINENT_RICH_CELLS*2,name)
    out={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":minland,
        "v379d_cells_before_environment_gate":len(rows),"matched_cells":len(matched),
        "matched_by_continent":dict(Counter(r["continent"] for r in matched)),
        "excluded_environment":qa,"region_supported_by_continent":dict(Counter(r["continent"] for r in supported)),
        "same_support_pre_adjustment":pre,"same_support_road_tree_main_adjustment":adj,
        "same_support_interactions":full,"prediction":cv,
        "within_region_interaction_by_continent":{c:regression(
            eligible_region_rows([r for r in matched if r["continent"]==c]),FULL,
            ("log_rail_km","log_road_density","rail_tree","road_tree"),False) for c in CONTINENTS}}
    print(json.dumps({"v383_panel":out},allow_nan=False),flush=True)
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("plants","rail","built","ne","climate","glonaf","road","tree","out","cell_counts"):
        p.add_argument("--"+k.replace("_","-"),dest=k,type=Path,required=True)
    a=p.parse_args()
    read=lambda f:json.loads(f.read_text())
    plant,rail,road,tree=map(read,(a.plants,a.rail,a.road,a.tree))
    verify_original_sources(plant,rail)
    if road["source_zip_MD5"]!=ROAD_MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not road["version"].startswith("v381b"):
        raise ValueError("Not verified original v381b road source")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:
        raise ValueError("Not original v380b source tree")
    g,listed,lqa=source_lists(a.glonaf)
    # Audit the publisher row-status distribution before replay, without
    # silently changing the v379d outcome definition.
    status=Counter(r["status"] for r in load(a.glonaf/"glonaf_flora2.csv",{"status"}))
    print(json.dumps({"original_GloNAF_flora_status_counts":status}),flush=True)
    panels,geoqa=design(plant,rail,a.built,a.ne,a.climate,g)
    client=source_client();canonical,_,_,_,sha=taxonomy_source(client)
    names,sha2=vascular_names(client,canonical)
    if sha!=sha2 or sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:
        raise ValueError("Original taxonomy replay mismatch")
    sets,sourceqa=replay(client,canonical,names,plant,panels)
    rows=rows_from_source(panels,sets,names,listed)
    # Store only coarse-cell aggregates, never raw occurrence/species/observer IDs.
    counts={"version":"v383-exact-v379d-coarse-cell-aggregate","source_replay":sourceqa,
        "panels":[{"width_m":w,"region_guard_degrees":p,"rows":rr} for (w,p),rr in rows.items()]}
    a.cell_counts.parent.mkdir(parents=True,exist_ok=True)
    a.cell_counts.write_text(json.dumps(counts,allow_nan=False))
    o=[]
    for (w,pad),rr in rows.items():
        rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        for land in LAND_THRESHOLDS:o.append(model_panel(rr,rd,tr,w,pad,land))
    out={"version":"v383-original-v379d-naturalized-fraction-road-tree-retest",
        "exact_v379d_reproduction":True,"original_source_replay":sourceqa,
        "glonaf_source_QA":lqa,"source_original_flora_status_counts":dict(status),
        "original_geo_QA":geoqa,"bootstrap_draws":BOOT,"bootstrap_seed":SEED,
        "panels":o,"limits":["observed GloNAF-documented species fraction, nonlisted not native",
        "same original coarse grid, not rail-local plants or direct invasion/establishment",
        "regional source coverage highly uneven and no direct precision geoprivacy audit",
        "historical roads/tree consensus, not current forest canopy; tree is full-cell fraction",
        "block bootstrap conditional on sampled region-support and current data, not causal uncertainty",
        "no individual source GPS, taxa IDs, observer IDs in output"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

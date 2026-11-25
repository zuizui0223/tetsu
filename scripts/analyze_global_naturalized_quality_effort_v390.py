#!/usr/bin/env python3
"""v390: GloNAF list-quality and common geography effort comparison."""
import argparse, hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from audit_global_GloNAF_rail2048_v378a import MD5 as GLO_MD5, load,unique_map,complete_vascular,norm
from analyze_global_naturalized_count_decomposition_v389 import (
    source_rows,fit,conditional,FIELDS,OUTCOMES,TARGETS,EXPECTED,BOOT)
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows,regression,ADJ
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
SEED=39020261125
PRIOR={
 (100000,.3,.5):(.06376788599199941,.0948168156573006,197),
 (100000,.3,.8):(.067639900673987,.10048558927487392,193),
 (100000,0.,.5):(.06610531802001234,.07944742470266788,375),
 (100000,0.,.8):(.06359517051768895,.08004512330214006,366),
 (200000,.3,.5):(.034794799602417235,None,18),
 (200000,.3,.8):(.059423576572170136,None,18),
 (200000,0.,.5):(.08623556116781521,None,28),
 (200000,0.,.8):(.09715623038132795,None,28)}
def list_quality(root):
    for f in ("glonaf_list.csv","glonaf_region.csv"):
        with (root/f).open("rb") as z:
            if hashlib.file_digest(z,"md5").hexdigest()!=GLO_MD5[f]:
                raise ValueError("Original GloNAF source changed: "+f)
    regions=unique_map(load(root/"glonaf_region.csv",{"id","code"}),"id")
    lists=unique_map(load(root/"glonaf_list.csv",{"id","region_id","taxon_group_id",
                                               "completeness_id","invasives_only","naturalized_explicit_id"}),"id")
    if len(regions)!=1343 or len(lists)!=2774:
        raise ValueError("Original source list or region rows changed")
    complete=[r for r in lists.values() if complete_vascular(r)]
    if len(complete)!=473:raise ValueError("Original complete vascular-list count changed")
    good={norm(r["region_id"]) for r in complete}
    explicit={norm(r["region_id"]) for r in complete if norm(r["naturalized_explicit_id"])=="1"}
    if not good<=set(regions):raise ValueError("Qualifying region absent from original table")
    return set(regions),good,explicit,{"original_list_rows":len(lists),"original_region_rows":len(regions),
        "complete_vascular_list_rows":len(complete),"regions_with_qualifying_list":len(good),
        "regions_with_qualifying_explicit_naturalized_list":len(explicit)}
def high(r):return r["observations"]>=200 and r["observers"]>=20
def spatial_key(r):return r["region"],r["continent"],r["block"]
def common_support(rows):
    """Common high/low effort geography, not an outcome-dependent gate."""
    g=defaultdict(list)
    for r in rows:g[spatial_key(r)].append(r)
    keys={k for k,v in g.items() if sum(high(r) for r in v)>=2 and sum(not high(r) for r in v)>=2}
    return [r for r in rows if spatial_key(r) in keys],keys
def arrays(rows):
    groups=defaultdict(list);blocks=defaultdict(list)
    for i,r in enumerate(rows):
        groups[spatial_key(r)].append(i)
        blocks[(r["continent"],r["block"])].append(i)
    X=np.asarray([[r[k] for k in FIELDS+OUTCOMES] for r in rows],float)
    return X,[np.array(v,dtype=int) for v in groups.values()],blocks
def shared_effort_comparison(rows):
    # Both arms occupy exactly the same region x block strata.
    # Fit separate arm-specific slopes with region x block FE.
    selected,keys=common_support(rows)
    h=[r for r in selected if high(r)];l=[r for r in selected if not high(r)]
    blocks=sorted({(r["continent"],r["block"]) for r in selected})
    out={"n_original_region_block_strata":len(keys),"n_original_regions":len({k[0] for k in keys}),
         "n_spatial_blocks":len(blocks),"n_high":len(h),"n_low":len(l),
         "high_by_continent":dict(Counter(r["continent"] for r in h)),
         "low_by_continent":dict(Counter(r["continent"] for r in l))}
    if len(blocks)<10 or len(keys)<10 or min(len(h),len(l))<30:
        out["status"]="insufficient shared geographic support";return out
    H,hg,hb=arrays(h);L,lg,lb=arrays(l)
    if set(hb)!=set(lb) or set(hb)!=set(blocks):raise ValueError("Shared effort geography mismatch")
    a=conditional(H,hg,block_weights(h),FIELDS)
    b=conditional(L,lg,block_weights(l),FIELDS)
    if a is None or b is None:
        out["status"]="insufficient residual transport variation";return out
    out["status"]="descriptive comparison on same region-block strata"
    out["high_coefficients"]=a;out["low_coefficients"]=b
    out["high_minus_low"]={k:{t:a[k][t]-b[k][t] for t in OUTCOMES} for k in TARGETS}
    bycont=defaultdict(list)
    for k in blocks:bycont[k[0]].append(k)
    rng=np.random.default_rng(SEED)
    draws={k:{t:[] for t in OUTCOMES} for k in TARGETS}
    for _ in range(BOOT):
        wh=np.zeros(len(h));wl=np.zeros(len(l))
        for bs in bycont.values():
            mult=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            for k,n in zip(bs,mult):
                wh[hb[k]]=n/len(hb[k]);wl[lb[k]]=n/len(lb[k])
        aa=conditional(H,hg,wh,FIELDS);bb=conditional(L,lg,wl,FIELDS)
        if aa is None or bb is None:continue
        for k in TARGETS:
            for t in OUTCOMES:draws[k][t].append(aa[k][t]-bb[k][t])
    out["paired_block_resampling_difference_95_descriptive"]={
        k:{t:{"valid":len(v),"low":float(np.quantile(v,.025)) if len(v)>=.95*BOOT else None,
              "high":float(np.quantile(v,.975)) if len(v)>=.95*BOOT else None}
           for t,v in ss.items()} for k,ss in draws.items()}
    return out
def run(cells,roads,tree,glonaf):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:
        raise ValueError("Not original verified v383 aggregate")
    if roads["source_zip_MD5"]!=MD5 or roads["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or not roads["version"].startswith("v381b"):
        raise ValueError("Not original v381b road")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:
        raise ValueError("Not original v380b tree")
    known,complete,explicit,qa=list_quality(glonaf)
    panels=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];rows=source_rows(p)
        if any(r["region"] not in known for r in rows):raise ValueError("Unknown original region")
        rd=make_lookup(roads["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        for gate in LAND_THRESHOLDS:
            m,_=enrich_cells(rows,rd,tr,w,gate)
            whole=eligible_region_rows(m);h=eligible_region_rows([r for r in m if high(r)])
            n,s,b=EXPECTED[w,pad,gate];d,dh,nh=PRIOR[w,pad,gate]
            old=regression(whole,ADJ,("log_rail_km",),False)
            a=fit(whole,FIELDS,False);ah=fit(h,FIELDS,False)
            if len(m)!=n or len(whole)!=s or abs(old["coefficients_pp"]["log_rail_km"]-b)>1e-9 or len(h)!=nh:
                raise ValueError("Original v383/v389 support mismatch")
            if abs(a["coefficients"]["log_rail_km"]["log1p_documented"]-d)>1e-9:
                raise ValueError("Original v389 coefficient mismatch")
            if dh is not None and (ah["coefficients"] is None or abs(ah["coefficients"]["log_rail_km"]["log1p_documented"]-dh)>1e-9):
                raise ValueError("Original v389 high effort coefficient mismatch")
            out={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":gate,
                 "n_matched":len(m),"n_whole_supported":len(whole),
                 "n_high_supported":len(h),"n_qualifying_list_cells_before_region_support":sum(r["region"] in complete for r in m),
                 "n_qualifying_explicit_list_cells_before_region_support":sum(r["region"] in explicit for r in m),
                 "verified_v389_replay":True,"models":[]}
            high_regions={r["region"] for r in h}
            high_keys={spatial_key(r) for r in h}
            subsets={
                "all":whole,"high_effort":h,
                "all_in_high_effort_supported_regions":eligible_region_rows([r for r in m if r["region"] in high_regions]),
                "all_in_high_effort_supported_region_blocks":eligible_region_rows([r for r in m if spatial_key(r) in high_keys]),
                "complete_list_regions_all":eligible_region_rows([r for r in m if r["region"] in complete]),
                "complete_list_regions_high":eligible_region_rows([r for r in m if r["region"] in complete and high(r)]),
                "no_qualifying_complete_list_all":eligible_region_rows([r for r in m if r["region"] not in complete]),
                "no_qualifying_complete_list_high":eligible_region_rows([r for r in m if r["region"] not in complete and high(r)])}
            for label,rr in subsets.items():
                f=fit(rr,FIELDS,True)
                out["models"].append({"support":label,"fit":f})
            out["common_region_block_high_low"]=shared_effort_comparison(whole)
            panels.append(out)
            print(json.dumps({"v390_panel":out},allow_nan=False),flush=True)
    return {"version":"v390-original-v389-list-quality-and-common-geography-effort",
         "source_list_quality_QA":qa,"panels":panels,"shared_effort_draws":BOOT,"shared_effort_seed":SEED,
         "limits":["complete/close-to-complete list criterion is a metadata quality indicator, NOT perfect detection or proof that nonlisted taxa are native",
         "posthoc same observational sample; shared region-block support is not random assignment, different cells remain",
         "unmatched or low-support geographies not extrapolated",
         "block bootstrap conditional descriptive stability; small blocks and correlated observations limit inference",
         "historical roads/tree, no rail-local or causal transport inference"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","glonaf","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.cells),read(a.road),read(a.tree),a.glonaf);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()

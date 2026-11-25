#!/usr/bin/env python3
"""v398: outcome-blind overlap graph support on the same v396 caliper."""
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from analyze_global_naturalized_count_decomposition_v389 import source_rows,OUTCOMES,BOOT,EXPECTED
from analyze_global_naturalized_quality_effort_v390 import spatial_key,high
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_global_naturalized_population_v392 import add_population
from analyze_global_naturalized_local_swap_v393 import supported_local
from analyze_global_naturalized_joint_population_local_v394 import POP_SHA,MATCH,LOCAL_EXPECT
from analyze_global_rail_presence_intensity_v395 import add_exposure,mixed
from analyze_global_rail_presence_matched_v396 import standardized,optimal_pairs,CALIPER,MIN_PAIRS,MIN_BLOCKS
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
SEED=39820261125
MIN_ESS=30
MAX_BALANCE=.20
EXPECTED_PAIRS={(100000,.3,.5):(15,2),(100000,.3,.8):(14,2),(100000,0.,.5):(33,11),(100000,0.,.8):(31,11),
                (200000,.3,.5):(0,0),(200000,.3,.8):(0,0),(200000,0.,.5):(0,0),(200000,0.,.8):(0,0)}
def graph(rows,Z):
    strata=defaultdict(list)
    for i,r in enumerate(rows):strata[spatial_key(r)].append(i)
    adj={}
    for ix in strata.values():
        a=[i for i in ix if rows[i]["rail_present"]==0]
        b=np.array([i for i in ix if rows[i]["rail_present"]==1],int)
        if not len(b):continue
        for i in a:
            ok=np.max(np.abs(Z[b]-Z[i]),axis=1)<=CALIPER+1e-12
            if np.any(ok):adj[i]=b[ok].tolist()
    return adj
def ess(w):
    w=np.asarray(w,float)
    return float(w.sum()**2/(w@w)) if len(w) and w@w>0 else 0.
def make_weights(rows,adj,blocks=None):
    # Each block has equal mass; each eligible zero within block equal
    # mass; each admissible positive neighbour of that zero equal share.
    g=defaultdict(list)
    for i in adj:g[(rows[i]["continent"],rows[i]["block"])].append(i)
    if blocks is None:blocks={b:1 for b in g}
    a=np.zeros(len(rows));b=np.zeros(len(rows))
    for k,m in blocks.items():
        if not m:continue
        ix=g[k]
        for i in ix:
            w=m/len(ix);a[i]+=w
            b[adj[i]]+=w/len(adj[i])
    if a.sum()<=0 or b.sum()<=0:raise ValueError("No support for graph weights")
    a/=a.sum();b/=b.sum()
    return a,b
def audit(rows,seed,expected):
    if not rows:
        if expected:raise ValueError("Previous match count changed")
        return {"status":"no mixed strata","n_pairs":0}
    Z=standardized(rows)
    opt=optimal_pairs(rows,Z)
    if len(opt)!=expected:raise ValueError("Exact v396 maximum-cardinality match count mismatch")
    adj=graph(rows,Z)
    pos={j for js in adj.values() for j in js}
    blocks=sorted({(rows[i]["continent"],rows[i]["block"]) for i in adj})
    o={"n_mixed_cells":len(rows),"n_zero":sum(r["rail_present"]==0 for r in rows),
       "n_positive":sum(r["rail_present"]==1 for r in rows),"n_pairs":len(opt),
       "n_zero_with_admissible_partner":len(adj),"n_positive_with_admissible_partner":len(pos),
       "n_admissible_edges":sum(map(len,adj.values())),"n_blocks":len(blocks),
       "supported_zeros_by_continent":dict(Counter(rows[i]["continent"] for i in adj))}
    if not adj:return dict(o,status="no admissible edges")
    a,b=make_weights(rows,adj)
    bal=np.abs((b-a)@Z)
    o.update(zero_ESS=ess(a),positive_ESS=ess(b),max_zero_weight=float(a.max()),max_positive_weight=float(b.max()),
             max_abs_weighted_standardized_difference=float(bal.max()),
             abs_weighted_standardized_difference=dict(zip(MATCH,map(float,bal))),
             zero_with_partner_fraction=len(adj)/o["n_zero"],
             n_zero_with_partner_not_in_optimal_pairing=len(adj)-len(opt))
    if min(len(adj),len(pos))<MIN_PAIRS or len(blocks)<MIN_BLOCKS or min(ess(a),ess(b))<MIN_ESS or bal.max()>MAX_BALANCE:
        return dict(o,status="insufficient unique cells, effective sample, block support or balance; no biological estimate")
    Y=np.array([[r[k] for k in OUTCOMES] for r in rows],float)
    contrast=(b-a)@Y
    if abs(contrast[0]-contrast[1]-contrast[2])>1e-9:raise ValueError("Nonadditive contrast")
    o["graph_weighted_difference"]=dict(zip(OUTCOMES,map(float,contrast)))
    bc=defaultdict(list)
    for k in blocks:bc[k[0]].append(k)
    rng=np.random.default_rng(seed);draws=[]
    for _ in range(BOOT):
        m={}
        for bs in bc.values():
            v=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
            m.update(zip(bs,map(int,v)))
        aa,bb=make_weights(rows,adj,m);draws.append((bb-aa)@Y)
    D=np.asarray(draws)
    o["block_percentile_95_descriptive"]={t:{"low":float(np.quantile(D[:,i],.025)),
       "high":float(np.quantile(D[:,i],.975)),"valid":len(D)} for i,t in enumerate(OUTCOMES)}
    o["status"]="estimated on caliper-overlap subset, descriptive only"
    return o
def run(cells,roads,tree,pop):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:raise ValueError("Wrong original cells")
    if roads["source_zip_MD5"]!=MD5 or roads["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong roads")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree")
    if pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong population")
    out=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];rr=source_rows(p)
        rd=make_lookup(roads["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
        for land in LAND_THRESHOLDS:
            k=w,pad,land;m,_=enrich_cells(rr,rd,tr,w,land);z,_=add_population(m,pp,rd)
            if len(m)!=EXPECTED[k][0] or len(z)!=len(m):raise ValueError("Original source support changed")
            panel={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":land,"fits":[]}
            for i,mode in enumerate(("all","high_effort")):
                r=eligible_region_rows(z if i==0 else [x for x in z if high(x)])
                local=[add_exposure(x) for x in supported_local(r)]
                if len(local)!=LOCAL_EXPECT[k][i][0]:raise ValueError("v396 local support changed")
                panel["fits"].append({"support":mode,"overlap":audit(mixed(local),SEED+w+int(pad*100)+int(land*10)+i,EXPECTED_PAIRS[k][i])})
            out.append(panel);print(json.dumps({"v398_panel":panel},allow_nan=False),flush=True)
    return {"version":"v398-original-v396-admissible-overlap-graph",
            "panels":out,"caliper":CALIPER,"matching_covariates":MATCH,"minimum_effective_n_each_arm":MIN_ESS,
            "minimum_blocks":MIN_BLOCKS,"maximum_mean_standardized_imbalance":MAX_BALANCE,
            "limits":["same fixed source, no new field observations or additional years",
                      "graph weighting targets supported mapped-zero cells with admissible positive neighbours, not global landscape",
                      "partner reuse is not independent replication; Kish ESS is weight concentration not spatial effective n",
                      "block resampling descriptive, no causal interpretation; unlisted not native"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","population","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.cells),read(a.road),read(a.tree),read(a.population));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()

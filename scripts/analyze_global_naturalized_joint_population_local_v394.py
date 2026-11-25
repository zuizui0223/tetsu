#!/usr/bin/env python3
"""v394 joint population and local region-block control of observed naturalized richness."""
import argparse,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from analyze_global_naturalized_count_decomposition_v389 import source_rows,fit,conditional,FIELDS,OUTCOMES,TARGETS,EXPECTED,BOOT
from analyze_global_naturalized_quality_effort_v390 import PRIOR,spatial_key,high
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import block_weights
from analyze_global_naturalized_population_v392 import add_population
from analyze_global_naturalized_local_swap_v393 import supported_local
from validate_global_rail_local_matched_placebo_v385 import FEATURES,CALIPER,swapped
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
POP_SHA="846c03076da12d9d998e6b9946a0774ca134ea76b0731d11483058a7f9185bbf"
FULL=FIELDS+("log_population_density",)
MATCH=FEATURES+("log_population_density",)
SEED=39420261125
DRAWS=199
# panel order: strict 50%, strict 80%, relaxed 50%, relaxed 80%.
POP_EXPECT={
 (100000,.3,.5):(.048715221507353475,.06901818740265908),
 (100000,.3,.8):(.052335183142607655,.07247223423034625),
 (100000,0.,.5):(.055925199639997404,.06799561222983505),
 (100000,0.,.8):(.05535887154591203,.0696031780650326),
 (200000,.3,.5):(.035554756881568596,None),
 (200000,.3,.8):(.06175016670768389,None),
 (200000,0.,.5):(.07703934515362282,None),
 (200000,0.,.8):(.08961114488017524,None)}
LOCAL_EXPECT={
 (100000,.3,.5):((354,.023570395765293985),(178,.09766881904617292)),
 (100000,.3,.8):((343,.026313841432427564),(172,.1022665322178685)),
 (100000,0.,.5):((576,.03805909452604064),(334,.06079733536086792)),
 (100000,0.,.8):((557,.03843211286531771),(325,.05726289907127122)),
 (200000,.3,.5):((35,None),(11,None)),
 (200000,.3,.8):((31,None),(11,None)),
 (200000,0.,.5):((53,None),(17,None)),
 (200000,0.,.8):((48,None),(17,None))}
def groups(rows):
    g=defaultdict(list);b=defaultdict(list)
    for i,r in enumerate(rows):
        g[spatial_key(r)].append(i);b[(r["continent"],r["block"])].append(i)
    return [np.array(v,int) for v in g.values()],b
def mat(rows,fields):
    return np.array([[r[k] for k in fields+OUTCOMES] for r in rows],float)
def supported(rows,fields):
    g,b=groups(rows)
    return len(rows)>=40 and len(g)>=10 and len(b)>=10 and len(rows)-len(g)-len(fields)>=10
def match(rows):
    cont=defaultdict(list);strata=defaultdict(list)
    for i,r in enumerate(rows):
        cont[r["continent"]].append(i);strata[spatial_key(r)].append(i)
    X=np.array([[r[k] for k in MATCH] for r in rows],float);Z=X.copy()
    for ix in cont.values():
        sd=X[ix].std(axis=0);sd[sd<1e-12]=1
        Z[ix]=(X[ix]-X[ix].mean(axis=0))/sd
    pairs=[]
    for ix in strata.values():
        candidates=[]
        for n,i in enumerate(ix):
            for j in ix[n+1:]:
                d=np.abs(Z[i]-Z[j])
                if d.max()<=CALIPER:candidates.append((float(d@d),min(i,j),max(i,j)))
        used=set()
        for _,i,j in sorted(candidates):
            if i in used or j in used:continue
            used.update((i,j));pairs.append((i,j))
    return pairs
def summarize(draws):
    return {k:{t:{"valid":len(v),"low":float(np.quantile(v,.025)) if len(v)>=.95*BOOT else None,
        "high":float(np.quantile(v,.975)) if len(v)>=.95*BOOT else None} for t,v in d.items()} for k,d in draws.items()}
def newdraws():return {k:{t:[] for t in OUTCOMES} for k in TARGETS}
def compare(rows,expected,seed):
    r=supported_local(rows);g,b=groups(r)
    n,old_d=expected
    if len(r)!=n:raise ValueError("v393 original local support mismatch")
    out={"n_input":len(rows),"n_cells":len(r),"n_region_block_strata":len(g),"n_blocks":len(b),
         "n_regions":len({x["region"] for x in r}),"cells_by_continent":dict(Counter(x["continent"] for x in r))}
    if not supported(r,FIELDS):
        if old_d is not None:raise ValueError("Expected original v393 local support")
        out["status"]="insufficient original local support";return out
    X=mat(r,FIELDS);Z=mat(r,FULL);w=block_weights(r)
    old=conditional(X,g,w,FIELDS)
    if old is None or old_d is None or abs(old["log_rail_km"]["log1p_documented"]-old_d)>1e-9:
        raise ValueError("v393 original local coefficient not reproduced")
    out["same_support_local_before_population"]=old
    if not supported(r,FULL):
        out["status"]="insufficient additional-model support";return out
    new=conditional(Z,g,w,FULL)
    if new is None:out["status"]="insufficient residual variation";return out
    out["status"]="joint local and population adjustment estimated"
    out["same_support_local_after_population"]=new
    out["same_support_population_region_FE"]=fit(r,FULL,False)
    out["coefficient_change_after_minus_before"]={k:{t:new[k][t]-old[k][t] for t in OUTCOMES} for k in TARGETS}
    bycont=defaultdict(list)
    for k in sorted(b):bycont[k[0]].append(k)
    rng=np.random.default_rng(seed)
    draws={"before":newdraws(),"after":newdraws(),"paired_change":newdraws()}
    for _ in range(BOOT):
        bw=np.zeros(len(r))
        for keys in bycont.values():
            mul=rng.multinomial(len(keys),np.full(len(keys),1/len(keys)))
            for k,m in zip(keys,mul):
                ix=b[k];bw[ix]=m/len(ix)
        a=conditional(X,g,bw,FIELDS);c=conditional(Z,g,bw,FULL)
        if a is None or c is None:continue
        for k in TARGETS:
            for t in OUTCOMES:
                draws["before"][k][t].append(a[k][t])
                draws["after"][k][t].append(c[k][t])
                draws["paired_change"][k][t].append(c[k][t]-a[k][t])
    out["block_percentile_95_descriptive"]={k:summarize(d) for k,d in draws.items()}
    # Whole block deletion and independent continental direction.
    for label,units,keyfunc in (("leave_one_block_out",sorted(b),lambda x:(x["continent"],x["block"])),
                                ("leave_one_continent_out",sorted({x["continent"] for x in r}),lambda x:x["continent"])):
        vals=newdraws()
        for unit in units:
            sub=[x for x in r if keyfunc(x)!=unit]
            if not supported(sub,FULL):continue
            gg,_=groups(sub);co=conditional(mat(sub,FULL),gg,block_weights(sub),FULL)
            if co is None:continue
            for k in TARGETS:
                for t in OUTCOMES:vals[k][t].append(co[k][t])
        out[label]={"n_attempted":len(units),"ranges":{k:{t:{"valid":len(v),"min":min(v) if v else None,
              "max":max(v) if v else None,"positive":sum(x>0 for x in v)} for t,v in d.items()} for k,d in vals.items()}}
    pairs=match(r);mov=[(i,j) for i,j in pairs if r[i]["rail_km"]!=r[j]["rail_km"]]
    qa={"matched_cells":2*len(pairs),"movable_cells":2*len(mov),
        "movable_fraction":2*len(mov)/len(r),"movable_blocks":len({(r[i]["continent"],r[i]["block"]) for i,j in mov})}
    out["population_matched_swap_support"]=qa
    if qa["movable_fraction"]<.2 or qa["movable_blocks"]<10:
        out["placebo_status"]="insufficient mobility";return out
    vals={t:[] for t in OUTCOMES};col=FULL.index("log_rail_km")
    for _ in range(DRAWS):
        f,_=swapped(r,pairs,rng);P=Z.copy();P[:,col]=[x["log_rail_km"] for x in f]
        cf=conditional(P,g,w,FULL)
        if cf is not None:
            for t in OUTCOMES:vals[t].append(cf["log_rail_km"][t])
    out["placebo_status"]="technical matched swaps, not randomization inference"
    out["technical_rank"]={t:{"actual":new["log_rail_km"][t],"valid":len(v),
        "fake_mean":float(np.mean(v)) if v else None,
        "fake_at_least_actual_plus_one_fraction":(1+sum(x>=new["log_rail_km"][t]-1e-12 for x in v))/(len(v)+1) if len(v)>=.95*DRAWS else None}
        for t,v in vals.items()}
    return out
def run(cells,road,tree,pop):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:raise ValueError("Wrong original cells")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong roads")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree")
    if pop["version"]!="v391-original-JRC-GHS-POP-2020-1km-global-counts" or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong population original")
    panels=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];rr=source_rows(p)
        rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
        for land in LAND_THRESHOLDS:
            key=w,pad,land
            rows,_=enrich_cells(rr,rd,tr,w,land)
            m,qa=add_population(rows,pp,rd)
            n,k,_=EXPECTED[key]
            if len(rows)!=n or len(m)!=n:raise ValueError("v392 population support changed")
            out={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":land,"population_support_excluded":qa,"fits":[]}
            for i,mode in enumerate(("all","high_effort")):
                r=eligible_region_rows(m if i==0 else [x for x in m if high(x)])
                prior=fit(r,FULL,False);prior_d=POP_EXPECT[key][i]
                if prior_d is not None and (prior["coefficients"] is None or abs(prior["coefficients"]["log_rail_km"]["log1p_documented"]-prior_d)>1e-9):
                    raise ValueError("v392 original population-adjusted coefficient mismatch")
                if prior_d is None and prior["coefficients"] is not None:raise ValueError("Unexpected prior support")
                out["fits"].append({"support":mode,"joint_test":compare(r,LOCAL_EXPECT[key][i],SEED+w+int(pad*100)+int(land*10)+i)})
            panels.append(out);print(json.dumps({"v394_panel":out},allow_nan=False),flush=True)
    return {"version":"v394-original-joint-population-region-block-naturalized-richness",
            "population_SHA256":POP_SHA,"panels":panels,"bootstrap_draws":BOOT,"swap_draws":DRAWS,
            "limitations":["posthoc joint sensitivity, same observational records not independent replication",
            "residential population 2020 and road/tree historical, not traffic",
            "region-block FE does not remove all within-stratum confounding or detection bias",
            "technical swap ranks not causal p-values, exposure geometry not preserved",
            "unlisted is not native; cross-sectional observed richness not actual temporal increase"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","population","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.cells),read(a.road),read(a.tree),read(a.population));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()

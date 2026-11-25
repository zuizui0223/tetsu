#!/usr/bin/env python3
"""Environmental-frame audit of selection into citizen-science diversity analyses.
Zero counts mean no eligible records, not species absence. No raw scan."""
import argparse,json,math,hashlib
from pathlib import Path
from collections import Counter,defaultdict
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import GEO,extract,linear,interval,combine
from analyze_world_railstatus_observer_split_connectivity import adapt_pair,split_pair,build_pairs
from analyze_world_railstatus_plant_insect_bird_beta import attach_agb,join,EXPECTED
from scan_world_railstatus_observer_split_beta import MIN_FULL,MIN_SPLIT,SALT
from analyze_world_rail_status_richness_railstatus import CATS,LOG
from analyze_world_railstatus_cross_taxa import GROUPS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import continent_holdout
PB=("plants","birds")
FIELDS=GEO
OUTCOMES=("plant_record","bird_record","both_record","both_full","both_split_endpoint")
def make_frame(clim,rails):
    m={}
    for z in rails:
        k=b.key(z);v=float(z["source_real_original_Overture_normal_rail_km_clipped_into_cell"])
        if k in m or not math.isfinite(v) or v<=0:raise ValueError("Bad original rail census")
        m[k]=v
    return [dict(key=k,continent=z["continent"],block=z["block"],bio1=z["bio1"],bio1_sq=z["bio1"]**2,
                 abs_lat=z["abs_lat"],abs_lat_sq=z["abs_lat"]**2,rail_km=m.get(k,0.),log_rail_km=math.log1p(m.get(k,0.))) for k,z in sorted(clim.items())]
def classify(r):
    if r["active_normal_km"]>0:return "active_positive"
    if all(r[c+"_km"]==0 for c in CATS):return "all_four_zero"
    return "other_only"
def outcomes(row,counts,splits):
    k=row["key"];q=dict(row,rail_class=classify(row))
    for g in PB:
        n=counts[g].get(k,(0,0,0))
        if len(n)!=3 or min(n)<0 or n[1]>n[0] or n[2]>n[0]:raise ValueError("Invalid source counts")
        if k in splits[g] and splits[g][k][0]!=n:raise ValueError("Split count replay failed")
        q[g+"_record"]=float(n[0]>0)
        q[g+"_full"]=float(all(x>=y for x,y in zip(n,MIN_FULL)))
        q[g+"_split_endpoint"]=float(k in splits[g])
        if q[g+"_split_endpoint"]>q[g+"_full"]:raise ValueError("Non-nested support")
    q["plant_record"]=q["plants_record"];q["bird_record"]=q["birds_record"]
    for label,suffix in (("both_record","record"),("both_full","full"),("both_split_endpoint","split_endpoint")):
        q[label]=q["plants_"+suffix]*q["birds_"+suffix]
    return q
def funnel(rows):
    def report(rr):
        return {"n":len(rr),"counts":{k:int(sum(r[k] for r in rr)) for k in OUTCOMES},
                "fractions":{k:sum(r[k] for r in rr)/len(rr) for k in OUTCOMES} if rr else {}}
    return {"all":report(rows),"by_class":{c:report([r for r in rows if r["rail_class"]==c]) for c in ("active_positive","all_four_zero","other_only")},
        "by_continent":{c:{cl:report([r for r in rows if r["continent"]==c and r["rail_class"]==cl]) for cl in ("active_positive","all_four_zero","other_only")} for c in sorted({r["continent"] for r in rows})}}
def fit(rows,adjusted=True):
    rr=b.local_support(rows)
    fields=FIELDS if adjusted else tuple(LOG[c] for c in CATS)
    if len(set(fields))!=len(fields) or not {b.ACTIVE,b.INACTIVE}<=set(fields):raise ValueError("Invalid fields")
    out={}
    for target in OUTCOMES:
        f={t:linear(rr,fields,target,t) for t in (b.ACTIVE,b.INACTIVE)}
        if any(v is None for v in f.values()):out[target]={"status":"insufficient support or rank"};continue
        out[target]={"status":"estimated","active":interval(f[b.ACTIVE]),"inactive":interval(f[b.INACTIVE]),
          "active_minus_inactive":interval(combine(f,{b.ACTIVE:1.,b.INACTIVE:-1.}))}
    return {"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),"effects":out}
def native(rows):
    active=[r for r in rows if r["rail_class"]=="active_positive"]
    zero=[r for r in rows if r["rail_class"]=="all_four_zero"]
    ca=Counter((r["continent"],r["block"]) for r in active);cz=Counter((r["continent"],r["block"]) for r in zero)
    common=set(ca)&set(cz)
    # Equal block weight for each category; paired descriptive fraction difference.
    diffs={}
    for k in OUTCOMES:
        a=defaultdict(list);z=defaultdict(list)
        for r in active:a[r["continent"],r["block"]].append(r[k])
        for r in zero:z[r["continent"],r["block"]].append(r[k])
        vals=[sum(a[g])/len(a[g])-sum(z[g])/len(z[g]) for g in sorted(common)]
        diffs[k]={"mean_active_minus_zero":sum(vals)/len(vals) if vals else None,"positive_blocks":sum(v>0 for v in vals),"negative_blocks":sum(v<0 for v in vals)}
    return {"n_common_blocks":len(common),"n_active_cells":sum(ca[k] for k in common),"n_zero_cells":sum(cz[k] for k in common),"differences":diffs}
def main():
    p=argparse.ArgumentParser()
    for k in ("split","prior","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    src,prior,edge,agb,am,plant,rail,status,tree,road,pop=map(b.read,(a.split,a.prior,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
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
                   "funnel":funnel(rows),"same_block_raw":native(rows),"conditional":fit(rows),
                   "unadjusted_local":fit(rows,False),"AGB_supported_funnel":funnel(sub)}
                panels.append(z);print(json.dumps({"observation_support_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-observation-support-audit-v1","source_replay":checks,"source_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")},"panels":panels,
        "primary":"100km land80 independent environmental land-cell frame; conditional local block linear probability contrast active minus inactive mapped rail log-km coefficients for both plant and bird split-supported endpoints.",
        "limits":["no record is not species absence","eligible split endpoints depend on richness and observer effort and adjacency, not just observer participation","environmental complete-case frame is not all land","map status is not measured traffic, no causal observation bias or biodiversity effect","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

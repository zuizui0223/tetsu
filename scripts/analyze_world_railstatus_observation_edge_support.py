#!/usr/bin/env python3
"""Independent land-frame adjacent-pair observation support and mapped rail edges.
No occurrence scan, no ecological absences inferred. Local and endpoint FE
are descriptive nonrandom network contrasts, not effects of actual travel.
"""
import argparse,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import GEO,extract,combine,interval
from analyze_world_railstatus_observation_support import make_frame,attach_agb
from analyze_world_railstatus_plant_insect_bird_beta import score_fit,EXPECTED
from analyze_world_railstatus_observation_stages import outcomes,verify_replay,fit as stage_fit,funnel
from analyze_world_railstatus_coupling_climate import climate_extra,EXTRA as CLIMATE
from analyze_world_railstatus_observer_split_connectivity import CATS
from scan_world_railstatus_observer_split_beta import SALT,MIN_FULL,MIN_SPLIT
from analyze_world_railstatus_cross_taxa import GROUPS
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import fit_predict,metrics,CONTINENTS
PB=("plants","birds")
CELL=GEO+(b.BIO,)+CLIMATE
CROSS=tuple("cross_"+c for c in CATS)
ACTIVE="cross_active_normal";INACTIVE="cross_inactive_normal"
MEAN=tuple("mean_"+k for k in CELL)
DIFF=tuple("diff_"+k for k in CELL)
FIELDS=MEAN+DIFF+("horizontal",)+CROSS
ENDPOINT=DIFF+("horizontal",)+CROSS
OUTCOMES=("both_record","both_effort","both_full","both_split_pairs","plants_effort","birds_effort")

def pairkey(p):
    a,b=tuple(p["a"]),tuple(p["b"])
    if b not in ((a[0]+1,a[1]),(a[0],a[1]+1)):raise ValueError("Nonadjacent source pair")
    return a,b

def edge_map(edge,w):
    out={}
    for c in CATS:
        d={}
        for z in edge["edges"][str(w)][c]:
            k=pairkey(z);n=z["crossing_source_segments"]
            if k in d or n<=0 or int(n)!=n:raise ValueError("Invalid mapped edge")
            d[k]=int(n)
        out[c]=d
    return out

def split_maps(sp):
    out={}
    for g in GROUPS:
        # Verify source pair structure and split endpoint counts before outcome use.
        ex=extract(sp[g]);m={}
        for z in sp[g]:
            k=pairkey(z)
            if k in m:raise ValueError("Repeated original pair")
            m[k]=("split_similarity" in z)
            if m[k] and not set(k)<=set(ex):raise ValueError("Invalid split endpoint")
        out[g]={k for k,v in m.items() if v}
    return out

def make_pairs(rows,edges,split):
    by={r["key"]:r for r in rows}
    if len(by)!=len(rows):raise ValueError("Repeated environmental cell")
    out=[]
    for a,x in sorted(by.items()):
        for k in ((a[0]+1,a[1]),(a[0],a[1]+1)):
            y=by.get(k)
            if y is None or (x["continent"],x["block"])!=(y["continent"],y["block"]):continue
            r={"a":a,"b":k,"continent":x["continent"],"block":x["block"],"horizontal":float(k[0]!=a[0]),
               "both_active":min(x["active_normal_km"],y["active_normal_km"])>0,
               "both_all_four_zero":all(t[c+"_km"]==0 for t in (x,y) for c in CATS)}
            for f in CELL:
                r["mean_"+f]=(x[f]+y[f])/2;r["diff_"+f]=abs(x[f]-y[f])
            for c in CATS:
                n=edges[c].get((a,k),0)
                if n>0 and min(x[c+"_km"],y[c+"_km"])<=0:raise ValueError("Edge and endpoint km disagree")
                r["cross_"+c]=float(n>0)
            for lab in ("both_record","both_effort","both_full"):
                r[lab]=x[lab]*y[lab]
            for g,name in (("plants","plant"),("birds","bird")):
                r[g+"_effort"]=x[name+"_effort"]*y[name+"_effort"]
            for g in GROUPS:r[g+"_split_pair"]=float((a,k) in split[g])
            r["both_split_pairs"]=r["plants_split_pair"]*r["birds_split_pair"]
            if not r["both_record"]>=r["both_effort"]>=r["both_full"]>=r["both_split_pairs"]:raise ValueError("Non-nested pair support")
            out.append(r)
    return out

def describe(rows):
    groups={"all":rows,"active_crossing":[r for r in rows if r[ACTIVE]],
       "no_active_crossing":[r for r in rows if not r[ACTIVE]],
       "inactive_crossing":[r for r in rows if r[INACTIVE]],
       "no_inactive_crossing":[r for r in rows if not r[INACTIVE]],
       "both_all_four_zero":[r for r in rows if r["both_all_four_zero"]]}
    return {k:{"n":len(v),"n_blocks":len({(r["continent"],r["block"]) for r in v}),
       "rates":{t:sum(r[t] for r in v)/len(v) for t in OUTCOMES} if v else {}} for k,v in groups.items()}

def fit(rows,endpoint=False):
    fields=ENDPOINT if endpoint else FIELDS
    qa={"n_pairs":len(rows),"n_blocks":len({(r["continent"],r["block"]) for r in rows}),
        "n_endpoints":len({r[k] for r in rows for k in ("a","b")})}
    effects={};fits={}
    for t in OUTCOMES:
        f=score_fit(rows,fields,t,(ACTIVE,INACTIVE),endpoint)
        if f is None:effects[t]={"status":"insufficient support or rank"};continue
        z={"active":f[ACTIVE],"inactive":f[INACTIVE],"active_minus_inactive":combine(f,{ACTIVE:1.,INACTIVE:-1.})}
        fits[t]=z;effects[t]={"status":"estimated",**{k:interval(v) for k,v in z.items()}}
    if all(t in fits for t in ("plants_effort","birds_effort")):
        effects["plant_minus_bird_effort"]={"status":"estimated",**{k:interval(combine({g:fits[g][k] for g in ("plants_effort","birds_effort")},{"plants_effort":1.,"birds_effort":-1.})) for k in fits["plants_effort"]}}
    return dict(qa,effects=effects)

def predict(rows):
    base=tuple(k for k in FIELDS if k not in (ACTIVE,INACTIVE))
    out={}
    for t in ("both_effort","both_split_pairs"):
        folds=[]
        for c in CONTINENTS:
            tr=[r for r in rows if r["continent"]!=c];te=[r for r in rows if r["continent"]==c]
            if len(te)<100 or len(tr)<max(300,4*len(FIELDS)):continue
            specs={"both":FIELDS,"base":base,"no_active":tuple(k for k in FIELDS if k!=ACTIVE),"no_inactive":tuple(k for k in FIELDS if k!=INACTIVE)}
            m={k:metrics(te,fit_predict(tr,te,f,t,True),t)[0] for k,f in specs.items()}
            folds.append({"continent":c,"n_test":len(te),"Brier":m,"gain_both":m["base"]-m["both"],
                "gain_active":m["no_active"]-m["both"],"gain_inactive":m["no_inactive"]-m["both"]})
        terms=("gain_both","gain_active","gain_inactive")
        out[t]={"n_folds":len(folds),"adequate_folds":len(folds)>=4,"folds":folds,
                "mean_gain":{k:sum(r[k] for r in folds)/len(folds) for k in terms} if folds else {},
                "positive_folds":{k:sum(r[k]>0 for r in folds) for k in terms}}
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ("prior","split","stage_model","edges","agb","agb_model","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    prior,src,st,edge,agb,am,plant,rail,status,tree,road,pop=map(b.read,(a.prior,a.split,a.stage_model,a.edges,a.agb,a.agb_model,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if st["version"]!="railstatus-observation-stage-decomposition-v1" or any(st["source_SHA256"][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("split","prior","road","population","status")):raise ValueError("Wrong stage source")
    old={(p["width_m"],p["land"]):p for p in st["panels"]}
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Wrong split source")
    if any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")) or any(src["qa"][g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Frozen source mismatch")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in CATS):raise ValueError("Wrong mapped geometry")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or am["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or am["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Wrong biomass source")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9" or b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong environmental sources")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in b.GRID_WIDTHS:
            grid=prior["grids"][str(w)];sp=src["pairs"][str(w)]
            counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in PB}
            splits={g:extract(sp[g]) for g in PB};ps=split_maps(sp);ed=edge_map(edge,w)
            for g in PB:
                if any(counts[g].get(k)!=v[0] for k,v in splits[g].items()):raise ValueError("Split counts changed")
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            clim,_=b.source_climate_lookup(land,ras,scale,set(rd),w);ext,meta=climate_extra(a.climate,ras,clim,w)
            frame=make_frame(clim,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)])
            frame,_=b.inject_built(frame,built[str(w)]);frame,_=b.inject_status(frame,status,w)
            for gate in b.LAND_THRESHOLDS:
                rr,_=b.enrich_cells(frame,rd,tr,w,gate);rr,_=b.join_population(rr,pp,rd)
                r=[outcomes(x,counts,splits) for x in rr];f=funnel(r)
                if f["all"]["counts"]!=old[w,gate]["funnel"]["all"]["counts"] or f["all"]["n"]!=old[w,gate]["funnel"]["all"]["n"]:raise ValueError("Previous stage denominator changed")
                ag=attach_agb(rr,ag20,ag24,rd)
                # Check prior active-only pair census before applying added AGB/climate support.
                keys={x["key"]:x for x in rr}
                def retained(k):
                    x,y=(keys.get(q) for q in k)
                    return x is not None and y is not None and (x["continent"],x["block"])==(y["continent"],y["block"]) and min(x["active_normal_km"],y["active_normal_km"])>0
                replay=tuple(sum(retained(k) for k in ps[g]) for g in GROUPS)
                if replay!=EXPECTED[w,gate]:raise ValueError("Previous active pair counts changed")
                cells=[dict(x,**ext[x["key"]],**{b.BIO:ag[x["key"]][2024]}) for x in r if x["key"] in ag and x["key"] in ext]
                pairs=make_pairs(cells,ed,ps);active=[x for x in pairs if x["both_active"]]
                z={"width_m":w,"land":gate,"n_environmental_cells":len(cells),"original_active_pair_replay":dict(zip(GROUPS,replay)),
                   "all_pairs":describe(pairs),"active_endpoint_pairs":describe(active),
                   "primary_active_endpoints_block":fit(active),"active_endpoints_endpoint_FE":fit(active,True),
                   "all_pairs_block":fit(pairs),"all_pairs_endpoint_FE":fit(pairs,True),
                   "heldout_active_endpoints":predict(active),
                   "leave_continent_out":{c:fit([x for x in active if x["continent"]!=c]) for c in sorted({x["continent"] for x in active})}}
                panels.append(z);print(json.dumps({"observation_edge_support_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-observation-edge-support-v1","source_replay":checks,"source_SHA256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ("stage_model","split","prior","edges","status","agb","climate")},
        "primary":"100km land80, adjacent same-block pairs with positive active mapped km in both endpoints, direct active-minus-inactive actual crossing coefficient for both taxa effort sufficient at both endpoints, block FE plus fixed geography/mean/difference exposure and climate.",
        "panels":panels,"limits":["effort eligibility is recorded coverage, not true species absence or travel","crossings are mapped geometry not train movements","endpoint FE absorbs additive marginal endpoint tendencies only; binary joint selection is nonlinear","same-source descriptive conditional contrasts and shared continent folds","no causal rail or historic mechanism, no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

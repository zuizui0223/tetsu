"""Global mapped railway role: does status add functional information beyond rail amount?

Reuses exact 2025 aggregate supports; no new species scan, no closure history,
no forest-stratum search. HWI composition is separate from the existing HWI
connectivity evidence. No inference of causal operation/management effects.
"""
import argparse, collections, hashlib, json, math
from pathlib import Path
import numpy as np
import analyze_world_railstatus_disturbance_traits as d
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    continent_holdout, fit_predict, metrics, CONTINENTS)
b=d.b
AXES={"plants":("height_nonwoody",),"insects":("edge_affinity",),
      "birds":("HWI","grassland_vs_forest")}
PRIMARY={"height_nonwoody":200000,"edge_affinity":100000,"HWI":100000}
ROAD="log_road_density"
MIN_FOLDS=4
def source_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def fields(flex=True):
    base=d.GEO+d.CORE+d.CLIMATE+(d.EXTRA if flex else ())+(b.BIO,)
    if len(set(base))!=len(base):raise ValueError("duplicate frozen fields")
    env=tuple(k for k in base if k not in (b.ACTIVE,b.INACTIVE))
    if ROAD not in env:raise ValueError("road comparator missing")
    return {"environment":env,"total_rail":env+(d.TOTAL,),
        "total_plus_status":env+(d.TOTAL,d.SHARE),
        "separate_status_lengths":env+(b.ACTIVE,b.INACTIVE),
        "total_plus_status_without_road":tuple(k for k in env if k!=ROAD)+(d.TOTAL,d.SHARE)}
def heldout(rows,flex=True):
    """All models share exactly the same train/test rows and thresholds."""
    rr=b.with_epoch(b.local_support([r for r in rows if r["active_normal_km"]+r["inactive_normal_km"]>0]),2024)
    mods=fields(flex);folds=[]
    min_train=max(300,4*max(map(len,mods.values())))
    for continent in CONTINENTS:
        test=[r for r in rr if r["continent"]==continent]
        train=[r for r in rr if r["continent"]!=continent]
        if len(test)<50 or len(train)<min_train:continue
        loss={}
        for label,p in mods.items():
            pred=fit_predict(train,test,p,"trait_value",False)
            mse,mae=metrics(test,pred,"trait_value")
            loss[label]=float(mse)
        gain={"total_rail_beyond_environment":loss["environment"]-loss["total_rail"],
              "status_beyond_total":loss["total_rail"]-loss["total_plus_status"],
              "separate_lengths_beyond_total":loss["total_rail"]-loss["separate_status_lengths"],
              "road_beyond_total_and_status":loss["total_plus_status_without_road"]-loss["total_plus_status"]}
        folds.append({"continent":continent,"n_test":len(test),"n_train":len(train),
                      "loss":loss,"gain":gain})
    out={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),
         "n_folds":len(folds),"min_training":min_train,"min_test":50,
         "status":"estimated" if len(folds)>=MIN_FOLDS else "insufficient continental support",
         "folds":folds}
    if folds:
        comps={"total_rail_beyond_environment":"environment","status_beyond_total":"total_rail",
               "separate_lengths_beyond_total":"total_rail","road_beyond_total_and_status":"total_plus_status_without_road"}
        out["summary"]={}
        for k,base in comps.items():
            g=np.array([r["gain"][k] for r in folds])
            losses=np.array([r["loss"][base] for r in folds])
            out["summary"][k]={"mean_delta_MSE":float(g.mean()),
                "relative_mean_MSE_gain":float(g.mean()/losses.mean()) if losses.mean()>0 else None,
                "positive_folds":int(sum(g>0)),"n_folds":len(folds)}
    return out

def replay(rows,previous):
    for label,flex,share in (("primary",True,False),("linear_effort",False,False),("fixed_total_active_share",True,True)):
        new=d.estimate(rows,flex,share)[0];old=previous[label]
        if new["n_cells"]!=old["n_cells"] or new["status"]!=old["status"]:
            raise ValueError("Native mean support changed: "+label)
        if new["status"]=="estimated":
            for k,v in new["effects"].items():
                q=old["effects"][k]
                if abs(v["estimate"]-q["estimate"])>1e-9 or max(abs(x-y) for x,y in zip(v["interval95_approx"],q["interval95_approx"]))>1e-9:
                    raise ValueError("Native mean effect changed: "+label)

def contrast(rows):
    return {"separate_lengths":d.estimate(rows)[0],"fixed_total_share":d.estimate(rows,share=True)[0]}
def evaluate(ready,tm,axis,gate,prior):
    rr=d.attach(ready,tm,axis,gate)
    replay(rr,prior)
    both=[r for r in rr if min(r["active_normal_km"],r["inactive_normal_km"])>0]
    fam=d.attach(ready,tm,axis,gate,family=True)
    h={i:d.attach(ready,tm,axis,gate,side=i) for i in (0,1)}
    ks=set(r["key"] for r in h[0])&set(r["key"] for r in h[1])
    hs={i:d.estimate([r for r in h[i] if r["key"] in ks],share=True) for i in h}
    half={"native":[hs[i][0] for i in (0,1)]}
    if all(hs[i][1] is not None for i in hs):
        half["paired_share"]=d.interval(d.combine({i:hs[i][1]["active_share_at_fixed_total"] for i in hs},{0:.5,1:.5}))
    return {"native_exact_replay":True,"conditional":contrast(rr),
            "both_status_positive":contrast(both),
            "family_centered":contrast(fam),"observer_half_fixed_total":half,
            "heldout_native_positive_total":heldout(rr),
            "heldout_both_status_positive":heldout(both)},rr

def selftest():
    z=fields()
    assert b.ACTIVE not in z["total_rail"] and b.INACTIVE not in z["total_rail"]
    assert set(z["total_plus_status"])-set(z["total_rail"])=={d.SHARE}
    assert set(z["total_plus_status"])-set(z["total_plus_status_without_road"])=={ROAD}
    assert len(set(z["total_plus_status"]))==len(z["total_plus_status"])
    assert heldout([])["status"]=="insufficient continental support"

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true")
    for k in ("size","disturbance","size_model","disturbance_model","taxa","split","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+k,type=Path)
    a=p.parse_args();selftest()
    if a.self_test:return
    size,src,sm,dm,prior,split,agb,plant,rail,status,tree,road,pop=map(b.read,(
        a.size,a.disturbance,a.size_model,a.disturbance_model,a.taxa,a.split,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(prior,plant)
    if src["version"]!="railstatus-disturbance-trait-cells-v1" or not src["exact_full_and_split_replay"]:
        raise ValueError("Wrong disturbance source")
    if size["version"]!="railstatus-verified-scaled-plant-size-cells-v1" or not size["exact_full_split_target_and_life_replay"]:
        raise ValueError("Wrong plant size source")
    if src["source"]["source_SHA256"]!={"LifeForm":d.LIFE_SHA,"LepTraits":d.LEP_SHA,"AVONET":d.SOURCE_SHA,"size_archive":d.SIZE_SHA} or size["source_SHA256"]!={"LifeForm":d.LIFE_SHA,"size_archive":d.SIZE_SHA}:
        raise ValueError("Frozen trait source changed")
    if any(len({src[k],size[k],split[k]})!=1 for k in ("obs_sha256","taxa_sha256","split_seed")):
        raise ValueError("Observer partition mismatch")
    if any(prior[k]!=src[k] for k in ("obs_sha256","taxa_sha256")):
        raise ValueError("Frozen occurrence source mismatch")
    if size["prior_cells_sha256"]!=source_hash(a.disturbance) or dm["trait_cells_sha256"]!=source_hash(a.disturbance) or sm["trait_cells_sha256"]!=source_hash(a.size):
        raise ValueError("Prior analysis input digests differ")
    if dm["version"]!="railstatus-disturbance-trait-association-v1" or sm["version"]!="railstatus-verified-scaled-plant-size-association-v1":
        raise ValueError("Prior model versions differ")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1":
        raise ValueError("Wrong biomass source")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":
        raise ValueError("Frozen road/tree source differs")
    if b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:
        raise ValueError("Frozen geography mismatch")
    old={(x["width"],x["land"],x["coverage"],x["taxon"]):x for x in dm["panels"]}
    sz={(x["width"],x["land"],x["coverage"],x["taxon"]):x for x in sm["panels"]}
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);ras,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(ras)
        for w in d.WIDTHS:
            grid=prior["grids"][str(w)]
            keys={b.key(r) for g in d.GROUPS for r in grid[g]}
            clim,_=b.source_climate_lookup(land,ras,scale,keys,w)
            ext,_=d.climate_extra(a.climate,ras,clim,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
            rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)])
            pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)])
            ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            for g,axes in AXES.items():
                tm={tuple(r["key"]):r for r in (size if g=="plants" else src)["grids"][str(w)][g]}
                op=d.extract(split["pairs"][str(w)][g])
                if {k:(tuple(v["full"]),tuple(tuple(x) for x in v["halves"])) for k,v in tm.items()}!=op:
                    raise ValueError("Exact split counts changed")
                rows,_,_=b.original_cell_rows([b.convert(v) for v in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w)
                for landgate in (.8,.5):
                    rr,_=b.enrich_cells(rows,rd,tr,w,landgate)
                    rr,_=b.join_population(rr,pp,rd)
                    ready=[]
                    for r in rr:
                        k=r["key"]
                        if k not in tm or k not in ext or k not in ag20 or k not in ag24:continue
                        if tuple(tm[k]["full"])!=(r["observations"],r["observers"],r["species"]):
                            raise ValueError("Exact full count changed")
                        av,cv=ag20[k],ag24[k]
                        if min(av["valid_area_m2"],cv["valid_area_m2"])<.9*rd[k]["source_GRIP4_total_land_area_m2"]:continue
                        if any(not math.isfinite(z["agb_Mg_ha"]) or z["agb_Mg_ha"]<0 for z in (av,cv)):
                            raise ValueError("Invalid biomass")
                        ready.append(dict(r,agb2020=av["agb_Mg_ha"],agb2024=cv["agb_Mg_ha"],**ext[k]))
                    for gate in (.5,.25):
                        oldp=(sz if g=="plants" else old)[w,landgate,gate,g]
                        if len(ready)!=oldp["n_base"]:raise ValueError("Prior geographic support changed")
                        for axis in axes:
                            z,_=evaluate(ready,tm,axis,gate,oldp["traits"][axis])
                            z.update(width=w,land=landgate,coverage=gate,taxon=g,trait=axis,
                                    primary=PRIMARY.get(axis)==w and landgate==.8 and gate==.5)
                            panels.append(z)
                            # Compact results in the job log, full fits in artifact.
                            print(json.dumps({"global_rail_role_panel":{
                                "width":w,"land":landgate,"coverage":gate,"trait":axis,"primary":z["primary"],
                                "conditional":z["conditional"],"both_status_positive":z["both_status_positive"],
                                "family_centered":z["family_centered"],"observer_half_fixed_total":z["observer_half_fixed_total"].get("paired_share"),
                                "holdout":{k:{x:v for x,v in z[k].items() if x!="folds"} for k in ("heldout_native_positive_total","heldout_both_status_positive")}
                            }},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"global-rail-role-status-vs-total-v1","source_replay":checks,"panels":panels,
         "primary":"Height 200km, butterfly edge and HWI 100km; land80 coverage50. Compare total rail plus active share vs total rail on same positive-rail cells, heldout continents, original ridge and block weights. Separately report direct native and fixed total associations. Both-status-positive and all other panels sensitivities.",
         "limits":["Exploratory follow-up of existing trait signals, same 2025 iNaturalist occurrence data and static 2026 map",
          "Fixed total active share is a conditional contrast, not actual line closure or a causal equal-length replacement",
          "Different trait scales and species pools; cannot rank coefficients or merge them into a single effect",
          "Observed landscape trait composition is not trackside abundance, individual adaptation, actual dispersal or measured management",
          "HWI composition is distinct from historical low-versus-high HWI boundary similarity; previous conditional HWI mechanism unsupported",
          "Continental folds share training data and selected source coverage; prediction is not causal inference; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

"""Global common-cell joint trait signature of active vs inactive mapped rail.

Prespecified exploratory directional synthesis, not a new independent sample.
Reuse exact verified aggregate geometry and models; no raw scan.
"""
import argparse, collections, hashlib, json, math
from pathlib import Path
import numpy as np
from analyze_global_rail_role_status_vs_total import d,b,source_hash,replay
AXES={"plants":("height_nonwoody",),"insects":("edge_affinity",),"birds":("HWI",)}
ORDER=("height_nonwoody","edge_affinity","HWI")
SIGN={"height_nonwoody":-1.,"edge_affinity":1.,"HWI":-1.}
EFFORT=d.CORE+d.EXTRA
ROAD="log_road_density"
GEO=d.GEO+d.CLIMATE
if len(set(GEO))!=len(GEO):raise ValueError("duplicate geography fields")
ENV=tuple(k for k in GEO if k not in (b.ACTIVE,b.INACTIVE))+(b.BIO,)
JOINT=tuple(f"{a}__{k}" for a in ORDER for k in EFFORT)

def join(maps,mode="full",keys=None):
    by={a:{r["key"]:r for r in maps[a][mode]} for a in ORDER}
    common=set.intersection(*(set(v) for v in by.values()))
    if keys is not None:common &= keys
    out=[]
    for k in sorted(common):
        p=by[ORDER[0]][k]
        if p["active_normal_km"]+p["inactive_normal_km"]<=0:continue
        r=dict(p)
        for a in ORDER:
            q=by[a][k]
            for f in GEO+("continent","block","agb2020","agb2024",d.TOTAL,d.SHARE,
                          "active_normal_km","inactive_normal_km"):
                if p[f]!=q[f]:raise ValueError("Joint geographic source differs: "+f)
            r[a]=q["trait_value"]
            for f in EFFORT:r[f"{a}__{f}"]=q[f]
        out.append(r)
    return b.local_support(out)

def scales(rows):
    if not rows:return None
    w=b.block_weights(rows)
    ans={}
    for a in ORDER:
        x=np.asarray([r[a] for r in rows])
        m=float(np.average(x,weights=w))
        sd=float(np.sqrt(np.average((x-m)**2,weights=w)))
        if not math.isfinite(sd) or sd<1e-10:return None
        ans[a]={"mean":m,"sd":sd}
    return ans

def effect(rows,sc,share=True,flex=True,epoch=2024):
    rr=b.with_epoch(b.local_support(rows),epoch)
    groups=collections.defaultdict(list)
    for r in rr:
        groups[r["continent"],r["block"]].append(r)
    f=JOINT if flex else tuple(f"{a}__{k}" for a in ORDER for k in d.CORE)
    if share:
        terms=(d.TOTAL,d.SHARE)
    else:terms=(b.ACTIVE,b.INACTIVE)
    fields=ENV+f+terms
    if len(set(fields))!=len(fields):raise ValueError("duplicate joint model fields")
    var={k:sum(len({r[k] for r in v})>=2 for v in groups.values()) for k in terms}
    qa={"n_cells":len(rr),"n_blocks":len(groups),"n_varying_blocks":var,
        "continents":dict(collections.Counter(r["continent"] for r in rr))}
    if sc is None or len(rr)<100 or len(groups)<20 or min(var.values(),default=0)<20 or len(rr)-len(groups)-len(fields)<30:
        return dict(qa,status="insufficient common support"),None
    # Outcomes are standardized with the same primary full-support scale,
    # including all subsets; no data-dependent direction selection.
    for r in rr:
        for a in ORDER:r["z_"+a]=(r[a]-sc[a]["mean"])/sc[a]["sd"]
    allfits={}
    for name in (*terms,ROAD):
        ff={a:d.linear(rr,fields,"z_"+a,name) for a in ORDER}
        if any(x is None for x in ff.values()):
            return dict(qa,status="insufficient rank or target variation",failed_term=name),None
        allfits[name]=ff
    if share:
        target=d.SHARE
    else:
        target="active_minus_inactive"
        allfits[target]={a:d.combine({0:allfits[b.ACTIVE][a],1:allfits[b.INACTIVE][a]},{0:1.,1:-1.}) for a in ORDER}
    def describe(ff):
        oriented={a:d.combine({a:ff[a]},{a:SIGN[a]}) for a in ORDER}
        avg=d.combine(oriented,{a:1/3 for a in ORDER})
        original={a:d.interval(d.combine({a:ff[a]},{a:sc[a]["sd"]})) for a in ORDER}
        oo={a:d.interval(x) for a,x in oriented.items()}
        return {"raw_effects":original,"oriented_standardized":oo,
                "directional_mean_standardized":d.interval(avg),
                "all_three_point_directions":all(v["estimate"]>0 for v in oo.values()),
                "all_three_95_intervals_positive":all(v["interval95_approx"][0]>0 for v in oo.values())},oriented
    effects={k:describe(v)[0] for k,v in allfits.items()}
    return dict(qa,status="estimated",target=target,effects=effects),allfits[target]

def both(rows):return [r for r in rows if min(r["active_normal_km"],r["inactive_normal_km"])>0]

def panel(maps):
    rr=join(maps);sc=scales(rr)
    keys={r["key"] for r in rr}
    fam=join(maps,"family",keys)
    h0,h1=join(maps,"half0",keys),join(maps,"half1",keys)
    hk={r["key"] for r in h0}&{r["key"] for r in h1}
    hh=[b.local_support([r for r in h if r["key"] in hk]) for h in (h0,h1)]
    if {r["key"] for r in hh[0]}!={r["key"] for r in hh[1]}:
        raise ValueError("Half support differs")
    hp=[effect(h,sc) for h in hh]
    pair=None
    if all(x[1] is not None for x in hp):
        ff={a:d.combine({0:hp[0][1][a],1:hp[1][1][a]},{0:.5,1:.5}) for a in ORDER}
        oo={a:d.combine({a:ff[a]},{a:SIGN[a]}) for a in ORDER}
        pair={"oriented_standardized":{a:d.interval(x) for a,x in oo.items()},
              "directional_mean_standardized":d.interval(d.combine(oo,{a:1/3 for a in ORDER})),
              "all_three_95_intervals_positive":all(d.interval(x)["interval95_approx"][0]>0 for x in oo.values())}
    return {"native_n_by_trait":{a:len(maps[a]["full"]) for a in ORDER},
            "common_n":len(rr),"standardization":sc,"prior_native_replay":True,
            "fixed_total_primary":effect(rr,sc)[0],"separate_lengths_same_support":effect(rr,sc,False)[0],
            "both_positive_fixed_total":effect(both(rr),sc)[0],
            "both_positive_separate":effect(both(rr),sc,False)[0],
            "linear_effort":effect(rr,sc,flex=False)[0],
            "AGB2020":effect(rr,sc,epoch=2020)[0],
            "family_centered":effect(fam,sc)[0],
            "full_on_family_support":effect([r for r in rr if r["key"] in {x["key"] for x in fam}],sc)[0],
            "half_results":[p[0] for p in hp],"half_paired":pair,
            "full_on_half_support":effect([r for r in rr if r["key"] in {x["key"] for x in hh[0]}],sc)[0],
            "leave_continent_out":{c:effect([r for r in rr if r["continent"]!=c],sc)[0] for c in sorted({r["continent"] for r in rr})}}

def selftest():
    assert set(AXES)=={"plants","insects","birds"} and set(SIGN)==set(ORDER)
    assert scales([]) is None and effect([],None)[0]["status"]=="insufficient common support"
    assert len(set(ENV+JOINT+(d.TOTAL,d.SHARE)))==len(ENV+JOINT+(d.TOTAL,d.SHARE))

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
    panels=[]; common={}
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
                            native=d.attach(ready,tm,axis,gate)
                            replay(native,oldp["traits"][axis])
                            z=common.setdefault((w,landgate,gate),{})
                            z[axis]={"full":native,"family":d.attach(ready,tm,axis,gate,family=True),
                                     "half0":d.attach(ready,tm,axis,gate,side=0),
                                     "half1":d.attach(ready,tm,axis,gate,side=1)}
    finally:ras.close();mem.close()
    for (w,landgate,gate),z in sorted(common.items()):
        ans=panel(z)
        ans.update(width=w,land=landgate,coverage=gate,
                   primary=(w==100000 and landgate==.8 and gate==.5))
        panels.append(ans)
        print(json.dumps({"global_joint_role_panel":ans},allow_nan=False),flush=True)
    out={"version":"global-rail-joint-functional-signature-v1","source_replay":checks,
         "panels":panels,"direction":{"height_nonwoody":-1,"edge_affinity":1,"HWI":-1},
         "primary":"100km land80 coverage50, same three-trait supported cells with positive total rail, fixed-total active share; test individual directional associations and their conjunction. Equal-weight standardized directional average is a descriptive composite, not evidence all axes respond.",
         "limits":["Exploratory synthesis of previously discovered directions, not independent confirmation",
         "Each grid has species-equal means from different species pools, not individual trait adaptation or trackside sampling",
         "No interpretation of HWI mean as actual connectivity; separate HWI-guild similarity result retained in evidence",
         "Conditional status association, not causal operation or rail closure, all sources 2025 occurrence and mapped 2026 rail",
         "Observation effort and geographic support differ by taxon; joint support and both-positive subset can be selective; no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

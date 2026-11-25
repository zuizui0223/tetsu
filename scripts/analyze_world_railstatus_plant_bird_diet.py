#!/usr/bin/env python3
"""Trophic contrast of plant richness and woody biomass association.
Shared geography, same effort design in both guilds, observer-disjoint slopes.
"""
import argparse,json,math,hashlib
from collections import Counter
from pathlib import Path
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_bird_disjoint import (
    cells_from_pairs,candidate_support,pair_key,MIN_FULL,MIN_SPLIT,SALT,YEARS,GROUPS,
    GEO,EFF,FLEX,ORIENTATIONS,TERMS,EXPECTED,extract,attach,linear,interval,combine,one,temporal)
from scan_world_railstatus_bird_diet_cells import SOURCE_SHA,GUILDS,GATE
AXES=(("fruit_seed","invertivore"),("frugivore","invertivore"),("granivore","invertivore"))

def diet_attach(rows,records,gs):
    out=[]
    for r in rows:
        k=r["key"]
        if k not in records:raise ValueError("Missing exact replay diet cell")
        q=records[k]
        for i in (0,1):
            for name,count in zip(("observations","observers","species"),q["halves"][i]):
                if abs(r[f"birds_{i}_{name}"]-math.log1p(count))>1e-10:raise ValueError("Trophic source split replay failed")
        if not all(g in q["guilds"] for g in gs):continue
        z=dict(r)
        for j,g in enumerate(gs):
            for i,counts in enumerate(q["guilds"][g]):
                if any(v<m for v,m in zip(counts,GATE)):raise ValueError("Guild gate invalid")
                for name,count in zip(("observations","observers","species"),counts):
                    f=f"d{j}_{i}_{name}";z[f]=math.log1p(count)
                    if name!="species":z[f+"_sq"]=z[f]**2;z[f+"_cu"]=z[f]**3
                z[f"d{j}_{i}_effort_interaction"]=z[f"d{j}_{i}_observations"]*z[f"d{j}_{i}_observers"]
        out.append(z)
    return out

DEFF=tuple(f"d{j}_{i}_{k}" for j in (0,1) for i in (0,1) for k in ("observations","observers"))
DFLEX=tuple(f"{k}_{s}" for k in DEFF for s in ("sq","cu"))+tuple(f"d{j}_{i}_effort_interaction" for j in (0,1) for i in (0,1))
EFFECTS=("plant_richness","woody_AGB")

def one_diet(rows,epoch=2024,flex=False):
    rr=b.with_epoch(b.local_support(rows),epoch)
    fields=GEO+(b.BIO,)+EFF+DEFF+(FLEX+DFLEX if flex else ())
    if len(set(fields))!=len(fields):raise ValueError("Duplicate effort fields")
    qa={"n_cells":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),"continents":dict(Counter(r["continent"] for r in rr))}
    f={}
    for i,j in ORIENTATIONS:
        term=f"plants_{i}_species";pred=fields+(term,)
        # Same regressors for both guild outcomes. Direct difference
        # uses matched block scores; never compare two p-values.
        for metric,k in (("plant_richness",term),("woody_AGB",b.BIO)):
            vals={g:linear(rr,pred,f"d{g}_{j}_species",k) for g in (0,1)}
            if any(x is None for x in vals.values()):return dict(qa,status="insufficient support or rank"),None
            f[i,j,metric]={"guild0":vals[0],"guild1":vals[1],"guild0_minus_guild1":combine(vals,{0:1.,1:-1.})}
    result={}
    for metric in EFFECTS:
        result[metric]={}
        for c in ("guild0","guild1","guild0_minus_guild1"):
            z={ij:f[(*ij,metric)][c] for ij in ORIENTATIONS}
            result[metric][c]={label:combine(z,weight) for label,weight in TERMS.items()}
    return dict(qa,status="estimated",effects={m:{c:{s:interval(v) for s,v in z.items()} for c,z in q.items()} for m,q in result.items()}),result

def native(rows,records,gs,extended=False):
    rr=diet_attach(rows,records,gs);rr=b.local_support(rr)
    z={"groups":gs,"n_supported":len(rr),"primary":one_diet(rr)[0],"flexible_effort":one_diet(rr,flex=True)[0]}
    if extended:
        z["AGB2020_same_cells"]=one_diet(rr,2020)[0]
        z["high_effort"]=one_diet([r for r in rr if min(r["observations"],r["bird_observations"])>=200 and min(r["observers"],r["bird_observers"])>=20])[0]
        z["leave_continent_out"]={c:one_diet([r for r in rr if r["continent"]!=c])[0] for c in sorted({r["continent"] for r in rr})}
    return z

def annual_diet(by,records,gs):
    rr={y:diet_attach(by[y],records[y],gs) for y in YEARS}
    keys=set.intersection(*({r["key"] for r in v} for v in rr.values()))
    rr={y:b.local_support([r for r in v if r["key"] in keys]) for y,v in rr.items()}
    if any({r["key"] for r in rr[y]}!={r["key"] for r in rr["2025"]} for y in YEARS):raise ValueError("Frozen annual diet geography changed")
    out={"n_common_cells":len(rr["2025"])}
    for label,flex in (("linear",False),("flexible",True)):
        z={};f={}
        for y in YEARS:z[y],f[y]=one_diet(rr[y],flex=flex)
        q={"years":z}
        if all(f.values()):
            for name,ys in (("earlier_mean",YEARS[:2]),("all_year_mean",YEARS)):
                q[name]={m:{c:{s:interval(combine({y:f[y][m][c][s] for y in ys},{y:1/len(ys) for y in ys})) for s in TERMS} for c in f[ys[0]][m]} for m in EFFECTS}
        out[label]=q
    return out

def replay(new,old):
    if new["n_cells"]!=old["n_cells"] or new["status"]!=old["status"]:raise ValueError("Prior observer split support changed")
    if new["status"]=="estimated":
        for k,v in new["contrasts"].items():
            q=old["contrasts"][k]
            if abs(v["estimate"]-q["estimate"])>1e-9 or max(abs(a-b) for a,b in zip(v["interval95_approx"],q["interval95_approx"]))>1e-9:raise ValueError("Prior slope replay failed")

def main():
    p=argparse.ArgumentParser()
    for n in ("diet","disjoint","agb","agb_model","split","annual","taxa","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();diet,prior,agb,old,src,annual,taxa,plant,rail,status,tree,road,pop=map(b.read,(a.diet,a.disjoint,a.agb,a.agb_model,a.split,a.annual,a.taxa,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=b.verify_original_sources(plant,rail);b.verify_replay(taxa,plant)
    if diet['version']!='railstatus-avonet-trophic-disjoint-cell-counts-v1' or not diet['exact_full_and_split_cell_replay'] or diet['source']['source_SHA256']!=SOURCE_SHA or any(diet[k]!=src[k] for k in ('obs_sha256','taxa_sha256','split_seed')):raise ValueError('Wrong trophic source')
    if prior['version']!='railstatus-observer-disjoint-plant-bird-biomass-v1' or prior['biomass_source_sha256']!=agb['source_sha256'] or any(prior['source_artifact_sha256'][k]!=hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ('agb','agb_model','split','annual')):raise ValueError('Wrong previous disjoint model')
    prior_panels={(p['width_m'],p['land']):p for p in prior['panels']}
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or old["version"]!="railstatus-plant-bird-woody-biomass-association-v1" or old["biomass_source"]!={k:v for k,v in agb.items() if k!="grids"}:raise ValueError("Frozen biomass replay mismatch")
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or src["split_seed"]!=SALT.decode() or src["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or src["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Frozen observer split invalid")
    if annual["version"]!="railstatus-fixed-2025-observer-split-annual-v1" or annual["width_m"]!=200000 or annual["qa"]["rows"]!=288721938 or any(len({src[k],annual[k],taxa[k]})!=1 for k in ("obs_sha256","taxa_sha256")) or annual["split_seed"]!=src["split_seed"]:raise ValueError("Annual source mismatch")
    cand,_=candidate_support(src)
    if annual["candidate_pairs"]!={g:len(cand[g]) for g in GROUPS}:raise ValueError("Frozen annual candidates changed")
    for y in YEARS[:2]:
        for g in GROUPS:
            if not {pair_key(p) for p in annual["pairs"][y][g]}<=cand[g]:raise ValueError("Annual candidates changed")
    if tree["source_file_SHA256"]!=b.EXPECTED_SHA or road["source_zip_MD5"]!=b.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if b.sha_file(a.ne)!=b.ORIGINAL_NE_SHA or b.sha_file(a.climate)!=b.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=b.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen map/climate/pop mismatch")
    built,_=b.read_source_builtup(a.built);land,_=b.source_continents(a.ne);raster,mem,_=b.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=b.climate_scale(raster)
        for w in b.GRID_WIDTHS:
            grid=taxa["grids"][str(w)];yr=YEARS if w==200000 else ("2025",)
            ps={y:{g:src["pairs"][str(w)][g] if y=="2025" else annual["pairs"][y][g] for g in GROUPS} for y in yr}
            maps={y:{g:extract(ps[y][g]) for g in GROUPS} for y in yr}
            prior_counts={g:{b.key(x):(x["observations"],x["observers"],x["species"]) for x in grid[g]} for g in GROUPS}
            for g in GROUPS:
                if any(prior_counts[g].get(k)!=v[0] for k,v in maps["2025"][g].items()):raise ValueError("2025 exact source replay failed")
            keys={b.key(z) for g in GROUPS for z in grid[g]}|{k for y in yr for g in GROUPS for k in maps[y][g]}
            clim,_=b.source_climate_lookup(land,raster,scale,keys,w)
            tr=b.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=b.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=b.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=b.make_lookup(agb["grids"]["2020"][str(w)]);ag24=b.make_lookup(agb["grids"]["2024"][str(w)])
            geo={}
            for y in yr:
                geo[y]={}
                for g in GROUPS:
                    raw=[b.convert(z) for z in grid[g]] if y=="2025" else cells_from_pairs([p for p in ps[y][g] if "split_similarity" in p])
                    rows,_,_=b.original_cell_rows(raw,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                    rows,_=b.inject_built(rows,built[str(w)]);rows,_=b.inject_status(rows,status,w);geo[y][g]=rows
            for gate in b.LAND_THRESHOLDS:
                by={}
                for y in yr:
                    rows={}
                    for g in GROUPS:
                        m,_=b.enrich_cells(geo[y][g],rd,tr,w,gate);m,_=b.join_population(m,pp,rd)
                        rows[g]=[b.flexible(r) for r in m]
                    joined,_=b.join(rows,ag20,ag24,w,rd)
                    bc={r["key"]:r["species"] for r in rows["birds"]}
                    joined=[dict(r,bird_species=bc[r["key"]]) for r in joined]
                    if y=="2025":
                        n,expected=EXPECTED[w][gate];m=b.models(joined,2024)
                        q=m["bird_plant_richness_association"]["effects"]["plant_richness"]["bird_log_richness"]["estimate"]
                        if m["n"]!=n or abs(q-expected)>5.1e-6:raise ValueError("Prior full-count plant/bird coefficient not replayed")
                    by[y]=attach(joined,maps[y])
                oldp=prior_panels[w,gate]
                replay(one(by["2025"])[0],oldp["2025"]["primary"])
                replay(one(by["2025"],flex=True)[0],oldp["2025"]["flexible_effort"])
                records={y:{tuple(r["key"]):r for r in diet["grids"][str(w)][y]} for y in yr}
                z={"width_m":w,"land":gate,"2025":{gs[0]+"_vs_"+gs[1]:native(by["2025"],records["2025"],gs,gs==AXES[0]) for gs in AXES}}
                if w==200000:
                    prev=temporal(by)
                    for lab in ("linear","flexible"):
                        for y in YEARS:replay(prev[lab]["years"][y],oldp["common_annual"][lab]["years"][y])
                    z["common_annual"]=annual_diet(by,records,AXES[0])
                panels.append(z)
                def compact(q):
                    if q["status"]!="estimated":return q
                    return {"status":q["status"],"n_cells":q["n_cells"],"n_blocks":q["n_blocks"],"cross":{m:{c:v["cross"] for c,v in r.items()} for m,r in q["effects"].items()}}
                print(json.dumps({"diet_coupling_panel":{"width_m":w,"land":gate,"2025":{k:{j:(compact(v) if isinstance(v,dict) and "status" in v else {c:compact(q) for c,q in v.items()} if j=="leave_continent_out" else v) for j,v in r.items()} for k,r in z["2025"].items()},"common_annual":z.get("common_annual")}},allow_nan=False),flush=True)
    finally:raster.close();mem.close()
    result={"version":"railstatus-avonet-diet-plant-biomass-coupling-v1","source_replay":checks,"diet_source":diet["source"],"biomass_source_sha256":agb["source_sha256"],"panels":panels,
      "primary":"200km land80 2025, direct fruit+seed feeding minus terrestrial invertivore plant-richness slope, mean opposite-observer orientations, same cells and joint guild+full efforts, woody AGB and geographic controls, local block FE. Secondary woody AGB slope contrast, separate fruit/seed guilds, 100km, temporal retest and sensitivities.",
      "limits":["food-resource guilds from pinned AVONET derivative, not realized diets or interactions","plant richness not seed/fruit quantity or specific food-plant diversity","invertivores may depend indirectly on plants","guild-specific detectability, richness saturation, habitat and taxonomic composition can differ","same platform and geography, no causal inference","AGB not plant productivity or bird biomass","no eBird occurrences"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False))
if __name__=="__main__":main()

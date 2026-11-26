#!/usr/bin/env python3
"""Do actual railway crossings differentially connect disturbance/edge guilds?
No assertion of dispersal, identical paired design for both guilds.
"""
import argparse,json,hashlib,math
from pathlib import Path
from collections import Counter
import analyze_world_railstatus_observer_split_connectivity as s
import analyze_world_railstatus_plant_bird_biomass as b
from analyze_world_railstatus_plant_insect_bird_beta import score_fit,attach_agb
from analyze_world_railstatus_plant_bird_disjoint import combine,interval
from analyze_world_railstatus_coupling_climate import climate_extra,EXTRA
from scan_world_railstatus_disturbance_guild_pairs import key,AXES,GATE,WIDTHS,GROUPS
EFFECTS=(s.ACTIVE,s.INACTIVE)
ENV=("log_woody_AGB_center50",)+EXTRA
def join_guild(row,q,coverage,env):
    if q["min_half_coverage"]<coverage:return None
    z=dict(row);z["trait_coverage"]=q["min_half_coverage"]
    for f in ENV:
        a,c=env[row["a"]][f],env[row["b"]][f]
        z["mean_"+f]=(a+c)/2;z["diff_"+f]=abs(a-c)
    for i,g in enumerate(q["guilds"]):
        for f,v in s.symmetric_effort(g["split_counts_a"],g["split_counts_b"]).items():z[f"g{i}_"+f]=v
        for metric,name in (("jaccard","jaccard"),("simpson","simpson_similarity")):
            v=[g["split_similarity"][k][name] for k in ("01","10")]
            if any(x is None or not 0<=x<=1 for x in v):raise ValueError("Invalid guild similarity")
            z[f"g{i}_"+metric]=sum(v)/2
            z[f"g{i}_same_"+metric]=sum(g["split_similarity"][k][name] for k in ("00","11"))/2
    return z

def fields(endpoint,flex):
    f=(s.ENDPOINT if endpoint else s.PRED)+tuple(f"g{i}_"+k for i in (0,1) for k in s.EFF_PAIR if not endpoint or k.startswith("diff_"))
    f+=tuple(pre+k for k in ENV for pre in (("diff_",) if endpoint else ("mean_","diff_")))+("trait_coverage",)
    if flex:f+=tuple(k+"_sq" for k in f if "log_observations" in k or "log_observers" in k)
    return f
def estimate(rows,metric="jaccard",endpoint=False,flex=True,same=False):
    ff=fields(endpoint,flex)
    rr=[]
    for r in rows:
        q=dict(r)
        for k in ff:
            if k.endswith("_sq") and k not in q:q[k]=q[k[:-3]]**2
        rr.append(q)
    qa={"n_pairs":len(rr),"n_blocks":len({(r["continent"],r["block"]) for r in rr}),"n_active_cross":sum(r[s.ACTIVE]>0 for r in rr),"n_inactive_cross":sum(r[s.INACTIVE]>0 for r in rr)}
    fits={i:score_fit(rr,ff,f"g{i}_"+("same_" if same else "")+metric,EFFECTS,endpoint) for i in (0,1)}
    if any(v is None for v in fits.values()):return dict(qa,status="insufficient support or rank")
    for i in (0,1):
        fits[i]["active_minus_inactive"]=combine(fits[i],{s.ACTIVE:1.,s.INACTIVE:-1.})
    fits["guild1_minus_guild0"]={k:combine({i:fits[i][k] for i in (0,1)},{0:-1.,1:1.}) for k in fits[0]}
    return dict(qa,status="estimated",effects={str(i):{k:interval(v) for k,v in q.items()} for i,q in fits.items()})
def panel(rows):
    return {"primary":estimate(rows),"linear_effort":estimate(rows,flex=False),
            "endpoint_FE":estimate(rows,endpoint=True,flex=False),
            "simpson":estimate(rows,"simpson"),"same_observer_side":estimate(rows,same=True),
            "both_inactive_endpoints":estimate([r for r in rows if r["both_inactive"]]),
            "leave_continent_out":{c:estimate([r for r in rows if r["continent"]!=c]) for c in sorted({r["continent"] for r in rows})}}
def main():
    p=argparse.ArgumentParser()
    for k in ("guild","trait_cells","split","prior","edges","agb","plant","rail","status","built","tree","road","population","ne","climate","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();gs,src,prior,edge,agb,plant,rail,status,tree,road,pop=map(s.read,(a.guild,a.split,a.prior,a.edges,a.agb,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=s.verify_original_sources(plant,rail);s.verify_replay(prior,plant)
    if gs["version"]!="railstatus-disturbance-guild-pairs-v1" or not gs["exact_full_target_guild_count_and_pair_replay"] or gs["guild_gate"]!=list(GATE) or gs["trait_cells_sha256"]!=hashlib.sha256(a.trait_cells.read_bytes()).hexdigest() or any(gs[k]!=src[k] for k in ("obs_sha256","taxa_sha256","split_seed")):raise ValueError("Wrong guild source")
    if src["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not src["exact_grid_replay"] or not src["exact_plant_beta_replay"] or any(src[k]!=prior[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Wrong original split")
    if edge["version"]!="railstatus-direct-grid-boundary-v1" or not edge["source_verified"] or edge["release"]!=rail["source_Oveture_release"] or any(abs(edge["km"][c]-status["total_km"][c])>.01 for c in s.CATS):raise ValueError("Wrong rail geometry")
    if agb["version"]!="railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1" or tree["source_file_SHA256"]!=s.EXPECTED_SHA or road["source_zip_MD5"]!=s.MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong environment source")
    if s.sha_file(a.ne)!=s.ORIGINAL_NE_SHA or s.sha_file(a.climate)!=s.WORLDCLIM_SHA or pop["source_archive_SHA256"]!=s.POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen environment changed")
    built,_=s.read_source_builtup(a.built);land,_=s.source_continents(a.ne);ras,mem,_=s.worldclim_tif(a.climate)
    panels=[]
    try:
        scale,_=s.climate_scale(ras)
        for w in WIDTHS:
            grid=prior["grids"][str(w)];sp=src["pairs"][str(w)]
            keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for g in GROUPS for r in grid[g]}
            clim,_=s.source_climate_lookup(land,ras,scale,keys,w);ext,_=climate_extra(a.climate,ras,clim,w)
            tr=s.make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)]);rd=s.make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(w)]);pp=s.make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
            ag20=s.make_lookup(agb["grids"]["2020"][str(w)]);ag24=s.make_lookup(agb["grids"]["2024"][str(w)])
            el={c:{(tuple(z["a"]),tuple(z["b"])):z["crossing_source_segments"] for z in edge["edges"][str(w)][c]} for c in s.CATS}
            for g in GROUPS:
                raw={key(p):p for p in sp[g]};guild={key(p):p for p in gs["pairs"][str(w)][g]}
                if len(guild)!=len(gs["pairs"][str(w)][g]) or not set(guild)<=set(raw):raise ValueError("Guild geometry not frozen")
                rows,_,_=s.original_cell_rows([s.convert(z) for z in grid[g]],rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(w)],clim)
                rows,_=s.inject_built(rows,built[str(w)]);rows,_=s.inject_status(rows,status,w)
                for gate in (.8,.5):
                    rr,_=s.enrich_cells(rows,rd,tr,w,gate);rr,_=s.join_population(rr,pp,rd)
                    ag=attach_agb(rr,ag20,ag24,rd)
                    env={k:dict(ext[k],**{"log_woody_AGB_center50":v[2024]}) for k,v in ag.items() if k in ext}
                    base=s.build_pairs([s.adapt_pair(raw[k]) for k in guild],rr,el)
                    base=[s.split_pair(r,raw[r["a"],r["b"]]) for r in base if r["a"] in env and r["b"] in env]
                    for cov in (.5,.25):
                        z={"width":w,"land":gate,"coverage":cov,"taxon":g,"axes":{}}
                        for ax in AXES[g]:
                            rr0=[v for r in base if ax in guild[r["a"],r["b"]]["axes"] and (v:=join_guild(r,guild[r["a"],r["b"]]["axes"][ax],cov,env)) is not None]
                            z["axes"][ax]=panel(rr0)
                        panels.append(z);print(json.dumps({"disturbance_guild_panel":z},allow_nan=False),flush=True)
    finally:ras.close();mem.close()
    out={"version":"railstatus-disturbance-guild-crossing-v1","source_replay":checks,"guild_sha256":hashlib.sha256(a.guild.read_bytes()).hexdigest(),"panels":panels,
         "primary":"100km land80 each half coverage>=50%, direct (guild1-guild0) x (active-inactive) Jaccard, local block FE, all joint guild effort, regional climate and woody AGB. Plant annual vs perennial nonwoody; butterfly edge-associated vs avoidant, disturbance-associated vs avoidant secondary. 200km, coverage25, land50 and endpoint FE sensitivity.",
         "limits":["whole-landscape recorded overlap, not movement or rail flora","both endpoints contain active mapped normal railway as prior crossing design","status classified in 2026, not measured operation/maintenance in 2025","strict directional known guilds; missingness and pool-size confounding remain","no eBird"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

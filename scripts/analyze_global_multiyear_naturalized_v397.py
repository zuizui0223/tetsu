#!/usr/bin/env python3
"""v397: original global 2023-25 observations, temporal replication and expanded support.
2025 original observation-cell counts and naturalized aggregates must reproduce.
Only coarse cell aggregates and model summaries leave the runner.
"""
import argparse, json, math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from pyproj import Transformer
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    source_client,open_original_source,taxonomy_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,
    TAXA_ETAG,TAXA_BYTES,MAX_UNCERT_M,OBS_QUALITY)
from analyze_global_rail_documented_naturalized_v379 import source_lists,vascular_names,cell_region,grid_counts,grid_key
from analyze_original_world_railkm_plant_richness_Jaccard_v371 import (
    verify_original_sources,source_climate_lookup,original_cell_rows,block_weights)
from analyze_global_original_railkm_biodiversity_GHSL_control_v376c import read_source_builtup,inject_built,ORIGINAL_NE_SHA
from generate_global_equal_area_climate_land_grid_v351 import WORLDCLIM_SHA,source_continents,worldclim_tif,climate_scale,sha_file
from analyze_global_naturalized_count_decomposition_v389 import source_rows,OUTCOMES,fit,FIELDS,BOOT
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup
from analyze_global_naturalized_population_v392 import add_population
from analyze_global_naturalized_local_swap_v393 import supported_local
from analyze_global_naturalized_joint_population_local_v394 import FULL,POP_SHA,groups,mat,supported
from analyze_global_rail_presence_intensity_v395 import add_exposure,mixed,BASE,coefficients,model
from analyze_global_rail_presence_matched_v396 import standardized,optimal_pairs,compare,pair_weights,MIN_PAIRS,MIN_BLOCKS
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
YEARS=("2023","2024","2025")
W=100000
PADS=(.3,0.)
LAND=(.5,.8)
SEED=39720261125
EXPECTED_MATCH={(.3,.5):(15,2),(.3,.8):(14,2),(0.,.5):(33,11),(0.,.8):(31,11)}
EXPECTED_2025_PAIR_FRACTION={(0.,.5):1.2332559936662024,(0.,.8):1.4503846917211587}
def entry():return [0,set(),set()]
def scan(client,canonical,plant):
    data={y:defaultdict(entry) for y in YEARS}
    qa=Counter();annual=Counter()
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rows,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    try:
        required={"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}
        if not required<=set(rows.fieldnames or []):raise ValueError("Original observation schema changed")
        for r in rows:
            qa["all_rows"]+=1
            if qa["all_rows"]%10000000==0:
                print(json.dumps({"all_rows":qa["all_rows"],"eligible_by_year":annual}),flush=True)
            y=str(r["observed_on"])[:4]
            if y not in YEARS or not str(r["observed_on"]).startswith(y+"-"):continue
            if str(r["quality_grade"]).strip().lower().replace(" ","_").replace("-","_")!=OBS_QUALITY:continue
            try:
                tid=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"])
                lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if tid is None or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k is None:continue
            annual[y]+=1
            a=data[y][k];a[0]+=1;a[1].add(tid)
            if r["observer_id"]:a[2].add(int(r["observer_id"]))
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or raw.hash.hexdigest()!=plant["source_summary"]["source_iNat_original_2026_10_27_observations_gzip_SHA256_full"] or qa["all_rows"]!=288721938 or annual["2025"]!=7827936:
        raise ValueError("Full original source replay mismatch")
    old={grid_key(r):(r["source_year_2025_photo_linked_plants_observation_count"],r["source_2025_plant_species_richness"],r["source_2025_distinct_observers_with_photo_linked_plant_records"]) for r in grid_counts(plant,W)}
    new={k:(a[0],len(a[1]),len(a[2])) for k,a in data["2025"].items()}
    if old!=new:raise ValueError("Exact v369 2025 all-world coarse counts/observers/richness not reproduced")
    return data,{"all_rows":qa["all_rows"],"annual_eligible_original_plant_records":dict(annual),"original_2025_all_cell_counts_exact":True,
                 "source_observation_gzip_sha256":raw.hash.hexdigest()}
def pooled(data):
    out=defaultdict(entry)
    for y in YEARS:
        for k,a in data[y].items():
            b=out[k];b[0]+=a[0];b[1].update(a[1]);b[2].update(a[2])
    return out
def region_cells(geo,g):
    inv=Transformer.from_crs(6933,4326,always_xy=True)
    return {pad:{k:r for k in sorted(geo) if (r:=cell_region(k,W,g,inv,pad)) is not None} for pad in PADS}
def actual_annual_grid(data):
    # Adapter to the frozen v371 field names. Values are genuinely from
    # the selected observation period; no 2025 observations copied into other years.
    return [{"source_equal_area_grid_ix":k[0],"source_equal_area_grid_iy":k[1],
             "source_year_2025_photo_linked_plants_observation_count":a[0],
             "source_2025_distinct_observers_with_photo_linked_plant_records":len(a[2]),
             "source_2025_plant_species_richness":len(a[1])} for k,a in sorted(data.items())]
def fraction_rows(rows,data,vasc,region,listed):
    o=[]
    for r in rows:
        k=r["key"];rid=region.get(k)
        if rid is None or not listed[rid]:continue
        ids=data[k][1];n=sum(i in vasc for i in ids)
        if n<20:continue
        d=sum(i in vasc and vasc[i] in listed[rid] for i in ids)
        f=d/n;x=math.log1p(d);u=math.log1p(n-d)
        o.append(dict(r,region=rid,vascular_richness=n,documented_richness=d,known_fraction=f,
            log_documented_regional_species_pool=math.log1p(len(listed[rid])),
            log1p_documented=x,log1p_unlisted=u,log1p_count_contrast=x-u,known_fraction_pp=100*f))
    return o
def exact_2025(rows,reference,pad):
    a={r["key"]:r for r in rows};b={r["key"]:r for r in source_rows(reference)}
    if set(a)!=set(b):raise ValueError("v383 2025 geography or vascular gate changed")
    for k,r in a.items():
        for f in ("region","continent","block","observations","observers","species","vascular_richness","documented_richness"):
            if r[f]!=b[k][f]:raise ValueError("Exact v383 2025 cell field mismatch: "+f)
        for f in ("rail_km","known_fraction","log_documented_regional_species_pool"):
            if abs(r[f]-b[k][f])>1e-10:raise ValueError("Exact v383 numeric field mismatch: "+f)
def support(rows,period,mode):
    # Pooled support is a new cumulative observation endpoint, not
    # three independent samples; stricter pooled effort is separate.
    if mode=="base":return rows
    if period=="2023-2025":
        n,b=(150,15) if mode=="annualized_support" else (600,60)
    else:n,b=(200,20)
    return [r for r in rows if r["observations"]>=n and r["observers"]>=b]
def compute(rows,period,pad,land,mode):
    s=eligible_region_rows(support(rows,period,mode));local=[add_exposure(r) for r in supported_local(s)]
    mix=mixed(local);seed=SEED+int(pad*100)+int(land*10)
    m=compare(mix,seed)
    fields=BASE+("rail_present",)
    binary=model(mix,fields,("rail_present",),True,seed)
    return {"period":period,"mode":mode,"n_regional_supported":len(s),"n_local":len(local),
        "regional_count_model":fit(s,FULL,True),"binary_presence_local":binary,"matched_presence":m}
def fixed_pairs(rows25,high_effort):
    s=eligible_region_rows([r for r in rows25 if not high_effort or (r["observations"]>=200 and r["observers"]>=20)])
    m=mixed([add_exposure(r) for r in supported_local(s)])
    if not m:return []
    pairs=optimal_pairs(m,standardized(m))
    return [(m[i]["key"],m[j]["key"]) for i,j in pairs]
def fixed_temporal(pairs,annual,high_effort):
    maps={y:{r["key"]:r for r in annual[y]} for y in YEARS}
    def eligible(y,k):
        r=maps[y].get(k)
        return r is not None and (not high_effort or (r["observations"]>=200 and r["observers"]>=20))
    # Use the same 2025-selected pairs in both compared years. Do not
    # rematch with earlier outcomes or loosen support thresholds.
    o=[]
    for past in ("2024","2023"):
        keys=[(a,b) for a,b in pairs if all(eligible(y,k) for y in (past,"2025") for k in (a,b))]
        blocks=[(maps["2025"][a]["continent"],maps["2025"][a]["block"]) for a,b in keys]
        rec={"past_year":past,"original_2025_pairs":len(pairs),"both_years_supported_pairs":len(keys),
             "n_blocks":len(set(blocks)),"pairs_by_continent":dict(Counter(x[0] for x in blocks))}
        if len(keys)<MIN_PAIRS or len(set(blocks))<MIN_BLOCKS:
            rec["status"]="insufficient same-pair temporal support";o.append(rec);continue
        yy={}
        for y in (past,"2025"):
            yy[y]=np.asarray([[maps[y][b][t]-maps[y][a][t] for t in OUTCOMES] for a,b in keys],float)
        w=pair_weights(blocks)
        rec["same_pairs_means"]={y:dict(zip(OUTCOMES,map(float,np.average(a,axis=0,weights=w)))) for y,a in yy.items()}
        g=defaultdict(list);bc=defaultdict(list)
        for i,k in enumerate(blocks):g[k].append(i)
        for k in sorted(g):bc[k[0]].append(k)
        rng=np.random.default_rng(SEED+int(past));draws={past:[],"2025":[],"2025_minus_past":[]}
        for _ in range(BOOT):
            bw=np.zeros(len(keys))
            for bs in bc.values():
                mul=rng.multinomial(len(bs),np.full(len(bs),1/len(bs)))
                for k,n in zip(bs,mul):bw[g[k]]=n/len(g[k])
            a=np.average(yy[past],axis=0,weights=bw);b=np.average(yy["2025"],axis=0,weights=bw)
            draws[past].append(a);draws["2025"].append(b);draws["2025_minus_past"].append(b-a)
        rec["block_percentile_95_descriptive"]={y:{t:{"low":float(np.quantile(np.asarray(a)[:,j],.025)),"high":float(np.quantile(np.asarray(a)[:,j],.975))} for j,t in enumerate(OUTCOMES)} for y,a in draws.items()}
        rec["status"]="same geography different observation years, 2025 exposure matching"
        o.append(rec)
    return o
def main():
    p=argparse.ArgumentParser()
    for k in ("plants","rail","built","ne","climate","glonaf","road","tree","population","v383","out","aggregates"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text())
    plant,rail,road,tree,pop,old=map(read,(a.plants,a.rail,a.road,a.tree,a.population,a.v383))
    verify_original_sources(plant,rail)
    if old["version"]!="v383-exact-v379d-coarse-cell-aggregate" or old["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:raise ValueError("Wrong original v383")
    if road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Wrong roads")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:raise ValueError("Wrong tree")
    if pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Wrong pop")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA:raise ValueError("Wrong climate/map")
    g,listed,lqa=source_lists(a.glonaf)
    client=source_client();canonical,_,_,_,sha=taxonomy_source(client);vasc,sha2=vascular_names(client,canonical)
    if sha!=sha2 or sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Original taxonomy hash mismatch")
    data,qa=scan(client,canonical,plant);data["2023-2025"]=pooled(data)
    land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster)
        geo,_=source_climate_lookup(land,raster,scale,set(data["2023-2025"]),W)
        reg=region_cells(geo,g)
    finally:raster.close();mem.close()
    built,_=read_source_builtup(a.built)
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    references={p["region_guard_degrees"]:p for p in old["panels"] if p["width_m"]==W}
    per={};aggregate=[];design=[]
    for period,values in data.items():
        rows,qa0,_=original_cell_rows(actual_annual_grid(values),rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],geo)
        b,qa1=inject_built(rows,built[str(W)])
        if period=="2025" and (len(rows)!=3455 or len(b)!=3455):raise ValueError("Original 2025 full model support changed")
        for pad in PADS:
            r=fraction_rows(b,values,vasc,reg[pad],listed)
            if period=="2025":exact_2025(r,references[pad],pad)
            aggregate.append({"period":period,"region_guard_degrees":pad,"rows":r})
            design.append({"period":period,"pad":pad,"original_eligible_global_cells":len(rows),"with_builtup":len(b),
                           "assigned_eligible_vascular_cells":len(r),"cell_qa":qa0,"builtup_qa":qa1})
            for land in LAND:
                m,e=enrich_cells(r,rd,tr,W,land);z,f=add_population(m,pp,rd)
                per[period,pad,land]=z
    # Publish only coarse-cell species/effort counts and pre-existing covariates.
    a.aggregates.parent.mkdir(parents=True,exist_ok=True)
    a.aggregates.write_text(json.dumps({"version":"v397-original-2023-2025-coarse-aggregates","source_replay":qa,
        "panels":aggregate},allow_nan=False))
    out=[]
    for pad in PADS:
        for land in LAND:
            annual={y:per[y,pad,land] for y in YEARS}
            fixed=[]
            for h in (False,True):
                pairs=fixed_pairs(annual["2025"],h)
                if len(pairs)!=EXPECTED_MATCH[pad,land][int(h)]:raise ValueError("Original v396 matching support changed")
                fixed.append({"high_effort":h,"results":fixed_temporal(pairs,annual,h)})
            panel={"region_guard_degrees":pad,"min_land_fraction":land,"fixed_2025_pair_temporal_retest":fixed,"models":[]}
            for period in (*YEARS,"2023-2025"):
                modes=("base","high_effort") if period!="2023-2025" else ("base","annualized_support","high_effort")
                for mode in modes:
                    m=compute(per[period,pad,land],period,pad,land,mode)
                    if period=="2025" and mode=="base" and (pad,land) in EXPECTED_2025_PAIR_FRACTION:
                        v=m["matched_presence"]["matched_blockweighted_difference"]["known_fraction_pp"]
                        if abs(v-EXPECTED_2025_PAIR_FRACTION[pad,land])>1e-9:raise ValueError("Exact v396 2025 paired effect mismatch")
                    panel["models"].append(m)
            out.append(panel)
            # Print a concise source-backed result, retain full estimates in artifact.
            summary=[{"period":m["period"],"mode":m["mode"],"n_local":m["n_local"],
                "match_n":m["matched_presence"].get("n_pairs",0),"match_blocks":m["matched_presence"].get("n_blocks",0),
                "match_status":m["matched_presence"]["status"],
                "match_difference":m["matched_presence"].get("matched_blockweighted_difference"),
                "match_interval":m["matched_presence"].get("block_percentile_95_descriptive"),
                "binary_coef":m["binary_presence_local"].get("coefficients"),
                "binary_interval":m["binary_presence_local"].get("block_percentile_95_descriptive")} for m in panel["models"]]
            print(json.dumps({"v397_panel":{"pad":pad,"land":land,"models":summary,"fixed":fixed}},allow_nan=False),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({"version":"v397-original-three-year-temporal-and-support-expansion","source_replay":qa,
        "source_GloNAF":lqa,"design":design,"panels":out,
        "limitations":["historical observation years, frozen 2026 rail map; not railway opening history",
          "same geography and possible repeated observers/species, not independent ecological systems",
          "historical GloNAF checklist applied to all years, not year-specific naturalization state",
          "pooled species/observer unions are cumulative richness, not independent annual replicates",
          "no raw individual GPS or species/observer IDs persisted outside workflow",
          "mapped zero not confirmed absence, nonlisted not native, no causal or temporal-invasion claims"]},allow_nan=False))
if __name__=="__main__":main()

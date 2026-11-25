#!/usr/bin/env python3
"""Complementary support-form functional diversity test on the same verified
growthform coarse artifact. No new observation scan or source trait imputation.
"""
import argparse,json,math
from pathlib import Path
from analyze_world_railstatus_growthform import (
    FULL,POP,ACTIVE,INACTIVE,POP_SHA,join_population,inject_status,local_support,
    CATS,LOG,read,make_lookup,enrich_cells,read_source_builtup,ORIGINAL_NE_SHA,
    WORLDCLIM_SHA,EXPECTED_SHA,MD5,sha_file,source_continents,worldclim_tif,
    climate_scale,source_climate_lookup,original_cell_rows,verify_original_sources,
    inject_built,GRID,fit,flexible,EXTRA,SHA,W)
from analyze_world_railstatus_growthform_missing_bounds import evaluate as joint_bounds
CLASSES=("F","C","E","P","M","A","H")
N_DRAW=20
TARGETS=("support_diversity","free_vs_nonfree","nonfree_within")
SECONDARY=("rarefied_support_richness","free_standing_fraction")
FIELDS=FULL+(POP,"log_species")

def expected_present(n,k,m=N_DRAW):
    if not 0<=k<=n or not 0<=m<=n:raise ValueError("Invalid rarefaction")
    if k==0 or m==0:return 0.
    if n-k<m:return 1.
    # Exact hypergeometric absence probability, without replacement.
    return -math.expm1(math.fsum(math.log1p(-k/(n-i)) for i in range(m)))

def diversity(counts):
    if not set(counts)<=set(CLASSES):raise ValueError("Unknown original support class")
    if any(type(v)!=int or v<0 for v in counts.values()):raise ValueError("Invalid class counts")
    n=sum(counts.values())
    if n<N_DRAW:raise ValueError("Insufficient known species")
    f=counts.get("F",0);non=n-f
    denom=n*(n-1)
    # Probability of different growth support forms for two distinct known
    # species. Decompose exactly into free vs nonfree and within nonfree.
    between=2*f*non/denom
    within=(non*non-sum(counts.get(k,0)**2 for k in CLASSES if k!="F"))/denom
    gini=(n*n-sum(v*v for v in counts.values()))/denom
    if abs(gini-between-within)>1e-12:raise ValueError("Diversity decomposition failed")
    return dict(support_diversity=gini,free_vs_nonfree=between,nonfree_within=within,
                rarefied_support_richness=math.fsum(expected_present(n,counts.get(k,0)) for k in CLASSES),
                free_standing_fraction=f/n)

def attach(rows,traits,coverage):
    out=[]
    for r in rows:
        t=traits.get(r["key"])
        if t is None:raise ValueError("Missing original trait count cell")
        if any(t[k]!=r[k] for k in ("species","observations","observers")):raise ValueError("Original count mismatch")
        n=t["known_support"]
        if sum(t["support_counts"].values())!=n or not 0<=n<=r["species"]:raise ValueError("Support counts do not replay")
        if n<N_DRAW or n/r["species"]<coverage:continue
        out.append(dict(r,coverage=n/r["species"],known_support=n,log_known_support=math.log1p(n),
                        **diversity(t["support_counts"])))
    return out

def evaluate(rows):
    r=local_support(rows)
    return {"n_before_local":len(rows),
            "primary":fit(r,FIELDS,TARGETS),
            "secondary":fit(r,FIELDS,SECONDARY),
            "without_total_richness":fit(r,FULL+(POP,),TARGETS),
            "coverage_known_richness_flexible_effort":fit([flexible(x) for x in r],FIELDS+EXTRA+("coverage","log_known_support"),TARGETS),
            "high_effort":fit(local_support([x for x in r if x["observations"]>=200 and x["observers"]>=20]),FIELDS,TARGETS),
            "both_status_positive":fit(local_support([x for x in r if x["active_normal_km"]>0 and x["inactive_normal_km"]>0]),FIELDS,TARGETS),
            "leave_continent_out":{c:fit(local_support([x for x in r if x["continent"]!=c]),FIELDS,TARGETS)
                                    for c in sorted({x["continent"] for x in r})}}

def main():
    p=argparse.ArgumentParser()
    for n in ("traits","plant","rail","status","built","tree","road","population","ne","climate","out"):
        p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();src,plant,rail,status,tree,road,pop=map(read,(a.traits,a.plant,a.rail,a.status,a.tree,a.road,a.population))
    checks=verify_original_sources(plant,rail)
    if src["version"]!="railstatus-growthform-2025-coarse-v1" or src["source"]["commit"]!=SHA or src["source_replay"]["rows"]!=288721938 or src["source_replay"]["plants"]!=7827936:
        raise ValueError("Wrong original growthform source")
    if src["source"]["taxonomy_SHA256"]!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Taxonomy SHA mismatch")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or road["source_zip_MD5"]!=MD5 or road["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":raise ValueError("Frozen road/tree mismatch")
    if sha_file(a.ne)!=ORIGINAL_NE_SHA or sha_file(a.climate)!=WORLDCLIM_SHA or pop["source_archive_SHA256"]!=POP_SHA or pop["relative_population_sum_difference"]>.01:raise ValueError("Frozen climate/map/pop mismatch")
    if not status["version"].startswith("railstatus-v1"):raise ValueError("Wrong railway status")
    traits={tuple(r["key"]):r for r in src["cells"]}
    if len(traits)!=src["n_cells"]:raise ValueError("Duplicate original trait cells")
    original=plant["source_grids"][str(W)][GRID]
    keys={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]) for r in original}
    built,_=read_source_builtup(a.built);land,_=source_continents(a.ne);raster,mem,_=worldclim_tif(a.climate)
    try:
        scale,_=climate_scale(raster);clim,_=source_climate_lookup(land,raster,scale,keys,W)
    finally:raster.close();mem.close()
    rows,_,_=original_cell_rows(original,rail["source_global_gridded_rail_exposure_km_by_cell_width"][str(W)],clim)
    rows,_=inject_built(rows,built[str(W)]);rows,worst=inject_status(rows,status,W)
    tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(W)])
    rd=make_lookup(road["source_global_GRIP4_road_density_by_grid"][str(W)])
    pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(W)])
    panels=[]
    for land in (.5,.8):
        m,_=enrich_cells(rows,rd,tr,W,land);m,_=join_population(m,pp,rd)
        for coverage in (.5,.25):
            r=attach(m,traits,coverage)
            # Same cells with identifiable woodiness as a secondary matched
            # support sensitivity, not a new primary result.
            joint=[x for x in r if traits[x["key"]]["known_woodiness"]>=N_DRAW and
                   traits[x["key"]]["known_woodiness"]/x["species"]>=coverage]
            panel={"land":land,"coverage":coverage,"base_n":len(m),
                   "support_diversity":evaluate(r),
                   "common_woodiness_support":evaluate(joint),
                   "joint_richness_and_missingness_bounds":joint_bounds(joint,traits)}
            panels.append(panel);print(json.dumps({"growthform_diversity_panel":panel},allow_nan=False),flush=True)
    out={"version":"railstatus-growthform-support-diversity-v1","source_replay":checks,
         "growthform_source":src["source"],"original_coverage":src["coverage"],
         "rail_km_replay_error":worst,"panels":panels,
         "primary":"100km, land>=50%, >=20 source-known support forms and >=50% coverage; active-minus-inactive slope on within-cell Gini-Simpson support-form diversity; negative predicted by selective filtering",
         "secondary":"exact decomposition into free/nonfree and within nonfree contributions; expected category richness at 20 known species and free-standing fraction; land/coverage/effort/joint-woodiness sensitivity",
         "limits":["categorical support forms only, not multidimensional functional diversity or among-cell convergence",
                   "source-known taxa, no imputation; not a within-family support-form null",
                   "whole landscape records, not actual railway microhabitat",
                   "species equally weighted, not abundance",
                   "posthoc complementary endpoint; approximate cluster intervals and correlated sensitivities"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()

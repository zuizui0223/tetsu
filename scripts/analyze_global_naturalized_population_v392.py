#!/usr/bin/env python3
"""v392: population-adjusted documented species count contrast using exact v389 aggregates."""
import argparse,json,math
from pathlib import Path
from collections import Counter
from analyze_global_naturalized_count_decomposition_v389 import source_rows,fit,FIELDS,EXPECTED
from analyze_global_naturalized_quality_effort_v390 import PRIOR
from analyze_global_rail_naturalized_road_tree_v383 import eligible_region_rows,regression,ADJ
from analyze_global_rail_road_tree_richness_v382 import enrich_cells,make_lookup,LAND_THRESHOLDS
from aggregate_global_tree_consensus_EarthEnv_v380a import EXPECTED_SHA
from aggregate_original_GRIP4_road_density_v381a import MD5
def add_population(rows,pop,roads):
    matched=[];qa=Counter()
    for r in rows:
        p=pop.get(r["key"]);rd=roads[r["key"]]
        land=float(rd["source_GRIP4_total_land_area_m2"])
        if p is None:qa["no_population_source"]+=1;continue
        n=float(p["source_GHS_POP2020_people"]);a=float(p["source_GHS_POP2020_approx_valid_pixel_area_m2"])
        if not math.isfinite(n) or n<0 or land<=0 or not math.isfinite(a) or a<0:
            raise ValueError("Invalid original population support")
        # This is an aggregate area-adequacy check, NOT exact land-mask overlap.
        if a<.90*land:qa["population_source_area_less_than_90pct_GRIP_land_area"]+=1;continue
        matched.append(dict(r,log_population_density=math.log1p(n/(land/1e6))))
    return matched,dict(qa)
def run(cells,roads,tree,pop):
    if cells["version"]!="v383-exact-v379d-coarse-cell-aggregate" or cells["source_replay"]!={"all_rows":288721938,"all_original_eligible_plant_observations":7827936,"all_vascular_observations":7726031}:
        raise ValueError("Not original v383 aggregate")
    if roads["source_zip_MD5"]!=MD5 or roads["source_zip_SHA256"]!="782fe4b609befbb294528d2b7abce212df4d8f883c07775134fe9a6372dae7d9":
        raise ValueError("Not original GRIP4")
    if tree["source_file_SHA256"]!=EXPECTED_SHA or tree["source_pixel_QA"]["joint_valid_source_pixels_in_grid"]!=756864000:
        raise ValueError("Not original EarthEnv")
    if pop["version"]!="v391-original-JRC-GHS-POP-2020-1km-global-counts" or pop["relative_population_sum_difference"]>.01 or not 7e9<pop["source_full_original_population_sum"]<8.5e9:
        raise ValueError("Not valid original GHS-POP")
    panels=[]
    for p in cells["panels"]:
        w,pad=p["width_m"],p["region_guard_degrees"];rr=source_rows(p)
        rd=make_lookup(roads["source_global_GRIP4_road_density_by_grid"][str(w)])
        tr=make_lookup(tree["source_global_100km_200km_tree_consensus"][str(w)])
        pp=make_lookup(pop["source_global_GHS_POP2020_by_grid"][str(w)])
        for g in LAND_THRESHOLDS:
            rows,_=enrich_cells(rr,rd,tr,w,g);old=eligible_region_rows(rows)
            n,m,prior=EXPECTED[w,pad,g]
            d=fit(old,FIELDS,False)
            if len(rows)!=n or len(old)!=m or abs(d["coefficients"]["log_rail_km"]["log1p_documented"]-PRIOR[w,pad,g][0])>1e-9:
                raise ValueError("v389 original count coefficient not reproduced")
            joined,qa=add_population(rows,pp,rd)
            out={"width_m":w,"region_guard_degrees":pad,"min_land_fraction":g,
                 "before_population_support":len(rows),"after_population_support":len(joined),"population_support_excluded":qa,"fits":[]}
            for mode in ("all","high_effort"):
                r=joined if mode=="all" else [x for x in joined if x["observations"]>=200 and x["observers"]>=20]
                r=eligible_region_rows(r)
                out["fits"].append({"support":mode,"same_support_before_population":fit(r,FIELDS,True),
                    "same_support_after_population":fit(r,FIELDS+("log_population_density",),True),
                    "original_fraction_conditional_total_after_population":regression(
                        r,ADJ+("log_population_density",),("log_rail_km","log_road_density","log_population_density"),True)})
            panels.append(out);print(json.dumps({"v392_panel":out},allow_nan=False),flush=True)
    return {"version":"v392-v389-original-population-adjusted-documented-counts",
            "source_population_SHA256":pop["source_archive_SHA256"],"panels":panels,
            "limits":["population 2020 modeled residential counts, not 2025 travel or rail traffic",
            "aggregate coverage check not exact land-mask overlap, historical road/tree",
            "observed documented naturalized richness, unlisted is not native",
            "descriptive cross-sectional same-data sensitivity, no causal dispersal or temporal change"]}
def main():
    p=argparse.ArgumentParser()
    for k in ("cells","road","tree","population","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    o=run(read(a.cells),read(a.road),read(a.tree),read(a.population));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
if __name__=="__main__":main()

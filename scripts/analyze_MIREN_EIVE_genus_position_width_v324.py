#!/usr/bin/env python3
"""v324: MIREN European mountain ROAD thermal-position/breadth genus partition.

Actual MIREN Plot1 roadside vs Plot3 52-102m interior, both 100m2.
This is an external road ecology comparator; NOT a railway field sample.
"""
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path

from analyze_original_MIREN_EIVE_road_vs_interior_v323 import (
    read_miren, select_first_complete, analyze_panel)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Slovenia_thermal_generalist_section_v319 import section_test, section_drop

SYSTEMS=3
MIN_SPECIES=5
MIN_COVERAGE=0.5
MIN_ROAD_SITES=4
MIN_REGION_ROADS=2

def mean(vals):
    return sum(vals)/len(vals) if vals else None

def eligible_roads(pairs,eive):
    candidates=defaultdict(list)
    for p in pairs:
        bs={s for s in p["roadside"] if s in eive and eive[s]["systems"]>=SYSTEMS}
        bt={s for s in p["inland"] if s in eive and eive[s]["systems"]>=SYSTEMS}
        if (len(bs)<MIN_SPECIES or len(bt)<MIN_SPECIES
            or len(bs)/len(p["roadside"])<MIN_COVERAGE
            or len(bt)/len(p["inland"])<MIN_COVERAGE):
            continue
        candidates[(p["region"],p["road"])].append(
            {"site":p["site"],"year":p["year"],"roadside":bs,"inland":bt})
    enough_roads={k:v for k,v in candidates.items() if len(v)>=MIN_ROAD_SITES}
    counts=defaultdict(int)
    for reg,road in enough_roads:
        counts[reg]+=1
    selected={k:v for k,v in enough_roads.items() if counts[k[0]]>=MIN_REGION_ROADS}
    return selected

def partition(selected,eive,metric):
    regions=sorted({r for r,_ in selected})
    region_species=defaultdict(set)
    for (reg,road),sites in selected.items():
        for p in sites:
            region_species[reg].update(p["roadside"]|p["inland"])
    region_genus=defaultdict(dict)
    region_groups=defaultdict(lambda:defaultdict(list))
    for reg,taxa in region_species.items():
        for s in taxa:
            region_groups[reg][s.split()[0]].append(s)
        for gen,names in region_groups[reg].items():
            region_genus[reg][gen]=mean([eive[s][metric] for s in names])

    allroad=defaultdict(list)
    siteinfo={}
    for (reg,road),pairs in selected.items():
        for p in pairs:
            a=p["roadside"];b=p["inland"]
            full=mean([eive[s][metric] for s in a])-mean([eive[s][metric] for s in b])
            between=(mean([region_genus[reg][s.split()[0]] for s in a])-
                     mean([region_genus[reg][s.split()[0]] for s in b]))
            within=full-between
            if abs(full-between-within)>1e-12:
                raise AssertionError("MIREN genus CWM additivity violated")
            allroad[(reg,road)].append((full,between,within))
            siteinfo[(reg,road,p["site"])]=(a,b)
    road_effects={k:[mean([r[i] for r in v]) for i in range(3)]
                  for k,v in allroad.items()}
    regional={i:defaultdict(list) for i in range(3)}
    for (reg,road),deltas in road_effects.items():
        for i,x in enumerate(deltas):
            regional[i][reg].append(x)
    effects=[mean([mean(v) for v in regional[i].values()]) for i in range(3)]
    if not math.isclose(effects[0],sum(effects[1:]),rel_tol=0,abs_tol=1e-12):
        raise AssertionError("Road/region equal-weight genus decomposition violated")

    # Independently computed source main analysis must match directly.
    allpairs=[]
    for (reg,road),rows in selected.items():
        for x in rows:
            allpairs.append({"region":reg,"road":road,"site":x["site"],
                "year":x["year"],"roadside":x["roadside"],"inland":x["inland"]})
    reference=analyze_panel(allpairs,eive,metric,SYSTEMS,MIN_COVERAGE,MIN_SPECIES)
    if reference["supported_region_codes"]!=len(regions):
        raise ValueError("v323 same-source region support mismatch")
    if not math.isclose(reference["region_equal_mean_roadside_minus_inland"],
                        effects[0],rel_tol=0,abs_tol=1e-10):
        raise ValueError("v323 source reference effect mismatch")

    # Species-level coefficients preserve region -> road -> site hierarchy.
    coef=defaultdict(lambda:defaultdict(float))
    for (reg,road),sites in selected.items():
        nroads=len({rd for rr,rd in selected if rr==reg})
        weight=1/(len(regions)*nroads*len(sites))
        for s in sites:
            for species in s["roadside"]:
                coef[reg][species]+=weight/len(s["roadside"])
            for species in s["inland"]:
                coef[reg][species]-=weight/len(s["inland"])
    actual=0
    for reg,taxon_coefs in coef.items():
        for species,c in taxon_coefs.items():
            actual+=c*(eive[species][metric]-region_genus[reg][species.split()[0]])
    if not math.isclose(actual,effects[2],rel_tol=0,abs_tol=1e-10):
        raise AssertionError("Source-road hierarchical permutation weights mismatch")

    return {
        "metric":metric,"regions":len(regions),
        "road_ids":len(selected),"source_first_visit_sites":len(siteinfo),
        "region_source_species_pools":{r:len(v) for r,v in sorted(region_species.items())},
        "region_genus_counts":{r:len(v) for r,v in sorted(region_groups.items())},
        "region_multi_species_genera":{
            r:sum(len(x)>=2 for x in gr.values()) for r,gr in sorted(region_groups.items())},
        "total_region_equal_roadside_minus_inland":effects[0],
        "between_genus_contribution":effects[1],
        "within_genus_species_residual_contribution":effects[2],
        "regional_contributions":{
            reg:{"total":mean(regional[0][reg]),
                 "between":mean(regional[1][reg]),
                 "within":mean(regional[2][reg])} for reg in regions},
        "component_region_direction":{
           name:{"sign_flip":section_test(regional[i]),
                 "leave_one_region_out":section_drop(regional[i])}
           for i,name in enumerate(("total","between_genus","within_genus"))},
        "_permutation":(coef,region_genus,region_groups,eive,metric)}

def genus_shuffle(result,draws=999,seed=324):
    coef,gmeans,groups,eive,metric=result["_permutation"]
    rng=random.Random(seed)
    stats=[]
    for _ in range(draws):
        val=0
        for reg,grp in sorted(groups.items()):
            for g,taxa in sorted(grp.items()):
                if len(taxa)<=1:continue
                taxa=sorted(taxa)
                items=[eive[t][metric] for t in taxa]
                rng.shuffle(items)
                for taxon,trait in zip(taxa,items):
                    val+=coef[reg].get(taxon,0)*(trait-gmeans[reg][g])
        stats.append(val)
    stats.sort()
    actual=result["within_genus_species_residual_contribution"]
    k=sum(abs(x)>=abs(actual)-1e-12 for x in stats)
    return {"null":"shuffle EIVE temperature values among species WITHIN GENUS and WITHIN each European region's observed species pool",
           "draws":draws,"seed":seed,"actual":actual,"null_mean":mean(stats),
           "null_95pct":[stats[int(.025*draws)],stats[int(.975*draws)]],
           "two_sided_fraction_plus_one":(k+1)/(draws+1),
           "caution":"Conditional trait-label permutation, not independent plots or a randomized roadside treatment"}

def analyze(miren_path,eive_path):
    source,sourceqa=read_miren(miren_path)
    eive,eiveqa=eive_table(eive_path)
    first=select_first_complete(source,2010)
    sel=eligible_roads(first,eive)
    if not sel:raise ValueError("No eligible European original paired mountain roads")
    panels=[]
    for metric in ("position","width"):
        p=partition(sel,eive,metric)
        p["within_genus_trait_permutation"]=genus_shuffle(p)
        del p["_permutation"]
        panels.append(p)
    return {"version":"v324","context":"FOUR EUROPEAN MIREN ROADS REGIONS ONLY, NOT RAILWAY",
      "MIREN_Zenodo_DOI":"10.5281/zenodo.15749536",
      "EIVE_Zenodo_DOI":"10.5281/zenodo.7534792",
      "source_miren_qa":sourceqa,"source_eive_qa":eiveqa,
      "panels":panels,
      "limitations":[
        "Four MIREN region codes from one survey program are not four independent publications",
        "Roadside vs inland 100m2 surveys do not correspond to Slovenia railway shoulder vs toe zones",
        "Genus source pools are constructed within each region from observed MIREN plot taxa, not independent regional vegetation surveys",
        "Genus-centred EIVE components are descriptive; within can change with genus relative representation",
        "Year>=2010 first paired site used; original protocol versions can still differ across regions",
        "No actual ballast temperature, train airflow, seed-carriage mechanism or railway vs nonrail climate offset estimated"]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--miren-csv",required=True)
    ap.add_argument("--eive-xlsx",required=True)
    ap.add_argument("--output",required=True)
    x=ap.parse_args()
    result=analyze(x.miren_csv,x.eive_xlsx)
    dst=Path(x.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"version":"v324","panels":[{
       k:v.get(k) for k in ("metric","regions","road_ids","source_first_visit_sites",
          "total_region_equal_roadside_minus_inland","between_genus_contribution",
          "within_genus_species_residual_contribution",
          "within_genus_trait_permutation","regional_contributions")}
          for v in result["panels"]]},ensure_ascii=False),flush=True)

if __name__=="__main__": main()

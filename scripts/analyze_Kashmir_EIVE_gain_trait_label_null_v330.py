#!/usr/bin/env python3
"""v330 source-correct null for Kashmir recorded species thermal turnover.

Recorded taxon membership & 2014/2017 site-zone incidence are FIXED.
Shuffle independently published species EIVE values across source-matched
taxa *globally* or within original genus, using ONE consistent trait mapping
per permutation for all sites. This is a conditional source taxon-label null,
NOT field treatment randomization, causation or a confirmatory p-value.
"""
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import (
    species_names_from_publisher,cohort)
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_EIVE_loss_gain_identity_v329 import components

YEAR=(2014,2017)
REF="Interior Plot 2"
METRIC="position"
SYSTEMS=3
DRAWS=1999
SEED=330

def validated_source_taxa(plots,names,eive):
    source=set()
    for v in plots.values():
        if v is None:continue
        source.update(names[i] for i in v if names[i] in eive and
                      eive[names[i]]["systems"]>=SYSTEMS)
    return sorted(source)

def site_cohorts(plots,names,eive,reference,min_n):
    out=[]
    for site in cohort(plots):
        states={}
        ok=True
        for zone in ("Track",reference):
            years={}
            for year in YEAR:
                indices=plots.get((site,year,zone))
                if indices is None:
                    ok=False;break
                years[year]={names[i] for i in indices
                        if names[i] in eive and eive[names[i]]["systems"]>=SYSTEMS}
                if len(years[year])<min_n:ok=False
            if not ok:break
            if not years[2014]&years[2017]:
                ok=False;break
            states[zone]=years
        if ok:out.append({"source_site":site,
             "station_or_between":"between" if site.startswith("Between ") else "station",
             "states":states})
    return out

def turnover_stats(cohort_sites,values,reference=REF):
    if not cohort_sites:
        raise ValueError("Source cohort has no complete repeated spatial sites")
    x=[]
    for r in cohort_sites:
        state=r["states"]
        comp={}
        for zone in ("Track",reference):
            t14={s:values[s] for s in state[zone][2014]}
            t17={s:values[s] for s in state[zone][2017]}
            comp[zone]=components(t14,t17)
            if comp[zone] is None:
                raise AssertionError("No persistent anchor after fixed cohort selection")
        a=comp["Track"];b=comp[reference]
        x.append({
            "reference_gain_minus_track_gain":b["gain_component"]-a["gain_component"],
            "reference_loss_minus_track_loss":b["loss_component"]-a["loss_component"],
            "reference_total_minus_track_total":b["total_cwm_change"]-a["total_cwm_change"],
            "track_loss_component":a["loss_component"],
            "track_gain_component":a["gain_component"],
            "reference_gain_component":b["gain_component"],
            "reference_loss_component":b["loss_component"]})
    fields=x[0].keys()
    out={f:mean(p[f] for p in x) for f in fields}
    if not math.isclose(
        out["reference_total_minus_track_total"],
        out["reference_gain_minus_track_gain"]+out["reference_loss_minus_track_loss"],
        rel_tol=0,abs_tol=1e-10):
        raise AssertionError("Site-balanced temporal decomposition additive check fails")
    return out

def permuted_values(names,base,groups,rng,mode):
    assigned=dict(base)
    if mode=="global":
        vv=[base[s] for s in names]
        rng.shuffle(vv)
        assigned.update(zip(names,vv))
    elif mode=="within_genus":
        for genus,taxa in sorted(groups.items()):
            if len(taxa)<2:continue
            vv=[base[s] for s in taxa]
            rng.shuffle(vv)
            assigned.update(zip(taxa,vv))
    else:
        raise ValueError("Unknown original-source conditional null")
    return assigned

def null(cohort_sites,pool,base,mode,draws=DRAWS,seed=SEED,reference=REF):
    groups=defaultdict(list)
    for sp in pool:groups[sp.split()[0]].append(sp)
    multi={k:v for k,v in groups.items() if len(v)>1}
    obs=turnover_stats(cohort_sites,base,reference)
    keys=("reference_gain_minus_track_gain",
          "reference_total_minus_track_total")
    vec={f:[] for f in keys}
    rng=random.Random(seed+(0 if mode=="global" else 500))
    for _ in range(draws):
        assigned=permuted_values(pool,base,groups,rng,mode)
        st=turnover_stats(cohort_sites,assigned,reference)
        for k in keys:vec[k].append(st[k])
    results={}
    for key in keys:
        v=sorted(vec[key]);observed=obs[key]
        k=sum(abs(x)>=abs(observed)-1e-12 for x in v)
        results[key]={"observed":observed,
                      "null_mean":mean(v),
                      "null_95pct":[v[int(.025*draws)],v[int(.975*draws)]],
                      "technical_two_sided_fraction_plus_one":(k+1)/(draws+1),
                      "n_draws":draws}
    return {"mode":mode,"taxa_pool":len(pool),
            "genera":len(groups),
            "permutable_genera":len(multi),
            "taxa_in_permutable_genera":sum(map(len,multi.values())),
            "n_source_sites":len(cohort_sites),
            "source_site_strata":{
                x:sum(s["station_or_between"]==x for s in cohort_sites)
                for x in ("station","between")},
            "source_observed_components":obs,
            "null_statistics":results,
            "warning":"Conditioned on exact recorded taxon identities and fitted EIVE values, NOT independent regional climate replication"}

def analyze(kashmir,eivefile):
    plots,qa=parse_xlsx_matrix(Path(kashmir))
    names=species_names_from_publisher(kashmir)
    eive,eiveqa=eive_table(eivefile)
    pool=validated_source_taxa(plots,names,eive)
    base={s:float(eive[s][METRIC]) for s in pool}
    if len(pool)<30:raise ValueError("Too few true Kashmir source taxa with EIVE >=3")
    results=[]
    for minimum in (1,3):
        subset=site_cohorts(plots,names,eive,REF,minimum)
        for mode in ("global","within_genus"):
            results.append({"minimum_matched_EIVE_taxa_each_year_and_zone":minimum,
                            **null(subset,pool,base,mode)})
    return {"version":"v330","status":"source-authenticated exploratory taxon-label permutation",
        "original_kashmir":qa,"original_eive":eiveqa,
        "source_EIVE_min_contributing_systems":SYSTEMS,
        "source_years":YEAR,"source_reference":REF,
        "results":results,
        "limitations":[
            "Taxon label permutations do not randomize railway treatment or observed plant transport",
            "Global trait-label permutation breaks phylogenetic association; within-genus null has fewer movable taxa",
            "Reference interior zone is at most 150m from rail and affected by landscape processes",
            "Inference does not correct exploration across temporal filters, references and source trait systems",
            "2017 original railway track has sparse trait-known communities; min1 source panel deliberately descriptive",
            "EIVE European niche indicators are not independent Himalayan physiology measurements",
            "Repeated original source spatial sites belong to one geographically autocorrelated railway corridor",
            "Taxon disappearance is recorded absence, not mortality or dispersal",
            "No original field matrices, site coordinates or species-by-site source observations redistributed"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=analyze(a.kashmir_original,a.eive_original)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({"version":"v330","results":result["results"]},
         ensure_ascii=False),flush=True)
if __name__=="__main__":main()

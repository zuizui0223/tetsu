#!/usr/bin/env python3
"""v322: partition Slovenian railway BS-BT thermal CWM differences by genus.

Two independent publisher sources (Slovenian railway vegetation & European
EIVE thermal traits); both habitats are railway. Do not infer physical train
dispersal, measured temperatures, natural-habitat climate offsets, adaptation
or independent climate-region replication.
"""
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_thermal_generalist_section_v319 import (
    make_pairs, section_test, section_drop)
from analyze_Slovenia_EIVE_temperature_width_v321 import (
    eive_table, latin_map, compare)

SOURCE_EIVE_SYSTEMS=(1,3)
COVERAGE_THRESHOLDS=(0.0,0.5)
METRICS=("position","width")
NULL_DRAWS=1999

def mean(xs):
    return sum(xs)/len(xs) if xs else None

def usable_traits(pairs, names, eive, metric, minsystems):
    """Unique accepted names avoid duplicated EPPO aliases overweighting taxa."""
    observed=set().union(*(p["BS"]|p["BT"] for p in pairs))
    matched={}
    for code in observed:
        name=names.get(code)
        if not name or name not in eive or eive[name]["systems"] < minsystems:
            continue
        matched[name]=float(eive[name][metric])
    genuses=defaultdict(list)
    for name,value in matched.items():
        genuses[name.split()[0]].append(value)
    genus_mean={g:mean(vals) for g,vals in genuses.items()}
    groups=defaultdict(list)
    for name in sorted(matched):
        groups[name.split()[0]].append(name)
    return matched,genus_mean,groups

def observed_zone(codes, code_to_name, matched):
    found={code_to_name[c] for c in codes
           if c in code_to_name and code_to_name[c] in matched}
    return found, len(found)/len(codes) if codes else 0.0

def contrasts(pairs,names,eive,metric,minsystems,coverage):
    matched,gmean,groups=usable_traits(pairs,names,eive,metric,minsystems)
    by_site=defaultdict(list)
    # Each paired side contributes equally within a field site;
    # each site then contributes equally to the overall mean.
    for p in pairs:
        bs,cb=observed_zone(p["BS"],names,matched)
        bt,ct=observed_zone(p["BT"],names,matched)
        if not bs or not bt or min(cb,ct)+1e-12<coverage:
            continue
        get=lambda subset,metric_map:mean([metric_map[x] for x in subset])
        total=get(bs,matched)-get(bt,matched)
        between=(mean([gmean[x.split()[0]] for x in bs])
                -mean([gmean[x.split()[0]] for x in bt]))
        residual=total-between
        by_site[(str(p["section"]),str(p["site"]))].append(
            (total,between,residual,bs,bt))
    if not by_site:
        raise ValueError("No eligible field sites after coverage threshold")
    site_vals={}
    sections_by_component=[defaultdict(list) for _ in range(3)]
    for (sec,site),rows in sorted(by_site.items()):
        val=[mean([row[i] for row in rows]) for i in range(3)]
        if not math.isclose(val[0],val[1]+val[2],abs_tol=1e-12):
            raise AssertionError("Genus decomposition fails additivity")
        site_vals[(sec,site)]=val
        for k in range(3):
            sections_by_component[k][sec].append(val[k])
    allvals=[mean([v[k] for v in site_vals.values()]) for k in range(3)]
    assert math.isclose(allvals[0],allvals[1]+allvals[2],abs_tol=1e-12)

    # v321 independent-source CWM baseline must agree on identical raw taxa.
    original=compare(pairs,names,eive,metric,minsystems,coverage)
    if original["eligible_distinct_field_locations"]!=len(by_site):
        raise ValueError("v321 matched-site cross-check n mismatch")
    if not math.isclose(
        original["field_location_mean_BS_minus_BT"],allvals[0],
        rel_tol=0,abs_tol=1e-10):
        raise ValueError("v321 thermal CWM baseline changed: check taxon crosswalk")
    return {
       "metric":metric,"min_EIVE_source_systems":minsystems,
       "min_numeric_coverage_each_habitat":coverage,
       "matched_unique_eive_taxa":len(matched),
       "represented_genera":len(groups),
       "genera_with_at_least_two_species":sum(len(x)>1 for x in groups.values()),
       "species_in_multi_species_genera":sum(len(x) for x in groups.values() if len(x)>1),
       "eligible_track_sides":sum(map(len,by_site.values())),
       "eligible_field_sites":len(site_vals),
       "eligible_line_sections":len(sections_by_component[0]),
       "total_location_mean_delta":allvals[0],
       "between_genus_location_mean_delta":allvals[1],
       "within_genus_trait_residual_location_mean_delta":allvals[2],
       "section_direction":{
         k:{"exact_signed_flip":section_test(sections_by_component[i]),
            "leave_section_out":section_drop(sections_by_component[i])}
         for i,k in enumerate(("total","between_genus","within_genus_residual"))},
       "_raw":{"matched":matched,"gmean":gmean,"groups":groups,"sites":by_site,
               "section_effects":sections_by_component}}

def within_genus_permutation(result,n_draws=NULL_DRAWS,seed=322):
    src=result["_raw"]
    matched=src["matched"]
    gmean=src["gmean"]
    groups=src["groups"]
    sites=src["sites"]
    # Coefficient representation reproduces equal-site, equal-side CWMs.
    coeff=defaultdict(float)
    for site_rows in sites.values():
        for (_,_,_,bs,bt) in site_rows:
            weight=1/(len(sites)*len(site_rows))
            for taxon in bs:
                coeff[taxon]+=weight/len(bs)
            for taxon in bt:
                coeff[taxon]-=weight/len(bt)
    observed=sum(coeff[name]*(v-gmean[name.split()[0]])
                 for name,v in matched.items())
    if not math.isclose(
         observed,result["within_genus_trait_residual_location_mean_delta"],
         rel_tol=0,abs_tol=1e-10):
        raise AssertionError("Inconsistent within-genus linear coefficients")
    rng=random.Random(seed)
    null=[]
    items=[(genus,names) for genus,names in sorted(groups.items()) if len(names)>1]
    # Taxa belonging to singleton genera have zero residuals and are fixed.
    for _ in range(n_draws):
        stat=0.0
        for g,taxa in items:
            values=[matched[t] for t in taxa]
            rng.shuffle(values)
            baseline=gmean[g]
            stat+=sum(coeff[t]*(v-baseline) for t,v in zip(taxa,values))
        null.append(stat)
    null.sort()
    extreme=sum(abs(v)>=abs(observed)-1e-12 for v in null)
    return {"null_model":"permute EIVE values between source-observed taxa within the SAME genus",
        "statistic":"field-site equal-weight BS minus BT genus-centred EIVE mean",
        "n_draws":n_draws,"seed":seed,
        "actual_within_genus_component":observed,
        "conditional_null_mean":mean(null),
        "conditional_null_95pct":[null[int(.025*n_draws)],null[int(.975*n_draws)]],
        "conditional_two_sided_fraction_with_plus_one":(extreme+1)/(n_draws+1),
        "n_genera_permuted":len(items),
        "n_source_species_in_permutable_genera":sum(len(t) for _,t in items),
        "caveat":"Exploratory taxon-trait exchangeability null; not random field treatment or phylogenetic species-selection proof"}

def analyze(original,eivepath):
    pairs,_,source_qa,_=make_pairs(original)
    names=latin_map(original)
    eive,meta=eive_table(eivepath)
    panels=[]
    for metric in METRICS:
        for systems in SOURCE_EIVE_SYSTEMS:
            for coverage in COVERAGE_THRESHOLDS:
                stats=contrasts(pairs,names,eive,metric,systems,coverage)
                perms=within_genus_permutation(stats)
                stats.pop("_raw")
                stats["genus_trait_shuffle"]=perms
                panels.append(stats)
    return {"version":"v322","status":"REAL TWO ORIGINAL SOURCE REPLAY",
      "source_slovenia_DOI":"10.1007/s44473-026-00155-z",
      "source_EIVE_DOI":"10.5281/zenodo.7534792",
      "source_eive_metadata":meta,
      "source_rail_field_sites":len({p["site"] for p in pairs}),
      "source_rail_sections":len({p["section"] for p in pairs}),
      "source_side_pairs":len(pairs),
      "method":"Exact additive EIVE CWM partition using genus mean over source-observed matched taxa and within-genus species residual, then equal field-location BS-BT difference",
      "panels":panels,
      "limitations":[
        "Both BS and BT are active-railway zones, NOT railway versus background natural vegetation",
        "Genus-centred residual can be affected by changing genus relative abundances, even without strict within-genus replacement; do not claim this identifies only within-genus turnover",
        "Genus taxon pool is defined by species observed in these railway zones and is not the complete regional species pool",
        "Genus-level trait partition does not demonstrate genetic adaptation or phylogenetic causality",
        "Original vegetation survey is nonrandom and physical sampled area is not identical between zones",
        "Exact section sign-flips and within-genus trait permutations are exploratory conditional nulls, not randomized treatment p-values",
        "Many EIVE niche-width values are model-derived; evaluate >=3 source systems separately",
        "A macroclimatic railway-vs-nonrail offset cannot be estimated from this within-rail study"]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-slovenia-xlsx",required=True)
    ap.add_argument("--verified-v102",required=True)
    ap.add_argument("--eive-xlsx",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    original=extract_original_xlsx(a.source_slovenia_xlsx,a.verified_v102)
    out=analyze(original,a.eive_xlsx)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf8")
    compact=[{k:p.get(k) for k in (
         "metric","min_EIVE_source_systems","min_numeric_coverage_each_habitat",
         "matched_unique_eive_taxa","represented_genera",
         "genera_with_at_least_two_species","eligible_field_sites",
         "total_location_mean_delta","between_genus_location_mean_delta",
         "within_genus_trait_residual_location_mean_delta",
         "genus_trait_shuffle")} for p in out["panels"]]
    print(json.dumps({"version":out["version"],"panels":compact},
                     ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

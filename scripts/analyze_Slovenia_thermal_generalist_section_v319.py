#!/usr/bin/env python3
"""v319 source-backed Slovenian railway temperature indicator robustness.

Original source: Sabic et al. 2026 401279-byte publisher XLSX previously
SHA256-verified by the inherited v103 parser. BS and BT are both *railway*
habitats. No matched nonrail plot, causal thermal exposure or train-seed
transport inference. Source data never republished.
"""
import argparse
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_SHA256, extract_original_xlsx, validate_original_keys
)

def eiv(v):
    if isinstance(v,(float,int)) and not isinstance(v,bool) and math.isfinite(v):
        if not 1 <= v <= 12:
            raise ValueError("Source EIV_T outside plausible 1-12 scale")
        return float(v)
    if v is None or str(v).strip().lower() in ("", "x"):
        return None
    try:
        n=float(v)
    except (ValueError,TypeError):
        raise ValueError("Unknown EIV_T token; stop instead of coercing") from None
    if not math.isfinite(n) or not 1<=n<=12:
        raise ValueError("Invalid EIV_T numeric value")
    return n

def trait_class(v):
    if v is None or str(v).strip()=="":
        return "blank"
    if str(v).strip().lower()=="x":
        return "generalist_x"
    return "numeric" if eiv(v) is not None else "blank"

def make_pairs(original):
    geo,zones,incidences,traits,_,qa=validate_original_keys(original)
    vals=[r.get("EIV_T") for r in traits.values()]
    n_num=sum(trait_class(v)=="numeric" for v in vals)
    n_x=sum(trait_class(v)=="generalist_x" for v in vals)
    n_blank=sum(trait_class(v)=="blank" for v in vals)
    if (len(vals),n_num,n_x,n_blank)!=(312,142,135,35):
        raise ValueError("Original EIV_T source coverage drift")
    pairs=[]
    for (site, zid),z in zones.items():
        if z["area"]!="BS":
            continue
        if zid not in ("LBS","RBS"):
            raise ValueError("Unexpected original shoulder zone")
        toe=zid[0]+"BT"
        zz=zones.get((site,toe))
        if zz is None:
            continue
        if zz["area"]!="BT" or z.get("area_side")!=zz.get("area_side"):
            raise ValueError("Original site-side zone mismatch")
        pairs.append({"site":site, "side":zid[0],
                      "section":str(geo[site]["track_section_number"]),
                      "BS":incidences[(site,zid)], "BT":incidences[(site,toe)]})
    if len(pairs)!=173 or len({x["site"] for x in pairs})!=89:
        raise ValueError("Original 173 pairs / 89 field sites mismatch")
    return pairs,traits,qa,{"numeric":n_num,"generalist_x":n_x,"blank":n_blank}

def zone_measure(codes, traits, excluded=None):
    codes=(codes-{excluded}) if excluded else codes
    if not codes:
        return {"value":None,"coverage":0.,"x":0.,"blank":0.,"n":0}
    ev=[eiv(traits[c].get("EIV_T")) for c in codes]
    known=[v for v in ev if v is not None]
    counts=[trait_class(traits[c].get("EIV_T")) for c in codes]
    return {"value":sum(known)/len(known) if known else None,
            "coverage":len(known)/len(codes),
            "x":counts.count("generalist_x")/len(codes),
            "blank":counts.count("blank")/len(codes), "n":len(codes)}

def field_summary(pairs,traits,min_coverage=0.,excluded=None):
    bysite=defaultdict(list)
    for p in pairs:
        bs=zone_measure(p["BS"],traits,excluded)
        bt=zone_measure(p["BT"],traits,excluded)
        if bs["value"] is None or bt["value"] is None:
            continue
        if min(bs["coverage"],bt["coverage"])+1e-12 < min_coverage:
            continue
        bysite[(p["section"],p["site"])].append((
            bs["value"]-bt["value"],bs["x"]-bt["x"],
            bs["blank"]-bt["blank"],bs["coverage"]-bt["coverage"]))
    sites={k:[sum(v[j] for v in rows)/len(rows) for j in range(4)]
           for k,rows in bysite.items()}
    sections=defaultdict(list)
    for (section,_),values in sites.items():
        sections[section].append(values[0])
    vals=[x[0] for x in sites.values()]
    return {"n_sites":len(sites), "n_sides":sum(map(len,bysite.values())),
            "n_sections":len(sections),
            "mean_eiv_difference":sum(vals)/len(vals) if vals else None,
            "bs_greater_sites":sum(v>0 for v in vals),
            "bt_greater_sites":sum(v<0 for v in vals),
            "mean_generalist_x_fraction_difference":(
                sum(x[1] for x in sites.values())/len(sites) if sites else None),
            "mean_blank_fraction_difference":(
                sum(x[2] for x in sites.values())/len(sites) if sites else None),
            "mean_numeric_coverage_difference":(
                sum(x[3] for x in sites.values())/len(sites) if sites else None),
            "_sections":dict(sections),"_sites":sites}

def section_test(groups):
    """Exact section-level sign reversal: equal original *section* weight.

    This is an exploratory randomization diagnostic conditional on survey
    sections, not an independent-field experiment p-value.
    """
    effects=[sum(v)/len(v) for _,v in sorted(groups.items())]
    n=len(effects)
    if not n:
        return None
    if n>18:
        raise ValueError("Excessive section combinations")
    obs=sum(effects)/n
    k=0
    for signs in itertools.product((-1,1),repeat=n):
        test=sum(x*s for x,s in zip(effects,signs))/n
        if abs(test)>=abs(obs)-1e-12:
            k+=1
    return {"section_count":n,"positive_section_count":sum(x>0 for x in effects),
            "negative_section_count":sum(x<0 for x in effects),
            "equal_section_mean":obs, "signed_flip_two_sided_fraction":k/(2**n),
            "enumerated_section_sign_flips":2**n}

def section_drop(groups):
    if len(groups)<2:
        return None
    items=sorted(groups)
    means=[]
    for exclude in items:
        vals=[x for s,arr in groups.items() if s!=exclude for x in arr]
        means.append(sum(vals)/len(vals))
    return {"leave_one_section_out_min":min(means),
            "leave_one_section_out_max":max(means),
            "all_leave_one_section_out_positive":all(x>0 for x in means)}

def summarize(s):
    return {k:v for k,v in s.items() if not k.startswith("_")}

def analyze(original):
    pairs,traits,qa,categories=make_pairs(original)
    allsource=field_summary(pairs,traits,0)
    base=allsource["mean_eiv_difference"]
    if not math.isclose(base,0.1292075293885701,rel_tol=0,abs_tol=1e-10):
        raise ValueError("v318 original-source paired EIV baseline drift")
    out={"version":"v319","scope":"Slovenian RAILWAY shoulder vs toe; NO nonrail controls",
         "original_sha256":SOURCE_SHA256,"source_trait_categories":categories,
         "source_qa":qa, "n_original_pairs":len(pairs),
         "n_distinct_locations":len({p["site"] for p in pairs}),
         "sections":len({p["section"] for p in pairs}),
         "thresholds":[],"single_taxon_omit":[]}
    for thresh in (0,0.5,0.7,0.9):
        s=field_summary(pairs,traits,thresh)
        obj=summarize(s)
        obj["min_numeric_trait_coverage_both_zones"]=thresh
        obj["section_equal_weight_sign_flip"]=section_test(s["_sections"])
        obj["leave_one_section_out"]=section_drop(s["_sections"])
        out["thresholds"].append(obj)
    # Exclude one taxon with an actual numeric indicator across all plots,
    # not an arbitrary source "x" generalist or blank-coded taxon.
    covered=sorted(c for c,t in traits.items()
                   if trait_class(t.get("EIV_T"))=="numeric")
    for thresh in (0,0.5):
        estimates=[]
        for code in covered:
            s=field_summary(pairs,traits,thresh,excluded=code)
            estimates.append((code,s["mean_eiv_difference"],s["n_sites"]))
        valid=[(code,val,n) for code,val,n in estimates if val is not None]
        out["single_taxon_omit"].append({
            "threshold":thresh,"eligible_numeric_taxa_considered":len(covered),
            "results_with_estimable_mean":len(valid),
            "min_delta":min(x[1] for x in valid) if valid else None,
            "max_delta":max(x[1] for x in valid) if valid else None,
            "number_nonpositive":sum(x[1]<=0 for x in valid),
            "min_sites_after_omission":min(x[2] for x in valid) if valid else None,
            "max_sites_after_omission":max(x[2] for x in valid) if valid else None})
    out["limits"]=[
        "Both BS and BT are railway habitats; not a railway vs natural habitat test",
        "EIV_T x identifies temperature generalists, and should NOT be recoded to zero or an imputed numeric warmth preference",
        "Exact sign flips assume exchangeability of section-level contrast directions under a no-direction null; observational design is not randomized",
        "Original railway line sections are not independent climatic regions",
        "Taxon omission is exploratory sensitivity, NOT source attribution or a causal mechanism",
        "No regional climate predictor, railway microclimate sensor, train carriage or nonrail plot in source"]
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-xlsx",required=True)
    ap.add_argument("--verified-v102",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.verified_v102)
    out=analyze(original)
    p=Path(a.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    # Print compact output, never source rows
    print(json.dumps({"thresholds":out["thresholds"],
        "single_taxon_omit":out["single_taxon_omit"],
        "source_trait_categories":out["source_trait_categories"]},
         ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

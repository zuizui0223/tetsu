#!/usr/bin/env python3
"""v11.4 REAL global MIREN ROAD flora source: species weighting benchmark.

Original Haider et al. (2026) GEB 10.1111/geb.70284 and actual
MIREN 2007–2023 Zenodo v3 10.5281/zenodo.15749536.
Not a railway dataset. Test general statistical property of guild
composition when each original plant species has equal weight versus
a site's source-recorded occupancy-opportunity weight, comparing
regional Alien vs Native labels at matched equal-area (100 m2)
roadside (Plot1) and interior (Plot3) quadrats. Species labels are
region-specific; no invented within-plot presences/absences, no
railway wind-diaspore category or direct train vector.
Early MIREN protocols often measured aliens more fully than natives,
so 2017+ is a mandatory separate secondary survey-quality control.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
from statistics import mean,median

from analyze_MIREN_standardized_road_edge_vegetation_v99 import (
    SOURCE_MD5,DOI,PAPER_DOI,source_integrity,native_status_safe,
    plot_species_richness,matched_siteyears,SITE,SITE_YEAR
)

KINDS=("Alien","Native")
MIN_SPECIES_PER_GUILD=3
MIN_SITES_PER_REGION=10
BOOT=5000
SEED=20261114


def species_preference_from_site_rows(site_taxa_rows):
    """Original one site-year x species distinct two original plot classes.

    Input rows may repeat a taxon in both plots. A one-sided observation
    is not a temporal colonization, only an observed source plot contrast.
    """
    unions=defaultdict(lambda:set())
    status={}
    for x in site_taxa_rows:
        if x["status"] not in KINDS:continue
        if x["plot"] not in ("1","3"):
            raise ValueError("Not genuine same-site field Plot1/Plot3 source")
        k=(x["region"],x["site"],x["taxon"])
        key=(x["region"],x["taxon"])
        if key in status and status[key]!=x["status"]:
            raise ValueError("Original region-specific taxon Alien/Native status conflicted")
        status[key]=x["status"]
        unions[k].add(x["plot"])
    counts=defaultdict(lambda:[0,0,0])
    for (region,site,taxon),plots in unions.items():
        v=counts[(region,taxon)]
        if plots=={"1"}:v[0]+=1
        elif plots=={"3"}:v[1]+=1
        elif plots=={"1","3"}:v[2]+=1
        else:raise ValueError("Invalid original paired vegetation source plot set")
    byreg=defaultdict(list)
    for (region,taxon),v in sorted(counts.items()):
        present=sum(v)
        if present<=0:raise ValueError("Plot taxon with no observed site")
        byreg[region].append({
            "region":region,"taxon":taxon,"guild":status[(region,taxon)],
            "roadside_only":v[0],"interior_only":v[1],"both":v[2],
            "site_occupancy":present,
            "roadside_minus_interior_preference":(v[0]-v[1])/present})
    return byreg


def region_contrast(rows,min_sites_per_taxon=1):
    """Comparative regional Alien-vs-Native preference two different weights."""
    grp={s:[] for s in KINDS}
    for r in rows:
        if r["site_occupancy"]>=min_sites_per_taxon:
            grp[r["guild"]].append(r)
    counts={kind:len(grp[kind]) for kind in KINDS}
    if min(counts.values())<MIN_SPECIES_PER_GUILD:
        return {"identified":False,"n_taxa_by_source_status":counts}
    means={}
    for kind,items in grp.items():
        means[kind]={
            "n_species":len(items),
            "species_equal_mean":mean(x["roadside_minus_interior_preference"] for x in items),
            "frequency_weighted_mean":sum(
                x["roadside_minus_interior_preference"]*x["site_occupancy"] for x in items
                )/sum(x["site_occupancy"] for x in items),
        }
    equal=means["Alien"]["species_equal_mean"]-means["Native"]["species_equal_mean"]
    weighted=means["Alien"]["frequency_weighted_mean"]-means["Native"]["frequency_weighted_mean"]
    return {
        "identified":True,"n_taxa_by_source_status":counts,
        "n_comparable_source_species":sum(counts.values()),
        "alien_minus_native_species_equal_preference":equal,
        "alien_minus_native_occurrence_weighted_preference":weighted,
        "opposite_signs":equal*weighted<0,
        "means_by_status":means,
    }


def original_source_plots(raw,cohort):
    """Use first confirmed matched 1/3 field survey for EACH road-site."""
    dedup,_=native_status_safe(raw)
    ranks=plot_species_richness(dedup)
    _,paired,_=matched_siteyears(ranks)
    if cohort=="2017plus":
        paired=paired[paired["Year"]>=2017].copy()
    elif cohort!="firstall":
        raise ValueError("Unrecognized MIREN original source cohort")
    selected=(paired.sort_values("Year")
              .drop_duplicates(list(SITE),keep="first"))[list(SITE_YEAR)].copy()
    matched=dedup[dedup["Plot"].isin(["1","3"]) & dedup["Status"].isin(KINDS)]
    matched=matched.merge(selected,on=list(SITE_YEAR),how="inner",validate="many_to_one")
    if matched.duplicated(list(SITE_YEAR)+["Plot","Accepted.Name.MIREN"]).any():
        raise ValueError("Original repeated taxon-plot source records not removed")
    ns=selected.groupby("Region").size().to_dict()
    records=[]
    for v in matched[list(SITE_YEAR)+["Plot","Accepted.Name.MIREN","Status"]].itertuples(index=False,name=None):
        region,road,site,year,plot,taxon,status=v
        # Site ID must include ROAD plus distinct site source-year so no collision.
        records.append({"region":str(region),
                        "site":(str(road),str(site),int(year)),
                        "plot":str(plot),"taxon":str(taxon),"status":str(status)})
    return records,{"n_original_matched_unique_road_sites":len(selected),
                    "n_real_taxon_times_plot_source_rows":len(records),
                    "n_original_region_codes_with_matched_data":len(ns),
                    "n_original_site_sampled_by_region":ns,
                    "original_native_status_source_coverage_limited_pre2012":True}


def summarize_regions(byreg,meta):
    regions=[]
    for region,rows in sorted(byreg.items()):
        ns=meta["n_original_site_sampled_by_region"].get(region,0)
        d=region_contrast(rows)
        d.update({"source_region_code":region,
                  "n_original_matched_source_road_sites":ns})
        for cutoff in (3,5,10):
            x=region_contrast(rows,cutoff)
            d["sensitivity_min_taxon_sites_"+str(cutoff)]={
                "identified":x["identified"],
                "species_equal":x.get("alien_minus_native_species_equal_preference"),
                "occurrence_weighted":x.get("alien_minus_native_occurrence_weighted_preference"),
                "n_taxa_by_guild":x["n_taxa_by_source_status"]}
        regions.append(d)
    return regions


def quantile(xs,p):
    xs=sorted(xs)
    if not xs:return None
    q=(len(xs)-1)*p
    i=int(math.floor(q));j=int(math.ceil(q))
    return xs[i] if i==j else xs[i]+(xs[j]-xs[i])*(q-i)


def macro_summary(regions):
    selected=[r for r in regions if r["identified"] and
              r["n_original_matched_source_road_sites"]>=MIN_SITES_PER_REGION]
    if not selected:
        return {"identified":False,"n_original_regions_eligible":0}
    a=[r["alien_minus_native_species_equal_preference"] for r in selected]
    b=[r["alien_minus_native_occurrence_weighted_preference"] for r in selected]
    orig=(mean(a),mean(b))
    rng=random.Random(SEED)
    ea=[];eb=[];nrev=[]
    for _ in range(BOOT):
        chosen=[rng.choice(selected) for _ in selected]
        eq=mean(x["alien_minus_native_species_equal_preference"] for x in chosen)
        oc=mean(x["alien_minus_native_occurrence_weighted_preference"] for x in chosen)
        ea.append(eq);eb.append(oc);nrev.append(eq*oc<0)
    return {
        "identified":True,
        "n_original_MIREN_road_region_codes_eligible_10plus_sites":len(selected),
        "n_region_codes_alien_relative_roadside_preference_positive_species_equal":
            sum(x>0 for x in a),
        "n_region_codes_alien_relative_roadside_preference_positive_frequency_weighted":
            sum(x>0 for x in b),
        "n_region_codes_with_local_weighting_sign_reversal":sum(x*y<0 for x,y in zip(a,b)),
        "region_equal_species_equal_alien_minus_native_preference":orig[0],
        "region_equal_occupancy_weighted_alien_minus_native_preference":orig[1],
        "original_macro_mean_signs_reversed":orig[0]*orig[1]<0,
        "region_resample_source_species_equal_95pct":[quantile(ea,.025),quantile(ea,.975)],
        "region_resample_source_frequency_weighted_95pct":[quantile(eb,.025),quantile(eb,.975)],
        "fraction_region_resample_gives_mean_sign_reversal":mean(nrev),
        "n_region_bootstrap_draws":BOOT,
        "MIREN_codes_not_claimed_independent_randomized_countries":True,
    }


def run(raw,receipt99):
    if (receipt99.get("complete") is not True or
        receipt99.get("source_original_MIREN_Zenodo_DOI")!=DOI or
        receipt99.get("source",{}).get("independently_reverified_original_Zenodo_MD5")!=SOURCE_MD5):
        raise ValueError("Original source v9.9 MIREN genuine source receipt missing")
    result={"complete":True,"version":"v11.4","original_MIREN_Zenodo_DOI":DOI,
            "original_MIREN_source_MD5":SOURCE_MD5,"original_data_paper_DOI":PAPER_DOI,
            "n_original_vegetation_source_rows":len(raw),
            "source_record_type":"MIREN ROAD 100m2 equal-area plots, NOT RAIL",
            "cohorts":{}}
    for cohort in ("firstall","2017plus"):
        observations,meta=original_source_plots(raw,cohort)
        byreg=species_preference_from_site_rows(observations)
        regional=summarize_regions(byreg,meta)
        macro=macro_summary(regional)
        result["cohorts"][cohort]={
            "source_metadata":meta,
            "n_region_taxon_units_with_real_native_or_alien_status":
                sum(len(v) for v in byreg.values()),
            "n_region_codes_computed":len(byreg),
            "macro_region_equal_result":macro,
            "all_regions_including_missing_source_qualification":regional,
        }
    first=result["cohorts"]["firstall"]["source_metadata"]
    if first["n_original_matched_unique_road_sites"]!=1005:
        raise ValueError("Original MIREN v9.9 1005 independently surveyed road sites missing")
    return result


def main():
    p=argparse.ArgumentParser()
    for n in ("source-csv","source-audit-v98","verified-v99","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    raw,md5=source_integrity(a.source_csv,a.source_audit_v98)
    prior=json.loads(Path(a.verified_v99).read_text(encoding="utf-8"))
    result=run(raw,prior)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v11.4",
        "source_original_MIREN_rows":result["n_original_vegetation_source_rows"],
        "cohorts":{k:{
            "matched_original_road_sites":v["source_metadata"]["n_original_matched_unique_road_sites"],
            "n_true_plant_region_taxon_units":v["n_region_taxon_units_with_real_native_or_alien_status"],
            "macro":v["macro_region_equal_result"],
            "region_preview":v["all_regions_including_missing_source_qualification"][:6],
        }for k,v in result["cohorts"].items()}
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()

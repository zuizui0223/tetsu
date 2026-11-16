#!/usr/bin/env python3
"""v22.4: REAL paired Slovenian ballast shoulder/toe functional breadth.

Original Šabić et al 2026 author XLSX, 89 locations / 173 same-side pairs.
Original LEDA canopy height (m) and seed mass (mg) are species medians,
not measured railway individuals or evidence of within-species adaptation.

Primary: within-community trace of the population covariance matrix of
log(height), log(seed mass), each standardized by unique species SD across
this genuine 89-location species pool. BS-minus-BT < 0 = narrower BS.
Unit: match side -> equal location -> equal 15 original sections.

Null: within SAME location AND side, keep shared species in BOTH zones,
keep exactly n exclusive BS and BT species, permute assignment ONLY of
exclusive source taxa. This fixes each zone's trait-known richness, overlap,
and the local actually observed species pool. It does not equal a separate
nonrail reference community or randomized treatment of actual railway sites.

Sensitivity: richness-rarefy both sides to identical min(nBS,nBT) species.
Technical random draws are not independent field experiments or empirical
95% confidence intervals; section bootstrap is descriptive uncertainty.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import (
    read_leda,binomial,audit_si
)
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx,validate_original_keys
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping

HEIGHT_SHA="cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92"
SEED_SHA="34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef"
PRIMARY_DRAWS=999
RARE_DRAWS=199
BOOT_DRAWS=3999
SEED=20261116
DIMENSIONS=("joint_trace","height_var","seed_var","height_mean","seed_mean")

def quart(v,p):
    v=sorted(v)
    if not v:return None
    k=(len(v)-1)*p
    lo=int(k);hi=min(lo+1,len(v)-1)
    return v[lo]*(1-(k-lo))+v[hi]*(k-lo)

def ztransform(taxa,raw_height,raw_seed):
    observed={taxon:(raw_height[latin],raw_seed[latin])
              for taxon,latin in taxa.items() if latin in raw_height and latin in raw_seed}
    if len(observed)<50:raise ValueError("No representative source joint traits")
    ln={t:(math.log(h),math.log(s)) for t,(h,s) in observed.items()}
    mu=[mean(v[k] for v in ln.values()) for k in (0,1)]
    sd=[math.sqrt(mean((v[k]-mu[k])**2 for v in ln.values())) for k in (0,1)]
    if min(sd)<=0:raise ValueError("No original trait population spread")
    return {t:((v[0]-mu[0])/sd[0],(v[1]-mu[1])/sd[1]) for t,v in ln.items()},{
        "n_unique_source_taxon_codes_joint":len(observed),
        "log_metric_mean_unweighted_source_species":mu,
        "log_metric_sd_unweighted_source_species":sd,
        "original_traits_as_median_across_LEDA_rows":True}

def variance(xs):
    if len(xs)<2:return None
    mu=sum(xs)/len(xs)
    return sum((x-mu)**2 for x in xs)/len(xs)

def components(codes,values):
    if len(codes)<2:raise ValueError("Cannot measure trait breadth for <2 source species")
    a=[values[t][0] for t in codes]
    b=[values[t][1] for t in codes]
    x=variance(a);y=variance(b)
    return dict(joint_trace=x+y,height_var=x,seed_var=y,
                height_mean=mean(a),seed_mean=mean(b))

def difference(a,b,values):
    x=components(a,values);y=components(b,values)
    return {k:x[k]-y[k] for k in DIMENSIONS}

def create_sites(original,values,min_taxa=5):
    loc,zone,sp,traits,codebook,integrity=validate_original_keys(original)
    sites=defaultdict(list);all_p=[]
    for (lid,zid),one in zone.items():
        if one["area"]!="BS" or zid not in ("LBS","RBS"):continue
        other=zone.get((lid,zid[0]+"BT"))
        if other is None or other["area"]!="BT":continue
        if one.get("area_side")!=other.get("area_side"):
            raise ValueError("True source track-side mismatch")
        a=set(sp.get((lid,zid),set()))
        b=set(sp.get((lid,zid[0]+"BT"),set()))
        if not a or not b:continue
        all_p.append((lid,zid,a,b))
        x=a&set(values);y=b&set(values)
        if len(x)<min_taxa or len(y)<min_taxa:continue
        geo=loc[lid]
        sites[str(lid)].append({
            "section":str(geo["track_section_number"]),
            "track_type":str(geo["track_type"]),
            "original_year_serial":geo.get("date"),
            "original_total_BS_n":len(a),"original_total_BT_n":len(b),
            "nBS":len(x),"nBT":len(y),
            "n_shared":len(x&y),"BS":x,"BT":y,
            "n_known_missing_BS":len(a-x),"n_known_missing_BT":len(b-y)
        })
    if len(all_p)!=173 or len({p[0] for p in all_p})!=89:
        raise ValueError("Original 89 sites 173 sides source pair cohort changed")
    if len(sites)<60 or len({p["section"] for rows in sites.values() for p in rows})<10:
        raise ValueError("Fail closed: original source sampling cannot support primary paired trait breadth")
    if any(len({p["section"] for p in rows})!=1 for rows in sites.values()):
        raise ValueError("Site is not nested in exactly one original rail section")
    return sites,all_p

def site_effects(site_pairs,values,draw=None,rg=None,rarefy=False):
    # Return one equal-weight effect per real site; sides are within-site repeats.
    out={}
    for site,prs in sorted(site_pairs.items()):
        v=[]
        for p in prs:
            a,b=p["BS"],p["BT"]
            if draw is not None:
                shared=a&b
                only=sorted(a^b)
                k=len(a-shared)
                # Preserve shared observations, both zone sizes, local species pool.
                assigned=rg.sample(only,k)
                a=shared|set(assigned);b=shared|(set(only)-set(assigned))
                assert len(a)==p["nBS"] and len(b)==p["nBT"] and len(a&b)==p["n_shared"]
            if rarefy:
                k=min(len(a),len(b))
                a=set(rg.sample(sorted(a),k))
                b=set(rg.sample(sorted(b),k))
            v.append(difference(a,b,values))
        out[site]={"section":prs[0]["section"],
                   "delta":{k:mean(q[k] for q in v) for k in DIMENSIONS},
                   "n_matched_sides":len(prs),
                   "mean_BS_species":mean(q["nBS"] for q in prs),
                   "mean_BT_species":mean(q["nBT"] for q in prs),
                   "mean_shared_species":mean(q["n_shared"] for q in prs)}
    return out

def section_estimates(locations):
    sect=defaultdict(list)
    for x in locations.values():sect[x["section"]].append(x["delta"])
    by={k:{m:mean(e[m] for e in x) for m in DIMENSIONS} for k,x in sect.items()}
    return by,{m:mean(v[m] for v in by.values()) for m in DIMENSIONS}

def finite_summary(vals,observed=None):
    res={"mean":mean(vals),"q025":quart(vals,.025),"q975":quart(vals,.975),"draws":len(vals)}
    if observed is not None:
        res["observed_minus_mean"]=observed-res["mean"]
        res["observed_less_equal_count"]=sum(x<=observed for x in vals)
        res["lower_tail_rank_plusone"]=(1+sum(x<=observed for x in vals))/(1+len(vals))
        res["centered_two_sided_rank_plusone"]=(1+sum(abs(x-res["mean"])>=abs(observed-res["mean"]) for x in vals))/(1+len(vals))
    return res

def evaluate(original,raw_height,raw_seed,draws=PRIMARY_DRAWS,raredraws=RARE_DRAWS,bootdraws=BOOT_DRAWS):
    source_qc=audit_si_source(original,raw_height,raw_seed)
    loc,zone,sp,traits,cb,integrity=validate_original_keys(original)
    names,codebook=original_slovenian_EPPO_codebook_binomial_mapping(original)
    n=all_original_matched_taxa(zone,sp)
    exact={t:binomial(names.get(t)) for t in n}
    values,transform=ztransform(exact,raw_height,raw_seed)
    if len(values)!=224 or source_qc["joint_height_seed_mass"]["n_source_taxa_joined"]!=224:
        raise ValueError("v22.0 genuine 224 source taxa LEDA original-breadth support failed")
    sites,pairs=create_sites(original,values)
    actual=site_effects(sites,values)
    bysection,ob=section_estimates(actual)
    use={m:[] for m in DIMENSIONS}
    rg=random.Random(SEED)
    for i in range(draws):
        _,score=section_estimates(site_effects(sites,values,draw=i,rg=rg))
        for m in DIMENSIONS:use[m].append(score[m])
    null={m:finite_summary(use[m],ob[m]) for m in DIMENSIONS}
    rare={m:[] for m in DIMENSIONS}
    rr=random.Random(SEED+91)
    for _ in range(raredraws):
        _,score=section_estimates(site_effects(sites,values,rg=rr,rarefy=True))
        for m in DIMENSIONS:rare[m].append(score[m])
    rarefaction={m:finite_summary(rare[m]) for m in DIMENSIONS}
    brng=random.Random(SEED+222)
    sections=sorted(bysection)
    boot={m:[] for m in DIMENSIONS}
    for _ in range(bootdraws):
        sampled=[bysection[sections[brng.randrange(len(sections))]] for j in sections]
        for m in DIMENSIONS:boot[m].append(mean(x[m] for x in sampled))
    bstrap={m:finite_summary(boot[m]) for m in DIMENSIONS}
    per_section={k:{"n_locations":sum(x["section"]==k for x in actual.values()),
                    **v} for k,v in sorted(bysection.items())}
    cov=[{k:mean(p[k] for p in v) for k in ("nBS","nBT","original_total_BS_n","original_total_BT_n")} for v in sites.values()]
    known_site_count=sum(len(v) for v in sites.values())
    assert known_site_count>=139 and len(actual)>=82
    return {
       "version":"v22.4",
       "status":"SOURCE_VERIFIED_SI_BS_BT_CONTINUOUS_TRAIT_SPREAD",
       "source_real_89_field_sites":89,
       "source_real_173_matched_track_sides":173,
       "source_original_15_line_sections":15,
       "source_original_plant_DOI":"10.1007/s44473-026-00155-z",
       "source_genuine_joint_tracked_taxa":len(values),
       "source_original_LEDA_trait_provenance":source_qc,
       "metric_definition":"population variance log height + population variance log seed mass, each rescaled to z by unique species SD across source 89-location flora; species incidence unweighted; BS-BT",
       "unit_order":"match side, equal location, equal original railway section",
       "sample_support":{
         "qualifying_matched_side_pairs":known_site_count,
         "qualifying_original_sites":len(actual),
         "qualifying_original_sections":len(sections),
         "all_original_matched_sites":89,
         "n_actual_sites_per_original_section":{k:v["n_locations"] for k,v in per_section.items()},
         "mean_BS_joint_species_across_locations":mean(c["nBS"] for c in cov),
         "mean_BT_joint_species_across_locations":mean(c["nBT"] for c in cov),
         "mean_original_BS_all_species":mean(c["original_total_BS_n"] for c in cov),
         "mean_original_BT_all_species":mean(c["original_total_BT_n"] for c in cov),
         "mean_paired_shared_species":mean(x["mean_shared_species"] for x in actual.values()),
       },
       "trait_standardization":transform,
       "observed_section_equal_BS_minus_BT":ob,
       "observed_site_equal_BS_minus_BT":{m:mean(x["delta"][m] for x in actual.values()) for m in DIMENSIONS},
       "observed_by_original_section":per_section,
       "sections_BS_spread_narrower_than_BT":{
           m:sum(x[m]<0 for x in bysection.values()) for m in ("joint_trace","height_var","seed_var")},
       "conditional_local_pool_shuffle_999":null,
       "richness_equalization_rarefaction":rarefaction,
       "section_cluster_bootstrap_DESCRIPTIVE_NOT_FIELD_CI":bstrap,
       "null_preserves":"each exact same-side observed local species pool, overlap both-zone taxa, unique zone trait-known richness, not independent regional flora",
       "limitations":[
         "different BS and BT source sampling area/zone width and imperfect botanical detectability",
         "LEDA species medians cannot estimate in situ plant size plasticity",
         "the section bootstrap is ONE observational study not independent geographic experiments",
         "local species union is already selected by railway habitat, never independent natural habitat pool",
         "randomized species assignment conditional on local pool and shared species is technical null, not randomized railway treatment",
         "no measured freight/passenger departures, official line main-vs-branch or train-induced wind exposure",
         "no measured railway adaptation alleles or genetics",
         "do not derive an artificial statistical p-value from 999 technical null draws"
       ]}
def all_original_matched_taxa(zones,species):
    all_codes=set()
    for (lid,zid),row in zones.items():
        if row["area"]!="BS" or zid not in ("LBS","RBS"):continue
        if (lid,zid[0]+"BT") not in zones:continue
        all_codes.update(species.get((lid,zid),set()))
        all_codes.update(species.get((lid,zid[0]+"BT"),set()))
    return all_codes

def audit_si_source(original,height,seed):
    # Recheck v22.0 source coverage result BEFORE deriving phenotypic spreads.
    from audit_SI_CZ_original_continuous_trait_coverage_v220 import summarize
    loc,zones,sp,traits,book,integrity=validate_original_keys(original)
    names,codebook=original_slovenian_EPPO_codebook_binomial_mapping(original)
    all_codes=all_original_matched_taxa(zones,sp)
    # 285 reported SI source-taxonomy cohort in original v22.0.
    if len(all_codes)!=285:raise ValueError("v22.0 original 285 matched author taxa changed")
    units=[]
    for (lid,zid),r in zones.items():
        if r["area"]!="BS" or zid not in ("LBS","RBS"):continue
        if (lid,zid[0]+"BT") not in zones:continue
        units.extend([{"codes":sp.get((lid,zid),set())},{"codes":sp.get((lid,zid[0]+"BT"),set())}])
    stat,_=summarize({t:names.get(t) for t in all_codes},units,
                     {"height_m":height,"seed_mass_mg":seed})
    if stat["joint_height_seed_mass"]["n_source_taxa_joined"]!=224:
        raise ValueError("v22.0 genuine original LEDA coverage 224 joint taxa not reproducible")
    return {"cohort_source_taxa_n":285,
            "height_joined_n":stat["traits"]["height_m"]["n_source_taxa_joined"],
            "seed_joined_n":stat["traits"]["seed_mass_mg"]["n_source_taxa_joined"],
            "joint_joined_n":224,
            "joint_source_trait_known_incidence_fraction":
                stat["joint_height_seed_mass"]["fraction_incidence_joined"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-xlsx",required=True)
    p.add_argument("--source-audit",required=True)
    p.add_argument("--height",required=True)
    p.add_argument("--seed",required=True)
    p.add_argument("--output",required=True)
    p.add_argument("--draws",type=int,default=PRIMARY_DRAWS)
    p.add_argument("--raredraws",type=int,default=RARE_DRAWS)
    p.add_argument("--bootdraws",type=int,default=BOOT_DRAWS)
    a=p.parse_args()
    source=extract_original_xlsx(a.source_xlsx,a.source_audit)
    if hashlib.sha256(Path(a.height).read_bytes()).hexdigest()!=HEIGHT_SHA:
        raise ValueError("Original LEDA source height SHA drift, not validated")
    if hashlib.sha256(Path(a.seed).read_bytes()).hexdigest()!=SEED_SHA:
        raise ValueError("Original LEDA source seed SHA drift, not validated")
    height,hqa=read_leda(a.height,"height_m")
    seed,sqa=read_leda(a.seed,"seed_mass_mg")
    ans=evaluate(source,height,seed,draws=a.draws,raredraws=a.raredraws,bootdraws=a.bootdraws)
    ans["LEDA_source_SHA"]={"height":HEIGHT_SHA,"seed":SEED_SHA}
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ans,ensure_ascii=False,indent=2),encoding="utf-8")
    summary={k:ans[k] for k in ("status","sample_support",
        "observed_section_equal_BS_minus_BT","sections_BS_spread_narrower_than_BT",
        "conditional_local_pool_shuffle_999","richness_equalization_rarefaction",
        "section_cluster_bootstrap_DESCRIPTIVE_NOT_FIELD_CI")}
    print("V224_BREADTH_SUMMARY="+json.dumps(summary,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
